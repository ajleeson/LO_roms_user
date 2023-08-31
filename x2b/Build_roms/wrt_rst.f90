      MODULE wrt_rst_mod
!
!svn $Id: wrt_rst.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes requested model fields into the RESTART output   !
!  file using the standard NetCDF library or the Parallel-IO (PIO)     !
!  library.                                                            !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE nf_fwrite4d_mod, ONLY : nf_fwrite4d
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_rst
      PRIVATE :: wrt_rst_nf90
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_rst (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_rst.F"
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (RST(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_rst_nf90 (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj)
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) RST(ng)%IOtype
  10        FORMAT (' WRT_RST - Illegal output type, io_type = ',i0)
          END IF
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 96, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE wrt_rst
!
!***********************************************************************
      SUBROUTINE wrt_rst_nf90 (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, i, itrc, status
      integer :: ntmp(1)
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_rst.F"//", wrt_rst_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 131, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      RST(ng)%Rindex=RST(ng)%Rindex+1
      Fcount=RST(ng)%Fcount
      RST(ng)%Nrec(Fcount)=RST(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        RST(ng)%Rindex=MOD(RST(ng)%Rindex-1,2)+1
      END IF
!
!  Write out time-stepping indices.
!
      ntmp(1)=1+MOD((iic(ng)-1)-ntstart(ng),2)
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'nstp',             &
     &                      ntmp, (/RST(ng)%Rindex/), (/1/),            &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 165, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'nrhs',             &
     &                      ntmp, (/RST(ng)%Rindex/), (/1/),            &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 170, MyFile)) RETURN
      ntmp(1)=3-ntmp(1)
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'nnew',             &
     &                      ntmp, (/RST(ng)%Rindex/), (/1/),            &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 176, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'kstp',             &
     &                      kstp(ng:), (/RST(ng)%Rindex/), (/1/),       &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 181, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'krhs',             &
     &                      krhs(ng:), (/RST(ng)%Rindex/), (/1/),       &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 186, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, iNLM, RST(ng)%name, 'knew',             &
     &                      knew(ng:), (/RST(ng)%Rindex/), (/1/),       &
     &                      ncid = RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 191, MyFile)) RETURN
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, RST(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/RST(ng)%Rindex/), (/1/),                  &
     &                      ncid = RST(ng)%ncid,                        &
     &                      varid = RST(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, 201, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      scale=1.0_dp
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % zeta)
      IF (FoundError(status, nf90_noerr, 347, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out RHS of free-surface equation.
!
      scale=1.0_dp
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRzet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % rzeta)
      IF (FoundError(status, nf90_noerr, 368, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRzet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % ubar,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 402, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out RHS of 2D momentum equation in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRu2d),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % rubar,                             &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 424, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu2d)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 3, scale,               &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 457, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out RHS of 2D momentum equation in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRv2d),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, 2, scale,               &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % rvbar,                             &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 480, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv2d)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % u,                                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 513, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out RHS of 3D momentum equation in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRu3d),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % ru,                                &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 536, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu3d)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % v,                                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 568, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out RHS of 3D momentum equation in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRv3d),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % rv,                                &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 591, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv3d)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Tid(itrc),   &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 1, 2, scale,   &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % t(:,:,:,:,itrc))
        IF (FoundError(status, nf90_noerr, 623, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_dp
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idDano),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % rho)
      IF (FoundError(status, nf90_noerr, 644, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvis),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 737, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTdif),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 759, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSdif),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 782, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out turbulent kinetic energy.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idMtke),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % tke,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 807, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtke)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define turbulent kinetic energy time length scale.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite4d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idMtls),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), 1, 2, scale,     &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % gls,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, 828, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtls)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define vertical mixing turbulent length scale.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmLS),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Lscale)
      IF (FoundError(status, nf90_noerr, 848, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmLS)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define turbulent kinetic energy vertical diffusion coefficient.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmKK),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akk)
      IF (FoundError(status, nf90_noerr, 868, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmKK)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define turbulent length scale vertical diffusion coefficient.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVmKP),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akp)
      IF (FoundError(status, nf90_noerr, 889, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVmKP)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, RST(ng)%name, RST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 1138, MyFile)) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), RST(ng)%Rindex
!
  10  FORMAT (/,' WRT_RST_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into restart NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_RST_NF90     - wrote re-start', t40,              &
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
!
      RETURN
      END SUBROUTINE wrt_rst_nf90
      END MODULE wrt_rst_mod
