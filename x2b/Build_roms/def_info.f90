      MODULE def_info_mod
!
!svn $Id: def_info.F 1122 2022-04-13 19:50:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!  On input the NetCDF dimesions IDs is an interger vector as follows: !
!                                                                      !
!    DimIDs( 1) => XI-dimension at RHO-points                          !
!    DimIDs( 2) => XI-dimension at U-points                            !
!    DimIDs( 3) => XI-dimension at V-points                            !
!    DimIDs( 4) => XI-dimension at PSI-points                          !
!    DimIDs( 5) => ETA-dimension at RHO-points                         !
!    DimIDs( 6) => ETA-dimension at U-points                           !
!    DimIDs( 7) => ETA-dimension at V-points                           !
!    DimIDs( 8) => ETA-dimension at PSI-points                         !
!    DimIDs( 9) => S-dimension at RHO-points                           !
!    DimIDs(10) => S-dimension at W-points                             !
!    DimIDs(11) => Number of tracers dimension                         !
!    DimIDs(12) => Unlimited time record dimension                     !
!    DimIDs(13) => Number of stations dimension                        !
!    DimIDs(14) => Boundary dimension                                  !
!    DimIDs(15) => Number of floats dimension                          !
!    DimIDs(16) => Number sediment bed layers dimension                !
!    DimIDs(17) => Dimension 2D water RHO-points                       !
!    DimIDs(18) => Dimension 2D water U-points                         !
!    DimIDs(19) => Dimension 2D water V-points                         !
!    DimIDs(20) => Dimension 3D water RHO-points                       !
!    DimIDs(21) => Dimension 3D water U-points                         !
!    DimIDs(23) => Dimension 3D water W-points                         !
!    DimIDs(24) => Dimension sediment bed water points                 !
!    DimIDs(25) => Number of EcoSim phytoplankton groups               !
!    DimIDs(26) => Number of EcoSim bacteria groups                    !
!    DimIDs(27) => Number of EcoSim DOM groups                         !
!    DimIDs(28) => Number of EcoSim fecal groups                       !
!    DimIDs(29) => Number of state variables                           !
!    DimIDs(30) => Number of 3D variables time levels (2)              !
!    DimIDs(31) => Number of 2D variables time levels (3)              !
!    DimIDs(32) => Number of sediment tracers                          !
!    DimIDs(33) => Number of light spectral bands.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_strings
!
      USE def_dim_mod, ONLY : def_dim
      USE def_var_mod, ONLY : def_var
      USE lbc_mod,     ONLY : lbc_putatt
      USE strings_mod, ONLY : FoundError, join_string
      USE tadv_mod,    ONLY : tadv_putatt
!
      implicit none
!
      INTERFACE def_info
        MODULE PROCEDURE def_info_nf90
      END INTERFACE def_info
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_info_nf90 (ng, model, ncid, ncname, DimIDs)
!***********************************************************************
!                                                                      !
!  This routine defines information variables for the requested NetCDF !
!  file using the standard NetCDF-3 or NetCDF-4 library.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF filename (character)                             !
!     DimIDs   NetCDF dimensions IDs (integer vector of size nDimID)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE mod_netcdf
!
      USE distribute_mod, ONLY : mp_bcasti
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      integer, intent(in) :: DimIDs(nDimID)
!
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: srdim, stadim, status, swdim, trcdim, usrdim
      integer :: ibuffer(2)
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=11 )    :: bryatt, clmatt, frcatt
      character (len=50 )    :: tiling
      character (len=80 )    :: type
      character (len=512)    :: bio_file
      character (len=4096)   :: string
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_info.F"//", def_info_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.GST(ng)%ncid)) THEN
        status=def_dim(ng, model, ncid, ncname, 'Nuser',                &
     &                 Nuser, usrdim)
        IF (FoundError(exit_flag, NoError, 198, MyFile)) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define history global attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',       &
     &                                      TRIM( version),             &
     &                                      TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version',              &
     &                               TRIM(version)
        END IF
!
!  Set tile decomposition global attribute.
!
        WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'file',                &
     &                        TRIM(ncname))
          IF (FoundError(status, nf90_noerr, 238, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define NetCDF format type.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'format',              &
     &                        'netCDF-4/HDF5 file')
          IF (FoundError(status, nf90_noerr, 257, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
        type='CF-1.4, SGRID-0.3'
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'Conventions',         &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 273, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file type global attribute.
!
        IF (ncid.eq.ADM(ng)%ncid) THEN
          type='ROMS/TOMS adjoint history file'
        ELSE IF (ncid.eq.AVG(ng)%ncid) THEN
          type='ROMS/TOMS nonlinear model averages file'
        ELSE IF (ncid.eq.DIA(ng)%ncid) THEN
          type='ROMS/TOMS diagnostics file'
        ELSE IF (ncid.eq.FLT(ng)%ncid) THEN
          type='ROMS/TOMS floats file'
        ELSE IF (ncid.eq.ERR(ng)%ncid) THEN
          type='ROMS/TOMS posterior analysis error covariance matrix'
        ELSE IF (ncid.eq.GST(ng)%ncid) THEN
          type='ROMS/TOMS GST check pointing restart file'
        ELSE IF (ncid.eq.HAR(ng)%ncid) THEN
          type='ROMS/TOMS Least-squared Detiding Harmonics file'
        ELSE IF (ncid.eq.HSS(ng)%ncid) THEN
          type='ROMS/TOMS 4D-Var Hessian eigenvectors file'
        ELSE IF (ncid.eq.HIS(ng)%ncid) THEN
          type='ROMS/TOMS history file'
        ELSE IF (ncid.eq.ITL(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear model initial file'
        ELSE IF (ncid.eq.LCZ(ng)%ncid) THEN
          type='ROMS/TOMS 4D-Var Lanczos vectors file'
        ELSE IF (ncid.eq.LZE(ng)%ncid) THEN
          type='ROMS/TOMS 4D-Var Evolved Lanczos vectors file'
        ELSE IF (ncid.eq.NRM(1,ng)%ncid) THEN
          type='ROMS/TOMS initial conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(2,ng)%ncid) THEN
          type='ROMS/TOMS model error covariance norm file'
        ELSE IF (ncid.eq.NRM(3,ng)%ncid) THEN
         type='ROMS/TOMS boundary conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(4,ng)%ncid) THEN
          type='ROMS/TOMS surface forcing error covariance norm file'
        ELSE IF (ncid.eq.QCK(ng)%ncid) THEN
          type='ROMS/TOMS quicksave file'
        ELSE IF (ncid.eq.RST(ng)%ncid) THEN
          type='ROMS/TOMS restart file'
        ELSE IF (ncid.eq.STA(ng)%ncid) THEN
          type='ROMS/TOMS station file'
        ELSE IF (ncid.eq.TLF(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear impulse forcing file'
        ELSE IF (ncid.eq.TLM(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear history file'
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'type',                &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 346, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define other global attributes to NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'title',               &
     &                        TRIM(title))
          IF (FoundError(status, nf90_noerr, 372, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'var_info',            &
     &                        TRIM(varname))
          IF (FoundError(status, nf90_noerr, 382, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'var_info', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'rst_file',            &
     &                        TRIM(RST(ng)%name))
          IF (FoundError(status, nf90_noerr, 417, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'his_base',        &
     &                            TRIM(HIS(ng)%base))
            ELSE
              status=nf90_put_att(ncid, nf90_global, 'his_file',        &
     &                            TRIM(HIS(ng)%name))
            END IF
            IF (FoundError(status, nf90_noerr, 433, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'grd_file',            &
     &                        TRIM(GRD(ng)%name))
          IF (FoundError(status, nf90_noerr, 532, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'grd_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'ini_file',            &
     &                        TRIM(INI(ng)%name))
          IF (FoundError(status, nf90_noerr, 545, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'ini_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LuvSrc(ng).or.LwSrc(ng).or.(ANY(LtracerSrc(:,ng)))) THEN
            status=nf90_put_att(ncid, nf90_global, 'river_file',        &
     &                          TRIM(SSF(ng)%name))
            IF (FoundError(status, nf90_noerr, 688, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'river_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LprocessTides(ng)) THEN
            status=nf90_put_att(ncid, nf90_global, 'tide_file',         &
     &                          TRIM(TIDE(ng)%name))
            IF (FoundError(status, nf90_noerr, 702, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'tide_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO i=1,nFfiles(ng)
            CALL join_string (FRC(i,ng)%files, FRC(i,ng)%Nfiles,        &
     &                        string, lstr)
            WRITE (frcatt,30) 'frc_file_', i
            status=nf90_put_att(ncid, nf90_global, frcatt,              &
     &                          string(1:lstr))
            IF (FoundError(status, nf90_noerr, 719, MyFile)) THEN
              IF (Master) WRITE (stdout,20) TRIM(frcatt), TRIM(ncname)
              exit_flag=3
              ioerror=status
              EXIT
            END IF
          END DO
        END IF
        IF (ObcData(ng)) THEN
          DO i=1,nBCfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (BRY(i,ng)%files, BRY(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (bryatt,30) 'bry_file_', i
              status=nf90_put_att(ncid, nf90_global, bryatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 737, MyFile)) THEN
                IF (Master) WRITE (stdout,20) TRIM(bryatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
        IF (Lclimatology(ng)) THEN
          DO i=1,nCLMfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (CLM(i,ng)%files, CLM(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (clmatt,30) 'clm_file_', i
              status=nf90_put_att(ncid, nf90_global, clmatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 756, MyFile)) THEN
                IF (Master) WRITE (stdout,20) TRIM(clmatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
        IF (Lnudging(ng)) THEN
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncid, nf90_global, 'nud_file',          &
     &                        TRIM(NUD(ng)%name))
            IF (FoundError(status, nf90_noerr, 771, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'nud_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'script_file',         &
     &                        TRIM(Iname))
          IF (FoundError(status, nf90_noerr, 810, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'script_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'bpar_file',           &
     &                        TRIM(bparnam))
          IF (FoundError(status, nf90_noerr, 834, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'bpar_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  NLM tracer advection scheme.
!
        IF (exit_flag.eq.NoError) THEN
          CALL tadv_putatt (ng, ncid, ncname, 'NLM_TADV',               &
     &                      Hadvection, Vadvection, status)
          IF (FoundError(status, nf90_noerr, 874, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'NLM_TADV', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  NLM Lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'NLM_LBC', LBC, status)
          IF (FoundError(status, nf90_noerr, 898, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  SVN repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_url',             &
     &                        TRIM(svn_url))
          IF (FoundError(status, nf90_noerr, 950, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_rev',             &
     &                        TRIM(svn_rev))
          IF (FoundError(status, nf90_noerr, 961, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'svn_rev', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'code_dir',            &
     &                        TRIM(Rdir))
          IF (FoundError(status, nf90_noerr, 977, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_dir',          &
     &                        TRIM(Hdir))
          IF (FoundError(status, nf90_noerr, 989, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_file',         &
     &                        TRIM(Hfile))
          IF (FoundError(status, nf90_noerr, 1001, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Attributes describing platform and compiler
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'os',                  &
     &                        TRIM(my_os))
          IF (FoundError(status, nf90_noerr, 1016, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'cpu',                 &
     &                        TRIM(my_cpu))
          IF (FoundError(status, nf90_noerr, 1026, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_system',     &
     &                        TRIM(my_fort))
          IF (FoundError(status, nf90_noerr, 1036, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_system',            &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_command',    &
     &                        TRIM(my_fc))
          IF (FoundError(status, nf90_noerr, 1047, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_command',           &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          lstr=INDEX(my_fflags, 'free')-2
          IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
          status=nf90_put_att(ncid, nf90_global, 'compiler_flags',      &
     &                        my_fflags(1:lstr))
          IF (FoundError(status, nf90_noerr, 1060, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Tiling and history attributes.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tiling',              &
     &                        TRIM(tiling))
          IF (FoundError(status, nf90_noerr, 1072, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'history',             &
     &                        TRIM(history))
          IF (FoundError(status, nf90_noerr, 1082, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Analytical header files used.
!
        IF (exit_flag.eq.NoError) THEN
          CALL join_string (ANANAME, SIZE(ANANAME), string, lstr)
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'ana_file',          &
     &                          string(1:lstr))
            IF (FoundError(status, nf90_noerr, 1096, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Biology model header file used.
!
        IF (exit_flag.eq.NoError) THEN
          DO i=1,512
            bio_file(i:i)='-'
          END DO
          status=nf90_put_att(ncid, nf90_global, 'bio_file',            &
     &                        bio_file)
          IF (FoundError(status, nf90_noerr, 1115, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'bio_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Activated CPP options.
!
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(Coptions)-1
          status=nf90_put_att(ncid, nf90_global, 'CPP_options',         &
     &                        TRIM(Coptions(1:lstr)))
          IF (FoundError(status, nf90_noerr, 1131, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (FoundError(exit_flag, NoError, 1147, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1168, MyFile)) RETURN
      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1175, MyFile)) RETURN
      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1183, MyFile)) RETURN
      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1191, MyFile)) RETURN
      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      WRITE (Vinfo( 3),'(a,a)') 'days since ', TRIM(Rclock%string)
      Vinfo( 4)=TRIM(Rclock%calendar)
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1200, MyFile)) RETURN
      Vinfo( 1)='shuffle'
      Vinfo( 2)='NetCDF-4/HDF5 file format shuffle filer flag'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1226, MyFile)) RETURN
      Vinfo( 1)='deflate'
      Vinfo( 2)='NetCDF-4/HDF5 file format deflate filer flag'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1233, MyFile)) RETURN
      Vinfo( 1)='deflate_level'
      Vinfo( 2)='NetCDF-4/HDF5 file format deflate level parameter'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1240, MyFile)) RETURN
      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1248, MyFile)) RETURN
      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1256, MyFile)) RETURN
      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1266, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1465, MyFile)) RETURN
      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1472, MyFile)) RETURN
      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1479, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1666, MyFile)) RETURN
      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1674, MyFile)) RETURN
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1684, MyFile)) RETURN
      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1693, MyFile)) RETURN
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1750, MyFile)) RETURN
      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1757, MyFile)) RETURN
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1766, MyFile)) RETURN
      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1774, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1785, MyFile)) RETURN
      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1792, MyFile)) RETURN
      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1799, MyFile)) RETURN
      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1806, MyFile)) RETURN
      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1813, MyFile)) RETURN
      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1820, MyFile)) RETURN
      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1827, MyFile)) RETURN
      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1834, MyFile)) RETURN
      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1841, MyFile)) RETURN
      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1848, MyFile)) RETURN
      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1855, MyFile)) RETURN
      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1862, MyFile)) RETURN
      Vinfo( 1)='Charnok_alpha'
      Vinfo( 2)='Charnok factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1869, MyFile)) RETURN
      Vinfo( 1)='Zos_hsig_alpha'
      Vinfo( 2)='wave amplitude factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1876, MyFile)) RETURN
      Vinfo( 1)='sz_alpha'
      Vinfo( 2)='surface flux from wave dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1883, MyFile)) RETURN
      Vinfo( 1)='CrgBan_cw'
      Vinfo( 2)='surface flux due to Craig and Banner wave breaking'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1890, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1902, MyFile)) RETURN
      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1910, MyFile)) RETURN
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1919, MyFile)) RETURN
      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1927, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1941, MyFile)) RETURN
        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1949, MyFile)) RETURN
        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1957, MyFile)) RETURN
        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1965, MyFile)) RETURN
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1974, MyFile)) RETURN
        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1982, MyFile)) RETURN
        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1990, MyFile)) RETURN
        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 1998, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2011, MyFile)) RETURN
!
!  Various parameters.
!
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2076, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      Vinfo( 1)='LuvSrc'
      Vinfo( 2)='momentum point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2089, MyFile)) RETURN
      Vinfo( 1)='LwSrc'
      Vinfo( 2)='mass point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2098, MyFile)) RETURN
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      Vinfo( 1)='LtracerSrc'
      Vinfo( 2)='tracer point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2112, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      Vinfo( 1)='LsshCLM'
      Vinfo( 2)='sea surface height climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2124, MyFile)) RETURN
      Vinfo( 1)='Lm2CLM'
      Vinfo( 2)='2D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2133, MyFile)) RETURN
      Vinfo( 1)='Lm3CLM'
      Vinfo( 2)='3D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2143, MyFile)) RETURN
      Vinfo( 1)='LtracerCLM'
      Vinfo( 2)='tracer climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2152, MyFile)) RETURN
!
!  Logical switches for nudging of climatology fields.
!
      Vinfo( 1)='LnudgeM2CLM'
      Vinfo( 2)='2D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2164, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeM3CLM'
      Vinfo( 2)='3D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2174, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeTCLM'
      Vinfo( 2)='tracer climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2183, MyFile)) RETURN
!
!  Define Fennel et al. (2006) ecosystem model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 23, MyFile)) RETURN
      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation by seawater'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 31, MyFile)) RETURN
      Vinfo( 1)='AttChl'
      Vinfo( 2)='light attenuation by chlorophyll'
      Vinfo( 3)='meter-2 milligram_Chl-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 39, MyFile)) RETURN
      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 46, MyFile)) RETURN
      Vinfo( 1)='Vp0'
      Vinfo( 2)='Eppley temperature-limited growth parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 53, MyFile)) RETURN
      Vinfo( 1)='I_thNH4'
      Vinfo( 2)='radiation threshold for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 61, MyFile)) RETURN
      Vinfo( 1)='D_p5NH4'
      Vinfo( 2)='half-saturation radiation for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 69, MyFile)) RETURN
      Vinfo( 1)='NitriR'
      Vinfo( 2)='nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 77, MyFile)) RETURN
      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 85, MyFile)) RETURN
      Vinfo( 1)='K_NH4'
      Vinfo( 2)='inverse half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 93, MyFile)) RETURN
      Vinfo( 1)='K_PO4'
      Vinfo( 2)='inverse half-saturation for phytoplankton PO4 uptake'
      Vinfo( 3)='meter3 millimole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 101, MyFile)) RETURN
      Vinfo( 1)='K_Phy'
      Vinfo( 2)='zooplankton half-saturation constant for ingestion'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 109, MyFile)) RETURN
      Vinfo( 1)='Chl2C_m'
      Vinfo( 2)='maximum chlorophyll to carbon ratio'
      Vinfo( 3)='milligram_chl milligram_carbon-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 117, MyFile)) RETURN
      Vinfo( 1)='ChlMin'
      Vinfo( 2)='minimum chlorophyll threshold'
      Vinfo( 3)='milligram_chl meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 125, MyFile)) RETURN
      Vinfo( 1)='PhyCN'
      Vinfo( 2)='phytoplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 133, MyFile)) RETURN
      Vinfo( 1)='R_P2N'
      Vinfo( 2)='phytoplankton Phosphate:Nitrogen ratio'
      Vinfo( 3)='mole_P mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 141, MyFile)) RETURN
      Vinfo( 1)='PhyIP'
      Vinfo( 2)='phytoplankton, NH4 inhibition parameter'
      Vinfo( 3)='millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 149, MyFile)) RETURN
      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton, initial slope of P-I curve'
      Vinfo( 3)='milligram_C milligram_Chl-1 watt-1 meter2 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 157, MyFile)) RETURN
      Vinfo( 1)='PhyMin'
      Vinfo( 2)='minimum phytoplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 165, MyFile)) RETURN
      Vinfo( 1)='PhyMR'
      Vinfo( 2)='phytoplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 173, MyFile)) RETURN
      Vinfo( 1)='ZooAE_N'
      Vinfo( 2)='zooplankton mitrogen assimilation efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 180, MyFile)) RETURN
      Vinfo( 1)='ZooBM'
      Vinfo( 2)='rate for zooplankton basal metabolism'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 188, MyFile)) RETURN
      Vinfo( 1)='ZooCN'
      Vinfo( 2)='zooplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 196, MyFile)) RETURN
      Vinfo( 1)='ZooER'
      Vinfo( 2)='zooplankton specific excretion rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 204, MyFile)) RETURN
      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton maximum growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 212, MyFile)) RETURN
      Vinfo( 1)='ZooMin'
      Vinfo( 2)='minimum zooplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 220, MyFile)) RETURN
      Vinfo( 1)='ZooMR'
      Vinfo( 2)='zooplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 228, MyFile)) RETURN
      Vinfo( 1)='LDeRRN'
      Vinfo( 2)='rate of large detritus nitrogen re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 236, MyFile)) RETURN
      Vinfo( 1)='LDeRRC'
      Vinfo( 2)='rate of large detritus carbon re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 244, MyFile)) RETURN
      Vinfo( 1)='CoagR'
      Vinfo( 2)='coagulation rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 252, MyFile)) RETURN
      Vinfo( 1)='SDeRRN'
      Vinfo( 2)='remineralization rate for small detritus nitrogen'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 260, MyFile)) RETURN
      Vinfo( 1)='SDeRRC'
      Vinfo( 2)='remineralization rate for small detritus carbon'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 268, MyFile)) RETURN
      Vinfo( 1)='RDeRRN'
      Vinfo( 2)='remineralization rate for river detritus nitrogen'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 276, MyFile)) RETURN
      Vinfo( 1)='RDeRRC'
      Vinfo( 2)='remineralization rate for river detritus carbon'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 284, MyFile)) RETURN
      Vinfo( 1)='wPhy'
      Vinfo( 2)='vertical sinking velocity for phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 292, MyFile)) RETURN
      Vinfo( 1)='wLDet'
      Vinfo( 2)='vertical sinking velocity for large detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 300, MyFile)) RETURN
      Vinfo( 1)='wSDet'
      Vinfo( 2)='vertical sinking velocity for small detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 308, MyFile)) RETURN
      Vinfo( 1)='pCO2air'
      Vinfo( 2)='partial pressure of CO2 in the air'
      Vinfo( 3)='parts per million by volume'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 316, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2876, MyFile)) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2886, MyFile)) RETURN
      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2894, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2904, MyFile)) RETURN
      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2911, MyFile)) RETURN
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2918, MyFile)) RETURN
      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2925, MyFile)) RETURN
      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2933, MyFile)) RETURN
      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2941, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids.
!
      Vinfo( 1)='grid'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2949, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2973, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2998, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3011, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3024, MyFile)) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 3037, MyFile)) RETURN
      END IF
      IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bathymetry'
        Vinfo(21)='sea_floor_depth'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3078, MyFile)) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3082, MyFile)) RETURN
        END IF
!
!  Coriolis Parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='Coriolis parameter'
          Vinfo(21)='coriolis_parameter'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3098, MyFile)) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm'
          Vinfo(21)='inverse_grid_x_spacing'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3113, MyFile)) RETURN
          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn'
          Vinfo(21)='inverse_grid_y_spacing'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3124, MyFile)) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='longitude'
          Vinfo(21)='grid_longitude_at_cell_center'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3138, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3142, MyFile)) RETURN
          END IF
          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='latitude'
          Vinfo(21)='grid_latitude_at_cell_center'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (FoundError(exit_flag, NoError, 3153, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3157, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Xr'
          Vinfo(21)='grid_x_location_at_cell_center'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3168, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3172, MyFile)) RETURN
          END IF
          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Yr'
          Vinfo(21)='grid_y_location_at_cell_center'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3183, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3187, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_u'
          Vinfo( 2)='longitude of U-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='longitude'
          Vinfo(21)='grid_longitude_at_cell_y_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3202, MyFile)) RETURN
          END IF
          Vinfo( 1)='lat_u'
          Vinfo( 2)='latitude of U-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='latitude'
          Vinfo(21)='grid_latitude_at_cell_y_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3213, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_u'
          Vinfo( 2)='x-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Xu'
          Vinfo(21)='grid_x_location_at_cell_y_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3224, MyFile)) RETURN
          END IF
          Vinfo( 1)='y_u'
          Vinfo( 2)='y-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Yu'
          Vinfo(21)='grid_y_location_at_cell_y_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3235, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_v'
          Vinfo( 2)='longitude of V-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='longitude'
          Vinfo(21)='grid_longitude_at_cell_x_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3250, MyFile)) RETURN
          END IF
          Vinfo( 1)='lat_v'
          Vinfo( 2)='latitude of V-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='latitude'
          Vinfo(21)='grid_latitude_at_cell_x_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3261, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_v'
          Vinfo( 2)='x-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Xv'
          Vinfo(21)='grid_x_location_at_cell_x_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3272, MyFile)) RETURN
          END IF
          Vinfo( 1)='y_v'
          Vinfo( 2)='y-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Yv'
          Vinfo(21)='grid_y_location_at_cell_x_edges'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3283, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_psi'
          Vinfo( 2)='longitude of PSI-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='longitude'
          Vinfo(21)='grid_longitude_at_cell_corners'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3298, MyFile)) RETURN
          END IF
          Vinfo( 1)='lat_psi'
          Vinfo( 2)='latitude of PSI-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='latitude'
          Vinfo(21)='grid_latitude_at_cell_corners'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3309, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_psi'
          Vinfo( 2)='x-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Xp'
          Vinfo(21)='grid_x_location_at_cell_corners'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3320, MyFile)) RETURN
          END IF
          Vinfo( 1)='y_psi'
          Vinfo( 2)='y-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='Yp'
          Vinfo(21)='grid_y_location_at_cell_corners'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3331, MyFile)) RETURN
          END IF
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='mask_rho'
          Vinfo( 2)='mask on RHO-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)='land_sea_mask_at_cell_center'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3370, MyFile)) RETURN
          Vinfo( 1)='mask_u'
          Vinfo( 2)='mask on U-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)='land_sea_mask_at_cell_y_edges'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3381, MyFile)) RETURN
          Vinfo( 1)='mask_v'
          Vinfo( 2)='mask on V-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)='land_sea_mask_at_cell_x_edges'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3392, MyFile)) RETURN
          Vinfo( 1)='mask_psi'
          Vinfo( 2)='mask on psi-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)='land_sea_mask_at_cell_corners'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3403, MyFile)) RETURN
        END IF
      END IF
!
  10  FORMAT (i3.3,'x',i3.3)
  20  FORMAT (/,' DEF_INFO_NF90 - error while creating global'          &
     &        ' attribute: ',a,/,17x,a)
  30  FORMAT (a,i2.2)
!
      RETURN
      END SUBROUTINE def_info_nf90
      END MODULE def_info_mod
