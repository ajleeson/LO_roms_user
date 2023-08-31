      MODULE wrt_info_mod
!
!svn $Id: wrt_info.F 1113 2022-03-01 20:31:29Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_grid
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_sources
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      INTERFACE wrt_info
        MODULE PROCEDURE wrt_info_nf90
      END INTERFACE wrt_info
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_info_nf90 (ng, model, ncid, ncname)
!***********************************************************************
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file using the standard NetCDF-3 or NetCDF-4 library.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF filename (string)                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE distribute_mod,  ONLY : mp_bcasti
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid
      integer, dimension(2) :: ibuffer
      integer :: ifield = 0
!
      real(dp) :: scale
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_info.F"//", wrt_info_nf90"
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 132, MyFile)) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 139, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 144, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 149, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 154, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 159, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'shuffle',               &
     &                      shuffle, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 177, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'deflate',               &
     &                      deflate, (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 182, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'deflate_level',         &
     &                      deflate_level, (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 187, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 193, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 198, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 203, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Falpha',                &
     &                      Falpha, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 344, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fbeta',                 &
     &                      Fbeta, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 349, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fgamma',                &
     &                      Fgamma, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 354, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
!
!  Background vertical mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Akt_bak',               &
     &                      Akt_bak(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 490, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akv_bak',               &
     &                      Akv_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 495, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akk_bak',               &
     &                      Akk_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 501, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akp_bak',               &
     &                      Akp_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 506, MyFile)) RETURN
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 548, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 553, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zob',                   &
     &                      Zob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 559, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zos',                   &
     &                      Zos(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 564, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_p',                 &
     &                      gls_p(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 574, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_m',                 &
     &                      gls_m(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 579, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_n',                 &
     &                      gls_n(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 584, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',              &
     &                      gls_cmu0(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 589, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c1',                &
     &                      gls_c1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 594, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c2',                &
     &                      gls_c2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 599, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3m',               &
     &                      gls_c3m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 604, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3p',               &
     &                      gls_c3p(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 609, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigk',              &
     &                      gls_sigk(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 614, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigp',              &
     &                      gls_sigp(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 619, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',              &
     &                      gls_Kmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 624, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',              &
     &                      gls_Pmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 629, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',         &
     &                      charnok_alpha(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 634, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',        &
     &                      zos_hsig_alpha(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 639, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'sz_alpha',              &
     &                      sz_alpha(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 644, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',             &
     &                      crgban_cw(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 649, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 657, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 662, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudg',                &
     &                      M3nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 668, MyFile)) RETURN
      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudg',                 &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 676, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',            &
     &                        FSobc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 687, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',           &
     &                        FSobc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 692, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',            &
     &                        M2obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 697, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',           &
     &                        M2obc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 702, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_in',             &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 713, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_out',            &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 723, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_in',            &
     &                        M3obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 728, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_out',           &
     &                        M3obc_out(ng,:), (/1/), (/4/),            &
     &                      ncid = ncid)
        IF (FoundError(exit_flag, NoError, 733, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 743, MyFile)) RETURN
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 791, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSrc',                &
     &                      LuvSrc(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 800, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LwSrc',                 &
     &                      LwSrc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 805, MyFile)) RETURN
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',            &
     &                      LtracerSrc(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 814, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LsshCLM',               &
     &                      LsshCLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 822, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',                &
     &                      Lm2CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 827, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',                &
     &                      Lm3CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 833, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',            &
     &                      LtracerCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 838, MyFile)) RETURN
!
!  Logical switches for nudging climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',           &
     &                      LnudgeM2CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 846, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',           &
     &                      LnudgeM3CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 852, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',            &
     &                      LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 857, MyFile)) RETURN
!
!  Write out Fennel et al. (2006) ecosystem model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 21, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttSW',                 &
     &                      AttSW(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 26, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttChl',                &
     &                      AttChl(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 31, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PARfrac',               &
     &                      PARfrac(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 36, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Vp0',                   &
     &                      Vp0(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 41, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'I_thNH4',               &
     &                      I_thNH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 46, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'D_p5NH4',               &
     &                      D_p5NH4(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 51, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'NitriR',                &
     &                      NitriR(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 56, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_NO3',                 &
     &                      K_NO3(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 61, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_NH4',                 &
     &                      K_NH4(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 66, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_PO4',                 &
     &                      K_PO4(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 71, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'K_Phy',                 &
     &                      K_Phy(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 76, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Chl2C_m',               &
     &                      Chl2C_m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 81, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ChlMin',                &
     &                      ChlMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 86, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyCN',                 &
     &                      PhyCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 91, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'R_P2N',                 &
     &                      R_P2N(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 96, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIP',                 &
     &                      PhyIP(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 101, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyIS',                 &
     &                      PhyIS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 106, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMin',                &
     &                      PhyMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 111, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PhyMR',                 &
     &                      PhyMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 116, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooAE_N',               &
     &                      ZooAE_N(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 121, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooBM',                 &
     &                      ZooBM(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 126, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooCN',                 &
     &                      ZooCN(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 131, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooER',                 &
     &                      ZooER(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 136, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooGR',                 &
     &                      ZooGR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 141, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMin',                &
     &                      ZooMin(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 146, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'ZooMR',                 &
     &                      ZooMR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 151, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRN',                &
     &                      LDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 156, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'LDeRRC',                &
     &                      LDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 161, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'CoagR',                 &
     &                      CoagR(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 166, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRN',                &
     &                      SDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 171, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'SDeRRC',                &
     &                      SDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 176, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'RDeRRN',                &
     &                      RDeRRN(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 181, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'RDeRRC',                &
     &                      RDeRRC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 186, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wPhy',                  &
     &                      wPhy(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 191, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wLDet',                 &
     &                      wLDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 196, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'wSDet',                 &
     &                      wSDet(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 201, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'pCO2air',               &
     &                      pCO2air(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 206, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'spherical',             &
     &                      spherical, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1324, MyFile)) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1331, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1336, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'Vtransform',            &
     &                      Vtransform(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1345, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'Vstretching',           &
     &                      Vstretching(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1350, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_s',               &
     &                      theta_s(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1355, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_b',               &
     &                      theta_b(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1360, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Tcline',                &
     &                      Tcline(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1365, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'hc',                    &
     &                      hc(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1370, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'grid',                  &
     &                      (/1/), (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1379, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL netcdf_put_fvar (ng, model, ncname, 's_rho',                 &
     &                      SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1386, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 's_w',                   &
     &                      SCALARS(ng)%sc_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1391, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_r',                  &
     &                      SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1398, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_w',                  &
     &                      SCALARS(ng)%Cs_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1403, MyFile)) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 1412, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Write out grid tiled variables.
!-----------------------------------------------------------------------
!
      GRID_VARS : IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF (ncid.ne.STA(ng)%ncid) THEN
            IF (find_string(var_name, n_var, 'h', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1456, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'f', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1495, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'pm', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1521, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'pn', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1545, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lon_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1573, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lon_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lat_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1611, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lat_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'x_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1651, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'x_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'y_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1689, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'y_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lon_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1731, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lon_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lat_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1756, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lat_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'x_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1783, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'x_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'y_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1808, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'y_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lon_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1837, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lat_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1862, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lat_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'x_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1889, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'y_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1914, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'y_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lon_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1943, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lon_psi', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_psi', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'lat_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1968, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'lat_psi', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'x_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1995, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'x_psi', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_psi', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, 'y_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2020, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) 'y_psi', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_psi', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'mask_rho', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2088, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'mask_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'mask_u', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2110, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'mask_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'mask_v', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2132, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'mask_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, 'mask_psi', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2154, MyFile)) THEN
                IF (Master) WRITE (stdout,10) 'mask_psi', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 2333, MyFile)) RETURN
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
  10  FORMAT (/,' WRT_INFO_NF90 - error while writing variable: ',a,/,  &
     &        17x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO_NF90 - error while inquiring ID for',        &
     &        ' variable: ',a,/,17x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO_NF90 - unable to synchronize to disk file:', &
     &        /,17x,a)
!
      RETURN
      END SUBROUTINE wrt_info_nf90
      END MODULE wrt_info_mod
