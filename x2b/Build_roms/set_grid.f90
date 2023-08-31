      SUBROUTINE set_grid (ng, model)
!
!svn $Id: set_grid.F 1099 2022-01-06 21:01:01Z arango $
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine sets application grid and associated variables and     !
!  parameters. It called only once during the initialization stage.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE analytical_mod
      USE distribute_mod,       ONLY : mp_bcasti
      USE get_grid_mod,         ONLY : get_grid
      USE get_nudgcoef_mod,     ONLY : get_nudgcoef
      USE metrics_mod,          ONLY : metrics
      USE strings_mod,          ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/set_grid.F"
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL get_grid (ng, MyRank, model)
!$OMP END MASTER
      CALL mp_bcasti (ng, model, exit_flag)
!$OMP BARRIER
      IF (FoundError(exit_flag, NoError, 85, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Set vertical terrain-following coordinate transformation function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL set_scoord (ng)
!$OMP END MASTER
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL set_weights (ng)
!$OMP END MASTER
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL metrics (ng, tile, model)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
      IF (Lnudging(ng)) THEN
!$OMP MASTER
        CALL get_nudgcoef (ng, MyRank, model)
!$OMP END MASTER
        CALL mp_bcasti (ng, model, exit_flag)
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, 173, MyFile)) RETURN
      END IF
!
      RETURN
      END SUBROUTINE set_grid
