      SUBROUTINE ntimesteps (model, RunInterval, nl, Nsteps, Rsteps)
!
!svn $Id: ntimestep.F 1099 2022-01-06 21:01:01Z arango $
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine set the number of time-steps to compute. In nesting    !
!  applications,  the number of time-steps depends on nesting layer    !
!  configuration.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model        Calling model identifier (integer)                  !
!     RunInterval  Time interval (seconds) to advance forward or       !
!                    backwards the primitive equations (scalar)        !
!     nl           Nesting layer number (integer)                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     nl           Updated nesting layer number (integer)              !
!     Nsteps       Number of time-steps to solve for all grids in      !
!                    current nesting layer "nl" (integer)              !
!     Rsteps       Number of time-steps to complete RunInterval        !
!                    time window (integer)                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in)    :: model
      integer, intent(inout) :: nl
      integer, intent(out)   :: Nsteps, Rsteps
!
      real(dp), intent(in)   :: RunInterval
!
!  Local variable declarations.
!
      integer :: extra, gn, ig, il, ng, ngm1
      integer, dimension(Ngrids) :: WindowSteps, my_Nsteps
!
!=======================================================================
!  Set number of time steps to execute for grid "ng".
!=======================================================================
!
!  Initialize.
!
      WindowSteps=0
      my_Nsteps=0
      Nsteps=0
      Rsteps=0
!
!  If appropriate, advance the nested layer counter.
!
      nl=nl+1
!
     IF ((0.lt.nl).and.(nl.le.NestLayers)) THEN
!
!  In composite and mosaic grids or all grids in the same nesting
!  layer, it is assumed that the donor and receiver grids have the
!  same time-step size. This is done to avoid the time interpolation
!  between donor and receiver grids. Only spatial interpolations
!  are possible in the current nesting design.
!
!  In grid refinement, it is assumed that the donor and receiver grids
!  are an interger factor of the grid size and time-step size.
!
        WindowSteps=0
!
!  Loop over all grids in current layer nesting layer.
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
!
!  Here, extra=1, indicates that the RunInterval is the same as
!  simulation interval.
!
          extra=1
!
!  Determine number of steps in time-interval window.
!
          WindowSteps(ig)=INT((RunInterval+0.5_r8*dt(ng))/dt(ng))
!
!  Advancing model forward: Nonlinear, tangent linear, and representer
!  models.
!
          IF ((model.eq.iNLM).or.                                       &
     &        (model.eq.iTLM).or.                                       &
     &        (model.eq.iRPM)) THEN
            my_Nsteps(ig)=MAX(my_Nsteps(ig), WindowSteps(ig)+extra)
            step_counter(ng)=step_counter(ng)+WindowSteps(ig)+extra
          END IF
        END DO
!
!  Insure that the steps per time-window are the same.
!
        IF (GridsInLayer(nl).gt.1) THEN
          DO ig=2,GridsInLayer(nl)
            IF (WindowSteps(ig).ne.WindowSteps(ig-1)) THEN
              ngm1=GridNumber(ig-1,nl)
              ng  =GridNumber(ig  ,nl)
              IF (Master) THEN
                WRITE (stdout,10) nl, ngm1, dt(ngm1), ng, dt(ng)
  10            FORMAT (/,' NTIMESTEPS - timestep size are not the ',   &
     &                    ' same in nesting layer: ',i2,                &
     &                  2(/,14x,'Grid ',i2.2,3x,'dt = ',f11.3))
              END IF
              exit_flag=5
            END IF
          END DO
        END IF
!
!  Set number of time-steps to execute. Choose minimum values.
!
        Nsteps=my_Nsteps(1)
        Rsteps=WindowSteps(1)+extra
        DO ig=2,GridsInLayer(nl)
          Nsteps=MIN(Nsteps, my_Nsteps(ig))
          Rsteps=MIN(Rsteps, WindowSteps(ig)+extra)
        END DO
      END IF
      RETURN
      END SUBROUTINE ntimesteps
