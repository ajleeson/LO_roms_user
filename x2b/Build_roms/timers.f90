      RECURSIVE SUBROUTINE wclock_on (ng, model, region, line, routine)
!
!svn $Id: timers.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine turns on wall clock to meassure the elapsed time in    !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     region     Profiling reagion number (integer)                    !
!     line       Calling model routine line (integer)                  !
!     routine    Calling model routine (string)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      USE distribute_mod, ONLY : mp_barrier
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, region, line
      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      integer :: iregion, MyModel, NSUB
      integer :: my_getpid
      integer :: MyCOMM, nPETs, PETrank
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
!
!-----------------------------------------------------------------------
!  Initialize timing for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
      MyCOMM=OCN_COMM_WORLD
      nPETs=numthreads
      PETrank=MyRank
      NSUB=1
!
!  Insure that MyModel is not zero.
!
      MyModel=MAX(1,model)
!
!  Start the wall CPU clock for specified region, model, and grid.
!
      Cstr(region,MyModel,ng)=my_wtime(wtime)
!
!  If region zero, indicating first call from main driver, initialize
!  time profiling arrays and set process ID.
!
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,MyModel,ng)=0.0_r8
          Csum(iregion,MyModel,ng)=0.0_r8
        END DO
        proc(1,MyModel,ng)=1
        proc(0,MyModel,ng)=my_getpid()
!$OMP CRITICAL (START_WCLOCK)
        IF (ng.eq.1) THEN
          CALL mp_barrier (ng, model, MyCOMM)
          WRITE (stdout,10) ' Node #', PETrank,                         &
     &                      ' (pid=',proc(0,MyModel,ng),') is active.'
          CALL my_flush (stdout)
        END IF
 10     FORMAT (a,i5,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.NSUB) thread_count=0
!$OMP END CRITICAL (START_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_on
!
      RECURSIVE SUBROUTINE wclock_off (ng, model, region, line, routine)
!
!=======================================================================
!                                                                      !
!  This routine turns off wall clock to meassure the elapsed time in   !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     region     Profiling region number (integer)                     !
!     line       Calling model routine line (integer)                  !
!     routine    Calling model routine (string)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      USE distribute_mod, ONLY : mp_barrier, mp_reduce
      USE distribute_mod, ONLY : mp_collect
      USE strings_mod,    ONLY : uppercase
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region, line
      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      integer :: ig, imodel, iregion, MyModel, NSUB
      integer :: MyCOMM, nPETs, PETrank
      integer :: my_threadnum
      real(r8) :: percent, sumcpu, sumper, sumsum, total
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
      real(r8) :: TendMin, TendMax
      real(r8), parameter :: Tspv = 0.0_r8
      real(r8), allocatable :: Tend(:)
      real(r8), dimension(0:Nregion) :: rbuffer
      character (len= 3), dimension(0:Nregion) :: op_handle
      character (len=14), dimension(4) :: label
!
!-----------------------------------------------------------------------
!  Compute elapsed wall time for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
      MyCOMM=OCN_COMM_WORLD
      nPETs=numthreads
      PETrank=MyRank
      NSUB=1
!
!  Insure that MyModel is not zero.
!
      MyModel=MAX(1,model)
!
!  Compute elapsed CPU time (seconds) for each profile region, except
!  for region zero which is called by the main driver before the
!  simulatiom is stopped.
!
      IF (region.ne.0) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
      END IF
!
!-----------------------------------------------------------------------
!  If simulation is compleated, compute and report elapsed CPU time for
!  all regions.
!-----------------------------------------------------------------------
!
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.1)) THEN
!
!  Computed elapsed wall time for the driver, region=0.  Since it is
!  called only once, "MyModel" will have a value and the other models
!  will be zero.
!
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
        DO imodel=1,4
          proc(1,imodel,ng)=0
        END DO
!$OMP CRITICAL (FINALIZE_WCLOCK)
!
!  Report elapsed time (seconds) for each CPU.  We get the same time
!  time for all nested grids.
!
        IF (ng.eq.1) THEN
          CALL mp_barrier (ng, model, MyCOMM)
          IF (.not.allocated(Tend)) THEN
            allocate ( Tend(nPETs) )
            Tend=0.0_r8
          END IF
          Tend(PETrank+1)=Cend(region,MyModel,ng)
          CALL mp_collect (ng, model, nPETs, Tspv, Tend, MyCOMM)
          TendMin=MINVAL(Tend)
          TendMax=MAXVAL(Tend)
          WRITE (stdout,10) ' Node   #', PETrank,                       &
     &                      ' CPU:', Tend(PETrank+1)
          CALL my_flush (stdout)
 10       FORMAT (a,i5,a,f12.3)
        END IF
!
!  Sum the elapsed time for each profile region by model.
!
        thread_count=thread_count+1
        DO imodel=1,4
          DO iregion=0,Nregion
            Csum(iregion,imodel,ng)=Csum(iregion,imodel,ng)+            &
     &                              Cend(iregion,imodel,ng)
          END DO
        END DO
!
!  Compute total elapsed CPU wall time between all parallel processes.
!
        IF (thread_count.eq.NSUB) THEN
          thread_count=0
          op_handle(0:Nregion)='SUM'      ! Gather all values using a
          DO imodel=1,4                   ! reduced sum between nodes
            DO iregion=0,Nregion
              rbuffer(iregion)=Csum(iregion,imodel,ng)
            END DO
            CALL mp_reduce (ng, MyModel, Nregion+1, rbuffer(0:),        &
     &                      op_handle(0:), MyCOMM)
            DO iregion=0,Nregion
              Csum(iregion,imodel,ng)=rbuffer(iregion)
            END DO
          END DO
          IF (Master) THEN
            IF (ng.eq.1) THEN             ! Same for all nested grids
              total_cpu=total_cpu+Csum(region,model,ng)
            END IF
            DO imodel=1,4
              total_model(imodel)=0.0_r8
              DO iregion=1,Nregion
                total_model(imodel)=total_model(imodel)+                &
     &                              Csum(iregion,imodel,ng)
              END DO
            END DO
            IF (ng.eq.1) THEN
              WRITE (stdout,20) ' Total:', total_cpu
 20           FORMAT (a,t18,f14.3)
              IF (numthreads.gt.1) THEN
                WRITE (stdout,20) ' Average:', total_cpu/numthreads
                WRITE (stdout,20) ' Minimum:', TendMin
                WRITE (stdout,20) ' Maximum:', TendMax
              END IF
            END IF
          END IF
          IF (allocated(Tend)) deallocate (Tend)
!
!  Report profiling times.
!
          label(iNLM)='Nonlinear     '
          label(iTLM)='Tangent linear'
          label(iRPM)='Representer   '
          label(iADM)='Adjoint       '
          DO imodel=1,4
            IF (Master.and.(total_model(imodel).gt.0.0_r8)) THEN
              WRITE (stdout,30) TRIM(label(imodel)),                    &
     &                          'model elapsed CPU time profile, Grid:',&
     &                          ng
 30           FORMAT (/,1x,a,1x,a,1x,i2.2/)
            END IF
            sumcpu=0.0_r8
            sumper=0.0_r8
            DO iregion=1,Mregion-1
              IF (Master.and.(Csum(iregion,imodel,ng).gt.0.0_r8)) THEN
                percent=100.0_r8*Csum(iregion, imodel,ng)/total_cpu
                WRITE (stdout,40) Pregion(iregion),                     &
     &                            Csum(iregion,imodel,ng), percent
                sumcpu=sumcpu+Csum(iregion,imodel,ng)
                sumper=sumper+percent
              END IF
            END DO
            Ctotal=Ctotal+sumcpu
 40         FORMAT (2x,a,t53,f14.3,2x,'(',f7.4,' %)')
            IF (Master.and.(total_model(imodel).gt.0.0_r8)) THEN
              WRITE (stdout,50) sumcpu, sumper
 50           FORMAT (t47,'Total:',f14.3,2x,f8.4,' %')
            END IF
          END DO
!
!  Sometimes the profiling does not fully accounts for all the CPU
!  spend outside of the ROMS kernels. In data assimilation algorithms,
!  there is a lot of CPU expend outside the kernels. A separated
!  profiling is reported bellow.
!
          IF (Master.and.(ng.eq.Ngrids)) THEN
            percent=100.0_r8*Ctotal/total_cpu
            WRITE (stdout,60) Ctotal, percent,                          &
     &                        total_cpu-Ctotal, 100.0_r8-percent
 60         FORMAT (/,2x,                                               &
     &             'Unique kernel(s) regions profiled ................',&
     &               f14.3,2x,f8.4,' %'/,2x,                            &
     &             'Residual, non-profiled code ......................',&
     &               f14.3,2x,f8.4,' %'/)
            WRITE (stdout,70) total_cpu
 70         FORMAT (/,' All percentages are with respect to',           &
     &                ' total time =',5x,f12.3,/)
          END IF
!
!  Report elapsed time for message passage communications.
!
          total=0.0_r8
          DO iregion=Mregion,Fregion-1
            DO imodel=1,4
              total=total+Csum(iregion,imodel,ng)
            END DO
          END DO
          IF (Master.and.(total.gt.0.0_r8)) THEN
            WRITE (stdout,30) uppercase('mpi'),                         &
     &                        'communications profile, Grid:', ng
          END IF
          IF (total.gt.0.0_r8) THEN
            sumper=0.0_r8
            sumsum=0.0_r8
            DO iregion=Mregion,Fregion-1
              sumcpu=0.0_r8
              DO imodel=1,4
                sumcpu=sumcpu+Csum(iregion,imodel,ng)
              END DO
              IF (Master.and.(sumcpu.gt.0.0_r8)) THEN
                percent=100.0_r8*sumcpu/total_cpu
                WRITE (stdout,40) Pregion(iregion), sumcpu, percent
                sumsum=sumsum+sumcpu
                sumper=sumper+percent
              END IF
            END DO
            IF (Master.and.(total.gt.0.0_r8)) THEN
              WRITE (stdout,50) sumsum, sumper
            END IF
          END IF
        END IF
!$OMP END CRITICAL (FINALIZE_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_off
