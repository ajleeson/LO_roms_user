      MODULE stiffness_mod
!
!svn $Id: stiffness.F 1099 2022-01-06 21:01:01Z arango $
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine surveys the 3D grid in order to determine maximum      !
!  grid stiffness ratio:                                               !
!                                                                      !
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)              !
!      r_x = ---------------------------------------------             !
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)              !
!                                                                      !
!  This is done for diagnostic purposes and it does not affect the     !
!  computations.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: stiffness
      CONTAINS
!
!***********************************************************************
      SUBROUTINE stiffness (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL stiffness_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % omn,                              &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng)% zeta)
      RETURN
      END SUBROUTINE stiffness
!
!***********************************************************************
      SUBROUTINE stiffness_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           rmask, umask, vmask,                   &
     &                           h, omn,                                &
     &                           Hz, z_w,                               &
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_reduce
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      real(r8) :: cff, ratio
      real(r8) :: my_rx0, my_rx1
      real(r8) :: my_volume0, my_volume1, my_volume2
      real(r8), dimension(5) :: rbuffer
      character (len=3), dimension(5) :: op_handle
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      my_rx0=0.0_r8
      my_rx1=0.0_r8
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (umask(i,j).gt.0.0_r8) THEN
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i-1,j,0))/            &
     &                            (z_w(i,j,0)+z_w(i-1,j,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i-1,j,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i-1,j,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))))
            END DO
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (vmask(i,j).gt.0.0_r8) THEN
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i,j-1,0))/            &
     &                            (z_w(i,j,0)+z_w(i,j-1,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i,j-1,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i,j-1,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))))
            END DO
          END IF
        END DO
      END DO
!
!-------------------------------------------------------------------------
!  Compute initial basin volume and grid cell minimum and maximum volumes.
!-------------------------------------------------------------------------
!
      my_volume0=0.0_r8
      my_volume1=1.0E+20_r8
      my_volume2=0.0_r8
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (rmask(i,j).gt.0.0_r8) THEN
              cff=omn(i,j)*Hz(i,j,k)
              my_volume0=my_volume0+cff
              my_volume1=MIN(my_volume1,cff)
              my_volume2=MAX(my_volume2,cff)
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute global values.
!-----------------------------------------------------------------------
!
      NSUB=1                             ! distributed-memory
!$OMP CRITICAL (R_FACTOR)
      TotVolume(ng)=TotVolume(ng)+my_volume0
      MinVolume(ng)=MIN(MinVolume(ng),my_volume1)
      MaxVolume(ng)=MAX(MaxVolume(ng),my_volume2)
      rx0(ng)=MAX(rx0(ng),my_rx0)
      rx1(ng)=MAX(rx1(ng),my_rx1)
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        rbuffer(1)=TotVolume(ng)
        rbuffer(2)=MinVolume(ng)
        rbuffer(3)=MaxVolume(ng)
        rbuffer(4)=rx0(ng)
        rbuffer(5)=rx1(ng)
        op_handle(1)='SUM'
        op_handle(2)='MIN'
        op_handle(3)='MAX'
        op_handle(4)='MAX'
        op_handle(5)='MAX'
        CALL mp_reduce (ng, model, 5, rbuffer, op_handle)
        TotVolume(ng)=rbuffer(1)
        MinVolume(ng)=rbuffer(2)
        MaxVolume(ng)=rbuffer(3)
        rx0(ng)=rbuffer(4)
        rx1(ng)=rbuffer(5)
        IF (Master.and.LwrtInfo(ng)) THEN
          WRITE (stdout,10) ng
  10      FORMAT (/,' Basin information for Grid ',i2.2,':',/)
          WRITE (stdout,20) rx0(ng), rx1(ng)
  20      FORMAT (' Maximum grid stiffness ratios:  rx0 = ',1pe14.6,    &
     &            ' (Beckmann and Haidvogel)',/,t34,'rx1 = ',1pe14.6,   &
     &            ' (Haney)')
          IF (MinVolume(ng).ne.0.0_r8) THEN
            ratio=MaxVolume(ng)/MinVolume(ng)
          ELSE
            ratio=0.0_r8
          END IF
          WRITE (stdout,30) TotVolume(ng), MinVolume(ng),               &
     &                      MaxVolume(ng), ratio
  30      FORMAT (/,' Initial domain volumes:  TotVolume = ',1p,e17.10, &
     &            0p,' m3',/,t26,'MinCellVol = ',1p,e17.10,0p,' m3',    &
     &            /,t26, 'MaxCellVol = ',1p,e17.10,0p,' m3',            &
     &            /,t26, '   Max/Min = ',1p,e17.10,0p,/)
        END IF
      END IF
!$OMP END CRITICAL (R_FACTOR)
      RETURN
      END SUBROUTINE stiffness_tile
      END MODULE stiffness_mod
