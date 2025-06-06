      MODULE set_3dfld_mod
!
!svn $Id: set_3dfld.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates requested 3D field from snapshots    !
!  of input data.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_3dfld_tile (ng, tile, model, ifield,               &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           Finp, Fout, update,                    &
     &                           SetBC)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE dateclock_mod,   ONLY : time_string
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetBC
      logical, intent(out) :: update
      integer, intent(in) :: ng, tile, model, ifield
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(in) :: Finp(LBi:,LBj:,LBk:,:)
      real(r8), intent(out) :: Fout(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: LapplyBC, Lgrided, Lonerec
!
      integer :: Tindex, gtype, i, it1, it2, j, k
!
      real(dp) :: SecScale, fac, fac1, fac2
      real(r8) :: Fval
!
      character (len=22) :: DateString1, DateString2
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
!--------------------------------------------------------------------
!  Set-up requested field for current tile.
!--------------------------------------------------------------------
!
!  Set switch to apply boundary conditions.
!
      IF (PRESENT(SetBC)) THEN
        LapplyBC=SetBC
      ELSE
        LapplyBC=.TRUE.
      END IF
!
!  Get requested field information from global storage.
!
      Lgrided=Linfo(1,ifield,ng)
      Lonerec=Linfo(3,ifield,ng)
      gtype  =Iinfo(1,ifield,ng)
      Tindex =Iinfo(8,ifield,ng)
      update=.TRUE.
!
!  Set linear, time interpolation factors. Fractional seconds are
!  rounded to the nearest milliseconds integer towards zero in the
!  time interpolation weights.
!
      SecScale=1000.0_dp              ! seconds to milliseconds
      it1=3-Tindex
      it2=Tindex
      fac1=ANINT((Tintrp(it2,ifield,ng)-time(ng))*SecScale,dp)
      fac2=ANINT((time(ng)-Tintrp(it1,ifield,ng))*SecScale,dp)
!
!  Load time-invariant data. Time interpolation is not necessary.
!
      IF (Lonerec) THEN
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Finp(i,j,k,Tindex)
              END DO
            END DO
          END DO
        ELSE
          Fval=Fpoint(Tindex,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Time-interpolate from gridded or point data.
!
      ELSE IF (((fac1*fac2).ge.0.0_dp).and.                             &
     &        ((fac1+fac2).gt.0.0_dp)) THEN
        fac=1.0_dp/(fac1+fac2)
        fac1=fac*fac1                             ! nondimensional
        fac2=fac*fac2                             ! nondimensional
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=fac1*Finp(i,j,k,it1)+fac2*Finp(i,j,k,it2)
              END DO
            END DO
          END DO
        ELSE
          Fval=fac1*Fpoint(it1,ifield,ng)+fac2*Fpoint(it2,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)+dt(ng)).gt.Tintrp(it2,ifield,ng)) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) synchro_flag(ng)=.TRUE.
        END IF
!
!  Unable to set-up requested field.  Activate error flag to quit.
!
      ELSE
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tdays(ng),         &
     &                        Finfo(1,ifield,ng), Finfo(2,ifield,ng),   &
     &                        Finfo(3,ifield,ng), Finfo(4,ifield,ng),   &
     &                        Tintrp(it1,ifield,ng)*sec2day,            &
     &                        Tintrp(it2,ifield,ng)*sec2day,            &
     &                        fac1*sec2day/SecScale,                    &
     &                        fac2*sec2day/SecScale
          END IF
  10      FORMAT (/,' SET_3DFLD  - current model time',                 &
     &            ' exceeds ending value for variable: ',a,             &
     &            /,14x,'TDAYS     = ',f15.4,                           &
     &            /,14x,'Data Tmin = ',f15.4,2x,'Data Tmax = ',f15.4,   &
     &            /,14x,'Data Tstr = ',f15.4,2x,'Data Tend = ',f15.4,   &
     &            /,14x,'TINTRP1   = ',f15.4,2x,'TINTRP2   = ',f15.4,   &
     &            /,14x,'FAC1      = ',f15.4,2x,'FAC2      = ',f15.4)
          exit_flag=2
          update=.FALSE.
        END IF
      END IF
!
!  Exchange boundary data.
!
      IF (update) THEN
        IF (LapplyBC.and.(EWperiodic(ng).or.NSperiodic(ng))) THEN
          IF (gtype.eq.r3dvar) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.u3dvar) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.v3dvar) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          ELSE IF (gtype.eq.w3dvar) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fout)
          END IF
        END IF
        IF (.not.LapplyBC) THEN
          CALL mp_exchange3d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        NghostPoints,                             &
     &                        .FALSE., .FALSE.,                         &
     &                        Fout)
        ELSE
          CALL mp_exchange3d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        Fout)
        END IF
      END IF
      RETURN
      END SUBROUTINE set_3dfld_tile
      END MODULE set_3dfld_mod
