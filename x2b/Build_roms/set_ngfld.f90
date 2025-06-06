      SUBROUTINE set_ngfld (ng, model, ifield,                          &
     &                      LBi, UBi, UBj, Istr, Iend, Jrec,            &
     &                      Finp, Fout, update)
!
!svn $Id: set_ngfld.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates requested non-grided field from      !
!  time snapshots of input data.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     LBi        Finp/Fout 1st dimension lower-bound value.            !
!     UBi        Finp/Fout 1st dimension upper-bound value.            !
!     UBj        Finp/Fout 2nd dimension upper-bound value, if any.    !
!                  Otherwise, a value of one is expected.              !
!     Istr       Starting location to process in the 1st dimension.    !
!     Iend       Ending location to process in the 1st dimension.      !
!     Jrec       Number of records to process in the 2nd dimenision,   !
!                  if any, Otherwise, a value of one is expected.      !
!     Finp       Latest two-snapshopts of field to interpolate.        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated field.                                   !
!     update     Switch indicating successful interpolation.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update
      integer, intent(in) :: ng, model, ifield
      integer, intent(in) :: LBi, UBi, UBj, Istr, Iend, Jrec
      real(r8), intent(in) :: Finp(LBi:UBi,UBj,2)
      real(r8), intent(out) :: Fout(LBi:UBi,UBj)
!
!  Local variable declarations.
!
      logical :: Lonerec
      integer :: Tindex, i, it1, it2, j
      real(dp) :: SecScale, fac, fac1, fac2
!
!----------------------------------------------------------------------
!  Set up requested field from data snapshots.
!----------------------------------------------------------------------
!
!  Get requested field information from global storage.
!
      Lonerec=Linfo(3,ifield,ng)
      Tindex=Iinfo(8,ifield,ng)
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
        DO j=1,Jrec
          DO i=Istr,Iend
            Fout(i,j)=Finp(i,j,Tindex)
          END DO
        END DO
!
!  Time-interpolate.
!
      ELSE IF (((fac1*fac2).ge.0.0_dp).and.(fac1+fac2).gt.0.0_dp) THEN
        fac=1.0_dp/(fac1+fac2)
        fac1=fac*fac1                             ! nondimensional
        fac2=fac*fac2                             ! nondimensional
        DO j=1,Jrec
          DO i=Istr,Iend
            Fout(i,j)=fac1*Finp(i,j,it1)+fac2*Finp(i,j,it2)
          END DO
        END DO
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)+dt(ng)).gt.Tintrp(it2,ifield,ng)) THEN
          synchro_flag(ng)=.TRUE.
        END IF
!
!  Unable to interpolate field.  Activate error flag to quit.
!
      ELSE
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tdays(ng),           &
     &                      Finfo(1,ifield,ng), Finfo(2,ifield,ng),     &
     &                      Finfo(3,ifield,ng), Finfo(4,ifield,ng),     &
     &                      Tintrp(it1,ifield,ng)*sec2day,              &
     &                      Tintrp(it2,ifield,ng)*sec2day,              &
     &                      fac1*sec2day/SecScale,                      &
     &                      fac2*sec2day/SecScale
        END IF
  10    FORMAT (/,' SET_NGFLD  - current model time',                   &
     &          ' exceeds ending value for variable: ',a,               &
     &          /,14x,'TDAYS     = ',f15.4,                             &
     &          /,14x,'Data Tmin = ',f15.4,2x,'Data Tmax = ',f15.4,     &
     &          /,14x,'Data Tstr = ',f15.4,2x,'Data Tend = ',f15.4,     &
     &          /,14x,'TINTRP1   = ',f15.4,2x,'TINTRP2   = ',f15.4,     &
     &          /,14x,'FAC1      = ',f15.4,2x,'FAC2      = ',f15.4)
        exit_flag=2
        update=.FALSE.
      END IF
      RETURN
      END SUBROUTINE set_ngfld
