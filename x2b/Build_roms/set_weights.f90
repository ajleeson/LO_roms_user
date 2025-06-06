      SUBROUTINE set_weights (ng)
!
!svn $Id: set_weights.F 1099 2022-01-06 21:01:01Z arango $
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine sets the weigth functions for the time averaging of    !
!  2D fields over all short time-steps.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, j, iter
      real(dp) :: gamma, scale
      real(dp) :: cff1, cff2
      real(r16) :: wsum, shift, cff
!
!=======================================================================
!  Compute time-averaging filter for barotropic fields.
!=======================================================================
!
!  Initialize both sets of weights to zero.
!
      nfast(ng)=0
      DO i=1,2*ndtfast(ng)
        weight(1,i,ng)=0.0_dp
        weight(2,i,ng)=0.0_dp
      END DO
!
!-----------------------------------------------------------------------
!  Power-law shape filters.
!-----------------------------------------------------------------------
!
!  The power-law shape filters are given by:
!
!     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
!
!  where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
!  normalization are chosen to yield the correct zeroth-order
!  (normalization), first-order (consistency), and second-order moments,
!  resulting in overall second-order temporal accuracy for time-averaged
!  barotropic motions resolved by baroclinic time step. There parameters
!  are set in "mod_scalars".
!
      scale=(Falpha+1.0_dp)*(Falpha+Fbeta+1.0_dp)/                      &
     &      ((Falpha+2.0_dp)*(Falpha+Fbeta+2.0_dp)*REAL(ndtfast(ng),dp))
!
!  Find center of gravity of the primary weighting shape function and
!  iteratively adjust "scale" to place the  centroid exactly at
!  "ndtfast".
!
      gamma=Fgamma*MAX(0.0_dp, 1.0_dp-10.0_dp/REAL(ndtfast(ng),dp))
      DO iter=1,16
        nfast(ng)=0
        DO i=1,2*ndtfast(ng)
          cff=scale*REAL(i,dp)
          weight(1,i,ng)=cff**Falpha-cff**(Falpha+Fbeta)-gamma*cff
          IF (weight(1,i,ng).gt.0.0_dp) nfast(ng)=i
          IF ((nfast(ng).gt.0).and.(weight(1,i,ng).lt.0.0_dp)) THEN
            weight(1,i,ng)=0.0_dp
          END IF
        END DO
        wsum=0.0_r16
        shift=0.0_r16
        DO i=1,nfast(ng)
          wsum=wsum+weight(1,i,ng)
          shift=shift+weight(1,i,ng)*REAL(i,dp)
        END DO
        scale=scale*shift/(wsum*REAL(ndtfast(ng),dp))
      END DO
!
!-----------------------------------------------------------------------
!  Post-processing of primary weights.
!-----------------------------------------------------------------------
!
!  Although it is assumed that the initial settings of the primary
!  weights has its center of gravity "reasonably close" to NDTFAST,
!  it may be not so according to the discrete rules of integration.
!  The following procedure is designed to put the center of gravity
!  exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
!  applying basically an upstream advection of weights to eliminate
!  the mismatch iteratively. Once this procedure is complete primary
!  weights are normalized.
!
!  Find center of gravity of the primary weights and subsequently
!  calculate the mismatch to be compensated.
!
      DO iter=1,ndtfast(ng)
        wsum=0.0_r16
        shift=0.0_r16
        DO i=1,nfast(ng)
          wsum=wsum+weight(1,i,ng)
          shift=shift+REAL(i,dp)*weight(1,i,ng)
        END DO
        shift=shift/wsum
        cff=REAL(ndtfast(ng),dp)-shift
!
!  Apply advection step using either whole, or fractional shifts.
!  Notice that none of the four loops here is reversible.
!
        IF (cff.gt.1.0_r16) THEN
          nfast(ng)=nfast(ng)+1
          DO i=nfast(ng),2,-1
            weight(1,i,ng)=weight(1,i-1,ng)
          END DO
          weight(1,1,ng)=0.0_dp
        ELSE IF (cff.gt.0.0_r16) THEN
          wsum=1.0_r16-cff
          DO i=nfast(ng),2,-1
            weight(1,i,ng)=wsum*weight(1,i,ng)+cff*weight(1,i-1,ng)
          END DO
          weight(1,1,ng)=wsum*weight(1,1,ng)
        ELSE IF (cff.lt.-1.0_r16) THEN
          nfast(ng)=nfast(ng)-1
          DO i=1,nfast(ng),+1
            weight(1,i,ng)=weight(1,i+1,ng)
          END DO
          weight(1,nfast(ng)+1,ng)=0.0_dp
        ELSE IF (cff.lt.0.0_r16) THEN
          wsum=1.0_r16+cff
          DO i=1,nfast(ng)-1,+1
            weight(1,i,ng)=wsum*weight(1,i,ng)-cff*weight(1,i+1,ng)
          END DO
          weight(1,nfast(ng),ng)=wsum*weight(1,nfast(ng),ng)
        END IF
      END DO
!
!  Set SECONDARY weights assuming that backward Euler time step is used
!  for free surface.  Notice that array weight(2,i,ng) is assumed to
!  have all-zero status at entry in this segment of code.
!
      DO j=1,nfast(ng)
        cff=weight(1,j,ng)
        DO i=1,j
          weight(2,i,ng)=weight(2,i,ng)+cff
        END DO
      END DO
!
!  Normalize both set of weights.
!
      wsum=0.0_r16
      cff=0.0_r16
      DO i=1,nfast(ng)
        wsum=wsum+weight(1,i,ng)
        cff=cff+weight(2,i,ng)
      END DO
      wsum=1.0_r16/wsum
      cff=1.0_r16/cff
      DO i=1,nfast(ng)
        weight(1,i,ng)=wsum*weight(1,i,ng)
        weight(2,i,ng)=cff*weight(2,i,ng)
      END DO
!
!  Report weights.
!
      IF (Master.and.LwrtInfo(ng)) THEN
        WRITE (stdout,10) ng, ndtfast(ng), nfast(ng)
        cff=0.0_r16
        cff1=0.0_dp
        cff2=0.0_dp
        wsum=0.0_r16
        shift=0.0_r16
        DO i=1,nfast(ng)
          cff=cff+weight(1,i,ng)
          cff1=cff1+weight(1,i,ng)*REAL(i,dp)
          cff2=cff2+weight(1,i,ng)*REAL(i*i,dp)
          wsum=wsum+weight(2,i,ng)
          shift=shift+weight(2,i,ng)*(REAL(i,dp)-0.5_dp)
          WRITE (stdout,20) i, weight(1,i,ng), weight(2,i,ng), cff, wsum
        END DO
        cff1=cff1/REAL(ndtfast(ng),dp)
        cff2=cff2/(REAL(ndtfast(ng),dp)*REAL(ndtfast(ng),dp))
        shift=shift/REAL(ndtfast(ng),dp)
        WRITE (stdout,30) ndtfast(ng), nfast(ng),                       &
     &                    REAL(nfast(ng),dp)/REAL(ndtfast(ng),dp)
        WRITE (stdout,40) cff1, cff2, shift, cff, wsum, Fgamma, gamma
        IF (cff2.lt.1.0001_dp) WRITE (stdout,50)
      END IF
!
  10  FORMAT (/,' Time Splitting Weights for Grid ',i2.2,               &
     &        ':',4x,'ndtfast = ',i3,4x,'nfast = ',i3,/,                &
     &          ' ==================================',/,                &
     &        /,4x,'Primary',12x,'Secondary',12x,                       &
     &        'Accumulated to Current Step',/)
  20  FORMAT (i3,4f19.16)
  30  FORMAT (/,1x,'ndtfast, nfast = ',2i4,3x,'nfast/ndtfast = ',f8.5)
  40  FORMAT (/,1x,'Centers of gravity and integrals ',                 &
     &        '(values must be 1, 1, approx 1/2, 1, 1):',/,             &
     &        /,3x,5F15.12,/,/,1x,'Power filter parameters, ',          &
     &        'Fgamma, gamma = ', f8.5,2x,f8.5)
  50  FORMAT (/,' WARNING: unstable weights, reduce parameter',         &
     &          ' Fgamma in mod_scalars.F',/)
      RETURN
      END SUBROUTINE set_weights
