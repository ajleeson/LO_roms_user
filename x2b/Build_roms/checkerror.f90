      FUNCTION CheckError (ErrorFlag, source, routine, message)
!
!svn $Id: checkerror.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This function checks an error flag returned by a particular source   !
! code and prints the appropriate error message, if any.               !
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
      integer, intent(in) :: ErrorFlag
      character (len=*), intent(in) :: source
      character (len=*), intent(in) :: routine
      character (len=*), intent(in) :: message
!
!  Local variable declarations.
!
      logical :: CheckError
!
!-----------------------------------------------------------------------
!  Check error flag.
!-----------------------------------------------------------------------
!
!  First, initialize to no error.
!
      CheckError=.FALSE.
!
!  Check error flags according to source.
!
      SELECT CASE (TRIM(ADJUSTL(source)))
        CASE ('ROMS')
          IF (ErrorFLAG.ne.NoError) THEN
            CheckError=.TRUE.
          END IF
        CASE DEFAULT
          CheckError=.FALSE.
      END SELECT
!
!  Report error message.
!
      IF (CheckError) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(ADJUSTL(routine)),                     &
     &                      TRIM(ADJUSTL(source)),                      &
     &                      TRIM(ADJUSTL(message))
        END IF
 10     FORMAT(/,1x,a,' - ',a,': Error while ',a,'.'/)
      END IF
      RETURN
      END FUNCTION CheckError
