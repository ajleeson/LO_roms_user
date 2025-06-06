      MODULE close_io_mod
!
!svn $Id: close_io.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module closes input and output files using either the standard !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!  During initialization, the input input files need to be in closed   !
!  state to facilitate multi-file processing. This is important in     !
!  iterative algorithms that run the model kernels repetitevely.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE dateclock_mod, ONLY : get_date
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
      PUBLIC :: close_file
      PUBLIC :: close_inp
      PUBLIC :: close_out
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE close_file (ng, model, S, ncname, Lupdate)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      TYPE(T_IO), intent(inout) :: S
!
      logical, intent(in), optional :: Lupdate
!
      character (len=*), intent(in), optional :: ncname
!
!  Local variable declarations.
!
      integer :: ClosedState = -1
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/close_io.F"//", close_file_nf90"
!
!-----------------------------------------------------------------------
!  Close specified NetCDF file.
!-----------------------------------------------------------------------
!
      SELECT CASE (S%IOtype)
        CASE (io_nf90)
          IF (S%ncid.ne.ClosedState) THEN
            CALL netcdf_close (ng, model, S%ncid,                       &
     &                         TRIM(ncname), Lupdate)
            S%ncid=ClosedState
          END IF
      END SELECT
      IF (FoundError(exit_flag, NoError, 85, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE close_file
!
!***********************************************************************
      SUBROUTINE close_inp (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: Fcount, i
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/close_io.F"//", close_inp"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  If multi-file input fields, close several input files.
!-----------------------------------------------------------------------
!
!  Skip if configuration error.
!
      IF ((exit_flag.eq.5).or.(exit_flag.eq.6)) RETURN
!
!  If appropriate, close input forcing files and set several parameter
!  to closed state.
!
      DO i=1,nFfiles(ng)
        IF ((FRC(i,ng)%Nfiles.gt.0).and.(FRC(i,ng)%ncid.ne.-1)) THEN
          Fcount=FRC(i,ng)%Fcount
          CALL close_file (ng, model, FRC(i,ng),                        &
     &                     FRC(i,ng)%files(Fcount), .FALSE.)
          IF (FoundError(exit_flag, NoError, 125, MyFile)) RETURN
          FRCids=-1
          FRCncid=-1
          Fcount=1
          FRC(i,ng)%Fcount=Fcount
          FRC(i,ng)%name=TRIM(FRC(i,ng)%files(Fcount))
        END IF
      END DO
!
!  If appropriate, close boundary files.
!
      IF (ObcData(ng)) THEN
        DO i=1,nBCfiles(ng)
          IF ((BRY(i,ng)%Nfiles.gt.0).and.(BRY(i,ng)%ncid.ne.-1)) THEN
            Fcount=BRY(i,ng)%Fcount
            CALL close_file (ng, model, BRY(i,ng),                      &
     &                       BRY(i,ng)%files(Fcount), .FALSE.)
            IF (FoundError(exit_flag, NoError, 143, MyFile)) RETURN
            BRYids=-1
            BRYncid=-1
            Fcount=1
            BRY(i,ng)%Fcount=Fcount
            BRY(i,ng)%name=TRIM(BRY(i,ng)%files(Fcount))
          END IF
        END DO
      END IF
!
!  If appropriate, close climatology files.
!
      IF (CLM_FILE(ng)) THEN
        DO i=1,nCLMfiles(ng)
          IF ((CLM(i,ng)%Nfiles.gt.0).and.(CLM(i,ng)%ncid.ne.-1)) THEN
            Fcount=CLM(i,ng)%Fcount
            CALL close_file (ng, model, CLM(i,ng),                      &
     &                       CLM(i,ng)%files(Fcount), .FALSE.)
            IF (FoundError(exit_flag, NoError, 161, MyFile)) RETURN
            CLMids=-1
            CLMncid=-1
            Fcount=1
            CLM(i,ng)%Fcount=Fcount
            CLM(i,ng)%name=TRIM(CLM(i,ng)%files(Fcount))
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE close_inp
!
      SUBROUTINE close_out
!
!=======================================================================
!                                                                      !
! This subroutine flushes and closes all output files.                 !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      logical :: First
!
      integer :: Fcount, MyError, i, ivalue, ng
!
      character (len=256) :: ana_string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/close_io.F"//", close_out"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
!  Skip if configuration error.
!
      IF ((exit_flag.eq.5).or.(exit_flag.eq.6)) RETURN
!
!  Close output NetCDF files.
!
      DO ng=1,Ngrids
        CALL close_file (ng, iNLM, RST(ng), RST(ng)%name)
        CALL close_file (ng, iNLM, HIS(ng), HIS(ng)%name)
        CALL close_file (ng, iNLM, QCK(ng), QCK(ng)%name)
!
!  Report number of time records written.
!
        IF (Master) THEN
          WRITE (stdout,10) ng
          IF (associated(HIS(ng)%Nrec)) THEN
            IF (ANY(HIS(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'HISTORY', SUM(HIS(ng)%Nrec)
            END IF
          END IF
          IF (associated(RST(ng)%Nrec)) THEN
            Fcount=RST(ng)%load
            IF (RST(ng)%Nrec(Fcount).gt.0) THEN
              IF (LcycleRST(ng)) THEN
                IF (RST(ng)%Nrec(Fcount).gt.1) THEN
                  RST(ng)%Nrec(Fcount)=2
                ELSE
                  RST(ng)%Nrec(Fcount)=1
                END IF
              END IF
              WRITE (stdout,20) 'RESTART', RST(ng)%Nrec(Fcount)
            END IF
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report analytical header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,39
          ana_string=TRIM(ANANAME(i))
          ivalue=ICHAR(ana_string(1:1))
          IF (ivalue.ge.47) THEN          ! decimal value for characters
            IF (exit_flag.ne.5) THEN
              IF (First) THEN
                First=.FALSE.
                WRITE (stdout,30) ' Analytical header files used:'
              END IF
              WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(ANANAME(i)))
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Report biology model header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,4
          ana_string=TRIM(BIONAME(i))
          ivalue=ICHAR(ana_string(1:1))
          IF (ivalue.ge.47) THEN          ! decimal value for characters
            IF (exit_flag.ne.5) THEN
              IF (First) THEN
                First=.FALSE.
                WRITE (stdout,30) ' Biology model header files used:'
              END IF
              WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(BIONAME(i)))
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (Master.and.                                                   &
     &    (FoundError(exit_flag, NoError, 401, MyFile))) THEN
        WRITE (stdout,40) Rerror(exit_flag), exit_flag
      END IF
      IF (blowup.ne.0) THEN
        IF (Master) WRITE (stdout,50) TRIM(blowup_string)
      ELSE IF (exit_flag.eq.NoError) THEN
        CALL get_date (date_str)
        IF (Master) WRITE (stdout,60) TRIM(date_str)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,70) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,80) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,90)
      ELSE IF (exit_flag.eq.5) THEN
        IF (Master) WRITE (stdout,100)
      ELSE IF (exit_flag.eq.6) THEN
        IF (Master) WRITE (stdout,110)
      ELSE IF (exit_flag.eq.7) THEN
        IF (Master) WRITE (stdout,120)
      ELSE IF (exit_flag.eq.8) THEN
        IF (Master) WRITE (stdout,130)
      END IF
!
 10   FORMAT (/,' ROMS/TOMS - Output NetCDF summary for Grid ',         &
     &        i2.2,':')
 20   FORMAT (13x,'number of time records written in ',                 &
     &        a,' file = ',i0)
 30   FORMAT (/,a,/)
 40   FORMAT (/,a,i0,/)
 50   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.',/,               &
     &          ' REASON: ',a)
 60   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: I/O related problem.')
100   FORMAT (/,' ERROR: Illegal model configuration.')
110   FORMAT (/,' ERROR: Illegal domain partition.')
120   FORMAT (/,' ERROR: Illegal input parameter.')
130   FORMAT (/,' ERROR: Fatal algorithm result.')
!
      RETURN
      END SUBROUTINE close_out
      END MODULE close_io_mod
