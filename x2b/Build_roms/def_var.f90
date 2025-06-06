      MODULE def_var_mod
!
!svn $Id: def_var.F 1120 2022-04-08 19:14:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines the requested NetCDF variable.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     Vtype        NetCDF variable type (integer)                      !
!     nVdim        Number of variable dimensions (integer; 0=scalar)   !
!     Vdim         Dimensions IDs for this variable (integer vector)   !
!     Aval         Attribute values (real vector):                     !
!                    Aval(1)   =>  Add offset value                    !
!                    Aval(2)   =>  Valid minimum value                 !
!                    Aval(3)   =>  Valid maximum value                 !
!                    Aval(4)   =>  Missing value                       !
!                    Aval(5)   =>  C-grid variable type                !
!                    Aval(6)   =>  Fill value                          !
!     Vinfo        Variable information (character array):             !
!                    Vinfo( 1) =>  Variable name                       !
!                    Vinfo( 2) =>  Variable "longname" attribute       !
!                    Vinfo( 3) =>  Variable "units" attribute          !
!                    Vinfo( 4) =>  Variable "calendar" attribute       !
!                    Vinfo( 5) =>  Variable "valid_min" attribute      !
!                    Vinfo( 6) =>  Variable "valid_max" attribute      !
!                    Vinfo( 7) =>  Variable "option_T" attribute       !
!                    Vinfo( 8) =>  Variable "option_F" attribute       !
!                    Vinfo( 9) =>  Variable "option_0" attribute       !
!                    Vinfo(10) =>  Variable "option_1" attribute       !
!                    Vinfo(11) =>  Variable "negative_value" attribute !
!                    Vinfo(12) =>  Variable "positive_value" attribute !
!                    Vinfo(13) =>  Variable "cycle" attribute          !
!                    Vinfo(14) =>  Variable "field" attribute          !
!                    Vinfo(15) =>  Variable "positions" attribute      !
!                    Vinfo(16) =>  Variable "time" attribute           !
!                    Vinfo(17) =>  Variable "missing_value" attribute  !
!                    Vinfo(18) =>  Variable "add_offset" attribute     !
!                    Vinfo(19) =>  Variable "size_class" attribute     !
!                    Vinfo(20) =>  Variable "water_points" attribute   !
!                    Vinfo(21) =>  Variable "standard_name" attribute  !
!                    Vinfo(22) =>  Variable "coordinates" attribute    !
!                    Vinfo(23) =>  Variable "formula_terms" attribute  !
!                    Vinfo(24) =>  Variable "_FillValue" attribute     !
!                    Vinfo(25) =>  Variable "positive" attribute       !
!     ncname       NetCDF file name                                    !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional)                                        !
!     SetParAccess Logical switch to set parallel I/O access flag to   !
!                    either collective or independent (optional)       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     def_var      Error flag (integer)                                !
!     Vid          NetCDF variable ID (integer)                        !
!                                                                      !
!  Notice that arrays "Aval" and "Vinfo" is destroyed on output to     !
!  facilitate the definition of the next variable.                     !
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
      INTERFACE def_var
        MODULE PROCEDURE def_var_nf90
      END INTERFACE def_var
!
      CONTAINS
!
!***********************************************************************
      FUNCTION def_var_nf90 (ng, model, ncid, Vid,                      &
     &                       Vtype, nVdim, Vdim,                        &
     &                       Aval, Vinfo, ncname,                       &
     &                       SetFillVal, SetParAccess) RESULT (status)
!***********************************************************************
!
      USE mod_netcdf
!
      USE distribute_mod, ONLY : mp_bcasti
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
      logical, intent(in), optional :: SetParAccess
!
      integer, intent (in) :: ng, model, ncid, Vtype, nVdim
      integer, dimension(:), intent (in) :: Vdim
      integer, intent (out) :: Vid
!
      real(r8), dimension(:), intent(inout) :: Aval
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(inout) :: Vinfo(25)
!
!  Local variable declarations.
!
      logical :: LandFill
      integer :: i, j, latt
      integer :: status
!
      character (len=  5) location
      character (len=160) text
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_var.F"//", def_var_nf90"
!
!-----------------------------------------------------------------------
!  Define requested variable and its attributes.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
!
      IF (OutThread) THEN
!
!  Define variable.
!
        IF (exit_flag.eq.NoError) THEN
          IF (LEN_TRIM(Vinfo(1)).gt.0) THEN
            IF ((nVdim.eq.1).and.(Vdim(1).eq.0)) THEN
              status=nf90_def_var(ncid, TRIM(Vinfo(1)), Vtype,          &
     &                            varid = Vid)
            ELSE
              status=nf90_def_var(ncid, TRIM(Vinfo(1)), Vtype,          &
     &                            Vdim(1:nVdim), Vid)
            END IF
            IF (FoundError(status, nf90_noerr, 164, MyFile)) THEN
              IF (Master) WRITE (stdout,10) ng, TRIM(Vinfo(1)),         &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Define deflate (file compresion) parameters. Notice that deflation
!  cannot be used in parallel I/O for writing data. This is because
!  the compression makes it impossible for the  library to exactly
!  map the data to the disk location. However, deflated data can be
!  read with parallel I/O
!
        IF (exit_flag.eq.NoError) THEN
          IF (LEN_TRIM(Vinfo(1)).gt.0) THEN
            IF ((nVdim.gt.1).and.(Vdim(1).ne.0)) THEN
              status=nf90_def_var_deflate (ncid, Vid, shuffle, deflate, &
     &                                     deflate_level)
              IF (FoundError(status, nf90_noerr, 186, MyFile)) THEN
                IF (Master) WRITE (stdout,20) TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
!  Define special attributes for SGRID conventions variable "grid".
!
        IF (exit_flag.eq.NoError) THEN
          IF (TRIM(Vinfo(1)).eq.'grid') THEN
            status=nf90_put_att(ncid, Vid, 'cf_role',                   &
     &                          'grid_topology')
            IF (FoundError(status, nf90_noerr, 203, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'cf_role',                  &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'topology_dimension',        &
     &                          (/2/))
            IF (FoundError(status, nf90_noerr, 213, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'topology_dimension',       &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'node_dimensions',           &
     &                          'xi_psi eta_psi')
            IF (FoundError(status, nf90_noerr, 223, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'node_dimensions',          &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_rho: xi_psi (padding: both) '//                    &
     &           'eta_rho: eta_psi (padding: both)'
            status=nf90_put_att(ncid, Vid, 'face_dimensions',           &
     &                          TRIM(text))
            IF (FoundError(status, nf90_noerr, 235, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'face_dimensions',          &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_u: xi_psi eta_u: eta_psi (padding: both)'
            status=nf90_put_att(ncid, Vid, 'edge1_dimensions',          &
     &                          TRIM(text))
            IF (FoundError(status, nf90_noerr, 246, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'edge1_dimensions',         &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            text='xi_v: xi_psi (padding: both) eta_v: eta_psi'
            status=nf90_put_att(ncid, Vid, 'edge2_dimensions',          &
     &                          TRIM(text))
            IF (FoundError(status, nf90_noerr, 257, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'edge2_dimensions',         &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            IF (spherical) THEN
              status=nf90_put_att(ncid, Vid, 'node_coordinates',        &
     &                            'lon_psi lat_psi')
              IF (FoundError(status, nf90_noerr, 268, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'node_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'face_coordinates',        &
     &                            'lon_rho lat_rho')
              IF (FoundError(status, nf90_noerr, 278, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'face_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'edge1_coordinates',       &
     &                            'lon_u lat_u')
              IF (FoundError(status, nf90_noerr, 288, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'edge1_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'edge2_coordinates',       &
     &                            'lon_v lat_v')
              IF (FoundError(status, nf90_noerr, 298, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'edge2_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              status=nf90_put_att(ncid, Vid, 'node_coordinates',        &
     &                            'x_psi y_psi')
              IF (FoundError(status, nf90_noerr, 308, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'node_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'face_coordinates',        &
     &                            'x_rho y_rho')
              IF (FoundError(status, nf90_noerr, 318, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'face_coordinates',       &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'edge1_coordinates',       &
     &                            'x_u y_u')
              IF (FoundError(status, nf90_noerr, 328, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'edge1_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
!
              status=nf90_put_att(ncid, Vid, 'edge2_coordinates',       &
     &                            'x_v y_v')
              IF (FoundError(status, nf90_noerr, 338, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'edge2_coordinates',      &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
!
            status=nf90_put_att(ncid, Vid, 'vertical_dimensions',       &
     &                          's_rho: s_w (padding: none)')
            IF (FoundError(status, nf90_noerr, 350, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'vertical_dimensions',      &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "standard_name" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(21))
          IF (latt.gt.0.and.(Vinfo(21)(1:6).ne.'nulval')) THEN
            status=nf90_put_att(ncid, Vid, 'standard_name',             &
     &                          Vinfo(21)(1:latt))
            IF (FoundError(status, nf90_noerr, 368, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'standard_name',            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Define "long_name" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          IF (LEN_TRIM(Vinfo(2)).gt.0) THEN
            vinfo(2)=TRIM(ADJUSTL(vinfo(2)))
            latt=LEN_TRIM(Vinfo(2))
            status=nf90_put_att(ncid, Vid, 'long_name',                 &
     &                          Vinfo(2)(1:latt))
            IF (FoundError(status, nf90_noerr, 386, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'long_name',                &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Define "size_class" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(19))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'size_class',                &
     &                          Vinfo(19)(1:latt))
            IF (FoundError(status, nf90_noerr, 403, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'size_class',               &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "units" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(3))
          IF (latt.gt.0) THEN
            IF (TRIM(Vinfo(3)).ne.'nondimensional') THEN
              status=nf90_put_att(ncid, Vid, 'units',                   &
     &                            Vinfo(3)(1:latt))
              IF (FoundError(status, nf90_noerr, 421, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'units',                  &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END IF
        END IF
!
!  If applicable, define "calendar" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(4))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'calendar',                  &
     &                          Vinfo(4)(1:latt))
            IF (FoundError(status, nf90_noerr, 439, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'calendar',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "valid_min" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(5))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            INT(Aval(2)))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            REAL(Aval(2),r4))
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(5)),            &
     &                            Aval(2))
            END IF
            IF (FoundError(status, nf90_noerr, 466, MyFile)) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(5)),             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(2)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "valid_max" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(6))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            INT(Aval(3)))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            REAL(Aval(3),r4))
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(6)),            &
     &                            Aval(3))
            END IF
            IF (FoundError(status, nf90_noerr, 494, MyFile)) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(6)),             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(3)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "flag_values" and "flag_meanings" attributes
!  for logical variables.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((LEN_TRIM(Vinfo(7)).gt.0).and.                            &
     &        (LEN_TRIM(Vinfo(8)).gt.0)) THEN
            text='T, F'
            latt=LEN_TRIM(text)
            status=nf90_put_att(ncid, Vid, 'flag_values',               &
     &                          text(1:latt))
            IF (status.eq.nf90_noerr) THEN
              text=TRIM(Vinfo(7))//' '//TRIM(Vinfo(8))
              latt=LEN_TRIM(text)
              status=nf90_put_att(ncid, Vid, 'flag_meanings',           &
     &                            text(1:latt))
              IF (FoundError(status, nf90_noerr, 520, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'flag_meanings (T/F)',    &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,30) 'flag_values (T/F)',        &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "flag_values" and "flag_meanings" attributes
!  for integer and floating point variables.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((LEN_TRIM(Vinfo( 9)).gt.0).and.                           &
     &        (LEN_TRIM(Vinfo(10)).gt.0)) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0,1/))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0.0_r4, 1.0_r4/))
            ELSE
              status=nf90_put_att(ncid, Vid, 'flag_values',             &
     &                            (/0.0_r8, 1.0_r8/))
            END IF
            IF (status.eq.nf90_noerr) THEN
              text=TRIM(Vinfo(9))//' '//TRIM(Vinfo(10))
              latt=LEN_TRIM(text)
              status=nf90_put_att(ncid, Vid, 'flag_meanings',           &
     &                            text(1:latt))
              IF (FoundError(status, nf90_noerr, 560, MyFile)) THEN
                IF (Master) WRITE (stdout,30) 'flag_meanings',          &
     &                                        TRIM(Vinfo(1)),           &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,30) 'flag_values',              &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "negative_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(11))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'negative_value',            &
     &                          Vinfo(11)(1:latt))
            IF (FoundError(status, nf90_noerr, 584, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'negative_value',           &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "positive_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(12))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positive_value',            &
     &                          Vinfo(12)(1:latt))
            IF (FoundError(status, nf90_noerr, 601, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'positive_value',           &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define CF "positive" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(25))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positive',                  &
     &                          Vinfo(25)(1:latt))
            IF (FoundError(status, nf90_noerr, 618, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'positive',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "cycle" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(13))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'cycle',                     &
     &                          Vinfo(13)(1:latt))
            IF (FoundError(status, nf90_noerr, 635, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'cycle',                    &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "positions" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(15))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'positions',                 &
     &                          Vinfo(15)(1:latt))
            IF (FoundError(status, nf90_noerr, 652, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'positions',                &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "time" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(16))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'time',                      &
     &                          Vinfo(16)(1:latt))
            IF (FoundError(status, nf90_noerr, 669, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'time',                     &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "cell_methods" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          IF (((ABS(INT(Aval(5))).le.4).and.(nVdim.gt.2)).or.           &
     &        ((ABS(INT(Aval(5))).gt.4).and.(nVdim.gt.3))) THEN
            IF ((ncid.eq.AVG(ng)%ncid).or.                              &
     &          (ncid.eq.DIA(ng)%ncid)) THEN
              text='ocean_time: mean'
            ELSE
              text='ocean_time: point'
            END IF
            latt=LEN_TRIM(text)
            status=nf90_put_att(ncid, Vid, 'cell_methods',              &
     &                          text(1:latt))
            IF (FoundError(status, nf90_noerr, 693, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'cell_methods',             &
     &                                      text(1:latt),               &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "missing_value" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(17))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            INT(Aval(4)))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            REAL(Aval(4),r4))
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(17)),           &
     &                            Aval(4))
            END IF
            IF (FoundError(status, nf90_noerr, 720, MyFile)) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(17)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(4)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "add_offset" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(18))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            INT(Aval(1)))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            REAL(Aval(1),r4))
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(18)),           &
     &                            Aval(1))
            END IF
            IF (FoundError(status, nf90_noerr, 748, MyFile)) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(18)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(1)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "_FillValue" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(24))
          IF (latt.gt.0) THEN
            IF (Vtype.eq.nf90_int) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            INT(Aval(6)))
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            REAL(Aval(6),r4))
            ELSE
              status=nf90_put_att(ncid, Vid, TRIM(Vinfo(24)),           &
     &                            Aval(6))
            END IF
            IF (FoundError(status, nf90_noerr, 776, MyFile)) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vinfo(24)),            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(6)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "water_points" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(20))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'water_points',              &
     &                          Vinfo(20)(1:latt))
            IF (FoundError(status, nf90_noerr, 794, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'water_points',             &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "coordinates" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(22))
          IF (latt.gt.0) THEN
            IF (spherical) THEN
              IF (INT(Aval(5)).eq.r2dvar) THEN
                text='lon_rho lat_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.r3dvar) THEN
                text='lon_rho lat_rho s_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.w3dvar) THEN
                text='lon_rho lat_rho s_w'
                location='face'
              ELSE IF (INT(Aval(5)).eq.b3dvar) THEN
                text='lon_rho lat_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.u2dvar) THEN
                text='lon_u lat_u'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.u3dvar) THEN
                text='lon_u lat_u s_rho'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.-u3dvar) THEN
                text='lon_u lat_u s_w'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.v2dvar) THEN
                text='lon_v lat_v'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.v3dvar) THEN
                text='lon_v lat_v s_rho'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.-v3dvar) THEN
                text='lon_v lat_v s_w'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.p2dvar) THEN
                text='lon_psi lat_psi'
                location='node'
              ELSE IF (INT(Aval(5)).eq.p3dvar) THEN
                text='lon_psi lat_psi s_rho'
                location='node'
              ELSE IF (INT(Aval(5)).eq.l3dvar) THEN
                text='lon_rho lat_rho light'
                location='face'
              ELSE IF (INT(Aval(5)).eq.l4dvar) THEN
                text='lon_rho lat_rho s_rho light'
                location='face'
              END IF
            ELSE
              IF (INT(Aval(5)).eq.r2dvar) THEN
                text='x_rho y_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.r3dvar) THEN
                text='x_rho y_rho s_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.w3dvar) THEN
                text='x_rho y_rho s_w'
                location='face'
              ELSE IF (INT(Aval(5)).eq.b3dvar) THEN
                text='x_rho y_rho'
                location='face'
              ELSE IF (INT(Aval(5)).eq.u2dvar) THEN
                text='x_u y_u'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.u3dvar) THEN
                text='x_u y_u s_rho'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.-u3dvar) THEN
                text='x_u y_u s_w'
                location='edge1'
              ELSE IF (INT(Aval(5)).eq.v2dvar) THEN
                text='x_v y_v'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.v3dvar) THEN
                text='x_v y_v s_rho'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.-v3dvar) THEN
                text='x_v y_v s_w'
                location='edge2'
              ELSE IF (INT(Aval(5)).eq.p2dvar) THEN
                text='x_psi y_psi'
                location='node'
              ELSE IF (INT(Aval(5)).eq.p3dvar) THEN
                text='x_psi y_psi s_rho'
                location='node'
              ELSE IF (INT(Aval(5)).eq.l3dvar) THEN
                text='x_rho y_rho Nbands'
                location='face'
              ELSE IF (INT(Aval(5)).eq.l4dvar) THEN
                text='x_rho y_rho s_rho Nbands'
                location='face'
              END IF
            END IF
            status=nf90_put_att(ncid, Vid, 'grid',                      &
     &                          'grid')
            IF (FoundError(status, nf90_noerr, 901, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'grid',                     &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            status=nf90_put_att(ncid, Vid, 'location',                  &
     &                          TRIM(location))
            IF (FoundError(status, nf90_noerr, 911, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'location',                 &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
!
            latt=LEN_TRIM(text)
            IF (((ABS(INT(Aval(5))).le.4).and.(nVdim.gt.2)).or.         &
     &          ((ABS(INT(Aval(5))).gt.4).and.(nVdim.gt.3))) THEN
              text=text(1:latt)//' ocean_time'
              latt=LEN_TRIM(text)
            END IF
            status=nf90_put_att(ncid, Vid, TRIM(Vinfo(22)),             &
     &                          text(1:latt))
            IF (FoundError(status, nf90_noerr, 927, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'coordinates',              &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            Aval(5)=0.0_r8
          END IF
        END IF
!
!  If applicable, define "formula_terms" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(23))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'formula_terms',             &
     &                          Vinfo(23)(1:latt))
            IF (FoundError(status, nf90_noerr, 945, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'formula_terms',            &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If applicable, define "field" attribute (always last).
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(14))
          IF (latt.gt.0) THEN
            status=nf90_put_att(ncid, Vid, 'field',                     &
     &                          Vinfo(14)(1:latt))
            IF (FoundError(status, nf90_noerr, 962, MyFile)) THEN
              IF (Master) WRITE (stdout,30) 'field',                    &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  If land/sea masking, define "_FillValue" attribute since masking
!  areas are overwritten with special value (spval) during output.
!  Notice that the coordinate attribute is used to check which
!  variables need the "_FillValue" attribute.
!
        IF (exit_flag.eq.NoError) THEN
          latt=LEN_TRIM(Vinfo(22))
          IF (PRESENT(SetFillVal)) THEN
            LandFill=SetFillVal
          ELSE
            LandFill=(latt.gt.0).and.(nVdim.gt.2)
          END IF
          IF (LandFill) THEN
            IF (Vtype.eq.nf90_double) THEN
              status=nf90_put_att(ncid, Vid, '_FillValue',              &
     &                            spval)
            ELSE IF (Vtype.eq.nf90_float) THEN
              status=nf90_put_att(ncid, Vid, '_FillValue',              &
     &                            REAL(spval,r4))
            END IF
            IF (FoundError(status, nf90_noerr, 996, MyFile)) THEN
              IF (Master) WRITE (stdout,30) '_FillValue',               &
     &                                      TRIM(Vinfo(1)),             &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
      END IF
!
!  Clean information variables.
!
      DO i=1,SIZE(Vinfo)
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
!
!  Broadcast error flag.
!
      CALL mp_bcasti (ng, model, exit_flag)
!
 10   FORMAT (/,' DEF_VAR_NF90 - Grid ',i2.2,                           &
     &          ', unable to define variable: ',a,/,                    &
     &        16x,'in NetCDF file: ',a)
 20   FORMAT (/,' DEF_VAR_NF90 - error while setting deflate',          &
     &        ' parameters for variable: ',a,/,16x,'in NetCDF file: ',a)
 30   FORMAT (/,'DEF_VAR_NF90 - error while defining attribute: ',a,    &
     &        ' for variable: ',a,/,16x,'in NetCDF file: ',a)
!
      RETURN
      END FUNCTION def_var_nf90
      END MODULE def_var_mod
