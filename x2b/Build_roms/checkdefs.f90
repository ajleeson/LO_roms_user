      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 1122 2022-04-13 19:50:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: iatms, ibbl, ibiology, idriver
      integer :: ivelHadv, ivelVadv, ivmix, nearshore
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Initialze counters
!-----------------------------------------------------------------------
!
      iatms = 0
      ibbl = 0
      ibiology = 0
      idriver = 0
      ivelHadv = 0
      ivelVadv = 0
      ivmix = 0
      nearshore = 0
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t27,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
      IF (Master) WRITE (stdout,20) 'ADD_FSOBC',                        &
     &   'Adding tidal elevation to processed OBC data'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_FSOBC,'
!
      IF (Master) WRITE (stdout,20) 'ADD_M2OBC',                        &
     &   'Adding tidal currents to processed OBC data'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_M2OBC,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BPFLUX',                       &
     &   'Analytical bottom passive tracers fluxes'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BPFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_SPFLUX',                       &
     &   'Analytical surface passive tracer fluxes'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SPFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
!
      IF (Master) WRITE (stdout,20) 'BIO_FENNEL',                       &
     &   'Fennel et al. (2006) nitrogen-based model'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BIO_FENNEL,'
      ibiology=ibiology+1
!
      IF (Master) WRITE (stdout,20) 'BIO_SEDIMENT',                     &
     &   'Restore fallen particulate material to the nutrient pool'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' BIO_SEDIMENT,'
      IF (Master) WRITE (stdout,20) 'BOUNDARY_ALLREDUCE',               &
     &   'Using mpi_allreduce in mp_boundary routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' BOUNDARY_ALLREDUCE,'
!
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES',                      &
     &   'Surface bulk fluxes parameterization'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
!
      IF (Master) WRITE (stdout,20) 'CANUTO_A',                         &
     &   'Canuto A-stability function formulation'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CANUTO_A,'
!
      IF (Master) WRITE (stdout,20) 'CARBON',                           &
     &   'Add Carbon constituents'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CARBON,'
!
      IF (Master) WRITE (stdout,20) 'COLLECT_ALLREDUCE',                &
     &   'Using mpi_allreduce in mp_collect routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' COLLECT_ALLGATHER,'
!
      IF (Master) WRITE (stdout,20) 'DEFLATE',                          &
     &   'Setting compression in output NetCDF-4/HDF5 files'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DEFLATE,'
!
      IF (Master) WRITE (stdout,20) 'DENITRIFICATION',                  &
     &   'Add denitrification processes'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' DENITRIFICATION,'
!
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002)'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
!
      IF (Master) WRITE (stdout,20) 'DOUBLE_PRECISION',                 &
     &   'Double precision arithmetic numerical kernel.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
!
      IF (Master) WRITE (stdout,20) 'EMINUSP',                          &
     &   'Compute Salt Flux using E-P'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' EMINUSP,'
!
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
!
      IF (Master) WRITE (stdout,20) 'HDF5',                             &
     &   'Creating NetCDF-4/HDF5 format files'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' HDF5,'
!
      IF (Master) WRITE (stdout,20) 'LONGWAVE_OUT',                     &
     &   'Compute outgoing longwave radiation internally'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LONGWAVE_OUT,'
!
      IF (Master) WRITE (stdout,20) 'MASKING',                          &
     &   'Land/Sea masking'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
!
      IF (Master) WRITE (stdout,20) 'MPI',                              &
     &   'MPI distributed-memory configuration'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
!
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
!
      IF (Master) WRITE (stdout,20) 'NONLIN_EOS',                       &
     &   'Nonlinear Equation of State for seawater'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NONLIN_EOS,'
!
      IF (Master) WRITE (stdout,20) 'N2S2_HORAVG',                      &
     &   'Horizontal smoothing of buoyancy and shear'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
!
      IF (Master) WRITE (stdout,20) '!OCMIP_OXYGEN_SC',                 &
     &   'Schmidt number from Wanninkhof (1992)'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !OCMIP_OXYGEN_SC,'
!
      IF (Master) WRITE (stdout,20) 'OXYGEN',                           &
     &   'Add oxygen dynamics'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' OXYGEN,'
!
      IF (Master) WRITE (stdout,20) 'PCO2AIR_SECULAR',                  &
     &   'CO2 gas exchange time-dependent evolution'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' PCO2AIR_SECULAR,'
!
      IF (Master) WRITE (stdout,20) 'PERFECT_RESTART',                  &
     &   'Processing perfect restart variables'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' PERFECT_RESTART,'
!
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
!
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
!
      IF (Master) WRITE (stdout,20) '!RW14_OXYGEN_SC',                  &
     &   'Schmidt number from Wanninkhof (1992)'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' !RW14_OXYGEN_SC,'
!
      IF (Master) WRITE (stdout,20) 'K_GSCHEME',                        &
     &   'Third-order upstream advection of TKE fields'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' K_GSCHEME,'
!
      IF (Master) WRITE (stdout,20) 'RADIATION_2D',                     &
     &   'Use tangential phase speed in radiation conditions'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' RADIATION_2D,'
!
      IF (Master) WRITE (stdout,20) 'REDUCE_ALLREDUCE',                 &
     &   'Using mpi_allreduce in mp_reduce routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' REDUCE_ALLREDUCE,'
!
      IF (Master) WRITE (stdout,20) 'RI_SPLINES',                       &
     &   'Parabolic Spline Reconstruction for Richardson Number'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RI_SPLINES,'
!
      IF (Master) WRITE (stdout,20) '!RST_SINGLE',                      &
     &   'Double precision fields in restart NetCDF file'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
!
      IF (Master) WRITE (stdout,20) 'SALINITY',                         &
     &   'Using salinity'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
!
      IF (Master) WRITE (stdout,20) 'SOLAR_SOURCE',                     &
     &   'Solar Radiation Source Term'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SOLAR_SOURCE,'
!
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
!
      IF (Master) WRITE (stdout,20) 'SPHERICAL',                        &
     &   'Spherical grid configuration'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SPHERICAL,'
!
      IF (Master) WRITE (stdout,20) 'SSH_TIDES',                        &
     &   'Adding tidal elevation to SSH climatology'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SSH_TIDES,'
!
      IF (Master) WRITE (stdout,20) 'TALK_NONCONSERV',                  &
     &   'Alkalinity is affected by changes in nitrate or ammonium'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TALK_NONCONSERV,'
!
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
!
      IF (Master) WRITE (stdout,20) 'UV_COR',                           &
     &   'Coriolis term'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
!
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4VADVECTION,'
      ivelVadv=ivelVadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_QDRAG',                         &
     &   'Quadratic bottom stress'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_QDRAG,'
      ibbl=ibbl+1
!
      IF (Master) WRITE (stdout,20) 'UV_TIDES',                         &
     &   'Add tidal currents to 2D momentum climatologies'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_TIDES,'
!
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
!
      RETURN
      END SUBROUTINE checkdefs
