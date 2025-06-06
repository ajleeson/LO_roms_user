      MODULE distribute_mod
!
!svn $Id: distribute.F 1110 2022-02-27 21:37:32Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used for distrubuted-memory communications       !
!  between parallel nodes:                                             !
!                                                                      !
!  mp_aggregate2d    aggregates 2D tiled data into a 2D global array   !
!  mp_aggregate3d    aggregates 3D tiled data into a 3D global array   !
!  mp_barrier        barrier sychronization                            !
!  mp_bcastf         broadcasts floating point variables               !
!  mp_bcasti         broadcasts integer variables                      !
!  mp_bcastl         broadcasts logical variables                      !
!  mp_bcasts         broadcasts character variables                    !
!  mp_bcast_struc    broadcats NetCDF IDs of an IO_TYPE structure      !
!  mp_boundary       exchanges boundary data between tiles             !
!  mp_assemblef_1d   assembles 1D floating point array from tiles      !
!  mp_assemblef_2d   assembles 2D floating point array from tiles      !
!  mp_assemblef_3d   assembles 3D floating point array from tiles      !
!  mp_assemblei_1d   assembles 1D integer array from tiles             !
!  mp_assemblei_2d   assembles 2D integer array from tiles             !
!  mp_collect_f      collects 1D floating point array from tiles       !
!  mp_collect_i      collects 1D integer array from tiles              !
!  mp_dump           writes 2D and 3D tiles arrays for debugging       !
!  mp_gather2d       collects a 2D tiled array for output purposes     !
!  mp_gather3d       collects a 3D tiled array for output purposes     !
!  mp_gather_state   collects state vector for unpacking of variables  !
!  mp_ncread1d       reads  in  1D state array from NetCDF file        !
!  mp_ncread2d       reads  in  2D state array from NetCDF file        !
!  mp_ncwrite1d      writes out 1D state array into NetCDF file        !
!  mp_ncwrite2d      writes out 2D state array into NetCDF file        !
!  mp_reduce         global reduction operations                       !
!  mp_reduce2        global reduction operations (MINLOC, MAXLOC)      !
!  mp_scatter2d      scatters input data to a 2D tiled array           !
!  mp_scatter3d      scatters input data to a 3D tiled array           !
!  mp_scatter_state  scatters global data for packing of state vector  !
!                                                                      !
!  Notice that the tile halo exchange can be found in "mp_exchange.F"  !
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
      INTERFACE mp_assemble
        MODULE PROCEDURE mp_assemblef_1d
        MODULE PROCEDURE mp_assemblef_2d
        MODULE PROCEDURE mp_assemblef_3d
        MODULE PROCEDURE mp_assemblei_1d
        MODULE PROCEDURE mp_assemblei_2d
      END INTERFACE mp_assemble
!
      INTERFACE mp_bcastf
        MODULE PROCEDURE mp_bcastf_0d
        MODULE PROCEDURE mp_bcastf_1d
        MODULE PROCEDURE mp_bcastf_2d
        MODULE PROCEDURE mp_bcastf_3d
        MODULE PROCEDURE mp_bcastf_4d
      END INTERFACE mp_bcastf
!
      INTERFACE mp_bcastl
        MODULE PROCEDURE mp_bcastl_0d
        MODULE PROCEDURE mp_bcastl_1d
        MODULE PROCEDURE mp_bcastl_2d
      END INTERFACE mp_bcastl
!
      INTERFACE mp_bcasti
        MODULE PROCEDURE mp_bcasti_0d
        MODULE PROCEDURE mp_bcasti_1d
        MODULE PROCEDURE mp_bcasti_2d
      END INTERFACE mp_bcasti
!
      INTERFACE mp_bcasts
        MODULE PROCEDURE mp_bcasts_0d
        MODULE PROCEDURE mp_bcasts_1d
        MODULE PROCEDURE mp_bcasts_2d
        MODULE PROCEDURE mp_bcasts_3d
      END INTERFACE mp_bcasts
!
      INTERFACE mp_collect
        MODULE PROCEDURE mp_collect_f
        MODULE PROCEDURE mp_collect_i
      END INTERFACE mp_collect
!
      INTERFACE mp_reduce
        MODULE PROCEDURE mp_reduce_i8    ! integer reduction
        MODULE PROCEDURE mp_reduce_0d
        MODULE PROCEDURE mp_reduce_1d
      END INTERFACE mp_reduce
!
      CONTAINS
!
      SUBROUTINE mp_barrier (ng, model, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine blocks the caller until all group members have called  !
!  it.                                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
!  Local variable declarations.
!
      integer ::  MyCOMM, MyError
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_barrier"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 72, 146, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize all distribute-memory nodes in the group.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier (MyCOMM, MyError)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 72, 174, MyFile)
!
      RETURN
      END SUBROUTINE mp_barrier
!
      SUBROUTINE mp_bcastf_0d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a floating-point scalar variable to all     !
!  processors in the communicator. It is called by all the members     !
!  in the group.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (real).                         !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastf_0d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 404, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=1
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 442, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastf_0d
!
      SUBROUTINE mp_bcastf_1d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D floating-point, non-tiled, array       !
!  to all processors in the communicator. It is called by all the      !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (real).                         !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastf_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 492, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 531, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastf_1d
!
      SUBROUTINE mp_bcastf_2d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D floating-point, non-tiled, array       !
!  to all processors in the communicator. It is called by all the      !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (real).                         !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: Asize(2)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastf_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 583, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 624, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastf_2d
!
      SUBROUTINE mp_bcastf_3d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 3D floating-point, non-tiled, array       !
!  to all processors in the communicator. It is called by all the      !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          3D array to broadcast (real).                         !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 3D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(inout) :: A(:,:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: Asize(3)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastf_3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 676, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Npts=Asize(1)*Asize(2)*Asize(3)
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_3D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 718, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastf_3d
!
      SUBROUTINE mp_bcastf_4d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 4D floating-point, non-tiled, array       !
!  to all processors in the communicator. It is called by all the      !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          4D array to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 4D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(inout) :: A(:,:,:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: Asize(4)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastf_4d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 769, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Asize(4)=UBOUND(A, DIM=4)
      Npts=Asize(1)*Asize(2)*Asize(3)*Asize(4)
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_4D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 812, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastf_4d
!
      SUBROUTINE mp_bcasti_0d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts an integer scalar variable to all           !
!  processors in the communicator.  It is called by all the            !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (integer).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
      integer, intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasti_0d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_on (ng, model, 64, 863, MyFile)
      END IF
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=1
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_off (ng, model, 64, 903, MyFile)
      END IF
!
      RETURN
      END SUBROUTINE mp_bcasti_0d
!
      SUBROUTINE mp_bcasti_1d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D non-tiled, integer array to all        !
!  processors in the communicator. It is called by all the             !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (integer).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
      integer, intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasti_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 954, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 993, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcasti_1d
!
      SUBROUTINE mp_bcasti_2d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D non-tiled, integer array to all        !
!  processors in the communicator. It is called by all the             !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (integer).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
      integer, intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: Asize(2)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasti_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1044, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1085, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcasti_2d
!
      SUBROUTINE mp_bcastl_0d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a logical scalar variable to all            !
!  processors in the communicator. It is called by all the             !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (logical).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      logical, intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastl_0d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1135, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=1
      CALL mpi_bcast (A, Npts, MPI_LOGICAL, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1173, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastl_0d
!
      SUBROUTINE mp_bcastl_1d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D nontiled, logical array to all         !
!  processors in the communicator. It is called by all the             !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (logical).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      logical, intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastl_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1223, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)
      CALL mpi_bcast (A, Npts, MPI_LOGICAL, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1262, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastl_1d
!
      SUBROUTINE mp_bcastl_2d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D non-tiled, logical array to all        !
!  processors in the communicator. It is called by all the             !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (logical).                      !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      logical, intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: Asize(2)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcastl_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1313, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)
      CALL mpi_bcast (A, Npts, MPI_LOGICAL, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1354, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcastl_2d
!
      SUBROUTINE mp_bcasts_0d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a string scalar variable to all processors  !
!  in the communicator. It is called by all the members in the group.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (string).                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Nchars, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasts_0d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_on (ng, model, 64, 1404, MyFile)
      END IF
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Nchars=LEN(A)
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_off (ng, model, 64, 1444, MyFile)
      END IF
!
      RETURN
      END SUBROUTINE mp_bcasts_0d
!
      SUBROUTINE mp_bcasts_1d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D string array to all processors in the  !
!  communicator. It is called by all the members in the group.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (string).                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Asize, Lstr, MyCOMM, MyError, Nchars, Serror
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasts_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1494, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize=UBOUND(A, DIM=1)
      Nchars=LEN(A(1))*Asize
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1534, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcasts_1d
!
      SUBROUTINE mp_bcasts_2d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D string array to all processors in the  !
!  communicator. It is called by all the members in the group.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (string).                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Nchars, Serror
      integer :: Asize(2)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasts_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1584, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Nchars=LEN(A(1,1))*Asize(1)*Asize(2)
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1625, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcasts_2d
!
      SUBROUTINE mp_bcasts_3d (ng, model, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 3D string array to all processors in the  !
!  communicator. It is called by all the members in the group.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          3D array to broadcast (string).                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 3D array.                                 !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(inout) :: A(:,:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Nchars, Serror
      integer :: Asize(3)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcasts_3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1675, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Nchars=LEN(A(1,1,1))*Asize(1)*Asize(2)*Asize(3)
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_3D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1717, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcasts_3d
!
      SUBROUTINE mp_bcast_struc (ng, model, S, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts the NetCDF IDs of a TYPE_IO structure to    !
!  all processors in the communicator. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     S          ROMS I/O structure, TYPE(T_IO).                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     S          Broadcasted ROMS I/O structure.                       !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: InpComm
!
      TYPE(T_IO), intent(inout) :: S(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyCOMM, MyError, Nchars, Npts, Serror
      integer :: ibuffer(5)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_bcast_struc"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 64, 1768, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast variables in structure.
!-----------------------------------------------------------------------
!
!  Structure scalar integer variables.
!
      ibuffer(1)=S(ng)%Nfiles
      ibuffer(2)=S(ng)%Fcount
      ibuffer(3)=S(ng)%load
      ibuffer(4)=S(ng)%Rindex
      ibuffer(5)=S(ng)%ncid
!
      Npts=5
      CALL mpi_bcast (ibuffer, Npts, MPI_INTEGER, MyMaster,             &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCAST_STRUC - error during ',a,' call, Node = ', &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      ELSE
        S(ng)%Nfiles=ibuffer(1)
        S(ng)%Fcount=ibuffer(2)
        S(ng)%load  =ibuffer(3)
        S(ng)%Rindex=ibuffer(4)
        S(ng)%ncid  =ibuffer(5)
      END IF
!
!  Variables IDs.
!
      Npts=UBOUND(S(ng)%Vid, DIM=1)
      CALL mpi_bcast (S(ng)%Vid, Npts, MPI_INTEGER, MyMaster,           &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Tracer variables IDs.
!
      Npts=UBOUND(S(ng)%Tid, DIM=1)
      CALL mpi_bcast (S(ng)%Tid, Npts, MPI_INTEGER, MyMaster,           &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Structure Filenames.
!
      Nchars=LEN(S(ng)%head)
      CALL mpi_bcast (S(ng)%head, Nchars, MPI_BYTE, MyMaster,           &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
      Nchars=LEN(S(ng)%base)
      CALL mpi_bcast (S(ng)%base, Nchars, MPI_BYTE, MyMaster,           &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
      Nchars=LEN(S(ng)%name)
      CALL mpi_bcast (S(ng)%name, Nchars, MPI_BYTE, MyMaster,           &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
      Nchars=LEN(S(ng)%files(1))*S(ng)%Nfiles
      CALL mpi_bcast (S(ng)%files, Nchars, MPI_BYTE, MyMaster,          &
     &                MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 64, 1897, MyFile)
!
      RETURN
      END SUBROUTINE mp_bcast_struc
!
      SUBROUTINE mp_boundary (ng, model, Imin, Imax,                    &
     &                        LBi, UBi, LBk, UBk,                       &
     &                        update, A)
!
!***********************************************************************
!                                                                      !
!  This routine exchanges boundary arrays between tiles.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Imin       Starting tile index.                                  !
!     Imax       Ending   tile index.                                  !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     update     Switch activated by the node that updated the         !
!                  boundary data.                                      !
!     A          Boundary array (1D or 2D) to process.                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Updated boundary array (1D or 2D).                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: update
!
      integer, intent(in) :: ng, model, Imin, Imax
      integer, intent(in) :: LBi, UBi, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:UBi,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Ilen, Ioff, Lstr, MyError, Nnodes, Npts, Serror
      integer :: i, ik, k, kc, rank
!
      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1)) :: Asend
      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1)) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_boundary"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 68, 1970, MyFile)
!
!-----------------------------------------------------------------------
!  Pack boundary data.  Zero-out boundary array except points updated
!  by the appropriate node, so sum reduction can be perfomed during
!  unpacking.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Asend)+                   &
     &                                   SIZE(Arecv))*KIND(A),r8))
!
!  Initialize buffer to the full range so unpacking is correct with
!  summation.  This also allows even exchange of segments with
!  communication routine "mpi_allgather".
!
      Ilen=UBi-LBi+1
      Ioff=1-LBi
      Npts=Ilen*(UBk-LBk+1)
      DO i=1,Npts
        Asend(i)=0.0_r8
      END DO
!
!  If a boundary tile, load boundary data.
!
      IF (update) THEN
        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            ik=i+Ioff+kc
            Asend(ik)=A(i,k)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
 10     FORMAT (/,' MP_BOUNDARY - error during ',a,' call, Node = ',    &
     &          i3.3,' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data: reduction sum.
!-----------------------------------------------------------------------
!
      ik=0
      DO k=LBk,UBk
        DO i=LBi,UBi
          ik=ik+1
          A(i,k)=Arecv(ik)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 68, 2072, MyFile)
!
      RETURN
      END SUBROUTINE mp_boundary
!
      SUBROUTINE mp_assemblef_1d (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 1D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          1D array to collect.                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 1D array.                                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(in) :: Aspv
      real(r8), intent(inout) :: A(:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, MyNpts, Nnodes, Serror
      integer :: i, rank, request
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      real(r8), dimension(Npts) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_assemblef_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 70, 2139, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(Npts*KIND(A),r8))
!
      MyNpts=UBOUND(A, DIM=1)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!
!  Coppy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MP_FLOAT, MPI_SUM,            &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 70, 2301, MyFile)
!
 10   FORMAT (/,' MP_ASSEMBLEF_1D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_1D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_1D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)
!
      RETURN
      END SUBROUTINE mp_assemblef_1d
!
      SUBROUTINE mp_assemblef_2d (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 2D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          2D array to collect.                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 2D array.                                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(in) :: Aspv
      real(r8), intent(inout) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, MyNpts, Nnodes, Serror
      integer :: i, rank, request
      integer :: Asize(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      real(r8), dimension(Npts) :: Arecv, Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_assemblef_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 70, 2381, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(2*Npts*KIND(A),r8))
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      MyNpts=Asize(1)*Asize(2)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Reshape input 2D data into 1D array to facilitate communications.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data into output 2D array.
!
      A=RESHAPE(Arecv, Asize)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 70, 2557, MyFile)
!
 10   FORMAT (/,' MP_ASSEMBLEF_2D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_2D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_2D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)
!
      RETURN
      END SUBROUTINE mp_assemblef_2d
!
      SUBROUTINE mp_assemblef_3d (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 3D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          3D array to collect.                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 3D array.                                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
      real(r8), intent(in) :: Aspv
      real(r8), intent(inout) :: A(:,:,:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, MyNpts, Nnodes, Serror
      integer :: i, rank, request
      integer :: Asize(3)
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      real(r8), dimension(Npts) :: Arecv, Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_assemblef_3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 70, 2637, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(2*Npts*KIND(A),r8))
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      MyNpts=Asize(1)*Asize(2)*Asize(3)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Reshape input 3D data into 1D array to facilitate communications.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data into output 3D array.
!
      A=RESHAPE(Arecv, Asize)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 70, 2814, MyFile)
!
 10   FORMAT (/,' MP_ASSEMBLEF_3D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_3D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_3D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)
!
      RETURN
      END SUBROUTINE mp_assemblef_3d
!
      SUBROUTINE mp_assemblei_1d (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 1D integer array from all members in the   !
!  group.  The collection of data from all nodes is achieved as a      !
!  reduction sum.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          1D array to collect.                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 1D array.                                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
      integer, intent(in) :: Aspv
      integer, intent(inout) :: A(:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, MyNpts, Nnodes, Serror
      integer :: i, rank, request
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(Npts) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_assemblei_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 70, 2890, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(Npts*KIND(A),r8))
!
      MyNpts=UBOUND(A, DIM=1)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MPI_INTEGER, MPI_SUM,         &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 70, 3052, MyFile)
!
 10   FORMAT (/,' MP_ASSEMBLEI_1D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEI_1D - illegal special value, Aspv = ',i4, &
     &        /,19x,'a zero value is needed for global reduction.')
 30   FORMAT (/,' MP_ASSEMBLEI_1D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)
!
      RETURN
      END SUBROUTINE mp_assemblei_1d
!
      SUBROUTINE mp_assemblei_2d (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 2D integer array from all members in the   !
!  group.  The collection of data from all nodes is achieved as a      !
!  reduction sum.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          2D array to collect.                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 2D array.                                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
      integer, intent(in) :: Aspv
      integer, intent(inout) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: Lstr,  MyCOMM, MyError, MyNpts, Nnodes, Serror
      integer :: i, rank, request
      integer :: Asize(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(Npts) :: Arecv, Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_assemblei_2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 70, 3131, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(2*Npts*KIND(A),r8))
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      MyNpts=Asize(1)*Asize(2)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Reshape input 2D data into 1D array to facilitate communications.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MPI_INTEGER, MPI_SUM,     &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data.
!
      A=RESHAPE(Arecv, Asize)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 70, 3308, MyFile)
!
 10   FORMAT (/,' MP_ASSEMBLEI_2D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEI_2D - illegal special value, Aspv = ',i4, &
     &        /,19x,'a zero value is needed for global reduction.')
 30   FORMAT (/,' MP_ASSEMBLEI_2D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)
!
      RETURN
      END SUBROUTINE mp_assemblei_2d
!
      SUBROUTINE mp_collect_f (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 1D floating-point array from all members    !
!  in the group. Then, it packs distributed data by removing the       !
!  special values. This routine is used when extracting station        !
!  data from tiled arrays.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points.                      !
!     Aspv       Special value indicating no data.  This implies that  !
!                  desired data is tile unbouded.                      !
!     A          Collected data.                                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in), optional :: InpComm
!
      real(r8), intent(in) :: Aspv
      real(r8), intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Nnodes, Serror
      integer :: i, rank, request
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      real(r8), dimension(Npts) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_collect_f"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 69, 3383, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(Npts*KIND(A),r8))
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MP_FLOAT, MPI_SUM,            &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
 10   FORMAT (/,' MP_COLLECT_F - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,14x,a)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 69, 3527, MyFile)
!
      RETURN
      END SUBROUTINE mp_collect_f
!
      SUBROUTINE mp_collect_i (ng, model, Npts, Aspv, A, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 1D integer array from all members in        !
!  the group. Then, it packs distributed data by removing the          !
!  special values. This routine is used when extracting station        !
!  data from tiled arrays.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points.                      !
!     Aspv       Special value indicating no data.  This implies that  !
!                  desired data is tile unbouded.                      !
!     A          Collected data.                                       !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts
      integer, intent(in) :: Aspv
      integer, intent(in), optional :: InpComm
      integer, intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Nnodes, Serror
      integer :: i, rank, request
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(Npts) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_collect_i"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 69, 3594, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(Npts*KIND(A),r8))
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MPI_INTEGER, MPI_SUM,         &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
 10   FORMAT (/,' MP_COLLECT_I - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,14x,a)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 69, 3738, MyFile)
!
      RETURN
      END SUBROUTINE mp_collect_i
!
      SUBROUTINE mp_gather2d (ng, model, LBi, UBi, LBj, UBj,            &
     &                        tindex, gtype, Ascl,                      &
     &                        Amask,                                    &
     &                        A, Npts, Awrk, SetFillVal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 2D tiled, floating-point array from each    !
!  spawned node and stores it into one dimensional global array. It    !
!  is used to collect and  pack output data.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     tindex     Time record index to process.                         !
!     gtype      C-grid type. If negative and Land-Sea is available,   !
!                  only water-points processed.                        !
!     Ascl       Factor to scale field before writing.                 !
!     Amask      Land/Sea mask, if any.                                !
!     A          2D tiled, floating-point array to process.            !
!     SetFillVal Logical switch to set fill value in land areas        !
!                  (optional).                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts       Number of points processed in Awrk.                   !
!     Awrk       Collected data from each node packed into 1D array    !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(out) :: Npts
!
      real(dp), intent(in) :: Ascl
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: A(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
      logical :: LandFill
      integer :: Cgrid, ghost, rank
      integer :: Io, Ie, Jo, Je, Ioff, Joff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen
      integer :: Lstr, MyError, MyType, Serror, Srequest
      integer :: i, ic, j, jc, np
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: MySize
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus
!
      real(r8), dimension(TileSize(ng)) :: Asend
      real(r8), dimension(TileSize(ng),                                 &
     &                    NtileI(ng)*NtileJ(ng)-1) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_gather2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 66, 3834, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Asend)+                   &
     &                                   SIZE(Arecv))*KIND(A),r8))
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Npts=Ilen*Jlen
!
!  Set physical, non-overlapping (no ghost-points) ranges according to
!  tile rank.
!
      ghost=0
!
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT                              ! RHO-points
          Cgrid=2
      END SELECT
      Imin=BOUNDS(ng) % Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng) % Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,MyRank)
!
!  Compute size of distributed buffers.
!
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        MySize(rank)=(BOUNDS(ng) % Imax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Imin(Cgrid,ghost,rank)+1)*           &
     &               (BOUNDS(ng) % Jmax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Jmin(Cgrid,ghost,rank)+1)
      END DO
!
!  Initialize local arrays to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Collect requested array data.
!-----------------------------------------------------------------------
!
!  Pack and scale input data.
!
      np=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          np=np+1
          Asend(np)=A(i,j)*Ascl
        END DO
      END DO
!
!  If overwriting Land/Sea mask or processing water-points only, flag
!  land-points with special value.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
      IF (gtype.lt.0) THEN
        np=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            IF (Amask(i,j).eq.0.0_r8) THEN
              Asend(np)=spval
            END IF
          END DO
        END DO
      ELSE IF (LandFill) THEN
        np=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            IF (Amask(i,j).eq.0.0_r8) THEN
              Asend(np)=spval
            END IF
          END DO
        END DO
      END IF
!
!  If master processor, unpack the send buffer since there is not
!  need to distribute.
!
      IF (MyRank.eq.MyMaster) THEN
        np=0
        DO j=Jmin,Jmax
          jc=(j-Joff)*Ilen
          DO i=Imin,Imax
            np=np+1
            ic=i+Ioff+jc
            Awrk(ic)=Asend(np)
          END DO
        END DO
      END IF
!
!  Send, receive, and unpack data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), MySize(rank), MP_FLOAT, rank,  &
     &                    rank+5, OCN_COMM_WORLD, Rrequest(rank),       &
     &                    MyError)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_GATHER2D - error during ',a,' call, Node = ',&
     &              i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF
          np=0
          Imin=BOUNDS(ng) % Imin(Cgrid,ghost,rank)
          Imax=BOUNDS(ng) % Imax(Cgrid,ghost,rank)
          Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,rank)
          Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,rank)
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen
            DO i=Imin,Imax
              np=np+1
              ic=i+Ioff+jc
              Awrk(ic)=Arecv(np,rank)
            END DO
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, MySize(MyRank), MP_FLOAT, MyMaster,      &
     &                  MyRank+5, OCN_COMM_WORLD, Srequest, MyError)
        CALL mpi_wait (Srequest, Sstatus, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
! If pocessing only water-points, remove land points and repack.
!
      IF ((MyRank.eq.MyMaster).and.(gtype.lt.0)) THEN
        ic=0
        np=Ilen*Jlen
        DO i=1,np
          IF (Awrk(i).lt.spval) THEN
            ic=ic+1
            Awrk(ic)=Awrk(i)
          END IF
        END DO
        Npts=ic
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 66, 4065, MyFile)
!
      RETURN
      END SUBROUTINE mp_gather2d
!
      SUBROUTINE mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,  &
     &                        tindex, gtype, Ascl,                      &
     &                        Amask,                                    &
     &                        A, Npts, Awrk, SetFillVal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 3D tiled, floating-point array from each    !
!  spawned node and stores it into one dimensional global array. It    !
!  is used to collect and  pack output data.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     tindex     Time record index to process.                         !
!     gtype      C-grid type. If negative and Land-Sea is available,   !
!                  only water-points processed.                        !
!     Ascl       Factor to scale field before writing.                 !
!     Amask      Land/Sea mask, if any.                                !
!     A          3D tiled, floating-point array to process.            !
!     SetFillVal Logical switch to set fill value in land areas        !
!                  (optional).                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts       Number of points processed in Awrk.                   !
!     Awrk       Collected data from each node packed into 1D array    !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(out) :: Npts
!
      real(dp), intent(in) :: Ascl
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
      logical :: LandFill
!
      integer :: Cgrid, ghost, rank
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: Lstr, MyError, MyType, Serror, Srequest
      integer :: i, ic, j, jc, k, kc, np
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: MySize
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus
!
      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)) :: Asend
      real(r8), dimension(TileSize(ng)*(UBk-LBk+1),                     &
     &                    NtileI(ng)*NtileJ(ng)-1) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_gather3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 66, 4164, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Asend)+                   &
     &                                   SIZE(Arecv))*KIND(A),r8))
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
      Npts=IJlen*Klen
!
!  Set tile physical, non-overlapping (no ghost-points) ranges according
!  to tile rank.
!
      ghost=0
!
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT                              ! RHO-points
          Cgrid=2
      END SELECT
      Imin=BOUNDS(ng) % Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng) % Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,MyRank)
!
!  Compute size of distributed buffers.
!
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        MySize(rank)=(BOUNDS(ng) % Imax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Imin(Cgrid,ghost,rank)+1)*           &
     &               (BOUNDS(ng) % Jmax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Jmin(Cgrid,ghost,rank)+1)*           &
     &               (UBk-LBk+1)
      END DO
!
!  Initialize local arrays to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Collect requested array data.
!-----------------------------------------------------------------------
!
!  Pack and scale input data.
!
      np=0
      DO k=LBk,UBk
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            Asend(np)=A(i,j,k)*Ascl
          END DO
        END DO
      END DO
!
!  If overwriting Land/Sea mask or processing water-points only, flag
!  land-points with special value.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
      IF (gtype.lt.0) THEN
        np=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              np=np+1
              IF (Amask(i,j).eq.0.0_r8) THEN
                Asend(np)=spval
              END IF
            END DO
          END DO
        END DO
      ELSE IF (LandFill) THEN
        np=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              np=np+1
              IF (Amask(i,j).eq.0.0_r8) THEN
                Asend(np)=spval
              END IF
            END DO
          END DO
        END DO
      END IF
!
!  If master processor, unpack the send buffer since there is not
!  need to distribute.
!
      IF (MyRank.eq.MyMaster) THEN
        np=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen+kc
            DO i=Imin,Imax
              np=np+1
              ic=i+Ioff+jc
              Awrk(ic)=Asend(np)
            END DO
          END DO
        END DO
      END IF
!
!  Send, receive, and unpack data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), MySize(rank), MP_FLOAT, rank,  &
     &                    rank+5, OCN_COMM_WORLD, Rrequest(rank),       &
     &                    MyError)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_GATHER3D - error during ',a,' call, Node = ',&
     &              i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF
          np=0
          Imin=BOUNDS(ng) % Imin(Cgrid,ghost,rank)
          Imax=BOUNDS(ng) % Imax(Cgrid,ghost,rank)
          Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,rank)
          Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,rank)
          DO k=LBk,UBk
            kc=(k-Koff)*IJlen
            DO j=Jmin,Jmax
              jc=(j-Joff)*Ilen+kc
              DO i=Imin,Imax
                np=np+1
                ic=i+Ioff+jc
                Awrk(ic)=Arecv(np,rank)
              END DO
            END DO
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, MySize(MyRank), MP_FLOAT, MyMaster,      &
     &                  MyRank+5, OCN_COMM_WORLD, Srequest, MyError)
        CALL mpi_wait (Srequest, Sstatus, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
! If pocessing only water-points, remove land points and repack.
!
      IF ((MyRank.eq.MyMaster).and.(gtype.lt.0)) THEN
        ic=0
        np=IJlen*Klen
        DO i=1,np
          IF (Awrk(i).lt.spval) THEN
            ic=ic+1
            Awrk(ic)=Awrk(i)
          END IF
        END DO
        Npts=ic
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 66, 4416, MyFile)
!
      RETURN
      END SUBROUTINE mp_gather3d
!
      SUBROUTINE mp_gather_state (ng, model, Mstr, Mend, Asize,         &
     &                            A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine gathers (threaded to global) state data to all nodes   !
!  in the group. This  routine  is used to unpack the state data for   !
!  the GST analysis propagators.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Mstr       Threaded array lower bound.                           !
!     Mend       Threaded array upper bound.                           !
!     Asize      Size of the full state.                               !
!     A          Threaded 1D array process.                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       Collected data from each node packed into 1D full     !
!                  state array.                                        !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend, Asize
!
      real(r8), intent(in)  :: A(Mstr:Mend)
      real(r8), intent(out) :: Awrk(Asize)
!
!  Local variable declarations.
!
      integer :: LB, Lstr, MyError, Serror
      integer :: i, np, rank, request
      integer :: my_bounds(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(2,0:NtileI(ng)*NtileJ(ng)-1) :: Abounds
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_gather_state"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 66, 4475, MyFile)
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Collect data lower and upper bound dimensions.
!
      np=2
      my_bounds(1)=Mstr
      my_bounds(2)=Mend
      CALL mpi_allgather (my_bounds, np, MPI_INTEGER, Abounds, np,      &
     &                    MPI_INTEGER, OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 10     FORMAT (/,' MP_GATHER_STATE - error during ',a,                 &
     &          ' call, Node = ',i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, loop over other nodes and receive the data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          np=Abounds(2,rank)-Abounds(1,rank)+1
          LB=Abounds(1,rank)
          CALL mpi_irecv (Awrk(LB:), np, MP_FLOAT, rank, rank+5,        &
     &                    OCN_COMM_WORLD, request, MyError)
          CALL mpi_wait (request, status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
        END DO
!
!  Load master node contribution.
!
        DO i=Mstr,Mend
          Awrk(i)=A(i)
        END DO
!
!  Otherwise, send data to master node.
!
      ELSE
        np=Mend-Mstr+1
        CALL mpi_isend (A(Mstr:), np, MP_FLOAT, MyMaster, MyRank+5,     &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast collected data to all nodes.
!
      CALL mpi_bcast (Awrk, Asize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,  &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 66, 4558, MyFile)
!
      RETURN
      END SUBROUTINE mp_gather_state
!
      FUNCTION mp_ncread1d (ng, model, ncid, ncvname, ncname,           &
     &                      ncrec, LB1, UB1, Ascale, A)                 &
     &                     RESULT (io_error)
!
!***********************************************************************
!                                                                      !
!  This function reads floating point 1D state array from specified    !
!  NetCDF file and scatters it to the other nodes.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvname      NetCDF variable name.                               !
!     ncname       NetCDF file name.                                   !
!     ncrec        NetCDF record index to write. If negative, it       !
!                    assumes that the variable is recordless.          !
!     LB1          First-dimension Lower bound.                        !
!     UB1          First-dimension Upper bound.                        !
!     Ascale       Factor to scale field after reading (real).         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            Field to read in (real).                            !
!     io_error     Error flag (integer).                               !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1
!
      real(r8), intent(in) :: Ascale
      real(r8), intent(out) :: A(LB1:UB1)
!
      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: io_error
      integer :: ibuffer(2), my_bounds(2), start(1), total(1)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize
!
      real(r8), allocatable :: Asend(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_ncread1d_nf90"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 67, 4634, MyFile)
!
!-----------------------------------------------------------------------
!  Read requested NetCDF file and scatter it to all nodes.
!-----------------------------------------------------------------------
!
      io_error=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      np=2
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      CALL mpi_allgather (my_bounds, np, MPI_INTEGER,                   &
     &                    Asize, np, MPI_INTEGER,                       &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If not master node, receive data from master node.
!
      IF (MyRank.ne.MyMaster) THEN
        np=UB1-LB1+1
        CALL mpi_irecv (A(LB1:), np, MP_FLOAT, MyMaster, MyRank+5,      &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_IRECV', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
!
!  Scale recieved (read) data.
!
        DO i=LB1,UB1
          A(i)=A(i)*Ascale
        END DO
!
!  Otherwise, if master node allocate the send buffer.
!
      ELSE
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=Asize(2,rank)-Asize(1,rank)+1
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Asend)) THEN
          allocate (Asend(Npts))
        END IF
!
!  If master node, loop over all nodes and read buffers to send.
!
        io_error=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (io_error.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
          exit_flag=2
          ioerror=io_error
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO rank=0,NtileI(ng)*NtileJ(ng)-1
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            io_error=nf90_get_var(ncid, varid,  Asend, start, total)
            IF (io_error.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
              exit_flag=2
              ioerror=io_error
              EXIT
            END IF
!
!  Send buffer to all nodes, except itself.
!
            IF (rank.eq.MyMaster) THEN
              np=0
              DO i=LB1,UB1
                np=np+1
                A(i)=Asend(np)*Ascale
              END DO
            ELSE
              np=Asize(2,rank)-Asize(1,rank)+1
              CALL mpi_isend (Asend, np, MP_FLOAT, rank, rank+5,        &
     &                        OCN_COMM_WORLD, request, MyError)
              CALL mpi_wait (request, status, MyError)
              IF (MyError.ne.MPI_SUCCESS) THEN
                CALL mpi_error_string (MyError, string, Lstr, Serror)
                Lstr=LEN_TRIM(string)
                WRITE (stdout,10) 'MPI_ISEND', rank, MyError,           &
     &                            string(1:Lstr)
                exit_flag=2
                RETURN
              END IF
            END IF
          END DO
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Asend)*KIND(A),r8))
!
!  Deallocate send buffer.
!
      IF (allocated(Asend).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Asend)
      END IF
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 67, 4763, MyFile)
!
 10   FORMAT (/,' MP_NCREAD1D - error during ',a,' call, Node = ',i0,   &
     &        ' Error = ',i0,/,15x,a)
 20   FORMAT (/,' MP_NCREAD1D - error while inquiring ID for',          &
     &        ' variable: ',a,/,15x,'in file: ',a)
 30   FORMAT (/,' MP_NCREAD1D - error while reading variable: ',a,      &
     &        /,15x,'in file: ',a)
!
      RETURN
      END FUNCTION mp_ncread1d
!
      FUNCTION mp_ncread2d (ng, model, ncid, ncvname, ncname,           &
     &                      ncrec, LB1, UB1, LB2, UB2, Ascale, A)       &
     &                     RESULT (io_error)
!
!***********************************************************************
!                                                                      !
!  This function reads floating point 2D state array from specified    !
!  NetCDF file and scatters it to the other nodes.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvname      NetCDF variable name.                               !
!     ncname       NetCDF file name.                                   !
!     ncrec        NetCDF record index to write. If negative, it       !
!                    assumes that the variable is recordless.          !
!     LB1          First-dimension Lower bound.                        !
!     UB1          First-dimension Upper bound.                        !
!     LB2          Second-dimension Lower bound.                       !
!     UB2          Second-dimension Upper bound.                       !
!     Ascale       Factor to scale field after reading (real).         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            Field to read in (real).                            !
!     io_error     Error flag (integer).                               !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1, LB2, UB2
!
      real(r8), intent(in) :: Ascale
      real(r8), intent(out) :: A(LB1:UB1,LB2:UB2)
!
      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: io_error
      integer :: ibuffer(2), my_bounds(4), start(2), total(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize
!
      real(r8), allocatable :: Asend(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_ncread2d_nf90"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 67, 4848, MyFile)
!
!-----------------------------------------------------------------------
!  Read requested NetCDF file and scatter it to all nodes.
!-----------------------------------------------------------------------
!
      io_error=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      np=4
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      my_bounds(3)=LB2
      my_bounds(4)=UB2
      CALL mpi_allgather (my_bounds, np, MPI_INTEGER,                   &
     &                    Asize, np, MPI_INTEGER,                       &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If not master node, receive data from master node.
!
      IF (MyRank.ne.MyMaster) THEN
        np=(UB1-LB1+1)*(UB2-LB2+1)
        CALL mpi_irecv (A(LB1,LB2), np, MP_FLOAT, MyMaster, MyRank+5,   &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_IRECV', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
!
!  Scale recieved (read) data.
!
        DO j=LB2,UB2
          DO i=LB1,UB1
            A(i,j)=A(i,j)*Ascale
          END DO
        END DO
!
!  Otherwise, if master node allocate the send buffer.
!
      ELSE
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)*                           &
     &       (Asize(4,rank)-Asize(3,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Asend)) THEN
          allocate (Asend(Npts))
        END IF
!
!  If master node, loop over all nodes and read buffers to send.
!
        io_error=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (io_error.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
          exit_flag=2
          ioerror=io_error
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO rank=0,NtileI(ng)*NtileJ(ng)-1
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            start(2)=Asize(3,rank)
            total(2)=Asize(4,rank)-Asize(3,rank)+1
            io_error=nf90_get_var(ncid, varid,  Asend, start, total)
            IF (io_error.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
              exit_flag=2
              ioerror=io_error
              EXIT
            END IF
!
!  Send buffer to all nodes, except itself.
!
            IF (rank.eq.MyMaster) THEN
              np=0
              DO j=LB2,UB2
                DO i=LB1,UB1
                  np=np+1
                  A(i,j)=Asend(np)*Ascale
                END DO
              END DO
            ELSE
              np=(Asize(2,rank)-Asize(1,rank)+1)*                       &
     &           (Asize(4,rank)-Asize(3,rank)+1)
              CALL mpi_isend (Asend, np, MP_FLOAT, rank, rank+5,        &
     &                        OCN_COMM_WORLD, request, MyError)
              CALL mpi_wait (request, status, MyError)
              IF (MyError.ne.MPI_SUCCESS) THEN
                CALL mpi_error_string (MyError, string, Lstr, Serror)
                Lstr=LEN_TRIM(string)
                WRITE (stdout,10) 'MPI_ISEND', rank, MyError,           &
     &                            string(1:Lstr)
                exit_flag=2
                RETURN
              END IF
            END IF
          END DO
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Asend)*KIND(A),r8))
!
!  Deallocate send buffer.
!
      IF (allocated(Asend).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Asend)
      END IF
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 67, 4987, MyFile)
!
 10   FORMAT (/,' MP_NCREAD2D - error during ',a,' call, Node = ',i0,   &
     &        ' Error = ',i0,/,15x,a)
 20   FORMAT (/,' MP_NCREAD2D - error while inquiring ID for',          &
     &        ' variable: ',a,/,15x,'in file: ',a)
 30   FORMAT (/,' MP_NCREAD2D - error while reading variable: ',a,      &
     &        /,15x,'in file: ',a)
!
      RETURN
      END FUNCTION mp_ncread2d
!
      FUNCTION mp_ncwrite1d (ng, model, ncid, ncvname, ncname,          &
     &                       ncrec, LB1, UB1, Ascale, A)                &
     &                      RESULT (io_error)
!
!***********************************************************************
!                                                                      !
!  This function collects floating point 1D state array data from the  !
!  other nodes and writes it into specified NetCDF file.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number.                                !
!     model         Calling model identifier.                          !
!     ncid          NetCDF file ID.                                    !
!     ncvname       NetCDF variable name.                              !
!     ncname        NetCDF file name.                                  !
!     ncrec         NetCDF record index to write. If negative, it      !
!                     assumes that the variable is recordless.         !
!     LB1           First-dimension Lower bound.                       !
!     UB1           First-dimension Upper bound.                       !
!     Ascale        Factor to scale field before writing (real).       !
!     A             Field to write out (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     io_error      Error flag (integer).                              !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1
!
      real(r8), intent(in) :: Ascale
      real(r8), intent(in) :: A(LB1:UB1)
      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: io_error
      integer :: ibuffer(2), my_bounds(2), start(1), total(1)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(2,0:NtileI(ng)*NtileJ(ng)-1) :: Asize
!
      real(r8), allocatable :: Arecv(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_ncwrite1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 66, 5070, MyFile)
!
!-----------------------------------------------------------------------
!  Collect and write data into requested NetCDF file.
!-----------------------------------------------------------------------
!
      io_error=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      np=2
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      CALL mpi_allgather (my_bounds, np, MPI_INTEGER,                   &
     &                    Asize, np, MPI_INTEGER,                       &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, allocate the receive buffer.
!
      IF (MyRank.eq.MyMaster) THEN
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  Write out master node contribution.
!
        start(1)=LB1
        total(1)=UB1-LB1+1
        np=0
        DO i=LB1,UB1
          np=np+1
          Arecv(np)=A(i)
        END DO
        io_error=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (io_error.eq.nf90_noerr) THEN
          io_error=nf90_put_var(ncid, varid, Arecv, start, total)
          IF (io_error.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
            exit_flag=3
            ioerror=io_error
          END IF
        ELSE
          WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
          exit_flag=3
          ioerror=io_error
        END IF
!
!  If master node, loop over other nodes and receive the data.
!
        IF (exit_flag.eq.NoError) THEN
          DO rank=1,NtileI(ng)*NtileJ(ng)-1
            np=Asize(2,rank)-Asize(1,rank)+1
            CALL mpi_irecv (Arecv, np, MP_FLOAT, rank, rank+5,          &
     &                      OCN_COMM_WORLD, request, MyError)
            CALL mpi_wait (request, status, MyError)
            IF (MyError.ne.MPI_SUCCESS) THEN
              CALL mpi_error_string (MyError, string, Lstr, Serror)
              Lstr=LEN_TRIM(string)
              WRITE (stdout,10) 'MPI_IRECV', rank, MyError,             &
     &                          string(1:Lstr)
              exit_flag=3
              RETURN
            END IF
!
!  Write out data into NetCDF file.
!
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            DO i=1,np
              Arecv(i)=Arecv(i)*Ascale
            END DO
            io_error=nf90_put_var(ncid, varid, Arecv, start, total)
            IF (io_error.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
              exit_flag=3
              ioerror=io_error
              EXIT
            END IF
          END DO
        END IF
!
!  Otherwise, send data to master node.
!
      ELSE
        np=UB1-LB1+1
        CALL mpi_isend (A(LB1:), np, MP_FLOAT, MyMaster, MyRank+5,      &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Arecv)*KIND(A),r8))
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 66, 5205, MyFile)
!
 10   FORMAT (/,' MP_NCWRITE1D - error during ',a,' call, Node = ',i0,  &
     &        ' Error = ',i0,/,21x,a)
 20   FORMAT (/,' MP_NCWRITE1D - error while writing variable: ',a,     &
     &        /,16x,'in file: ',a)
 30   FORMAT (/,' MP_NCWRITE1D - error while inquiring ID for',         &
     &        ' variable: ',a,/,16x,'in file: ',a)
!
      RETURN
      END FUNCTION mp_ncwrite1d
!
      FUNCTION mp_ncwrite2d (ng, model, ncid, ncvname, ncname,          &
     &                       ncrec, LB1, UB1, LB2, UB2, Ascale, A)      &
     &                      RESULT (io_error)
!
!***********************************************************************
!                                                                      !
!  This function collects floating point 2D state array data from the  !
!  other nodes and writes it into specified NetCDF file.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number.                                !
!     model         Calling model identifier.                          !
!     ncid          NetCDF file ID.                                    !
!     ncvname       NetCDF variable name.                              !
!     ncname        NetCDF file name.                                  !
!     ncrec         NetCDF record index to write. If negative, it      !
!                     assumes that the variable is recordless.         !
!     LB1           First-dimension Lower bound.                       !
!     UB1           First-dimension Upper bound.                       !
!     LB2           Second-dimension Lower bound.                      !
!     UB2           Second-dimension Upper bound.                      !
!     Ascale        Factor to scale field before writing (real).       !
!     A             Field to write out (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     io_error      Error flag (integer).                              !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1, LB2, UB2
!
      real(r8), intent(in) :: Ascale
      real(r8), intent(in) :: A(LB1:UB1,LB2:UB2)
!
      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: io_error
      integer :: ibuffer(2), my_bounds(4), start(2), total(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize
!
      real(r8), allocatable :: Arecv(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_ncwrite2d_nf90"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 66, 5290, MyFile)
!
!-----------------------------------------------------------------------
!  Collect and write data into requested NetCDF file.
!-----------------------------------------------------------------------
!
      io_error=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      np=4
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      my_bounds(3)=LB2
      my_bounds(4)=UB2
      CALL mpi_allgather (my_bounds, np, MPI_INTEGER,                   &
     &                    Asize, np, MPI_INTEGER,                       &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, allocate the receive buffer.
!
      IF (MyRank.eq.MyMaster) THEN
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)*                           &
     &       (Asize(4,rank)-Asize(3,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  Write out master node contribution.
!
        start(1)=LB1
        total(1)=UB1-LB1+1
        start(2)=LB2
        total(2)=UB2-LB2+1
        np=0
        DO j=LB2,UB2
          DO i=LB1,UB1
            np=np+1
            Arecv(np)=A(i,j)
          END DO
        END DO
        io_error=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (io_error.eq.nf90_noerr) THEN
          io_error=nf90_put_var(ncid, varid, Arecv, start, total)
          IF (io_error.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
            exit_flag=3
            ioerror=io_error
          END IF
        ELSE
          WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
          exit_flag=3
          ioerror=io_error
        END IF
!
!  If master node, loop over other nodes and receive the data.
!
        IF (exit_flag.eq.NoError) THEN
          DO rank=1,NtileI(ng)*NtileJ(ng)-1
            np=(Asize(2,rank)-Asize(1,rank)+1)*                         &
     &         (Asize(4,rank)-Asize(3,rank)+1)
            CALL mpi_irecv (Arecv, np, MP_FLOAT, rank, rank+5,          &
     &                      OCN_COMM_WORLD, request, MyError)
            CALL mpi_wait (request, status, MyError)
            IF (MyError.ne.MPI_SUCCESS) THEN
              CALL mpi_error_string (MyError, string, Lstr, Serror)
              Lstr=LEN_TRIM(string)
              WRITE (stdout,10) 'MPI_IRECV', rank, MyError,             &
     &                          string(1:Lstr)
              exit_flag=3
              RETURN
            END IF
!
!  Write out data into NetCDF file.
!
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            start(2)=Asize(3,rank)
            total(2)=Asize(4,rank)-Asize(3,rank)+1
            DO i=1,np
              Arecv(i)=Arecv(i)*Ascale
            END DO
            io_error=nf90_put_var(ncid, varid, Arecv, start, total)
            IF (io_error.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
              exit_flag=3
              ioerror=io_error
              EXIT
            END IF
          END DO
        END IF
!
!  Otherwise, send data to master node.
!
      ELSE
        np=(UB1-LB1+1)*(UB2-LB2+1)
        CALL mpi_isend (A(LB1:,LB2:), np, MP_FLOAT, MyMaster, MyRank+5, &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Arecv)*KIND(A),r8))
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 66, 5435, MyFile)
!
 10   FORMAT (/,' MP_NCWRITE2D - error during ',a,' call, Node = ',i0,  &
     &        ' Error = ',i0,/,21x,a)
 20   FORMAT (/,' MP_NCWRITE2D - error while writing variable: ',a,     &
     &        /,16x,'in file: ',a)
 30   FORMAT (/,' MP_NCWRITE2D - error while inquiring ID for',         &
     &        ' variable: ',a,/,16x,'in file: ',a)
!
      RETURN
      END FUNCTION mp_ncwrite2d
!
      SUBROUTINE mp_reduce_i8 (ng, model, Asize, A, handle_op, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine collects and reduces requested 1D integer array        !
!  variables from all nodes in the group.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     handle_op  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Asize
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(in) :: handle_op(Asize)
!
      integer(i8b), intent(inout) :: A(Asize)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Serror
      integer :: handle, i, rank, request
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus
!
      integer(i8b), dimension(Asize,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
      integer(i8b), dimension(Asize) :: Areduce
      integer(i8b), dimension(Asize) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_reduce_1di"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 65, 5507, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Arecv)+                   &
     &                                   2*Asize)*KIND(A),r8))
!
!  Pack data to reduce.
!
      DO i=1,Asize
        Asend(i)=A(i)
      END DO
!
!  Collect and reduce.
!
      DO i=1,Asize
        IF (handle_op(i)(1:3).eq.'MIN') THEN
          handle=MPI_MIN
        ELSE IF (handle_op(i)(1:3).eq.'MAX') THEN
          handle=MPI_MAX
        ELSE IF (handle_op(i)(1:3).eq.'SUM') THEN
          handle=MPI_SUM
        END IF
        CALL mpi_allreduce (Asend(i), Areduce(i), 1, MPI_INTEGER,       &
     &                      handle, MyCOMM, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,           &
     &                      string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END DO
 10   FORMAT (/,' MP_REDUCE_I8 - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,16x,a)
!
!  Unpack.
!
      DO i=1,Asize
        A(i)=Areduce(i)
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 65, 5652, MyFile)
!
      RETURN
      END SUBROUTINE mp_reduce_i8
!
      SUBROUTINE mp_reduce_0d (ng, model, Asize, A, handle_op, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine collects and reduces requested variables from all      !
!  nodes in the group.  Then,  it broadcasts reduced variables to      !
!  all nodes in the group.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     handle_op  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Asize
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(in) :: handle_op
!
      real(r8), intent(inout) :: A
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Npts, Serror
      integer :: handle, i, rank, request
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus
!
      real(r8) :: Areduce
      real(r8) :: Asend
      real(r8), dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_reduce_0d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 65, 6127, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Pack data to reduce.
!
      Asend=A
      Npts=1
!
!  Collect and reduce.
!
      IF (handle_op(1:3).eq.'MIN') THEN
        handle=MPI_MIN
      ELSE IF (handle_op(1:3).eq.'MAX') THEN
        handle=MPI_MAX
      ELSE IF (handle_op(1:3).eq.'SUM') THEN
        handle=MPI_SUM
      END IF
      CALL mpi_allreduce (Asend, Areduce, Npts, MP_FLOAT, handle,       &
     &                    MyCOMM, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
 10   FORMAT (/,' MP_REDUCE_0D - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,16x,a)
!
!  Unpack.
!
      A=Areduce
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 65, 6256, MyFile)
!
      RETURN
      END SUBROUTINE mp_reduce_0d
!
      SUBROUTINE mp_reduce_1d (ng, model, Asize, A, handle_op, InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine collects and reduces requested variables from all      !
!  nodes in the group.  Then,  it broadcasts reduced variables to      !
!  all nodes in the group.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     handle_op  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Asize
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(in) :: handle_op(Asize)
!
      real(r8), intent(inout) :: A(Asize)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Serror
      integer :: handle, i, rank, request
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus
!
      real(r8), dimension(Asize,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
      real(r8), dimension(Asize) :: Areduce
      real(r8), dimension(Asize) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_reduce_1d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 65, 6322, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Arecv)+                   &
     &                                   2*Asize)*KIND(A),r8))
!
!  Pack data to reduce.
!
      DO i=1,Asize
        Asend(i)=A(i)
      END DO
!
!  Collect and reduce.
!
      DO i=1,Asize
        IF (handle_op(i)(1:3).eq.'MIN') THEN
          handle=MPI_MIN
        ELSE IF (handle_op(i)(1:3).eq.'MAX') THEN
          handle=MPI_MAX
        ELSE IF (handle_op(i)(1:3).eq.'SUM') THEN
          handle=MPI_SUM
        END IF
        CALL mpi_allreduce (Asend(i), Areduce(i), 1, MP_FLOAT, handle,  &
     &                      MyCOMM, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,           &
     &                      string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END DO
 10   FORMAT (/,' MP_REDUCE_1D - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,16x,a)
!
!  Unpack.
!
      DO i=1,Asize
        A(i)=Areduce(i)
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 65, 6466, MyFile)
!
      RETURN
      END SUBROUTINE mp_reduce_1d
!
      SUBROUTINE mp_reduce2 (ng, model, Isize, Jsize, A, handle_op,     &
     &                       InpComm)
!
!***********************************************************************
!                                                                      !
!  This routine computes the global minimum/maximum and its associated !
!  qualifiers like: location and/or other scalar components.  Then, it !
!  it broadcasts reduced variables to all nodes in the group.          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Isize      Size of I-dimension: the minimum/maximum to reduce    !
!                  is in location A(1,:) and qualifiers A(2:Isize,:).  !
!     Jsize      Size of J-dimension: number of different sets of      !
!                  minimum and/or maximum to process.                  !
!     A          Matrix of variables and qualifiers to reduce.         !
!     handle_op  Reduction operation handle (string) of size Jsize.    !
!                  The following  reduction operations are supported:  !
!                  'MINLOC', 'MAXLOC'                                  !
!     InpComm    Communicator handle (integer, OPTIONAL).              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Matrix of reduced variables and qualifiers.           !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Isize, Jsize
      integer, intent(in), optional :: InpComm
!
      character (len=*), intent(in) :: handle_op(Jsize)
!
      real(r8), intent(inout) :: A(Isize,Jsize)
!
!  Local variable declarations.
!
      integer :: Lstr, MyCOMM, MyError, Serror
      integer :: handle, i, j
!
      real(r8), dimension(2,Isize) :: Areduce
      real(r8), dimension(2,Isize) :: Asend
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_reduce2"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 65, 6530, MyFile)
!
!-----------------------------------------------------------------------
!  Set distributed-memory communicator handle (context ID).
!-----------------------------------------------------------------------
!
      IF (PRESENT(InpComm)) THEN
        MyCOMM=InpComm
      ELSE
        MyCOMM=OCN_COMM_WORLD
      END IF
!
!-----------------------------------------------------------------------
!  Reduce requested variables and qualifiers.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Areduce)+                 &
     &                                   SIZE(Asend))*KIND(A),r8))
!
!  Pack and reduce.
!
      DO j=1,Jsize
        DO i=1,Isize
          Asend(1,i)=A(1,j)
          Asend(2,i)=A(i,j)
        END DO
        IF (handle_op(j)(1:6).eq.'MINLOC') THEN
          handle=MPI_MINLOC
        ELSE IF (handle_op(j)(1:6).eq.'MAXLOC') THEN
          handle=MPI_MAXLOC
        END IF
        CALL mpi_allreduce (Asend, Areduce, Isize,                      &
     &                      MPI_2DOUBLE_PRECISION,                      &
     &                      handle, MyCOMM, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,           &
     &                      string(1:Lstr)
 10       FORMAT (/,' MP_REDUCE2 - error during ',a,' call, Node = ',   &
     &            i3.3,' Error = ',i3,/,16x,a)
          exit_flag=2
          RETURN
        END IF
!
!  Unpack.
!
        A(1,j)=Areduce(1,1)
        DO i=2,Isize
          A(i,j)=Areduce(2,i)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 65, 6598, MyFile)
!
      RETURN
      END SUBROUTINE mp_reduce2
!
      SUBROUTINE mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,           &
     &                         Nghost, gtype, Amin, Amax,               &
     &                         Npts, A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts input global data, packed as 1D real array, !
!  to each spawned mpi node.  Because this routine is also used by the !
!  adjoint model,  the ghost-points in the halo region are NOT updated !
!  in the ouput tile array (Awrk).  It is used by the  master node  to !
!  scatter input global data to each tiled node.                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     gtype      C-grid type. If negative and Land-Sea mask is         !
!                  available, only water-points are processed.         !
!     Amin       Input array minimum value.                            !
!     Amax       Input array maximum value.                            !
!     NWpts      Number of water points.                               !
!     IJ_water   IJ-indices for water points.                          !
!     Npts       Number of points to processes in A.                   !
!     A          Input global data from each node packed into 1D array !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!     Npts       Number of points to processes in A.                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       2D tiled, floating-point array.                       !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost, gtype, Npts
!
      real(r8), intent(inout) :: Amin, Amax
      real(r8), intent(inout) :: A(Npts+2)
      real(r8), intent(out) :: Awrk(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: Io, Ie, Jo, Je, Ioff, Joff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen
      integer :: Lstr, MyError, MySize, MyType, Serror, ghost
      integer :: i, ic, ij, j, jc, mc, nc
!
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_scatter2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 67, 6683, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Arecv)*KIND(A),r8))
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      IJlen=Ilen*Jlen
!
!  Set physical, non-overlapping (Nghost=0) or overlapping (Nghost>0)
!  ranges according to tile rank.
!
      IF (Nghost.eq.0) THEN
        ghost=0                                   ! non-overlapping
      ELSE
        ghost=1                                   ! overlapping
      END IF
!
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Imin=BOUNDS(ng) % Imin(1,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(1,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(1,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(1,ghost,MyRank)
        CASE (r2dvar, r3dvar)
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
        CASE (u2dvar, u3dvar)
          Imin=BOUNDS(ng) % Imin(3,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(3,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(3,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(3,ghost,MyRank)
        CASE (v2dvar, v3dvar)
          Imin=BOUNDS(ng) % Imin(4,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(4,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(4,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(4,ghost,MyRank)
        CASE DEFAULT                              ! RHO-points
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
      END SELECT
!
!  Size of broadcast buffer.
!
      IF (gtype.gt.0) THEN
        MySize=IJlen
      ELSE
        MySize=Npts
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Scatter requested array data.
!-----------------------------------------------------------------------
!
!  If master processor, append minimum and maximum values to the end of
!  the buffer.
!
      IF (MyRank.eq.MyMaster) Then
        A(MySize+1)=Amin
        A(MySize+2)=Amax
      END IF
      MySize=MySize+2
!
!  Broadcast data to all processes in the group, itself included.
!
      CALL mpi_bcast (A, MySize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
         Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_SCATTER2D - error during ',a,' call, Node = ',   &
     &          i3.3, ' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If water points only, fill land points.
!
      IF (gtype.gt.0) THEN
        DO nc=1,MySize-2
          Arecv(nc)=A(nc)
        END DO
      END IF
!
!  Unpack data buffer.
!
      DO j=Jmin,Jmax
        jc=(j-Joff)*Ilen
        DO i=Imin,Imax
          ic=i+Ioff+jc
          Awrk(i,j)=Arecv(ic)
        END DO
      END DO
      Amin=A(MySize-1)
      Amax=A(MySize)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 67, 6865, MyFile)
!
      RETURN
      END SUBROUTINE mp_scatter2d
!
      SUBROUTINE mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk, &
     &                         Nghost, gtype, Amin, Amax,               &
     &                         Npts, A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts input global data, packed as 1D real array, !
!  to each spawned mpi node.  Because this routine is also used by the !
!  adjoint model,  the ghost-points in the halo region are NOT updated !
!  in the ouput tile array (Awrk).  It is used by the  master node  to !
!  scatter input global data to each tiled node.                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     gtype      C-grid type. If negative and Land-Sea mask is         !
!                  available, only water-points are processed.         !
!     Amin       Input array minimum value.                            !
!     Amax       Input array maximum value.                            !
!     NWpts      Number of water points.                               !
!     IJ_water   IJ-indices for water points.                          !
!     Npts       Number of points to processes in A.                   !
!     A          Input global data from each node packed into 1D array !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!     Npts       Number of points to processes in A.                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       3D tiled, floating-point array.                       !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost, gtype, Npts
!
      real(r8), intent(inout) :: Amin, Amax
      real(r8), intent(inout) :: A(Npts+2)
      real(r8), intent(out) :: Awrk(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: Lstr, MyError, MySize, MyType, Serror, ghost
      integer :: i, ic, ij, j, jc, k, kc, mc, nc
!
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_scatter3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 67, 6952, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(SIZE(Arecv)*KIND(A),r8))
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
!
!  Set physical, non-overlapping (Nghost=0) or overlapping (Nghost>0)
!  ranges according to tile rank.
!
      IF (Nghost.eq.0) THEN
        ghost=0                                   ! non-overlapping
      ELSE
        ghost=1                                   ! overlapping
      END IF
!
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Imin=BOUNDS(ng) % Imin(1,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(1,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(1,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(1,ghost,MyRank)
        CASE (r2dvar, r3dvar)
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
        CASE (u2dvar, u3dvar)
          Imin=BOUNDS(ng) % Imin(3,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(3,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(3,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(3,ghost,MyRank)
        CASE (v2dvar, v3dvar)
          Imin=BOUNDS(ng) % Imin(4,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(4,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(4,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(4,ghost,MyRank)
        CASE DEFAULT                              ! RHO-points
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
      END SELECT
!
!  Size of broadcast buffer.
!
      IF (gtype.gt.0) THEN
        MySize=IJlen*Klen
      ELSE
        MySize=Npts
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Scatter requested array data.
!-----------------------------------------------------------------------
!
!  If master processor, append minimum and maximum values to the end of
!  the buffer.
!
      IF (MyRank.eq.MyMaster) Then
        A(MySize+1)=Amin
        A(MySize+2)=Amax
      END IF
      MySize=MySize+2
!
!  Broadcast data to all processes in the group, itself included.
!
      CALL mpi_bcast (A, MySize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
         Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_SCATTER3D - error during ',a,' call, Node = ',   &
     &          i3.3, ' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If water points only, fill land points.
!
      IF (gtype.gt.0) THEN
        DO nc=1,MySize-2
          Arecv(nc)=A(nc)
        END DO
      END IF
!
!  Unpack data buffer.
!
      DO k=LBk,UBk
        kc=(k-Koff)*IJlen
        DO j=Jmin,Jmax
          jc=(j-Joff)*Ilen+kc
          DO i=Imin,Imax
            ic=i+Ioff+jc
            Awrk(i,j,k)=Arecv(ic)
          END DO
        END DO
      END DO
      Amin=A(MySize-1)
      Amax=A(MySize)
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 67, 7147, MyFile)
!
      RETURN
      END SUBROUTINE mp_scatter3d
!
      SUBROUTINE mp_scatter_state (ng, model, Mstr, Mend, Asize,        &
     &                             A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine scatters (global to threaded) state data to all nodes  !
!  in the group. Before this can be done, the global data needs to be  !
!  collected from all the  nodes  by the master.  This is achieved by  !
!  summing the input values at each point.  This  routine  is used to  !
!  pack the state data for the GST analysis propagators.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Mstr       Threaded array lower bound.                           !
!     Mend       Threaded array upper bound.                           !
!     Asize      Size of array A.                                      !
!     A          Threaded 1D array process.                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data from all nodes.                        !
!     Awrk       Threaded block of data.                               !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend, Asize
!
      real(r8), intent(inout)  :: A(Asize)
      real(r8), intent(out) :: Awrk(Mstr:Mend)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Serror
      integer :: i, rank, request
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      real(r8), allocatable :: Arecv(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_scatter_state"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 67, 7211, MyFile)
!
!-----------------------------------------------------------------------
!  Collect data blocks from all nodes and scatter the data to all nodes.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL(Asize*KIND(A),r8))
!
!  All nodes have distinct pieces of the data and zero everywhere else.
!  So the strategy here is for the master node to receive the data from
!  the other nodes (excluding itself) and accumulate the sum at each
!  point. Then, the master node broadcast (itself included) its copy of
!  the accumlated data to other the nodes in the group. After this, each
!  node loads only the required block of the data into output array.
!
!  Notice that only the master node allocates the recieving buffer
!  (Arecv). It also receives only buffer at the time to avoid having
!  a very large communication array.  So here memory is more important
!  than time.
!
      IF (MyRank.eq.MyMaster) THEN
!
!  If master node, allocate and receive buffer.
!
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Asize))
        END IF
!
!  If master node, loop over other nodes to receive and accumulate the
!  data.
!
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv, Asize, MP_FLOAT, rank, rank+5,         &
     &                    OCN_COMM_WORLD, Rrequest(rank), MyError)
          CALL mpi_wait (Rrequest(rank), status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_SCATTER_STATE - error during ',a,            &
     &              ' call, Node = ', i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Asize
            A(i)=A(i)+Arecv(i)
          END DO
        END DO
!
!  Otherwise, send data to master node.
!
      ELSE
        CALL mpi_isend (A, Asize, MP_FLOAT, MyMaster, MyRank+5,         &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast accumulated (full) data to all nodes.
!
      CALL mpi_bcast (A, Asize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,     &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load appropriate data block into output array.
!
      DO i=Mstr,Mend
        Awrk(i)=A(i)
      END DO
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 67, 7308, MyFile)
!
      RETURN
      END SUBROUTINE mp_scatter_state
!
      SUBROUTINE mp_dump (ng, tile, gtype,                              &
     &                    ILB, IUB, JLB, JUB, KLB, KUB, A, name)
!
!***********************************************************************
!                                                                      !
!  This routine is used to debug distributed-memory communications.    !
!  It writes field into an ASCII file for further post-processing.     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype
      integer, intent(in) :: ILB, IUB, JLB, JUB, KLB, KUB
!
      real(r8), intent(in) :: A(ILB:IUB,JLB:JUB,KLB:KUB)
!
      character (len=*) :: name
!
!  Local variable declarations.
!
      common /counter/ nc
      integer :: nc
!
      logical, save :: first = .TRUE.
!
      integer :: Imin, Imax, Ioff, Jmin, Jmax, Joff
      integer :: unit
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_dump"
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
!------------------------------------------------------------------------
!  Write out requested field.
!------------------------------------------------------------------------
!
      IF (first) THEN
        nc=0
        first=.FALSE.
      END IF
      nc=nc+1
      IF (Master) THEN
        WRITE (10,'(a,i3.3,a,a)') 'file ', nc, ': ', TRIM(name)
        CALL my_flush (10)
      END IF
!
!  Write out field including ghost-points.
!
      Imin=0
      Imax=Lm(ng)+1
      IF (EWperiodic(ng)) THEN
        Ioff=3
      ELSE
        Ioff=1
      END IF
      Jmin=0
      Jmax=Mm(ng)+1
      IF (NSperiodic(ng)) THEN
        Joff=3
      ELSE
        Joff=1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=1
      END IF
      unit=(MyRank+1)*1000+nc
      WRITE (unit,*) ILB, IUB, JLB, JUB, KLB, KUB,                      &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(ILB:IUB,JLB:JUB,KLB:KUB)
      CALL my_flush (unit)
!
!  Write out non-overlapping field.
!
      Imin=IstrR
      Imax=IendR
      IF (EWperiodic(ng)) THEN
        Ioff=2
      ELSE
        Ioff=1
      END IF
      Jmin=JstrR
      Jmax=JendR
      IF (NSperiodic(ng)) THEN
        Joff=2
      ELSE
        Joff=1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=Istr
        Ioff=Ioff-1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=Jstr
        Joff=Joff-1
      END IF
      unit=(MyRank+1)*10000+nc
      WRITE (unit,*) Imin, Imax, Jmin, Jmax, KLB, KUB,                  &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(Imin:Imax,Jmin:Jmax,KLB:KUB)
      CALL my_flush (unit)
      RETURN
      END SUBROUTINE mp_dump
!
      SUBROUTINE mp_aggregate2d (ng, model, gtype,                      &
     &                           LBiT, UBiT, LBjT, UBjT,                &
     &                           LBiG, UBiG, LBjG, UBjG,                &
     &                           Atiled, Aglobal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 2D tiled, floating-point array from each    !
!  spawned node and stores it into 2D global array. If nesting, the    !
!  global array contains the contact points data.                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     gtype      C-grid type.                                          !
!     LBiT       Tiled  array, I-dimension Lower bound.                !
!     UBiT       Tiled  array, I-dimension Upper bound.                !
!     LBjT       Tiled  array, J-dimension Lower bound.                !
!     UBjT       Tiled  array, J-dimension Upper bound.                !
!     LBiG       Global array, I-dimension Lower bound.                !
!     UBiG       Global array, I-dimension Upper bound.                !
!     LBjG       Global array, J-dimension Lower bound.                !
!     UBjG       Global array, J-dimension Upper bound.                !
!     Atiled     2D tiled, floating-point array to process.            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aglobal    2D global array, all tiles are aggregated.            !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype
      integer, intent(in) :: LBiT, UBiT, LBjT, UBjT
      integer, intent(in) :: LBiG, UBiG, LBjG, UBjG
!
      real(r8), intent(in)  :: Atiled(LBiT:UBiT,LBjT:UBjT)
      real(r8), intent(out) :: Aglobal(LBiG:UBiG,LBjG:UBjG)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyType, Nnodes, Npts, Serror
      integer :: i, j, np, rank
      integer,  dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: my_bounds
!
      real(r8), dimension(TileSize(ng)) :: Asend
      real(r8), dimension(TileSize(ng)*                                 &
     &                    NtileI(ng)*NtileJ(ng)) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_aggregate2d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 71, 7496, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Asend)+                   &
     &                                   SIZE(Aglobal)+                 &
     &                                   SIZE(Arecv))*KIND(Asend),r8))
!
!  Number of nodes in the group.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
!
!  Set starting and ending indices to process including contact points
!  (if nesting) according to the staggered C-grid classification.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
        CASE (r2dvar, r3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (u2dvar, u3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (v2dvar, v3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
      END SELECT
!
!  Determine the maximum number of points to process between all tiles.
!  In collective communications, the amount of data sent must be equal
!  to the amount of data received.
!
      Npts=0
      DO rank=0,Nnodes
        np=(my_bounds(2,rank)-my_bounds(1,rank)+1)*                     &
     &     (my_bounds(4,rank)-my_bounds(3,rank)+1)
        Npts=MAX(Npts, np)
      END DO
      IF (Npts.gt.TileSize(ng)) THEN
        IF (Master) THEN
          WRITE (stdout,10) ' TileSize = ', TileSize(ng), Npts
 10       FORMAT (/,' MP_AGGREGATE2D - communication buffer to small,', &
     &            a, 2i8)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Initialize local arrays to facilitate collective communicatios.
!  This also avoid denormalized values, which facilitates debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Pack tile data.
!-----------------------------------------------------------------------
!
      np=0
      DO j=my_bounds(3,MyRank),my_bounds(4,MyRank)
        DO i=my_bounds(1,MyRank),my_bounds(2,MyRank)
          np=np+1
          Asend(np)=Atiled(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Aggregate data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,20) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 20     FORMAT (/,' MP_AGGREGATE2D - error during ',a,' call, Node = ', &
     &          i3.3,' Error = ',i3,/,18x,a)
        exit_flag=5
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data into a global 2D array.
!-----------------------------------------------------------------------
!
      DO rank=0,Nnodes
        np=rank*Npts
        DO j=my_bounds(3,rank),my_bounds(4,rank)
          DO i=my_bounds(1,rank),my_bounds(2,rank)
            np=np+1
            Aglobal(i,j)=Arecv(np)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 71, 7626, MyFile)
!
      RETURN
      END SUBROUTINE mp_aggregate2d
!
      SUBROUTINE mp_aggregate3d (ng, model, gtype,                      &
     &                           LBiT, UBiT, LBjT, UBjT,                &
     &                           LBiG, UBiG, LBjG, UBjG,                &
     &                           LBk,  UBk,                             &
     &                           Atiled, Aglobal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 3D tiled, floating-point array from each    !
!  spawned node and stores it into 3D global array. If nesting, the    !
!  global array contains the contact points data.                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     gtype      C-grid type.                                          !
!     LBiT       Tiled  array, I-dimension Lower bound.                !
!     UBiT       Tiled  array, I-dimension Upper bound.                !
!     LBjT       Tiled  array, J-dimension Lower bound.                !
!     UBjT       Tiled  array, J-dimension Upper bound.                !
!     LBkT       Tiled  array, K-dimension Lower bound.                !
!     UBkT       Tiled  array, K-dimension Upper bound.                !
!     LBiG       Global array, I-dimension Lower bound.                !
!     UBiG       Global array, I-dimension Upper bound.                !
!     LBjG       Global array, J-dimension Lower bound.                !
!     UBjG       Global array, J-dimension Upper bound.                !
!     LBkG       Global array, K-dimension Lower bound.                !
!     UBkG       Global array, K-dimension Upper bound.                !
!     Atiled     3D tiled, floating-point array to process.            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aglobal    3D global array, all tiles are aggregated.            !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype
      integer, intent(in) :: LBiT, UBiT, LBjT, UBjT
      integer, intent(in) :: LBiG, UBiG, LBjG, UBjG
      integer, intent(in) :: LBk,  UBk
!
      real(r8), intent(in)  :: Atiled(LBiT:UBiT,LBjT:UBjT,LBk:UBk)
      real(r8), intent(out) :: Aglobal(LBiG:UBiG,LBjG:UBjG,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Klen, Lstr, MyError, MyType, Nnodes, Npts, Serror
      integer :: i, j, k, np, rank
      integer,  dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: my_bounds
!
      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)) :: Asend
      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)*                     &
     &                    NtileI(ng)*NtileJ(ng)) :: Arecv
!
      character (len=MPI_MAX_ERROR_STRING) :: string
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/distribute.F"//", mp_aggregate3d"
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 71, 7702, MyFile)
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Maximum automatic buffer memory size in bytes.
!
      BmemMax(ng)=MAX(BmemMax(ng), REAL((SIZE(Asend)+                   &
     &                                   SIZE(Aglobal)+                 &
     &                                   SIZE(Arecv))*KIND(Asend),r8))
!
!  Number of nodes in the group.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
!
!  Set starting and ending indices to process including contact points
!  (if nesting) according to the staggered C-grid classification.
!
      MyType=ABS(gtype)
      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
        CASE (r2dvar, r3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (u2dvar, u3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (v2dvar, v3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
      END SELECT
      Klen=UBk-LBk+1
!
!  Determine the maximum number of points to process between all tiles.
!  In collective communications, the amount of data sent must be equal
!  to the amount of data received.
!
      Npts=0
      DO rank=0,Nnodes
        np=(my_bounds(2,rank)-my_bounds(1,rank)+1)*                     &
     &     (my_bounds(4,rank)-my_bounds(3,rank)+1)*                     &
     &     Klen
        Npts=MAX(Npts, np)
      END DO
      IF (Npts.gt.TileSize(ng)*Klen) THEN
        IF (Master) THEN
          WRITE (stdout,10) ' TileSize = ', TileSize(ng)*Klen, Npts
 10       FORMAT (/,' MP_AGGREGATE3D - communication buffer to small,', &
     &            a, 2i8)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Initialize local arrays to facilitate collective communicatios.
!  This also avoid denormalized values, which facilitates debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Pack tile data.
!-----------------------------------------------------------------------
!
      np=0
      DO k=LBk,UBk
        DO j=my_bounds(3,MyRank),my_bounds(4,MyRank)
          DO i=my_bounds(1,MyRank),my_bounds(2,MyRank)
            np=np+1
            Asend(np)=Atiled(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Aggregate data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,20) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 20     FORMAT (/,' MP_AGGREGATE3D - error during ',a,' call, Node = ', &
     &          i3.3,' Error = ',i3,/,18x,a)
        exit_flag=5
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data into a global 2D array.
!-----------------------------------------------------------------------
!
      DO rank=0,Nnodes
        np=rank*Npts
        DO k=LBk,UBk
          DO j=my_bounds(3,rank),my_bounds(4,rank)
            DO i=my_bounds(1,rank),my_bounds(2,rank)
              np=np+1
              Aglobal(i,j,k)=Arecv(np)
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 71, 7838, MyFile)
!
      RETURN
      END SUBROUTINE mp_aggregate3d
      END MODULE distribute_mod
