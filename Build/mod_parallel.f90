      MODULE mod_parallel
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
      include 'mpif.h'
!
!  Switch to identify master processor. In serial and shared-memory
!  applications it is always true.
!
      logical :: Master
!
!  Switch to identify which thread is processing input/output files.
!  In distributed-memory applications, this thread can be the master
!  thread or all threads in case of parallel output. In serial and
!  shared-memory applications it is always true.
!
      logical :: InpThread
      logical :: OutThread
!
!  Number of shared-memory parallel threads or distributed-memory
!  parallel nodes.
!
      integer :: numthreads
!
!  First and last tile to process in a tiled application.
!
      integer, allocatable :: first_tile(:)
      integer, allocatable :: last_tile(:)
!$OMP THREADPRIVATE (first_tile, last_tile)
!
!  Parallel nodes assined to the ocean model.
!
      integer :: peOCN_frst            ! first ocean parallel node
      integer :: peOCN_last            ! last  ocean parallel node
!
!  Parallel threads/nodes counters used in critical parallel regions.
!
      integer :: tile_count = 0
      integer :: block_count = 0
      integer :: thread_count = 0
!
!  Profiling variables as function of parallel thread:
!
!    proc          Parallel process ID.
!    Cstr          Starting time for program region.
!    Cend          Ending time for program region.
!    Csum          Accumulated time for progam region.
!
      integer, allocatable :: proc(:,:,:)
      real(r8), allocatable :: Cstr(:,:,:)
      real(r8), allocatable :: Cend(:,:,:)
      real(r8), allocatable :: Csum(:,:,:)
!$OMP THREADPRIVATE (proc)
!$OMP THREADPRIVATE (Cstr, Cend)
!
!  Switch manage time clock in "mp_bcasts". During initialization is
!  set to .FALSE. because the profiling variables cannot be allocated
!  and initialized before the "Ngrids" parameter is known.
!
      logical :: Lwclock = .FALSE.
!
!  Distributed-memory master process.
!
      integer, parameter :: MyMaster = 0
!
!  Rank of the parallel local process.
!
      integer :: MyRank = 0
      integer :: MyThread = 0
!$OMP THREADPRIVATE (MyThread)
!
!  Ocean model 1 group communicator handle.
!
      integer :: OCN_COMM_WORLD
!
!  Set mpi_info opaque object handle.
!
      integer :: MP_INFO = MPI_INFO_NULL
!
!  Type of message-passage floating point bindings.
!
      integer, parameter :: MP_FLOAT = MPI_DOUBLE_PRECISION
!
      CONTAINS
!
      SUBROUTINE allocate_parallel
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_strings, ONLY: Nregion
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
!  First and last tile to process in a tiled application.
!
      allocate ( first_tile(Ngrids) )
      allocate ( last_tile (Ngrids) )
!
!  Time profiling variables.
!
      allocate ( proc(0:1,4,Ngrids) )
      proc(0:1,1:4,1:Ngrids)=0
      allocate ( Cstr(0:Nregion,4,Ngrids) )
      Cstr(0:Nregion,1:4,1:Ngrids)=0.0_r8
      allocate ( Cend(0:Nregion,4,Ngrids) )
      Cend(0:Nregion,1:4,1:Ngrids)=0.0_r8
!$OMP END PARALLEL
      allocate ( Csum(0:Nregion,4,Ngrids) )
      Csum(0:Nregion,1:4,1:Ngrids)=0.0_r8
!
! Activate wall clock switch used only in "mp_bcasts". This switch
! is set to .FALSE. during initialization before calling "inp_par.F"
! because the above profiling variables are allocated and initialized
! after the value of "Ngrids" is known.
!
      Lwclock=.TRUE.
      RETURN
      END SUBROUTINE allocate_parallel
      SUBROUTINE initialize_parallel
!
!=======================================================================
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: i
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Initialize shared-memory (OpenMP) or serial configuration.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (1) configuration.
!-----------------------------------------------------------------------
!
!  Get the number of processes in the group associated with the world
!  communicator.
!
      CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' ROMS/TOMS - Unable to inquire number of',           &
     &            ' processors in the group.')
        exit_flag=6
        RETURN
      END IF
!
!  Identify master, input and output threads.
!
      Master=.FALSE.
      InpThread=.FALSE.
      OutThread=.FALSE.
      IF (MyRank.eq.MyMaster) THEN
        Master=.TRUE.
        InpThread=.TRUE.
        OutThread=.TRUE.
      END IF
      RETURN
      END SUBROUTINE initialize_parallel
      END MODULE mod_parallel
