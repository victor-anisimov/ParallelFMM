! Parallel Fast Multipole Method
! Copyright (C) 2001-2026  Victor Anisimov
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

#define MEM_SIMULATOR
#undef  MEM_SIMULATOR

program start
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  !
  use timer
  use common_data, only: M2Lene, nRanks, nSplitRanks, szSubComm, rmSize, &
#ifdef MPI
                         MPI_COMM_RMGROUP, MPI_COMM_NODE, MPI_COMM_NZERO, &
#endif
                         M2Lrot, nNodes, nodeSize, prtRank, memUsed, memMax, nodeRank
  implicit none
  !
  integer :: iam
#ifdef MPI
  integer :: mpiError, color
#endif

#ifdef MPI
  call MPI_Init(mpiError)

  ! create the global (world) communicator MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)

  ! create the node-level sub-communicators MPI_COMM_NODE, which combine the ranks on single node
  ! determing the number of ranks per node, nodeSize
  call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_NODE, mpiError)
  call MPI_Comm_rank(MPI_COMM_NODE, nodeRank, mpiError)
  call MPI_Comm_size(MPI_COMM_NODE, nodeSize, mpiError)

  ! determine the number of allocated nodes
  nNodes = nRanks / nodeSize

  ! default number of ranks per node in a single pass in FMM Plan
  nSplitRanks = nodeSize

  ! create a sub-communicator MPI_COMM_NZERO, which connects ranks 0 of each node
  if(nNodes > 1) then
    if(nodeRank == 0) then
      color = 1
    else
      color = MPI_UNDEFINED
    endif
    call MPI_Comm_split(MPI_COMM_WORLD, color, 0, MPI_COMM_NZERO, mpiError)
  endif
#else
  nRanks = 1
  nNodes = 1
  rmSize = 1
  nodeRank = 0
  nodeSize = 1
#endif

  call TimerInit(1)

  ! initialize memory usage counters; use double precision to process integer values
  memUsed = 0.0d0
  memMax  = 0.0d0

#ifdef MEM_SIMULATOR
  ! parameters for memory-usage simulator (10T)
  nodeSize = 12
  nNodes = 2732
  !nNodes = 2731
  nRanks = nNodes * nodeSize
#endif

  ! read input data including prtRank, szSubComm, and nSplitRanks
  call readinp()

#ifdef MEM_SIMULATOR
  ! parameters for the memory-usage simulator (10T)
  szSubComm = 2
#endif

#ifdef MPI
  ! create sub-communicator MPI_COMM_RMGROUP to handle the computation of regular multipoles
  color = iam / (nodeSize * szSubComm)
  call MPI_Comm_split(MPI_COMM_WORLD, color, 0, MPI_COMM_RMGROUP, mpiError)
  rmSize = nodeSize * szSubComm
#endif

#ifdef MEM_SIMULATOR
  ! parameters for the memory-usage simulator (10T)
  rmSize = 1366 * nodeSize
  !rmSize = 2731 * nodeSize
  nSplitRanks = 12
#endif

  ! reset timers for rotation-based multipole translations
  call TimerRotReset

  ! Test the correctness of selected compilation flags for M2L computation
  M2Lene = 0
#ifdef M2LSTD
  M2Lene = M2Lene * 10 + 1
#endif
#ifdef M2LROT
  M2Lene = M2Lene * 10 + 2
#endif
  select case (M2Lene)
    case (1)
      if(iam == prtRank) write(*,'(/a)') "M2L compute mode in energy: standard O(N4)"
    case (2)
      if(iam == prtRank) write(*,'(/a)') "M2L compute mode in energy: rotation O(N3)"
    case default
      if(iam == prtRank) then
        write(*,'(/a)') "Error: Please specify either -DM2LSTD or -DM2LROT in compiler options"
        write(*,'(a,i10)') "You specified: ", M2Lene
      endif
    error stop
  endselect

  M2Lrot = 0
#ifdef KSTD
  M2Lrot = M2Lrot * 10 + 1
#endif
#ifdef MSTD
  M2Lrot = M2Lrot * 10 + 2
#endif
#ifdef KSCA
  M2Lrot = M2Lrot * 10 + 3
#endif
#ifdef MSCA
  M2Lrot = M2Lrot * 10 + 4
#endif
#ifdef MKSCA
  M2Lrot = M2Lrot * 10 + 5
#endif
  select case (M2Lrot)
    case (1)
      if(iam == prtRank) write(*,'(a)') "M2Lrot compute mode: k-standard  rks"
    case (2)
      if(iam == prtRank) write(*,'(a)') "M2Lrot compute mode: m-standard  rms"
    case (3)
      if(iam == prtRank) write(*,'(a)') "M2Lrot compute mode: k-scaled    rkk"
    case (4)
      if(iam == prtRank) write(*,'(a)') "M2Lrot compute mode: m-scaled    rmm"
    case (5)
      if(iam == prtRank) write(*,'(a)') "M2Lrot compute mode: mbar-ktild  rmk"
    case default
      if(iam == prtRank) then
        write(*,'(a)') &
        "Error: M2L rotation routine requires either -DKSTD or -DMSTD or -DKSCA or -DMSCA or -DMKSCA in compiler options"
        write(*,'(a,i10)') "You specified: ", M2Lrot
      endif
      error stop
  endselect

  ! start FMM computation in prg.F90
  call fmm()

#ifdef MPI
  call MPI_Comm_Free(MPI_COMM_NODE, mpiError)
  if(nNodes > 1 .and. MPI_COMM_NZERO /= MPI_COMM_NULL) call MPI_Comm_Free(MPI_COMM_NZERO, mpiError)
  call MPI_Comm_Free(MPI_COMM_RMGROUP, mpiError)
  call MPI_Finalize(mpiError)
#endif
end program start
