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

! Decide if we do direct sum computation, define - YES, undef - NO
!#define DO_DIRECT_SUM
#undef  DO_DIRECT_SUM

! getLocalFromGlobal() uses more memory but is much faster
#define USE_GETLOCAL_FROMGLOBAL
! findloc() uses less memory but is very slow
!#undef  USE_GETLOCAL_FROMGLOBAL 

! monitor min/max values of various parameters across all MPI ranks
!#define DEBUG_TIMER
#undef  DEBUG_TIMER

! allow mod(PrimitiveEdge,CubeEdge) > 0
#define FLEX_SIDELENGTH
!#undef  FLEX_SIDELENGTH

! enable or disable the report of binding of threads to cpu ids
!#define DEBUG_THREADS
#undef  DEBUG_THREADS

! maximize the performance by disabling the non-essential diagnostics
! one should keep the diagnostics enabled in the day-to-day runs
!#define MAX_PERF
#undef MAX_PERF

module fmm_plan
  implicit none
  save                
  public

  ! box that carries regular (R) multipole moments
  type R_box
!    integer            :: iLevel                            ! 4   hierarchy level
!    integer            :: proc                                ! 4   MPI rank
!    integer (kind = 8) :: parent                              ! 8   parent box index
!    integer (kind = 8) :: idx                               ! 8   box index in R-list
    integer            :: xyz(3)                            ! 4 + 4 + 4
    integer (kind = 8) :: atomStart,  atomEnd !,  nAtoms    ! 8 + 8
!    integer (kind = 8) :: childStart, childEnd !, nChildren   ! 8 + 8 + 8
!    integer (kind = 8) :: nearStart,  nearEnd,  nNear         ! 8 + 8 + 8
!    integer (kind = 8) :: farStart,   farEnd,   nFar          ! 8 + 8 + 8
    ! TODO rename expStart to expLocal
    integer (kind = 8) :: expStart                          ! 8   starting index for box local expansion
    integer (kind = 8) :: expRemote                         ! 8   starting index for box remote expansion
!    double complex, allocatable :: chargeExpansion(:)         ! 16 * (nSzTri+1)
!    double complex, allocatable :: potentialExpansion(:)      ! 16 * (nSzTri+1)
    ! allocated space: 12+16+16 = 44 bytes
  end type R_box

  ! box that carries irregular (S) multipole moments
  type S_box
    integer (kind = 1) :: iLevel                            ! 1   hierarchy level
    integer (kind = 1) :: group                             ! 1   box group
    integer (kind = 8) :: parent                            ! 8   parent box index
!    integer (kind = 8) :: idx                               ! 8   box index in B-list
    integer            :: xyz(3)                            ! 4 + 4 + 4
    integer (kind = 8) :: expStart                          ! 8   starting index for box expansion
    integer (kind = 8) :: atomStart,  atomEnd !,  nAtoms    ! 8 + 8
!    integer (kind = 8) :: childStart, childEnd !, nChildren   ! 8 + 8 + 8
    integer (kind = 8) :: nearStart,  nearEnd !,  nNear     ! 8 + 8
    integer (kind = 8) :: farStart,   farEnd !,   nFar      ! 8 + 8
!    double complex, allocatable :: chargeExpansion(:)         ! 16 * (nSzTri+1)
!    double complex, allocatable :: potentialExpansion(:)      ! 16 * (nSzTri+1)
    ! allocated space: 1+1+8+12+8+16+16+16 = 78 bytes
  end type S_box

  ! binary tree structure
  type Tier
    !integer (kind = 8) :: bStart, bEnd        !  8 + 8 contiguous box ids across all MPI ranks
    double precision   :: edge  !, reciprocal    !  8 + 8 = 32
  end type Tier

  type (R_box), allocatable :: rList(:)
  type (S_box), allocatable :: sList(:)
  type (Tier),  allocatable :: treeLevel(:)

  ! far-field interaction, RMG parallelization, global domain decomposition, and the lowest levels
  integer :: ffLevel, rmLevel, ddLevel, loLevel
  integer (kind=8) :: ffLevel_8, rmLevel_8, ddLevel_8, loLevel_8
  integer (kind=4), allocatable :: segmentSBox(:)

  integer(kind=8), parameter :: nBits = 8   ! number of bits in a byte (kind=1)
  integer(kind=8) :: sizeBytes              ! size of the logical array in bytes

  type Logical_Array
    integer(kind=1), allocatable :: lArray(:)
  contains
    private
    final :: Logical_Array_free   ! destructor
  end type Logical_Array

contains
!
! Constructor
!
function Logical_Array_constructor(sizeBits) result(newInstance)
  implicit none
  integer(kind=8), intent(in) :: sizeBits
  type(Logical_Array) :: newInstance
  !
  ! set the array size
  sizeBytes = sizeBits / nBits + merge(1_8, 0_8, mod(sizeBits,nBits) > 0)
  !
  ! allocate the array in bytes
  allocate(newInstance%lArray(0 : sizeBytes - 1_8))
  !
  ! zero-initialize the array; set all array elements to .false.
  newInstance%lArray = 0
end function Logical_Array_constructor
!
! Destructor
!
subroutine Logical_Array_free(this)
  type(Logical_Array), intent(inout) :: this
  if( allocated(this%lArray) ) then
    deallocate(this%lArray)
  endif
end subroutine Logical_Array_free
!
! Create FMM Plan
!
subroutine CreatePlan()
  !$ use omp_lib
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: nL, ws, nAtomsTotal, aC, aQ, bC, bQ, cC, cQ, cubeEdge, &
            edgeReciprocal, w, cell, nAtoms, nBoxesTotal, nSzTri, vertexMin, vertexMax, &
            bBoxesStart, sBoxLLstart, fBoxes, nearNeighbors, nNearNeighbors, &
            nBoxesSList, maxBoxesSList, rankAtoms, dispAtoms, &
            rankBoxes, dispBoxes, farNeighbors, nFarNeighbors, scaleToUnitSize, &
#ifdef MPI
            MPI_COMM_RMGROUP, MPI_COMM_NODE, MPI_COMM_PNODE, MPI_COMM_NZERO, MPI_COMM_PZERO, &
            MPI_COMM_PGROUP, &
#endif
            nodeRank, memUsed, memMax, memMaxMyNode, memMaxAllNodes, lowestLevel, nThreads, &
            rBoxBottomStart, rBoxBottomEnd, intCubeEdge, intPrimitiveEdge, &
            nLLBoxes_rmRank, segmentRank, segmentNBoxes, segmentGaps, nBoxSegments, &
            nBoxesRList, nRemoteFFBoxes, nParticlesRank, nParticlesRMRank, &
            lComputeDirectSum, particleRanges, nRepLayers, nRepLayers1, nRepLayers2, &
            nReplicas, primitiveCellEdge, leadRank, prtRank, nBoxesBefore_loLevel, &
            nFloatingBoxes, sBoxesTopNum, sBoxesTopDispl, iAtomsBegin, iAtomsEnd, &
            minTopRtreeSize, minTopRtreeRanks, minTopRtreeGroup, expShift, numBoxesSNRF, &
            nNodes, nodeSize, nRanks, rmSize, nSplitRanks, nExtraMultipoleBoxes, debugEachMPIrank, &
            lMultipleDataSets, nDataSets, atomRange, datasetid, nGroups, groupid, &
            nSelections, nSelectedLLboxesRank, ddLevel_input
  use common_constants, only: izero, ione, invLog8
  use timer
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_f_pointer, c_associated
  use allocator  
!$ use omp_lib
  implicit none
  !
  ! Declare the interface for POSIX fsync function
  interface
    function fsync (fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int) :: fsync
    end function fsync
  end interface
  !
#ifdef MPI
  integer            :: mpiError, rankSender, rankReceiver, prtRankSaved, prtNode, color
#endif
  integer            :: pnodeRank, znodeRank, pgroupRank
  integer            :: iam, iRank, jRank, iRankOld, rankPrint=0, rmRank, rmRank1
  integer            :: planSize
  integer            :: allocStat, i,j,k,m,n, nLostAtoms, a, nxb
  integer (kind = 8) :: iAtom, jAtom, kAtom, nextAtom, sumAtoms, globalIndexStart, atomStart, ij, ijk
  integer            :: ix, iy, iz,  jx, jy, jz
  double precision, allocatable   :: tmpCoor(:,:), tmpQ(:)
  integer (kind = 4), allocatable :: nAtomsCube(:), boxRank(:), backupBuffer(:)
  double precision                :: xCenter, yCenter, zCenter, halfCubeEdge
  integer            :: iLevel, iLevelOld, eLevel, nNF, nFF
  integer (kind = 8) :: iBox, jBox, kBox, pBox, tBox, pStart, nA, nB, nC, nD, nF, nP, nLBoxes, i0
  integer (kind = 8) :: aBox, bBox, cBox, bb, gSize, nBoxesBefore, bStart, bEnd, bBefore, bIndex
  integer (kind = 8) :: nChildCubes, iAtomStart, jAtomStart, kAtomStart, mAtomStart
  integer            :: nChildBoxes, nParentEdge, nChildEdge, nParentXY, nChildXY
  double precision   :: edge, reciprocal, memLocal, cm(3)
  integer            :: tjx,tjy,tjz, pjx,pjy,pjz, nbEdge, nbEdgeMAxID, nbStart, nbParent
  integer            :: nbParentEdge_init, nThreadsLocal

  ! TODO try kind=4 with the largest problem size
  integer (kind = 8) :: nbParentEdge, nbParentCube

  integer            :: fax,fay,faz, fbx,fby,fbz
  integer            :: nax,nay,naz, nbx,nby,nbz, nbParentEdgeMaxID, aLevel, sLevel
#ifdef USE_GETLOCAL_FROMGLOBAL
  integer (kind = 4), allocatable :: getLocalFromGlobal(:)
#endif
  integer (kind = 8), allocatable :: getGlobalFromLocal(:)
  integer (kind = 4), allocatable :: boxesRank(:)
  integer (kind = 8), allocatable :: nAtomsR(:), nAtomsNS(:), nAtomsSNR(:)
  logical (kind = 1), allocatable :: btypeNS(:), btypeF(:)
  double precision   :: ffCoeff
  integer            :: indexfb, fBox, xyz(3), nFloatingBoxesEstimated
  integer (kind = 8) :: xBox, yBox, zBox
  integer            :: mBox
  integer (kind = 8) :: nParticlesBuffer, lastLLbox, boxesStart
  integer (kind = 8) :: nParticlesSlist
  integer            :: nBoxesBefore_rmLevel, nBoxesBefore_ffLevel, nBoxesBefore_iLevel
  double precision, allocatable :: bufferC(:,:), bufferQ(:)
  double precision   :: translation(3), sideRatio
  integer, allocatable :: nSegmentsThread(:)  ! number of box segments per OpenMP thread
  integer, allocatable :: nBoxesThread(:)     ! number of boxes handled by each OpenMP thread
  integer, allocatable :: nCubesThread(:)     ! number of boxes handled by each OpenMP thread
  integer              :: iThread, jThread, kThread, bsBox, bsNat, bsRid, sBox, ddBox
  integer (kind=8)     :: iBoxesStart, iBoxesEnd
  integer (kind=8), allocatable :: idx(:), nearNeighbors8(:), farNeighbors8(:)
  integer (kind=8)     :: expStart, expRemote, tmpVal
  integer              :: minRx,minRy,minRz, maxRx,maxRy,maxRz
  integer              :: minNx,minNy,minNz, maxNx,maxNy,maxNz
  integer              :: rBoxTopStart, rBoxTopEnd
  integer (kind=4)     :: xloCell, yloCell, zloCell, xupCell, yupCell, zupCell, lastLLBoxAxis
  integer (kind=4)     :: idx_4, cellSize1, cellSize2, cellSize3
  integer (kind=8)     :: num1_8, num2_8, k_8, nSzTri1_8
  integer (kind=8)     :: nNearNeighborsMin, nNearNeighborsMax, nFarNeighborsMin, nFarNeighborsMax
  type(Logical_Array)  :: boxTypeNS, boxTypeF
  integer (kind=8)     :: nBoxesTotal_div8
  integer (kind=1)     :: dummy
  integer (kind=4)     :: ijkCube, xCube, yCube, zCube
  integer (kind=1)     :: replicaUsed
  integer (kind=4)     :: memFreeNode, memFreeNodeMin
#if defined(DEBUG_TIMER) && defined(MPI)
  integer (kind=8)     :: nReplicasUsedRank, nReplicasUsedRankMin, nReplicasUsedRankMax, nReplicasUsedRankAve
  integer (kind=8)     :: nReplicasRank, nReplicasRankMin, nReplicasRankMax, nReplicasRankAve
  integer (kind=8)     :: numBoxesSNRFMin, numBoxesSNRFMax, numBoxesSNRFAve
  integer (kind=8)     :: nParticlesRMRankMin, nParticlesRMRankMax, nParticlesRMRankAve
  integer (kind=8)     :: nParticlesBufferMin, nParticlesBufferMax, nParticlesBufferAve
  integer (kind=4)     :: nReplicasUsed, nReplicasNeeded

  nReplicasUsed = 0
  nReplicasNeeded = 0
#endif

#ifndef FLEX_SIDELENGTH
  if(mod(intPrimitiveEdge,intCubeEdge) > 0) then
    write(*,'(a)') "Error: mod(intPrimitiveEdge,intCubeEdge) > 0; Recompile with FLEX_SIDELENGTH enabled in plan.F90"
    call exitPrg(100)
  endif
#endif

  numBoxesSNRF = 0     ! number of S+N+R+F boxes
  nParticlesSlist = 0  ! number of particles in LL boxes of S-tree assigned to the rank in ddLevel domain decomposition
  nSelectedLLboxesRank = 0  ! number of selected LL-boxes in S-link in the rank

#ifdef MPI
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_rank(MPI_COMM_RMGROUP, rmRank, mpiError)

  ! RM Group (RMG) is an MPI sub-communicator. Each RMG computes regular multipoles (RM)
  ! of all boxes in a replicated fashion. Each RMG is identical to other RMGs
  ! in details of RM computation. However, each RMG serves the needs of the unique set
  ! of MPI ranks constituting that sub-communicator. Motivation behind this model is
  ! following. Generally speaking, each MPI rank needs access to regular multipoles of
  ! all boxes existing in the system. Distributing RM computation among all ranks would
  ! lead to too much comminication when the ranks need to retrieve the RM data stored in
  ! remote ranks. Such FMM computation would not scale due to communication bottleneck.
  ! Since RM computation is relatively fast, the solution is to split the global set of
  ! MPI ranks into a set of sub-communicators each computing RM data in a replicated
  ! fashion. This model confines the exchange of RM data to within sub-communicators.
  ! This parallelization model eliminates global communications from the FMM computation.

  ! the global MPI space is split into RMGs
  ! each rank iam from MPI world is also an rmRank
  ! rmRank is a rank in RMG

#else
  iam    = 0
  rmRank = 0
  nodeRank = 0
  pnodeRank = 0
  znodeRank = 0
  pgroupRank = 0
#endif
  rmRank1 = rmRank + 1  !   rmRank1 is a unit-based index (1...rmSize) for rmRank (0...rmSize-1)
  !
#if defined (TIMER) || defined (DEBUG_TIMER)
  call TimerInit(2)  ! measures total time spent in this routine
  call TimerInit(3)  ! provides time break down
#endif

#ifdef DEBUG_THREADS
  ! get the cpu id of the running process
  !$omp parallel
    if(iam == prtRank) write(*,'(a,i0,a,i0)') "Thread ", omp_get_thread_num(), " runs on cpu ", getcpuid()
  !$omp end parallel
#endif

  ! the number of ranks in the node that run simultaneously
  planSize = nodeSize  ! default = all ranks in the node

  ! zero-initialize dummy variables
  m = 0
  n = 0

  ! total number of solid harmionics stored in triangular form using zero-based index (0:nSzTri)
  nSzTri = (nL+1)*(nL+2)/2 - 1   ! used in allocate(0:nSzTri), therefore it needs - 1
  nSzTri1_8 = nSzTri + 1         ! int*8
  !
  edgeReciprocal = 1.0d0 / cubeEdge

  sideRatio = cubeEdge / primitiveCellEdge  ! sidelength ratio

  ! hierarchy depth
  loLevel = lowestLevel  ! lowest hierarchy level, LL

  ! ffLevel is the top-most hierarchy level where far field boxes exist, FF
  ffLevel = min0(ceiling(log(dble(ws+2)) / log(2.0d0)), lowestLevel)

  ! The present code employs 2-level parallelization

  ! Both parallelizaion steps start from the same basic principle
  ! Find a hierarchy level that has the number of boxes closely matching the number of ranks
  ! The number of ranks should not exceed the number of boxes at that hierarchy level
  ! The algorith to find iLevel for nRanks:  8^(iLevel-1) < nRanks <= 8^iLevel

  ! 1st level of parallelization performs domain decomposition (DD) across all nRanks
  ! this parallelization happens at ddLevel and deals with far-field interactions
  ! find the upmost ddLevel that has the number of boxes >= nRanks
  ! theoretical min value of ddLevel is ffLevel and the max value is loLevel
  ! going beyond ddLevel=11 requires additional int4->int8 conversion       max         max           min     min
  ddLevel = min0( min0( max0( max0(ceiling(log(dble(nRanks))/log(8.0d0)), ffLevel), ddLevel_input), loLevel), 11)

  ! 2nd level of parallelization generates regular multipoles (RM) for global MPI ranks
  ! this parallelization happens at rmLevel within RMG accross rmSize ranks
  ! rmLevel resides between ffLevel and ddLevel: ffLevel <= rmLevel <= ddLevel
  rmLevel = min0( max0(ceiling(log(dble(rmSize))/log(8.0d0)), ffLevel), ddLevel)

  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(/a,i12)')  "Far-field level,            ffLevel: ", ffLevel  ! FF-level
    write(*,'(a,i12)')   "RMG parallelization level,  rmLevel: ", rmLevel  ! RM-level
    write(*,'(a,i12)')   "Domain decomposition level, ddLevel: ", ddLevel  ! DD-level
    write(*,'(a,i12)')   "The lowest hierarchy level, loLevel: ", loLevel  ! LL-level
  endif

  ! number of boxes before (not including) the hierarchy level
  nBoxesBefore_ffLevel = (8 ** ffLevel - 1) / 7
  nBoxesBefore_rmLevel = (8 ** rmLevel - 1) / 7
  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
  bBoxesStart = nBoxesBefore_loLevel

  ! local index of the last LL box for loLevel
  lastLLbox = 8_8 ** loLevel - 1_8

  ! zero-based grid index of the last loLevel box along a coordinate axis
  lastLLBoxAxis = 2 ** loLevel - 1

  nbParentEdge_init = 2 ** (loLevel - 1)  ! number of boxes along one side of the cube at the parent level

  ! Construct an R-tree by distributing rmLevel boxes among the ranks in RMG sub-communicator
  ! The R-tree begins from the boxes at rmLevel down to loLevel through parent-child relationship
  ! minTopRtreeRanks of RMG ranks (beginning from rmRank 0) will get the minimal number of boxes
  ! The remaining ranks in RMG will get +1 box due to uneven domain decomposition
  !
  ! the minimal number of boxes at rmLevel per rank in RMG
  minTopRtreeSize  = 8**rmLevel / rmSize
  ! the number of ranks in RMG holding the minimal number of boxes
  minTopRtreeRanks = rmSize - 8**rmLevel + minTopRtreeSize*rmSize 
  ! the total number of rmLevel boxes in the minTopRtreeRanks group of RMG ranks
  minTopRtreeGroup = minTopRtreeSize * minTopRtreeRanks

  ! starting and ending global indices of boxes at the top of the R-tree at rmLevel for rmRank
  rBoxTopStart = nBoxesBefore_rmLevel + minTopRtreeSize*rmRank &
               + merge((rmRank-minTopRtreeRanks), 0, rmRank > minTopRtreeRanks)
  rBoxTopEnd   = rBoxTopStart + minTopRtreeSize + merge(0, -1, rmRank >= minTopRtreeRanks)

  ! Compute starting and ending global box indices for loLevel boxes in R-tree for rmRank
  ! It will be used to distribute particles across the ranks in RMG
  ! R-boxes at loLevel hold particles
  ! All R-boxes will be used for computation of regular multipole expansions
  ! The rmRank computes multipole expansions for the boxes of R-list only
  k_8 = 8 ** (loLevel - rmLevel)   ! the number of loLevel boxes per RM parent box
  rBoxBottomStart = nBoxesBefore_loLevel + minTopRtreeSize*k_8*rmRank &
            + merge((rmRank-minTopRtreeRanks)*k_8, 0_8, rmRank > minTopRtreeRanks)
  rBoxBottomEnd   = rBoxBottomStart + minTopRtreeSize*k_8 + merge(k_8, 0_8, rmRank >= minTopRtreeRanks) - 1

  ! number of LL boxes of R-tree to be used in FMM Step 1 for regular multipole computation
  nLLBoxes_rmRank = rBoxBottomEnd - rBoxBottomStart + 1  ! number of LL boxes assigned to rmRank

  ! assign (dummy) starting values for min/max boundaries
  ! lower-left (min) and upper-right (max) grid boundaries of the R-tree boxes at LoLevel on rmRank
  minRx = 9999999
  minRy = 9999999
  minRz = 9999999
  maxRx = 0
  maxRy = 0
  maxRz = 0

  ! lower-left (min) and upper-right (max) grid boundaries of the near-field boxes at loLevel on rmRank
  ! these boundaries automatically include S-boxes
  minNx = 9999999
  minNy = 9999999
  minNz = 9999999
  maxNx = 0
  maxNy = 0
  maxNz = 0

  nParticlesRMRank = 0  ! total number of particles in R boxes assigned to rmRank in RMG

  nParticlesBuffer = 0  ! total number of particle in a tmp buffer, which includes S N R boxes

  ! allocate memory for the binary tree structure
  memLocal = dble( (loLevel+1)*32 )   ! memLocal tracks the amount of memory used in CreatePlan routine
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan treeLevel", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(treeLevel(0:loLevel), stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  ! total number of boxes
  ! Sum of the power series Sn = 8^0 + 8^1 + 8^2 + ... + 8^n = sum[k=0 to n](8^k) = (8^(n+1) - 1)/(8 - 1)
  loLevel_8 = loLevel
  nBoxesTotal = (8_8 ** (loLevel_8 + 1_8) - 1_8) / 7_8

  !!  iLevel  maxCubes             nCubes                 nCubes !!
  !!          along edge           in the level           total  !!
  !!    0       1 =  2^0              1 =    1^3              1  !!
  !!    1       2 =  2^1              8 =    2^3              9  !!
  !!    2       4 =  2^2             64 =    4^3             73  !!
  !!    3       8 =  2^3            512 =    8^3            585  !!
  !!    4      16 =  2^4           4096 =   16^3           4681  !!
  !!    5      32 =  2^5          32768 =   32^3          37449  !!
  !!    6      64 =  2^6         262144 =   64^3         299593  !!
  !!    7     128 =  2^7        2097152 =  128^3        2396745  !!
  !!    8     256 =  2^8       16777216 =  256^3       19173961  !!
  !!    9     512 =  2^9      134217728 =  512^3      153391689  !!
  !!   10    1024 = 2^10     1073741824 = 1024^3     1227133513  !!
  !!   11    2048 = 2^11     8589934592 = 2048^3     9817068105  !!

  ! find the closest number from above to nBoxesTotal divisible by 8
  nBoxesTotal_div8 = (nBoxesTotal / 8_8 + merge(1_8, 0_8, mod(nBoxesTotal,8_8) > 0_8)) * 8_8

  ! B-list of boxes (bStart, bEnd, nBoxes) includes all boxes before domain decomposition
  ! S-list of boxes (sStart, sEnd, sBoxesTopNum) includes the boxes assigned to the current MPI rank

  ! level 0 holds the master box
  treeLevel(0)%edge       =  cubeEdge * dble(2_8 ** loLevel_8)
  !treeLevel(0)%reciprocal =  1.0d0 / treeLevel(0)%edge
  !treeLevel(0)%bStart     =  0 ! starts from the master box
  !treeLevel(0)%bEnd       =  0 ! ends with the master box

  ! initialize the remaining hierarchy levels in the binary tree structure
  do iLevel = 1, loLevel
    ! child's box edge is twice smaller than the parent's box edge
    treeLevel(iLevel)%edge       = treeLevel(iLevel-1)%edge * 0.5d0
    !treeLevel(iLevel)%reciprocal = treeLevel(iLevel-1)%reciprocal * 2.0d0

    !treeLevel(iLevel)%bStart     = (8 ** iLevel - 1) / 7             ! stating index of the box in B-list
    !treeLevel(iLevel)%bEnd       = (8 ** (iLevel+1) - 1) / 7 - 1     ! ending index of the box in B-list
  enddo

  if(iam == prtRank .or. debugEachMPIrank) then
    ! at this point, nB is the number of boxes along the edge at the lowest hierarchy level
    write(*,'(/a,i8)') "Well-separation:      ", ws
    write(*,'( a,i0)') "Theoretical number of boxes: ", nBoxesTotal   ! including the master box
    nB = 2 ** loLevel
    write(*,'(a,8x,i0," * ",i0," * ",i0)') "Boxes along the cube edge:", nB, nB, nB
  endif

  if(loLevel < 2) then
    ! minimal initialization to test PBC implementation 

    ! allocate memory for boxes
    memLocal = memLocal + dble(2*44) + dble(nAtomsTotal)*4*8
    memMax = dmax1(memMax, memUsed)
    if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan bC bQ", memLocal, memUsed, memMax, planSize)
    allocate(rList(0:1), bC(3,nAtomsTotal), bQ(nAtomsTotal), stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)

    ! assemble the box hierarchy in compact form (include only populated boxes)

    ! level = 0; master box; iBox = 0
    rList(0)%xyz(1)     =  0
    rList(0)%xyz(2)     =  0
    rList(0)%xyz(3)     =  0
    rList(0)%atomStart  =  1
    rList(0)%atomEnd    =  nAtomsTotal

!    TODO This part is broken and needs fixing
!    do iAtom = 1, nAtomsTotal
!      bC(1,iAtom) = coor(1,iAtom)
!      bC(2,iAtom) = coor(2,iAtom)
!      bC(3,iAtom) = coor(3,iAtom)
!      bQ(iAtom) = q(iAtom)
!    enddo
    return
  endif

  ! from this point down we have a regular FMM case

#ifdef MPI
  if(nSplitRanks < nodeSize) then
    ! purpose: reduce the memory usage in the node
    ! restrict the number of simultaneously running MPI ranks per node to nSplitRanks in the FMM Plan
    ! first pass includes nodeRank(s) in the range: 0 ... nSplitRanks-1 per node
    ! second pass includes nodeRank(s) in the range: nSplitRanks ... 2*nSplitRanks-1 < nodeSize
    ! and so on
    !
    ! connect nSplitRanks in the pass accross all nodes into a sub-communicator MPI_COMM_PGROUP
    color = nodeRank / nSplitRanks
    call MPI_Comm_Split(MPI_COMM_WORLD, color, 0, MPI_COMM_PGROUP, mpiError)
    call MPI_Comm_rank(MPI_COMM_PGROUP, pgroupRank, mpiError)
    !
    ! split each single-node sub-communicator MPI_COMM_NODE into smaller sub-communicators MPI_COMM_PNODE
    ! MPI_COMM_PNODE combines ranks in the node included in a pass
    ! each nSplitRanks ranks in the node form a separate sub-communicator
    ! all nodes in the entire node-allocation undergo the same decomposition
    color = nodeRank / nSplitRanks
    call MPI_Comm_Split(MPI_COMM_NODE, color, nodeRank, MPI_COMM_PNODE, mpiError)
    !
    ! connect the ranks pnodeRank = 0 from the same pass across all nodes
    ! create a sub-communicator MPI_COMM_PZERO, which connects ranks 0 (pnodeRank) in each pass
    ! distinguish each pass in the node one from another by color = nodeRank
    ! ranks with the same color represent pnodeRank = 0 from the same pass accross all nodes
    call MPI_Comm_rank(MPI_COMM_PNODE, pnodeRank, mpiError)
    if(pnodeRank == 0) then
      color = nodeRank  ! distinguish each pass in the node one from another by color (nodeRank)
    else
      color = MPI_UNDEFINED
    endif
    call MPI_Comm_Split(MPI_COMM_WORLD, color, iam, MPI_COMM_PZERO, mpiError)
    if(MPI_COMM_PZERO /= MPI_COMM_NULL) then
      call MPI_Comm_rank(MPI_COMM_PZERO, znodeRank, mpiError)
    else
      znodeRank = -1  ! the value is undefined for non-existing sub-communicator
    endif
    !
    ! assign prtRank to each pass; for this to work well, it is best to set prtRank=0
    prtRankSaved = prtRank        ! save the original value of prtRank in prtRankSaved
    prtNode = prtRank / nodeSize  ! determine the node which has prtRank
    prtRank = prtNode * nodeSize  ! nodeRank == 0 on the node prtNode becomes prtRank for the first pass
    if(iam == prtRank .or. debugEachMPIrank) write(*,'(/a,i0, a,i0,a)') "FMM Plan pass 1 out of ", &
                                                  ceiling(dble(nodeSize) / dble(nSplitRanks)), &
                                                  " using ", nSplitRanks, " ranks per node"

    ! assign the number of ranks in the node to each pass
    planSize = nSplitRanks        ! the first pass includes nSplitRanks

    ! the nodeRank(s) in the range 0 ... nSplitRanks-1 continue execution without waiting for a wake-up signal
    ! the remainder of ranks will wait until the ranks from the previous phase (rankSender) finish their work
    if(nodeRank >= nSplitRanks) then
      ! suspend the nodeRank if it is outside the range until it gets awaken by rankSender
      rankSender = nodeRank - nSplitRanks  ! rankSender notifies nodeRank (iam) to exit the wait mode and go ahead
      call MPI_Recv(i, 1, MPI_INTEGER, rankSender, rankSender, MPI_COMM_NODE, MPI_STATUS_IGNORE, mpiError)
      !
      ! if we got here, the nodeRank has received a wake-up call from the rankSender
      !
      ! assign prtRank to the pass number nodeRank / nSplitRanks
      prtRank = prtNode * nodeSize + (nodeRank / nSplitRanks) * nSplitRanks
      if(iam == prtRank .or. debugEachMPIrank) &
              write(*,'(/a,i0, a,i0, a,i0,a)') "FMM Plan pass ", nodeRank / nSplitRanks + 1, &
              " out of ", ceiling(dble(nodeSize) / dble(nSplitRanks)),                                  &
              " using ", nSplitRanks, " ranks per node"
      !
      ! assign the number of ranks to the pass number nodeRank / nSplitRanks
      planSize = merge(nSplitRanks,                                   &  ! every pass except the last one
                       nodeSize - (nodeRank/nSplitRanks)*nSplitRanks, &  ! last pass
                       nodeRank / nSplitRanks + 1 < ceiling(dble(nodeSize) / dble(nSplitRanks)) )
      !
      ! reset the break-down timer for the just awaken MPI rank
      call TimerInit(3)  ! provides time break down
    endif
  else
    ! execute the rest of the FMM Plan in a single pass
    if(iam == prtRank .or. debugEachMPIrank) &
            write(*,'(/a,i0,a)') "FMM Plan single pass using ", nSplitRanks, " ranks per node"
  endif
#endif

#ifndef MAX_PERF
  ! Report free (host) memory on the compute node
  memFreeNode = stat_allmemory() / 1048576_8  ! 1024*1024, reported in bytes and converted to Mb
  memFreeNodeMin = memFreeNode
#ifdef MPI
  if(nNodes > 1 .and. .not. debugEachMPIrank) then
    if(nSplitRanks == nodeSize) then
      if(MPI_COMM_NZERO /= MPI_COMM_NULL) &
      call MPI_Reduce(memFreeNode, memFreeNodeMin, 1, MPI_INTEGER4, MPI_MIN, leadRank, MPI_COMM_NZERO, mpiError)
    else 
      if(MPI_COMM_PZERO /= MPI_COMM_NULL) &
      call MPI_Reduce(memFreeNode, memFreeNodeMin, 1, MPI_INTEGER4, MPI_MIN, leadRank, MPI_COMM_PZERO, mpiError)
    endif
  endif
#endif
  if(iam == prtRank .or. debugEachMPIrank) write(*,'(a,i7,a,i7,a)') &
  "Free Host Memory in the nodes              Local:", memFreeNode, " Mb    Min:", memFreeNodeMin, " Mb"
#endif

  ! Estimate the max number of floating boxes per rank (needed for the case when ddLevel /= ffLevel)
  ! Floating boxes are parents, parents of parents, etc of ddLevel boxes 
  ! The parental leneage ends at ffLevel
  ! Here, we estimate a number of floating boxes owned by the rank
  ! Assume that each ddLevel box in the rank has a unique parent
  ! Parallelization applies to the boxes at ddLevel accross the ranks in MPI World
  nFloatingBoxesEstimated = ((8_8 ** ddLevel) / nRanks + 1) * (ddLevel - ffLevel)

  ! allocate memory to track the number of boxes and particles per rank
  ! allocate memory for boxes
  ! add an extra box at the end(+1) since the boxes are indexed from 0:nBoxesTotal
  ! allocate memory for particle coordinated and weights
  ! sBoxesTopNum is a number of ddLevel boxes in the S-tree for each MPI rank
  ! sBoxesTopDispl is a displacement in the number of ddLevel boxes in the S-tree accross MPI ranks
  ! fBoxes is a list of floating boxes for the rank iam
  ! particleRanges(1:2,0:nLLBoxes_rmRank-1) determines beginning and ending particle index in 
  !   the particle coordinate array cC for LL boxes of R-tree to be used in FMM Step 1.
  !   The arrays cC, cQ, and particleRanges will be deallocated after FMM Step 1.
  memLocal = memLocal + dble(nRanks*(4+8)) + dble(nRanks*8*4) + dble((nFloatingBoxesEstimated+1)*8) &
           + dble(rmSize*4) &
           + dble(nBoxesTotal_div8)*2 &             ! btypeNS + btypeF
           + dble(nBoxesTotal_div8/nBits)*2 + &     ! boxTypeNS + boxTypeF
#ifdef USE_GETLOCAL_FROMGLOBAL
           ! account for memory usage of getLocalFromGlobal
           dble(nBoxesTotal)*4 + &
#endif
           dble(nLLBoxes_rmRank*8*2)
  memMax = dmax1(memMax, (memUsed + memLocal))
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan + nBoxesTotal", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(sBoxesTopNum(0:nRanks-1), sBoxesTopDispl(0:nRanks-1),       &
           rankAtoms(0:nRanks-1), dispAtoms(0:nRanks-1), &
           rankBoxes(0:nRanks-1), dispBoxes(0:nRanks-1), &
           fBoxes(nFloatingBoxesEstimated+1), & ! need extra 1 to avoid out-of-boundary in fBoxes(nFloatingBoxes+1) = zBox
           nExtraMultipoleBoxes(0:rmSize-1),             &
           btypeNS(0:nBoxesTotal_div8-1),                &
           btypeF(0:nBoxesTotal_div8-1),                 &
#ifdef USE_GETLOCAL_FROMGLOBAL
           getLocalFromGlobal(0:nBoxesTotal-1),          &
#endif
           particleRanges(1:2,0:nLLBoxes_rmRank-1),      &
           stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  boxTypeNS = Logical_Array_constructor(nBoxesTotal_div8)  ! allocate nBoxesTotal_div8/nBits bytes for N,S boxes
  boxTypeF  = Logical_Array_constructor(nBoxesTotal_div8)  ! allocate nBoxesTotal_div8/nBits bytes for F boxes

  btypeNS = .false.   ! includes the S and N boxes at the loLevel only; boxes that carry particles
  btypeF  = .false.   ! includes F boxes from ffLevel..loLevel
#ifdef USE_GETLOCAL_FROMGLOBAL
  getLocalFromGlobal = 0  ! array to convert global box index to a local box index (local to MPI rank)
#endif

  !purpose: Determine the number of E-boxes carrying multipole moments (when rmLevel /= ffLevel) per rank

  ! Extra multipole (E) boxes are parents, parents of parents, etc of rmLevel boxes 
  ! Only one box from ffLevel when ffLevel < rmLevel becomes assigned to rmRank
  ! Boxes from ffLevel become assigned to rmRanks in sequential order
  ! Boxes between ffLevel and rmLevel form E-tree starting from the top box at ffLevel and ending at rmLevel-1
  ! Each E-tree becomes assigned to a separate sequential rmRank
  ! nExtraMultipoleBoxes(rmRank) reserves space for E-boxes (and their children) in boxRegExpansion
  ! The lowest child E-boxes hold a copy of data from rmLevel and are temporary placeholders in boxRegExpansion
  ! TODO replace this trivial array with a conditional test whenever it is needed
  nExtraMultipoleBoxes(0:rmSize-1) = 0  ! the number of E-boxes in E-tree

  if(rmLevel > ffLevel) then
    ! each rmRank computes their extra multipole boxes from ffLevel...rmLevel-1
    ! E-tree reserves 9 storage units (boxes) in boxRegExpansion when rmLevel = ffLevel+1
    ! 9 = 1 (for parent from ffLevel) + 8 (for its children holding a copy of the data from rmLevel)
    i = (8 ** ffLevel) - 1  ! the box max index in ffLevel (it never exceeds rmSize)
    nExtraMultipoleBoxes(0:i) = (8 ** (rmLevel-ffLevel+1) - 1) / 7  ! the number of E-boxes in E-tree
  endif

  !if(iam == prtRank) then
  !  do i = 0, rmRank-1
  !    write(*,'(2i3)') i, nExtraMultipoleBoxes(i)
  !  enddo
  !endif
  !stop

  ! Floating boxes are part of the S-list. They expand the S-tree from its top ddLevel up to
  !   ffLevel through a parental lineage. E-boxes are the sources of far-field potential
  !   for the floating boxes located above rmLevel.

  rankAtoms = 0  ! number of particles assinged to an MPI rank
  dispAtoms = 0  ! displacement in the number of particles for the ranks; rank 0 has displacement 0

  !purpose: Start 1st-level of parallelization by distributing ddLevel boxes across global MPI ranks

  ! start building the S-tree, which is associated with ddLevel
  ! each MPI rank knows the S-tree structure on other ranks from the data stored in sBoxesTopNum and sBoxesTopDispl
  ! each global rank iRank will handle a unique list of ddLevel boxes
  ! from its ddLevel information, the iRank can find children, children of children, etc for computation
  k_8 = max0(1_8, 8_8**ddLevel / nRanks)  ! the minimal number of boxes at ddLevel per rank
  j = 8_8**ddLevel - k_8 * nRanks       ! the total number of excess boxes in the uneven workload
  i = nRanks - j                    ! (iRank >= i) is the rank index to assign the extra box
  do iRank = 0, nRanks-1
    sBoxesTopNum(iRank) = k_8 + merge(1, 0, (iRank >= i))                 ! the number of boxes at ddLevel assigned to iRank
    sBoxesTopDispl(iRank) = k_8 * iRank + merge(iRank-i, 0, (iRank > i))  ! box starting local index at ddLevel for iRank
  enddo

!  ! number of boxes from ffLevel to ddLevel-1 (floating boxes)
!  nFloatingBoxes = merge((8**ffLevel) * (8**(ddLevel-ffLevel)-1)/7, 0, (ddLevel > ffLevel))
!  ! minimal number of floating (above ddLevel) boxes per rank
!  m = merge(max0(1, nFloatingBoxes / nRanks), 0, (ddLevel > ffLevel)) 
!  ! number of excess (wrap around) floating boxes per rank
!  n = merge(nFloatingBoxes - m * nRanks, 0, (nFloatingBoxes > nRanks)) 
!  do iRank = 0, rmSize-1
!    fBoxes(iRank) = merge(m, 0, (iRank < nFloatingBoxes)) &     ! floating boxes from ddLevel-1 and up
!                  + merge(1, 0, (iRank < n))                    ! the remainder of floating boxes wrapped around
!  enddo

  nFloatingBoxes = 0  ! there are no floating boxes in a serial (non-MPI) mode
#ifdef MPI
  ! floating box is a parent box (goes up to ffLevel) of ddLevel boxes on the top of S-tree
  ! floating box belongs to S-list and participates in far-field interactions

  !purpose: Create a list of floating boxes for the rank iam

  ! not using OpenMP parallelization due to sequential population of the array fBoxes
  ! for ddLevel > 10, sBoxesTopDispl must be int8, sBoxesTopNum can stay int4
  ! the box indices xBox, yBox, zBox, and the array fBoxes must be int8 for ddLevel > 10
  do xBox = sBoxesTopDispl(iam), sBoxesTopDispl(iam) + sBoxesTopNum(iam) - 1
    yBox = xBox  ! local index of the ddLevel box in B-list

    ! find the global index of all parents of yBox up to and including FF level
    do iLevel = ddLevel-1, ffLevel, -1
      ! get the local index of the parent
      yBox = yBox / 8
      ! get the global index of the parent
      zBox = yBox + (8**iLevel - 1) / 7
      ! need extra +1 in the array allocation to avoid out-of-boundary in the assignment
      fBoxes(nFloatingBoxes+1) = zBox
      ! check the presence of the newly found zBox in the array fBoxes
      j = 1
      do i = 1, nFloatingBoxes
        if(fBoxes(i) == zBox) j = 0
      enddo
      nFloatingBoxes = nFloatingBoxes + j
    enddo
    call quicksort8(fBoxes(1:nFloatingBoxes))
  enddo

  if(nFloatingBoxes > nFloatingBoxesEstimated) then
    write(*,'(3(a,i0))') "Error in FMM Plan: The actual value of nFloatingBoxes ", &
          nFloatingBoxes, " exceeds the estimation ", nFloatingBoxesEstimated, " on rank ", iam
    call MPI_Abort(MPI_COMM_WORLD, -1, mpiError)
  endif
#endif

  ! nBoxesRList is the number of boxes assigned to rmRank in RMG domain decomposition
  ! it includes:
  ! +1 is for the master box
  ! + nExtraMultipoleBoxes to account for rmLevel > ffLevel
  ! + the boxes of R-tree (from rmLevel down to loLevel)
  ! nBoxesRList determines the size of boxRegExpansion
  nBoxesRList = 1 + nExtraMultipoleBoxes(rmRank) &        ! E-boxes (from E-tree)
              + (minTopRtreeSize + merge(1, 0, rmRank >= minTopRtreeRanks)) &
              * (8**(loLevel-rmLevel+1) - 1) / 7

  ! Array boxRegExpansion(0:(nSzTri+1)*(nBoxesRList+nRemoteFFBoxes)-1) holds regular multipole
  !   expansions of the boxes used in rmRank. Each box occupies (nSzTri + 1) array elements 
  !   in boxRegExpansion. The boxes (their regular multipole expansions) follow the order:
  ! 1) box with index 0 - virtual master box holds no multipole expansion (empty placeholder)
  ! 2) followed by nExtraMultipoleBoxes 
  ! 3) boxes with index (nExtraMultipoleBoxes+1 : nBoxesRList-1) represent boxes of R-tree
  !      their index value goes from top to bottom of the R-tree
  ! 4) remote far-field (FF) boxes occupy the remainder of boxRegExpansion

  ! can stay int4 with ddLevel=11 since sBoxesTopNum can fit in int4
  ! nBoxesSList is a total number of boxes in the S-list: 1 master, floating boxes, S-tree boxes
  nBoxesSList = 1 + nFloatingBoxes + sBoxesTopNum(iam) * (8 ** (loLevel - ddLevel + 1) - 1) / 7

  ! can stay int4 with ddLevel=11 since sBoxesTopNum can fit in int4
  ! max number of boxes per MPI rank (last rank carries the highest workload)
  ! we do not know nFloatingBoxes for the rank nRanks-1 so we approximate maxBoxesSList
  maxBoxesSList = 1 + nFloatingBoxes + sBoxesTopNum(nRanks-1) * (8 ** (loLevel - ddLevel + 1) - 1) / 7

  ! allocate memory for boxes, S-list associated with MPI rank
  ! gather those boxes in array sList; possibly, do not need a separate array
  memLocal = memLocal + dble(nBoxesSList*78) + dble(nBoxesSList*8) + &
#ifdef DO_DIRECT_SUM
             dble(((max0(8_8**loLevel_8,nBoxesSList*2_8)+1_8))*4_8*3_8)
#else
             dble((nBoxesSList*2_8+1_8)*4_8)
#endif
  memMax = dmax1(memMax, (memUsed + memLocal))
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan + nBoxesSList", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(sList(0:nBoxesSList-1), idx(0:nBoxesSList-1),     &
#ifdef DO_DIRECT_SUM
           nAtomsCube(0_8:max0(8_8**loLevel_8,nBoxesSList*2_8)),     &
           boxRank(0_8:max0(8_8**loLevel_8,nBoxesSList*2_8)),        &
           backupBuffer(0_8:max0(8_8**loLevel_8,nBoxesSList*2_8)),   &
#else
           boxRank(0:nBoxesSList*2),                                 &
#endif
           stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  ! initialize the master box in S-list
  sList(0)%iLevel   = 0
  sList(0)%group    = 0
  sList(0)%parent   = 0
 !sList(0)%idx      = 0
  sList(0)%xyz(1)   = 0
  sList(0)%xyz(2)   = 0
  sList(0)%xyz(3)   = 0
  sList(0)%expStart = 0
  idx(0)            = 0

  ! sBoxLLstart is a starting local S-list index of S-boxes at loLevel
  ! It shows how many boxes in the S-list are before those at loLevel
  ! The S-list includes a master box, floating boxes, and the S-tree boxes
  ! The local index appears after renumbering the boxes in S-list assigned to an MPI rank
  sBoxLLstart = 1 + nFloatingBoxes + sBoxesTopNum(iam) * (8 ** (loLevel - ddLevel) - 1) / 7

  ! Global starting index of S-boxes at loLevel
  num1_8 = 8_8 ** (loLevel - ddLevel)
  iBoxesStart = sBoxLLstart + nBoxesBefore_loLevel + sBoxesTopDispl(iam) * num1_8 &
              - ((num1_8 - 1_8)/7_8 * sBoxesTopNum(iam) + 1_8) - nFloatingBoxes

  ! Global index of the last S-box at loLevel
  iBoxesEnd = iBoxesStart + nBoxesSList - sBoxLLstart - 1

  nChildEdge = 2 ** loLevel            ! number of boxes along the edge at the lowest level
  nChildXY   = nChildEdge * nChildEdge ! number of boxes in XY-plane at the lowest level

#ifdef DO_DIRECT_SUM
  nAtomsCube = 0                       ! zero-initialize the array
#endif

  ! box types:
  ! bit values for S N M F boxes:
  ! box    bits   numerical value
  ! S box: 1000   1  boxes from S-list
  ! N box: 0100   2  near-field boxes (carry particles) at loLevel to those from S-tree
  ! R box: 0010   4  multipole boxes (carry particles) at loLevel from R-tree
  ! F box: 0001   8  far-field boxes to those from S-list
  bsBox =  4  ! box type
  bsNat =  0  ! presently not used
  bsRid = 28  ! rank id and (later) local box id after box renumbering

  ! we know the box ranges but they do not have assigned particles yet
  ! TODO we already know the box ranges (the ranges are fixed) so the following is unnecessary
  ! TODO we can use a smaller array just for LL boxes
  ! mark the boxes which carry particles for rmRank in node-based domain decomposition

  nNearNeighbors = 0  ! number of boxes in the near-neighbor list in the MPI rank
  nFarNeighbors  = 0  ! number of boxes in the far-neighbor list in the MPI rank

  ! starting global index of iBox for S-boxes at ddLevel
  globalIndexStart = (8_8**ddLevel - 1_8) / 7_8 + sBoxesTopDispl(iam) - 1_8

  !$omp parallel

  !purpose: Determine the min-max grid values of the R-tree top boxes at rmLevel

  ! these boundaries will later be projected to loLevel
  !$omp do private(jBox,ix,iy,iz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz) &
  !$omp reduction(min:minRx,minRy,minRz) reduction(max:maxRx,maxRy,maxRz)
  do iBox = rBoxTopStart, rBoxTopEnd
    ! determine ix, iy, iz of iBox located on ffLevel from the box index
    jBox = iBox - nBoxesBefore_rmLevel  ! jBox is a local index of iBox for rmLevel
    iz = 0
    iy = 0
    ix = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = 2 ** (rmLevel - 1)
    do aLevel = 1, rmLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / nbParentCube        ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      iz = iz + tjz * nbParentEdge
      iy = iy + tjy * nbParentEdge
      ix = ix + tjx * nbParentEdge

      jBox = mod(jBox,nbParentCube)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    minRx = min0(minRx, ix)
    minRy = min0(minRy, iy)
    minRz = min0(minRz, iz)

    maxRx = max0(maxRx, ix)
    maxRy = max0(maxRy, iy)
    maxRz = max0(maxRz, iz)
  enddo

  !purpose: Initialize floating S-boxes (ffLevel..ddLevel-1)

  ! sList(sBox)%expStart points to boxIrrExpansion
  ! boxIrrExpansion begins with multipole expansions of floating boxes
  !$omp do private(sBox,iLevel,jBox,ix,iy,iz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz,i)
  do sBox = 1, nFloatingBoxes  ! loop over floating S-boxes at the ffLevel..ddLevel-1
    ! iBox = fBoxes(sBox)  ! B-box
    ! get iLevel from iBox by using the inequality iLevel <= log(iBox*7+1)/log8 < iLevel+1
    iLevel = floor(dlog(dble(fBoxes(sBox)*7_8 + 1_8)) * invLog8)  ! iLevel = loLevel when ddLevel = loLevel

    ! get ix, iy, iz from the box index
    ! must use int8 math when iLevel = ddLevel > 10
    jBox = fBoxes(sBox) - (8_8**iLevel - 1_8) / 7_8   ! get the local index of iBox
    iz = 0
    iy = 0
    ix = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = 2 ** (iLevel - 1)
    do aLevel = 1, iLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / nbParentCube        ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      iz = iz + tjz * nbParentEdge
      iy = iy + tjy * nbParentEdge
      ix = ix + tjx * nbParentEdge

      jBox = mod(jBox,nbParentCube)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    jBox = fBoxes(sBox)        ! global index of the box in B-list
    !box(jBox)%idx = sBox      ! let the B-box know its local S-index (sBox)

    ! find the parent of sBox in S-list
    pBox =  (8**(iLevel-1)-1)/7 + (jBox - (8_8**iLevel - 1_8)/7_8)/8_8  ! the parent box index in B-list
    do i = 1, nFloatingBoxes
      if(pBox == fBoxes(i)) then
        pBox = i
        exit
      endif
    enddo

    ! the boxes at ffLevel do not need to know their actual parent
    sList(sBox)%iLevel   = iLevel
    sList(sBox)%group    = 0
    sList(sBox)%parent   = merge(0_8, pBox, iLevel == ffLevel)
   !sList(sBox)%idx      = jBox
    sList(sBox)%xyz(1)   = ix
    sList(sBox)%xyz(2)   = iy
    sList(sBox)%xyz(3)   = iz
    sList(sBox)%expStart = sBox * (nSzTri + 1)  ! using int*4 is fine since it is the beginnig of the array
    idx(sBox)            = jBox  ! save the global box index for sBox
  enddo

  !purpose: Initialize S-boxes at ddLevel

  ! create a consecutive list of boxes (S-list) associated with MPI rank
  ! the MPI ranks do not need boxes above ffLevel (ddLevel)
  ! S-list consistes of the master box, floating boxes, boxes at ddLevel, and child boxes down in the hierarchy
  ! the total number of boxes above ddLevel, plus the shift due to domain decomposition, minus one since ddBox=1
  !$omp do private(jBox,ix,iy,iz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz,sBox,i)
  do ddBox = 1, sBoxesTopNum(iam)  ! loop over S-boxes at the ddLevel
    ! get ix, iy, iz from the box, ddBox
    jBox = sBoxesTopDispl(iam) - 1 + ddBox   ! get the local index of the box at ddLevel for the MPI rank iam
    iz = 0
    iy = 0
    ix = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = 2 ** (ddLevel - 1)
    do aLevel = 1, ddLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / nbParentCube        ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      iz = iz + tjz * nbParentEdge
      iy = iy + tjy * nbParentEdge
      ix = ix + tjx * nbParentEdge
    
      jBox = mod(jBox,nbParentCube)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    !sBox = ddBox
    sBox = nFloatingBoxes + ddBox       ! shift ddLevel boxes in the S-list index by adding floating boxes
    jBox = globalIndexStart+ddBox       ! global index of sBox in B-list
    !box(jBox)%idx = sBox               ! let the B-box (globalIndexStart+ddBox) know its S-index (sBox)

    ! find the parent of sBox in S-list among floating boxes
    ! boxes at ffLevel do not need to know their parent
    !pBox = (jBox - (8**ddLevel-1)/7) / (8**(ddLevel-1)) + (8**(ddLevel-1)-1)/7 ! the parent box index in B-list
    ! TODO use pBox4
    pBox =  (8**(ddLevel-1)-1)/7 + (jBox - (8_8**ddLevel - 1_8)/7_8)/8_8  ! the parent box index in B-list
    do i = 1, nFloatingBoxes
      if(pBox == fBoxes(i)) then
        pBox = i   ! i is a local box index in S-list
        exit
      endif
    enddo

    ! the boxes at ffLevel do not need to know their actual parent
    sList(sBox)%iLevel   = ddLevel
    sList(sBox)%group    = 0
    sList(sBox)%parent   = merge(0_8, pBox, ddLevel == ffLevel)
   !sList(sBox)%idx      = jBox
    sList(sBox)%xyz(1)   = ix
    sList(sBox)%xyz(2)   = iy
    sList(sBox)%xyz(3)   = iz
    sList(sBox)%expStart = sBox * (nSzTri + 1)
    idx(sBox)            = jBox

    !if(iam == prtRank) write(*,'(a,3i4,4x,20i4)') "debug: sBox, idx, pBox,  fBoxes: ", sBox, sList(sBox)%idx, pBox, fBoxes
  enddo

  !purpose: iInitialize the rest of the S-list starting from ddLevel+1 and down

  if(ddLevel < loLevel) then
  !$omp do private(kBox,slevel,iLevel,pBox,iBox,jBox,ix,iy,iz,nbParentEdge,aLevel,nbParentCube,aBox,tjx,tjy,tjz)
  do kBox = 1 + sBoxesTopNum(iam), nBoxesSList-1 - nFloatingBoxes
    ! get iLevel from kBox
    ! sLevel is a local hierarchy level for kBox
    sLevel = floor(dlog(dble(((kBox - 1) / sBoxesTopNum(iam)) * 7 + 1)) * invLog8)
    iLevel = sLevel + ddLevel  ! global hierarchy level for kBox

    ! get pBox from kBox and sLevel
    pBox = (8_8**(sLevel-1) - 1_8) / 7_8 * sBoxesTopNum(iam) + 1_8  &
         + (kBox - 1_8 - (8_8**sLevel - 1_8) / 7_8 * sBoxesTopNum(iam)) / 8_8 + nFloatingBoxes

    ! get the global index iBox of kBox
    iBox = (8_8**iLevel - 1_8) / 7_8 + sBoxesTopDispl(iam) * 8_8**sLevel + kBox  &
         -((8_8**sLevel - 1_8) / 7_8 * sBoxesTopNum(iam) + 1_8)

    ! get ix, iy, iz from iBox
    jBox = iBox - (8_8**iLevel - 1_8) / 7_8   ! get the local index of iBox for iLevel
    iz = 0
    iy = 0
    ix = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = 2 ** (iLevel - 1)
    do aLevel = 1, iLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      aBox = jBox / nbParentCube        ! the parent box index, aBox that holds jBox (iBox)
      tjz =  aBox / 4
      tjy = (aBox - tjz * 4) / 2
      tjx =  aBox - tjz * 4 - tjy * 2
      iz = iz + tjz * nbParentEdge
      iy = iy + tjy * nbParentEdge
      ix = ix + tjx * nbParentEdge

      jBox = mod(jBox,nbParentCube)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    jBox = kBox + nFloatingBoxes
    sList(jBox)%iLevel   = iLevel
    sList(jBox)%group    = 0
    sList(jBox)%parent   = pBox
   !sList(jBox)%idx      = iBox
    sList(jBox)%xyz(1)   = ix
    sList(jBox)%xyz(2)   = iy
    sList(jBox)%xyz(3)   = iz
    sList(jBox)%expStart = jBox * (nSzTri + 1)
    idx(jBox)            = iBox

!   sList(kBox)%childStart = (8**(aLevel+1))/7 * sBoxesTopNum(iam) + 1 + (kBox - ((8**aLevel - 1)/7 * sBoxesTopNum(iam) + 1)) * 8
!   sList(kBox)%childEnd   = sList(kBox)%childStart + 7

    !box(iBox)%idx = jBox

    ! it is more convenient to mark the LL boxes in R-list together with NF boxes
  enddo
  endif

  !
  !purpose: assign each LL-box in S-list a selection group if there are any (default group 0 means no selection)
  !purpose: compute the number of selected LL-boxes in S-list in the rank

  if(nSelections > 0) then
    !$omp do private(ix,iy,iz, ijk) reduction(+:nSelectedLLboxesRank)
    do sBox = sBoxLLstart, nBoxesSList-1 ! Loop over boxes of the lowest level
      !
      ! convert box coordinates to the replica coordinates
      ix = dble(sList(sBox)%xyz(1)) * sideRatio
      iy = dble(sList(sBox)%xyz(2)) * sideRatio
      iz = dble(sList(sBox)%xyz(3)) * sideRatio

      ! check if the LL-box coordinates match a selected replica coordinates
      do ijk = 0, nReplicas-1
        if(groupid(ijk) > 0) then
          k = ijk / nRepLayers2
          j = (ijk - k*nRepLayers2) / nRepLayers1
          i = mod(ijk, nRepLayers1)
          if(i==ix .and. j==iy .and. k==iz) then
            sList(sBox)%group = groupid(ijk)
            nSelectedLLboxesRank = nSelectedLLboxesRank + 1
            exit
          endif
        endif
      enddo
      !
    enddo
  endif  ! if(nSelections > 0)

  !TODO purpose: perform partial reduction for selected boxes in each rank on the host
  !TODO purpose: perform final MPI array reduction on rank 0

  !purpose: Count the number of near- and far-neighbor boxes for each box in S-list in MPI rank

  ! each MPI rank obtains a separate count
  ! the counting applies to floating boxes as well
  ! count the near-neighbor boxes at loLevel only to account for particles
  !$omp do schedule(static,1) &
  !$omp reduction(+:nNearNeighbors) reduction(+:nFarNeighbors) &
  !$omp private(iLevel, nbEdgeMaxID, nbParentEdgeMaxID, ix,iy,iz, pBox) &
  !$omp private(fax,fay,faz, fbx,fby,fbz) &
  !$omp private(nax,nay,naz, nbx,nby,nbz, nNF, nFF, jx,jy,jz)
  do sBox = 1, nBoxesSList-1
    iLevel = sList(sBox)%iLevel
    nbEdgeMaxID = 2**iLevel - 1  ! max box id along the edge at the given level of hierarchy
    nbParentEdgeMaxID = 2**(iLevel - 1) - 1  ! max box id along the edge on parent level
    !
    ix = sList(sBox)%xyz(1)
    iy = sList(sBox)%xyz(2)
    iz = sList(sBox)%xyz(3)

    ! dynamically build the far-field box list
    pBox = sList(sBox)%parent  ! parent box of sBox

    ! far-field limits
    fax = max0(sList(pBox)%xyz(1) - ws, 0) * 2
    fay = max0(sList(pBox)%xyz(2) - ws, 0) * 2
    faz = max0(sList(pBox)%xyz(3) - ws, 0) * 2

    fbx = min0(sList(pBox)%xyz(1) + ws, nbParentEdgeMaxID) * 2 + 1
    fby = min0(sList(pBox)%xyz(2) + ws, nbParentEdgeMaxID) * 2 + 1
    fbz = min0(sList(pBox)%xyz(3) + ws, nbParentEdgeMaxID) * 2 + 1

    ! near-field limits
    nax = max0(ix-ws,0)
    nay = max0(iy-ws,0)
    naz = max0(iz-ws,0)

    nbx = min0(ix+ws,nbEdgeMaxID)
    nby = min0(iy+ws,nbEdgeMaxID)
    nbz = min0(iz+ws,nbEdgeMaxID)

    nNF = 0  ! local count of near-field boxes for sBox
    nFF = 0  ! local count of  far-field boxes for sBox
    do jz = faz, fbz
      do jy = fay, fby
        do jx = fax, fbx
          if(jx.ge.nax .and. jx.le.nbx .and. &
             jy.ge.nay .and. jy.le.nby .and. &
             jz.ge.naz .and. jz.le.nbz ) then
            ! account for near-neighbor boxes at loLevel only
            if(iLevel == loLevel) nNF = nNF + 1 ! near-neighbour includes self-count
          else
            nFF = nFF + 1                       !  far-neighbour
          endif
        enddo
      enddo
    enddo
    nNF = max0(0, nNF - 1)                 ! exclude self-count
    boxRank(sBox*2) = nNF                  ! the number of near-neighbors for sBox
    nNearNeighbors  = nNearNeighbors + nNF ! total number of boxes in the near-neighbor list per rank

    boxRank(sBox*2+1)= nFF                 ! the number of  far-neighbors for sBox
    nFarNeighbors    = nFarNeighbors + nFF ! total number of boxes in the far-neighbor list per rank
  enddo  !sBox

  !$omp end parallel

#ifdef DEBUG_TIMER
  ! print selected box coordinates
  if(iam == prtRank .and. nGroups > 0) then
    do sBox = sBoxLLstart, nBoxesSList-1 ! Loop over boxes of the lowest level
      if(sList(sBox)%group > 0) write(*,'(a,3i6,4x,a,i0)') &
        "dbg: LL-box: ", sList(sBox)%xyz, "is in group: ", sList(sBox)%group
    enddo
  endif

  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif 

  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 1:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  !purpose: Project the min-max boundaries of R-tree at rmLevel to those at loLevel

  ! it must be done outside of the parallel region
  n = 2**(loLevel-rmLevel)
  minRx = minRx * n
  minRy = minRy * n
  minRz = minRz * n

  maxRx = (maxRx + 1) * n - 1
  maxRy = (maxRy + 1) * n - 1
  maxRz = (maxRz + 1) * n - 1

  !purpose: Sequentially numerate the starting and ending indices for near- and far-neighbor lists for S-boxes

  nC = 0
  nF = 0
  do sBox = 1, nBoxesSList-1    ! sBox is a local int*4 index for sList
    sList(sBox)%nearStart = nC + 1
    nC = nC + boxRank(sBox*2)     ! boxRank(sBox*2) holds the number of near-field boxes for sBox
    sList(sBox)%nearEnd   = nC

    sList(sBox)%farStart  = nF + 1
    nF = nF + boxRank(sBox*2+1)   ! boxRank(sBox*2+1) holds the number of far-field boxes for sBox
    sList(sBox)%farEnd    = nF
  enddo 

#if defined (DEBUG) && defined (DETAILED1)
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,i0)') &
      "nBoxesTotal             = ", nBoxesTotal,                 &
      "maxBoxesSList           = ", maxBoxesSList,               &
     ("nBoxesSList (rank X)    = ", sBoxesTopDispl(k), k=0, nRanks-1),  &
      "loLevel boxes           = ", 8_8 ** loLevel,              &
      "ffLevel boxes           = ", 8 ** ffLevel
    !write(*,'(/a)') "Domain decomposition. All done"
  endif
  !call MPI_Finalize(mpiError)
  !error stop
#endif

  ! TODO the actual threshold is unknown
  ! TODO test if we still need it
  if(nAtomsTotal >= 380352) then
    nThreadsLocal = nThreads
  else
    nThreadsLocal = 1  ! needed for !$omp atomic update for small systems
  endif

  ! allocate memory
  memLocal = memLocal + dble(nNearNeighbors)*(8) + dble(nFarNeighbors)*(8) + dble(nThreadsLocal*4*3)
  memMax = dmax1(memMax, (memUsed + memLocal))
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan + nNeighbors", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(                                    &
!          nearNeighbors8(nNearNeighbors),     &
!          farNeighbors8(nFarNeighbors),       &
           nearNeighbors(nNearNeighbors),      &
           farNeighbors(nFarNeighbors),        &
           nSegmentsThread(0:nThreadsLocal-1), &
           nBoxesThread(0:nThreadsLocal),      &  ! need an extra element to simplify renumbering
           nCubesThread(0:nThreadsLocal-1),    &
           stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  !$omp parallel num_threads(nThreadsLocal)

  !purpose: Build the near-field and far-field lists; mark the unique boxes in btypeNS and btypeF

  ! TODO get min (first) max (last) box for each iLevel: use min/max array reduction for int8
  ! create a far-neighbor list; produce a separate list for each MPI rank
  !!$omp do schedule(static,1) &
  !$omp do &
  !$omp private(bBox, iLevel, nbEdge, nbParentEdgeMaxID, n, ix,iy,iz, pBox) &
  !$omp private(fax,fay,faz, fbx,fby,fbz, num1_8, dummy) &
  !$omp private(nax,nay,naz, nbx,nby,nbz, nC, nF, jx,jy,jz) &
  !$omp private(jBox, tjx,tjy,tjz, nbParent, aLevel, pjx,pjy,pjz, boxesStart) &
  !$omp reduction(min:minNx,minNy,minNz) reduction(max:maxNx,maxNy,maxNz)
  do sBox = 1, nBoxesSList-1
    !bBox = sList(sBox)%idx     ! convert S-box to B-box
    bBox = idx(sBox)          ! convert S-box to B-box
    iLevel = sList(sBox)%iLevel
!    nBoxesBefore = (8**iLevel - 1) / 7  ! number of boxes before iLevel
    nbEdge = 2 ** iLevel      ! number of boxes along the edge at the given level of hierarchy
!    nbEdgeMaxID = nbEdge  - 1  ! max box id along the edge at the given level of hierarchy
    nbParentEdgeMaxID = 2**(iLevel - 1) - 1  ! max box id along the edge on parent level
!    n = 8**(iLevel-ffLevel)   ! the number of boxes at iLevel spawned from the ffLevel parent

    ix = sList(sBox)%xyz(1)
    iy = sList(sBox)%xyz(2)
    iz = sList(sBox)%xyz(3)

    ! dynamically build the far-field box list
    pBox = sList(sBox)%parent  ! parent box of sBox

    ! far-field limits
    fax = max0(sList(pBox)%xyz(1) - ws, 0) * 2
    fay = max0(sList(pBox)%xyz(2) - ws, 0) * 2
    faz = max0(sList(pBox)%xyz(3) - ws, 0) * 2

    fbx = min0(sList(pBox)%xyz(1) + ws, nbParentEdgeMaxID) * 2 + 1
    fby = min0(sList(pBox)%xyz(2) + ws, nbParentEdgeMaxID) * 2 + 1
    fbz = min0(sList(pBox)%xyz(3) + ws, nbParentEdgeMaxID) * 2 + 1

    ! near-field limits
    nax = max0(ix-ws,0)
    nay = max0(iy-ws,0)
    naz = max0(iz-ws,0)

    nbx = min0(ix+ws,nbEdge-1)
    nby = min0(iy+ws,nbEdge-1)
    nbz = min0(iz+ws,nbEdge-1)

    ! the local index for far-neighbor boxes MUST BE int*8 (confirmed in a single-node test)
    nC = sList(sBox)%nearStart ! nC is a starting local index of near-boxes for sBox
    nF = sList(sBox)%farStart  ! nF is a starting local index of  far-boxes for sBox
    boxesStart = (8_8 ** iLevel - 1_8) / 7_8  ! starting number for the boxes at this level
    do jz = faz, fbz
      do jy = fay, fby
        do jx = fax, fbx
          ! find jBox based on grid coordinates jx,jy,jz for hierarchy level iLevel
          jBox = boxesStart
          tjx = jx
          tjy = jy
          tjz = jz
          ! count the number of boxes at the increasing granularity
          nbParent = nbEdge  ! number of boxes along the edge at the parent level
          do aLevel = 1, iLevel  ! although this loop appears going down the hierarchy, it actually goes upward
            tjx = mod(tjx,nbParent)
            tjy = mod(tjy,nbParent)
            tjz = mod(tjz,nbParent)
            nbParent = nbParent / 2
            pjx = tjx / nbParent
            pjy = tjy / nbParent
            pjz = tjz / nbParent
            num1_8 = nbParent
            jBox = jBox + (pjz * 4 + pjy * 2 + pjx) * num1_8 * num1_8 * num1_8
          enddo

          ! the above block may be replaced by a call to function
          !jBox = idBox(jx, jy, jz, iLevel)

          if(jx.ge.nax .and. jx.le.nbx .and. &
             jy.ge.nay .and. jy.le.nby .and. &
             jz.ge.naz .and. jz.le.nbz ) then
            if(iLevel == loLevel) then
              ! check the R-NF-group grid boundaries and expand if necessary to accomodate jBox (jx,jy,jz)
              minNx = min0(minNx, jx)
              minNy = min0(minNy, jy)
              minNz = min0(minNz, jz)

              maxNx = max0(maxNx, jx)
              maxNy = max0(maxNy, jy)
              maxNz = max0(maxNz, jz)

              ! need to know the NF boxes at the LL only (those that carry particles)
              ! include S-boxes at LL only (those that carry particles)
              if(bBox /= jBox) then
                !nearNeighbors8(nC) = jBox
                nearNeighbors(nC) = jBox
                nC = nC + 1
              endif
              ! if jBox == bBox then jBox is S-box  (descendent of domain decomposition at ddLevel in R-tree)
              ! if jBox /= bBox then jBox is NF-box (is in near-field for the S-box)
              ! no need to distinguish these two types since we need to access particles that they hold
              btypeNS(jBox) = .true.
            endif
          else
            ! need to know the FF boxes at FF-level and down
            !farNeighbors8(nF) = jBox
            farNeighbors(nF) = jBox
            nF = nF + 1
            btypeF(jBox) = .true.  ! mark the FF-box
          endif
        enddo
      enddo
    enddo
  enddo  !sBox

  ! byte-wise compress & copy btypeNS to boxTypeNS and btypeF to boxTypeF
  ! do it here outside of assignment of btypeNS and btypeF to avoid racing condition
  ! multiple threads in fortran cannot concurrently set (different) bits of the same byte
  ! it is because fortran does not support bit-wise memory pointers
  !$omp do
  do iBox = 0, nBoxesTotal_div8-1, nBits
    boxTypeNS%lArray(iBox/nBits) =       &
      + merge(  1, 0, btypeNS(iBox))     &
      + merge(  2, 0, btypeNS(iBox + 1)) &
      + merge(  4, 0, btypeNS(iBox + 2)) &
      + merge(  8, 0, btypeNS(iBox + 3)) &
      + merge( 16, 0, btypeNS(iBox + 4)) &
      + merge( 32, 0, btypeNS(iBox + 5)) &
      + merge( 64, 0, btypeNS(iBox + 6)) &
      + merge(128, 0, btypeNS(iBox + 7))
    boxTypeF%lArray(iBox/nBits) =        &
      + merge(  1, 0, btypeF(iBox))      &
      + merge(  2, 0, btypeF(iBox  + 1)) &
      + merge(  4, 0, btypeF(iBox  + 2)) &
      + merge(  8, 0, btypeF(iBox  + 3)) &
      + merge( 16, 0, btypeF(iBox  + 4)) &
      + merge( 32, 0, btypeF(iBox  + 5)) &
      + merge( 64, 0, btypeF(iBox  + 6)) &
      + merge(128, 0, btypeF(iBox  + 7))
  enddo

  !! compare btypeNS to boxTypeNS and and btypeF to boxTypeF
  !!$omp do
  !do iBox = 0, nBoxesTotal-1
  !  if(btypeNS(iBox) .neqv. btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits))) then
  !    write(*,'(a,i0)') "nBoxesTotal = " , nBoxesTotal
  !    write(*,'(a,i0,4x,l4,2x,l4,8x,i0)') "Error for iBox = ", iBox, &
  !       btypeNS(iBox), btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)), iam
  !    stop "Error in boxTypeNS"
  !  endif
  !  if(btypeF(iBox) .neqv. btest(boxTypeF%lArray(iBox/nBits), mod(iBox,nBits))) then
  !    write(*,'(a,i0)') "nBoxesTotal = " , nBoxesTotal
  !    write(*,'(a,i0,4x,l4,2x,l4,8x,i0)') "Error for iBox = ", iBox, &
  !       btypeF(iBox), btest(boxTypeF%lArray(iBox/nBits), mod(iBox,nBits)), iam
  !    stop "Error in boxTypeF"
  !  endif
  !enddo

  ! purpose: sum-up S + N + R + F box types used in the rank (the S + N + R are from LL only; carry particles)

  ! numBoxesSNRF is the total number of S + N + R + F boxes per MPI rank
  ! include S + N + R boxes at loLevel only (those tht carry particles); ignore boxes at other levels
  ! nBoxesThread() is the number of S + N + R + F boxes per OpenMP thread on the MPI rank
  ! the data distribution over OpenMP threads in nBoxesThread() corresponds to the thread-iterator(iBox) 
  !   assignment taking place in the do-loop over limits: iBox = nBoxesBefore_ffLevel, nBoxesTotal-1
  ! TODO should limit the use of nBoxesThread to a single parallel level

  numBoxesSNRF = 0
  nBoxesThread(0) = 0
  nBoxesThread(omp_get_thread_num()+1) = 0
  !$omp do reduction(+:numBoxesSNRF)
  do iBox = nBoxesBefore_ffLevel, nBoxesTotal-1
    if((iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd) &
!            .or. btypeNS(iBox) &
            .or. btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
!            .or. btypeF(iBox)) then
            .or. btest(boxTypeF%lArray(iBox/nBits), mod(iBox,nBits))) then
      numBoxesSNRF = numBoxesSNRF + 1
      nBoxesThread(omp_get_thread_num()+1) = numBoxesSNRF  ! the number of boxes per thread shifted in thread number
    endif
  enddo

  ! test for potential int4 overflow
  if(numBoxesSNRF <= 0) then
    write(*,'(2(a,i0))') "Error: numBoxesSNRF = ", numBoxesSNRF, " for rank ", iam
    call exitPrg(100)
  endif

  ! S-box is a box from S-tree constructed by domain decomposition of ddLevel for MPI World
  !   The S-tree spans from ddLevel to loLevel
  !   Each rank in MPI World carries a unique S-tree
  !   S-box carries particles when it is at loLevel to compute interaction with near-field boxes
  ! N-box is a loLevel box (carries particles), which is in near-field position to a loLevel S-box
  ! R-box is a loLevel box (carries particles) from R-tree in RMG sub-communicator (rmRank)
  ! F-box is a box positioned in far-field to S-boxes ranging from ffLevel to loLevel
  ! E-box is a non-typical F-box located in between rmLevel-1 and ffLevel when rmLevel /= ffLevel

  ! The actual difference between E-boxes (non-typical F-boxes) and other (typical) F-boxes is
  !   that particles necessary for computation of multipole expansion of E-boxes reside in two
  !   or more RMG ranks, whereas for F-boxes, all particles needed for computation of their
  !   multipole expansions reside in the same rank.

  !$omp end parallel

  memLocal = memLocal - dble(nBoxesTotal_div8 * 2_8)  ! btypeNS and btypeF are no longer needed
  deallocate(btypeNS, btypeF)

  ! purpose: replace global index of S + N + R + F boxes with local box index in MPI rank to save memory

  ! The MPI rank needs to know only its S + N + R + F boxes out of the total (huge) pool of the boxes
  ! Presently, nBoxesThread() holds the number of S + N + R + F boxes per OpenMP thread 
  !   shifted (+1) in thread number in the MPI rank
  ! Need to consecutively sum up the number of boxes per thread accross the threads
  ! Each thread obtains a starting local box index for the boxes to be handled by that thread
  ! This is necessary to renumber S + N + R + F boxes from 0 to numBoxesSNRF-1

  ! Presently, the summation must be done outside of the parallel clause
  do i = 2, nThreadsLocal
    ! remember that the sums are shifted to the right in thread index
    nBoxesThread(i) = nBoxesThread(i) + nBoxesThread(i-1)
  enddo

#if defined(DEBUG_TIMER) && defined(MPI)
  numBoxesSNRFMin = numBoxesSNRF
  numBoxesSNRFMax = numBoxesSNRF
  numBoxesSNRFAve = numBoxesSNRF
  if(debugEachMPIrank) then
    if(nSplitRanks == nodeSize) then
      call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    else
      call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
    endif
  else
    num1_8 = numBoxesSNRF
    if(nSplitRanks == nodeSize) then
      call MPI_Reduce(num1_8, numBoxesSNRFMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(num1_8, numBoxesSNRFMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(num1_8, numBoxesSNRFAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    else
      call MPI_Reduce(num1_8, numBoxesSNRFMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(num1_8, numBoxesSNRFMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(num1_8, numBoxesSNRFAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_PGROUP, mpiError)
    endif
    numBoxesSNRFAve = numBoxesSNRFAve / nRanks
  endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(2(a,i12))')  "dbg: numBoxesSNRF for    prtRank: ", numBoxesSNRF,    "   rank: ", iam
    write(*,'(3(a,i12))')  "dbg: numBoxesSNRF per rank   min: ", numBoxesSNRFMin, &
                "    max: ", numBoxesSNRFMax,  "    ave: ", numBoxesSNRFAve
  endif
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
  !call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 2:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  ! nAtomsSNR(0:numBoxesSNRF) holds the number of particles in S-, N-, and R-boxes.
  ! nAtomsNS(0:numBoxesSNRF)  holds the number of particles in N- and S-boxes. It is a copy of nAtomsSNR.
  !                              Its purpose is to reorder particles in S-boxes.
  ! nAtomsR(0:numBoxesSNRF)   holds the number of particles in R-boxes at loLevel. It is a copy of nAtomsSNR.
  !                              Its purpose is to reorder particles in R-boxes for computation of regular
  !                              multipoles of the boxes.

  ! allocate memory for local (renumbered) boxes on the rank
  memLocal = memLocal + dble(numBoxesSNRF)*8 + dble(numBoxesSNRF)*4 + dble(numBoxesSNRF)*8*3
  memMax = dmax1(memMax, memUsed+memLocal)
  if(iam == prtRank .or. debugEachMPIrank) &
          call MemPrint("Plan + getGlobalFromLocal", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(getGlobalFromLocal(0:numBoxesSNRF), boxesRank(0:numBoxesSNRF), &
          nAtomsR(0:numBoxesSNRF), nAtomsNS(0:numBoxesSNRF), nAtomsSNR(0:numBoxesSNRF), stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  ! S N R boxes carry particles at loLevel
  getGlobalFromLocal = 0  ! convert the local (for the MPI rank) box index to the global box index
  nAtomsSNR    = 0  ! zero the counter of number of particles in the boxes
  numBoxesSNRF = numBoxesSNRF + 1  ! add the master box
  boxesRank = 0

  !$omp parallel num_threads(nThreadsLocal)

  ! purpose: renumber R + N + S + F boxes handled in the MPI rank to produce a local index

  ! start from ffLevel boxes and go down to loLevel boxes
  ! renumbering involves all selected boxes
  ! this is a loop over the total number of boxes in the system

  ! Each thread works on a separate consecutive (based on global index) range of R + N + S + F boxes
  ! Presently, nBoxesThread holds a starting value of the local box index for each thread, so
  !   each thread can renumber its boxes independently from the work done in other threads
  ! OpenMP splits the do-loop index range into nThreadsLocal continous segments
  ! Each OpenMP thread handles a separate index segment

  !$omp do
  do iBox = nBoxesBefore_ffLevel, nBoxesTotal-1   ! go through all boxes; do not need boxes above ffLevel
    if((iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd) &
!            .or. btypeNS(iBox) &
            .or. btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
!            .or. btypeF(iBox)) then
            .or. btest(boxTypeF%lArray(iBox/nBits), mod(iBox,nBits))) then
      ! this iBox is S N R F type box in the MPI rank
      ! increment the local box index in each thread
      nBoxesThread(omp_get_thread_num()) = nBoxesThread(omp_get_thread_num()) + 1
      ! each selected box gets a sequential index local to the MPI rank
#ifdef USE_GETLOCAL_FROMGLOBAL
      getLocalFromGlobal(iBox) = nBoxesThread(omp_get_thread_num()) ! provides global index -> local index conversion
#endif
      getGlobalFromLocal(nBoxesThread(omp_get_thread_num())) = iBox ! provides local index -> global index conversion 
    endif
  enddo

  !purpose: Assign rank id to FF boxes in boxesRank and count remote FF boxes. Do not need btypeF after that.

  ! Arrange the retrieval of FF boxes from the remote RMG ranks in segments
  ! FF boxes are those positioned in far-field to S-boxes
  ! consecutive FF boxes residing on the same remore rank (in my RMG) form a segment
  ! it is more performant to retrieve the remote boxes in segments rather than individually
  ! compute the number of remote FF boxes and the number of their segments per thread
  nBoxesThread(omp_get_thread_num()) = 0
  nSegmentsThread(omp_get_thread_num()) = 0
  jBox = 0
  jRank = 0
  nRemoteFFBoxes = 0 ! number of unique remote FF boxes needed for the rank
  !
  !purpose: Assign rank id to FF boxes
  !
  iLevelOld = 0
  j = 0
  !$omp do firstprivate(iLevelOld) private(iBox, iLevel, nBoxesBefore, n, fBox) &
  !$omp reduction(+:nRemoteFFBoxes) firstprivate(jBox,jRank) private(iRank)
!  do iBox = nBoxesBefore_ffLevel, nBoxesTotal-1
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)  ! convert the local box index to the global box index
!    if(btypeF(iBox)) then
    if(btest(boxTypeF%lArray(iBox/nBits), mod(iBox,nBits))) then
      ! this is an F-box (E-box) which brings far-field potential to the boxes of S-tree
      ! F-box (E-box) can be local or remote

      ! get iLevel from iBox using the inequality iLevel <= log(iBox*7+1)/log8 < iLevel+1
      iLevel = floor(dlog(dble(iBox*7+1)) * invLog8)

      if(iLevel > iLevelOld) then
        iLevelOld = iLevel
        nBoxesBefore = (8_8 ** iLevel - 1_8) / 7_8  ! the total number of boxes before iLevel in B-list
        ! n is the number of boxes at iLevel spawned from a single rmLevel parent
        ! n = 1 for iLevel == (ffLevel ... rmLevel)
        n = 8**(max0(0, iLevel-rmLevel))
      endif

      ! determine which rank in RMG owns the F-box, iBox
      ! do not need to know the MPI World rank id since the communication will use RMG sub-communicator
      ! the RMG rank id of this F-box is the same for all MPI World ranks
      ! each RMG sub-communicator has exactly the same distribution of all boxes of R-list
      ! each R-list is replicated by RMG sub-communicators

      ! the F-box can either be local or remote to the rank iam "==" rmRank
      ! the F-box can be either a typical (F-box) or an extra multipole box (E-box)
      ! the only presently relevant difference between F-box and E-box is that determination of their
      !   RMG rank id requires a different technique

      if(iLevel >= rmLevel) then
        ! this is an F-box
        ! F-box is always someone's R-box in RMG
        ! Find which rank in RMG has the F-box as R-box (which has multipole moments)
        ! (iBox - nBoxesBefore) is a local index of iBox at iLevel
        ! fBox is a local index of the box at rmLevel where R-tree begins
        fBox = (iBox - nBoxesBefore) / n   ! for iLevel == rmLevel we have n = 1
        ! save the RMG rank (1...rmSize) (offset by 1) that owns the FF box needed for the global rank iam
        boxesRank(a) = 1 + merge(fBox/minTopRtreeSize, &
                                 minTopRtreeRanks+(fBox-minTopRtreeGroup)/(minTopRtreeSize+1), &
                                 fBox < minTopRtreeGroup)
      else
        ! this is an E-box (extra multipole box) due to ffLevel <= iLevel < rmLevel
        ! the RMG rank id of the E-box is the local box id of its parent (or self) at ffLevel
        ! (iBox - nBoxesBefore) is the local index of iBox at iLevel; maximum one box at ffLevel per rmRank
        ! (iBox - nBoxesBefore)/8**(iLevel-ffLevel) is the local index of a parent of iBox at ffLevel
        boxesRank(a) = 1 + (iBox - nBoxesBefore)/8**(iLevel-ffLevel)  ! +1 to shift the rank id (1...rmSize)
      endif

      iRank = boxesRank(a)   ! R and F boxes have MPI rank id

      ! if the F-box resides in a remote MPI rank, arrange its retrieval
      ! the local F-box is already available to the rmRank
      ! create/update segments of consecutive FF boxes residing on remote RMG ranks
      ! count the number (segment length) of consecutive FF boxes in the remote segment
      if(iRank /= rmRank1) then
        ! iBox is a remote FF box that resides on a remote rank iRank in my RMG
        nRemoteFFBoxes = nRemoteFFBoxes + 1
        nBoxesThread(omp_get_thread_num()) = nRemoteFFBoxes  ! the number of FF boxes in a segment per thread
        ! each OpenMP thread handles a consecutive segment of the array
        ! for this to function efficiently, OpenMP parallelization of the loop do a = 1, numBoxesSNRF-1
        ! should be static with chunk size = numBoxesSNRF / nThreadsLocal
        !   so that each thread handles a maximally long continous portion of box indices, a in the do-loop
        ! TODO consider experimenting with chunk size
        ! the arrangements of remote FF boxes in fewer segments will improve the MPI retrieval time
        ! TODO track the number of segments across all MPI ranks
        if(iBox == jBox+1 .and. iRank == jRank) then
          ! no gaps, continue the current segment
          jBox = iBox  ! jBox is the latest box in the current segment
        else
          ! start a new segment
          jBox = iBox    ! the new segment starts from iBox
          jRank = iRank  ! the new segment has the value iRank
          ! the number of segments per thread
          nSegmentsThread(omp_get_thread_num()) = nSegmentsThread(omp_get_thread_num()) + 1
        endif
      endif
    else
      ! if the box is not F-type, then it is local to rmRank
      ! this is a box of type S, N, or R (all are local to rmRank)
      boxesRank(a) = rmRank1  ! the iBox resides in rmRank (saved as rmRank1)
    endif
  enddo

  !purpose: Count the number of particles held in SNR boxes; use btypeNS to narrow down the search

  ! count the number of particles in the boxes
  ! do not save the replicated particle coordinates
  ! TODO use inverted coordinates to make the enclosed loop more efficient
#ifdef DO_DIRECT_SUM
  !$omp do private(i,j,k, translation, iAtom, ix,iy,iz, nC)
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the input box
    translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the input box
    translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the input box
    !do iAtom = 1, nAtoms
    do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
      !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
      ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
      iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
      iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)
      nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and R-list
      !$omp atomic update
      nAtomsCube(nC) = nAtomsCube(nC) + 1
    enddo
  enddo
#else

#if defined(DEBUG_TIMER) && defined(MPI)
  ! TODO add support for multiple datasets
  ! Count the number of replicas needed per rank in the debug mode
  ! Determine how many replicas overlap with R- and NS-projected boxes in the rank
  ! This is an upper-bound estimate since not all overlapping replicas will be needed
  nReplicasRank = 0
  !$omp do private(i,j,k, xloCell, yloCell, zloCell, xupCell, yupCell, zupCell) reduction(+:nReplicasRank)
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with N or R spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      nReplicasRank = nReplicasRank + 1_8
    endif
  enddo
  !$omp end parallel

  ! test for potential int32 overflow
  if(nReplicasRank <= 0) then
    write(*,'(2(a,i0))') "Error: nReplicasRank = ", nReplicasRank, " for rank ", iam
    call exitPrg(100)
  endif

  nReplicasRankMin = nReplicasRank
  nReplicasRankMax = nReplicasRank
  nReplicasRankAve = nReplicasRank
  if(debugEachMPIrank) then
    if(nSplitRanks == nodeSize) then
      call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    else
      call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
    endif
  else
    if(nSplitRanks == nodeSize) then
      call MPI_Reduce(nReplicasRank, nReplicasRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nReplicasRank, nReplicasRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nReplicasRank, nReplicasRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    else
      call MPI_Reduce(nReplicasRank, nReplicasRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nReplicasRank, nReplicasRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nReplicasRank, nReplicasRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_PGROUP, mpiError)
    endif
    nReplicasRankAve = nReplicasRankAve / nRanks
  endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(3(a,i12))')  "dbg: nReplicas for rank  prtRank: ", nReplicasRank,    &
                "   rank: ", iam,          "  total: ", nReplicas
    write(*,'(3(a,i12))')  "dbg: nReplicas per rank      min: ", nReplicasRankMin, &
                "    max: ", nReplicasRankMax, "    ave: ", nReplicasRankAve
  endif
  !$omp parallel num_threads(nThreadsLocal)
#endif ! DEBUG_TIMER

#ifdef FLEX_SIDELENGTH
  ! FLEX_SIDELENGTH allows mod(PrimitiveEdge,CubeEdge) > 0
  ! with this inequality, threads should not share adjacent replicas since it creates a racing condition
  ! threads processing adjacent replicas will be writing to the same LL boxes
  ! to prevent it, threads should not be allowed to work on adjacent replicas when mod(PrimitiveEdge,CubeEdge) > 0
  ! as a quick fix, spacing the index k and assigning it to threads guarantees threads writing to different LL boxes
  ! it requires to deal with two loops, one for even and the other for odd values of index k
  ! to keep a single loop instead of two, one would need to generate ijk with non-adjacent indices for the threads
  !$omp do private(i,j, ij,ijk, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, idx_4)
  do k = 0, nRepLayers, 2
  do ij = 0, nRepLayers2 - 1  ! if this is uncommented, make sure ijk is declared as a private variable
!  do ijk = 0, nReplicas - 1
!    k = ijk / nRepLayers2     ! if this is uncommented, make sure k is declared as a private variable
!    j = (ijk - k*nRepLayers2) / nRepLayers1
    ijk = k * nRepLayers2 + ij
    j = ij / nRepLayers1
    i = mod(ij, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with N or R spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the primitive cell
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the primitive cell
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the primitive cell
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif  ! MPI
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          ! each OpenMP thread handles a unique replica
          ! each replica consists of a unique set of boxes so the OpenMP atomic directive is unnecessary
          !!$omp atomic update
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          ! need this to prevent threads writing to the same LL-box
          !!$omp atomic update 
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! count the number of paricles in iBox
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo ! ijk
  enddo ! k

  ! make all threads wait here until the above do-loop over the even values of index k is complete
  !$omp barrier

  !$omp do private(i,j, ij,ijk, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, idx_4)
  do k = 1, nRepLayers, 2
  do ij = 0, nRepLayers2 - 1
!  do ijk = 0, nReplicas - 1
!    k = ijk / nRepLayers2
!    j = (ijk - k*nRepLayers2) / nRepLayers1
    ijk = k * nRepLayers2 + ij
    j = ij / nRepLayers1
    i = mod(ij, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with N or R spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the primitive cell
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the primitive cell
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the primitive cell
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif  ! MPI
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          ! each OpenMP thread handles a unique replica
          ! each replica consists of a unique set of boxes so the OpenMP atomic directive is unnecessary
          !!$omp atomic update
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          ! need this to prevent threads writing to the same LL-box
          !!$omp atomic update
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! count the number of paricles in iBox
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo ! ijk
  enddo ! k

#else   ! not FLEX_SIDELENGTH
  ! no racing condition happens when mod(PrimitiveEdge,CubeEdge) == 0

  !$omp do private(i,j,k, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, idx_4)
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with N or R spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the primitive cell
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the primitive cell
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the primitive cell
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          ! each OpenMP thread handles a unique replica
          ! each replica consists of a unique set of boxes so the OpenMP atomic directive is unnecessary
          !!$omp atomic update
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          !!$omp atomic update
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! count the number of paricles in iBox
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo ! ijk

#endif  ! test for FLEX_SIDELENGTH

#endif  ! not DO_DIRECT_SUM

#ifdef DEBUG_TIMER
  ! Count the number of SNR particles in the rank
  !$omp end parallel

!  if(iam == prtRank) then
!    write(*,'(/a)') "Box xyz:"
!    do a = 1, 30
!      iBox = a
!      iLevel = floor(dlog(dble(iBox*7+1)) * invLog8)
!
!      jBox = iBox - (8 ** iLevel - 1) / 7  ! jBox is a local index of iBox for ffLevel
!      iz = 0
!      iy = 0
!      ix = 0
!      ! the number of primitive boxes in the parent box along the edge
!      nbParentEdge = 2 ** (iLevel - 1)
!      do aLevel = 1, iLevel
!        ! the number of primitive (lowest-level) boxes in the parent box
!        nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge
!
!        pBox = jBox / nbParentCube        ! the parent box index that holds jBox (iBox)
!        tjz =  pBox / 4
!        tjy = (pBox - tjz * 4) / 2
!        tjx =  pBox - tjz * 4 - tjy * 2
!        iz = iz + tjz * nbParentEdge
!        iy = iy + tjy * nbParentEdge
!        ix = ix + tjx * nbParentEdge
!
!        jBox = mod(jBox,nbParentCube)
!        nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
!      enddo
!
!      write(*,'(8i8,2l8)') iam, a, iBox, iLevel, ix,iy,iz
!    enddo
!    write(*,'(a)') ""
!
!    write(*,'(/a)') "S-list:"
!    do a = 1, 20
!      iBox = idx(a)          ! convert S-box to B-box
!      iLevel = sList(a)%iLevel
!      write(*,'(8i8,2l8)') iam, a, iBox, iLevel, sList(a)%xyz(1), sList(a)%xyz(2), sList(a)%xyz(3)
!    enddo
!
!    write(*,'(/a)') "S N R F:"
!    do a = 1, 50
!      iBox = getGlobalFromLocal(a)
!      iLevel = floor(dlog(dble(iBox*7+1)) * invLog8)
!
!      jBox = iBox - (8 ** iLevel - 1) / 7  ! jBox is a local index of iBox for ffLevel
!      iz = 0
!      iy = 0
!      ix = 0
!      ! the number of primitive boxes in the parent box along the edge
!      nbParentEdge = 2 ** (iLevel - 1)
!      do aLevel = 1, iLevel
!        ! the number of primitive (lowest-level) boxes in the parent box
!        nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge
!
!        pBox = jBox / nbParentCube        ! the parent box index that holds jBox (iBox)
!        tjz =  pBox / 4
!        tjy = (pBox - tjz * 4) / 2
!        tjx =  pBox - tjz * 4 - tjy * 2
!        iz = iz + tjz * nbParentEdge
!        iy = iy + tjy * nbParentEdge
!        ix = ix + tjx * nbParentEdge
!
!        jBox = mod(jBox,nbParentCube)
!        nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
!      enddo
!
!      write(*,'(8i8,2l8)') iam, a, iBox, iLevel, ix,iy,iz, boxesRank(a), btypeF(iBox), btypeNS(iBox)
!    enddo
!
!    write(*,'(/a)') "End of debug printout"
!  endif

  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  nReplicasNeeded = max0(nReplicasNeeded, nReplicasUsed)
  nReplicasUsed = 0
  nReplicasUsedRank = nReplicasNeeded
  nReplicasUsedRankMin = nReplicasUsedRank
  nReplicasUsedRankMax = nReplicasUsedRank
  nReplicasUsedRankAve = nReplicasUsedRank
  if(debugEachMPIrank) then
    if(nSplitRanks == nodeSize) then
      call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    else
      call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
    endif
  else
    if(nSplitRanks == nodeSize) then
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    else
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nReplicasUsedRank, nReplicasUsedRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_PGROUP, mpiError)
    endif
    nReplicasUsedRankAve = nReplicasUsedRankAve / nRanks
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(3(a,i12))')  "dbg: nReplicasUsed  for  prtRank: ", nReplicasUsedRank,    &
                "   rank: ", iam,          " out of: ", nReplicasRank
    write(*,'(3(a,i12))')  "dbg: nReplicasUsed per rank  min: ", nReplicasUsedRankMin, &
                "    max: ", nReplicasUsedRankMax, "    ave: ", nReplicasUsedRankAve
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 3:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif

  call TimerInit(3)
  !$omp parallel num_threads(nThreadsLocal)
#endif
 
  !purpose: Count the total number of particles (nParticlesRMRank) in R boxes assigned to rmRank in RMG
  !         R-boxes appear at loLevel of the R-tree, which starts from rmLevel

#ifdef  DO_DIRECT_SUM
  !$omp do private(jBox,xyz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz) reduction(+:nParticlesRMRank)
  do iBox = rBoxBottomStart, rBoxBottomEnd

    ! get ix, iy, iz from iBox
    jBox = iBox - nBoxesBefore_loLevel   ! get the local index of iBox for loLevel
    xyz = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = nbParentEdge_init
    do aLevel = 1, loLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / (nbParentEdge * nbParentEdge * nbParentEdge)   ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      xyz(1) = xyz(1) + tjx * nbParentEdge
      xyz(2) = xyz(2) + tjy * nbParentEdge
      xyz(3) = xyz(3) + tjz * nbParentEdge

      jBox = mod(jBox, nbParentEdge * nbParentEdge * nbParentEdge)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    !nC = box(iBox)%xyz(1) &
    !   + box(iBox)%xyz(2) * nChildEdge &
    !   + box(iBox)%xyz(3) * nChildXY
    !nParticlesRMRank = nParticlesRMRank + nAtomsCube(nC)
    nParticlesRMRank = nParticlesRMRank + nAtomsCube(xyz(1) + xyz(2)*nChildEdge + xyz(3)*nChildXY)
  enddo
#else
  !$omp do reduction(+:nParticlesRMRank)
  do iBox = rBoxBottomStart, rBoxBottomEnd
#ifdef USE_GETLOCAL_FROMGLOBAL
    nParticlesRMRank = nParticlesRMRank + nAtomsSNR(getLocalFromGlobal(iBox))
#else
    nParticlesRMRank = nParticlesRMRank + nAtomsSNR(findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1)
#endif
  enddo
#endif

  !purpose: Count the total number of particles (nParticlesBuffer) in the tmp buffer: S + N + R boxes

#ifdef  DO_DIRECT_SUM
  !$omp do private(iBox,jBox,xyz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz) reduction(+:nParticlesBuffer)
  !do iBox = bBoxesStart, nBoxesTotal-1
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)
!    if(btypeNS(iBox) &
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
            .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
      ! get ix, iy, iz from iBox
      jBox = iBox - nBoxesBefore_loLevel   ! get the local index of iBox for iLevel
      xyz = 0
      ! the number of primitive boxes in the parent box along the edge
      nbParentEdge = nbParentEdge_init
      do aLevel = 1, loLevel
        ! the number of primitive (lowest-level) boxes in the parent box
        nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge
      
        pBox = jBox / (nbParentEdge * nbParentEdge * nbParentEdge)   ! the parent box index that holds jBox (iBox)
        tjz =  pBox / 4
        tjy = (pBox - tjz * 4) / 2
        tjx =  pBox - tjz * 4 - tjy * 2
        xyz(1) = xyz(1) + tjx * nbParentEdge
        xyz(2) = xyz(2) + tjy * nbParentEdge
        xyz(3) = xyz(3) + tjz * nbParentEdge
      
        jBox = mod(jBox, nbParentEdge * nbParentEdge * nbParentEdge)
        nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
      enddo
      
      nParticlesBuffer = nParticlesBuffer + &
                nAtomsCube(xyz(1) + xyz(2)*nChildEdge + xyz(3)*nChildXY)
    endif
  enddo
#else
  !$omp do private(iBox) reduction(+:nParticlesBuffer)
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)
!    if(btypeNS(iBox) &
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
            .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
      nParticlesBuffer = nParticlesBuffer + nAtomsSNR(a)
    endif
  enddo
#endif

  !purpose: Count the number of S + N + R + F boxes per OpenMP thread; skip the master box

  ! TODO Question: Do we need S-boxes above loLevel ?
  ! TODO If yes, they can be added as local F-boxes having the rank id rmRank1
  nCubesThread(omp_get_thread_num()) = 0
  !$omp do
  do a = 1, numBoxesSNRF-1
    nCubesThread(omp_get_thread_num()) = nCubesThread(omp_get_thread_num()) + 1
  enddo

  !$omp end parallel

  !purpose: Make a backup copy of the array that holds the number of particles in S N R boxes

#ifdef DO_DIRECT_SUM
  boxRank = nAtomsCube
  backupBuffer = nAtomsCube
#else
  nAtomsR  = nAtomsSNR
  nAtomsNS = nAtomsSNR
#endif

  !$omp parallel sections
  !!$omp section
  !! atom displacements; starting position of the first atom in the rank
  !  do iRank = 1, nRanks-1
  !    dispAtoms(iRank) = dispAtoms(iRank-1) + rankAtoms(iRank-1)
  !  enddo
  !$omp section

  !purpose: Count the number of segments and sum up the starting address per OpenMP thread

    nBoxSegments = 0
    j = 0  ! old value
    k = 0  ! running sum
    m = 0  ! running sum
    do i = 0, nThreadsLocal-1
      nBoxSegments = nBoxSegments + nSegmentsThread(i)
      j = nSegmentsThread(i)  ! save the value
      nSegmentsThread(i) = k  ! write the sum of previous elements
      k = k + j           ! accumulate the sum of previous elements
      j = nBoxesThread(i) ! save the value
      nBoxesThread(i) = m ! write the sum of previous elements
      m = m + j           ! accumulate the sum of previous elements
    enddo
  !$omp section

  !purpose: Determine the starting index of R+S+N+F boxes for each thread

    !kThread = 0  ! running sum; start from the master box
    kThread = 1  ! running sum; skip the master box
    do iThread = 0, nThreadsLocal-1
      jThread = nCubesThread(iThread)  ! save the value
      nCubesThread(iThread) = kThread  ! write the sum of previous elements
      kThread = kThread + jThread      ! accumulate the sum of previous elements
    enddo
  !$omp end parallel sections

#if defined(DEBUG_TIMER) && defined(MPI)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
  nParticlesRMRankMin = nParticlesRMRank
  nParticlesRMRankMax = nParticlesRMRank
  nParticlesRMRankAve = nParticlesRMRank
  nParticlesBufferMin = nParticlesBuffer
  nParticlesBufferMax = nParticlesBuffer
  nParticlesBufferAve = nParticlesBuffer
  if(debugEachMPIrank) then
    if(nSplitRanks == nodeSize) then
      call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    else
      call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
    endif
   !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
   !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
  else
    if(nSplitRanks == nodeSize) then
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    else
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nParticlesRMRank, nParticlesRMRankAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_PGROUP, mpiError)
      call MPI_Reduce(nParticlesBuffer, nParticlesBufferAve, 1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_PGROUP, mpiError)
    endif
    nParticlesRMRankAve = nParticlesRMRankAve / nRanks
    nParticlesBufferAve = nParticlesBufferAve / nRanks
  endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(2(a,i12))')    "dbg: nParticlesRMRank    prtRank: ", nParticlesRMRank,    "   rank: ", iam
    write(*,'(3(a,i12))')    "dbg: nParticlesRMRank        min: ", nParticlesRMRankMin, &
               "    max: ", nParticlesRMRankMax, "    ave: ", nParticlesRMRankAve
    write(*,'(2(a,i12))')    "dbg: nParticlesBuffer    prtRank: ", nParticlesBuffer,    "   rank: ", iam
    write(*,'(3(a,i12))')    "dbg: nParticlesBuffer        min: ", nParticlesBufferMin, &
               "    max: ", nParticlesBufferMax, "    ave: ", nParticlesBufferAve
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 4:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  ! allocate memory
  memLocal = memLocal + dble(nBoxSegments)*(4+4) &
                      + dble(nBoxSegments)*(4+4) &
                      + dble(numBoxesSNRF)*44 
  memMax = dmax1(memMax, (memUsed + memLocal))
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan + Segments", memLocal, memUsed+memLocal, memMax, planSize)
  allocate(segmentRank(nBoxSegments), segmentNBoxes(nBoxSegments), &
           segmentSBox(nBoxSegments), segmentGaps(nBoxSegments), &
           rList(0:numBoxesSNRF-1), stat=allocStat)
  if(allocStat /= 0) then
    if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan", memLocal, memUsed+memLocal, memMax, planSize)
    call exitPrg(2)
  endif
  segmentRank   = 0
  segmentSBox   = 0
  segmentNBoxes = 0 
  segmentGaps   = 0
  ! initialize the master box
  rList(0)%xyz(1)     =  0
  rList(0)%xyz(2)     =  0
  rList(0)%xyz(3)     =  0
  rList(0)%atomStart  =  1
  rList(0)%atomEnd    =  nAtomsTotal

  !$omp parallel num_threads(nThreadsLocal)

  !purpose: Initialize boxes in the expanded R-list: S + N + R + F
  !purpose: Identify segments of remote boxes (stored in other ranks in RMG) needed for FF computation
  ! TODO find a better name for "expanded R-list"
  
  j     = 0
  k     = 0
  jBox  = 0
  jRank = 0
  iLevelOld = 0
  !$omp do firstprivate(j,k,jBox,jRank,iLevelOld) private(num1_8) &
  !$omp private(iBox,iRank,iLevel,eLevel,nBoxesBefore,i,n,m,bBefore,bStart,bEnd,expStart,expRemote) &
  !$omp private(aBox,xyz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz,nbParentEdge_init,nxb)
!  do iBox = nBoxesBefore_ffLevel, nBoxesTotal-1
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)
    iRank = boxesRank(a) - 1 ! iRank is the RMG rank id; adjust for the shift in rank id
    ! TODO consider marking only remote boxes to avoid checking iRank /= rmRank1
    ! get iLevel from iBox using the inequality iLevel <= log(iBox*7+1)/log8 < iLevel+1
    iLevel = floor(dlog(dble(iBox*7+1)) * invLog8)  ! ffLevel <= iLevel <= loLevel

    if(iLevel > iLevelOld) then
      iLevelOld = iLevel
      nBoxesBefore = (8_8 ** iLevel - 1_8) / 7_8  ! the total number of boxes before iLevel in B-tree
      ! n is the number of boxes at iLevel spawned from a single rmLevel parent
      n = 8**(max0(0, iLevel-rmLevel))    ! n = 1 for iLevel <= rmLevel
      num1_8 = n
      ! the total number of boxes before iLevel-rmLevel in R-tree
      m = (n - 1) / 7                     ! m = 0 if iLevel <= rmLevel
      nbParentEdge_init = 2 ** (iLevel - 1)
    endif

    if(iRank /= rmRank) then
      ! this is a remote F-box (not in rmRank); it is an R-list box in a remote rank
      ! the R-list boxes occupy ffLevel..loLevel
      ! each OpenMP thread handles a consecutive segment of the array
      if(iBox == jBox+1 .and. iRank == jRank) then
        ! no gaps, continue the current segment
        jBox = iBox  ! jBox is the latest box in the current segment
        segmentNBoxes(i) = segmentNBoxes(i) + 1  ! number of boxes in the segment
      else
        ! start a new segment
        jBox = iBox    ! the new segment starts from iBox
        jRank = iRank  ! the new segment has the value iRank
        j = j + 1  ! advance the number of isegments in the thread
        i = nSegmentsThread(omp_get_thread_num()) + j  ! address to write the segment
        segmentRank(i)   = iRank
        segmentSBox(i)   = a
        segmentNBoxes(i) = 1
        segmentGaps(i)   = k  ! the number of unneeded boxes before the segment starts
        k = 0  ! reset the gap-width counter

        ! Determine the local address of iBox in R-list in the remote rank
        ! the iBox resides at iLevel-rmLevel down from the top of R-tree, which is rmLevel
        ! based on RMG domain decomposition at rmlevel,
        !   determine the number of boxes in the R-tree down iLevel-rmLevel
        !     from the top (rmLevel) of R-tree on the remote iRank
        ! m is number of boxes expanded from a box at rmLevel down by iLevel-rmLevel levels
        ! mod(iRank,rmSize) is an RMG rank id for iRank
        ! bBefore is a local (within R-list) number of boxes before iLevel on the remote rank
        !   + 1 is for the master box
        !   + E-boxes
        !       the number of extra boxes in the R-list before iBox
        !       there are (8 ** (rmLevel-ffLevel+1) - 1)/7 E-boxes (which compose E-tree) in rmRank
        !       for rmLevel == ffLevel + 1, there are 9 E-boxes in the E-tree
        !          the number 9 comes from 1 parent box followed by 8 boxes for its children
        !          the parent box (from ffLevel) is the actual extra multipole box for iRank
        !       for rmLevel > ffLevel, only the boxes above rmLevel are the actual E-boxes
        !       the children boxes (from rmLevel) at the end of E-tree are temporary data placeholders
        !   + the boxes of R-tree (from rmLevel down to iLevel-1)
        ! iRank is an RMG rank id

        ! find the number of extra multipole boxes (E-boxes) located before (above) iLevel in iRank
        if(iLevel >= rmLevel) then
          nxb = nExtraMultipoleBoxes(iRank) ! TODO replace this array with a test for iRank
        else
          ! count the number of E-boxes in E-tree in the levels above iLevel
!delete   nxb = 0
!delete   do eLevel = ffLevel, min0(iLevel-1, rmLevel-1) ! loop over levels above iLevel
!delete     nxb = nxb + merge(9, 0, iRank < 8**eLevel)   ! +9 = 1 + 8
!delete   enddo
          ! at this point, iBox is an E-box in remote iRank; determine the size of E-tree above iLevel
          ! every box at ffLevel starts an E-tree; one E-tree per rmRank
          ! only first 8**ffLevel rmRanks have E-boxes ! TODO replace 8**ffLevel with a constant
          !nxb = merge((8 ** (iLevel-ffLevel) - 1)/7, 0, iRank < 8**ffLevel)  ! 0 for iLevel == ffLevel
          nxb = (8 ** (iLevel-ffLevel) - 1)/7  ! no need to check for iRank < 8**ffLevel)
        endif

        ! bBefore is the number of boxes above iLevel in boxRegExpansion in the remote rank iRank
        ! each iLevel above rmLevel brings 8**(iLevel-ffLevel) E-boxes
        ! TODO rename bBefore to rBefore
        bBefore = 1 &
            + nxb   &  ! E-boxes: 1 parent + 8 children for rmLevel = ffLevel + 1
            + (minTopRtreeSize + merge(1, 0, iRank >= minTopRtreeRanks)) * m  ! =0 if iLevel <= rmLevel

        ! Determine the global (B-tree) index of the first box in R-list at iLevel on the remote iRank
        ! bStart is the total number of boxes (in B-tree) before the first box in the R-list at iLevel
        ! in other words, bStart is a global index pointer to the first box in R-list at iLevel
        ! TODO check that it doesn't overflow
        if(iLevel >= rmLevel) then
          bStart = nBoxesBefore &  ! the number of boxes before iLevel
             + minTopRtreeSize * num1_8 * iRank &  ! R-tree boxes held by other ranks before iRank
             + merge((iRank-minTopRtreeRanks) * num1_8, 0_8, iRank > minTopRtreeRanks)
        else
          ! at this point, iBox is an E-box in the remote rank iRank
          bStart = nBoxesBefore &  ! the number of boxes before iLevel
             + iRank * 8**(iLevel-ffLevel)  ! E-boxes held by other ranks before iRank at iLevel
!delete   bStart = nBoxesBefore &  ! the number of boxes before iLevel
!delete      + min0(iRank, 8 ** iLevel)  ! E-boxes held by other ranks before iRank at iLevel
!delete   ! TODO this test consumes too much resources
!delete   !if(iRank >= 64)  stop "Error in Plan: iRank >=  64 for E-box" ! ws = 2
!delete   !if(iRank >= 512) stop "Error in Plan: iRank >= 512 for E-box" ! ws = 3
        endif
      endif ! new segment

      ! iBox is a remote F-box
      ! multipole expansions of iBox will be computed in a remote rank and need to be retrieved
      ! need to provide a local (boxRegExpansion) storage for retrieved multipole expansions of iBox
      ! expStart is a pointer in boxRegExpansion for local storage of multipole expansions of iBox
      ! multipole expansions of the remote F-boxes in boxRegExpansion go after those of R-list
      !
      ! TODO check that it doesn't overflow
      expStart = (nBoxesRList + nBoxesThread(omp_get_thread_num())) * nSzTri1_8
      nBoxesThread(omp_get_thread_num()) = nBoxesThread(omp_get_thread_num()) + 1

      ! expRemote is a pointer for remote F-box in boxRegExpansion in the remote iRank
      ! the F-box is remote to rmRank, and, it is local to iRank
      ! (iBox - bStart) is the local box index for remote R-tree at iLevel
      ! (bBefore + iBox - bStart) is its local address on the remote rank iRank
      expRemote = (bBefore + iBox - bStart) * nSzTri1_8
      !if(expRemote < 0) stop "Error in Plan: expRemote < 0"
    else
      ! this is a local box for rmRank
      ! count the number of boxes in the gap before a new segment starts
      k = k + 1

      ! iRank == rmRank
      ! determine the number of boxes in the local (to rmRank) R-list before iLevel
      ! + 1 is for the master box
      ! + extra multipole boxes when rmLevel > ffLevel
      ! + the boxes of R-tree (from rmLevel down to iLevel-1)

      ! find the number of E-boxes before (above) iLevel for iRank
      if(iLevel >= rmLevel) then
        nxb = nExtraMultipoleBoxes(rmRank)
      else
        ! count the number of E-boxes in E-tree in the levels above iLevel
        ! every box at ffLevel starts an E-tree; one E-tree per rmRank
        ! determine the size of E-tree above iLevel
        ! only first 8**ffLevel rmRanks have E-boxes
        ! TODO replace 8**ffLevel with a constant
        !nxb = merge((8 ** (iLevel-ffLevel) - 1)/7, 0, rmRank < 8**ffLevel)  ! 0 for iLevel == ffLevel
        nxb = (8 ** (iLevel-ffLevel) - 1)/7  ! no need to check for rmRank < 8**ffLevel)
!delete ! each eligible level carries one E-box
!delete nxb = 0
!delete do eLevel = ffLevel, min0(iLevel-1, rmLevel-1) ! loop over levels above iLevel
!delete   nxb = nxb + merge(9, 0, rmRank < 8**eLevel)   ! +9 = 1 + 8
!delete enddo
      endif

      ! TODO rename bBefore to rBefore
      ! bBefore is the number of boxes above iLevel in boxRegExpansion in the local rmRank
      bBefore = 1 &
          + nxb   &  ! number of E-boxes above iLevel
          + (minTopRtreeSize + merge(1, 0, rmRank >= minTopRtreeRanks)) * m

      ! determine the global (B-tree) starting and ending box indices for R-list at iLevel in rmRank
      ! TODO it may potentially overflow for larger problem sizes
      if(iLevel >= rmLevel) then
        bStart = nBoxesBefore &
           + minTopRtreeSize * num1_8 * rmRank &
           + merge((rmRank-minTopRtreeRanks) * num1_8, 0_8, rmRank > minTopRtreeRanks)
!delete bEnd = bStart + minTopRtreeSize * num1_8 + merge(num1_8, 0_8, rmRank >= minTopRtreeRanks) - 1
      else
        bStart = nBoxesBefore &  ! the number of boxes before iLevel
           + rmRank * 8**(iLevel-ffLevel)  ! E-boxes held by other ranks before iRank at iLevel
!delete bStart = iBox
!delete bEnd = iBox
      endif

!delete ! check if iBox is within the global box indices for R-list
!delete ! only the box, which is a part of R-list, gets non-zero expStart (a pointer in boxRegExpansion)
!delete expStart = merge((bBefore + iBox - bStart) * nSzTri1_8, 0_8, iBox >= bStart .and. iBox <= bEnd)

      ! the expStart of the boxes of R-tree is a pointer to boxRegExpansion
      ! the value of expStart for a box outside of R-tree is meaningless
      ! if iBox is outside of R-tree, its expStart will not be used in FMM computation
      ! it is unnecessary to zero-initialize expStart (no need to check for R-tree boundaries)
      ! (bBefore + iBox - bStart) is the local box index of iBox in R-list in rmRank
      expStart  = (bBefore + iBox - bStart) * nSzTri1_8
      expRemote = 0_8
    endif

    ! sequentially renumber R+S+N+F boxes; save new box id in place of rank id

    nCubesThread(omp_get_thread_num()) = nCubesThread(omp_get_thread_num()) + 1

    ! get xyz from iBox
    aBox = iBox - nBoxesBefore   ! get the local index of iBox for iLevel
    xyz = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = nbParentEdge_init
    do aLevel = 1, iLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = aBox / nbParentCube        ! the parent box index that holds aBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      xyz(1) = xyz(1) + tjx * nbParentEdge
      xyz(2) = xyz(2) + tjy * nbParentEdge
      xyz(3) = xyz(3) + tjz * nbParentEdge

      aBox = mod(aBox,nbParentCube)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    rList(a)%xyz       = xyz
    rList(a)%expStart  = expStart
    rList(a)%expRemote = expRemote
  enddo
  !$omp end parallel

    !do a = 1, numBoxesSNRF-1
    !  iBox = getGlobalFromLocal(a)
    !  iRank = boxesRank(a)  ! iRank is the MPI World rank id
    !  iLevel = floor(dlog(dble(iBox*7+1)) * invLog8)
    !  iRank = iRank - 1  ! udjust for the shift in rank id
    !  if(rList(a)%expRemote < 0) then
    !    write(*,'(a,7i8)') "vam ", a, iBox, iLevel, iam, rmRank, iRank, rList(a)%expRemote
    !  endif
    !enddo
    !stop

#ifdef DEBUG_TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 5:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  !purpose: Initialize starting and ending values of particle indices in B-boxes for S-list and near field

  ! modify nAtomsCube to represent sequential particle ordering in B-boxes
  ! it automatically provides sequential particle ordering in S-boxes due to the ordering of S-boxes
  ! TODO run as an omp section

  iAtomStart =  1  ! use the positive value for S-boxes and their NF boxes
  jAtomStart = -1  ! use the negative value for direct summation boxes
  nParticlesRank = 0
#ifdef  DO_DIRECT_SUM
  do iBox = bBoxesStart, nBoxesTotal-1
  !do a = 1, numBoxesSNRF-1
  !  iBox = getGlobalFromLocal(a)

    ! get ix, iy, iz from iBox
    jBox = iBox - nBoxesBefore_loLevel   ! get the local index of iBox for iLevel
    xyz = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = nbParentEdge_init
    do aLevel = 1, loLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / (nbParentEdge * nbParentEdge * nbParentEdge)   ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      xyz(1) = xyz(1) + tjx * nbParentEdge
      xyz(2) = xyz(2) + tjy * nbParentEdge
      xyz(3) = xyz(3) + tjz * nbParentEdge

      jBox = mod(jBox, nbParentEdge * nbParentEdge * nbParentEdge)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    nC = xyz(1) &
       + xyz(2) * nChildEdge &
       + xyz(3) * nChildXY

!    if(btypeNS(iBox)) then  ! S or N
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits))) then  ! S or N
      ! choose only the boxes which carry particles required for the rank iam: S-list and near field
#ifdef USE_GETLOCAL_FROMGLOBAL
      idx_4 = getLocalFromGlobal(iBox)
#else
      idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
      rList(idx_4)%atomStart = iAtomStart
      rList(idx_4)%atomEnd   = iAtomStart + nAtomsCube(nC) - 1
      if(iBox >= iBoxesStart .and. iBox <= iBoxesEnd) then
        nParticlesSlist = nParticlesSlist + nAtomsCube(nC)
        cBox = iBox - iBoxesStart + sBoxLLstart
        sList(cBox)%atomStart  = iAtomStart
        sList(cBox)%atomEnd    = iAtomStart + nAtomsCube(nC) - 1
      endif
      iAtomStart           = iAtomStart + nAtomsCube(nC)   ! particle starting index for the next B-box
      nParticlesRank =   nParticlesRank + nAtomsCube(nC)   ! count the number of particles for the rank iam
      nAtomsCube(nC)       = iAtomStart - nAtomsCube(nC)   ! nAtomsCube(nC) = box(iBox)%atomStart
    else  ! neither S nor N
      ! use the particles from these boxes for accuracy validation by performing direct summation
      i = nAtomsCube(nC)           ! nAtomsCube(nC) holds a positive value at this moment
      nAtomsCube(nC) = jAtomStart  ! save a negative index in its place
      jAtomStart = jAtomStart - i  ! increment the negative index
    endif
  enddo
#else
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)

!    if(btypeNS(iBox)) then  ! the box type is N or S
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits))) then  ! the box type is N or S
      rList(a)%atomStart = iAtomStart
      rList(a)%atomEnd   = iAtomStart + nAtomsNS(a) - 1
      if(iBox >= iBoxesStart .and. iBox <= iBoxesEnd) then
        ! S-boxes at LL only
        nParticlesSlist = nParticlesSlist + nAtomsNS(a)
        cBox = iBox - iBoxesStart + sBoxLLstart
        sList(cBox)%atomStart  = iAtomStart
        sList(cBox)%atomEnd    = iAtomStart + nAtomsNS(a) - 1
      endif
      iAtomStart     =     iAtomStart + nAtomsNS(a)   ! particle starting index for the next S-box
      nParticlesRank = nParticlesRank + nAtomsNS(a)   ! count the number of particles for the rank iam
      nAtomsNS(a)    =     iAtomStart - nAtomsNS(a)   ! nAtomsNS(a) = box(iBox)%atomStart
    endif
  enddo
#endif

  !purpose: Renumber particles for two separate groups: (R) and (S + N + R)

  kAtomStart =  1  ! particle index for R boxes
  mAtomStart =  1  ! particle index for S + N + R boxes
#ifdef  DO_DIRECT_SUM
  do iBox = bBoxesStart, nBoxesTotal-1

    ! get ix, iy, iz from iBox
    jBox = iBox - nBoxesBefore_loLevel   ! get the local index of iBox for iLevel
    xyz = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = nbParentEdge_init
    do aLevel = 1, loLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      pBox = jBox / (nbParentEdge * nbParentEdge * nbParentEdge)   ! the parent box index that holds jBox (iBox)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      xyz(1) = xyz(1) + tjx * nbParentEdge
      xyz(2) = xyz(2) + tjy * nbParentEdge
      xyz(3) = xyz(3) + tjz * nbParentEdge

      jBox = mod(jBox, nbParentEdge * nbParentEdge * nbParentEdge)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo

    nC = xyz(1) &
       + xyz(2) * nChildEdge &
       + xyz(3) * nChildXY

    ! TODO run in a separate loop in an omp section
    ! particle ordering in the DD coordinate array for computation of regular multipoles
    if(iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd) then
      particleRanges(1,iBox-rBoxBottomStart) = kAtomStart
      particleRanges(2,iBox-rBoxBottomStart) = kAtomStart + boxRank(nC) - 1
      kAtomStart                       = kAtomStart + boxRank(nC)   ! particle starting index for the next B-box
      boxRank(nC)                      = kAtomStart - boxRank(nC)   ! boxRank(nC) = box(iBox)%atomStart
    else
      boxRank(nC) = 0  ! TODO can we do without it?
    endif

    ! TODO run in a separate loop in an omp section
    ! particle ordering in the buffer coordinate array (S + N + D)
!    if(btypeNS(iBox) &
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
            .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
      ! choose the boxes which carry particles required for S-list + NF + DD
      mAtomStart       = mAtomStart + backupBuffer(nC)
      backupBuffer(nC) = mAtomStart - backupBuffer(nC)
    else  ! neither S nor NF nor DD
      backupBuffer(nC) = 0  ! TODO can we do without it?
    endif
  enddo
#else
  ! TODO run in an omp section
  ! This loop runs on a single thread due to sequential index computation
  do a = 1, numBoxesSNRF-1
    iBox = getGlobalFromLocal(a)

    ! particle renumbering for the coordinate array cC for regular multipole computation in FMM Step 1
    if(iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd) then
      ! at this moment, nAtoms(a) holds the number of particles in iBox
      particleRanges(1,iBox-rBoxBottomStart) = kAtomStart                   ! particle start index
      particleRanges(2,iBox-rBoxBottomStart) = kAtomStart + nAtomsR(a) - 1  ! particle end index
      kAtomStart                       = kAtomStart + nAtomsR(a)  ! particle start index for the next R-box
      nAtomsR(a)                       = kAtomStart - nAtomsR(a)  ! nAtomsR(a) == box(iBox)%atomStart
      ! at this moment, nAtoms(a) holds the particle start index for iBox
      ! the particle index is local to coordinate array cC holding these particles
      ! the FMM Step 1 needs only these R-boxes and their particles, so particles can be renumbered
    endif

    ! particle ordering in the coordinate buffer array (S + N + R)
    ! buffer array provides a temporary particle storage to assist copying particles to different arrays
!    if(btypeNS(iBox) &
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
            .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
      ! at this moment, nAtomsSNR(a) holds a number of particles in iBox of type SNR
      mAtomStart   = mAtomStart + nAtomsSNR(a)
      nAtomsSNR(a) = mAtomStart - nAtomsSNR(a)
      ! at this moment, nAtomsSNR(a) holds a sequential starting particle index for boxes of type SNR
    endif
  enddo
#endif

  !purpose: Determine the starting and ending particle indices for S-boxes in the rank

#ifdef  DO_DIRECT_SUM
  nC = sList(sBoxLLstart)%xyz(1) &
     + sList(sBoxLLstart)%xyz(2) * nChildEdge &
     + sList(sBoxLLstart)%xyz(3) * nChildXY
  iAtomsBegin = nAtomsCube(nC)
#else
#ifdef USE_GETLOCAL_FROMGLOBAL
  iAtomsBegin = nAtomsNS(getLocalFromGlobal(iBoxesStart))
#else
  ! findloc is horribly slow and it simply cannot be used
  ! TODO see if storing getGlobalFromLocal in HBM will make findloc performance acceptable
  iAtomsBegin = nAtomsNS(findloc(getGlobalFromLocal, iBoxesStart, dim=1, kind=4) - 1)
#endif
#endif
  iAtomsEnd = iAtomsBegin + nParticlesSlist - 1

  rankAtoms(iam) = nParticlesSlist

  ! allocate memory
  memLocal = memLocal &
#ifdef DO_DIRECT_SUM
                      + dble(nAtomsTotal)*8*4 &
#endif
                      + dble(nParticlesRank)*8*4     &
                      + dble(nParticlesRMRank)*8*4 &
                      + dble(nParticlesBuffer)*8*4
  memMax = dmax1(memMax, (memUsed + memLocal))
#if defined(DEBUG_TIMER) && defined(MPI)
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(2(a,i12))') "dbg: nParticlesRank      prtRank: ", nParticlesRank,  "   rank: ", iam
    write(*,'(2(a,i12))') "dbg: nParticlesSlist     prtRank: ", nParticlesSlist, "   rank: ", iam
  endif
#endif

#ifndef MAX_PERF
#ifdef MPI
  ! determine the max memory requirement for each node
  if(nSplitRanks < nodeSize) then
    call MPI_AllReduce(memMax, memMaxMyNode, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_PNODE, mpiError)
  else
    call MPI_AllReduce(memMax, memMaxMyNode, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NODE, mpiError)
  endif
  if(iam == prtRank .or. debugEachMPIrank) &
          call MemNodePrint("Plan+Particles  myNodeMax", memLocal, memUsed+memLocal, memMax, memMaxMyNode)

  ! determine the max memory requirement per node across all nodes
  if(nSplitRanks < nodeSize) then
    ! FMM Plan uses multiple passes
    if(MPI_COMM_PZERO /= MPI_COMM_NULL) then
      call MPI_Reduce(memMaxMyNode, memMaxAllNodes, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_PZERO, mpiError)
    endif
    if(znodeRank == 0) call MemNodePrint("Plan+Particles allNodeMax", memLocal, memUsed+memLocal, memMax, memMaxAllNodes)
  else
    ! FMM Plan uses single pass
    if(nNodes > 1 .and. MPI_COMM_NZERO /= MPI_COMM_NULL) then
      call MPI_Reduce(memMaxMyNode, memMaxAllNodes, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_NZERO, mpiError)
      if(iam == 0) call MemNodePrint("Plan+Particles allNodeMax", memLocal, memUsed+memLocal, memMax, memMaxAllNodes)
    endif
  endif
#endif
#endif  ! MAX_PERF
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan + Particles", memLocal, memUsed+memLocal, memMax, planSize)

  allocate( &
#ifdef DO_DIRECT_SUM
           aC(3,nAtomsTotal), aQ(nAtomsTotal), &
#endif
           bC(3,nParticlesRank), bQ(nParticlesRank), &
           cC(3,nParticlesRMRank), cQ(nParticlesRMRank), &
           bufferC(3,nParticlesBuffer), bufferQ(nParticlesBuffer), stat=allocStat)
  if(allocStat /= 0) then
    if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan", memLocal, memUsed+memLocal, memMax, planSize)
    call exitPrg(2)
  endif

#ifdef DEBUG_TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 6:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  !$omp parallel num_threads(nThreadsLocal)

  ! renumber boxes in NF-list
  !$omp do
  do iBox = 1, nNearNeighbors
#ifdef USE_GETLOCAL_FROMGLOBAL
    nearNeighbors(iBox) = getLocalFromGlobal(nearNeighbors(iBox))
#else
    nearNeighbors(iBox) = findloc(getGlobalFromLocal, nearNeighbors(iBox), dim=1, kind=4) - 1
#endif
  enddo

  ! renumber boxes in FF-list
  !$omp do
  do iBox = 1, nFarNeighbors
#ifdef USE_GETLOCAL_FROMGLOBAL
    farNeighbors(iBox) = getLocalFromGlobal(farNeighbors(iBox))
#else
    farNeighbors(iBox) = findloc(getGlobalFromLocal, farNeighbors(iBox), dim=1, kind=4) - 1
#endif
  enddo

#ifdef DEBUG_TIMER
  !$omp end parallel
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 7:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
  !$omp parallel num_threads(nThreadsLocal)
#endif

  !purpose: Extract S + N + R particles needed for the rank and sort them to boxes

  ! store the particles in a temporary buffer
  ! TODO use inverted coordinates to make the enclosed loop more efficient
#ifdef DO_DIRECT_SUM
  !$omp do private(i,j,k, translation, iAtom, ix,iy,iz, nC, nextAtom)
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the cube
    translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the cube
    translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the cube
    !do iAtom = 1, nAtoms
    do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
      !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
      ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
      iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
      iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)
      nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and S-list

      ! copy particles assigned to the MPI rank, iam
      !$omp critical
      if(backupBuffer(nC) > 0) then
        nextAtom = backupBuffer(nC)

        bufferC(1,nextAtom) = cell(1,iAtom) + translation(1) + vertexMin(1)
        bufferC(2,nextAtom) = cell(2,iAtom) + translation(2) + vertexMin(2)
        bufferC(3,nextAtom) = cell(3,iAtom) + translation(3) + vertexMin(3)
        bufferQ(nextAtom)   = w(iAtom)

        !!$omp atomic update
        backupBuffer(nC) = backupBuffer(nC) + 1       ! increment next particle index for this box
      endif
      !$omp end critical
    enddo
  enddo
#else ! not DO_DIRECT_SUM

#ifdef FLEX_SIDELENGTH
  !$omp do private(i,j, ij,ijk, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, idx_4)
  do k = 0, nRepLayers, 2
  do ij = 0, nRepLayers2 - 1  ! if this is uncommented, make sure ijk is declared as a private variable
!  do ijk = 0, nReplicas - 1
!    k = ijk / nRepLayers2     ! if this is uncommented, make sure k is declared as a private variable
!    j = (ijk - k*nRepLayers2) / nRepLayers1
    ijk = k * nRepLayers2 + ij
    j = ij / nRepLayers1
    i = mod(ij, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with S-NF or DD spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the cube
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the cube
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the cube
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

        !nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and S-list

        ! copy particles assigned to the MPI rank, iam
!        if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          !!$omp critical
          nextAtom = nAtomsSNR(idx_4)

          bufferC(1,nextAtom) = cell(1,iAtom) + translation(1) + vertexMin(1)
          bufferC(2,nextAtom) = cell(2,iAtom) + translation(2) + vertexMin(2)
          bufferC(3,nextAtom) = cell(3,iAtom) + translation(3) + vertexMin(3)
          bufferQ(nextAtom)   = w(iAtom)

          ! each replica consists of a unique set of boxes so the atomic is unnecessary
          !!$omp atomic update
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! increment next particle index for this box
          !!$omp end critical
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo ! ijk
  enddo ! k

  !$omp barrier

  !$omp do private(i,j, ij,ijk, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, idx_4)
  do k = 1, nRepLayers, 2
  do ij = 0, nRepLayers2 - 1  ! if this is uncommented, make sure ijk is declared as a private variable
!  do ijk = 0, nReplicas - 1
!    k = ijk / nRepLayers2     ! if this is uncommented, make sure k is declared as a private variable
!    j = (ijk - k*nRepLayers2) / nRepLayers1
    ijk = k * nRepLayers2 + ij
    j = ij / nRepLayers1
    i = mod(ij, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with S-NF or DD spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the cube
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the cube
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the cube
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

        !nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and S-list

        ! copy particles assigned to the MPI rank, iam
!        if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          !!$omp critical
          nextAtom = nAtomsSNR(idx_4)

          bufferC(1,nextAtom) = cell(1,iAtom) + translation(1) + vertexMin(1)
          bufferC(2,nextAtom) = cell(2,iAtom) + translation(2) + vertexMin(2)
          bufferC(3,nextAtom) = cell(3,iAtom) + translation(3) + vertexMin(3)
          bufferQ(nextAtom)   = w(iAtom)

          ! each replica consists of a unique set of boxes so the atomic is unnecessary
          !!$omp atomic update
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! increment next particle index for this box
          !!$omp end critical
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo ! ijk
  enddo ! k

#else  ! not FLEX_SIDELENGTH

  !$omp do private(i,j,k, translation, iAtom, ix,iy,iz, iBox, xloCell, yloCell, zloCell) &
#if defined(DEBUG_TIMER) && defined(MPI)
  !$omp reduction(+:nReplicasUsed) &
#endif
#ifdef MPI
  !$omp private(replicaUsed, cellSize1, cellSize2, cellSize3, ijkCube, xCube, yCube, zCube) &
#endif
  !$omp private(xupCell, yupCell, zupCell, nextAtom, idx_4) 
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    ! check if the ijk replica of the primitive cell overlaps with S-NF or DD spatial boundaries
    ! skip the replicas not needed on the MPI rank
    xloCell = (i * intPrimitiveEdge) / intCubeEdge
    yloCell = (j * intPrimitiveEdge) / intCubeEdge
    zloCell = (k * intPrimitiveEdge) / intCubeEdge
    xupCell = min0(lastLLBoxAxis, ((i+1) * intPrimitiveEdge) / intCubeEdge)
    yupCell = min0(lastLLBoxAxis, ((j+1) * intPrimitiveEdge) / intCubeEdge)
    zupCell = min0(lastLLBoxAxis, ((k+1) * intPrimitiveEdge) / intCubeEdge)
    if((((xloCell >= minRx .and. xloCell <= maxRx) .or. &
         (xupCell >= minRx .and. xupCell <= maxRx) .or. &
         (xloCell <= minRx .and. xupCell >= maxRx)      &
        ) .and. &
        ((yloCell >= minRy .and. yloCell <= maxRy) .or. &
         (yupCell >= minRy .and. yupCell <= maxRy) .or. &
         (yloCell <= minRy .and. yupCell >= maxRy)      &
        ) .and. &
        ((zloCell >= minRz .and. zloCell <= maxRz) .or. &
         (zupCell >= minRz .and. zupCell <= maxRz) .or. &
         (zloCell <= minRz .and. zupCell >= maxRz)      &
        ) &
       ) .or. &
       (((xloCell >= minNx .and. xloCell <= maxNx) .or. &
         (xupCell >= minNx .and. xupCell <= maxNx) .or. &
         (xloCell <= minNx .and. xupCell >= maxNx)      &
        ) .and. &
        ((yloCell >= minNy .and. yloCell <= maxNy) .or. &
         (yupCell >= minNy .and. yupCell <= maxNy) .or. &
         (yloCell <= minNy .and. yupCell >= maxNy)      &
        ) .and. &
        ((zloCell >= minNz .and. zloCell <= maxNz) .or. &
         (zupCell >= minNz .and. zupCell <= maxNz) .or. &
         (zloCell <= minNz .and. zupCell >= maxNz)      &
        ) &
       ) &
    ) then
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the cube
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the cube
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the cube
      !
#ifdef MPI
      cellSize1 = (xupCell-xloCell+1)              ! number of LL boxes along X-axis in the replica
      cellSize2 = (yupCell-yloCell+1) * cellSize1  ! number of LL boxes in XY-plane in the replica
      cellSize3 = (zupCell-zloCell+1) * cellSize2  ! number of LL boxes in XYZ-box in the replica
      ! test if the replica carries the required boxes
      replicaUsed = 0  ! begin from the assumption that the replica is not needed
      ! break the primitive cell (replica) down to individual loLevel cubes
      do ijkCube = 0, cellSize3 - 1
        ! local cube coordinates in the primitive cell
        zCube = ijkCube / cellSize2
        yCube = mod(ijkCube, cellSize2) / cellSize1
        xCube = ijkCube - zCube*cellSize2 - yCube*cellSize1
        ! global flat coordinates of the cubes after translation
        ix = xloCell + xCube
        iy = yloCell + yCube
        iz = zloCell + zCube
 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1
 ! TODO replace with  nBoxesBefore_loLevel = (8_8 ** loLevel - 1_8) / 7_8
 ! TODO check which one of the above two is correct
        !
        ! test if the box iBox is needed this rank
!       if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
          replicaUsed = 1  ! the replica is needed
          exit             ! skip the rest of the test
        endif
      enddo

      ! count the number of particles in the boxes
      if(replicaUsed == 1) then
#endif
      !do iAtom = 1, nAtoms
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
        ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
        iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
        iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

        !nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and S-list

        ! copy particles assigned to the MPI rank, iam
!        if(btypeNS(iBox) &
        if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits)) &
                .or. (iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd)) then
#ifdef USE_GETLOCAL_FROMGLOBAL
          idx_4 = getLocalFromGlobal(iBox)
#else
          idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
          nextAtom = nAtomsSNR(idx_4)

          bufferC(1,nextAtom) = cell(1,iAtom) + translation(1) + vertexMin(1)
          bufferC(2,nextAtom) = cell(2,iAtom) + translation(2) + vertexMin(2)
          bufferC(3,nextAtom) = cell(3,iAtom) + translation(3) + vertexMin(3)
          bufferQ(nextAtom)   = w(iAtom)

          ! each replica consists of a unique set of boxes so the atomic is unnecessary
          !!$omp atomic update
          nAtomsSNR(idx_4) = nAtomsSNR(idx_4) + 1  ! increment next particle index for this box
        endif
      enddo  ! do iAtom = 1, nAtoms
#ifdef MPI
      endif  ! if(replicaUsed == 1)
#endif
#if defined(DEBUG_TIMER) && defined(MPI)
      nReplicasUsed = nReplicasUsed + replicaUsed
#endif
    endif
  enddo

#endif ! not FLEX_SIDELENGTH
#endif ! not DO_DIRECT_SUM

#ifdef DEBUG_TIMER
  !$omp end parallel
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
#if defined(DEBUG_TIMER) && defined(MPI)
    nReplicasNeeded = max0(nReplicasNeeded, nReplicasUsed)
    nReplicasUsed = 0
    write(*,'(2(a,i12))')    "dbg: nReplicasNeeded            : ", nReplicasNeeded, " out of: ", nReplicasRank
#endif
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 8:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
  !$omp parallel num_threads(nThreadsLocal)
#endif

  !purpose: Sort particles to boxes

#ifdef DO_DIRECT_SUM
  ! direct summation currently requires all particles in the system to be available to every MPI rank
  ! sort particles to boxes so the particles will be arranged sequentially for each box
  ! TODO use inverted coordinates to make the enclosed loop more efficient
  !$omp do private(i,j,k, translation, iAtom, ix,iy,iz, nC, nextAtom, m)
  do ijk = 0, nReplicas - 1
    k = ijk / nRepLayers2
    j = (ijk - k*nRepLayers2) / nRepLayers1
    i = mod(ijk, nRepLayers1)
    translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)  ! replicate the cube
    translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)  ! replicate the cube
    translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)  ! replicate the cube
    !do iAtom = 1, nAtoms
    do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
      !coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
      ix = int((cell(1,iAtom) + translation(1)) * edgeReciprocal)
      iy = int((cell(2,iAtom) + translation(2)) * edgeReciprocal)
      iz = int((cell(3,iAtom) + translation(3)) * edgeReciprocal)
      nC = ix + iy*nChildEdge + iz*nChildXY  ! this (flat) index is different from that in B-list and S-list

      ! increment atom counter
      !$omp critical
      if(nAtomsCube(nC) > 0) then
        ! particles assigned to the MPI rank
        nextAtom = nAtomsCube(nC)
        m = 1
      elseif(nAtomsCube(nC) < 0) then
        ! particles needed for direct summation
        nextAtom = nParticlesRank - nAtomsCube(nC)
        m = -1
      endif

      ! copy the particles assigned to the MPI rank
      aC(1,nextAtom) = cell(1,iAtom) + translation(1) + vertexMin(1)
      aC(2,nextAtom) = cell(2,iAtom) + translation(2) + vertexMin(2)
      aC(3,nextAtom) = cell(3,iAtom) + translation(3) + vertexMin(3)
      aQ(nextAtom)   = w(iAtom)

      !!$omp atomic update
      nAtomsCube(nC) = nAtomsCube(nC) + m       ! increment next particle index for this box
      !$omp end critical
    enddo
  enddo
#else
  ! select S + N particles from the buffer containing S + N + R particles needed for the rank
  ! sort particles to boxes so the particles are alined with the ordering of boxes
  !$omp do private(ix,iy,iz, iBox, nextAtom, idx_4)
  do iAtom = 1, nParticlesBuffer

    ! find the box index for iAtom
    ix = int((bufferC(1,iAtom) - vertexMin(1)) * edgeReciprocal)
    iy = int((bufferC(2,iAtom) - vertexMin(2)) * edgeReciprocal)
    iz = int((bufferC(3,iAtom) - vertexMin(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

!    if(btypeNS(iBox)) then
    if(btest(boxTypeNS%lArray(iBox/nBits), mod(iBox,nBits))) then
#ifdef USE_GETLOCAL_FROMGLOBAL
      idx_4 = getLocalFromGlobal(iBox)
#else
      idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
      nextAtom = nAtomsNS(idx_4)

      bC(1,nextAtom) = bufferC(1,iAtom)
      bC(2,nextAtom) = bufferC(2,iAtom)
      bC(3,nextAtom) = bufferC(3,iAtom)
      bQ(nextAtom)   = bufferQ(iAtom)

      !$omp atomic update
      nAtomsNS(idx_4) = nAtomsNS(idx_4) + 1  ! increment next particle index for this box
    endif
  enddo
#endif

#ifdef DEBUG_TIMER
  !$omp end parallel
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan 9:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
  !$omp parallel num_threads(nThreadsLocal)
#endif

  !purpose: Select R particles from the buffer containing S + N + R particles needed for the rank

  ! sort particles to boxes so the particles will be aligned with the box order
#ifdef DO_DIRECT_SUM
  !$omp do private(nC, nextAtom)
  do iAtom = 1, nParticlesBuffer
    ! identify the enclosing box by its grid coordinates
    nC = int((bufferC(1,iAtom) - vertexMin(1)) * edgeReciprocal) &
       + int((bufferC(2,iAtom) - vertexMin(2)) * edgeReciprocal) * nChildEdge &
       + int((bufferC(3,iAtom) - vertexMin(3)) * edgeReciprocal) * nChildXY

    ! copy particles assigned to the MPI rank, iam
    !$omp critical
    if(boxRank(nC) > 0) then
      nextAtom = boxRank(nC)

      cC(1,nextAtom) = bufferC(1,iAtom)
      cC(2,nextAtom) = bufferC(2,iAtom)
      cC(3,nextAtom) = bufferC(3,iAtom)
      cQ(nextAtom)   = bufferQ(iAtom)

      !!$omp atomic update
      boxRank(nC) = boxRank(nC) + 1       ! increment next particle index for this box
    endif
    !$omp end critical
  enddo
#else
  !$omp do private(ix,iy,iz, iBox, nextAtom, idx_4)
  do iAtom = 1, nParticlesBuffer

    ! find the box index for iAtom
    ix = int((bufferC(1,iAtom) - vertexMin(1)) * edgeReciprocal)
    iy = int((bufferC(2,iAtom) - vertexMin(2)) * edgeReciprocal)
    iz = int((bufferC(3,iAtom) - vertexMin(3)) * edgeReciprocal)

 ! TODO find out what is the largest supported loLevel
 iBox =                                     (iz/2048_8)*34359738368_8                                 + & ! 12
                                        (mod(iz,2048_8)/1024_8)*4294967296_8                            + & ! 11
                                    (mod(mod(iz,2048_8),1024_8)/512_8)*536870912_8                          + & ! 10
                                (mod(mod(mod(iz,2048_8),1024_8),512_8)/256_8)*67108864_8                      + & !  9
                            (mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8)/128_8)*8388608_8                  + & !  8
                        (mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8)/64_8)*1048576_8              + & !  7
                    (mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*131072_8           + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*16384_8        + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*2048_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*256_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*32_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iz,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*4_8  + & !  2
                                            (iy/2048_8)*17179869184_8                                 + & ! 12
                                        (mod(iy,2048_8)/1024_8)*2147483648_8                            + & ! 11
                                    (mod(mod(iy,2048_8),1024_8)/512_8)*268435456_8                          + & ! 10
                                (mod(mod(mod(iy,2048_8),1024_8),512_8)/256_8)*33554432_8                      + & !  9
                            (mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8)/128_8)*4194304_8                  + & !  8
                        (mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8)/64_8)*524288_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*65536_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*8192_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*1024_8      + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*128_8    + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*16_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(iy,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)*2_8  + & !  2
                                            (ix/2048_8)*8589934592_8                                  + & ! 12
                                        (mod(ix,2048_8)/1024_8)*1073741824_8                              + & ! 11
                                    (mod(mod(ix,2048_8),1024_8)/512_8)*134217728_8                          + & ! 10
                                (mod(mod(mod(ix,2048_8),1024_8),512_8)/256_8)*16777216_8                      + & !  9
                            (mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8)/128_8)*2097152_8                  + & !  8
                        (mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8)/64_8)*262144_8               + & !  7
                    (mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8)/32_8)*32768_8            + & !  6
                (mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8)/16_8)*4096_8         + & !  5
            (mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8)/8_8)*512_8       + & !  4
        (mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8)/4_8)*64_8     + & !  3
    (mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8)/2_8)*8_8 + & !  2
 mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(mod(ix,2048_8),1024_8),512_8),256_8),128_8),64_8),32_8),16_8),8_8),4_8),2_8)     + & !  2
 nBoxesBefore_loLevel
 !(8_8**loLevel)/7_8  ! the total number of boxes for a system that has loLevel-1

    ! copy particles assigned to the MPI rank, iam
    if(iBox >= rBoxBottomStart .and. iBox <= rBoxBottomEnd) then
#ifdef USE_GETLOCAL_FROMGLOBAL
      idx_4 = getLocalFromGlobal(iBox)
#else
      idx_4 = findloc(getGlobalFromLocal, iBox, dim=1, kind=4) - 1  ! -1 needed for 0-based array
#endif
      nextAtom = nAtomsR(idx_4)
    
      cC(1,nextAtom) = bufferC(1,iAtom)
      cC(2,nextAtom) = bufferC(2,iAtom)
      cC(3,nextAtom) = bufferC(3,iAtom)
      cQ(nextAtom)   = bufferQ(iAtom) 
      
      !$omp atomic update
      nAtomsR(idx_4) = nAtomsR(idx_4) + 1  ! increment next particle index for this box
    endif
  enddo
#endif

  !$omp end parallel

#ifdef DEBUG_TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifdef MPI
  if(nSplitRanks == nodeSize) then
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  else
    call MPI_Barrier(MPI_COMM_PGROUP, mpiError)
  endif
!  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,2(f7.1,a))') "Timer: FMM Plan10:      min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))') "Timer: Plan Consumed       :", TimerCheck(2), "s"
  endif
  call TimerInit(3)
#endif

  ! set the starting address for regular expansions in R-list
  !expShift = box(rBoxBottomStart)%expStart
#ifdef USE_GETLOCAL_FROMGLOBAL
  expShift = rList(getLocalFromGlobal(rBoxBottomStart))%expStart
#else
  expShift = rList(findloc(getGlobalFromLocal, rBoxBottomStart, dim=1, kind=4) - 1)%expStart
#endif

#ifdef DO_DIRECT_SUM
  ! copy the particles assigned to the MPI rank
  bC(1:3,1:nParticlesRank) = aC(1:3,1:nParticlesRank)
  bQ(1:nParticlesRank) = aQ(1:nParticlesRank)
  lComputeDirectSum = .TRUE.  ! do not change this parameter anywhere else in the code
#else
  lComputeDirectSum = .FALSE. ! do not change this parameter anywhere else in the code
#endif

  ! deallocate memory
  memLocal = memLocal &
#ifdef DO_DIRECT_SUM
                      - dble((max0(8**loLevel,nBoxesSList*2)+1))*4 & ! nAtomsCube
                      - dble((max0(8**loLevel,nBoxesSList*2)+1))*4 & ! boxRank
                      - dble((max0(8**loLevel,nBoxesSList*2)+1))*4 & ! backupBuffer
#else
                      - dble(nBoxesSList*2+1)*4                   &  ! boxRank
#endif
                     !- dble(nNearNeighbors*8)                    &  ! nearNeighbors8
                     !- dble(nFarNeighbors*8)                     &  ! farNeighbors8
                      - dble(nBoxesSList*8)                       &  ! idx
                      - dble(numBoxesSNRF)*8                      &  ! getGlobalFromLocal
                      - dble(numBoxesSNRF)*4                      &  ! boxesRank
#ifdef USE_GETLOCAL_FROMGLOBAL
                      - dble(nBoxesTotal+1)*4                     &  ! getLocalFromGlobal
#endif
                      - dble(numBoxesSNRF)*4                      &  ! nAtomsR
                      - dble(numBoxesSNRF)*4                      &  ! nAtomsNS
                      - dble(numBoxesSNRF)*4                      &  ! nAtomsSNR
!                     - dble(nBoxesTotal_div8)                    &  ! btypeNS
!                     - dble(nBoxesTotal_div8)                    &  ! btypeF
                      - dble(nBoxesTotal_div8/nBits)              &  ! boxTypeNS%lArray
                      - dble(nBoxesTotal_div8/nBits)              &  ! boxTypeF%lArray
                      - dble(nParticlesBuffer)*8*3                &  ! bufferC
                      - dble(nParticlesBuffer)*8                  &  ! bufferQ
                      - dble(nThreadsLocal*4*3)                   &  ! nSegmentsThread, nBoxesThread, nCubesThread
                      - dble(nDataSets*16_8) - dble(nReplicas*2*2)   ! atomRange, datasetid, groupid
  memUsed = memUsed + memLocal
  if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan Final        ", memLocal, memUsed, memMax, planSize)

  deallocate( &
#ifdef DO_DIRECT_SUM
             nAtomsCube, boxRank, backupBuffer,           &
#else
             boxRank,                                     &
#endif
            !nearNeighbors8, farNeighbors8,               &
             idx, getGlobalFromLocal, boxesRank,          &
#ifdef USE_GETLOCAL_FROMGLOBAL
             getLocalFromGlobal,                          &
#endif
             nAtomsR, nAtomsNS, nAtomsSNR,                &
!            btypeNS,                                     &
!            btypeF,                                      &
             boxTypeNS%lArray,                            &
             boxTypeF%lArray,                             &
             bufferC, bufferQ,                            &
             nSegmentsThread, nBoxesThread, nCubesThread, &
             atomRange, datasetid, groupid,               &
             stat=allocStat)
  if(allocStat /= 0) then
    if(iam == prtRank .or. debugEachMPIrank) call MemPrint("Plan", memLocal, memUsed+memLocal, memMax, planSize)
    call exitPrg(2)
  endif

  nNearNeighborsMin = nNearNeighbors
  nNearNeighborsMax = nNearNeighbors
  nFarNeighborsMin  = nFarNeighbors
  nFarNeighborsMax  = nFarNeighbors
!#ifdef MPI
!  call MPI_Reduce(nNearNeighbors, nNearNeighborsMin, 1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  call MPI_Reduce(nNearNeighbors, nNearNeighborsMax, 1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
!  call MPI_Reduce(nFarNeighbors,  nFarNeighborsMin,  1, MPI_INTEGER8, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
!  call MPI_Reduce(nFarNeighbors,  nFarNeighborsMax,  1, MPI_INTEGER8, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
!#endif

  if(iam == prtRank .or. debugEachMPIrank) then
   !write(*,'(a,i12)')   "Far-field level,            ffLevel: ", ffLevel
   !write(*,'(a,i12)')   "Particle division level,    rmLevel: ", rmLevel
   !write(*,'(a,i12)')   "Domain decomposition level, ddLevel: ", ddLevel
   !write(*,'(a,i12)')   "The lowest hierarchy level, loLevel: ", loLevel
    write(*,'(a,i12)')   "nFloatingBoxes, estimated          : ", nFloatingBoxesEstimated
    write(*,'(a,i12)')   "nFloatingBoxes, actual             : ", nFloatingBoxes
    write(*,'(a,2i12)')  "nExtraMultipoleBoxes  for  rmRank  : ", nExtraMultipoleBoxes(rmRank), rmRank
    write(*,'(a,i12)')   "nBoxesTotal                        : ", nBoxesTotal
    write(*,'(a,i12)')   "nBoxesSList S-tree+floating boxes  : ", nBoxesSList
    write(*,'(a,i12)')   "numBoxesSNRF S+N+R+F boxes         : ", numBoxesSNRF
    write(*,'(a,i12)')   "nBoxesRList (+ master + extra_mult): ", nBoxesRList
    write(*,'(a,i12)')   "nBoxes loLevel total               : ", 8_8 ** loLevel
    write(*,'(a,i12)')   "nBoxes loLevel in S-list           : ", nBoxesSList - sBoxLLstart
    write(*,'(a,i12)')   "nRemoteFFBoxes per rank            : ", nRemoteFFBoxes
    write(*,'(a,i12)')   "nBoxSegments                       : ", nBoxSegments
   !write(*,'(a,32i8)')  "Segments starting address per core : ", nSegmentsThread
   !write(*,'(a,32i8)')  "Num. of boxes in segments per core : ", nBoxesThread
   !write(*,'(a,32i8)')  "Num. of DD+S+NF+FF boxes  per core : ", nCubesThread
    write(*,'(a,i12)')   "nReplicas of the unit cell         : ", nReplicas
    write(*,'(a,i12)')   "nParticlesSlist                    : ", nParticlesSlist
    write(*,'(a,i12)')   "nParticlesRank                     : ", nParticlesRank
    write(*,'(a,i12)')   "nParticlesRMRank                   : ", nParticlesRMRank
    write(*,'(a,i12)')   "nParticlesBuffer                   : ", nParticlesBuffer
    write(*,'(a,2(i12,a))') "nNearNeighbors for S-list      min : ", nNearNeighborsMin, "   max: ", &
                                                                  nNearNeighborsMax, " boxes"
    write(*,'(a,2(i12,a))') "nFarNeighbors  for S-list      min : ", nFarNeighborsMin,  "   max: ", &
                                                                  nFarNeighborsMax,  " boxes"
  endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,1(f7.1,a,i0))')  "Timer: FMM Plan pass time          :", TimerCheck(3), "s  for rank ", prtRank
  endif
  call TimerInit(3)
#endif

#ifdef MPI
  if(nSplitRanks < nodeSize) then
    ! the FMM Plan runs in multiple passes to reduce the memory usage in the node
    !
    rankReceiver = nodeRank + nSplitRanks  ! rankReceiver is the rank to be woken up
    if(rankReceiver < nodeSize) then
      ! nodeRank notifies rankReceiver to wake-up
      i = 0
      call MPI_Send(i, 1, MPI_INTEGER, rankReceiver, nodeRank, MPI_COMM_NODE, mpiError)
    endif
    !
    ! restore the original value of prtRank
    prtRank = prtRankSaved
  endif
#endif

#ifdef TIMER
#ifdef MPI
  ! wait for all MPI ranks to finish the FMM Plan and deallocate the temporary arrays in the node before continuing
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  !call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  !call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
  !
  ! free the sub-communicators
  if(nSplitRanks < nodeSize) then
    call MPI_Comm_Free(MPI_COMM_PNODE, mpiError)
    call MPI_Comm_Free(MPI_COMM_PGROUP, mpiError)
    if(MPI_COMM_PZERO /= MPI_COMM_NULL) call MPI_Comm_Free(MPI_COMM_PZERO, mpiError)
  endif
#endif
  if(iam == prtRank .or. debugEachMPIrank) then
   !write(*,'(/a,2(f7.1,a))') "Timer: Total FMM Plan   min:", TimeConsumedMin, "s   max:", TimeConsumedMax,"s"
    write(*,'(a,1(f7.1,a))')  "Timer: FMM Plan total for all ranks:", TimerCheck(2), "s"
  endif
  call TimerInit(2)
#endif

end subroutine CreatePlan


subroutine FreeCubes
  use common_data, only: aC, aQ, bC, bQ, cC, cQ, cell, w, lowestLevel, rankAtoms, dispAtoms, &
          rankBoxes, dispBoxes, fBoxes, nExtraMultipoleBoxes, boxIrrExpansion, nearNeighbors, &
          farNeighbors, segmentRank, segmentNBoxes, particleRanges, sBoxesTopNum, sBoxesTopDispl
  implicit none
  integer :: allocStat
  deallocate( &
             rList, bC, bQ, treeLevel, &
             sList, fBoxes, nExtraMultipoleBoxes, &
#ifdef DO_DIRECT_SUM
             aC, aQ, &
#endif
             cC, cQ, particleRanges, cell, w, sBoxesTopNum, sBoxesTopDispl, &
             rankAtoms, dispAtoms, rankBoxes, dispBoxes, nearNeighbors, farNeighbors, &
             segmentSBox, segmentRank, segmentNBoxes, &
             stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
  ! no memory-usage accounting is needed here because computation is over
  ! unless there is a memory leak
end subroutine FreeCubes


subroutine RunCubeEdgeSimulation(edgeSize,nIter,vertexMin,vertexMax,iam)
  use common_data, only: prtRank
!
! Show the effect of choosing the right box edge
!   on how many empty boxes will be added to satisfy 
!   the requirement of power-2 recursive division
!
! Smaller boxes are more preferrable for WS=2, whereas larger boxes
!   are better for WS=1
!
  implicit none
  double precision, intent(in) :: edgeSize, vertexMin(3), vertexMax(3)
  integer, intent(in)          :: nIter, iam

  integer          :: i,  nX, nY, nZ, n1, n2, depth
  double precision :: edge, log2, reciprocal

  if(iam == prtRank) then
    write(*,'(/5x,a)') "Simulation to show the effect of different cubeEdges"
  endif

  edge = edgeSize
  log2 = log(2.0d0)

  do i = 1, nIter
    reciprocal = 1.0d0 / edge

    ! determine number of boxes along each axis
    nX = int((vertexMax(1) - vertexMin(1)) * reciprocal) + 1
    nY = int((vertexMax(2) - vertexMin(2)) * reciprocal) + 1
    nZ = int((vertexMax(3) - vertexMin(3)) * reciprocal) + 1
    n1 = max0(nX, nY, nZ)

    ! determine number of boxes to satisfy the recursive division-by-powers-of-2 requirement
    depth = int(log(dble(n1)) / log2)
    depth = depth + min0(max0(n1 - 2**depth, 0), 1) ! if(n1 > 2**depth) depth=depth+1
    n2 = 2 ** depth

    if(iam == prtRank) then
      write(*,'(5x,a,f4.1,a,3i8, 6x,a,3i8, 6x,a,3i8)') &
      "edge = ", edge, "  boxes:", nX,nY,nZ, " power2:", n2,n2,n2, &
      "depth: ", depth,depth,depth
    endif

    edge = edge + 1.0d0
  enddo

end subroutine RunCubeEdgeSimulation


! the original quicksort code is from https://gist.github.com/1AdAstra1/6f7785373efe5bb6c254d2e20c78ccc4
! bug: the original algorithm assumes that the starting index of the array is 1 whereas it can be 0
recursive subroutine quicksort4(a)
  implicit none
  integer :: a(:)
  integer :: x, t
  integer :: first = 1, last
  integer :: i, j

  last = size(a, 1)
  if(last > 1) then
    x = a( (first+last) / 2 )
    i = first
    j = last

    do
       ! find a(i) >= x
       do while (a(i) < x)
          i=i+1
       end do
       ! find a(j) <= x
       do while (a(j) > x)
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do

    if (first < i - 1) call quicksort4(a(first : i - 1))
    if (j + 1 < last)  call quicksort4(a(j + 1 : last))
  endif
end subroutine quicksort4


recursive subroutine quicksort8(a)
  implicit none
  integer (kind=8) :: a(:)
  integer (kind=8) :: x, t
  integer (kind=8) :: first = 1, last
  integer (kind=8) :: i, j

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last

  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do

  if (first < i - 1) call quicksort8(a(first : i - 1))
  if (j + 1 < last)  call quicksort8(a(j + 1 : last))
end subroutine quicksort8


integer (kind=2) function compar8( a, b )
  implicit none
  integer (kind=8) a, b
  if ( a .lt. b ) compar8 = -1
  if ( a .eq. b ) compar8 = 0
  if ( a .gt. b ) compar8 = 1
  return
end function compar8


end module fmm_plan
