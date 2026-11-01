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

#undef  USE_REDUCE
!#define USE_REDUCE

!#undef  USE_HOSTRDMA
#define USE_HOSTRDMA

! maximize the performance by disabling the non-essential diagnostics
! one should keep the diagnostics enabled in the day-to-day runs
!#define MAX_PERF
#undef MAX_PERF

! this section (using aPotential and aForce) is unusable and requires a redesign
#undef  FETCH_POTFORCE
!#define FETCH_POTFORCE

! compute force
!#undef  COMP_FORCE
#define COMP_FORCE

! sum up force for selected groups of atoms and return force for each selected atom
#undef  FSEL
!#define FSEL

! sum up force for a single selected group of atoms and return a force sum only
#undef  FSUM
!#define FSUM

! measure runtime of all ranks in step 6
!#undef TIME_STEP6
#define TIME_STEP6

! dependencies between the defines
#ifdef FSEL
#define COMP_FORCE
#endif

#undef EXIT_STEP1
!#define EXIT_STEP1

subroutine ElectrostaticEnergy()
  !
  ! Compute energy, potential and force
  !
!$ use omp_lib
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: ws, nL, tri, nSzTri, aC, aQ, bC, bQ, cC, cQ, cubeEdge, vertexMin, bBoxesStart, sBoxLLstart, &
                         nearNeighbors, nNearNeighbors, nBoxesSList, maxBoxesSList, rankAtoms, dispAtoms, &
                         rankBoxes, dispBoxes, nDevicesNode, farNeighbors, nFarNeighbors, nBoxesTotal, &
                         nAtomsTotal, stdEnergy, memUsed, memMax, memMaxMyNode, memMaxAllNodes, &
                         maxNNperBox, maxFNperBox, fBoxes, &
#ifdef MPI
                         MPI_COMM_RMGROUP, MPI_COMM_NODE, MPI_COMM_NZERO, &
#endif
                         boxRegExpansion, boxIrrExpansion, rBoxBottomStart, rBoxBottomEnd, &
                         segmentRank, segmentNBoxes, nBoxSegments, nExtraMultipoleBoxes, &
                         nBoxesRList, nRemoteFFBoxes, nParticlesRank, nParticlesRMRank, lComputeDirectSum, &
                         particleRanges, nLLBoxes_rmRank, leadRank, prtRank, nBoxesBefore_loLevel, &
                         nFloatingBoxes, sBoxesTopNum, sBoxesTopDispl, iAtomsBegin, iAtomsEnd, &
                         minTopRtreeSize, minTopRtreeRanks, minTopRtreeGroup, expShift, numBoxesSNRF, nodeSize, &
                         nRanks, rmSize, nNodes, debugEachMPIrank, nGroups, nSelectedLLboxesRank, nAtomsGroupSelected
  use common_constants, only: eConv, pConv, fConv, zero, zeroCoor, mSign, iGB, invLog8
  use fmm_plan, only: treeLevel, ffLevel, rmLevel, ddLevel, loLevel, sList, rList, segmentSBox
  use timer
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_f_pointer, c_associated
  use allocator
  implicit none
  !
#ifdef MPI
  !integer      :: status(MPI_STATUS_SIZE)
  type(MPI_Comm):: MPI_COMM_XXGROUP    ! cross-communicator that connects the sub-groups
  integer       :: mpiError, dcmplx_size, targetCount, targetRank, color, xxPrint, xxRank, xxSize 
  logical (kind=1) :: xxComm
  integer (kind=8) :: nMaxMessages, nTotMessagesRMG, nMaxMessagesWorld
#ifdef CPLR_NVFTN
  integer       :: window
  integer       :: status
#else
  type(MPI_Win) :: window
  type(MPI_Status) :: status
#endif
  integer (kind=MPI_ADDRESS_KIND) :: buffer_size, targetDisplacement, sourceDisplacement  ! MPI_ADDRESS_KIND == 8
  integer, allocatable :: nMessages(:), messageDisplacement(:)
  integer (kind=MPI_ADDRESS_KIND), allocatable :: bufferSizes(:)
  integer (kind=MPI_ADDRESS_KIND), allocatable :: message(:,:), imessage(:,:)  ! TODO can possibly be int4
#endif
  !
  integer :: iam, actual_nRanks, nOMPjobs, iDevice, nThreads, iThread, rmRank, iRank, jRank
  integer :: allocStat, success, iLevel, aLevel, rLevel, nbParent, nbParentEdge, nbParentCube, nbParentEdge_init
  integer (kind=4) :: ix,iy,iz, jx,jy,jz, kx,ky,kz, tjx,tjy,tjz
  double precision :: xShift, yShift, zShift, bx,by,bz
  double precision :: iOrigin(3), jOrigin(3), pOrigin(3), zOrigin(3)
  !double complex, allocatable :: tmpExpansion(:) !, shiftedExpansion(:), oShift(:)
#ifdef M2LROT
  double complex, allocatable :: scrExpansion(:), expalp(:)
#endif
  integer (kind=4) :: memFreeNode, memFreeNodeMin
  integer (kind=8) :: n1_8, n2_8, n3_8, n4_8, n5_8, n6_8, n7_8

! add debug print after the specified kernel
!#define DEBUG_STEP1
!#define DEBUG_STEP2
!#define DEBUG_STEP3
!#define DEBUG_STEP4
!#define DEBUG_PRINT_ENERGY

! print one-sided Get messages to stdout
! TODO cannot be used with MPI debug mode: need to replace prtRank with leadRank
!#define ONESIDED_PRINT

! apply fence after nFence one-sided messages
#define USE_ONESIDED_FENCE

! use MPI_Put instead of MPI_Get
!#define USE_ONESIDED_PUT

! do not modify preprocessor directives bellow this line
#if defined(ONESIDED_PRINT) || defined(USE_ONESIDED_PUT)
  #define ONESIDED_COLLECT
#endif
#if defined(ONESIDED_PRINT) && defined(USE_ONESIDED_PUT)
  ! cannot use ONESIDED_PRINT together with USE_ONESIDED_PUT since ONESIDED_PRINT overwrites the message buffer
  #undef USE_ONESIDED_PUT
#endif
!#define ONESIDED_COLLECT

#if defined(DEBUG_STEP1) || defined(DEBUG_STEP2) || defined(DEBUG_STEP3) || defined(DEBUG_STEP4)
  #define DEBUG_PRINT_MULIPOLES
  ! debug: print first few elements of multipole expansion
  ! TODO: implement dynamic allocation for outexp and prtexp buffers to accomodate arbitrary problem sizes
  integer, parameter :: nSzPrt = 5, nSzPrt1 = nSzPrt+1, nMaxRanks = 64, nMaxBoxesRank = 500  ! 0..nSzPrt elements
  integer, parameter :: szOut = nSzPrt1*nMaxBoxesRank*nMaxRanks
  integer, parameter :: szPrt = nSzPrt1*nMaxBoxesRank
  integer            :: nbr(1:2,0:nMaxRanks)   ! box start and box end for each rank
  double complex     :: outexp(0:szOut)        ! save buffer
  double complex     :: prtexp(0:szPrt)        ! print buffer
#endif
  double complex     :: shiftedExpansion(0:nSzTri), oShift(0:nSzTri)
  double precision, allocatable :: aPotential(:), aForce(:,:), bPotential(:), bForce(:,:)
#ifdef FSEL
  integer (kind=8), allocatable :: nAtomsGroupRank(:), nAtomsGroupGlobal(:)          ! number of atoms in the group
  double precision, allocatable :: boxFsum(:,:), refFsum(:,:), rForce(:,:), hostBoxFsum(:,:), hostRefFsum(:,:)
  double precision, allocatable :: boxGroupFsumRank(:,:),   refGroupFsumRank(:,:)    ! size of nGroups
  double precision, allocatable :: boxGroupFsumGlobal(:,:), refGroupFsumGlobal(:,:)  ! size of nGroups
#endif
  double precision :: sump, energy, energyStd, esum, memLocal, refEnergy, refEnergy_sum, fmmEnergy
  double precision :: pot, potA, potB, rij(3), forceA(3), forceB(3), bpot, br2, brij(3)
#if defined (DEBUG) && defined (DETAILED)
  double precision :: energy1, energy2, energy3, error
#endif
#ifdef DYNAMICLIST
  double complex, allocatable :: potentialExpansion(:) !, irrExpansion(:)
  integer :: pjx,pjy,pjz, nbEdge, nbEdgeMAxID, nbStart, nbParent, nbParentEdgeMaxID
  integer :: fax,fay,faz, fbx,fby,fbz, nax,nay,naz, nbx,nby,nbz, aLevel
  integer, allocatable :: boxList(:)
  integer :: idBox
#endif
  double precision :: edgeReciprocal
  double complex, allocatable :: regHarmonics(:,:)
  !
  double precision :: dL, dm, sumf(3), aprod, bprod
  double complex   :: cprod
  !
  integer          :: L, m, Lm, nLx, aL, aL1, aL2, LL, i, pdigits, xyz(3)
  integer (kind=8) :: i_8, iAtom, jAtom, aAtom, nBoxesBefore
  integer (kind=8) :: iBox, jBox, pBox, aBox, cBox, ptrBox, expStart, cBoxShift
  double precision :: x,y,z,x2y2,xDx2y2,yDx2y2,r2,z2lm1
  !
  integer          :: iMember, nMembers, family(loLevel), i0, j, k, aJ, aJL, aLm, aLJ, LL1
  integer (kind=4) :: bBefore, cBefore, pBefore, nChildrenToProcess, childBoxesBegin, childBoxesEnd
  integer (kind=4) :: cNodeStart, srcBox, dstBox, childBox, parentBox, firstBox, lastBox, nextBox
  integer (kind=4) :: nBoxesInSegment, pLevel, eBox
  double precision :: rsum, mSm, mS, ir, ir2, ir3
  double complex   :: csum, xpy
#ifndef CPLR_NVFTN
  double complex, pointer :: ptr(:)
#endif
  integer          :: nLL, nLL1, jL, km, jkx
  double precision :: memItem, memTotal, memH2D, memD2H, memDevTotal, maxmemDevTotal
  integer          :: found
  logical            :: lrefEnergyAvailable
  integer            :: nErrors, nErrors_sum, nFence
  integer (kind=8)   :: nMinMessages
  integer, parameter :: nhist = 16   ! number of elements in the histogram
  integer            :: rPhistogram(nhist), rFhistogram(nhist)  ! relative % error
  integer            :: aPhistogram(nhist), aFhistogram(nhist)  ! absolute error
  double precision   :: difPotMax,     difPotRel,     difForceMax,     difForceRel
  double precision   :: difPotMax_sum, difPotRel_sum, difForceMax_sum, difForceRel_sum
  double precision   :: rmsdP,     rmsdF,     averP,     averF
  double precision   :: rmsdP_sum, rmsdF_sum, averP_sum, averF_sum
  double precision   :: relPotential,  relForce
  double precision   :: difPotential,  difForce
  double precision   :: dsPotential, dsForce(3)
  double precision, dimension(nhist), parameter :: err = &
          (/ 1.0d-1, 1.0d-2,  1.0d-3,  1.0d-4,  1.0d-5,  1.0d-6,  1.0d-7,  1.0d-8, &
             1.0d-9, 1.0d-10, 1.0d-11, 1.0d-12, 1.0d-13, 1.0d-14, 1.0d-15, 1.0d-16 /)
  !
  ! Buffers to compute regular multipole expansion of E-boxes
  ! Intel Fortran compiler has problems with dynamically allocated arrays in GPU kernels
  ! Works for up to nL=40 when nSzTri+1 = 861
  integer (kind=4), parameter :: maxSupportedAngularMomentum = 40
  integer (kind=4), parameter :: szBuffer = (maxSupportedAngularMomentum+1)*(maxSupportedAngularMomentum+2)/2
  double complex     :: parentRegExpansion(0:szBuffer-1) ! buffer for 1 parent regular expansion
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
  double complex, allocatable :: hostRegExpansion(:) ! host copy of boxRegExpansion
  integer (kind=8)            :: szHostRegExp
#endif
  integer (kind=8)            :: szBoxRegExp
#ifdef FSUM
  double precision            :: rankFSum(3), totalFSum(3)
  rankFSum = 0.0d0
#endif
  !
  ! do not allow nL > 30 to avoid buffer overflow for parentRegExpansion(0:szBuffer-1)
  if(nL > maxSupportedAngularMomentum .and. ffLevel < rmLevel) then
    if(iam == prtRank) then
      write(*,'(/a)') "Error in energy.F90: nL > maxSupportedAngularMomentum"
    endif
    stop 1
  endif
  !
#ifndef MAX_PERF
#ifdef MPI
  call MPI_Barrier(MPI_COMM_RMGROUP, mpiError)
#endif
#endif

#ifdef TIMER
  call TimerInit(2) ! measure time spent in FMM computation
#endif
  !
  rPhistogram = 0        ! zero initialize the histogram
  rFhistogram = 0        ! zero initialize the histogram
  aPhistogram = 0        ! zero initialize the histogram
  aFhistogram = 0        ! zero initialize the histogram
  !
#ifdef MPI
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, actual_nRanks, mpiError)  ! to recognize the memory-usage simulator mode
  call MPI_Comm_rank(MPI_COMM_RMGROUP, rmRank, mpiError)
#else
  iam    = 0
  nRanks = 1
  rmRank = 0
  rmSize = 1
  actual_nRanks = 1
#endif
  !
#ifdef DEBUG_PRINT_MULIPOLES
  ! TODO: Use these limits for dynamic allocation
  if(iam == prtRank .and. nMaxRanks < nRanks) error stop "Error: nMaxRanks < nRanks"
  if(nMaxBoxesRank*nMaxRanks < rBoxBottomEnd-rBoxBottomStart+1 .or. &
     nMaxBoxesRank*nMaxRanks < nBoxesSList-sBoxLLstart) error stop "Error: insufficient buffer size for debug print"
#endif
  !
  iThread  = 0
  nThreads = 1
!!$ call omp_set_num_threads(min0(omp_get_max_threads(), 8))  ! limit the max number of threads to 8
!$ nThreads = omp_get_max_threads()
#ifdef OFFL_GPU
  iDevice = mod(iam,nDevicesNode)
#endif

  ! multipole order, L
  if(iam == prtRank) then
    write(*,'(/a,i12)')    "Polynomial degree                  : ", nL
    write(*,'(a,i12)')     "Multipole size                     : ", nSzTri+1
    write(*,'( a,1x,i11)') "Number of OMP threads per rank     : ", nThreads
  endif

  ! track the memory allocation in this subroutine
  memLocal = 0.0d0

  ! memory allocation report (TODO: inaccurate, need to revise)
  if(iam == prtRank) &
  write(*,'(/a)') "Memory Allocation Report:"
  memTotal    = 0.0d0
  memH2D      = 0.0d0
  memD2H      = 0.0d0
  memDevTotal = 0.0d0

  ! box structure
  ! allocated and initialized on CPU and used read-only on GPU
  memItem  = dble(numBoxesSNRF) * 44 * iGB  ! box -> cube
  memTotal = memTotal + memItem
  memH2D   = memH2D   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "rList(type B-list 44 * numBoxesSNRF)                = ", memItem, " GB"

  ! box structure
  ! allocated and initialized on CPU and used read-only on GPU
  memItem  = dble(maxBoxesSList) * 80 * iGB
  memTotal = memTotal + memItem
  memH2D   = memH2D   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "sList(type Struct Tree 80 * maxBoxesSList)          = ", memItem, " GB"

  if(lComputeDirectSum) then
    ! needed only for accuracy validation by using direct summation on GPU
    ! all particle coordinates and weights (charges)
    ! should be present both on CPU and GPU
    memItem  = 4 * dble(nAtomsTotal) * 8 * iGB
    memTotal = memTotal + memItem
    if(iam == prtRank) write(*,'(a,f10.2,a)') &
                         "a[C|Q](4,nAtomsTotal)                               = ", memItem, " GB"
  endif

  ! rank assigned particle coordinates and weights (charges)
  ! should be present both on CPU and GPU
  memItem  = 4 * dble(nParticlesRank) * 8 * iGB
  memTotal = memTotal + memItem
  memH2D   = memH2D   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "b[C|Q](4,nParticlesRank)                            = ", memItem, " GB"

  ! rank assigned particle coordinates and weights (charges)
  ! should be present both on CPU and GPU
  ! the array cC works in pair with particleRanges(1:2,0:nLLBoxes_rmRank-1)
  memItem  = (dble(nParticlesRMRank)*4*8 + dble(nLLBoxes_rmRank)*2*4) * iGB
  !memTotal = memTotal + memItem
  !memH2D   = memH2D   + memItem
  !memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "c[C|Q](4,nParticlesRMRank)                          = ", memItem, " GB"

#if defined(FETCH_POTFORCE)
  ! potentian and force in the position of each particle
  ! should be present only on CPU on rank 0
  memItem = 4 * dble(nAtomsTotal) * 8 * iGB
  memTotal = memTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "a[Potential|Force](4,nAtomsTotal)                   = ", memItem, " GB"
#endif

#if defined(COMP_FORCE)
  ! potentian and force in the position of each particle
  ! allocated on GPU only
  ! could be present both on CPU and GPU
  memItem = 4 * dble(maxval(rankAtoms)) * 8 * iGB
  !memTotal = memTotal + memItem
  !memH2D   = memH2D   + memItem
  !memD2H   = memD2H   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "[bPotential|bForce](4,rankAtoms(iam))           GPU = ", memItem, " GB"
#else
  ! potentian in the position of each particle
  ! allocated on GPU only
  ! could be present both on CPU and GPU
  memItem = dble(maxval(rankAtoms)) * 8 * iGB
  !memTotal = memTotal + memItem
  !memH2D   = memH2D   + memItem
  !memD2H   = memD2H   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "bPotential(rankAtoms(iam))                      GPU = ", memItem, " GB"
#endif

  ! box regular harmonics
  ! allocated both on CPU (hostRegExp) and on GPU
  memItem  = dble(nSzTri+1)*dble(nBoxesRList+nRemoteFFBoxes) * 16 * iGB
  memTotal = memTotal + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "boxRegExp(nSzTri*(nBoxesRList+nRemoteFFBoxes))  GPU = ", memItem, " GB"

  ! sList irregular harmonics (boxIrrExpansion) on an MPI rank
  ! allocated on GPU only
  memItem  = dble(nSzTri+1)*dble(nBoxesSList+1) * 16 * iGB
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "boxIrrExpansion(0:nBoxesSList*(nSzTri+1))       GPU = ", memItem, " GB"

#ifdef FSEL
  ! boxFsum is the sum of atomic forces for each LL-box in S-list
  ! boxFsum is allocated in GPU only
  ! hostBoxFsum is allocated in CPU only
  memItem  = dble(nBoxesSList-sBoxLLstart) * 3 * 8 * iGB
  memTotal = memTotal + memItem        ! add to the CPU count
  memH2D   = memH2D   + memItem
  memD2H   = memD2H   + memItem
  memDevTotal = memDevTotal + memItem  ! add to the GPU count
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "boxFsum(1:3,0:nBoxesSList-sBoxLLstart-1)        GPU = ", memItem, " GB"
  ! boxGroupFsumRank and boxGroupFsumGlobal are allocated in CPU only
  memItem  = dble(2*nGroups) * 3 * 8 * iGB
  memTotal = memTotal + memItem        ! add to the CPU count
  ! nAtomsGroupRank and nAtomsGroupGlobal are allocated in CPU only
  memItem  = dble(2*(nGroups+1)) * 8 * iGB
  memTotal = memTotal + memItem        ! add to the CPU count

  if(lComputeDirectSum) then
  ! refFsum is the sum of reference atomic forces for each LL-box in S-list
  ! allocated on GPU only
  ! hostRefFsum is allocated on CPU only
  memItem  = dble(nBoxesSList-sBoxLLstart) * 3 * 8 * iGB
  memTotal = memTotal + memItem        ! add to the CPU count
  memH2D   = memH2D   + memItem
  memD2H   = memD2H   + memItem
  memDevTotal = memDevTotal + memItem  ! add to the GPU count
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "refFsum(1:3,0:nBoxesSList-sBoxLLstart-1)        GPU = ", memItem, " GB"
  ! refGroupFsumRank and refGroupFsumGlobal are  allocated in CPU only
  memItem  = dble(2*nGroups) * 3 * 8 * iGB
  memTotal = memTotal + memItem        ! add to the CPU count

  ! rForce is the reference atomic force for all LL-boxes in S-list
  ! need rForce to compute refFsum
  ! allocated on GPU only
  memItem = dble(maxval(rankAtoms)) * 3 * 8 * iGB
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "rForce(3,rankAtoms(iam))                        GPU = ", memItem, " GB"
  endif  ! lComputeDirectSum
#endif

  ! list of near-neighbor boxes
  ! initialized on CPU, used and GPU
  !memItem = dble(nNearNeighbors) * 4 * iGB
  memItem = dble(nNearNeighbors) * 8 * iGB
  memTotal = memTotal + memItem
  memH2D   = memH2D   + memItem
  memDevTotal = memDevTotal + memItem
  if(iam == prtRank) &
  write(*,'(a,f10.2,a)') "nearNeighbors(nNearNeighbors)                       = ", memItem, " GB"

  ! list of far-neighbor boxes
  ! initialized on CPU, used and GPU
  !memItem = dble(nFarNeighbors) * 4 * iGB
  memItem = dble(nFarNeighbors) * 8 * iGB
  memTotal = memTotal + memItem
  memH2D   = memH2D   + memItem
  memDevTotal = memDevTotal + memItem
  ! TODO can dynamically compute the FF box indices to avoid the host-device data transfer and memory storage
  if(iam == prtRank) then
    write(*,'(a,f10.2,a)') "farNeighbors(nFarNeighbors)                         = ", memItem, " GB"

    write(*,'(a,f10.2,a)') "host -> device data transfer    per device | total  = ", memH2D, " GB"
    write(*,'(a,f10.2,a)') "device -> host data transfer    per device | total  = ", memD2H, " GB"
  endif

#ifndef MAX_PERF
  maxmemDevTotal = memDevTotal
#ifdef MPI
  call MPI_Reduce(memDevTotal, maxmemDevTotal, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank) &
  write(*,'(a,f10.2,a,f10.2,a)') &
                           "device memory                   per device | total  = ", &
                           memDevTotal, " GB      max: ", maxmemDevTotal, " GB"

  ! Report free (host) memory on the compute node
  memFreeNode = stat_allmemory() / 1024_8 / 1024_8 ! reported in bytes and converted to Mb
  memFreeNodeMin = memFreeNode
#ifdef MPI
  if(nNodes > 1 .and. .not. debugEachMPIrank) then
    if(MPI_COMM_NZERO /= MPI_COMM_NULL) &
    call MPI_Reduce(memFreeNode, memFreeNodeMin, 1, MPI_INTEGER4, MPI_MIN, leadRank, MPI_COMM_NZERO, mpiError)
  endif
#endif
  if(iam == prtRank .or. debugEachMPIrank) write(*,'(a,i7,a,i7,a)') &
  "Free Host Memory in the nodes              Local:", memFreeNode, " Mb    Min:", memFreeNodeMin, " Mb"
#endif  ! MAX_PERF

  ! Compute the amount of memory to be allocated
  memLocal = memLocal &
#ifdef M2LROT
           + dble(nSzTri+1)*16         &           ! scrExpansion
           + dble(nL+1)*16             &           ! expalp
#endif
#ifdef FETCH_POTFORCE
           + dble(nAtomsTotal)*8       &           ! aPotential
           + dble(nAtomsTotal)*8*3     &           ! aForce
#endif
#ifdef DYNAMICLIST
           + dble(nSzTri+1)*16         &           ! potentialExpansion
           + dble(max0(maxNNperBox,maxFNperBox))*4 & ! boxList
#endif
!           + dble(rankAtoms(iam))*8    &           ! bPotential
!           + dble(rankAtoms(iam))*8*3  &           ! bForce
!
! Do not account for memory allocation of boxRegExpansion in the host if allocated in the device
#if !defined(USE_DEVPTR)
           + dble(nSzTri+1) * (dble(nBoxesRList) + dble(nRemoteFFBoxes)) * 16  &  ! boxRegExpansion
#endif
           + dble(szBuffer*16)  ! static host memory used by parentRegExpansion

  ! szBoxRegExp is the max index in boxRegExpansion: allocate(0:szBoxRegExp)
  ! allocate(boxRegExpansion(0:(nSzTri+1)*(nBoxesRList+nRemoteFFBoxes)-1))
  n1_8 = nSzTri+1
  n2_8 = nBoxesRList
  n3_8 = nRemoteFFBoxes
  szBoxRegExp = n1_8 * (n2_8 + n3_8) - 1_8

#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
  ! szHostRegExp is the max index in hostRegExpansion
  ! merge(1, 0, rmRank >= minTopRtreeRanks) is equivalent to merge(0, -1, rmRank >= minTopRtreeRanks) + 1
  ! the following (small memory allocation) works for E-boxes only
  !szHostRegExp = (1 + nExtraMultipoleBoxes(rmRank) + minTopRtreeSize + merge(1, 0, rmRank >= minTopRtreeRanks)) &
  !         * (nSzTri+1) - 1  ! allocate(0:szHostRegExp)
  szHostRegExp = szBoxRegExp  ! make a full copy of boxRegExpansion in the host memory
  memLocal = memLocal + (dble(1_8 + szHostRegExp)) * 16  ! include the master box in hostRegExpansion
#endif

  ! reset the value of memMax to track the memory usage in the current function
  memMax = memUsed + memLocal
#ifndef MAX_PERF
#ifdef MPI
  ! determine the max memory requirement for each node
  call MPI_AllReduce(memMax, memMaxMyNode, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NODE, mpiError)
  if(iam == prtRank)  call MemNodePrint("Energy  myNodeMax", memLocal, memUsed+memLocal, memMax, memMaxMyNode)

  ! determine the max memory requirement per node across all nodes
  if(nNodes > 1 .and. MPI_COMM_NZERO /= MPI_COMM_NULL) then
    call MPI_Reduce(memMaxMyNode, memMaxAllNodes, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_NZERO, mpiError)
    if(iam == 0) call MemNodePrint("Energy allNodeMax", memLocal, memUsed+memLocal, memMax, memMaxAllNodes)
  endif
#else
  if(iam == prtRank) call MemPrint("Energy", memLocal, memUsed+memLocal, memMax, nodeSize)
#endif
#endif  ! MAX_PERF

  ! allocate memory
  allocate( &  !oShift(0:nSzTri), shiftedExpansion(0:nSzTri), irrExpansion(0:nSzTri), tmpExpansion(0:nSzTri), &
#ifdef M2LROT
           scrExpansion(0:nSzTri),        &
           expalp(0:nL),                  &
#endif
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
           hostRegExpansion(0:szHostRegExp) , &  ! index 0 is for the master box
#endif
#ifdef FETCH_POTFORCE
           aPotential(nAtomsTotal),       &
           aForce(3,nAtomsTotal),         &
#endif
#ifdef DYNAMICLIST
           potentialExpansion(0:nSzTri),  &
           boxList(max0(maxNNperBox,maxFNperBox)), &   ! use the higher value
#endif
!           bPotential(rankAtoms(iam)),    &
!           bForce(3,rankAtoms(iam)),      &
           stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)

  ! boxRegExpansion holds regular multipole moments of nBoxesRList boxes composing
  !   R-list and computed in rmRank in RMG in Step 1 and then in Step 2. The R-list
  !   includes virtual master box, followed by extra multipole E-boxes, which appear
  !   when rmLevel > ffLevel, and continued by the boxes of R-tree. The R-tree box
  !   numbering begins from the rmLevel boxes assigned to rmRank in RMG and hierarchically
  !   (binary-tree in 3D) continues down to loLevel where it ends with R-boxes, which hold
  !   particles. Step 1 performs the computation of multipole moments of R-boxes. Step 2
  !   computes multipole moments of other boxes (virtual master box remains unassigned).
  ! In addition to the boxes of R-list, boxRegExpansion includes nRemoteFFBoxes
  !   boxes (located on other ranks in RMG) required for FF (far-field) computation
  !   in Step 3. The regular multipoles of these boxes must be retrieved from the
  !   remote ranks in RMG prior the start of Step 3.
  ! boxRegExpansion begins with nBoxesRList boxes of R-list followed by nRemoteFFBoxes.
  !   Each of these two groups of boxes is separately sequentially numbered from top
  !   to bottom of the hierarchy.

#ifdef USE_DEVPTR
  ! make the array boxRegExpansion resident (allocated) in GPU
  !!$omp allocate allocator(omp_target_device_mem_alloc)
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):boxRegExpansion)
#endif
  allocate(boxRegExpansion(0:szBoxRegExp), stat=allocStat)
  !allocate(boxRegExpansion(0:(nSzTri+1)*(nBoxesRList+nRemoteFFBoxes)-1), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
#ifdef OFFL_CPU
  ! for CPU-only computation
  boxRegExpansion = dcmplx(0.0d0, 0.0d0)
#endif

  energy     = 0.0d0
#ifdef FETCH_POTFORCE
  aPotential = 0.0d0
  aForce     = 0.0d0
#endif
!  bPotential = 0.0d0
!  bForce     = 0.0d0
  !
  edgeReciprocal = 1.0d0 / cubeEdge

#ifdef TIMER
  call TimerInit(3)
#endif

  rLevel = loLevel - ddLevel  ! rLevel is a local hierarchy level for cBox

  ! get cBox from iBox
  cBoxShift = nBoxesBefore_loLevel + sBoxesTopDispl(iam) * 8_8**rLevel &
         - ((8_8**rLevel - 1_8)/7_8 * sBoxesTopNum(iam) + 1_8) - nFloatingBoxes

#if OFFL_GPU
  !$omp target data map(to: vertexMin, nBoxesBefore_loLevel, rLevel, tri, nSzTri, iAtomsBegin, iAtomsEnd, nFloatingBoxes) &
  !$omp map(to: nBoxesRList, nRemoteFFBoxes, cBoxShift) &
  !$omp map(alloc: oShift) device(iDevice) &
  !$omp map(to: cC, cQ, particleRanges) &
#ifdef FSUM
  !$omp map(from: rankFSum) &
#endif
#ifdef USE_DEVPTR
  !$omp use_device_addr(boxRegExpansion)
#else
  !$omp map(to: boxRegExpansion)
#endif
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
  ! do not need a barrier here unless it is enforced
#ifndef MAX_PERF
#if defined (MPI) && defined (USE_REDUCE)
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
#endif
  if(iam == prtRank) then
     write(*,'(/a,2(f6.1,a))') "Data Transfer: Host to Device       min: ", &
                               TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif
  !
  ! FMM step 1    Compute multipole expansion for the lowest-level boxes
  !
#ifdef OFFL_GPU
  ! zero-initialize the array in GPU
#ifdef CPLR_NVFTN
  !$omp target teams distribute parallel do device(iDevice) 
#else
  !$omp target teams distribute parallel do has_device_addr(boxRegExpansion) device(iDevice) 
#endif
  do i = 0, (nSzTri+1)*(nBoxesRList+nRemoteFFBoxes)-1
    boxRegExpansion(i) = dcmplx(0.0d0, 0.0d0)
  enddo
#else
  ! CPU-only computation
  boxRegExpansion = dcmplx(0.0d0, 0.0d0)
#endif
  !
  !expShift = box(rBoxBottomStart)%expStart
  nbParentEdge_init = 2 ** (loLevel - 1)
  nbParentCube = nbParentEdge_init * nbParentEdge_init * nbParentEdge_init

  ! catch potential overflow
  n1_8 = nbParentEdge_init       ! int*4 math
  n2_8 = 2_8 ** (loLevel - 1_8)  ! int*8 math
  n3_8 = nbParentCube
  n4_8 = n2_8 * n2_8 * n2_8
  if(n1_8 /= n2_8 .or. n3_8 /= n4_8) then
    if(iam == prtRank)  write(*,'(/a)') "Error in ElectrostaticEnergy(): int*4 overflow with nbParentCube"
    error stop 1
  endif

  ! Compute regular multipole expansion of R-boxes at loLevel assigned to rmRank in RMG
  ! Each rmRank holds a unique R-tree of boxes
  ! R-tree begins at rmLevel from boxes assigned to rmRanks in RMG domain decomposition
  ! R-tree marks a unique set of R-boxes at loLevel based on tree expansion
  ! If rmLevel == ffLevel, the R-tree is complete and includes a complete parental
  !   lineage starting from the rmRank selected boxes at ffLevel down to those at loLevel
  ! If rmLevel /= ffLevel, the existing R-boxes at loLevel will form an incomplete set
  !   insufficient to compute regular multipoles of the boxes above rmLevel upto ffLevel
  !
#ifdef OFFL_CPU
  !$omp parallel do     &
#elif defined(OFFL_GPU) 

#ifdef USE_ARRAY_REDUCTION
  !$omp target teams distribute &
#else
  !$omp target teams distribute parallel do &
#endif

#ifdef DEBUG_STEP1
  !$omp map(from: outexp) &
#endif

#ifndef CPLR_NVFTN
  !$omp has_device_addr(boxRegExpansion) &
#endif
  !$omp map(to:nbParentEdge_init,expShift) device(iDevice) &
#endif
  !$omp default(shared) &
  !$omp private(jBox,xyz,nbParentEdge,aLevel,nbParentCube,pBox,tjx,tjy,tjz, iOrigin) &
#ifndef CPLR_NVFTN
  !$omp private(ptr) &
#endif
  !$omp private(expStart, iAtom, x,y,z, r2, oShift, L,z2lm1,aL,aL1,aL2,LL,m)
  do iBox = rBoxBottomStart, rBoxBottomEnd  ! loop over R-boxes assigned to rmRank and located at loLevel
    ! iBox is a global index of R-box from R-tree at loLevel
    !
    ! get ix, iy, iz from iBox
    jBox = iBox - nBoxesBefore_loLevel   ! get the local index of iBox for iLevel
    xyz = 0
    ! the number of primitive boxes in the parent box along the edge
    nbParentEdge = nbParentEdge_init
    do aLevel = 1, loLevel
      ! the number of primitive (lowest-level) boxes in the parent box
      nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

      ! obtain the parent box index that has jBox (iBox) as its child
      pBox = jBox / (nbParentEdge * nbParentEdge * nbParentEdge)
      tjz =  pBox / 4
      tjy = (pBox - tjz * 4) / 2
      tjx =  pBox - tjz * 4 - tjy * 2
      xyz(1) = xyz(1) + tjx * nbParentEdge
      xyz(2) = xyz(2) + tjy * nbParentEdge
      xyz(3) = xyz(3) + tjz * nbParentEdge

      jBox = mod(jBox, nbParentEdge * nbParentEdge * nbParentEdge)
      nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
    enddo
    !
    ! find center of the box
    iOrigin(1) = vertexMin(1) + (dble(xyz(1)) + 0.5d0) * cubeEdge
    iOrigin(2) = vertexMin(2) + (dble(xyz(2)) + 0.5d0) * cubeEdge
    iOrigin(3) = vertexMin(3) + (dble(xyz(3)) + 0.5d0) * cubeEdge
    !iOrigin(1) = vertexMin(1) + (dble(box(iBox)%xyz(1)) + 0.5d0) * cubeEdge
    !iOrigin(2) = vertexMin(2) + (dble(box(iBox)%xyz(2)) + 0.5d0) * cubeEdge
    !iOrigin(3) = vertexMin(3) + (dble(box(iBox)%xyz(3)) + 0.5d0) * cubeEdge

    ! use the knowledge that boxes are arranged sequentually without gaps
    ! expShift provides the address shift in boxRegExpansion
    ! it points to the place where to write expansions for D boxes
    ! the loLevel R-boxes reside at the end of the R-tree so expShift is not zero
    expStart = expShift + (iBox - rBoxBottomStart) * (nSzTri + 1)
#ifndef CPLR_NVFTN
    ptr(0:nSzTri) => boxRegExpansion(expStart:expStart+nSzTri)
#endif

    ! compute multipole expansion for each particle and sum it up in the box expansion
#if defined(USE_ARRAY_REDUCTION) && !defined(CPLR_NVFTN)
    !$omp parallel do reduction (+:ptr) &
    !$omp default(shared) &
    !$omp private(iAtom, x,y,z, r2, oShift, L,z2lm1,aL,aL1,aL2,LL,m)
#endif
    !do iAtom = box(iBox)%atomStart,  box(iBox)%atomEnd
    !do iAtom = box(iBox)%atomStart-box(rBoxBottomStart)%atomStart+1,  box(iBox)%atomEnd-box(rBoxBottomStart)%atomStart+1
    do iAtom = particleRanges(1,iBox-rBoxBottomStart), particleRanges(2,iBox-rBoxBottomStart)

      ! shift particle coordinates to the new coordinate system
      !x  = bC(1,iAtom) - iOrigin(1)
      !y  = bC(2,iAtom) - iOrigin(2)
      !z  = bC(3,iAtom) - iOrigin(3)
      x  = cC(1,iAtom) - iOrigin(1)
      y  = cC(2,iAtom) - iOrigin(2)
      z  = cC(3,iAtom) - iOrigin(3)
      r2 = x*x + y*y + z*z

      ! initialize first few regular harmonics
      ! (0,0) = (L,m)
      oShift(0) = dcmplx(1.0d0, 0.0d0)

      ! (1,0)
      oShift(1) = dcmplx(z, 0.0d0)
      ! (1,1)
      oShift(2) = dcmplx(-0.5d0*x, 0.5d0*y)

      ! obtain higher-degree harmonics from recurrence relations
      do L = 2, nL
        z2lm1 = z * dble(2*L-1)
        aL  = tri(L)
        aL1 = tri(L-1)
        aL2 = tri(L-2)
        LL  = L*L

        ! (L,0) ... (L,L-2)   use first recurrece relation
        do m = 0, L-2
          oShift(aL+m) = (z2lm1 * oShift(aL1+m) - r2 * oShift(aL2+m)) / dble(LL-m*m)
        enddo

        ! (L,L-1)             use second recurrence relation
        oShift(aL+L-1) = z * oShift(aL1+L-1)

        ! (L,L)               use third recurrence relation
        oShift(aL+L) = oShift(aL1+L-1) * oShift(2) / dble(L)
      enddo

      ! accumulate multipole expansion
      !boxRegExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) = &
      !boxRegExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) + oShift(0:nSzTri) * bQ(iAtom)
      !ptr(0:nSzTri) = ptr(0:nSzTri) + oShift(0:nSzTri) * bQ(iAtom)
#ifdef CPLR_NVFTN
      boxRegExpansion(expStart:expStart+nSzTri) = boxRegExpansion(expStart:expStart+nSzTri) + oShift(0:nSzTri) * cQ(iAtom)
#else
      ptr(0:nSzTri) = ptr(0:nSzTri) + oShift(0:nSzTri) * cQ(iAtom)
#endif
    enddo
#ifdef DEBUG_STEP1
    ! save first nSzPrt1 elements of multipole expansion of each LL box
    outexp((iBox-rBoxBottomStart)*nSzPrt1 : (iBox-rBoxBottomStart)*nSzPrt1 + nSzPrt) = ptr(0:nSzPrt)
#endif
  enddo

#ifdef OFFL_GPU
  !$omp target exit data map(delete: cC, cQ, particleRanges) device(iDevice)
  !$omp target enter data map(to: rList) device(iDevice)
#endif

#ifdef DEBUG_STEP1
  ! print first few elements of the regular multipole expansion of the lowest-level boxes
  call PrintExpansion(rBoxBottomStart, rBoxBottomEnd, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1, nbr, outexp, prtexp)
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 1: Timer: Multipole comp   min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER

#ifdef EXIT_STEP1
  ! deallocate memory and return
#ifdef OFFL_GPU
  ! deallocate target offloading
  !$omp end target data
#endif
  ! deallocate host memory
  success = 0
  deallocate( &  !oShift, shiftedExpansion, irrExpansion, tmpExpansion, &
#ifdef M2LROT
             scrExpansion, &
             expalp, &
#endif
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
             hostRegExpansion, &
#endif
#if defined(FETCH_POTFORCE)
             aPotential, &
             aForce, &
#endif
#ifdef DYNAMICLIST
             potentialExpansion, &
             boxList, &
#endif
             stat=allocStat)
  return
#endif ! EXIT_STEP1

!  ! broadcast the computed multipole expansions
!#if defined (MPI) && defined(USE_ALLGATHERV)
!  if(rmSize > 1) call MPI_Allgatherv(boxRegExpansion(dispBoxes(rmRank)), rankBoxes(rmRank), MPI_DOUBLE_COMPLEX, &
!                      boxRegExpansion, rankBoxes, dispBoxes, MPI_DOUBLE_COMPLEX, MPI_COMM_RMGROUP, mpiError)
!#endif
!#if defined (MPI) && !defined(USE_ALLGATHERV)
!  call MPI_Type_size(MPI_DOUBLE_COMPLEX, dcmplx_size)
!  buffer_size = ((nSzTri+1)*(nBoxesTotal+1) + 1) * dcmplx_size
!  call MPI_Win_Create(boxRegExpansion, buffer_size, dcmplx_size, MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
!  call MPI_Win_fence(0, window)
!  do targetRank = 0, rmSize-1
!    if(targetRank /= rmRank) then
!      call MPI_Put(boxRegExpansion(dispBoxes(rmRank)), rankBoxes(rmRank), MPI_DOUBLE_COMPLEX, &
!              targetRank, dispBoxes(rmRank), rankBoxes(rmRank), MPI_DOUBLE_COMPLEX, window, mpiError)
!    endif
!  enddo
!  call MPI_Win_fence(0, window)
!  call MPI_Win_free(window, mpiError)
!#endif

  !!boxRegExpansion(0:(nSzTri+1)*(nBoxesTotal+1))
  !iBox = treeLevel(loLevel)%bstart
  !call print_complex(boxRegExpansion(box(iBox)%expStart), nSzTri+1, "step 1")
  !iBox = treeLevel(loLevel)%bend
  !call print_complex(boxRegExpansion(box(iBox)%expStart), nSzTri+1, "step 1")

  zOrigin = 0.0d0
  !
  ! FMM step 2    M2M translation
  !
  ! Accumulate children regular multipole expansions in their respective parents
  ! Only boxes at the FF level and down need to have their multipole expansion

#ifdef OFFL_GPU
  !$omp target enter data map(alloc: shiftedExpansion) device(iDevice) &
  !$omp map(to: treeLevel, bC, bQ, farNeighbors, nearNeighbors) &
  !$omp map(to: sList, sBoxesTopNum, sBoxesTopDispl, nFloatingBoxes, minTopRtreeSize, minTopRtreeRanks)
#endif
  !TODO do not need minTopRtreeSize, minTopRtreeRanks in the device

  !do iLevel = loLevel-1, 0, -1        ! this is a general case, a bit more expensive
  !do iLevel = loLevel-1, ffLevel, -1  ! this is sufficient for non-periodic case
  do iLevel = loLevel-1, rmLevel, -1   ! need to handle (ffLevel..rmLevel-1) separately
    ! go up in the hierarchy level - traverse the R-tree
    ! perform RMG-style domain decomposition on the fly for each iLevel to generate R-type parent box ids
    ! iLevel is the parent hierarchy level for iBox (parent sums up the children's expansions)
    ! determine parent box ids for iLevel based on the known child box id at iLevel+1
    ! the first loop iteration determines parent boxes at loLevel+1 from R-boxes at loLevel
    ! the loops overwrite rBoxBottomStart and rBoxBottomEnd so the initial value becomes lost
    !
    ! cBefore applies to the child level in R-list (iLevel + 1)
    ! (cBefore + rBoxBottomStart) is the number of boxes in the R-list before the child level. For this
    !   reason, it is also the local address of the first box in the range rBoxBottomStart..rBoxBottomEnd
    ! cBefore - (subtracted now) rBoxBottomStart + (will be added later) rBoxBottomStart
    !   converts the global index (rBoxesStart..rBoxesEnd) to a local index for R-list
    ! the box index in R-list provides a pointer to boxRegExpansion(0:(nSzTri+1)*(nBoxesRList+nRemoteFFBoxes)-1)
    ! TODO cBefore should be int*8
    ! first next two lines are the number of boxes in the R-list before the child level
    cBefore = 1 + nExtraMultipoleBoxes(rmRank) &   ! iLevel >= rmLevel
      + (minTopRtreeSize + merge(1, 0, rmRank >= minTopRtreeRanks)) * (8**(iLevel+1-rmLevel)-1)/7 &
      - rBoxBottomStart  ! this is a convenience offset to simplify further computations

    ! get starting and ending box indices at the parent level
    ! global index at parent level
    rBoxBottomStart = (rBoxBottomStart - (8_8**(iLevel+1_8) - 1_8) / 7_8) / 8_8 + (8_8**iLevel - 1_8) / 7_8 
    rBoxBottomEnd   = (rBoxBottomEnd   - (8_8**(iLevel+1_8) - 1_8) / 7_8) / 8_8 + (8_8**iLevel - 1_8) / 7_8

    nbParentEdge_init = 2 ** (iLevel - 1)        ! number of boxes at the parent level
    nBoxesBefore = (8_8 ** iLevel - 1_8) / 7_8   ! number of boxes before iLevel
    ! R-list in addition to R-tree includes extra multipole boxes when ffLevel /= rmLevel
    ! Use R-tree to translate the regular multipoles upward from loLevel to rmLevel and accumulate
    ! bBefore has the same meaning as cBefore but applies to the parent level in R-list
    bBefore = 1 + nExtraMultipoleBoxes(rmRank) &
      + (minTopRtreeSize + merge(1, 0, rmRank >= minTopRtreeRanks)) * (8_8**(iLevel-rmLevel)-1_8)/7_8 &
      - rBoxBottomStart  ! this is a convenience offset to simplify further computations
#ifdef OFFL_CPU
    !$omp parallel do     &
#elif  OFFL_GPU
    !!$omp target teams distribute parallel do device(iDevice) thread_limit(2) has_device_addr(boxRegExpansion) &
    !$omp target teams distribute parallel do device(iDevice) &
#ifndef CPLR_NVFTN
    !$omp has_device_addr(boxRegExpansion) &
#endif
    !$omp map(to: nbParentEdge_init, nBoxesBefore, bBefore, cBefore, iLevel) &
    !$omp map(to: rBoxBottomStart, rBoxBottomEnd) &

#ifdef DEBUG_STEP2
    !$omp map(from: outexp) &
#endif

#endif
    !$omp default(shared) &
    !$omp private(expStart, aBox, xyz, nbParentEdge, aLevel, nbParentCube, iBox, tjx,tjy,tjz) &
    !$omp private(pOrigin, shiftedExpansion, oShift) &
    !$omp private(x,y,z, r2, L, z2lm1, aL, aL1, aL2, LL, m, i0, aLm, j, aJ, aLJ, k)
    do pBox = rBoxBottomStart, rBoxBottomEnd
    !do pBox = treeLevel(iLevel)%bStart, treeLevel(iLevel)%bEnd
      ! uses global box index
      expStart = (bBefore+pBox)*(nSzTri+1)  ! pointer to parent's regular multipole expansion

        ! get box grid coordinates ix, iy, iz from the global box index pBox
        aBox = pBox - nBoxesBefore   ! get the local index of pBox for iLevel
        xyz = 0
        ! the number of primitive boxes in the parent box along the edge
        nbParentEdge = nbParentEdge_init
        do aLevel = 1, iLevel
          ! the number of primitive (lowest-level) boxes in the parent box
          nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

          iBox = aBox / (nbParentEdge * nbParentEdge * nbParentEdge)   ! the parent box index that holds aBox (pBox)
          tjz =  iBox / 4
          tjy = (iBox - tjz * 4) / 2
          tjx =  iBox - tjz * 4 - tjy * 2
          xyz(1) = xyz(1) + tjx * nbParentEdge
          xyz(2) = xyz(2) + tjy * nbParentEdge
          xyz(3) = xyz(3) + tjz * nbParentEdge

          aBox = mod(aBox, nbParentEdge * nbParentEdge * nbParentEdge)
          nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
        enddo

        pOrigin(1) = (dble(xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
        pOrigin(2) = (dble(xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
        pOrigin(3) = (dble(xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge

        !pOrigin(1) = (dble(box(pBox)%xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
        !pOrigin(2) = (dble(box(pBox)%xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
        !pOrigin(3) = (dble(box(pBox)%xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge
  
        !box(pBox)%chargeExpansion = 0.0d0
  
        ! loop over 8 children
        xyz = xyz * 2
        aBox = (pBox - (8_8**iLevel - 1_8) / 7_8) * 8_8 + (8_8**(iLevel+1_8) - 1_8) / 7_8 
        do iBox = aBox, aBox + 7
        !do iBox = box(pBox)%childStart, box(pBox)%childEnd
          !
          ! Move the center from child to parent, from iBox to jBox (from - to)
          !iOrigin(1) = (dble(box(iBox)%xyz(1)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(1)
          !iOrigin(2) = (dble(box(iBox)%xyz(2)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(2)
          !iOrigin(3) = (dble(box(iBox)%xyz(3)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(3)
          !
          ! shift multipoles from child to parent and accumulate
#ifdef CPLR_NVFTN
          !shiftedExpansion = dcmplx(0.0d0, 0.0d0)
#else
          shiftedExpansion = dcmplx(0.0d0, 0.0d0)
#endif
#ifdef M2MSTD
          !!call RegularHarmonicsLocal(iOrigin(1), iOrigin(2), iOrigin(3), oShift)   ! compute helper expansion
          !!call M2M(box(iBox)%chargeExpansion, oShift, shiftedExpansion)
          !call RegularExpansion(iOrigin, 1.0d0, 1, zOrigin, tmpExpansion, oShift)

          ! shift particle coordinates to the new coordinate system
          x  = (dble(xyz(1) + mod(iBox-aBox,2))                  + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(1)
          y  = (dble(xyz(2) + (iBox-aBox - ((iBox-aBox)/4)*4)/2) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(2)
          z  = (dble(xyz(3) + (iBox-aBox)/4)                     + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(3)
          !x  = (dble(box(iBox)%xyz(1)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(1)
          !y  = (dble(box(iBox)%xyz(2)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(2)
          !z  = (dble(box(iBox)%xyz(3)) + 0.5d0) * treeLevel(iLevel+1)%edge - pOrigin(3)
          r2 = x*x + y*y + z*z

          ! initialize first few regular harmonics
          ! (0,0) = (L,m)
          oShift(0) = dcmplx(1.0d0, 0.0d0)

          ! (1,0)
          oShift(1) = dcmplx(z, 0.0d0)
          ! (1,1)
          oShift(2) = dcmplx(-0.5d0*x, 0.5d0*y)

          ! obtain higher-degree harmonics from recurrence relations
          do L = 2, nL
            z2lm1 = z * dble(2*L-1)
            aL  = tri(L)
            aL1 = tri(L-1)
            aL2 = tri(L-2)
            LL  = L*L

            ! (L,0) ... (L,L-2)   use first recurrece relation
            do m = 0, L-2
              oShift(aL+m) = (z2lm1 * oShift(aL1+m) - r2 * oShift(aL2+m)) / dble(LL-m*m)
            enddo

            ! (L,L-1)             use second recurrence relation
            oShift(aL+L-1) = z * oShift(aL1+L-1)

            ! (L,L)               use third recurrence relation
            oShift(aL+L) = oShift(aL1+L-1) * oShift(2) / dble(L)
          enddo

          !call M2Momp(box(iBox)%chargeExpansion, oShift, shiftedExpansion)
          !call M2Momp(boxRegExpansion(box(iBox)%expStart), oShift, shiftedExpansion)

          ! (0,0)
          !i0 = box(iBox)%expStart
          i0 = (cBefore + iBox)*(nSzTri + 1) ! i0 is the address of child expansion to be shifted
          shiftedExpansion(0) = boxRegExpansion(i0) * oShift(0)

          do L = 1, nL
            do m = 0, L
              aLm = tri(L) + m

              shiftedExpansion(aLm) = dcmplx(0.0d0, 0.0d0)

              ! j <= L-m
              do j = 0, L-m   !  0 <= m <= L-j
                aJ  = tri(j)
                aLJ = tri(L-j)

                ! k < 0;  m-k > 0;  m+|k| <= L-j
                do k = 1, min0(j,L-j-m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + dconjg(boxRegExpansion(i0+aJ+k)) * oShift(aLJ+m+k) * dble((-1)**k)
                enddo

                ! k = 0
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ) * oShift(aLJ+m)

                ! 0 < k <= m;  m-k >= 0
                do k = 1, min0(j,m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                enddo

                ! k > m;  m-k < 0
                do k = m+1, min0(j,L-j+m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) &
                          * dble((-1)**(k-m))
                enddo
              enddo

              ! j > L-m
              do j = L-m+1, L   !  m > L-j
                aJ  = tri(j)
                aLJ = tri(L-j)

                ! k < 0;  m+|k| > L-j   (L-j,m-k) = 0

                ! k = 0 ... m-L+j-1   (L-j,m-k) = 0

                ! k = m-L+j ... m
                do k = m-L+j, min0(m,j)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                enddo

                ! k = m+1 ... L-j+m
                do k = m+1, min0(L-j+m,j)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) &
                          * dble((-1)**(k-m))
                enddo
              enddo

            enddo  ! index m
          enddo  ! index L

#elif M2MROT
          !call M2Mrotate2(-iOrigin, expalp, box(iBox)%chargeExpansion, tmpExpansion, shiftedExpansion)
          call M2Mrotate2(-iOrigin, expalp, boxRegExpansion(box(iBox)%expStart), tmpExpansion, shiftedExpansion)
#else
          write(*,'(/a/a)') "Error: either -DM2MSTD or -DM2MROT have to be specified at compile time", &
                            "       to choose between the standard and rotation-based M2M routines" 
          error stop
#endif
          !box(pBox)%chargeExpansion = box(pBox)%chargeExpansion + shiftedExpansion
          !boxRegExpansion(box(pBox)%expStart:box(pBox)%expStart+nSzTri) = &
          !boxRegExpansion(box(pBox)%expStart:box(pBox)%expStart+nSzTri) + shiftedExpansion(:)
          boxRegExpansion(expStart:expStart+nSzTri) = &
          boxRegExpansion(expStart:expStart+nSzTri) + shiftedExpansion(:)
        enddo  ! loop over iBox
#ifdef DEBUG_STEP2
      ! save first nSzPrt1 elements of multipole expansion of boxes at ffLevel
      if(iLevel == ffLevel) &
      outexp((pBox-rBoxBottomStart)*nSzPrt1 : (pBox-rBoxBottomStart)*nSzPrt1 + nSzPrt) = &
         boxRegExpansion(expStart:expStart+nSzPrt)
#endif
    enddo   ! loop over pBox
  enddo  ! loop over iLevel

!  ! make sure all work in the sub-communicator is complete
!#ifdef MPI
!  call MPI_Barrier(MPI_COMM_RMGROUP, mpiError)
!#endif

#ifdef DEBUG_STEP2
  ! print first few elements of the regular multipole expansion of the boxes at ffLevel
  call PrintExpansion(rBoxBottomStart, rBoxBottomEnd, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1, nbr, outexp, prtexp)
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 2: Timer: M2M translation  min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    write(*,'(a,f6.1,a)')    "                   M2M d-mat comp      : ", timeM2Mdmat, "s"
    write(*,'(a,f6.1,a)')    "                   M2M z-translat      : ", timeM2Mdztr, "s"
    write(*,'(a,f6.1,a)')    "                   M2M d-mat rot       : ", timeM2Mdrot, "s"
  endif
  call TimerInit(3)
#endif ! TIMER

if(nRanks /= actual_nRanks) then
  if(iam == prtRank) then
    stop "End of the run in memory-usage simulator mode"
  else
    stop
  endif
endif


#ifdef MPI  ! begin the computation of E-boxes

    ! make a copy of GPU-resident boxRegExpansion on the host, hostRegExpansion
    ! copy device->host the entire data set once to avoid further data transfers
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
    !$omp target teams distribute parallel do device(iDevice) &
#ifndef CPLR_NVFTN
    !$omp has_device_addr(boxRegExpansion) &
#endif
    !$omp map(from: hostRegExpansion) &
    !$omp map(to: szHostRegExp)
    do i_8 = 0_8, szHostRegExp
      hostRegExpansion(i_8) = boxRegExpansion(i_8)
    enddo
#endif
  !
  ! FMM step 2e   M2M translation for E-boxes
  !
  ! Accumulate children regular multipole expansions in their respective parents
  ! Only boxes at the FF level and down need to have their multipole expansion
  ! E-boxes exist only in distributed memory (MPI) parallelizm mode as workload disbalance
  !
  !purpose: copy the multipole expansion of child boxes (at rmLevel) of E-boxes from device to host

  if(ffLevel < rmLevel) then
    ! iLevel is a parent level
    ! this is a special case of rmLevel == ffLevel+1
    ! children E-boxes always reside in rmLevel

    nChildrenToProcess = 0

    ! check if the rmRank has an extra multipole E-box
    ! one box from ffLevel per rmRank
    ! each box at ffLevel starts an E-tree; one E-tree per rmRank
    ! only first 8**ffLevel rmRanks get E-boxes
    if(rmRank < 8**ffLevel) then  ! copy local child E-boxes

      ! getting here implies that rmRank owns extra multipole E-boxes
      ! parentBox is an ffLevel parent box on the top of E-tree assigned to rmRank
      parentBox = rmRank  ! parentBox is a local box id at ffLevel

      ! some children of parentBox may be local to rmRank, whereas other children may be remote
      ! local children have their data (multipole expansions) stored in rmRank
      ! remote children have their data residing in remote ranks in RMG

      ! find the local children of parentBox
      ! determine the location of the children boxes

      ! local index (at rmLevel) of the first and last boxes of R-tree for rmRank
      firstBox = minTopRtreeSize * rmRank &
               + merge((rmRank-minTopRtreeRanks), 0, rmRank > minTopRtreeRanks)
      lastBox = firstBox + minTopRtreeSize + merge(0, -1, rmRank >= minTopRtreeRanks)

      ! child boxes occupy the bottom of the E-tree located at rmLevel
      ! child boxes in boxRegExpansion carry a copy of the data of respective boxes at rmLevel
      ! there are 8**(rmLevel-ffLevel) child boxes at rmLevel

      ! copy the data for child boxes local to rmRank from the rmLevel to the temporary buffer

      ! these are local box indices at rmLevel (child E-boxes local to rmRank)
      childBoxesBegin = max0(firstBox, parentBox * 8**(rmLevel-ffLevel)) ! cannot be less than firstBox
      childBoxesEnd   = min0(lastBox, (parentBox+1) * 8**(rmLevel-ffLevel) - 1) ! cannot exceed lastBox

      ! count the number of local child boxes
      nChildrenToProcess = max0(0, childBoxesEnd - childBoxesBegin + 1)

      ! get the data source box address for rmRank
      ! srcBox is the starting address of childBox in the buffer boxRegExpansion
      !   it is the number of boxes in boxRegExpansion before the current box
      ! the value of srcBox will be bugus (and will not be used) when nChildrenToProcess == 0
      srcBox = ( 1 &  ! master box
        + nExtraMultipoleBoxes(rmRank)  & ! all E-boxes including (child) buffer boxes
        + (childBoxesBegin - firstBox)) & ! internal offset of box id within R-tree in rmLevel
        * (nSzTri + 1)

      ! get the data destination box address (in the E-box buffer in boxRegExpansion) for rmRank
      ! the value of dstBox will be bugus (and will not be used) when nChildrenToProcess == 0
      dstBox = (1 &  ! master box
        + (8**(rmLevel-ffLevel) - 1) / 7  &  ! the number of (parent) E-boxes above rmLevel
        + (childBoxesBegin - parentBox*8**(rmLevel-ffLevel))) & ! child offset in the E-box storage
        * (nSzTri + 1)

      ! copy the child data locally in the rank rmRank if nChildrenToProcess > 0
      ! use the src data, which are already on the host
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
      do i = 0, (nSzTri + 1) * nChildrenToProcess - 1
        hostRegExpansion(dstBox + i) = hostRegExpansion(srcBox + i)
      enddo
#else
      do i = 0, (nSzTri + 1) * nChildrenToProcess - 1
        boxRegExpansion(dstBox + i) = boxRegExpansion(srcBox + i)
      enddo
#endif

      ! number of remote children E-boxes remaining to be processed
      nChildrenToProcess = 8**(rmLevel-ffLevel) - nChildrenToProcess

    endif  ! if(rmRank < 8**ffLevel)

    if(iam == prtRank) then
      if(rmRank < 8**ffLevel) then
        write(*,'(a,5i10)') "LocalE-boxes,iam,rmRank,nChildrenProcessed,srcBox,dstBox ", &
            iam,rmRank,8**(rmLevel-ffLevel)-nChildrenToProcess,srcBox,dstBox
      else
        write(*,'(a,5i10)') "LocalE-boxes,iam,rmRank,nChildrenProcessed,srcBox,dstBox ", &
            iam,rmRank,nChildrenToProcess,srcBox,dstBox
      endif
    endif

#ifdef TIMER
    timeConsumed = TimerCheck(3)
    timeConsumedMin = 0.0d0  ! there is no value in checking the min time
    timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#if defined (MPI) && defined (USE_REDUCE)
    ! TODO no need to synchronize the ranks outside RMGROUP
    call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
#endif
    if(iam == prtRank) then
      write(*,'(a,2(f6.1,a))') "FMM Step 2: Timer: M2M E-boxes loc  min: ", &
                               TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    endif
    call TimerInit(3)
#endif

    ! Copy the remaining data by using one-sided communication

    ! create one-sided window
    call MPI_Type_size(MPI_DOUBLE_COMPLEX, dcmplx_size, mpiError)
    n4_8 = dcmplx_size
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
    ! create the one-sided window in the host memory
    buffer_size = (szHostRegExp + 1_8) * n4_8
    call MPI_Win_Create(hostRegExpansion, buffer_size, dcmplx_size, &
                      MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
#else
    ! create the one-sided window in the device memory
    buffer_size = (szBoxRegExp + 1_8) * n4_8
    call MPI_Win_Create(boxRegExpansion, buffer_size, dcmplx_size, &
                      MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
    if(iam == prtRank) write(*,'(a)') "After MPI_Win_Create device"
#endif
    call MPI_Win_fence(0, window, mpiError)

    ! the following applies to those RMG ranks which have unfinished work
    if(nChildrenToProcess > 0) then
      nChildrenToProcess = 8**(rmLevel-ffLevel) - nChildrenToProcess
      ! find the remote children of parentBox

      ! these are child E-box indices local to rmLevel
      childBoxesBegin = parentBox * 8**(rmLevel-ffLevel)
      childBoxesEnd   = (parentBox+1) * 8**(rmLevel-ffLevel) - 1

      ! initialize the box id of the first child E-box at rmLevel
      childBox = childBoxesBegin

      do while (childBox <= childBoxesEnd)  ! loop over all child E-boxes
        ! start a new segment
        nBoxesInSegment = 1  ! initialize  the number of child E-boxes in the continuous segment

        ! get the rank id from the child E-box id at rmLevel
        iRank = merge(childBox/minTopRtreeSize, &
                      minTopRtreeRanks+(childBox-minTopRtreeGroup)/(minTopRtreeSize+1), &
                      childBox < minTopRtreeGroup)

        if(iRank /= rmRank) then
          ! this is a remote child E-box

          ! firstBox is the local index (within rmLevel) of the first R-tree box for iRank
          firstBox = minTopRtreeSize * iRank &
               + merge((iRank-minTopRtreeRanks), 0, iRank > minTopRtreeRanks)
         
          ! get the source box address for iRank
          ! srcBox is the starting address of childBox in the buffer boxRegExpansion
          !   it is the number of boxes in boxRegExpansion before the current box
          srcBox = (1 &  ! master box
                 + nExtraMultipoleBoxes(iRank) &  ! all floating boxes including (child) buffer boxes
                 + (childBox - firstBox))      &  ! internal offset of box id within R-tree in rmLevel
                 * (nSzTri+1)
         
          ! get the destination box address for rmRank
          dstBox = (1 &  ! master box
                 + (8**(rmLevel-ffLevel) - 1) / 7  &  ! the number of (parent) E-boxes above rmLevel
                 + (childBox - childBoxesBegin))   &  ! child offset in the E-box storage
                 * (nSzTri+1)

          ! check if the following child E-boxes form a continuous segment
          do nextBox = childBox+1, childBoxesEnd
            ! get the rank id from the child E-box id at rmLevel
            jRank = merge(nextBox/minTopRtreeSize, &
                    minTopRtreeRanks+(nextBox-minTopRtreeGroup)/(minTopRtreeSize+1), &
                    nextBox < minTopRtreeGroup)
            if(jRank == iRank) then
              nBoxesInSegment = nBoxesInSegment + 1
            else
              exit
            endif
          enddo
         
          targetDisplacement = dstBox                        ! save the data to this address
          targetCount        = (nSzTri+1) * nBoxesInSegment  ! size of the data
          targetRank         = iRank                         ! rank that holds the data
          sourceDisplacement = srcBox                        ! take the data from this address

          nChildrenToProcess = nChildrenToProcess + nBoxesInSegment

          ! save the data locally in boxRegExpansion(targetDisplacement)
          !call MPI_Get(boxRegExpansion(targetDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
          call MPI_Get(hostRegExpansion(targetDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
             sourceDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)

          !if(iam == prtRank) then
          !  write(*,'(a,5i10)') "AfterMPI_Get,iam,_iRank,nChildrenProcessed,srcBox,dstBox ", &
          !      iam,iRank,nChildrenToProcess,srcBox,dstBox
          !endif

        endif  ! iRank /= rmRank

        childBox = childBox + nBoxesInSegment  ! increment the index of child E-box
      enddo  ! do childBox = childBoxesBegin, childBoxesEnd

    endif  ! if(nChildrenToProcess > 0)

    if(iam == prtRank) then
      write(*,'(a,3i10)') "RemotE-boxes,iam,rmRank,nChildrenProcessed,srcBox,dstBox ", &
          iam,rmRank,nChildrenToProcess
    endif

    ! close the one-sided window
    call MPI_Win_fence(0, window, mpiError)
    call MPI_Win_free(window, mpiError)

#ifdef TIMER
    timeConsumed = TimerCheck(3)
    timeConsumedMin = 0.0d0
    timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#if defined (MPI) && defined (USE_REDUCE)
    call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
#endif
    if(iam == prtRank) then
      write(*,'(a,2(f6.1,a))') "FMM Step 2: Timer: M2M E-boxes rmt  min: ", &
                               TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    endif
    call TimerInit(3)
#endif

#ifdef TIMER_X
    timeConsumed = TimerCheck(3)
    timeConsumedMin = timeConsumed
    timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#if defined (MPI) && defined (USE_REDUCE)
    call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
    call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
#endif
    if(iam == prtRank) then
      write(*,'(a,2(f6.1,a))') "FMM Step 2: Timer: M2M E-boxes d2h  min: ", &
                               TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    endif
    call TimerInit(3)
#endif

    !purpose: perform M2M computation for E-boxes

    if(rmRank < 8**ffLevel) then
      ! this rmRank has extra multipole E-boxes
      ! parentBox accumulates children's multipole expansions

      ! pLevel is the parent level, which accumulates multipole expansion of its child boxes
      do pLevel = rmLevel-1, ffLevel, -1

      ! each rmRank deals with a group of parent boxes (one or more) at pLevel
      do eBox = 0, 8**(pLevel-ffLevel)-1  ! loop over 8**(pLevel-ffLevel) parent boxes

        ! eBox is a parent E-box group iterator at pLevel
        ! group iterator + local starting address of the group + global address shift for pLevel
        parentBox = eBox + rmRank * 8**(pLevel-ffLevel) + (8**pLevel-1)/7 ! global parent box id

        ! pBefore is the number of boxes before the parent E-box in boxRegExpansion
        pBefore = 1 &                         ! master box
           + (8**(pLevel-ffLevel) - 1) / 7  & ! number of boxes in E-tree above pLevel
           + eBox                             ! E-box iterator at pLevel

        ! childBox is the global id of the first child box located at the child level, pLevel+1
        ! it is also the total number of boxes before the first child box
        childBox = rmRank * 8**(pLevel+1-ffLevel) + (8**(pLevel+1)-1)/7 ! global id of the first child box

        ! cBefore is the number of boxes before the first child E-box in boxRegExpansion
        cBefore = 1 &                         ! master box
           + (8**(pLevel+1-ffLevel) - 1) / 7  ! number of boxes in E-tree above pLevel+1

        expStart = pBefore*(nSzTri+1)  ! pointer to parent's regular multipole expansion

        ! get box grid coordinates ix, iy, iz for parentBox
        aBox = rmRank   ! rmRank is the local index of parentBox for pLevel
        xyz = 0
        ! the number of primitive boxes in the parent box along the edge
        nbParentEdge = 2 ** (pLevel - 1)  ! number of boxes at the parent/parent level
        do aLevel = 1, pLevel
          ! the number of primitive (lowest-level) boxes in the parent box
          nbParentCube = nbParentEdge * nbParentEdge * nbParentEdge

          iBox = aBox / nbParentCube   ! the parent box index that holds aBox
          tjz =  iBox / 4
          tjy = (iBox - tjz * 4) / 2
          tjx =  iBox - tjz * 4 - tjy * 2
          xyz(1) = xyz(1) + tjx * nbParentEdge
          xyz(2) = xyz(2) + tjy * nbParentEdge
          xyz(3) = xyz(3) + tjz * nbParentEdge

          aBox = mod(aBox, nbParentCube)
          nbParentEdge = nbParentEdge / 2   ! go up to the grand-parent
        enddo

        pOrigin(1) = (dble(xyz(1)) + 0.5d0) * treeLevel(pLevel)%edge
        pOrigin(2) = (dble(xyz(2)) + 0.5d0) * treeLevel(pLevel)%edge
        pOrigin(3) = (dble(xyz(3)) + 0.5d0) * treeLevel(pLevel)%edge

        ! loop over 8 children
        xyz = xyz * 2
        aBox = childBox  ! it is a global box id 

        parentRegExpansion = dcmplx(0.0d0, 0.0d0)

        !$omp parallel do     &
        !$omp private(shiftedExpansion, x,y,z, r2, oShift, L, z2lm1, aL, aL1, aL2, LL, m) &
        !$omp private(i0, aLm, j, aJ, aLJ, k) reduction(+:parentRegExpansion)
        do childBox = 0, 7
          ! shift multipoles from child to parent and accumulate
          shiftedExpansion = dcmplx(0.0d0, 0.0d0)

          ! shift particle coordinates to the new coordinate system
          x  = (dble(xyz(1) + mod(childBox,2))               + 0.5d0) * treeLevel(pLevel+1)%edge - pOrigin(1)
          y  = (dble(xyz(2) + (childBox - (childBox/4)*4)/2) + 0.5d0) * treeLevel(pLevel+1)%edge - pOrigin(2)
          z  = (dble(xyz(3) + childBox/4)                    + 0.5d0) * treeLevel(pLevel+1)%edge - pOrigin(3)
          r2 = x*x + y*y + z*z

          ! initialize first few regular harmonics
          ! (0,0) = (L,m)
          oShift(0) = dcmplx(1.0d0, 0.0d0)

          ! (1,0)
          oShift(1) = dcmplx(z, 0.0d0)
          ! (1,1)
          oShift(2) = dcmplx(-0.5d0*x, 0.5d0*y)

          ! obtain higher-degree harmonics from recurrence relations
          do L = 2, nL
            z2lm1 = z * dble(2*L-1)
            aL  = tri(L)
            aL1 = tri(L-1)
            aL2 = tri(L-2)
            LL  = L*L

            ! (L,0) ... (L,L-2)   use first recurrece relation
            do m = 0, L-2
              oShift(aL+m) = (z2lm1 * oShift(aL1+m) - r2 * oShift(aL2+m)) / dble(LL-m*m)
            enddo

            ! (L,L-1)             use second recurrence relation
            oShift(aL+L-1) = z * oShift(aL1+L-1)

            ! (L,L)               use third recurrence relation
            oShift(aL+L) = oShift(aL1+L-1) * oShift(2) / dble(L)
          enddo

          ! (0,0)
          i0 = (cBefore + childBox)*(nSzTri + 1)  ! i0 is the address of child expansion
          !shiftedExpansion(0) = boxRegExpansion(i0) * oShift(0)
          shiftedExpansion(0) = hostRegExpansion(i0) * oShift(0)

          do L = 1, nL
            do m = 0, L
              aLm = tri(L) + m

              shiftedExpansion(aLm) = dcmplx(0.0d0, 0.0d0)

              ! j <= L-m
              do j = 0, L-m   !  0 <= m <= L-j
                aJ  = tri(j)
                aLJ = tri(L-j)

                ! k < 0;  m-k > 0;  m+|k| <= L-j
                do k = 1, min0(j,L-j-m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) &
                          + dconjg(hostRegExpansion(i0+aJ+k)) * oShift(aLJ+m+k) * dble((-1)**k)
                          !+ dconjg(boxRegExpansion(i0+aJ+k)) * oShift(aLJ+m+k) * dble((-1)**k)
                enddo

                ! k = 0
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + hostRegExpansion(i0+aJ) * oShift(aLJ+m)
                  !shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ) * oShift(aLJ+m)

                ! 0 < k <= m;  m-k >= 0
                do k = 1, min0(j,m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + hostRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                  !shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                enddo

                ! k > m;  m-k < 0
                do k = m+1, min0(j,L-j+m)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) &
                          + hostRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) * dble((-1)**(k-m))
                          !+ boxRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) * dble((-1)**(k-m))
                enddo
              enddo

              ! j > L-m
              do j = L-m+1, L   !  m > L-j
                aJ  = tri(j)
                aLJ = tri(L-j)

                ! k < 0;  m+|k| > L-j   (L-j,m-k) = 0

                ! k = 0 ... m-L+j-1   (L-j,m-k) = 0

                ! k = m-L+j ... m
                do k = m-L+j, min0(m,j)
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + hostRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                  !shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * oShift(aLJ+m-k)
                enddo

                ! k = m+1 ... L-j+m
                do k = m+1, min0(L-j+m,j)
                  !shiftedExpansion(aLm) = shiftedExpansion(aLm) + boxRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) &
                  shiftedExpansion(aLm) = shiftedExpansion(aLm) + hostRegExpansion(i0+aJ+k) * dconjg(oShift(aLJ+k-m)) &
                          * dble((-1)**(k-m))
                enddo
              enddo

            enddo  ! index m
          enddo  ! index L

          !boxRegExpansion(expStart:expStart+nSzTri) = &
          !boxRegExpansion(expStart:expStart+nSzTri) + shiftedExpansion(:)
          parentRegExpansion(0:nSzTri) = parentRegExpansion(0:nSzTri) + shiftedExpansion(0:nSzTri)
        enddo  ! loop over childBox
        hostRegExpansion(expStart:expStart+nSzTri) = parentRegExpansion(0:nSzTri)

      enddo  ! do eBox = 0, 8**(pLevel-ffLevel)-1

      enddo  ! do pLevel = rmLevel-1, ffLevel, -1

    endif  ! if(rmRank < 8**ffLevel)

#ifdef NEVER_RUNS
    ! create one-sided window
    call MPI_Type_size(MPI_DOUBLE_COMPLEX, dcmplx_size, mpiError)
    n1_8 = nSzTri
    n2_8 = nBoxesRList
    n3_8 = nRemoteFFBoxes
    n4_8 = dcmplx_size
    buffer_size = (n1_8+1_8)*(n2_8+n3_8) * n4_8
    call MPI_Win_Create(boxRegExpansion, buffer_size, dcmplx_size, &
                      MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
    call MPI_Win_fence(0, window, mpiError)

    sourceDisplacement = 0                 ! take the data from this address
    targetCount        = 2 * (nSzTri+1)    ! size of the data
    targetRank         = rmRank            ! rank that receives the data
    targetDisplacement = 0                 ! put the data ito this address

    ! save the data at the remote destination in boxRegExpansion(targetDisplacement)
    call MPI_Put(hostRegExpansion(sourceDisplacement), targetCount, MPI_DOUBLE_COMPLEX, rmRank, &
             targetDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)

    ! close the one-sided window
    call MPI_Win_fence(0, window, mpiError)
    call MPI_Win_free(window, mpiError)
#endif

#ifdef TIMER
#ifndef MAX_PERF
#ifdef USE_REDUCE
    timeConsumed = TimerCheck(3)
    timeConsumedMin = timeConsumed
    timeConsumedMax = timeConsumed
    call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
    call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    timeConsumed = TimerCheck(3)
    timeConsumedMin = timeConsumed
    timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MAX_PERF
    if(iam == prtRank) then
      write(*,'(a,2(f6.1,a))') "FMM Step 2: Timer: M2M E-boxes      min: ", &
                               TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    endif
    call TimerInit(3)
#endif ! TIMER

  endif ! ffLevel /= rmLevel
#endif ! end of computation of E-boxes

  !iBox = treeLevel(loLevel-1)%bstart
  !call print_complex(boxRegExpansion(box(iBox)%expStart), 5, "step 2")
!  i=0
!  do L = 0, 5
!    i = i + L
!    write(*,'(10f15.6)') (dble(box(0)%chargeExpansion(i+m)),m=0,L)
!    write(*,'(3x,10f15.6)') (aimag(box(0)%chargeExpansion(i+m)),m=0,L)
!  enddo
!  error stop

  ! Retrieve remote boxes needed for FF computation
  ! TODO: Will it shorten the communication time to invoke MPI_Get from each of 8 OpenMP cores per node?
#ifdef MPI
  if(nRanks > 1) then
    ! nFence = 500 works for nRepLayers=127 with  6-node sub-communicator in 18-node run
    ! nFence = 500 works for nRepLayers=127 with  7-node sub-communicator in 21-node run
    ! nFence = 500 works for nRepLayers=127 with  8-node sub-communicator in 24-node run
    ! nFence = 500 works for nRepLayers=127 with  9-node sub-communicator in 18-node run
    ! nFence = 500 works for nRepLayers=127 with 10-node sub-communicator in 20-node run
    ! Put and Get use fence differently in this code
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
    nFence = 1000  ! TODO test with 1000000
#else
    nFence = 1  !5  !2000  !5  !10  !50
#endif
#if defined (USE_ONESIDED_FENCE) && !defined(ONESIDED_COLLECT)
    ! important: each sub-communicator has a unique set of messages due to the difference in S-list between MPI ranks

    ! determine the min/max number of Get calls accross the ranks inside the sub-communicator
    call MPI_AllReduce(nBoxSegments, nMinMessages,    1, MPI_INTEGER8, MPI_MIN, MPI_COMM_RMGROUP, mpiError)  ! sub-comm
    call MPI_AllReduce(nBoxSegments, nMaxMessages,    1, MPI_INTEGER8, MPI_MAX, MPI_COMM_RMGROUP, mpiError)  ! sub-comm
    call MPI_AllReduce(nBoxSegments, nTotMessagesRMG, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_RMGROUP, mpiError)  ! sub-comm

    ! print the size of messages (in terms of number of boxes) for the sub-communicator which holds prtRank
    if(iam == prtRank) write(*,'(6x,a,4(i0,2x))') &
      "Min,max,total number of one-sided Get messages and nFence = ", nMinMessages, nMaxMessages, nTotMessagesRMG, nFence
#endif

    ! determine the buffer size in bytes for the one-sided window
    ! each rank in the sub-communicator has a unique value for buffer_size
    call MPI_Type_size(MPI_DOUBLE_COMPLEX, dcmplx_size, mpiError)
    n1_8 = nSzTri
    n2_8 = nBoxesRList
    n3_8 = nRemoteFFBoxes
    n4_8 = dcmplx_size
    buffer_size = (n1_8+1_8)*(n2_8+n3_8) * n4_8
    if(iam == prtRank) write(*,'(6x,a,i0,a,i0)') "Buffer size = ", buffer_size, " bytes for datatype ", dcmplx_size

#ifdef ONESIDED_COLLECT
    ! create a local message buffer (imessage) for the subcommunicator
    allocate(imessage(5,nBoxSegments),        &
             nMessages(0:rmSize-1),           &
             messageDisplacement(0:rmSize-1), &
             stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)

    ! populate the local message (buffer) info in the rank
    ! TODO replace the source and target rank ids by their global counterparts
    do i = 1, nBoxSegments
      imessage(1,i) = rmRank                          ! sourceRank
      imessage(2,i) = rList(segmentSBox(i))%expStart  ! sourceDisplacement
      imessage(3,i) = (nSzTri+1)*segmentNBoxes(i)     ! targetCount
      imessage(4,i) = segmentRank(i)                  ! targetRank
      imessage(5,i) = rList(segmentSBox(i))%expRemote ! targetDisplacement
    enddo

    ! gather the number of messages requested by the ranks in the sub-communicator
    ! each rank knows which data it needs and where to Get the required data from (including which rank holds those data)
    ! (the rank computes that information in the FMM Plan)
    ! however, the rank does not know beforehand how many Put messages it will send

    ! determine the number of Get calls accross the ranks in the sub-communicator
    call MPI_Allgather(nBoxSegments, 1, MPI_INTEGER, nMessages, 1, MPI_INTEGER, MPI_COMM_RMGROUP, mpiError)
    nMinMessages =  minval(nMessages(0:rmSize-1), dim = 1)
    nMaxMessages =  maxval(nMessages(0:rmSize-1), dim = 1)

    ! convert from the number of messages to the number of data points
    nMessages(:) = nMessages(:) * 5  ! each message has 5 data points

    ! compute displacements
    messageDisplacement(0) = 0
    messageDisplacement(1) = nMessages(0)
    do i = 2, rmSize-1
      messageDisplacement(i) = messageDisplacement(i-1) + nMessages(i-1)
    enddo

    ! get the total number of messages in the sub-communicator
    nTotMessagesRMG  = (messageDisplacement(rmSize-1) + nMessages(rmSize-1)) / 5

    ! set the max size for collective message buffer that receives the message data from other sub-communicators
    nMaxMessagesWorld = nTotMessagesRMG

#ifdef ONESIDED_PRINT
    if(nRanks / rmSize > 1) then
      ! this section appies to the case of multiple RM sub-communicators
      ! the rank that prints to stdout needs to receive message data from other sub-communicators
      ! determine the max value of nTotMessagesRMG across all RM sub-communicators so the printing rank allocates enough memory

      ! create an XX cross-communicator MPI_COMM_XXGROUP which connects the RM sub-communicators MPI_COMM_RMGROUP
      ! the XX cross-communicator includes only one rank from each RM sub-communicator
      ! the cross-connection happens at the slice which goes through prtRank
      ! each rank knows the absolute value of prtRank
      ! xxPrint is a zero-based index of the RM sub-communicator which holds prtRank
      xxPrint = leadRank / rmSize  ! the xxPrint (==prtRank) prints to stdout

      ! slice the RM sub-communicators along the line connecting the shift of prtRank in each sub-communicator
      ! xxComm is .T. if the MPI rank is in XX sub-comm or is .F. otherwise
      xxComm = ((iam - (iam / rmSize) * rmSize) == (prtRank - (prtRank / rmSize) * rmSize))

      if(xxComm) then
        color = 1
      else
        color = MPI_UNDEFINED
      endif
      call MPI_Comm_split(MPI_COMM_WORLD, color, xxRank, MPI_COMM_XXGROUP, mpiError)

      if(xxComm) then
        call MPI_Comm_rank(MPI_COMM_XXGROUP, xxRank, mpiError)
        call MPI_Comm_size(MPI_COMM_XXGROUP, xxSize, mpiError)

        ! let the prtRank == xxPrint know the max size of the message buffer
        call MPI_Reduce(nTotMessagesRMG, nMaxMessagesWorld, 1, MPI_INTEGER8, MPI_MAX, xxPrint, MPI_COMM_XXGROUP, mpiError)
      else
        xxRank = 0
        xxSize = 0  ! this rank is not a member of the cross-communicator
      endif
    endif
#endif

#if !defined(USE_ONESIDED_PUT)
    ! print the size of messages (in terms of number of boxes) for the sub-communicator which holds prtRank
    ! the use of lower-case "get" in the printout indicates a debug mode (collecting message data)
    if(iam == prtRank) write(*,'(6x,a,4i8)') &
      "Min,max,total number of one-sided get messages and nFence = ", nMinMessages, nMaxMessages, nMaxMessagesWorld, nFence
#endif

#ifdef ONESIDED_COLLECT
    ! create a collective message buffer to hold the messages exchanged between all ranks in the sub-communicator
    ! the message buffer should have enough size to hold the largest data set sent from any of the sub-communicators
    allocate(message(5,nMaxMessagesWorld), &
             bufferSizes(0:rmSize-1),   &
             stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)

    ! every rank gathers all messages from the sub-communicator they belong to
    call MPI_Allgatherv(imessage, nMessages(rmRank), MPI_INTEGER8, &
                message, nMessages, messageDisplacement, MPI_INTEGER8, &
                MPI_COMM_RMGROUP, mpiError)

    ! gather buffer sizes
    call MPI_Allgather(buffer_size, 1, MPI_INTEGER8, bufferSizes, 1, MPI_INTEGER8, MPI_COMM_RMGROUP, mpiError)
#endif

#ifdef ONESIDED_PRINT
    ! all sub-communicators must have the same node-width
    if(iam == prtRank) then
      write(*,'(a,2i8)') "Number of sub-communicators and their size: ", xxSize, rmSize

      ! print the content of the message buffer owned by prtRank before it gets overwritten
      write(*,'(a,i8)') "Sub-communicator: ", xxPrint

      do i = 0, rmSize-1
        write(*,'(a,3i12)') "rank nMessages: ", i, nMessages(i)
      enddo

      do i = 0, rmSize-1
        write(*,'(a,3i12)') "rank buffer_size,bytes: ", i, bufferSizes(i), bufferSizes(i) / dcmplx_size
      enddo

      write(*,'(a,i0)') "The total number of messages: ", nTotMessagesRMG

      do i = 1, nTotMessagesRMG
        write(*,'(a,5i10)') "rank src_displ trg_count trg_rank trg_displ: ", message(1:5,i)
      enddo
    endif

    ! case of multiple sub-communicators
    if(nRanks / rmSize > 1) then
      if(iam == prtRank) then
        do xxRank = 0, xxSize-1   ! loop over sub-communicators if there are more than one
          ! the rank xxPrint receives the data from other ranks in the cross-communicator
          ! populate the message buffer with the data received from the printing rank of each sub-communicator
          if(xxRank == xxPrint) cycle  ! skip this sub-communicator because its data have already been printed above

          write(*,'(a,i8)') "Sub-communicator: ", xxRank

          call MPI_Recv(nMessages, rmSize, MPI_INTEGER, xxRank, xxRank, MPI_COMM_XXGROUP, MPI_STATUS_IGNORE, mpiError)
          do i = 0, rmSize-1
            write(*,'(a,3i12)') "rank nMessages: ", i, nMessages(i)
          enddo

          call MPI_Recv(bufferSizes, rmSize, MPI_INTEGER8, xxRank, xxRank, MPI_COMM_XXGROUP, MPI_STATUS_IGNORE, mpiError)
          do i = 0, rmSize-1
            write(*,'(a,3i12)') "rank buffer_size,bytes: ", i, bufferSizes(i), bufferSizes(i) / dcmplx_size
          enddo

          call MPI_Recv(nTotMessagesRMG, 1, MPI_INTEGER8, xxRank, xxRank, MPI_COMM_XXGROUP, MPI_STATUS_IGNORE, mpiError)
          write(*,'(a,i0)') "The total number of messages: ", nTotMessagesRMG

          call MPI_Recv(message, nTotMessagesRMG*5, MPI_INTEGER8, xxRank, xxRank, MPI_COMM_XXGROUP, MPI_STATUS_IGNORE, mpiError)
          do i = 1, nTotMessagesRMG
            !write(*,'(a,5i10)') "rank src_displ trg_count trg_rank trg_displ: ", message(1:5,i)
            write(*,'(a,5i10)') "src_rank src_displ trg_count trg_rank trg_displ: ", &
                    message(1,i), &
                    message(2,i) / (nSzTri+1), &
                    message(3,i) / (nSzTri+1), &
                    message(4,i), &
                    message(5,i) / (nSzTri+1)
          enddo
        enddo

      elseif(xxComm) then

        ! ranks in the cross-communicator send their data to the rank xxPrint (prtRank)
        call MPI_Send(nMessages, rmSize, MPI_INTEGER, xxPrint, xxRank, MPI_COMM_XXGROUP, mpiError)
        call MPI_Send(bufferSizes, rmSize, MPI_INTEGER8, xxPrint, xxRank, MPI_COMM_XXGROUP, mpiError)
        call MPI_Send(nTotMessagesRMG, 1, MPI_INTEGER8, xxPrint, xxRank, MPI_COMM_XXGROUP, mpiError)
        call MPI_Send(message, nTotMessagesRMG*5, MPI_INTEGER8, xxPrint, xxRank, MPI_COMM_XXGROUP, mpiError)
      endif

      ! free the XX cross-communicator
      call MPI_Comm_Free(MPI_COMM_XXGROUP, mpiError)
    endif

    if(iam == prtRank) then
      write(*,'(a)') "End of messages"
    endif

    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    call MPI_Comm_Free(MPI_COMM_NODE, mpiError)
    call MPI_Comm_Free(MPI_COMM_RMGROUP, mpiError)
    call MPI_Finalize(mpiError)
    if(iam == prtRank) write(*,'(a)') "Graceful Termination"
    stop
#endif

#if defined(USE_ONESIDED_PUT)
    ! determine the min/max number of Put calls accross the ranks in each sub-communicator
    nMessages = 0
    do i = 1, nTotMessagesRMG
      iRank = message(4,i)  ! local rank id which initiates Put in the sub-communicator
      nMessages(iRank) = nMessages(iRank) + 1
    enddo
    nMinMessages =  minval(nMessages(0:rmSize-1), dim = 1)
    nMaxMessages =  maxval(nMessages(0:rmSize-1), dim = 1)
    ! TODO we do not know how many Put messages a rank will execute
    ! print the data for sub-communicator which includes prtRank
    if(iam == prtRank) write(*,'(6x,a,4i8)') &
      "Min,max,total number of one-sided Put messages and nFence = ", nMinMessages, nMaxMessages, nTotMessagesRMG, nFence
#endif
#endif

    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    call TimerInit(3)

#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
    ! steer one-sided communications to the host memory
    call MPI_Win_Create(hostRegExpansion, buffer_size, dcmplx_size, MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
#else
    ! steer one-sided communications to the device memory
    call MPI_Win_Create(boxRegExpansion, buffer_size, dcmplx_size, MPI_INFO_NULL, MPI_COMM_RMGROUP, window, mpiError)
#endif
    call MPI_Win_fence(0, window, mpiError)

#if defined(USE_ONESIDED_PUT)
    ! presently, each consecutive Put may go to a different rank in the sub-communicator
    ! TODO is it worth to consolidate all Put messages to the same destination?
    do i = 1, nTotMessagesRMG
      ! check every message to see if it outgoes from rmRank
      if(message(4,i) == rmRank) then
        targetRank         = message(1,i)  ! send the data to
        sourceDisplacement = message(2,i)
        targetCount        = message(3,i)
       !sourceRank         = message(4,i)  ! take the data from;    here, it is always rmRank
        targetDisplacement = message(5,i)

#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
        call MPI_Put(hostRegExpansion(targetDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
                 sourceDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)
#else
        ! save the data at the remote destination in boxRegExpansion(targetDisplacement)
        call MPI_Put(boxRegExpansion(targetDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
                 sourceDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)
#endif
      endif

#ifdef USE_ONESIDED_FENCE
      if(mod(i,nFence) == 0) call MPI_Win_fence(0, window, mpiError)
    enddo
    do i_8 = nTotMessagesRMG/nFence + 1_8, nMaxMessages/nFence
      call MPI_Win_fence(0, window, mpiError)
#endif
    enddo

#else
    !USE_ONESIDED_GET
    do i_8 = 1_8, nBoxSegments
      ! sanity check: first and last box in the segment must belong to the same hierarchy level
      !if(floor(dlog(dble(segmentSBox(i)*7+1)) * invLog8) /= &
      !   floor(dlog(dble((segmentSBox(i)+segmentNBoxes(i)-1)*7+1)) * invLog8)) error stop
  
      !sourceDisplacement = box(segmentSBox(i))%expStart
      !targetDisplacement = box(segmentSBox(i))%expRemote
  
#ifdef ONESIDED_COLLECT
     !sourceRank         = imessage(1,i)  ! save the data to;    here, it is always rmRank
      sourceDisplacement = imessage(2,i_8)
      targetCount        = imessage(3,i_8)
      targetRank         = imessage(4,i_8)  ! read the data from
      targetDisplacement = imessage(5,i_8)
#else
      sourceDisplacement = rList(segmentSBox(i_8))%expStart
      targetCount        = (nSzTri+1)*segmentNBoxes(i_8)
      targetRank         = segmentRank(i_8)
      targetDisplacement = rList(segmentSBox(i_8))%expRemote
#endif

#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
      ! perform one-sided communications in the host memory, host-to-host RDMA
      call MPI_Get(hostRegExpansion(sourceDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
                 targetDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)
#else
      ! perform one-sided communications in the device memory, device-to-device RDMA
      ! save the data locally in boxRegExpansion(sourceDisplacement)
      call MPI_Get(boxRegExpansion(sourceDisplacement), targetCount, MPI_DOUBLE_COMPLEX, targetRank, &
                 targetDisplacement, targetCount, MPI_DOUBLE_COMPLEX, window, mpiError)
#endif

#ifdef USE_ONESIDED_FENCE
      if(mod(i_8,nFence) == 0) call MPI_Win_fence(0, window, mpiError)
    enddo
    do i_8 = nBoxSegments/nFence + 1, nMaxMessages/nFence
      call MPI_Win_fence(0, window, mpiError)
#endif
    enddo
#endif

    call MPI_Win_fence(0, window, mpiError)
    call MPI_Win_free(window, mpiError)

#ifdef ONESIDED_COLLECT
    ! free the message buffers
    deallocate(imessage,            &
               nMessages,           &
               messageDisplacement, &
               message,             &
               bufferSizes,         &
               stat=allocStat)
#endif

#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
    ! copy hostRegExpansion to the device
    !$omp target teams distribute parallel do device(iDevice) &
#ifndef CPLR_NVFTN
    !$omp has_device_addr(boxRegExpansion) &
#endif
    !$omp map(to: hostRegExpansion) &
    !$omp map(to: szHostRegExp)
    do i_8 = 0_8, szHostRegExp
      boxRegExpansion(i_8) = hostRegExpansion(i_8)
    enddo
    !TODO deallocate hostRegExpansion in here
#endif

  endif
#endif  ! ifdef MPI

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "      Retrieval of remote FF boxes: min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER


#ifdef USE_DEVPTR
  ! make the array boxIrrExpansion (allocated) in GPU
  !!$omp allocate allocator(omp_target_device_mem_alloc)
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):boxIrrExpansion)
#endif
  allocate(boxIrrExpansion(0:nBoxesSList*(nSzTri+1)), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)

  !
  ! FMM step 3      Perform M2L translation
  !
  ! Compute far-field interaction
  ! Compute potential expansion at each level by using MPI parallelization
  !
#ifdef OFFL_GPU
  ! zero-initialize the array on GPU
#ifdef CPLR_NVFTN
  !$omp target teams distribute parallel do device(iDevice) 
#else
  !$omp target teams distribute parallel do has_device_addr(boxIrrExpansion) device(iDevice) 
#endif
  do i = 0, nBoxesSList*(nSzTri+1)
    boxIrrExpansion(i) = dcmplx(0.0d0, 0.0d0)
  enddo
#else
  ! CPU-only computation
  boxIrrExpansion = dcmplx(0.0d0, 0.0d0)
#endif
  ! 
#ifdef DYNAMICLIST
  do iLevel = loLevel, ffLevel, -1
    nbEdge = 2 ** iLevel      ! number of cubes along the edge at the given level of hierarchy
    nbEdgeMaxID = nbEdge - 1  ! max box id along the edge at the lowest level of hierarchy
    nbParentEdgeMaxID = 2 ** (iLevel - 1) - 1  ! max box id along the edge on parent level
    !
    ! loop over boxes at iLevel
    do iBox = treeLevel(iLevel)%bStart, treeLevel(iLevel)%bEnd
      !
      ! assign job to MPI tasks
      !
      if(box(iBox)%proc == iam) then
        ix = box(iBox)%xyz(1)
        iy = box(iBox)%xyz(2)
        iz = box(iBox)%xyz(3)

        ! far-field expansion will be accumulated
        potentialExpansion = dcmplx(0.0d0, 0.0d0)  ! zero-initialization is necessary for serial code

        ! dynamically build the far-field box list
        pBox = box(iBox)%parent  ! parent box of iBox

        ! far-field limits
        fax = max0(box(pBox)%xyz(1) - ws, 0) * 2  
        fay = max0(box(pBox)%xyz(2) - ws, 0) * 2
        faz = max0(box(pBox)%xyz(3) - ws, 0) * 2

        fbx = min0(box(pBox)%xyz(1) + ws, nbParentEdgeMaxID) * 2 + 1
        fby = min0(box(pBox)%xyz(2) + ws, nbParentEdgeMaxID) * 2 + 1
        fbz = min0(box(pBox)%xyz(3) + ws, nbParentEdgeMaxID) * 2 + 1

        ! near-field limits
        nax = max0(ix-ws,0)
        nay = max0(iy-ws,0)
        naz = max0(iz-ws,0)

        nbx = min0(ix+ws,nbEdgeMaxID)
        nby = min0(iy+ws,nbEdgeMaxID)
        nbz = min0(iz+ws,nbEdgeMaxID)

        i = 0
        do jz = faz, fbz
          do jy = fay, fby
            do jx = fax, fbx
              if(jx.ge.nax .and. jx.le.nbx .and. &
                 jy.ge.nay .and. jy.le.nby .and. &
                 jz.ge.naz .and. jz.le.nbz ) cycle

              ! find jBox based on grid coordinates jx,jy,jz for hierarchy level iLevel
              jBox = treeLevel(iLevel)%bstart   ! starting number for the boxes at this level
              tjx = jx
              tjy = jy
              tjz = jz
              nbParent = nbEdge  ! number of boxes along the edge at the parent level
              do aLevel = 1, iLevel  ! although this loop looks going down the hierarchy, it actually goes upward
                tjx = mod(tjx,nbParent)
                tjy = mod(tjy,nbParent)
                tjz = mod(tjz,nbParent)
                nbParent = nbParent / 2
                pjx = tjx / nbParent
                pjy = tjy / nbParent
                pjz = tjz / nbParent
                jBox = jBox + (pjz * 4 + pjy * 2 + pjx) * nbParent * nbParent * nbParent
              enddo

              ! the above block may be replaced by call to function
              !jBox = idBox(jx, jy, jz, iLevel)

              i = i + 1
              boxList(i) = jBox   ! far-field list
            enddo
          enddo
        enddo
        if(i /= box(iBox)%nFar) then
          write(*,'(a,2i5)') "Error in the number of FF boxes: ", i, box(iBox)%nFar
          write(*,'(a,3i5,a,3i5)') "Box: ", ix,iy,iz, "  error    box: ", iBox, pBox, iLevel
          error stop
        endif

        ! loop over far-field boxes stored in boxList
        !$omp parallel do reduction (+:potentialExpansion) &
        !$omp private(aBox, jBox,jx,jy,jz,x,y,z,iThread,shiftedExpansion) default(shared) 
        do i = 1, box(iBox)%nFar
          jBox = boxList(i)

          jx = box(jBox)%xyz(1)
          jy = box(jBox)%xyz(2)
          jz = box(jBox)%xyz(3)
          !
          ! find helper vector that moves the center of jBox to the center of iBox (from - to)
          iOrigin(1) = dble(jx-ix) * treeLevel(iLevel)%edge
          iOrigin(2) = dble(jy-iy) * treeLevel(iLevel)%edge
          iOrigin(3) = dble(jz-iz) * treeLevel(iLevel)%edge
          !
          ! compute shifted expansion M2L
          irrExpansion = dcmplx(0.0d0, 0.0d0)
#ifdef M2LSTD
          ! use standard O(N^4) procedure
          !call IrregularHarmonicsLocal(iOrigin(1), iOrigin(2), iOrigin(3), shiftedExpansion)
          call IrregularExpansion(iOrigin, 1.0d0, 1, zOrigin, tmpExpansion, shiftedExpansion)
          !call PotentialExpShiftInverted(shiftedExpansion, box(jBox)%chargeExpansion, irrExpansion)
          call PotentialExpShiftInverted(shiftedExpansion, boxRegExpansion(box(jBox)%expStart), irrExpansion)
#elif M2LROT
          ! use rotation-based translation
          !call M2Lrotate2(-iOrigin, expalp, box(jBox)%chargeExpansion, tmpExpansion, irrExpansion)
          call M2Lrotate2(-iOrigin, expalp, boxRegExpansion(box(jBox)%expStart), tmpExpansion, irrExpansion)
          !call CompareExpansions(irrExpansion, shiftedExpansion, 1.0d-8)
#else
          write(*,'(/a/a)') "Error: either -DM2LSTD or -DM2LROT have to be specified at compile time", &
                            "       to choose between the standard and rotation-based M2L routines" 
          error stop
#endif
          potentialExpansion = potentialExpansion + irrExpansion

          !write(*,'(/a,3i6,4x,3f8.3,4x,2i6)') &
          !   "vam: iLevel,iBox,jBox,txyz,njboxes,i= ", iLevel, iBox, jBox, x,y,z, box(iBox)%nFar, i
          !write(*,'(7x,a,3i4,4x,3i4)') "ix,iy,iz, jx,jy,jz = ", ix,iy,iz, jx,jy,jz
          !call CompareExpansions(oShift, irrExpansion, 1.0d-5)
        enddo

        !box(iBox)%potentialExpansion = box(iBox)%potentialExpansion + potentialExpansion
        boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) = &
        boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) + potentialExpansion(:)
      endif
      !
    enddo
  enddo
#else
  ! .not. DYNAMICLIST
  !do iLevel = loLevel, ffLevel, -1
  !do iLevel = ffLevel, loLevel
    ! loop over boxes at iLevel
#ifdef OFFL_CPU
    !$omp parallel do &
#elif  defined(OFFL_GPU) && defined(CPLR_NVFTN)
    !$omp target teams distribute parallel do device(iDevice) thread_limit(16) &
#elif  defined(OFFL_GPU) && !defined(CPLR_NVFTN)
    !$omp target teams distribute parallel do device(iDevice) thread_limit(16) &
    !$omp has_device_addr(boxRegExpansion) &
#ifdef USE_DEVPTR
    !$omp has_device_addr(boxIrrExpansion) &
#else
    !$omp map(to: boxIrrExpansion) &
#endif

#ifdef DEBUG_STEP3
    !$omp map(from: outexp) &
#endif

#endif
    !$omp private(iLevel, ix,iy,iz, shiftedExpansion) &
    !$omp private(aBox, jBox, jx,jy,jz, x,y,z, ir,ir2,ir3,xpy, oShift) &
    !$omp private(L, z2lm1, aL, aL1, aL2, LL1, m) &
    !$omp private(i0, rsum, j, aJ, aJL, mS, k, mSm, csum)
    !do iBox = treeLevel(iLevel)%bStart, treeLevel(iLevel)%bEnd
    !do iBox = treeLevel(ffLevel)%bStart, nBoxesTotal-1  ! skip the boxes which have no FF neighbors
    !do iBox = treeLevel(ffLevel)%bStart * merge(1, 0, (iam==0)) + iam * nBoxesSList, &
    !          min0(iam * nBoxesSList + nBoxesSList-1, nBoxesTotal-1)
    !do iBox = iam * nBoxesSList, min0(iam * nBoxesSList + nBoxesSList-1, nBoxesTotal-1)
    ! skip the upper-level boxes which have no FF neighbors
    !do cBox = treeLevel(ffLevel)%rStart, nBoxesSList-1  
    do cBox = 1, nBoxesSList-1  ! skip the master box 0; the boxes in S-list start from index 1 at ffLevel
      !if(box(iBox)%proc == iam) then
      !
      ! each iBox is to be processed by the assigned MPI rank, box(iBox)%proc
      ! iBox interacts with jBox 
      ! jBox has to have data on the same rank as that of iBox
      !
      !if(box(iBox)%proc /= iam) then
      !  write(*,'(a,3(2x,i0))') "Error: Process id mismatch", box(iBox)%proc, iam, iBox
      !  error stop
      !endif
      iLevel = sList(cBox)%iLevel
      ix = sList(cBox)%xyz(1)
      iy = sList(cBox)%xyz(2)
      iz = sList(cBox)%xyz(3)
      !iLevel = box(iBox)%iLevel
      !ix = box(iBox)%xyz(1)
      !iy = box(iBox)%xyz(2)
      !iz = box(iBox)%xyz(3)

      ! loop over far neighboors
      !potentialExpansion = 0.0d0  ! zero-initialization is necessary for serial code
#if defined(CPLR_NVFTN) && defined(OFFL_GPU)
      !shiftedExpansion = dcmplx(0.0d0, 0.0d0)  ! disable for nvfortran when offloading to GPU
#else
      shiftedExpansion = dcmplx(0.0d0, 0.0d0) 
#endif

      do aBox = sList(cBox)%farStart, sList(cBox)%farEnd
      !do aBox = box(iBox)%farStart, box(iBox)%farEnd
      !do aBox = box(iBox)%farStart, box(iBox)%farEnd
      !do aBox = box(iBox)%farStart*(nSzTri+1), (box(iBox)%farStart + box(iBox)%nFar - 1)*(nSzTri+1)+nSzTri

        ! find translation vector that moves the center of jBox to the center of iBox (from - to)
        jBox = farNeighbors(aBox)
        !jBox = farNeighbors(aBox/(nSzTri+1))
        !jx = box(jBox)%xyz(1) ! the far-neighbor boxes are outside of S-list, so have to use box not sList
        !jy = box(jBox)%xyz(2)
        !jz = box(jBox)%xyz(3)
        jx = rList(jBox)%xyz(1) ! the far-neighbor boxes are outside of S-list, so have to use box not sList
        jy = rList(jBox)%xyz(2)
        jz = rList(jBox)%xyz(3)
        !
#ifdef M2LSTD
        ! use the standard O(N^4) procedure

        x = dble(jx-ix) * treeLevel(iLevel)%edge
        y = dble(jy-iy) * treeLevel(iLevel)%edge
        z = dble(jz-iz) * treeLevel(iLevel)%edge
        ir  = 1.0d0 / sqrt( dmax1(x*x + y*y + z*z, 1.0d-10) )
        ir2 = ir * ir
        ir3 = ir * ir2
        xpy = dcmplx(-x * ir2, -y * ir2)

        ! initialize first few irregular harmonics
        ! (0,0)
        oShift(0) = dcmplx(ir, 0.0d0)

        ! (1,0)
        oShift(1) = dcmplx(z * ir3, 0.0d0)
        ! (1,1)
        oShift(2) = dcmplx(-x * ir3, -y * ir3)

        ! obtain higher-order harmonics from recurrence relations
        do L = 2, nL
          z2lm1 = z * dble(2*L-1)
          aL  = tri(L)
          aL1 = tri(L-1)
          aL2 = tri(L-2)
          LL1 = (L-1)*(L-1)

          ! (L,0) ... (L,L-2)   use first recurrence relation
          do m = 0, L-2
            oShift(aL+m) = (z2lm1 * oShift(aL1+m) - dble(LL1-m*m) * oShift(aL2+m)) * ir2
          enddo

          ! (L,L-1)             use second recurrence relation
          oShift(aL+L-1) = z2lm1 * oShift(aL1+L-1) * ir2

          ! (L,L)               use third recurrence relation
          oShift(aL+L) = dble(2*L-1) * xpy * oShift(aL1+L-1)
        enddo

#ifdef NEWSTEP3
        i0 = box(jBox)%expStart
        nLx = nL + 2
        do Lm = 0, nSzTri
        !Lm = mod(aBox,nSzTri+1)
          L = Lm / nLx
          m = mod(Lm, nLx)
          if(m > L) then
            L = nLx - L - 2
            m = nLx - m - 1
          endif

!        do L = 0, nL
          aL  = tri(L)
!          do m = 0, L

            csum = dcmplx(0.0d0, 0.0d0)

!           nLL  = nL-L
!           nLL1 = nL-L+1
!           do jkx = 0, nLL1*nLL1-1
!             jL = jkx / nLL1
!             km = mod(jkx, nLL1)
!             if(km > jL) then
!               jL = nLL - jL
!               km = km - nLL1
!             endif
!             j = jL + L
!             k = km + m
            do j = L, nL
              aJ  = tri(j)
              aJL = tri(j-L)
              mS  = mSign(mod(j-L,2))
              mSm = mSign(mod(j-L+m,2))

              do k = -j+L+m, j-L+m
                csum = csum + mSm * dconjg(boxRegExpansion(i0+aJL+iabs(k-m)) * oShift(aJ+iabs(k))) &
                                  * merge(1.0d0, 0.0d0, (k < 0)) &
                + mSign(mod(j-L+iabs(k-m),2)) * dconjg(boxRegExpansion(i0+aJL+iabs(k-m))) * oShift(aJ+iabs(k)) &
                                  * merge(1.0d0, 0.0d0, (k >= 0 .and. k < m)) &
                + mS * boxRegExpansion(i0+aJL+iabs(k-m)) * oShift(aJ+iabs(k)) &
                                  * merge(1.0d0, 0.0d0, (k >= m)) 
              enddo  ! index k

            enddo  ! index j or jkx

            shiftedExpansion(aL+m) = shiftedExpansion(aL+m) + csum
!          enddo  ! index m
!            shiftedExpansion(Lm) = shiftedExpansion(Lm) + csum
        enddo  ! index Lm
#else
        !i0 = box(jBox)%expStart  ! jBox is a far-neighbor box; must use structure box not sList
        i0 = rList(jBox)%expStart  ! jBox is a far-neighbor box; must use structure box not sList
        do L = 0, nL
            aL  = tri(L)

            ! m=0
            rsum = 0.0d0
            do j = L, nL
              aJ  = tri(j)
              aJL = i0 + tri(j-L)
              mS = mSign(mod(j-L,2))
              ! k = 0
                rsum = rsum + mS * dble( boxRegExpansion(aJL) * oShift(aJ) ) * 0.5d0

              ! k > 0;  m+k <= j-L
              do k = 1, j-L
                rsum = rsum + mS * dble( boxRegExpansion(aJL+k) * oShift(aJ+k) )
              enddo
            enddo
            shiftedExpansion(aL) = shiftedExpansion(aL) + rsum * 2.0d0

          do m = 1, L
            mSm = mSign(mod(m,2))

            csum = dcmplx(0.0d0)
            do j = L, min0(nL,L+m-1)
              aJ  = tri(j)
              aJL = i0 + tri(j-L)
              mS = mSign(mod(j-L,2))

              ! k < 0

              ! k = 0

              ! k > 0; k < m;  m-k <= j-L
              do k = max0(1,m-j+L), m-1
                csum = csum + dconjg(mS * boxRegExpansion(aJL+m-k)) * oShift(aJ+k) * mSign(mod(m-k,2))
              enddo

              ! k = m;  m <= j
                csum = csum + mS * boxRegExpansion(aJL) * oShift(aJ+m)

              ! k > m;  k-m <= j-L
              do k = m+1, j-L+m
                csum = csum + mS * boxRegExpansion(aJL+k-m) * oShift(aJ+k)
              enddo
            enddo

            do j = L+m, nL
              aJ  = tri(j)
              aJL = i0 + tri(j-L)
              mS = mSign(mod(j-L,2))

              ! k < 0;  m+k <= j-L
              do k = 1, j-L-m
                csum = csum + dconjg(mS * boxRegExpansion(aJL+m+k)) * dconjg(oShift(aJ+k)) * mSm
              enddo

              ! k = 0; j >= L+m
                csum = csum + dconjg(mS * boxRegExpansion(aJL+m)) * oShift(aJ) * mSm

              ! k >= 0; k < m;  m-k <= j-L
              do k = max0(1,m-j+L), m-1
                csum = csum + dconjg(mS * boxRegExpansion(aJL+m-k)) * oShift(aJ+k) * mSign(mod(m-k,2))
              enddo

              ! k = m;  m <= j
                csum = csum + mS * boxRegExpansion(aJL) * oShift(aJ+m)

              ! k > m;  k-m <= j-L
              do k = m+1, j-L+m
                csum = csum + mS * boxRegExpansion(aJL+k-m) * oShift(aJ+k)
              enddo
            enddo

            shiftedExpansion(aL+m) = shiftedExpansion(aL+m) + csum
          enddo
        enddo
#endif
#elif M2LROT
        ! use rotation-based translation
        !call M2Lrotate3(-iOrigin, expalp, box(jBox)%chargeExpansion, tmpExpansion, scrExpansion, shiftedExpansion)
        call M2Lrotate3(-iOrigin, expalp, boxRegExpansion(box(jBox)%expStart), tmpExpansion, scrExpansion, shiftedExpansion)
        !call CompareExpansions(shiftedExpansion, oShift, 1.0d-8)
#else
        write(*,'(/a/a)') "Error: either -DM2LSTD or -DM2LROT have to be specified at compile time", &
                          "       to choose between the standard and rotation-based M2L routines"
        error stop
#endif
      enddo  ! loop over aBox

      !!box(iBox)%potentialExpansion = box(iBox)%potentialExpansion + shiftedExpansion
      !boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) = &
      !boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) + shiftedExpansion(:)
      !boxIrrExpansion(sList(cBox)%expStart:sList(cBox)%expStart+nSzTri) = &
      !boxIrrExpansion(sList(cBox)%expStart:sList(cBox)%expStart+nSzTri) + shiftedExpansion(:)
      ! TODO: can accumulate directly into boxIrrExpansion
      boxIrrExpansion(sList(cBox)%expStart:sList(cBox)%expStart+nSzTri) = shiftedExpansion(:)

#ifdef DEBUG_STEP3
      ! save first nSzPrt1 elements of multipole expansion of boxes at ffLevel
      if(iLevel == ffLevel) outexp((cBox-1)*nSzPrt1 : (cBox-1)*nSzPrt1 + nSzPrt) = shiftedExpansion(0:nSzPrt)
#endif
      
      !endif  ! box(iBox)%proc == iam
    enddo  ! loop over cBox (iBox)
  !enddo  ! loop over iLevel
#endif

#ifdef DEBUG_STEP3
  ! print first few elements of the irregular multipole expansion of the boxes at ffLevel
  rBoxBottomStart = 1
  do cBox = 1, nBoxesSList-1  ! skip the master box 0; the boxes in S-list start from index 1 at ffLevel
    if(sList(cBox)%iLevel == ffLevel) rBoxBottomEnd = cBox
  enddo
  call PrintExpansion(rBoxBottomStart, rBoxBottomEnd, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1, nbr, outexp, prtexp)
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 3: Timer: M2L              min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    write(*,'(a,f6.1,a)')    "                   M2L d-mat comp      : ", timeM2Ldmat, "s"
    write(*,'(a,f6.1,a)')    "                   M2L z-translat      : ", timeM2Ldztr, "s"
    write(*,'(a,f6.1,a)')    "                   M2L d-mat rot       : ", timeM2Ldrot, "s"
  endif
  call TimerInit(3)
#endif ! TIMER


#if defined(OFFL_GPU) && defined(USE_DEVPTR)
  !$omp target exit data map(delete: boxRegExpansion, farNeighbors)
#else
  deallocate(boxRegExpansion, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
#endif

  !cBox = treeLevel(ffLevel)%rStart
  !call print_complex(boxIrrExpansion(sList(cBox)%expStart), 5, "step3")
  !
  ! FMM step 4      Perform L2L translation
  !
  ! Pass the irregular multipole expansion from parents to children
  ! This step is skipped when the Number of hierarchy levels, loLevel < 3
  !
  ! Phase 1: pass the irr expansion from the floating boxes down to ddLevel
  !            this phase involves MPI comminication between MPI ranks
  ! Phase 2: pass the irr expansion from ddLevel boxes down to the lowest level boxes
  !            no MPI communication between MPI ranks happens at this phase
  !
  nMembers = loLevel - ffLevel  ! go from ffLevel to loLevel
#ifdef OFFL_CPU
  !$omp parallel do     &
#elif  OFFL_GPU
  !!$omp target teams distribute parallel do device(iDevice) thread_limit(2) &
  !$omp target teams distribute parallel do device(iDevice) &
#ifndef CPLR_NVFTN
  !$omp has_device_addr(boxIrrExpansion) &
#endif

#ifdef DEBUG_STEP4
    !$omp map(from: outexp) &
#endif

#endif
  !$omp private(iOrigin,pBox,iMember,x,y,z,r2,oShift,L,z2lm1,aL,aL1,aL2,LL,m) &
  !$omp private(i0,rsum,j,aJ,aJL,k,shiftedExpansion,mSm,csum)
  !do iBox = treeLevel(ffLevel)%bStart * merge(1, 0, (iam==0)) + iam * nBoxesSList, &
  !          min0(iam * nBoxesSList + nBoxesSList-1, nBoxesTotal-1)
  !do iBox = iam * nBoxesSList, min0(iam * nBoxesSList + nBoxesSList-1, nBoxesTotal-1)
  do iBox = sBoxLLstart, nBoxesSList-1 ! Loop over boxes of the lowest level
    !
    ! find center of the box
    iOrigin(1) = (dble(sList(iBox)%xyz(1)) + 0.5d0) * cubeEdge
    iOrigin(2) = (dble(sList(iBox)%xyz(2)) + 0.5d0) * cubeEdge
    iOrigin(3) = (dble(sList(iBox)%xyz(3)) + 0.5d0) * cubeEdge

    !! build a genealogy tree for iBox
    !family(1) = box(iBox)%parent
    !do iMember = 2, nMembers
    !  family(iMember) = box( family(iMember-1) )%parent
    !enddo

    ! traverse the hierarchy levels from the bottom to the top and retrieve the parental lineage
    pBox = iBox
    do iMember = 1, nMembers
      !pBox = family(iMember)
      pbox = sList(pBox)%parent   ! without using family(:)

      ! determine translation vector: center of the child box <- center of the parent box
      x = iOrigin(1) - (dble( sList(pBox)%xyz(1) ) + 0.5d0) * treeLevel(loLevel - iMember)%edge
      y = iOrigin(2) - (dble( sList(pBox)%xyz(2) ) + 0.5d0) * treeLevel(loLevel - iMember)%edge
      z = iOrigin(3) - (dble( sList(pBox)%xyz(3) ) + 0.5d0) * treeLevel(loLevel - iMember)%edge

#ifdef L2LSTD
      ! use standard O(N^4) procedure
      ! use translation vector that moves the center of pBox to the center of iBox
      ! compute regular expansion, oShift(:)
      !call RegularExpansion(pOrigin, 1.0d0, 1, zOrigin, tmpExpansion, oShift)

      r2 = x*x + y*y + z*z

      ! initialize first few regular harmonics
      ! (0,0) = (L,m)
      oShift(0) = dcmplx(1.0d0, 0.0d0)

      ! (1,0)
      oShift(1) = dcmplx(z, 0.0d0)
      ! (1,1)
      oShift(2) = dcmplx(-0.5d0*x, 0.5d0*y)

      ! obtain higher-degree harmonics from recurrence relations
      do L = 2, nL
        z2lm1 = z * dble(2*L-1)
        aL  = tri(L)
        aL1 = tri(L-1)
        aL2 = tri(L-2)
        LL  = L*L

        ! (L,0) ... (L,L-2)   use first recurrece relation
        do m = 0, L-2
          oShift(aL+m) = (z2lm1 * oShift(aL1+m) - r2 * oShift(aL2+m)) / dble(LL-m*m)
        enddo

        ! (L,L-1)             use second recurrence relation
        oShift(aL+L-1) = z * oShift(aL1+L-1)

        ! (L,L)               use third recurrence relation
        oShift(aL+L) = oShift(aL1+L-1) * oShift(2) / dble(L)
      enddo

      !!call PotentialExpShift(box(pBox)%potentialExpansion, oShift, shiftedExpansion)
      !call PotentialExpShift(boxIrrExpansion(box(pBox)%expStart), oShift, shiftedExpansion)

      i0 = sList(pBox)%expStart
      do L = 0, nL
        aL  = tri(L)

        ! m=0; then both R and S share the same azimuthal number k
        !      that means we can use symmetry in the product of R and S
        !      and that their product is a real number
        ! oShift(j-L,k) * boxIrrExpansion((j,k)
          rsum = 0.0d0
          do j = L, nL
            aJ  = tri(j)
            aJL = tri(j-L)
            ! k = 0
              rsum = rsum + dble( oShift(aJL) * boxIrrExpansion(i0+aJ) ) * 0.5d0

            ! k > 0;  k-m <= j-L otherwise R=0; max(k) = nL-L = j-L
            do k = 1, j-L
              rsum = rsum + dble( oShift(aJL+k) * boxIrrExpansion(i0+aJ+k) )
            enddo
          enddo
          shiftedExpansion(aL) = rsum * 2.0d0

        do m = 1, L
          mSm = mSign(mod(m,2))
          csum = dcmplx(0.0d0)

          ! -j+L <= k-m <= j-L
          ! -j+L+m <= k <= j-L+m

          ! on the interval j = L...L+m-1 we have k > 0

          do j = L, min0(nL,L+m-1)
            aJ  = tri(j)
            aJL = tri(j-L)

            ! k < 0

            ! k = 0

            ! k > 0; k < m;  for R to be non-zero, m-k <= j-L, that is k >= m-j+L
            ! since max(j) = L+m-1, we get k >= 1, which is the initial condition
            do k = m-j+L, m-1
              csum = csum + dconjg(oShift(aJL+m-k)) * boxIrrExpansion(i0+aJ+k) * mSign(mod(m-k,2))
            enddo

            ! k = m;  m <= j
              csum = csum + oShift(aJL) * boxIrrExpansion(i0+aJ+m)

            ! k > m;  k-m <= j-L
            do k = m+1, j-L+m
              csum = csum + oShift(aJL+k-m) * boxIrrExpansion(i0+aJ+k)
            enddo
          enddo

          do j = L+m, nL
            aJ  = tri(j)
            aJL = tri(j-L)

            ! k < 0;  m+k <= j-L
            do k = 1, j-L-m
              csum = csum + dconjg(oShift(aJL+m+k) * boxIrrExpansion(i0+aJ+k)) * mSm
            enddo

            ! k = 0; j >= L+m
              csum = csum + dconjg(oShift(aJL+m)) * boxIrrExpansion(i0+aJ) * mSm

            ! k > 0; k < m;  m-k <= j-L
            do k = max0(1,m-j+L), m-1
              csum = csum + dconjg(oShift(aJL+m-k)) * boxIrrExpansion(i0+aJ+k) * mSign(mod(m-k,2))
            enddo

            ! k = m;  m <= j
              csum = csum + oShift(aJL) * boxIrrExpansion(i0+aJ+m)

            ! k > m;  k-m <= j-L
            do k = m+1, j-L+m
              csum = csum + oShift(aJL+k-m) * boxIrrExpansion(i0+aJ+k)
            enddo
          enddo

          shiftedExpansion(aL+m) = csum
        enddo

      enddo
#elif L2LROT
      ! use rotation-based translation
      !call L2Lrotate2(iOrigin, expalp, box(pBox)%potentialExpansion, tmpExpansion, shiftedExpansion)
      call L2Lrotate2(pOrigin, expalp, boxIrrExpansion(box(pBox)%expStart), tmpExpansion, shiftedExpansion)
      !call CompareExpansions(irrExpansion, shiftedExpansion, 1.0d-8)
#else
      write(*,'(/a/a)') "Error: either -DL2LSTD or -DL2LROT have to be specified at compile time", &
                        "       to choose between the standard and rotation-based L2L routines"
      error stop
#endif
      !box(iBox)%potentialExpansion = box(iBox)%potentialExpansion + shiftedExpansion
      !boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) = &
      !boxIrrExpansion(box(iBox)%expStart:box(iBox)%expStart+nSzTri) + shiftedExpansion(:)
      ! TODO: can accumulate directly into boxIrrExpansion
      boxIrrExpansion(sList(iBox)%expStart:sList(iBox)%expStart+nSzTri) = &
      boxIrrExpansion(sList(iBox)%expStart:sList(iBox)%expStart+nSzTri) + shiftedExpansion(:)

    enddo ! iMember chain for each iBox

#ifdef DEBUG_STEP4
      ! save first nSzPrt1 elements of multipole expansion of boxes at loLevel
      outexp((iBox-sBoxLLstart)*nSzPrt1 : (iBox-sBoxLLstart)*nSzPrt1 + nSzPrt) = &
                        boxIrrExpansion(sList(iBox)%expStart:sList(iBox)%expStart+nSzPrt)
#endif
  enddo ! iBox

#ifdef DEBUG_STEP4
  ! print first few elements of the irregular multipole expansion of the boxes at ffLevel
  call PrintExpansion(sBoxLLstart, nBoxesSList-1, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1, nbr, outexp, prtexp)
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 4: Timer: L2L              min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    write(*,'(a,f6.1,a)')    "                   L2L d-mat comp      : ", timeL2Ldmat, "s"
    write(*,'(a,f6.1,a)')    "                   L2L z-translat      : ", timeL2Ldztr, "s"
    write(*,'(a,f6.1,a)')    "                   L2L d-mat rot       : ", timeL2Ldrot, "s"
  endif
  call TimerInit(3)
#endif ! TIMER

  !iBox = treeLevel(loLevel)%rstart
  !call print_complex(boxIrrExpansion(sList(iBox)%expStart), 5, "step 4")
!  !
!  ! Compute electrostatic energy due to multipole expansions
!  ! This section is redundant; the same computation is done in the next section too
!  !
!  energy = 0.0d0
!  do iBox = treeLevel(loLevel)%bStart, treeLevel(loLevel)%bEnd
!    if(box(iBox)%proc == iam) then
!      sump = dble(box(iBox)%chargeExpansion(0) * box(iBox)%potentialExpansion(0)) * 0.5d0
!      i = 1
!      do L = 1, nL
!        sump = sump + dble(box(iBox)%chargeExpansion(i) * box(iBox)%potentialExpansion(i)) * 0.5d0
!        i = i + 1
!        do m = 1, L
!          sump = sump + dble(box(iBox)%chargeExpansion(i) * box(iBox)%potentialExpansion(i))
!          i = i + 1
!        enddo
!      enddo
!      energy = energy + sump
!    endif
!  enddo
!#ifdef MPI
!  call MPI_AllReduce(energy, sump, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiError)
!  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
!  energy = sump
!#endif
  !
  ! FMM step 5
  !
  ! Compute far-field potential and force in the position of particles 
  !   from irregular multipoles of FF accumulated in the boxes at the lowest hierarchy level
  !!$omp target data &
  !!$omp map(to:bC,bQ,tri,nearNeighbors,box,vertexMin,rij,brij,boxIrrExpansion) &
  !!$omp map(tofrom:aPotential,aForce)
  !call ElectricForceFF(loLevel, energy, aPotential, aForce, tmpExpansion, oShift)
  !
  ! Use irregular FF expansion accumulated in the box
  ! Works with non-periodic and periodic FMM
  ! All computed values are in Gauss units (charge in electrons, coordinates in Angstroms)
  ! Computation of energy includes correction against double counting
  !
  ! compute far-field potential and force in the position of particles
  ! obtain force by differentiation of the above computed FF potential (in multipole expansion)
  ! compute FF contribution to the energy of the system
  !

  ! Works with non-periodic and periodic FMM
  ! All computed values are in Gauss units (charge in electrons, coordinates in Angstroms)
  ! Computation of energy includes correction against double counting
  !
  ! compute far-field potential and force in the position of particles
  ! obtain force by differentiation of the above computed FF potential (in multipole expansion)
  ! compute FF contribution to the energy of the system
  !
  !!$omp target data map(to:bC,bQ,nearNeighbors,box,vertexMin) map(tofrom:aPotential,aForce)

!#ifdef OFFL_GPU
!  !$omp target enter data map(to: nearNeighbors) device(iDevice) &
!#ifdef FETCH_POTFORCE
!  !$omp map(tofrom: bPotential, bForce) 
!#else
!  !$omp map(to: bPotential, bForce) 
!#endif
!#endif

  !iAtomsBegin = dispAtoms(iam) + 1
  !iAtomsEnd   = dispAtoms(iam) + rankAtoms(iam)

  ! TODO think how to manage particles from NF boxes
  !iAtomsBegin = box(sList(sBoxLLstart)%idx)%atomStart
  !iAtomsEnd   = box(sList(nBoxesSList-1)%idx)%atomEnd

#ifdef COMP_FORCE
#ifdef USE_DEVPTR
  ! make the arrays bPotential and bForce resident (allocated) in GPU
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):bPotential,bForce)
#endif
  allocate(bPotential(rankAtoms(iam)), bForce(3,rankAtoms(iam)),  stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)

#ifndef USE_DEVPTR
  bPotential = 0.0d0
  bForce     = 0.0d0
#endif

#else
#ifdef USE_DEVPTR
  ! make the array bPotential resident (allocated) in GPU
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):bPotential)
#endif
  allocate(bPotential(rankAtoms(iam)),  stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)

#ifndef USE_DEVPTR
  bPotential = 0.0d0
#endif
#endif  ! ifdef COMP_FORCE

  !
  fmmEnergy = 0.0d0
#ifdef OFFL_CPU
  !$omp parallel do     &

#elif  OFFL_GPU
  !!$omp target teams distribute parallel do device(iDevice) thread_limit(2) &
  !$omp target teams distribute parallel do device(iDevice) &

#ifndef CPLR_NVFTN
  !$omp has_device_addr(boxIrrExpansion) &
  !$omp has_device_addr(bPotential) &
#ifdef COMP_FORCE
  !$omp has_device_addr(bForce) &
#endif
#endif

#if defined(FETCH_POTFORCE) && defined(COMP_FORCE)
  !$omp map(tofrom: bPotential, bForce) &
#endif
! !$omp map(to: bPotential, bForce) &

#ifdef DEBUG_PRINT_ENERGY
  !$omp map(tofrom: fmmEnergy) &
#endif

#endif

#ifdef DEBUG_PRINT_ENERGY
  !$omp reduction(+:fmmEnergy) &
#endif
#ifdef FSUM
  !$omp reduction(+:rankFSum) &
#endif
  !$omp private(ix,iy,iz,iBox,cBox,iOrigin,x,y,z,x2y2,xDx2y2,yDx2y2,r2,oShift,L,z2lm1,aL,aL1,aL2,LL,m) &
  !$omp private(sump,sumf,dL,i,aprod,bprod,dm,cprod)
  !do iAtom = sList(sBoxLLstart)%atomStart, sList(nBoxesSList-1)%atomEnd
  do iAtom = iAtomsBegin, iAtomsEnd

    ! find box indices for iAtom
    ix = int((bC(1,iAtom) - vertexMin(1)) * edgeReciprocal)
    iy = int((bC(2,iAtom) - vertexMin(2)) * edgeReciprocal)
    iz = int((bC(3,iAtom) - vertexMin(3)) * edgeReciprocal)

    ! box numbering rule: 
    !   each parent has 8 sequentially numbered child boxes
    !   the boxes at the child level are numbered in the sequence of their parent ordering
    ! for iLevel = 0,
    !   iBox = 0
    ! for iLevel = 1, 
    !   iBox = iz*4 + 
    !          iy*2 + 
    !          ix   + 1  (1 is for the master box at iLevel=0)
    ! for iLevel = 2, 
    !   each iz/2 parent-box carries 4 8-box cubes in 4x4x2-box-layer: 4*8=32 child boxes
    !   each mod(iz,2) single-box layer carries 4 boxes in xy-plane
    !   each iy/2 parent-box carries 2 8-box cubes in 4x2x2-box-slab: 2*8=16 child boxes
    !   each mod(iy,2) slab carries 2 boxes in x-axis
    !   each ix/2 parent-box encounters 8 child-boxes before itself
    !   each mod(ix,2) adds 0 or 1
    !   iz -> mod(iz,2)
    !   iBox = (iz/2)*32  + mod(iz,2)*4 + 
    !          (iy/2)*16  + mod(iy,2)*2 + 
    !          (ix/2)*8   + mod(ix,2)   + 9 (1+8 for boxes at iLevel=0,1)
    ! for iLevel = 3,
    !   iBox = (iz/4)*256 + (mod(iz,4)/2)*32 + mod(mod(iz,4),2)*4 +
    !          (iy/4)*128 + (mod(iy,4)/2)*16 + mod(mod(iy,4),2)*2 + 
    !          (ix/4)*64  + (mod(ix,4)/2)*8  + mod(mod(ix,4),2)   + 73 (1+8+64 for boxes at iLevel=0,1,2)
    ! for iLevel = 4,
    !   iBox = (iz/8)*2038 + (mod(iz,8)/4)*256 + (mod(mod(iz,8),4)/2)*32 + mod(mod(mod(iz,8),4),2)*4 +
    !          (iy/8)*1024 + (mod(iy,8)/4)*128 + (mod(mod(iy,8),4)/2)*16 + mod(mod(mod(iy,8),4),2)*2 + 
    !          (ix/8)*512  + (mod(ix,8)/4)*64  + (mod(mod(ix,8),4)/2)*8  + mod(mod(mod(ix,8),4),2)   + 585 (1+8+64i+512)
    ! supported values of nRepLayers = 0, 1, 3, 7, 15, 31, 63, and 127
    ! find which iBox the particle iAtom belongs to (works for cubic systems only)            loLevel

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
 !(8_8**loLevel - 1)/7_8  ! the total number of boxes for a system that has loLevel-1

!    nbParent = 2 ** loLevel  ! number of boxes along the edge at the parent level
!    iBox = treeLevel(loLevel)%bStart   ! starting number for the boxes at this level
!      ! count the number of boxes at the increasing granularity
!      do iLevel = 1, loLevel
!        ix = mod(ix,nbParent)
!        iy = mod(iy,nbParent)
!        iz = mod(iz,nbParent)
!        nbParent = nbParent / 2
!        jx = ix / nbParent
!        jy = iy / nbParent
!        jz = iz / nbParent
!        iBox = iBox + (jz * 4 + jy * 2 + jx) * nbParent * nbParent * nbParent
!      enddo

    ! get the cBox (S-list) from iBox (global B-list)
    !cBox = iBox - nBoxesBefore_loLevel - sBoxesTopDispl(iam) * 8_8**rLevel &
    !     + ((8_8**rLevel - 1_8)/7_8 * sBoxesTopNum(iam) + 1) + nFloatingBoxes
    cBox = iBox - cBoxShift

    ! find center of the box
    iOrigin(1) = vertexMin(1) + (dble(sList(cBox)%xyz(1)) + 0.5d0) * cubeEdge
    iOrigin(2) = vertexMin(2) + (dble(sList(cBox)%xyz(2)) + 0.5d0) * cubeEdge
    iOrigin(3) = vertexMin(3) + (dble(sList(cBox)%xyz(3)) + 0.5d0) * cubeEdge

    !! Not going to return the variable energy
    !!$omp parallel do reduction (+:energy)    &
    !!$omp parallel do     &
    !!$omp default(shared) &
    !!$omp private(iAtom,x,y,z,x2y2,xDx2y2,yDx2y2,sump,sumf,L,dL,aL1,i,aprod,bprod,m,dm,cprod,oShift)
    !do iAtom = box(iBox)%atomStart, box(iBox)%atomEnd

      x = bC(1,iAtom) - iOrigin(1)
      y = bC(2,iAtom) - iOrigin(2)
      z = bC(3,iAtom) - iOrigin(3)
      x2y2 = x*x + y*y

      ! modify the original coordinate by a very small increment to avoid division by zero
      if(x2y2 < zero) then
        bC(1,iAtom) = iOrigin(1) + dsign(zeroCoor,bC(1,iAtom))
        x = bC(1,iAtom) - iOrigin(1)
        x2y2 = x*x + y*y
      endif

      xDx2y2 = x / x2y2
      yDx2y2 = y / x2y2

      ! compute multipole expansion for a radius-vector iAtom
      !call RegularExpansion(bC(1,iAtom), 1.0d0, 1, iOrigin, tmpExpansion, oShift)

      ! shift particle coordinates to the new coordinate system
      x  = bC(1,iAtom) - iOrigin(1)
      y  = bC(2,iAtom) - iOrigin(2)
      z  = bC(3,iAtom) - iOrigin(3)
      r2 = x*x + y*y + z*z

      ! initialize first few regular harmonics
      ! (0,0) = (L,m)
      oShift(0) = dcmplx(1.0d0, 0.0d0)

      ! (1,0)
      oShift(1) = dcmplx(z, 0.0d0)
      ! (1,1)
      oShift(2) = dcmplx(-0.5d0*x, 0.5d0*y)

      ! obtain higher-degree harmonics from recurrence relations
      do L = 2, nL
        z2lm1 = z * dble(2*L-1)
        aL  = tri(L)
        aL1 = tri(L-1)
        aL2 = tri(L-2)
        LL  = L*L

        ! (L,0) ... (L,L-2)   use first recurrece relation
        do m = 0, L-2
          oShift(aL+m) = (z2lm1 * oShift(aL1+m) - r2 * oShift(aL2+m)) / dble(LL-m*m)
        enddo

        ! (L,L-1)             use second recurrence relation
        oShift(aL+L-1) = z * oShift(aL1+L-1)

        ! (L,L)               use third recurrence relation
        oShift(aL+L) = oShift(aL1+L-1) * oShift(2) / dble(L)
      enddo

      ! in the calculation of electrostatic potential, there is no need to use an 
      ! additional expansion for chargeless multipoles; one can use charged multipoles
      ! in the calculation and divide the final values by atom charge

      ! L = 0 and m = 0
      ! sump is sum for potential; sumf is sum for force
      !sump = dble(oShift(0)) * dble(box(iBox)%potentialExpansion(0)) * 0.5d0  ! correct for default double counting
      !sump = dble(oShift(0)) * dble(boxIrrExpansion(box(iBox)%expStart+0)) * 0.5d0  ! correct for default double counting
      sump = dble(oShift(0)) * dble(boxIrrExpansion(sList(cBox)%expStart+0)) * 0.5d0  ! correct for default double counting
      sumf = 0.0d0

      ! L > 0
      do L = 1, nL
        dL  = dble(L)
        aL1 = tri(L-1)
        i   = tri(L)

        ! m = 0
        ! potential; correct for default double counting
        !sump = sump + dble(oShift(i) * box(iBox)%potentialExpansion(i)) * 0.5d0
        !sump = sump + dble(oShift(i) * boxIrrExpansion(box(iBox)%expStart+i)) * 0.5d0
        sump = sump + dble(oShift(i) * boxIrrExpansion(sList(cBox)%expStart+i)) * 0.5d0

        ! force; divide by 2 so that index m can run only half-range (positive values)
        !aprod = dble( oShift(i)   * boxIrrExpansion(box(iBox)%expStart+i) ) * 0.5d0
        aprod = dble( oShift(i)   * boxIrrExpansion(sList(cBox)%expStart+i) ) * 0.5d0
        !bprod = dble( oShift(aL1) * boxIrrExpansion(box(iBox)%expStart+i) ) * 0.5d0
        bprod = dble( oShift(aL1) * boxIrrExpansion(sList(cBox)%expStart+i) ) * 0.5d0
        sumf(1) = sumf(1) + dL * xDx2y2 * aprod - z * xDx2y2 * bprod
        sumf(2) = sumf(2) + dL * yDx2y2 * aprod - z * yDx2y2 * bprod
        sumf(3) = sumf(3) + bprod

        ! do only half-range (only positive values); multiply by 2 to get the effect of the negative-range
        ! 0 < m < L
        do m = 1, L-1
          i = i + 1
          dm = dble(m)
          !bprod = dble( oShift(aL1+m) * boxIrrExpansion(box(iBox)%expStart+i) )
          bprod = dble( oShift(aL1+m) * boxIrrExpansion(sList(cBox)%expStart+i) )
          !cprod = oShift(i) * boxIrrExpansion(box(iBox)%expStart+i)
          cprod = oShift(i) * boxIrrExpansion(sList(cBox)%expStart+i)

          ! no need to divide the potential contribution by 2; doing only half-range of m does that
          !sump = sump + dble(oShift(i) * boxIrrExpansion(box(iBox)%expStart+i))
          sump = sump + dble(oShift(i) * boxIrrExpansion(sList(cBox)%expStart+i))
          ! due to m doing only half-range, force will have to be scaled by factor of 2 in the end
          sumf(1) = sumf(1) + dble( dcmplx(dL*xDx2y2,  dm*yDx2y2) * cprod ) - z * xDx2y2 * bprod
          sumf(2) = sumf(2) + dble( dcmplx(dL*yDx2y2, -dm*xDx2y2) * cprod ) - z * yDx2y2 * bprod
          sumf(3) = sumf(3) + bprod
        enddo

        ! m = L
        i = i + 1
        !cprod = oShift(i) * boxIrrExpansion(box(iBox)%expStart+i)
        cprod = oShift(i) * boxIrrExpansion(sList(cBox)%expStart+i)

        !sump = sump + dble(oShift(i) * boxIrrExpansion(box(iBox)%expStart+i))
        sump = sump + dble(oShift(i) * boxIrrExpansion(sList(cBox)%expStart+i))
        sumf(1) = sumf(1) + dble( dcmplx(dL*xDx2y2,  dL*yDx2y2) * cprod )
        sumf(2) = sumf(2) + dble( dcmplx(dL*yDx2y2, -dL*xDx2y2) * cprod )
        !sumf(3) contribution is zero for m = L
      enddo

      !energy = sump * bQ(iAtom)      ! sump is already corrected for double counting
      !aPotential(iAtom) = aPotential(iAtom) + sump * 2.0d0       ! contribution from all particles
      !aForce(1:3,iAtom) = aForce(1:3,iAtom) - sumf(1:3) * 2.0d0  ! compensate for half summation in index m
      bPotential(iAtom-iAtomsBegin+1) = bPotential(iAtom-iAtomsBegin+1) + sump * 2.0d0
#ifdef COMP_FORCE
      bForce(1:3,iAtom-iAtomsBegin+1) = bForce(1:3,iAtom-iAtomsBegin+1) - sumf(1:3) * 2.0d0
#endif
#ifdef FSUM
      if(sList(cBox)%group > 0) rankFSum(1:3) = rankFSum(1:3) - sumf(1:3) * 2.0d0 * bQ(iAtom)
#endif
    !enddo
#ifdef DEBUG_PRINT_ENERGY
    fmmEnergy = fmmEnergy + sump * bQ(iAtom)  ! not going to return energy
#endif
  enddo

#ifndef MAX_PERF
#ifdef DEBUG_PRINT_ENERGY
#ifdef MPI
call MPI_Reduce(fmmEnergy, esum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  fmmEnergy = esum
#endif
  energy = fmmEnergy * 0.5d0 * eConv
  if(iam == prtRank) write(*,'(a,23x,es26.18,a)') "mul =", energy, " Step 5"
#endif
#endif ! MAX_PERF

#if defined(OFFL_GPU) && defined(USE_DEVPTR)
  !$omp target exit data map(delete: boxIrrExpansion)
#else
  deallocate(boxIrrExpansion, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
#endif

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 5: Timer: far-field  pot:  min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER

  !call print_double(bPotential, bForce, 5, "step 5")
  !
  ! FMM step 6
  !
  ! Compute near-field potential and force in the position of particles
  ! Compute NF contribution to the energy of the system
  !
#ifdef DYNAMIC_NNLIST
  nbEdge = 2**loLevel       ! number of cubes along the edge at the lowest level of hierarchy
  nbEdgeMaxID = nbEdge - 1  ! max box id along the edge at the lowest level of hierarchy
  sump = 0.0d0
  do iBox = iam + treeLevel(loLevel)%bStart, treeLevel(loLevel)%bEnd, nRanks

    ! compute the near neighbor list
    ix = box(iBox)%xyz(1)
    iy = box(iBox)%xyz(2)
    iz = box(iBox)%xyz(3)

    do jz = max0(iz-ws,0), min0(iz+ws,nbEdgeMaxID)
      do jy = max0(iy-ws,0), min0(iy+ws,nbEdgeMaxID)
        do jx = max0(ix-ws,0), min0(ix+ws,nbEdgeMaxID)
          if(jx==ix .and. jy==iy .and. jz==iz) cycle     ! skip the interaction of iBox with itself

          ! find jBox based on grid coordinates jx,jy,jz for the lowest hierarchy level
          jBox = treeLevel(loLevel)%bstart   ! starting number for sequential enumeration of boxes at this level
          tjx = jx
          tjy = jy
          tjz = jz
          nbParent = nbEdge  ! number of boxes along the edge at the parent level
          do aLevel = 1, loLevel  ! although this loop looks going down the hierarchy, it actually goes upward
            tjx = mod(tjx,nbParent)
            tjy = mod(tjy,nbParent)
            tjz = mod(tjz,nbParent)
            nbParent = nbParent / 2
            pjx = tjx / nbParent
            pjy = tjy / nbParent
            pjz = tjz / nbParent
            jBox = jBox + (pjz * 4 + pjy * 2 + pjx) * nbParent * nbParent * nbParent
          enddo

          ! the above block may be replaced by call to function
          !jBox = idBox(jx, jy, jz, loLevel)

          !! perform sanity check to make sure that jBox has the expected grid coordinates jx,jy,jz
          !if(box(jBox)%xyz(1) /= jx .or. box(jBox)%xyz(2) /= jy .or. box(jBox)%xyz(3) /= jz) &
          !  write(*,'(a,2(3x,3i3))') "error jBox: x,y,z  jx,jy,jz  ", &
          !  box(jBox)%xyz(1), box(jBox)%xyz(2), box(jBox)%xyz(3), jx,jy,jz

          ! box-box interaction
          esum = 0.0d0
          !$omp parallel do reduction (+:esum) &
          !$omp private(potA, forceA, rij, r2, pot) default(shared) 
          do iAtom = box(iBox)%atomStart, box(iBox)%atomEnd
            potA = 0.0d0
            forceA = 0.0d0
            do jAtom = box(jBox)%atomStart, box(jBox)%atomEnd
              rij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
              r2       = rij(1)**2 + rij(2)**2 + rij(3)**2
              pot      = bQ(jAtom) / sqrt(r2)
              potA     = potA + pot
              forceA(1:3)= forceA(1:3) + pot / r2 * rij(1:3)
              !potA = potA + bQ(jAtom) / sqrt( (bC(1,iAtom)-bC(1,jAtom))**2 + &
              !                                (bC(2,iAtom)-bC(2,jAtom))**2 + &
              !                                (bC(3,iAtom)-bC(3,jAtom))**2 )
            enddo
            aPotential(iAtom) = aPotential(iAtom) + potA          ! potential on particle iAtom
            aForce(1:3,iAtom) = aForce(1:3,iAtom) + forceA(1:3)   ! force on particle iAtom
            esum = esum + potA * bQ(iAtom)
          enddo
          sump = sump + esum

        enddo
      enddo
    enddo

  enddo

#else
  ! .not. DYNAMIC_LIST
  fmmEnergy = 0.0d0
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU
  !$omp target teams distribute parallel do device(iDevice) &

#ifdef DEBUG_PRINT_ENERGY
  !$omp map(tofrom: fmmEnergy) &
#endif

#ifndef CPLR_NVFTN
  !$omp has_device_addr(bPotential) &
#ifdef COMP_FORCE
  !$omp has_device_addr(bForce) &
#endif
#endif

#endif

#ifdef DEBUG_PRINT_ENERGY
  !$omp reduction(+:fmmEnergy) &
#endif
#ifdef FSUM
  !$omp reduction(+:rankFSum) &
#endif
  !$omp private(ix,iy,iz, iBox, cBox, potA, forceA, ptrBox, jBox, jAtom, rij, r2, pot) default(shared) 
  !do iAtom = sList(sBoxLLstart)%atomStart, sList(nBoxesSList-1)%atomEnd
  do iAtom = iAtomsBegin, iAtomsEnd

    ! find box indices for iAtom
    ix = int((bC(1,iAtom) - vertexMin(1)) * edgeReciprocal)
    iy = int((bC(2,iAtom) - vertexMin(2)) * edgeReciprocal)
    iz = int((bC(3,iAtom) - vertexMin(3)) * edgeReciprocal)

    ! for nRepLayers = 0, iBox = 4 * iz + 2 * iy + ix + 9
    ! supported values of nRepLayers = 0, 1, 3, 7, 15, 31, 63, and 127
    ! find which iBox the particle iAtom belongs to (works for cubic systems only)            loLevel

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
 !(8_8**loLevel - 1)/7_8  ! the total number of boxes for a system that has loLevel-1

    ! get the cBox (S-list) from iBox (global B-list)
    !cBox = iBox - nBoxesBefore_loLevel - sBoxesTopDispl(iam) * 8_8**rLevel &
    !     + ((8_8**rLevel - 1_8)/7_8 * sBoxesTopNum(iam) + 1_8) + nFloatingBoxes
    cBox = iBox - cBoxShift

    !cBox = box(iBox)%idx  ! get the R-index from B-index
      ! iBox-jBox interaction
      potA = 0.0d0
      forceA = 0.0d0
      do ptrBox = sList(cBox)%nearStart, sList(cBox)%nearEnd  ! loop over neighbor boxes for cBox
        jBox = nearNeighbors(ptrBox)
        !do jAtom = box(jBox)%atomStart, box(jBox)%atomEnd
        do jAtom = rList(jBox)%atomStart, rList(jBox)%atomEnd
          rij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
          r2       = rij(1)**2 + rij(2)**2 + rij(3)**2
          pot      = bQ(jAtom) / sqrt(r2)
          potA     = potA + pot
          forceA(1:3)= forceA(1:3) + pot / r2 * rij(1:3)
          !potA = potA + bQ(jAtom) / sqrt( (bC(1,iAtom)-bC(1,jAtom))**2 + &
          !                                (bC(2,iAtom)-bC(2,jAtom))**2 + &
          !                                (bC(3,iAtom)-bC(3,jAtom))**2 )
        enddo
      enddo
      !aPotential(iAtom) = aPotential(iAtom) + potA         ! potential on particle iAtom
      !aForce(1:3,iAtom) = aForce(1:3,iAtom) + forceA(1:3)  ! force on particle iAtom
      bPotential(iAtom-iAtomsBegin+1) = bPotential(iAtom-iAtomsBegin+1) + potA
#ifdef COMP_FORCE
      bForce(1:3,iAtom-iAtomsBegin+1) = bForce(1:3,iAtom-iAtomsBegin+1) + forceA(1:3)
#endif
#ifdef FSUM
      if(sList(cBox)%group > 0) rankFSum(1:3) = rankFSum(1:3) + forceA(1:3) * bQ(iAtom)
#endif
    ! compute electrostatic energy
#ifdef DEBUG_PRINT_ENERGY
    fmmEnergy = fmmEnergy + bQ(iAtom) * bPotential(iAtom-iAtomsBegin+1)
#endif
  enddo ! iAtom
#endif

#ifndef MAX_PERF
#ifdef DEBUG_PRINT_ENERGY
#ifdef MPI
  call MPI_Reduce(fmmEnergy, esum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  fmmEnergy = esum
#endif
  energy = fmmEnergy * 0.5d0 * eConv
  if(iam == prtRank) write(*,'(a,23x,es26.18,a)') "mul =", energy, " Step 6"
#endif
#endif ! MAX_PERF

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
  timeConsumedSum = timeConsumed * nRanks
#ifndef MAX_PERF
#ifdef MPI
#if defined(USE_REDUCE) || defined(TIME_STEP6)
  if(.not.debugEachMPIrank) then
    call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
    call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
    call MPI_Reduce(timeConsumed, timeConsumedSum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  endif
#else
    call MPI_Barrier(MPI_COMM_WORLD, mpiError)
    timeConsumed = TimerCheck(3)
    timeConsumedMin = timeConsumed
    timeConsumedMax = timeConsumed
    timeConsumedSum = timeConsumed * nRanks
#endif ! USE_REDUCE || TIME_STEP6
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f7.1,a),e18.10,a)') "FMM Step 6: Timer: near-field pot:  min:", &
                             TimeConsumedMin, "s   max:", TimeConsumedMax,"s   sum:", TimeConsumedSum, "s"
  endif
  call TimerInit(3)
#endif ! TIMER

  !call print_double(bPotential, bForce, 5, "Step 6")
  !
  ! FMM step 7
  !
  ! Compute box self interaction: self potential and self force
  ! Work with the boxes at the lowest hierarchy level
  ! Deal with explicit particle-particle interaction
  ! Compute total electrostatic energy from potential
  !
  fmmEnergy = 0.0d0
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU
  !$omp target teams distribute parallel do device(iDevice) &
#endif

#ifndef CPLR_NVFTN
  !$omp has_device_addr(bPotential) &
#ifdef COMP_FORCE
  !$omp has_device_addr(bForce) &
#endif
#endif

  !$omp reduction(+:fmmEnergy) &
#ifdef FSUM
  !$omp reduction(+:rankFSum) &
#endif
  !$omp private(ix,iy,iz, iBox, cBox, potA, forceA, jAtom, rij, r2, pot, potB, forceB, brij, br2, bpot, aAtom) 
  !do iAtom = sList(sBoxLLstart)%atomStart, sList(nBoxesSList-1)%atomEnd
  do iAtom = iAtomsBegin, iAtomsEnd   ! (local to the rank) atom indices of S-boxes (LL) in the (bC) pool of (R,S,N) particles

    ! find the global box index for iAtom from  atom coordinates
    ix = int((bC(1,iAtom) - vertexMin(1)) * edgeReciprocal)
    iy = int((bC(2,iAtom) - vertexMin(2)) * edgeReciprocal)
    iz = int((bC(3,iAtom) - vertexMin(3)) * edgeReciprocal)

    ! for nRepLayers = 0, iBox = 4 * iz + 2 * iy + ix + 9
    ! supported values of nRepLayers = 0, 1, 3, 7, 15, 31, 63, and 127
    ! find which iBox the particle iAtom belongs to (works for cubic systems only)            loLevel
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
 !(8_8**loLevel - 1)/7_8  ! the total number of boxes for a system that has loLevel-1

    ! get the cBox (S-list) from iBox (global B-list)
    !cBox = iBox - nBoxesBefore_loLevel - sBoxesTopDispl(iam) * 8_8**rLevel &
    !     + ((8_8**rLevel - 1_8)/7_8 * sBoxesTopNum(iam) + 1_8) + nFloatingBoxes
    cBox = iBox - cBoxShift

    potA = 0.0d0
    forceA = 0.0d0
    !!$omp parallel do reduction(+:potA,forceA) private(rij, r2, pot)
    ! TODO replace with sList(iBox)%atomStart
    !do jAtom = box(iBox)%atomStart, iAtom-1
    do jAtom = sList(cBox)%atomStart, iAtom-1
      rij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
      r2       = rij(1)**2 + rij(2)**2 + rij(3)**2
      pot      = bQ(jAtom) / sqrt(r2)
      potA     = potA + pot
      forceA(1:3)= forceA(1:3) + pot / r2 * rij(1:3)
      !potA = potA + bQ(jAtom) / sqrt( (bC(1,iAtom)-bC(1,jAtom))**2 + &
      !                                (bC(2,iAtom)-bC(2,jAtom))**2 + &
      !                                (bC(3,iAtom)-bC(3,jAtom))**2 )
    enddo

    potB = 0.0d0
    forceB = 0.0d0
    !!$omp parallel do reduction(+:potB,forceB) private(rij, r2, pot)
    !do jAtom = iAtom+1, box(iBox)%atomEnd
    do jAtom = iAtom+1, sList(cBox)%atomEnd
      brij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
      br2       = brij(1)**2 + brij(2)**2 + brij(3)**2
      bpot      = bQ(jAtom) / sqrt(br2)
      potB     = potB + bpot
      forceB(1:3)= forceB(1:3) + bpot / br2 * brij(1:3)
      !potB = potB + bQ(jAtom) / sqrt( (bC(1,iAtom)-bC(1,jAtom))**2 + &
      !                                (bC(2,iAtom)-bC(2,jAtom))**2 + &
      !                                (bC(3,iAtom)-bC(3,jAtom))**2 )
    enddo

    ! accumulate potential and force in position of particle iAtom
    !aPotential(iAtom) = aPotential(iAtom) + (potA + potB)
    !aForce(1:3,iAtom) = bQ(iAtom) * (aForce(1:3,iAtom) + forceA(1:3) + forceB(1:3))
    aAtom = iAtom - iAtomsBegin + 1
    bPotential(aAtom) = bPotential(aAtom) + (potA + potB)
#ifdef COMP_FORCE
    bForce(1:3,aAtom) = bQ(iAtom) * (bForce(1:3,aAtom) + forceA(1:3) + forceB(1:3))
#endif
#ifdef FSUM
    if(sList(cBox)%group > 0) rankFSum(1:3) = rankFSum(1:3) + (forceA(1:3) + forceB(1:3)) * bQ(iAtom)
#endif

    ! compute electrostatic energy
    !sump = sump + bQ(iAtom) * aPotential(iAtom)
    fmmEnergy = fmmEnergy + bQ(iAtom) * bPotential(aAtom)
  enddo  ! loop over iAtom

#ifdef MPI
  call MPI_Reduce(fmmEnergy, esum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  fmmEnergy = esum
#endif
  energy = fmmEnergy * 0.5d0 * eConv

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "FMM Step 7: Timer: box-self energy: min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER

  call flush(6)

  timeConsumed = TimerCheck(2)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  !call MPI_Barrier(MPI_COMM_WORLD, mpiError)  ! we already have a barrier above for fmmEnergy
  timeConsumed = TimerCheck(2)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(/a)') "Total Energy"
    write(*,'(a,23x,es26.18,a,2(f9.1,a))') "mul =", energy, "    FMM time = min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    write(*,'(a)') "--- The FMM computation has concluded ---"
  endif
  call TimerInit(2)
  call TimerInit(3)

#ifdef FSUM
#if defined(OFFL_GPU)
  !$omp target exit data map(from: rankFSum)
#endif
  totalFSum = rankFSum
#ifdef MPI
  call MPI_Reduce(rankFSum, totalFSum, 3, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
#endif
  if(iam == leadRank) then
    write(*,'(/a)') "FMM: Electric force acting on a selected group of boxes, in nano-Newtons:"
      write(*,'(a,4x,a,3E22.12,4x,a,E22.12)') &
        "Group: 1", "Force: ", totalFSum(1:3) * fConv, &  ! in nano-Newtons
        "|F|: ", sqrt( totalFSum(1)**2 + totalFSum(2)**2 + totalFSum(3)**2 ) * fConv
  endif ! if(iam == leadRank)
#endif

#ifdef FSEL
  !purpose: compute molecular properties of the system

  ! this is a data analysis section
  ! sum up the the force for each LL-box of S-List from the atomic forces
  ! the resulting array boxFsum stays in the device (not returned)

#ifdef USE_DEVPTR
  ! boxFsum holds the force (in solid-body approximation) acting on each LL box
  ! boxFsum provides a mechanism for checking the accuracy of computation
  ! make the array boxFsum resident (allocated) in GPU
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):boxFsum)
#endif
  allocate(boxFsum(1:3,0:nBoxesSList-sBoxLLstart-1), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
#ifndef USE_DEVPTR
  ! this applies when the array is not GPU-resident
  boxFsum = 0.0d0
#endif

  ! return a copy of boxFsum from the device to host in array hostBoxFsum
  ! TODO skip this if nGroup == 0
  allocate(hostBoxFsum(1:3,0:nBoxesSList-sBoxLLstart-1), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
  hostBoxFsum = 0.0d0

  ! compute boxFsum(:) from bForce(:) on GPU by summing up atomic contributions
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU

#ifdef CPLR_NVFTN
  !$omp target teams distribute parallel do device(iDevice) &
#else
  !$omp target teams distribute parallel do has_device_addr(bForce) has_device_addr(boxFsum) device(iDevice) &
#endif
  !TODO do not need hostBoxFsum if nGroups == 0
  !$omp map(from: hostBoxFsum) &

#endif
  !$omp private(aBox, aAtom)
  do cBox = sBoxLLstart, nBoxesSList-1  ! Loop over S-list LL-boxes; cBox is an S-List index
    aBox = cBox - sBoxLLstart  ! aBox is a local (0-based) index for LL-boxes of S-List
    !
    ! sum up the force from all atoms of cBox (which is also aBox)
    boxFsum(1:3,aBox) = 0.0d0
    do aAtom = sList(cBox)%atomStart - iAtomsBegin + 1, sList(cBox)%atomEnd - iAtomsBegin + 1
      boxFsum(1:3,aBox) = boxFsum(1:3,aBox) + bForce(1:3,aAtom)
    enddo
    hostBoxFsum(1:3,aBox) = boxFsum(1:3,aBox)
  enddo

  !purpose: contribute the rank portion to the force-sum of selected (global) groups of LL-boxes
  ! TODO we can skip the computation of boxFsum (above) when nSelectedLLboxesRank == 0 and lComputeDirectSum == .F.

  ! nGroups is the global number of selected groups in the entire system
  ! boxGroupFsumRank holds the force sumed up over the selected LL-boxes in the rank
  ! TODO allocate boxGroupFsumRank(1:3,0:nGroups) to skip checking sList(cBox)%group > 0
  if(nGroups > 0) then
    allocate(boxGroupFsumRank(1:3,nGroups), stat=allocStat)
    success = allocStat
    if(success /= 0) call exitPrg(2)
    boxGroupFsumRank = 0.0d0

    ! nAtomsGroupRank(0:nGroups) is the number of atoms in each group in the rank
    allocate(nAtomsGroupRank(0:nGroups), nAtomsGroupGlobal(0:nGroups), stat=allocStat)
    success = allocStat
    if(success /= 0) call exitPrg(2)
    nAtomsGroupRank   = 0
    nAtomsGroupGlobal = 0

    ! nSelectedLLboxesRank is the number of selected LL-boxes in the rank
    ! Group 0 is everything that is not selected

    if(nSelectedLLboxesRank > 0) then  ! Warning: Disable this check to sum up the number of atoms in group 0
      ! accumulate the rank contribution to the box group in boxGroupFsumRank
      do cBox = sBoxLLstart, nBoxesSList-1  ! Loop over S-list LL-boxes; cBox is an S-List index
        aBox = cBox - sBoxLLstart  ! aBox is a local (0-based) index for LL-boxes of S-List
        !
        !TODO allocate boxGroupFsumRank(1:3,0:nGroups) to skip checking sList(cBox)%group > 0
        if(sList(cBox)%group > 0) &
        boxGroupFsumRank(1:3,sList(cBox)%group) = boxGroupFsumRank(1:3,sList(cBox)%group) + hostBoxFsum(1:3,aBox)
        !
        ! compute the number of atoms in each group
        nAtomsGroupRank(sList(cBox)%group) = nAtomsGroupRank(sList(cBox)%group) + &
          sList(cBox)%atomEnd - sList(cBox)%atomStart + 1
      enddo
    else
      boxGroupFsumRank = 0.0d0  ! need it for global reduction; TODO skip this if having a sub-communicator
    endif

    !purpose: reduce boxGroupFsumRank across all ranks
    !
    ! allocate the destination array in leadRank
    allocate(boxGroupFsumGlobal(1:3,nGroups), stat=allocStat)
    success = allocStat
    if(success /= 0) call exitPrg(2)
    boxGroupFsumGlobal = 0.0d0

    call MPI_Reduce(boxGroupFsumRank, boxGroupFsumGlobal, nGroups*3, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    call MPI_Reduce(nAtomsGroupRank,  nAtomsGroupGlobal,  nGroups+1, MPI_INTEGER8, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)

    if(iam == leadRank) then
      write(*,'(/a)') "The number of atoms in each group:"
      do i = 1, nGroups  ! Warning: To print from i=0, need to disable the check if(nSelectedLLboxesRank > 0) above
        if(nAtomsGroupSelected(i) - nAtomsGroupGlobal(i) == 0) then  ! sanity check
          write(*,'(a,i0,4x,i15)') "Group: ", i, nAtomsGroupGlobal(i)
        else
          write(*,'(a,i0,4x,2i15)')"Error: nAtomsGroupGlobal(i) and nAtomsGroupSelected(i) do not match", &
          i, nAtomsGroupGlobal(i), nAtomsGroupSelected(i)
        endif
      enddo

      write(*,'(/a)') "FMM: Electric force acting on the selected groups of boxes, in nano-Newtons:"
      do i = 1, nGroups
        write(*,'(a,i0,4x,a,3E22.12,4x,a,E22.12)') &
          "Group: ", i, "Force: ", boxGroupFsumGlobal(1:3,i) * fConv, &  ! in nano-Newtons
          "|F|: ", sqrt( &
            boxGroupFsumGlobal(1,i)**2 + boxGroupFsumGlobal(2,i)**2 + boxGroupFsumGlobal(3,i)**2  &
            ) * fConv
      enddo
    endif ! if(iam == leadRank)

#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(/a,2(f6.1,a))') "Phys prop:  Timer: force-sum boxes: min: ", &
                              TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER

  endif  ! if(nGroups > 0)

#endif  ! FSEL

  ! step non-FMM 

  call TimerInit(3)  ! reset the break-down timer for std (reference) computation

  if(nAtomsTotal ==           5943) then          !  0 wbox40 direct summation
    refEnergy = -4.146567836684986833E+05
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==       7926) then          !  0 wbox44 direct summation
    refEnergy = -5.532615329520769883E+05
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==      47544) then          !  1 wbox40 direct summation
    refEnergy = -3.317412143273170572E+06
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==      63408) then          !  1 wbox44 direct summation
    refEnergy = -4.426240151916897856E+06
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==     380352) then          !  3 wbox40 direct summation
    refEnergy = -2.654005647719660401E+07
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==     507264) then          !  3 wbox44 direct summation
    refEnergy = -3.541086834999790788E+07
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==    3042816) then          !  7 wbox40 direct summation
    refEnergy = -2.123236267029418945D+08
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==    4058112) then          !  7 wbox44 direct summation
    refEnergy = -2.832913418500415087E+08
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==   24342528) then          ! 15 wbox40 direct summation
    refEnergy = -1.698601821201063871D+09
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==   32464896) then          ! 15 wbox44 direct summation
    refEnergy = -2.266349322607834339E+09
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==  194740224) then          ! 31 wbox40 direct summation
    refEnergy = -1.358886578240979576E+10
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==  259719168) then          ! 31 wbox44 direct summation      48 nodes 4046 sec
    refEnergy = -1.813087058869335175E+10
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 1557921792) then          ! 63 wbox40
    refEnergy = -1.087111307264077148E+11         ! from FMM computation with cube=10 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 2077753344) then          ! 63 wbox44
    refEnergy = -1.450472714123668823E+11         ! from FMM computation with cube=22 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 12463374336_8) then       ! 127 wbox40
    refEnergy = -8.696898623406467285E+11         ! from FMM computation with cube=20 nL=30 ws=2
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 16622026752_8) then       ! 127 wbox44
    refEnergy = -1.160379402257002686E+12         ! from FMM computation with cube=22 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 99706994688_8) then       ! 255 wbox40
    refEnergy = -6.957522161154004883E+12         ! from FMM computation with cube=10 nL=30 ws=3 sc=256 on 256 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 132976214016_8) then      ! 255 wbox44
    refEnergy = -9.283040147882595703E+12         ! from FMM computation with cube=11 nL=30 ws=3 sc=256 on 256 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 797655957504_8) then      ! 511 wbox40
    refEnergy = -5.566019032968338281E+13         ! from FMM computation with cube=20 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 1063809712128_8) then     ! 511 wbox44
    refEnergy = -7.426434090985431250E+13         ! from FMM computation with cube=22 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 8510477697024_8) then     ! 1023 wbox44
    refEnergy = -5.941148061920701250E+14         ! from FMM computation with cube=22 nL=30 ws=3
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==21285857918976_8) then     ! 2047 wbox30
    refEnergy = -1.482080710327435750E+15         ! from FMM computation with cube=30 nL=30 ws=2 on 5462 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 36616830) then            ! aon4wbox replication 0
    refEnergy = -2.544089061151780128E+09         ! from direct summation 12 nodes sub-comm 6
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==     292934640_8) then     ! aon4wbox replication 1
    refEnergy = -2.035264383158514404E+10         ! from direct summation 96 nodes sub-comm 24
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==    2343477120_8) then     ! aon4wbox replication 3
    refEnergy = -1.628208773400063477E+11         ! from FMM computation with cube=22 nL=30 ws=3 sc=1   on   1 node
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==   18747816960_8) then     ! aon4wbox replication 7
    refEnergy = -1.302565929493822266E+12         ! from FMM computation with cube=44 nL=30 ws=3 sc=32  on  64 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==  149982535680_8) then     ! aon4wbox replication 15
   !refEnergy = -1.042052308995316602E+13         ! from FMM computation with cube=11 nL=30 ws=3 sc=64  on 512 nodes
    refEnergy = -1.042052308995306055E+13         ! from FMM computation with cube=22 nL=30 ws=3 sc=32  on 128 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 1199860285440_8) then     ! aon4wbox replication 31
   !refEnergy = -8.336416736327840625E+13         ! from FMM computation with cube=11 nL=14 ws=3
    refEnergy = -8.336416736325871875E+13         ! from FMM computation with cube=22 nL=30 ws=3 sc=512 on 512 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==        4844847_8) then    ! aon10t replication 0
    refEnergy = -3.115430542882155180E+08         ! from direct summation 1 node 56 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==       38758776_8) then    ! aon10t replication 1
    refEnergy = -2.492305059979503632E+09         ! from direct summation 4 nodes 1009 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==      130810869_8) then    ! aon10t replication 2
    refEnergy = -8.411481388859910965E+09         ! from direct summation 12 nodes 4292 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==      310070208_8) then    ! aon10t replication 3
    refEnergy = -1.993826755983322525E+10         ! from direct summation 80 nodes 3882 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==     2480561664_8) then    ! aon10t replication 7
    refEnergy = -1.595054211016271057E+11         ! from FMM computation with cube=22 nL=30 ws=3 sc=22 on 44 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==    19844493312_8) then    ! aon10t replication 15
    refEnergy = -1.276040440982466797E+12         ! from FMM computation with cube=22 nL=30 ws=3 sc=22 on 44 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==   158755946496_8) then    ! aon10t replication 31
    refEnergy = -1.020831172768275781E+13         ! from FMM computation with cube=22 nL=30 ws=3 sc=22 on 88 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal ==  1270047571968_8) then    ! aon10t replication 63
    refEnergy = -8.166644646848773437E+13         ! from FMM computation with cube=22 nL=30 ws=3 sc=342 on 684 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 10160380575744_8) then    ! aon10t replication 127
   !refEnergy = -6.533313820842202500E+14         ! from FMM computation with cube=22 nL=29 ws=3 sc=2731 on 2731 nodes
    refEnergy = -6.533313820842210000E+14         ! from FMM computation with cube=22 nL=30 ws=3 sc=2731 on 2731 nodes
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 3004287_8) then           ! aon320 rep=0
    refEnergy = -1.840954312663380802E+08         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 24034296_8) then          ! aon320 rep=1
    refEnergy = -1.472710699437767506E+09         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 81115749_8) then          ! aon320 rep=2
    refEnergy = -4.970331150663422585E+09         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 192274368_8) then         ! aon320 rep=3
    refEnergy = -1.178144258446423149E+10         ! from direct summation, 64 nodes for 1976.4 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 1538194944_8) then        ! aon320 rep=7
    refEnergy = -9.425051490392872620E+10         ! from direct summation, 1024 nodes for 10865.5 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 12305559552_8) then       ! aon320 rep=15
    refEnergy = -7.539999276312036133E+11         ! cu=40 nL=30 ws=2, 32 nodes for 444.3 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 41531263488_8) then       ! aon320 rep=23
    refEnergy = -2.544745005642558594E+12         ! cu=30 nL=30 ws=2 nr=8 nn=16 nr=8 sc=16 dd=8 took 2034.7 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 98444476416_8) then       ! aon320 rep=31
    refEnergy = -6.031982520031993164E+12         ! cu=40 nL=30 ws=2, 32 nodes for 3354.5 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 332250107904_8) then      ! aon320 rep=47
    refEnergy = -2.035792194471541797E+13         ! cu=30 nL=30 ws=2 nr=8 nn=128 sc=128 dd=9 took 2090.9 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 2658000863232_8) then     ! aon320 rep=95
    refEnergy = -1.628632229674037500E+14         ! cu=30 nL=30 ws=2 nr=8 nn=1024 sc=1024 dd=9 took 2095.4 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 6300446490624_8) then     ! aon320 rep=127
    refEnergy = -3.860460677047483125E+14         ! cu=40 nL=30 ws=2, 8192 nodes for 1486.1 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 21264006905856_8) then    ! aon320 rep=191
    refEnergy = -1.302905173279482500E+15         ! cu=30 nL=30 ws=2 nr=8 nn=8192 dd=9 took 2135.8 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 50403571924992_8) then    ! aon320 rep=255 50T
    refEnergy = -3.088367456426904000E+15         ! cu=40 nL=30 ws=2, 8192 nodes for 7251.5 sec with ddLevel=9
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 79504225418904_8) then    ! aon320 rep=297 79.5T
    refEnergy = -4.873834950759904000E+15         ! cu=46.5625 nL=20 ws=1, 8192 nodes for 5384.9 sec with ddLevel=10
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 91244185883136_8) then    ! aon320 rep=311 91T
    refEnergy = -5.583702230863062000E+15         ! cu=48.75 nL=30 ws=1 nn=10000 sc=10000 dd=11 pp=2 for 5773.8 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 100301870845176_8) then   ! aon320 rep=321 100T
    refEnergy = -6.143183033875440000E+15         ! cu=50.3125 nL=30 ws=1 nn=9216 sc=9216 dd=11 pp=1 for 8468.2 sec
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 3196149_8) then           ! wbox320 rep=0
    refEnergy = -2.243943094447086751E+08         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 25569192_8) then          ! wbox320 rep=1
    refEnergy = -1.795184396953133821E+09         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 86296023_8) then          ! wbox320 rep=2
    refEnergy = -6.058780700692528725E+09         ! from direct summation
    lrefEnergyAvailable = .TRUE.
  elseif(nAtomsTotal == 204553536_8) then         ! wbox320 rep=3
    refEnergy = -1.436159360647168541E+10         ! from direct summation, 32 nodes for 4837.5 sec
    lrefEnergyAvailable = .TRUE.
  else
    lrefEnergyAvailable = .FALSE.  ! the decision to compute direct sum depends on availability of coordinates, aC
  endif

  ! Error and accuracy analysis
  pdigits     = 0
  nErrors     = 0

  ! decide either to skip or compute the direct sum
  if(lComputeDirectSum) then

  ! use rank 0 to print the results of direct summation
  if(iam == prtRank) then
    write(*,'(/a)') "--- The next (non-FMM) step is validation of accuracy of the computation ---"
  endif

  ! Perform direct sum computation to validate the accuracy of FMM computation
  refEnergy   = 0.0d0
  difPotMax   = 0.0d0
  difPotRel   = 0.0d0
  difForceMax = 0.0d0
  difForceRel = 0.0d0
  rmsdP       = 0.0d0
  rmsdF       = 0.0d0
  averP       = 0.0d0
  averF       = 0.0d0

#if defined(OFFL_GPU)
  !$omp target enter data map(to: aC, aQ)
#endif

#ifdef FSEL
#ifdef USE_DEVPTR
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):rForce)
#endif
  allocate(rForce(1:3,rankAtoms(iam)), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
#ifndef USE_DEVPTR
  ! this applies when the array is not GPU-resident
  refForce = 0.0d0
#endif
#endif  ! FSEL

  ! Split the external loop in the double summation over nRanks groups
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU
  !$omp target teams distribute parallel do device(iDevice)      &
  !$omp map(tofrom: difPotMax,difPotRel,difForceMax,difForceRel) &
  !$omp map(tofrom: refEnergy,nErrors,averP,averF,rmsdP,rmsdF)   &
#endif

#ifndef CPLR_NVFTN
  !$omp has_device_addr(bPotential) &
#ifdef COMP_FORCE
  !$omp has_device_addr(bForce)     &
#endif
#ifdef FSEL
  !$omp has_device_addr(rForce)     &
#endif
#endif

  !$omp private(iAtom,aAtom,forceA,dsPotential,dsForce,jAtom,rij,r2,pot) &
  !$omp private(difPotential,difForce,relPotential,relForce) &
  !$omp reduction(max:difPotMax,difPotRel,difForceMax,difForceRel) &
  !$omp reduction(+:refEnergy,nErrors,averP,averF,rmsdP,rmsdF) 
  do iAtom = iAtomsBegin, iAtomsEnd

    aAtom = iAtom - iAtomsBegin + 1   ! aAtom = iAtom - nAtomsRank * iam   is a local index of the particle
    forceA      = 0.0d0
    dsPotential = 0.0d0
    dsForce     = 0.0d0

    do jAtom = 1, iAtom-1
      rij(1:3)    = aC(1:3,iAtom) - aC(1:3,jAtom)
      r2          = rij(1)**2 + rij(2)**2 + rij(3)**2
      pot         = aQ(jAtom) / sqrt(r2)
      forceA(1:3) = forceA(1:3) + pot / r2 * rij(1:3)
      dsPotential = dsPotential + pot
    enddo

    ! skip jAtom == iAtom since a charged particle does not interact with its own charge

    do jAtom = iAtom+1, nAtomsTotal
      rij(1:3)    = aC(1:3,iAtom) - aC(1:3,jAtom)
      r2          = rij(1)**2 + rij(2)**2 + rij(3)**2
      pot         = aQ(jAtom) / sqrt(r2)
      forceA(1:3) = forceA(1:3) + pot / r2 * rij(1:3)
      dsPotential = dsPotential + pot
    enddo

    dsForce(1:3)  = aQ(iAtom) * forceA(1:3)
#ifdef FSEL
    rForce(1:3,iAtom - iAtomsBegin + 1) = dsForce(1:3)  ! save the reference atomic force for further analysis
#endif
    refEnergy = refEnergy + aQ(iAtom) * dsPotential

    ! test the correctness of potential
    if(dabs(bPotential(aAtom) - dsPotential) > 0.01d0) nErrors = nErrors + 1

    ! determine absolute unsigned difference for each data point
    difPotential = dabs(bPotential(aAtom) - dsPotential)
    difForce = sqrt((bForce(1,aAtom)-dsForce(1))**2  &
                  + (bForce(2,aAtom)-dsForce(2))**2  &
                  + (bForce(3,aAtom)-dsForce(3))**2)
    averP = averP + difPotential
    averF = averF + difForce
    rmsdP = rmsdP + difPotential * difPotential
    rmsdF = rmsdF + difForce * difForce

    ! determine relative fractional (%) unsigned difference for each data point
    relPotential = dabs( (bPotential(aAtom) - dsPotential) / dsPotential )
    relForce = sqrt(((bForce(1,aAtom)-dsForce(1))**2  &
                   + (bForce(2,aAtom)-dsForce(2))**2  &
                   + (bForce(3,aAtom)-dsForce(3))**2) &
                  / ((dsForce(1)**2 + dsForce(2)**2 + dsForce(3)**2)))

    difPotMax   = dmax1(difPotMax, dabs(bPotential(aAtom) - dsPotential))  ! unsigned max error in absolute value
    difPotRel   = dmax1(difPotRel, relPotential)                           ! unsigned max error in relative value
    difForceMax = dmax1(difForceMax, sqrt((bForce(1,aAtom)-dsForce(1))**2 &
                                         +(bForce(2,aAtom)-dsForce(2))**2 &
                                         +(bForce(3,aAtom)-dsForce(3))**2) )
    difForceRel = dmax1(difForceRel, relForce)
  enddo

  if(nErrors /= 0) write(*,"(/2(a,i0))") "Error:  nErrors = ", nErrors, " in potential for rank ", iam

#ifdef MPI
  call MPI_Reduce(nErrors, nErrors_sum, 1, MPI_INTEGER, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  nErrors = nErrors_sum

  call MPI_Reduce(refEnergy, refEnergy_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  refEnergy = refEnergy_sum
#endif
  refEnergy = refEnergy * 0.5d0 * eConv
  lrefEnergyAvailable = .TRUE.
  endif  ! if(lComputeDirectSum)

  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF

  if(iam == prtRank .and. lrefEnergyAvailable) then
    write(*,'(/a)') "Reference value"
    write(*,'(a,23x,es26.18,a,2(f9.1,a))') "std =", refEnergy, "     DS time = min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
    call show_residual(energy, refEnergy, pdigits)
  endif
  call TimerInit(3)

  if(lComputeDirectSum) then
#ifdef MPI
    call MPI_Reduce(averP, averP_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    averP = averP_sum
    call MPI_Reduce(averF, averF_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    averF = averF_sum

    call MPI_Reduce(rmsdP, rmsdP_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    rmsdP = rmsdP_sum
    call MPI_Reduce(rmsdF, rmsdF_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    rmsdF = rmsdF_sum

    call MPI_Reduce(difPotMax, difPotMax_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difPotMax = difPotMax_sum
    call MPI_Reduce(difPotRel, difPotRel_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difPotRel = difPotRel_sum

    call MPI_Reduce(difForceMax, difForceMax_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difForceMax = difForceMax_sum
    call MPI_Reduce(difForceRel, difForceRel_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difForceRel = difForceRel_sum
#endif

    ! compute RMSD
    rmsdP = dsqrt(rmsdP / dble(nAtomsTotal))
    rmsdF = dsqrt(rmsdF / dble(nAtomsTotal))

    ! compute average error
    averP = averP / dble(nAtomsTotal)
    averF = averF / dble(nAtomsTotal)

    ! print the error analysis
    if(iam == prtRank) then
      ! print the histogramm for relative error
      !write(*,'(/a)') "Relative error distribution normalized for number of particles (% of particles in each bin)"
      !do i = 1, nhist
      !  write(*,'(a,i2,f20.16,2f8.2,4x,2i20)') "hist: ", i, err(i), dble(rPhistogram(i)*100)/dble(nAtomsTotal), &
      !                                         dble(rFhistogram(i)*100)/dble(nAtomsTotal), rPhistogram(i),rFhistogram(i)
      !enddo
      write(*,'(/a)') "Potential and Force at the atomic centers:"
      write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "potential maxAbsErr maxRelErr = ", &
           difPotMax, int(dabs(dlog10(difPotMax))), difPotRel, int(dabs(dlog10(difPotRel)))
      write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "force     maxAbsErr maxRelErr = ", &
           difForceMax, int(dabs(dlog10(difForceMax))), difForceRel, int(dabs(dlog10(difForceRel)))

      ! print the histogramm for absolute error
      !write(*,'(/a)') "Absolute unsigned error distribution normalized for number of particles (% of particles in each bin)"
      !do i = 1, nhist
      !  write(*,'(a,i2,f20.16,2f8.2,4x,2i20)') "hist: ", i, err(i), dble(aPhistogram(i)*100)/dble(nAtomsTotal), &
      !                                         dble(aFhistogram(i)*100)/dble(nAtomsTotal), aPhistogram(i),aFhistogram(i)
      !enddo
      write(*,'(/a,2x,2(f20.14," (",i0,")",1x))') "potential averErr rmsd        = ", &
           averP, int(dabs(dlog10(averP))), rmsdP, int(dabs(dlog10(rmsdP)))
      write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "force     averErr rmsd        = ", &
           averF, int(dabs(dlog10(averF))), rmsdF, int(dabs(dlog10(rmsdF)))
    endif
  endif  !   if(lComputeDirectSum)

  ! TODO need to test this part
  ! TODO the present implementation requires the handling of all particles by each rank that is not scalable
  ! TODO instead, rank 0 should gather the number of particles owned by each rank
#if defined(MPI) && defined(FETCH_POTFORCE) && !defined(BCAST)
  ! retrieve potential and force to the rank 0
  call MPI_Gatherv(bPotential, rankAtoms(iam), MPI_DOUBLE_PRECISION, &
                  aPotential, rankAtoms, dispAtoms, MPI_DOUBLE_PRECISION, &
                  leadRank, MPI_COMM_WORLD, mpiError)

  if(mpiError /= MPI_SUCCESS) call exitPrg(9)

  rankAtoms = rankAtoms * 3
  dispAtoms = dispAtoms * 3
  call MPI_Gatherv(bForce, rankAtoms(iam), MPI_DOUBLE_PRECISION, &
                  aForce, rankAtoms, dispAtoms, MPI_DOUBLE_PRECISION, &
                  leadRank, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  rankAtoms = rankAtoms / 3
  dispAtoms = dispAtoms / 3
#endif

  ! TODO need to test this part
#if defined(MPI) && defined(FETCH_POTFORCE) && defined(BCAST)
  ! copy the data from the temporary buffer to the right place in the destination array
  aPotential(iAtomsBegin : iAtomsEnd) = bPotential(1:rankAtoms(iam))
  aForce(1:3,iAtomsBegin : iAtomsEnd) = bForce(1:3,1:rankAtoms(iam))

  ! broadcast the sections of the destination array; each rank will hold a copy of the potential and force
  iAtomsBegin = 1
  do i = 0, nRanks-1
    call MPI_Bcast(aPotential(iAtomsBegin), rankAtoms(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
    call MPI_Bcast(aForce(1,iAtomsBegin), rankAtoms(i)*3, MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(6)
    iAtomsBegin = iAtomsBegin + rankAtoms(i)
  enddo
#endif

  ! TODO need to test this part
#if defined(MPI) && defined(FETCH_POTFORCE) && defined(TIMER)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#if defined (MPI) && defined (USE_REDUCE)
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#endif
#endif
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "Pot & Force Timer:  MPI Collective: min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif

  ! TODO need to test this part
#if defined(FETCH_POTFORCE)
  ! TODO: all unit conversions should be done by each rank prior sending the data to rank 0
  aPotential = aPotential * pConv  ! in Volts
  aForce = aForce * fConv          ! in nano-Newtons
#else
  ! TODO: bPotential and bForce are GPU-resident and do not exist in the host memory
  ! TODO: temporarily disabling the conversion
  ! TODO: do not forget to apply this conversion before the use (otherwise, in electronic units)
  !bPotential = bPotential * pConv  ! in Volts
  !bForce = bForce * fConv          ! in nano-Newtons
#endif

  ! Perform conventional calculation to compare its result to that of FMM computation
  ! Potential and force computation goes in parallel
  ! Only rank 0 preforms histogram analysis
!  call Demo_energy_pairwise_test(energy, aPotential, aForce, pdigits, hError)

  !
  ! compute lattice sums
  ! p-p: particle-particle interaction
  ! m-m: multipole-multipole interaction
  !
  !call latticesums2(energy)  ! p-p in near field
  !call latticesums4(energy)  ! m-m in near-field
  !call latticesums9(energy)  ! m-m in far-field
  !
  !call LatticeSumsPPinNF(energy)  ! particle-particle interaction in near-field of unit cell with its images
  !call LatticeSumsMMinNF(energy, shiftedExpansion, potentialExpansion)  ! m-m in near-field
  !call LatticeSumsMMinNFserial(energy)  ! m-m in near-field; serial implementation
  !call latticesums9(energy)  ! m-m in far-field
  !
!  write(*,'(/a)') "Periodic Boundary Condition"
!  call LatticeSumFF(energy, aPotential, aForce)  ! m-m in far-field
  !call latticesums10()  ! m-m in far-field
  !
  ! TODO need to test this part
#if defined (DEBUG) && defined (DETAILED)
  !
  ! Activate this section when dealing with PBC
  !
  call TimerInit(3)
  energyStd = 0.0d0
  kzLoop: do kz = -ws, ws            ! generate images along z-axis
    zShift = treeLevel(0)%edge * dble(kz)
    kyLoop: do ky = -ws, ws          ! generate images along y-axis
      yShift = treeLevel(0)%edge * dble(ky)
      kxLoop: do kx = -ws, ws        ! generate images along x-axis
        if(kx==0 .and. ky==0 .and. kz==0) cycle kxLoop   ! skip the unit-cell
        xShift = treeLevel(0)%edge * dble(kx)
        do jAtom = 1, nAtomsTotal    ! loop over atoms in the replicated image box
          bx = bC(1,jAtom) + xShift
          by = bC(2,jAtom) + yShift
          bz = bC(3,jAtom) + zShift
          do iAtom = 1, nAtomsTotal  ! loop over atoms in the central (unit-cell) box
            energyStd = energyStd + bQ(iAtom) * bQ(jAtom) / sqrt( (bC(1,iAtom)-bx)**2 + &
                                                                  (bC(2,iAtom)-by)**2 + &
                                                                  (bC(3,iAtom)-bz)**2 ) 
          enddo
        enddo
      enddo kxLoop
    enddo kyLoop
  enddo kzLoop
  energyStd = energyStd * eConv
  if(iam == prtRank) then
    write(*,'(a,f24.10,a,f9.1,a)') "std = ", energyStd, " Std pbc near-field time = ",TimerCheck(3),"s"
  endif
#endif


#ifdef FSEL
  if(lComputeDirectSum) then
  !
  !purpose: compute refFsum for each LL-box in S-list

#ifdef USE_DEVPTR
  !$omp allocators allocate (allocator(omp_target_device_mem_alloc):refFsum)
#endif
  allocate(refFsum(1:3,0:nBoxesSList-sBoxLLstart-1), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
#ifndef USE_DEVPTR
  ! this applies when the array is not GPU-resident
  refFsum = 0.0d0
#endif

  ! return a copy of refFsum from the device to host in array hostRefFsum
  allocate(hostRefFsum(1:3,0:nBoxesSList-sBoxLLstart-1), stat=allocStat)
  success = allocStat
  if(success /= 0) call exitPrg(2)
  hostRefFsum = 0.0d0
  !
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU

#ifdef CPLR_NVFTN
  !$omp target teams distribute parallel do device(iDevice) &
#else
  !$omp target teams distribute parallel do has_device_addr(rForce) has_device_addr(refFsum) device(iDevice) &
#endif
  ! TODO do this only when nGroups > 0
  !$omp map(from: hostRefFsum) &

#endif
  !$omp private(aBox, aAtom)
  do cBox = sBoxLLstart, nBoxesSList-1  ! Loop over S-list LL-boxes; cBox is an S-List index
    aBox = cBox - sBoxLLstart  ! aBox is a local (0-based) index for LL-boxes of S-List
    !
    ! sum up the force from all atoms of cBox (which is also aBox)
    refFsum(1:3,aBox) = 0.0d0
    do aAtom = sList(cBox)%atomStart - iAtomsBegin + 1, sList(cBox)%atomEnd - iAtomsBegin + 1
      refFsum(1:3,aBox) = refFsum(1:3,aBox) + rForce(1:3,aAtom)
    enddo
    hostRefFsum(1:3,aBox) = refFsum(1:3,aBox)
  enddo

  !purpose: compare boxFsum against refFsum for each LL-box in S-list
  !
  averF       = 0.0d0
  rmsdF       = 0.0d0
  difForceMax = 0.0d0
  difForceRel = 0.0d0

#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU

#ifdef CPLR_NVFTN
  !$omp target teams distribute parallel do device(iDevice) &
#else
  !$omp target teams distribute parallel do has_device_addr(boxFsum) has_device_addr(refFsum) device(iDevice) &
#endif

  !$omp map(tofrom: averF,rmsdF,difForceMax,difForceRel)    &
#endif

  !$omp private(difForce,relForce) &
  !$omp reduction(+:averF,rmsdF)   &
  !$omp reduction(max:difForceMax,difForceRel)
  do aBox = 0, nBoxesSList - sBoxLLstart - 1  ! aBox is a local (0-based) index for LL-boxes of S-List
    !
    ! determine absolute unsigned difference for each data point
    difForce = sqrt((boxFsum(1,aBox) - refFsum(1,aBox))**2  &
                  + (boxFsum(2,aBox) - refFsum(2,aBox))**2  &
                  + (boxFsum(3,aBox) - refFsum(3,aBox))**2)
    averF = averF + difForce
    rmsdF = rmsdF + difForce * difForce

    ! determine relative fractional (%) unsigned difference for each data point
    relForce = sqrt(((boxFsum(1,aBox) - refFsum(1,aBox))**2  &
                   + (boxFsum(2,aBox) - refFsum(2,aBox))**2  &
                   + (boxFsum(3,aBox) - refFsum(3,aBox))**2) &
                  / ((refFsum(1,aBox)**2 + refFsum(2,aBox)**2 + refFsum(3,aBox)**2)))

    difForceMax = dmax1(difForceMax, sqrt((boxFsum(1,aBox) - refFsum(1,aBox))**2  &
                                        + (boxFsum(2,aBox) - refFsum(2,aBox))**2  &
                                        + (boxFsum(3,aBox) - refFsum(3,aBox))**2) )
    difForceRel = dmax1(difForceRel, relForce)
  enddo

#ifdef MPI
    call MPI_Reduce(averF, averF_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    averF = averF_sum

    call MPI_Reduce(rmsdF, rmsdF_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    rmsdF = rmsdF_sum

    call MPI_Reduce(difForceMax, difForceMax_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difForceMax = difForceMax_sum

    call MPI_Reduce(difForceRel, difForceRel_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    difForceRel = difForceRel_sum
#endif

    ! compute RMSD
    rmsdF = dsqrt(rmsdF / dble(nAtomsTotal))

    ! compute average error
    averF = averF / dble(nAtomsTotal)

    ! print the error analysis
    if(iam == prtRank) then
      write(*,'(/a)') "Force acting on LL boxes in the rigid-body approximation:"

      write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "force     maxAbsErr maxRelErr = ", &
           difForceMax, int(dabs(dlog10(difForceMax))), difForceRel, int(dabs(dlog10(difForceRel)))

      write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "force     averErr rmsd        = ", &
           averF, int(dabs(dlog10(averF))), rmsdF, int(dabs(dlog10(rmsdF)))
    endif

  !purpose: contribute the rank portion to the force-sum of selected (global) groups of LL-boxes
  ! TODO we can skip the computation of refFsum (above) when nSelectedLLboxesRank == 0 and lComputeDirectSum == .F.

  if(nGroups > 0) then
    allocate(refGroupFsumRank(1:3,nGroups), stat=allocStat)
    success = allocStat
    if(success /= 0) call exitPrg(2)
    refGroupFsumRank = 0.0d0

    if(nSelectedLLboxesRank > 0) then
      ! accumulate the rank contribution for the box group into refGroupFsumRank
      do cBox = sBoxLLstart, nBoxesSList-1  ! Loop over S-list LL-boxes; cBox is an S-List index
        aBox = cBox - sBoxLLstart  ! aBox is a local (0-based) index for LL-boxes of S-List
        if(sList(cBox)%group > 0) &
        refGroupFsumRank(1:3,sList(cBox)%group) = refGroupFsumRank(1:3,sList(cBox)%group) + hostRefFsum(1:3,aBox)
      enddo
    else
      refGroupFsumRank = 0.0d0  ! need it for global reduction; TODO skip this if having a sub-communicator
    endif

    !purpose: reduce refGroupFsumRank across all ranks
    !
    ! allocate the destination array in leadRank
    allocate(refGroupFsumGlobal(1:3,nGroups), stat=allocStat)
    success = allocStat
    if(success /= 0) call exitPrg(2)
    refGroupFsumGlobal = 0.0d0

    call MPI_Reduce(refGroupFsumRank, refGroupFsumGlobal, 3*nGroups, MPI_DOUBLE_PRECISION, MPI_SUM, leadRank, MPI_COMM_WORLD, mpiError)

    if(iam == leadRank) then
      write(*,'(/a)') "Direct Summation: Electric force acting on selected groups of boxes:"
      do i = 1, nGroups
        write(*,'(a,i0,4x,a,3E22.12)') "Group: ", i, "Force: ", refGroupFsumGlobal(1:3,i) * fConv  ! in nano-Newtons
      enddo

      write(*,'(/a)') "Error Analysis of group force: Relative (FMM - DS) and Absolute |FMM - DS|:"
      do i = 1, nGroups
        write(*,'(a,i0,4x,a,3E22.12,4x,E22.12)') "Group: ", i, "Force: ", &
        (boxGroupFsumGlobal(1:3,i) - refGroupFsumGlobal(1:3,i)) * fConv,  & ! in nano-Newtons
        sqrt( &
                (boxGroupFsumGlobal(1,i) - refGroupFsumGlobal(1,i))**2 &
              + (boxGroupFsumGlobal(2,i) - refGroupFsumGlobal(2,i))**2 &
              + (boxGroupFsumGlobal(3,i) - refGroupFsumGlobal(3,i))**2 &
            ) * fConv
      enddo
    endif

  endif  ! if(nGroups > 0)

  endif  !   if(lComputeDirectSum)
#endif  ! FSEL

  !
  ! deallocate device memory
#ifdef COMP_FORCE
#if defined(OFFL_GPU) && defined(USE_DEVPTR)
  !$omp target exit data map(delete: bPotential, bForce)
#else
  deallocate(bPotential, bForce, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
#endif

#else
#if defined(OFFL_GPU) && defined(USE_DEVPTR)
  !$omp target exit data map(delete: bPotential)
#else
  deallocate(bPotential, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
#endif
#endif  ! ifdef COMP_FORCE

  !
  ! deallocate arrays involved in the computation of force-sum
#ifdef FSEL

#if defined(OFFL_GPU) && defined(USE_DEVPTR)
  !$omp target exit data map(delete: boxFsum)
#else
  deallocate(boxFsum, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
#endif

  deallocate(hostBoxFsum, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  if(nGroups > 0) then
    deallocate(boxGroupFsumRank, boxGroupFsumGlobal, nAtomsGroupRank, nAtomsGroupGlobal, stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)
  endif

  if(lComputeDirectSum) then
#if defined(OFFL_GPU) && defined(USE_DEVPTR)
    !$omp target exit data map(delete: refFsum, rForce)
#else
    deallocate(refFsum, rForce, stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)
#endif

    deallocate(hostRefFsum, stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)

    if(nGroups > 0) then
      deallocate(refGroupFsumRank, refGroupFsumGlobal, stat=allocStat)
      if(allocStat /= 0) call exitPrg(2)
    endif

  endif  ! if(lComputeDirectSum)
 
#endif  !FSEL

  if(lComputeDirectSum) then
#ifdef TIMER
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#ifndef MAX_PERF
#ifdef MPI
#ifdef USE_REDUCE
  call MPI_Reduce(timeConsumed, timeConsumedMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, leadRank, MPI_COMM_WORLD, mpiError)
  call MPI_Reduce(timeConsumed, timeConsumedMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, leadRank, MPI_COMM_WORLD, mpiError)
#else
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
  timeConsumed = TimerCheck(3)
  timeConsumedMin = timeConsumed
  timeConsumedMax = timeConsumed
#endif ! USE_REDUCE
#endif ! MPI
#endif ! MAX_PERF
  if(iam == prtRank) then
    write(*,'(a,2(f6.1,a))') "DirectSum:  Timer:  error analysis: min: ", &
                             TimeConsumedMin, "s   max: ", TimeConsumedMax,"s"
  endif
  call TimerInit(3)
#endif ! TIMER
  endif  ! if(lComputeDirectSum)

  !
  ! deallocate target offloading
#if defined(OFFL_GPU) && !defined(EXIT_STEP1)
  !$omp end target data
#endif

  !
  ! deallocate host memory
  success = 0
  deallocate( &  !oShift, shiftedExpansion, irrExpansion, tmpExpansion, &
#ifdef M2LROT
             scrExpansion, &
             expalp, &
#endif
#if defined(FETCH_POTFORCE)
             aPotential, &
             aForce, &
#endif
             bPotential, &
#ifdef COMP_FORCE
             bForce, &
#endif
#ifdef DYNAMICLIST
             potentialExpansion, &
             boxList, &
#endif
#if defined(USE_HOSTRDMA) && defined(OFFL_GPU)
             hostRegExpansion, &
#endif
             nAtomsGroupSelected, &  ! this array is always allocated
             boxIrrExpansion, &
             stat=allocStat)
  success = success + allocStat
  if(success /= 0) call exitPrg(2)
  !
  if(lComputeDirectSum .and. iam == prtRank .and. (pdigits < 8 .or. nErrors /= 0)) then
    write(*,'(/a,2l5)') "Error: computation produced numerically incorrect result: pdigits < 8 .or. nErrors /= 0", &
            (pdigits < 8), (nErrors /= 0)
    error stop 1
  endif
  !
end subroutine ElectrostaticEnergy


subroutine PrintExpansion(nodeStart, nodeEnd, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1, nbr, outexp, prtexp)
  !
  ! print first few elements of the multipole expansion of boxes
  !
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: prtRank, leadRank
  implicit none
  !
#ifdef MPI 
  integer          :: mpiError
  type(MPI_Status) :: status
#endif
  integer, intent(in) :: nodeStart, nodeEnd, iam, nRanks, szOut, szPrt, nSzPrt, nSzPrt1
  double complex, intent(in)  :: outexp(0:szOut)
  double complex, intent(out) :: prtexp(0:szPrt)
  integer,        intent(out) :: nbr(1:2, 0:nRanks)
  integer :: iRank, iBox, jBox, i, k

  if(nRanks == 1) then
    ! serial computation
    write(*,'(a,4i8)') "rank nodeStart nodeEnd nBoxes: ", iam, nodeStart, nodeEnd, nodeEnd-nodeStart+1
    do iBox = nodeStart, nodeEnd
      do i = 0, nSzPrt
        write(*,'(f20.8,2x,f20.8)') dble(outexp((iBox-nodeStart)*(nSzPrt1)+i)), &
                                   aimag(outexp((iBox-nodeStart)*(nSzPrt1)+i))
      enddo  
    enddo 
    ! TODO Not sure if the following is right
    nbr = 0
    prtexp = 0
#ifdef MPI
  else
    ! parallel computation using multiple MPI ranks
    do iRank = 0, nRanks-1
      if(iam == iRank .and. iam /= leadRank) then
        call MPI_Send(nodeStart, 1, MPI_INTEGER, leadRank, iRank, MPI_COMM_WORLD, mpiError)
      elseif(iam /= iRank .and. iam == leadRank) then 
        call MPI_Recv(iBox, 1, MPI_INTEGER, iRank, iRank, MPI_COMM_WORLD, status, mpiError)
        nbr(1,iRank) = iBox
      elseif(iam == iRank .and. iam == leadRank) then
        nbr(1,iRank) = nodeStart
      endif
    enddo
    
    do iRank = 0, nRanks-1
      if(iam == iRank .and. iam /= leadRank) then
        call MPI_Send(nodeEnd, 1, MPI_INTEGER, leadRank, iRank, MPI_COMM_WORLD, mpiError)
      elseif(iam /= iRank .and. iam == leadRank) then
        call MPI_Recv(jBox, 1, MPI_INTEGER, iRank, iRank, MPI_COMM_WORLD, status, mpiError)
        nbr(2,iRank) = jBox
      elseif(iam == iRank .and. iam == leadRank) then
        nbr(2,iRank) = nodeEnd
      endif
    enddo

    if(iam == leadRank) then
      do iRank = 0, nRanks-1
        write(*,'(a,4i8)') "rank nodeStart nodeEnd nBoxes: ", &
             iRank, nbr(1,iRank), nbr(2,iRank), nbr(2,iRank)-nbr(1,iRank)+1
      enddo
    endif

    do iRank = 0, nRanks-1
      if(iam == iRank .and. iam /= leadRank) then
        k = (nodeEnd-nodeStart+1)*nSzPrt1
        call MPI_Send(outexp, k, MPI_DOUBLE_COMPLEX, leadRank, iRank, MPI_COMM_WORLD, mpiError)
      elseif(iam /= iRank .and. iam == leadRank) then
        k = (nbr(2,iRank)-nbr(1,iRank)+1)*nSzPrt1
        call MPI_Recv(prtexp, k, MPI_DOUBLE_COMPLEX, iRank, iRank, MPI_COMM_WORLD, status, mpiError)
        do i = 0, k - 1
          write(*,'(f20.8,2x,f20.8)') dble(prtexp(i)), aimag(prtexp(i))
        enddo
      elseif(iam == iRank .and. iam == leadRank) then
        k = (nodeEnd-nodeStart+1)*nSzPrt1
        do i = 0, k - 1
          write(*,'(f20.8,2x,f20.8)') dble(outexp(i)), aimag(outexp(i))
        enddo
      endif
    enddo
#endif
  endif
end subroutine PrintExpansion


subroutine ElectricForceFF(iLevel, fmmEnergy, bPotential, bForce, tmpExpansion, rExpansion)
  !
  ! Use irregular FF expansion accumulated in the box
  ! Works with non-periodic and periodic FMM
  ! All computed values are in Gauss units (charge in electrons, coordinates in Angstroms)
  ! Computation of energy includes correction against double counting
  !
  use common_data, only: nL, tri, nSzTri, bC, bQ, vertexMin, nAtomsTotal, boxIrrExpansion, prtRank, leadRank
  use common_constants, only: zero, zeroCoor
  use fmm_plan, only: treeLevel, rList
  use timer
!$ use omp_lib
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  integer, intent(in)             :: iLevel      ! hierarchy level for the box carrying particles
  double precision, intent(inout) :: fmmEnergy   ! output; add FF energy to the initial value of variable
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri), rExpansion(0:nSzTri)        ! scratch arrays
  double precision, intent(inout) :: bPotential(nAtomsTotal), bForce(3,nAtomsTotal)  ! output
  !
  integer :: iam, nRanks
  integer :: i, L, aL1, m, iAtom, iBox
  double precision :: iOrigin(3)
  !double precision :: energy
  double precision :: x,y,z, x2y2, xDx2y2, yDx2y2, dL, dm, sump, sumf(3), aprod, bprod
  double precision :: pot, potA, potB, r2, rij(3), forceA(3), forceB(3)
  double complex   :: cprod
  integer :: boxesStart, boxesEnd, xyz(3)
  !
#ifdef MPI
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
#endif
  !
  ! compute far-field potential and force in the position of particles
  ! obtain force by differentiation of the above computed FF potential (in multipole expansion)
  ! compute FF contribution to the energy of the system
  !
  xyz = 0  ! dummy initialization to stop the compiler from issuing warning about uninitialized variable
  stop "Error in ElectricForceFF"
  !
  boxesStart = (8**(iLevel+1) - 1) / 7
  boxesEnd   = boxesStart + 8**iLevel - 1
  do iBox = boxesStart, boxesEnd
  !do iBox = treeLevel(iLevel)%bStart, treeLevel(iLevel)%bEnd
#ifdef MPI0
    if(iam == box(iBox)%proc) then
#endif
      ! find center of the box
      !iOrigin(1) = vertexMin(1) + (dble(box(iBox)%xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
      !iOrigin(2) = vertexMin(2) + (dble(box(iBox)%xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
      !iOrigin(3) = vertexMin(3) + (dble(box(iBox)%xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge
      iOrigin(1) = vertexMin(1) + (dble(xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
      iOrigin(2) = vertexMin(2) + (dble(xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
      iOrigin(3) = vertexMin(3) + (dble(xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge

      !! Not going to return the variable energy
      !!$omp parallel do reduction (+:energy)    &
      !$omp parallel do     &
      !$omp default(shared) &
      !$omp private(iAtom,x,y,z,x2y2,xDx2y2,yDx2y2,sump,sumf,L,dL,aL1,i,aprod,bprod,m,dm,cprod,tmpExpansion,rExpansion)
      do iAtom = rList(iBox)%atomStart, rList(iBox)%atomEnd

        x = bC(1,iAtom) - iOrigin(1)
        y = bC(2,iAtom) - iOrigin(2)
        z = bC(3,iAtom) - iOrigin(3)
        x2y2 = x*x + y*y
   
        ! modify the original coordinate by a very small increment to avoid division by zero
        if(x2y2 < zero) then
          bC(1,iAtom) = iOrigin(1) + dsign(zeroCoor,bC(1,iAtom))
          x = bC(1,iAtom) - iOrigin(1)
          x2y2 = x*x + y*y
        endif

        xDx2y2 = x / x2y2
        yDx2y2 = y / x2y2

        ! compute multipole expansion for radius-vector iAtom
        call RegularExpansion(bC(1,iAtom), 1.0d0, 1, iOrigin, tmpExpansion, rExpansion)

        ! in the calculation of electrostatic potential, there is no need to use an 
        ! additional expansion for chargeless multipoles; one can use charged multipoles
        ! in the calculation and divide the final values by atom charge

        ! L = 0 and m = 0
        ! sump is sum for potential; sumf is sum for force
        !sump = dble(rExpansion(0) * box(iBox)%potentialExpansion(0)) * 0.5d0  ! correct for default double counting
        sump = dble(rExpansion(0) * boxIrrExpansion(rList(iBox)%expStart+0)) * 0.5d0  ! correct for default double counting
        sumf = 0.0d0

        ! L > 0
        do L = 1, nL
          dL  = dble(L)
          aL1 = tri(L-1)
          i   = tri(L)

          ! m = 0
          ! potential; correct for default double counting
          !sump = sump + dble(rExpansion(i) * box(iBox)%potentialExpansion(i)) * 0.5d0
          sump = sump + dble(rExpansion(i) * boxIrrExpansion(rList(iBox)%expStart+i)) * 0.5d0

          ! force; divide by 2 so that index m can run only half-range (positive values)
          !aprod = dble( rExpansion(i)   * box(iBox)%potentialExpansion(i) ) * 0.5d0
          !bprod = dble( rExpansion(aL1) * box(iBox)%potentialExpansion(i) ) * 0.5d0
          aprod = dble( rExpansion(i)   * boxIrrExpansion(rList(iBox)%expStart+i) ) * 0.5d0
          bprod = dble( rExpansion(aL1) * boxIrrExpansion(rList(iBox)%expStart+i) ) * 0.5d0
          sumf(1) = sumf(1) + dL * xDx2y2 * aprod - z * xDx2y2 * bprod
          sumf(2) = sumf(2) + dL * yDx2y2 * aprod - z * yDx2y2 * bprod
          sumf(3) = sumf(3) + bprod

          ! do only half-range (only positive values); multiply by 2 to get the effect of the negative-range
          ! 0 < m < L
          do m = 1, L-1
            i = i + 1
            dm = dble(m)
            !bprod = dble( rExpansion(aL1+m) * box(iBox)%potentialExpansion(i) )
            !cprod = rExpansion(i) * box(iBox)%potentialExpansion(i)
            bprod = dble( rExpansion(aL1+m) * boxIrrExpansion(rList(iBox)%expStart+i) )
            cprod = rExpansion(i) * boxIrrExpansion(rList(iBox)%expStart+i)

            ! no need to divide the potential contribution by 2; doing only half-range of m does that
            !sump = sump + dble(rExpansion(i) * box(iBox)%potentialExpansion(i))
            sump = sump + dble(rExpansion(i) * boxIrrExpansion(rList(iBox)%expStart+i))
            ! due to m doing only half-range, force will have to be scaled by factor of 2 in the end
            sumf(1) = sumf(1) + dble( dcmplx(dL*xDx2y2,  dm*yDx2y2) * cprod ) - z * xDx2y2 * bprod
            sumf(2) = sumf(2) + dble( dcmplx(dL*yDx2y2, -dm*xDx2y2) * cprod ) - z * yDx2y2 * bprod
            sumf(3) = sumf(3) + bprod
          enddo

          ! m = L
          i = i + 1
          !cprod = rExpansion(i) * box(iBox)%potentialExpansion(i)
          cprod = rExpansion(i) * boxIrrExpansion(rList(iBox)%expStart+i)

          !sump = sump + dble(rExpansion(i) * box(iBox)%potentialExpansion(i))
          sump = sump + dble(rExpansion(i) * boxIrrExpansion(rList(iBox)%expStart+i))
          sumf(1) = sumf(1) + dble( dcmplx(dL*xDx2y2,  dL*yDx2y2) * cprod )
          sumf(2) = sumf(2) + dble( dcmplx(dL*yDx2y2, -dL*xDx2y2) * cprod )
          !sumf(3) contribution is zero for m = L
        enddo

        !energy = sump * bQ(iAtom)      ! sump is already corrected for double counting
        bPotential(iAtom) = bPotential(iAtom) + sump * 2.0d0       ! contribution from all particles
        bForce(1:3,iAtom) = bForce(1:3,iAtom) - sumf(1:3) * 2.0d0  ! compensate for half summation in index m
      enddo
      !fmmEnergy = fmmEnergy + energy   ! not going to return energy
#ifdef MPI0
    endif  ! iam == box(iBox)%proc
#endif
  enddo

end subroutine ElectricForceFF


subroutine DemoElectricForceFMM(oxyz, coor, q, Lexpansion, energy, potential, force, &
                                       nL, nSzTri, tri, tmpExpansion, Rexpansion)
  ! Compute potential and force due to far field in the position of particle (coor, q)
  ! Use irregular multipole expansion (Lexpansion) of remote charges accumulated in the box holding that particle
  ! All computed values are in Gauss units (charge in electrons, coordinates in Angstroms)
  ! Computation of energy includes correction against double counting
  !
  implicit none
  !
  double precision, intent(in)    :: oxyz(3)   ! center of the box in global coordinate system
  double precision, intent(in)    :: coor(3)   ! particle coordinates in global coordinate system
  double precision, intent(in)    :: q         ! particle charge
  double complex, intent(in)         :: Lexpansion(0:nSzTri) ! irregular expansion of remote charges
  double precision, intent(inout) :: energy    ! energy due to far field
  double precision, intent(inout) :: potential ! potential on the particle
  double precision, intent(inout) :: force(3)  ! force on the particle
  integer,          intent(in)    :: nL        ! expansion truncation
  integer,          intent(in)    :: nSzTri    ! triangular size of expansion array
  integer,          intent(in)    :: tri(0:nL) ! triangular pointers to expansion elements
  double complex, intent(inout)      :: tmpExpansion(0:nSzTri) ! scratch array
  double complex, intent(inout)      :: Rexpansion(0:nSzTri)   ! scratch array
  !
  integer          :: i, L, aL1, m
  double precision :: x,y,z, x2y2, xDx2y2, yDx2y2, LxDx2y2, LyDx2y2, zxDx2y2, zyDx2y2
  double precision :: dL, dm, sump, sumf(3), aprod, bprod
  double precision :: pot, potA, potB, r2, rij(3), forceA(3), forceB(3)
  double complex      :: cprod
  double precision, parameter :: zeroCoor = 1.0d-8, zero = 1.0d-14
  !
  ! determine particle coordinates in the local coordinate system of the box
  x = coor(1) - oxyz(1)
  y = coor(2) - oxyz(2)
  z = coor(3) - oxyz(3)

  x2y2 = x*x + y*y
 
  ! replace the particle coordinate by a very small increment to avoid division by zero
  if(x2y2 < zero) then
    x = dsign(zeroCoor,x)
    x2y2 = x*x + y*y
  endif

  xDx2y2 = x / x2y2
  yDx2y2 = y / x2y2
  zxDx2y2 = z * xDx2y2
  zyDx2y2 = z * yDx2y2

  ! compute chargless multipole expansion for radius-vector of the particle
  call DemoRegularExpansion(coor, 1.0d0, 1, oxyz, nL, nSzTri, tri, tmpExpansion, Rexpansion)

  ! L = 0 and m = 0
  ! sump is sum for potential; sumf is sum for force
  ! divide by 2 so that the entire sum over m can be multiplied by 2 in the end

  sump = dble(Rexpansion(0) * Lexpansion(0)) * 0.5d0 
  sumf = 0.0d0

  ! L >= 1
  do L = 1, nL
    dL  = dble(L)
    aL1 = tri(L-1)
    i   = tri(L)

    LxDx2y2 = dL * xDx2y2
    LyDx2y2 = dL * yDx2y2

    ! m = 0
    ! potential
    sump  = sump + dble(Rexpansion(i) * Lexpansion(i)) * 0.5d0   ! Re[ R(l,0) * L(l,0) ]

    ! force
    aprod = dble( Rexpansion(i)   * Lexpansion(i) ) * 0.5d0      ! Re[ R(l,0) * L(l,0) ]
    bprod = dble( Rexpansion(aL1) * Lexpansion(i) ) * 0.5d0      ! Re[ R(l-1,0) * L(l,0) ]
    sumf(1) = sumf(1) + LxDx2y2 * aprod - zxDx2y2 * bprod
    sumf(2) = sumf(2) + LyDx2y2 * aprod - zyDx2y2 * bprod
    sumf(3) = sumf(3) + bprod

    ! do half-range in index m (only positive values) and then multiply by 2 in the end
    ! 1 <= m <= L-1
    do m = 1, L-1
      i = i + 1

      ! potential
      sump = sump + dble(Rexpansion(i) * Lexpansion(i))  ! Re[ R(l,m) * L(l,m) ]

      ! force
      dm = dble(m)
      bprod = dble( Rexpansion(aL1+m) * Lexpansion(i) )  ! Re[ R(l-1,m) * L(l,m) ] 
      cprod = Rexpansion(i) * Lexpansion(i)              ! complex R(l,m) * L(l,m)
      sumf(1) = sumf(1) + dble( dcmplx(LxDx2y2,  dm*yDx2y2) * cprod ) - z * xDx2y2 * bprod
      sumf(2) = sumf(2) + dble( dcmplx(LyDx2y2, -dm*xDx2y2) * cprod ) - z * yDx2y2 * bprod
      sumf(3) = sumf(3) + bprod
    enddo

    ! m = L
    i = i + 1

    ! potential
    sump = sump + dble(Rexpansion(i) * Lexpansion(i))    ! Re[ R(l,l) * L(l,l) ]

    ! force
    cprod = Rexpansion(i) * Lexpansion(i)                ! complex R(l,l) * L(l,l)
    sumf(1) = sumf(1) + dble( dcmplx(LxDx2y2,  LyDx2y2) * cprod )
    sumf(2) = sumf(2) + dble( dcmplx(LyDx2y2, -LxDx2y2) * cprod )
    !sumf(3) contribution is zero for m = l
  enddo

  energy = energy + sump * q  ! skipping multiplication of sump by 2 corrects against double counting
  potential = potential + sump * 2.0d0         ! multiply potential and force terms by 2
  force(1:3) = force(1:3) - sumf(1:3) * 2.0d0  ! to compensate for half-summation in index m

end subroutine DemoElectricForceFMM


subroutine DemoElectricForceFF(iLevel, fmmEnergy, bPotential, bForce, tmpExpansion, rExpansion)
  !
  ! Use irregular FF expansion accumulated in the box
  ! Works with non-periodic and periodic FMM
  ! All computed values are in Gauss units (charge in electrons, coordinates in Angstroms)
  ! Computation of energy includes correction against double counting
  !
  use common_data, only: nL, tri, nSzTri, bC, bQ, vertexMin, nAtomsTotal, boxIrrExpansion, prtRank, leadRank
  use common_constants, only: zero, zeroCoor
  use fmm_plan, only: treeLevel, rList
  use timer
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  integer, intent(in)             :: iLevel      ! hierarchy level for the box carrying particles
  double precision, intent(inout) :: fmmEnergy   ! output; add FF energy to the initial value of variable
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri), rExpansion(0:nSzTri)        ! scratch arrays
  double precision, intent(inout) :: bPotential(nAtomsTotal), bForce(3,nAtomsTotal)  ! output
  !
  integer :: iam, nRanks, iAtom, iBox
  double precision :: oxyz(3)
  integer :: boxesStart, boxesEnd, xyz(3)
  !
#ifdef MPI
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
#else
  iam = 0
  nRanks = 1
#endif

#ifdef TIMER
  call TimerInit(2)
#endif
  !
  ! compute far-field potential and force in the position of particles
  ! obtain force by differentiation of the above computed FF potential (in multipole expansion)
  ! compute FF contribution to the energy of the system
  !
  xyz = 0  ! dummy initialization to stop the compiler from issuing warning about uninitialized variable
  stop "Error in DemoElectricForceFF"
  !
  boxesStart = (8**(iLevel+1) - 1) / 7
  boxesEnd   = boxesStart + 8**iLevel - 1
  do iBox = boxesStart, boxesEnd
  !do iBox = treeLevel(iLevel)%bStart, treeLevel(iLevel)%bEnd
    ! find center of the box
    !oxyz(1) = vertexMin(1) + (dble(box(iBox)%xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
    !oxyz(2) = vertexMin(2) + (dble(box(iBox)%xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
    !oxyz(3) = vertexMin(3) + (dble(box(iBox)%xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge
    oxyz(1) = vertexMin(1) + (dble(xyz(1)) + 0.5d0) * treeLevel(iLevel)%edge
    oxyz(2) = vertexMin(2) + (dble(xyz(2)) + 0.5d0) * treeLevel(iLevel)%edge
    oxyz(3) = vertexMin(3) + (dble(xyz(3)) + 0.5d0) * treeLevel(iLevel)%edge

    do iAtom = rList(iBox)%atomStart, rList(iBox)%atomEnd
      ! add up potential and force on the particle, and energy contribution
      !call DemoElectricForceFMM(oxyz, bC(1,iAtom), bQ(iAtom), box(iBox)%potentialExpansion, &
      !     fmmEnergy, bPotential(iAtom), bForce(1,iAtom), nL, nSzTri, tri, tmpExpansion, rExpansion)
      call DemoElectricForceFMM(oxyz, bC(1,iAtom), bQ(iAtom), boxIrrExpansion(rList(iBox)%expStart), &
           fmmEnergy, bPotential(iAtom), bForce(1,iAtom), nL, nSzTri, tri, tmpExpansion, rExpansion)
    enddo
  enddo

#ifdef TIMER
  if(iam == prtRank) then
    write(*,'(a,f9.1,a)') "Timer: FMM far-field potential and force: ", TimerCheck(2), "s"
  endif
#endif

end subroutine DemoElectricForceFF


subroutine print_complex(a,n,str)
  use common_data, only: prtRank, leadRank
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  implicit none
  integer, intent(in)           :: n
  double complex, intent(in)    :: a(n)
  character (len=6), intent(in) :: str
  integer :: i, iam, nRanks
#ifdef MPI
  integer :: mpiError, status(MPI_STATUS_SIZE)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
#else
  iam = 0
  nRanks = 1
#endif
  if(iam == prtRank) then
    do i = 1, n
      write(*,'(a10,4x,i8,2x,2f16.8)') str, i, dble(a(i)), aimag(a(i))
    enddo
  endif
end subroutine print_complex


subroutine print_double(a,b,n,str)
  use common_data, only: prtRank, leadRank
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  implicit none
  integer, intent(in)           :: n
  double precision, intent(in)  :: a(n), b(n)
  character (len=6), intent(in) :: str
  integer :: i, iam, nRanks
#ifdef MPI
  integer :: mpiError, status(MPI_STATUS_SIZE)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
#else
  iam = prtRank
  nRanks = 1
#endif
  if(iam == prtRank) then
    do i = 1, n
      write(*,'(a10,4x,i8,2x,2f16.8)') str, i, a(i), b(i)
    enddo
  endif
end subroutine print_double


