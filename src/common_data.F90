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

module common_data
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  implicit none
  save

  integer (kind=8)              :: nAtoms       ! number of particles in the unit cell
  integer (kind=8)              :: nAtomsTotal  ! number of particles in the system
  integer (kind=8)              :: nAtomsRank   ! number of particles in the distributed array
  integer                       :: nRepLayers   ! number of replication layers
  integer                       :: nRepLayers1  ! nReplayers + 1
  integer                       :: nRepLayers2  ! nReplayers1 * nReplayers1
  integer (kind=8)              :: nReplicas    ! the total number of replicas of unit cell
  integer (kind=8)              :: nBoxesTotal  ! number of boxes in the system
  integer                       :: nBoxesSList  ! total number of boxes in the S-list
  integer                       :: maxBoxesSlist ! max number of boxes in the S-list across all ranks
  integer (kind=8)              :: bBoxesStart  ! starting index of B-boxes at the lowest level
  integer                       :: sBoxLLstart  ! local (for S-list) starting index of S-boxes at loLevel
  integer (kind=8)              :: rBoxBottomStart  ! starting global index of R-boxes at loLevel
  integer (kind=8)              :: rBoxBottomEnd    ! ending global index R-boxes at loLevel
  integer (kind=8)              :: nLLBoxes_rmRank  ! number of LL boxes assigned to rmRank in RMG
  integer (kind=8), allocatable :: rankAtoms(:) ! number of atoms per MPI rank
  integer (kind=8), allocatable :: dispAtoms(:) ! displacement, starting atom position in the rank
  integer (kind=8), allocatable :: rankBoxes(:) ! number of lowest-level boxes per MPI rank
  integer (kind=8), allocatable :: dispBoxes(:) ! displacement, starting position for box multipoles
  integer (kind=8), allocatable :: fBoxes(:)    ! number of floating boxes from ffLevel to ddLevel-1 per rank 
  integer                       :: nFloatingBoxes   ! number of floating boxes involved in far-field interactions
  integer (kind=4), allocatable :: nExtraMultipoleBoxes(:)  ! extra multipole (E) boxes per rank in RMG
  integer                       :: nDevicesNode ! number of devices per node
  integer                       :: nThreads     ! number of OpenMP threads per MPI rank
  integer                       :: ws           ! separation parameter
  real (kind=8), allocatable    :: cell(:,:)    ! input unit cell particle Cartesian coordinates in original order
  real (kind=8), allocatable    :: w(:)         ! input unit cell particle point charges (weight) in original order
  double precision, allocatable :: aC(:,:)      ! all particle Cartesian coordinates sorted to FMM boxes
  double precision, allocatable :: aQ(:)        ! all particle point charges sorted to FMM boxes
  double precision, allocatable :: bC(:,:)      ! rank-assigned particle Cartesian coordinates sorted to FMM boxes
  double precision, allocatable :: bQ(:)        ! rank-assigned particle point charges sorted to FMM boxes
  double precision, allocatable :: cC(:,:)      ! rmRank-assigned particle Cartesian coordinates sorted to FMM boxes
  double precision, allocatable :: cQ(:)        ! rmRank-assigned particle point charges sorted to FMM boxes
  double precision              :: vertexMin(3), vertexMax(3)
  double precision              :: cubeEdge, edgeReciprocal
  double precision              :: unitCellEdge          ! PBC box size
  double precision              :: primitiveCellEdge     ! input box size
  integer (kind=4)              :: intCubeEdge, intPrimitiveEdge  ! integer form of cubeEdge and unitCellEdge
  logical                       :: scaleToUnitSize  ! yes/no to scale down the PBC box to unit size
  integer                       :: nL           ! multipole order
  integer                       :: nSzTri       ! number of elements in lower triangle for l=nL
  integer                       :: nSzWset      ! number of elements in Wigner set for l=nL
  integer, allocatable          :: tri(:)       ! triangular pointers to multipole terms with m=0
  double complex, allocatable   :: rExpan(:)    ! regular hamonics multipole expansion
  double complex, allocatable   :: sExpan(:)    ! irregular harmonics multipole expansion (0:nSzTri)
  double complex, target, allocatable :: boxRegExpansion(:) ! regular expansion (0:(nSzTri+1)*(nBoxesRList+nFFBoxes)-1)
  double complex, allocatable   :: boxIrrExpansion(:) ! irregular expansion (0:numBoxesSList*(nSzTri+1))
  double precision, allocatable :: d(:)         ! Wigner set of d-matrices (1:nSzWset)
  integer                       :: M2Lene, M2Lrot   ! rotation-mode flags
  integer                       :: ncoeff       ! number of elements in arrays coefflkm()
  double precision, allocatable :: coefflkm(:)  ! sq[(l-k)!(l+k)! / (l-m)!(l+m)!]
  double precision, allocatable :: zonal(:)     ! zonal (real) harmonics
 !integer (kind=8), allocatable :: nearNeighbors(:), farNeighbors(:)
 !integer (kind=4), allocatable :: nearNeighbors(:), farNeighbors(:)
  integer (kind=8), allocatable :: nearNeighbors(:), farNeighbors(:)
  integer (kind=8)              :: nNearNeighbors, nFarNeighbors
  integer                       :: lowestLevel  ! index of the lowest hierarchy level
  integer                       :: ddLevel_input ! the value of ddLevel provided in the command-line input
  integer                       :: maxNNperBox, maxFNperBox   ! maximum number of near and far neighbors per box
#ifdef MPI
#ifdef CPLR_NVFTN
  integer                       :: MPI_COMM_RMGROUP    ! MPI communicator local for each sub-communicator
  integer                       :: MPI_COMM_NODE       ! MPI communicator for a single node
  integer                       :: MPI_COMM_PNODE      ! MPI communicator for a single pass in the FMM Plan
  integer                       :: MPI_COMM_NZERO      ! MPI communicator for ranks 0 on each node
  integer                       :: MPI_COMM_PZERO      ! MPI communicator for ranks 0 in each pass in FMM Plan
  integer                       :: MPI_COMM_PGROUP     ! MPI communicator for the ranks in Plan pass
#else
  type(MPI_Comm)                :: MPI_COMM_RMGROUP    ! MPI communicator local for each sub-communicator
  type(MPI_Comm)                :: MPI_COMM_NODE       ! MPI communicator for a single node
  type(MPI_Comm)                :: MPI_COMM_PNODE      ! MPI communicator for a single pass in the FMM Plan
  type(MPI_Comm)                :: MPI_COMM_NZERO      ! MPI communicator for ranks 0 on each node
  type(MPI_Comm)                :: MPI_COMM_PZERO      ! MPI communicator for ranks 0 in each pass in FMM Plan
  type(MPI_Comm)                :: MPI_COMM_PGROUP     ! MPI communicator for the ranks in Plan pass
#endif
#endif
  integer                       :: nNodes              ! the number of allocated nodes
  integer                       :: nodeSize            ! the number of ranks per node
  integer, allocatable          :: segmentRank(:)      ! MPI rank id for the segment
  integer, allocatable          :: segmentNBoxes(:)    ! number of boxes in the segment
  integer, allocatable          :: segmentGaps(:)      ! number of boxes in the segment
  integer (kind=8)              :: nBoxSegments        ! number of segments
  integer                       :: nBoxesRList         ! total number of boxes in R-list
  integer                       :: nRemoteFFBoxes      ! number of unique remote FF boxes needed for the rank
  integer (kind=8)              :: nParticlesRank      ! number of particles asigned to the rank iam (system-wide)
  integer (kind=8)              :: nParticlesRMRank    ! number of particles asigned to the rank rmRank (group-based)
  logical                       :: lComputeDirectSum   ! flag to skip or compute the direct sum
  integer, allocatable          :: particleRanges(:,:) ! particle start and end indices for the box
  integer (kind=4)              :: prtRank             ! MPI rank assigned to print function
  integer (kind=4)              :: leadRank            ! lead MPI rank
  integer (kind=4)              :: nodeRank            ! MPI rank of a single-node sub-communicator
  integer (kind=8)              :: nBoxesBefore_loLevel ! number of boxes before loLevel in the global index
  integer (kind=4), allocatable :: sBoxesTopNum(:)     ! number of boxes at the top of S-tree at ddLevel
  integer (kind=8), allocatable :: sBoxesTopDispl(:)   ! diplacemnt in the number of boxes at the top of S-tree
  integer (kind=8)              :: iAtomsBegin, iAtomsEnd ! starting and ending particle indices on the rank

  ! Construction of R-tree begins from distributing the boxes at rmLevel among the ranks in RMG sub-communicator
  integer                       :: minTopRtreeSize        ! minimum number of boxes at rmLevel per rank in RMG
  integer                       :: minTopRtreeRanks       ! the number of ranks that get the minimal number of boxes
  integer                       :: minTopRtreeGroup       ! the total number of boxes in the minTopRtreeRanks group

  integer (kind=8)              :: expShift               ! starting address of regular expansions for node-based DD
  integer                       :: numBoxesSNRF           ! number of S N R F boxes
  integer                       :: szSubComm              ! size of sub-communicator in terms of number of nodes
  integer                       :: rmSize                 ! number of ranks in sub-communicator
  integer                       :: nRanks                 ! total number of ranks
  integer                       :: nSplitRanks            ! number of ranks per node used in a single pass in the FMM Plan

  logical (kind=1)              :: lxyzq                  ! indicate if the input file is in xyzq format 

  ! Multiple input datasets
  logical (kind=1)              :: lMultipleDataSets      ! use multiple datasets
  integer (kind=4)              :: nDataSets              ! number of datasets
  integer (kind=4)              :: nGroups                ! number of box groups
  integer (kind=4)              :: nSelections            ! number of box selections`
  integer (kind=8), allocatable :: nAtomsGroupSelected(:) ! number of atoms in the group
  integer (kind=8), allocatable :: atomRange(:,:)         ! starting and ending atom indices for each dataset
  integer (kind=2), allocatable :: datasetid(:)           ! dataset id for each replica tuple, ijk
  integer (kind=2), allocatable :: groupid(:)             ! group id for each replica tuple, ijk
  character (len=512), allocatable :: dataname(:)         ! dataset filenames
  integer (kind=8)              :: nSelectedLLboxesRank   ! number of selected LL-boxes in S-list in the rank

  ! test parameters
  double precision              :: stdEnergy
  double precision              :: memUsed, memMax, memMaxMyNode, memMaxAllNodes
  logical (kind=1)              :: debugEachMPIrank

end module common_data
