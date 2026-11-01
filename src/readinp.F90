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

! convert the input file to xyzq format
#define OUTPUT_XYZQ_BINFORMAT
#undef  OUTPUT_XYZQ_BINFORMAT

! generate xyz file for visualization
#define OUTPUT_VIS
#undef  OUTPUT_VIS

subroutine readinp()

#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif

  use common_data, only: nAtomsTotal, w, cell, nAtoms, nAtomsRank, unitCellEdge, scaleToUnitSize, &
                         cubeEdge, nL, ws, memUsed, memMax, lowestLevel, vertexMin, vertexMax, &
                         nRepLayers, nRepLayers1, nRepLayers2, nReplicas, primitiveCellEdge, &
#ifdef MPI
                         MPI_COMM_NODE, MPI_COMM_NZERO, &
#endif
                         prtRank, leadRank, szSubComm, nNodes, nodeSize, nRanks, nSplitRanks, lxyzq, &
                         debugEachMPIrank, intCubeEdge, intPrimitiveEdge, &
                         lMultipleDataSets, nDataSets, atomRange, datasetid, dataname, &
                         nGroups, groupid, nSelections, nAtomsGroupSelected, ddLevel_input
  use common_constants, only: machineZero, TYPE_FILE, TYPE_INT, TYPE_DBLE
  use timer
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_f_pointer, c_associated
  use allocator
  implicit none

  character (len=512) :: fname_psf, fname_pdb, fname_idx, fname_xyzq
  integer             :: nArgs, iam, nodeZeroRank, nodeMinMem
  integer, external   :: iargc
  character (len=80)  :: str
  integer, parameter  :: inpUnit = 10, idxUnit = 11
  integer (kind=4)    :: i, j, k, n, io, allocStat, iArg, idx, dsid, grpid
  integer (kind=8)    :: iAtom, ijk
 !integer (kind=8)    :: atomStart
  real (kind=8)       :: qTotal
  real (kind=8)       :: atomCharge, atomMass
 !double precision    :: reciprocal, diff(3)
  ! debug options to replicate the unit cell
  double precision    :: memLocal
  double precision    :: translation(3)
  ! argument place holders
  character (len=512) :: tmp_file
  integer (kind=4)    :: tmp_int, nSkip, nWrite
  double precision    :: tmp_dble
  integer (kind=8)    :: tmp_int8
  !character (len=2),allocatable :: atName(:) ! atom names
  character (len=2)   :: atName
  integer, parameter  :: pos_prtRank   = 8    ! position of prtRank in the program arguments list
  integer, parameter  :: pos_szSubComm = 9    ! position of szSubComm in the program arguments list
  integer, parameter  :: pos_ddLevel  = 11    ! position of ddLevel in the program arguments list
  character (len=512) :: out_string
  logical (kind=1)    :: psfExtendedFormat    ! distinguish between regular and extended formats of PSF file
  integer (kind=8)    :: prtRank_8, nRanks_8  ! int*8 placeholders for prtRank and nRanks
  integer (kind=4)    :: commonMultiplier
  integer (kind=4)    :: memFreeNode, memFreeNodeMin
  integer (kind=8)    :: nReplicasIdx         ! the number of replicas specified in the IDX file

#ifdef MPI
  integer :: mpiError, status(MPI_STATUS_SIZE), resultlen
  character (len=MPI_MAX_PROCESSOR_NAME) :: hostname
  integer (kind=4), dimension(0:1) :: nodeData, nodeDataMin

  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
  call MPI_Get_processor_name(hostname, resultlen, mpiError)
#else
  iam = 0
  nRanks = 1
  nodeMinMem = 0
  nodeZeroRank = 0
#endif

#ifdef TIMER
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
#endif
  call TimerInit(2)
#endif

  ! set the default values
  prtRank = 0      ! use rank 0 for stdout
  leadRank = 0     ! use rank 0 as the lead
  szSubComm = 1    ! use single-node subcommunicator
  lxyzq = .false.  ! .false. indicates that the input file format is in PSF & PDB format
  debugEachMPIrank  = .false.
  lMultipleDataSets = .false.
  nGroups = 0      ! no box groups
  nSelections = 0  ! no box selections
  ddLevel_input = 0 ! use the default value of ddLevel

  ! number of arguments in the command-line options
  nArgs = command_argument_count()

  ! each MPI rank reads the value of prtRank from program arguments
  if(nArgs >= pos_prtRank) then
    call getarg(pos_prtRank, out_string)

    ! process the case when the value of prtRank provided in the command line is too big (type int*8)
    read(out_string,*) prtRank_8

    ! negative value of prtRank initiates the debug mode
    debugEachMPIrank = (prtRank_8 < 0_8)

    ! MPI does not support int*8 for the rank id and the number of ranks
    nRanks_8 = nRanks  ! convert int4 value to int8

    ! test for  0 <= prtRank <= nRanks-1  and type cast from int8 to int4
    prtRank = merge(abs(prtRank_8), nRanks_8 - 1_8, abs(prtRank_8) < nRanks_8)

    ! assign the lead rank
    if(debugEachMPIrank) then
      ! debug run
      leadRank = 0   ! use rank 0 as the lead
      prtRank = iam  ! make each MPI rank print to stdout
    else
      ! normal run
      leadRank = prtRank  ! set the lead rank equal to the print rank
    endif
  endif

  ! each MPI rank reads the value of szSubComm from program arguments
  if(nArgs >= pos_szSubComm) then
    call getarg(pos_szSubComm, out_string)
    read(out_string,*) szSubComm
    szSubComm = min0(szSubComm, nNodes)  ! do not exceed the total number of nodes
    if(szSubComm <= 0) szSubComm = 1
    if(mod(nNodes, szSubComm) > 0 .and. iam == prtRank) then
      write(*,'(/a,i0,a,i0)') "Error: Sub-communicator size ", szSubComm, &
                              " is not a divisor of the total number of nodes ", nNodes
      call ExitPrg(11)
    endif
  endif

  ! each MPI rank reads the value of ddLevel from program arguments
  if(nArgs >= pos_ddLevel) then
    call getarg(pos_ddLevel, out_string)
    if(len_trim(out_string) > 2 .and. iam == prtRank) then
      write(*,'(/a,a,a)') "Error: ddLevel_input  ", trim(out_string), "  is invalid"
      call ExitPrg(11)
    endif
    read(out_string,*) ddLevel_input
    ddLevel_input = min0(ddLevel_input, 11)  ! do not exceed the max allowed value
    if(ddLevel_input < 0) ddLevel_input = 0
  endif

  ! check the available host memory on the nodes
  memFreeNode = stat_allmemory() / 1024_8 / 1024_8 ! reported in bytes and converted to Mb
  memFreeNodeMin = memFreeNode
  nodeMinMem = 0

#ifdef MPI
  if(nNodes > 1 .and. MPI_COMM_NZERO /= MPI_COMM_NULL) then
    call MPI_Comm_rank(MPI_COMM_NZERO, nodeZeroRank, mpiError)
    nodeData(0) = memFreeNode
    nodeData(1) = nodeZeroRank
    if(debugEachMPIrank) then
      call MPI_AllReduce(nodeData, nodeDataMin, 1, MPI_2INTEGER, MPI_MINLOC, MPI_COMM_NZERO, mpiError)
    else
      call MPI_Reduce(nodeData, nodeDataMin, 1, MPI_2INTEGER, MPI_MINLOC, 0, MPI_COMM_NZERO, mpiError)
    endif
    memFreeNodeMin = nodeDataMin(0)
    nodeMinMem     = nodeDataMin(1)
  endif
#endif

  if(iam == prtRank .or. debugEachMPIrank) then
    write(*,'(a,i0,a,a)') "MPI rank ", iam, " runs on the host: ", trim(hostname)
    write(*,'(a,i7,a,i7,a)') &
      "Beginning Free Host Memory in the nodes    Local:", memFreeNode, " Mb    Min:", memFreeNodeMin, " Mb"
    if(MPI_COMM_NZERO /= MPI_COMM_NULL .and. nodeZeroRank == 0) &
      write(*,'(a,i0,a,i0)') "The min memory of ", memFreeNodeMin, " MB is in the node: ", nodeMinMem
  endif

  ! Each MPI rank reads the program arguments

  iArg = 1  ! start reading from first argument

  ! get FMM cube edge
  call ReadNextArg(iArg, nArgs, TYPE_DBLE, tmp_file, tmp_int, cubeEdge)
  if(cubeEdge < machineZero) cubeEdge = 10.0d0

  ! get polynomial order
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, nL, tmp_dble)
  if(nL == 0) nL = 20

  ! get well-separation parameter
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, ws, tmp_dble)
  if(ws == 0) ws = 2

  ! get input files
  call ReadNextArg(iArg, nArgs, TYPE_FILE, fname_psf, tmp_int, tmp_dble)
  call ReadNextArg(iArg, nArgs, TYPE_FILE, fname_pdb, tmp_int, tmp_dble)
  fname_xyzq = fname_psf

  ! get unit-cell scaling
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, i, tmp_dble)
  scaleToUnitSize = i .eq. 1

  ! get cell replication parameters (used for code testing)
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, nRepLayers, tmp_dble)
  if(nRepLayers < 0) then
    nRepLayers = abs(nRepLayers)  ! the multiple datasets mode uses this value of nReplayers
    lMultipleDataSets = .true.
  endif
  nRepLayers1 = nRepLayers + 1
  nRepLayers2 = nRepLayers1 * nRepLayers1
  nReplicas = (nRepLayers + 1_8) ** 3_8  ! requires int*8 for nRepLayers > 1023

  ! get prtRank
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, tmp_int, tmp_dble)

  ! get szCommSize
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, tmp_int, tmp_dble)

  ! get nSplitRanks
  i = 0
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, i, tmp_dble)
  nSplitRanks = i

  ! get ddLevel_input
  call ReadNextArg(iArg, nArgs, TYPE_INT, tmp_file, tmp_int, tmp_dble)
 
  if(iam == prtRank) then
    write(*,'(a,18x,a)') "Input PSF: ", trim(fname_psf)
    write(*,'(a,18x,a)') "Input PDB: ", trim(fname_pdb)
    write(*,'(a,9x,i0)') "Replication layers: ", nRepLayers
    write(*,'(a,9x,l1)') "Multiple datasets : ", lMultipleDataSets
  endif

  nDataSets = 1  ! single input dataset is the default
  if(lMultipleDataSets) then
    if(iam == leadRank) then
      ! only the lead rank reads the IDX file
      ! compose IDX file name
      write(fname_idx,'(a,a)') fname_pdb(1 : max0(0, index(fname_pdb,".") - 1) ), ".idx"
      ! open IDX file
      open (unit=idxUnit, file=fname_idx, status='old', form='formatted')
      ! read the number of datasets
      read(idxUnit,*) nDataSets, nReplicasIdx, nGroups
      write(*,'(a,9x,i0)') "nReplicasIdx:       ", nReplicasIdx
    endif
#ifdef MPI
    ! broadcast the number of datasets
    if(nRanks > 1) then
      call MPI_Bcast(nDataSets, 1, MPI_INTEGER4, leadRank, MPI_COMM_WORLD, mpiError)
      if(mpiError /= MPI_SUCCESS) call exitPrg(6)
      call MPI_Bcast(nGroups,   1, MPI_INTEGER4, leadRank, MPI_COMM_WORLD, mpiError)
      if(mpiError /= MPI_SUCCESS) call exitPrg(6)
    endif
#endif
  endif

  ! allocate storage for datasets
  if(iam == leadRank) then
    ! only the lead rank, which reads input data, allocates dataname
    ! each use of dataname must be protected with a test for leadRank
    allocate (dataname(nDataSets), atomRange(2,nDataSets), datasetid(0:nReplicas-1), &
              groupid(0:nReplicas-1), nAtomsGroupSelected(0:nGroups), stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)
    memLocal = dble(nDataSets) * dble(512+8*2) + dble(nReplicas * 2*2)
    memMax = dmax1(memMax, (memUsed + memLocal))
    !
    dataname(1) = fname_pdb ! assume a single dataset
  else
    allocate (atomRange(2,nDataSets), datasetid(0:nReplicas-1), &
              groupid(0:nReplicas-1), nAtomsGroupSelected(0:nGroups), stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)
    memLocal = dble(nDataSets) * dble(8*2) + dble(nReplicas * 2*2)
    memMax = dmax1(memMax, (memUsed + memLocal))
  endif

  ! set the default values for datasetid and groupid
  datasetid = 1  ! all replicas are presumed to be created from the same default dataset 1
  groupid   = 0  ! no selected groups

  ! number of atoms in each group
  nAtomsGroupSelected = 0

  ! assign a dataset to each replica
  if(lMultipleDataSets) then
    if(iam == leadRank) then
      ! read the dataset names
      do i = 1, nDataSets
        read(idxUnit,'(a80)') dataname(i)  ! file name of the input data set, i
      enddo
      !
      ! read the dataset id for each replica i,j,k (index to the left changes faster)
      do ijk = 1_8, nReplicasIdx
        read(idxUnit,*) dsid, grpid, i, j, k
        idx = i + j * nRepLayers1 + k * nRepLayers2
        datasetid(idx) = dsid
        groupid(idx)   = grpid
      enddo
      close(idxUnit)
    endif
#ifdef MPI
    ! broadcast the datasetid
    if(nRanks > 1) then
      call MPI_Bcast(datasetid, nReplicas, MPI_INTEGER2, leadRank, MPI_COMM_WORLD, mpiError)
      call MPI_Bcast(groupid, nReplicas, MPI_INTEGER2, leadRank, MPI_COMM_WORLD, mpiError)
      if(mpiError /= MPI_SUCCESS) call exitPrg(6)
    endif
#endif
  else
    ! single dataset
    datasetid = 1  ! all replicas refer to the dataset 1 in this case
    groupid   = 0  ! no groups
  endif

  nSelections = count(groupid > 0)  ! number of box selections

  if(nGroups > nSelections) then
    write(*,'(a,2i10)') "Error: nnGroups > nSelections ", nGroups, nSelections
    call exitPrg(11)
  endif

  ! sum up the number of atoms in all datasets to store coordinates of all datasets in one array
  nAtoms = 0

  if(iam == leadRank) then
    if(lMultipleDataSets) then
      !
      ! Multiple input datasets
      !
      ! read the number of atoms in each dataset
      ! requires XYZQ data format
      do i = 1, nDataSets
        ! open XYZQ file
        open(unit=idxUnit, file=trim(dataname(i)), form='unformatted', access='stream', status='old')

        ! read the number of atoms
        read(idxUnit) atomRange(2,i)  ! number of atoms in the dataset, i
        close(idxUnit)

        ! determine atom ranges in the particle buffer, cell(1:nAtoms)
        atomRange(1,i) = nAtoms + 1               ! starting atom number for dataset i
        atomRange(2,i) = nAtoms + atomRange(2,i)  ! ending atom number for dataset i
        nAtoms = atomRange(2,i)                   ! sum up the number of atoms
      enddo
    else
      !
      ! Single input data set
      !
      if(lxyzq) then
        ! the input file is in XYZQ format
        ! open XYZQ file
        open(unit=inpUnit, file=trim(fname_xyzq), form='unformatted', access='stream', status='old')

        ! get the number of atoms
        nAtoms = 0
        read(inpUnit) nAtoms

      else
        ! the input file is in PSF & PDB format
        ! open PSF file
        open (unit=inpUnit, file=fname_psf, status='old', form='formatted')

        ! determine the PSF format (regular or extended) by reading the first line
        read(inpUnit,'(a80)',iostat=io) str
        if(io == 0) then
          if(index(str,"EXT") > 0) then
            psfExtendedFormat = .true.
          else
            psfExtendedFormat = .false.
          endif
        else
          write(*,'(a)') "Cannot read first line in PSF file"
          close(inpUnit)
          call exitPrg(0)
        endif

        ! get the number of atoms
        nAtoms = 0
        do while (nAtoms == 0)
          read(inpUnit,'(a80)',iostat=io) str
          if(io == 0) then
            if(index(str,"!NATOM") > 0) then
              read(str,*) nAtoms
              exit
            endif
          else
            write(*,'(a)') "Cannot find string !NATOM in PSF file"
            close(inpUnit)
            call exitPrg(0)
          endif
        enddo
      endif ! lxyz

      atomRange(1,1) = 1       ! starting atom number for dataset 1
      atomRange(2,1) = nAtoms  ! ending atom number for dataset 1
    endif ! test for lMultipleDataSets
  endif ! test for iam == leadRank

#ifdef MPI
  if(nRanks > 1) then
    ! broadcast the starting and ending atomic indices for each dataset
    call MPI_Bcast(atomRange, 2*nDataSets, MPI_INTEGER8, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(6)
  endif
#endif
  nAtoms = atomRange(2,nDataSets)

  ! count the total number of atoms including replicated cells
  nAtomsTotal = 0_8
  do i = 1, nDataSets
    ! number of atoms multiplied by the number of replicas in the particular dataset
    nAtomsTotal = nAtomsTotal + (atomRange(2,i) - atomRange(1,i) + 1_8) * count(datasetid == i)
  enddo
  nAtomsRank = nAtomsTotal / nRanks

  if(iam == prtRank) then
    write(*,'(a,9x,i0)') "Number of datasets :", nDataSets
    do i = 1, nDataSets
      if(iam == leadRank) then
        write(*,'(/a,i6,3x, a)')"Dataset id name    :", i, trim(dataname(i))  ! leadRank
      else
        write(*,'(a,i6,3x, a)') "Dataset id name    :", i, trim(fname_pdb)    ! not leadRanks
      endif
      write(*,'(a,i6,3x,i0)')   "Dataset id nAtoms  :", i, atomRange(2,i) - atomRange(1,i) + 1_8
      write(*,'(a,i6,2(3x,i0))')"Dataset atom ranges:", i, atomRange(1,i), atomRange(2,i)
      write(*,'(a,i6,3x,i0)')   "Dataset id replicas:", i, count(datasetid == i)
    enddo

    write(*,'(/a,9x,i0)')"Total replicas     :", nReplicas
    write(*,'(a,9x,i0)') "Particles per rank :", nAtomsRank
    write(*,'(a,9x,i0)') "Number of particles:", nAtomsTotal
    write(*,'(a,9x,i0)') "Number of boxgroups:", nGroups
    write(*,'(a,9x,i0)') "Box selections     :", nSelections

!    if(nSelections > 0 .and. nGroups > 0) then 
!      write(*,'(/a)') "Selections:"
!      do n = 1, min0(nSelections,8)
!        do ijk = 0, nReplicas-1
!          if(groupid(ijk) == n) then
!            k = ijk / nRepLayers2
!            j = (ijk - k*nRepLayers2) / nRepLayers1
!            i = mod(ijk, nRepLayers1)
!            write(*,'(a,i5,2(4x,a,i0),4x,a,3(2x,i0))') &
!            "Group: ", n, "dataset: ", datasetid(ijk), "replica: ", ijk+1, " i,j,k ",i,j,k
!          endif
!        enddo
!      enddo
!    endif

  endif   ! iam

  ! compute the number of atoms in each group; need it for sanity check

  nAtomsGroupSelected(0) = nAtomsTotal  ! group 0 includes all atoms in the saystem by default
  if(nSelections > 0) then 
    do ijk = 0, nReplicas-1
      if(groupid(ijk) > 0) then
        k = ijk / nRepLayers2
        j = (ijk - k*nRepLayers2) / nRepLayers1
        i = mod(ijk, nRepLayers1)
        !
        ! update the number of particles in groupid(ijk)
        nAtomsGroupSelected(groupid(ijk)) = nAtomsGroupSelected(groupid(ijk)) + &
               (atomRange(2,datasetid(ijk)) - atomRange(1,datasetid(ijk)) + 1_8)
        !
        ! update the number of particles in group 0
        nAtomsGroupSelected(0) = nAtomsGroupSelected(0) - &
               (atomRange(2,datasetid(ijk)) - atomRange(1,datasetid(ijk)) + 1_8)
      endif
    enddo
  endif

  ! allocate arrays
  allocate (w(nAtoms), cell(3,nAtoms), stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
  memLocal = memLocal + dble(nAtoms*4*8)
  memMax = dmax1(memMax, (memUsed + memLocal))

  ! read atomic charges
  qTotal = 0.0d0
  if(iam == leadRank) then
    if(lMultipleDataSets) then
      do i = 1, nDataSets
        ! the input file is in XYZQ format
        ! read atomic coordinates and charges
        open(unit=inpUnit, file=trim(dataname(i)), form='unformatted', access='stream', status='old')
        read(inpUnit) tmp_int8  ! nAtoms
        do iAtom = atomRange(1,i), atomRange(2,i)
          read(inpUnit) cell(1,iAtom), cell(2,iAtom), cell(3,iAtom), w(iAtom)
          qTotal = qTotal + w(iAtom)
        enddo
        close(inpUnit)
      enddo
    else ! single dataset
      if(lxyzq) then
        ! the input file is in XYZQ format
        ! read atomic coordinates and charges
        do iAtom = 1, nAtoms
          read(inpUnit) cell(1,iAtom), cell(2,iAtom), cell(3,iAtom), w(iAtom)
          qTotal = qTotal + w(iAtom)
        enddo
      else
        ! the input file is in PSF & PDB format
        ! read atomic charges
        qTotal = 0.0d0
        if(psfExtendedFormat) then
          do iAtom = 1, nAtoms
            read(inpUnit,'(53x,d11.6,3x,d11.4)') atomCharge, atomMass
            w(iAtom) = atomCharge
            qTotal = qTotal + w(iAtom)
          enddo
        else
          do iAtom = 1, nAtoms
            read(inpUnit,'(33x,d11.6,3x,d11.4)') atomCharge, atomMass
            w(iAtom) = atomCharge
            qTotal = qTotal + atomCharge
          enddo
        endif
      endif  ! not lxyzq
      close(inpUnit)
    endif  ! not lMultipleDataSets

    !if(abs(qTotal) > machineZero) then
    !  ! useful for development purposes when reading a portion of the system charges
    !  write(*,'(a)') "Compensating charge has been aded to the last particle"
    !  w(nAtoms) = w(nAtoms) - qTotal
    !  qTotal = sum(w, 1)
    !endif

    write(*,'(/a,3x,f24.10)') "Total charge: ", qTotal
    !if(abs(qTotal) > machineZero) call exitPrg(5)  ! Non-zero charge
  endif  ! iam

#ifdef MPI
  if(nRanks > 1) then
    ! broadcast the charges
    call MPI_Bcast(w, nAtoms, MPI_DOUBLE_PRECISION, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(6)
  endif
#endif

  ! read atomic coordinates
  if(iam == leadRank .and. .not. lxyzq .and. .not. lMultipleDataSets) then
    ! single dataset
    ! the input file is in PSF & PDB format
    ! open PDB file
    open (unit=inpUnit, file=fname_pdb, status='old', form='formatted')
 
    iAtom = 0
    j = 0
    do while (iAtom < nAtoms)
      j = j + 1
      read(inpUnit,'(a80)',iostat=io) str
      if(io == 0) then
        if(index(str,"ATOM") == 1 .or. index(str,"HETATM") == 1) then
          iAtom = iAtom + 1
          read(str,'(30x,3f8.3,22x,a2)') cell(1,iAtom), cell(2,iAtom), cell(3,iAtom), atName
        endif
      else
        write(*,'(a,i0)') "Failed reading atomic coordinates; line: ", j
        write(*,*) str
        close(inpUnit)
        call exitPrg(0)
      endif
    enddo
    close(inpUnit)
  endif  ! iam

#ifdef MPI
  if(nRanks > 1) then
    ! broadcast the coordinates
    call MPI_Bcast(cell, 3*nAtoms, MPI_DOUBLE_PRECISION, leadRank, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(6)
  endif
#endif

  ! each rank performs the following computation for each dataset
  ! all datasets must have identical dimensions
  ! TODO compare the dimensions of the datasets
  do i = 1, nDataSets
    ! determine the unit-cell dimensions of the original coordinate set
    vertexMin = minval(cell(1:3, atomRange(1,i):atomRange(2,i)), dim = 2) - 0.1d-8
   !vertexMax = maxval(cell(1:3, atomRange(1,i):atomRange(2,i)), dim = 2) + 0.1d-8

    ! shift particle coordinates
    do iAtom = atomRange(1,i), atomRange(2,i)
      cell(1:3,iAtom) = cell(1:3,iAtom) - vertexMin(1:3)
    enddo

    ! determine the unit-cell dimensions of the shifted coordinate set
   !vertexMin = minval(cell(1:3, atomRange(1,i):atomRange(2,i)), dim = 2) - 0.1d-8
    vertexMin = 0.0d0
    vertexMax = maxval(cell(1:3, atomRange(1,i):atomRange(2,i)), dim = 2) + 0.1d-8

    ! determine the max edge length of the input data cell
    primitiveCellEdge = dmax1(vertexMax(1)-vertexMin(1), vertexMax(2)-vertexMin(2), vertexMax(3)-vertexMin(3))

    ! lowest hierarchy level
    unitCellEdge = primitiveCellEdge * dble(nRepLayers1)
    lowestLevel = ceiling( log(unitCellEdge/cubeEdge) / log(2.0d0) )  

    ! round up the cell edge
    unitCellEdge = cubeEdge * dble (2 ** lowestLevel)
    primitiveCellEdge = unitCellEdge / dble(nRepLayers1)

    ! convert cubeEdge and primitiveCellEdge to an int value
    ! use commonMultiplier to handle the decimal value
    commonMultiplier = 0
    do while((dabs(dble(                                                      &
                   nint(cubeEdge * 10.0d0 ** commonMultiplier)                &
                    ) - cubeEdge * 10.0d0 ** commonMultiplier) > 1.0d-4) .or. &
            (dabs(dble(                                                       &
                   nint(primitiveCellEdge * 10.0d0 ** commonMultiplier)       &
                    ) - primitiveCellEdge * 10.0d0 ** commonMultiplier) > 1.0d-4) )
      commonMultiplier = commonMultiplier + 1
    enddo
    intCubeEdge = nint(cubeEdge * 10.0d0 ** commonMultiplier)
    intPrimitiveEdge = nint(primitiveCellEdge * 10.0d0 ** commonMultiplier)

    if(iam == prtRank) then
      write(*,'(/a,i0,a)') "Dataset ", i, ":"
      write(*,'(2(a,3f12.3,3x))') "Input vertexMin:", vertexMin, "vertexMax:", vertexMax
      write(*,'(  a,3f12.3    )') "Initial   edges:", &
         vertexMax(1)-vertexMin(1), vertexMax(2)-vertexMin(2), vertexMax(3)-vertexMin(3)
      write(*,'(a,2x,f16.5,a,i8)') "Input  box edge:     ", primitiveCellEdge, " Angstroms   int value ", intPrimitiveEdge
      write(*,'(a,2x,f16.5,a,i8)') "CubeEdge:            ", cubeEdge, " Angstroms   int value ", intCubeEdge
      write(*,'(a,2x,f16.5,a)') "System box edge:     ", unitCellEdge, " Angstroms"
      write(*,'(a,2x,i0)')      "Number of hierarchy levels:", lowestLevel + 1
    endif

    ! translate the lowest-left box corner to the center of coordinate system
  enddo

  if(iam == prtRank) then
    write(*,'(/a,10x,l1)') "Unit cell scaling: ", scaleToUnitSize
    write(*,'(a)') "Auto-adjustment of cubeEdge to the system size is disabled to test PBC implementation"
  endif

  ! uncomment the next line to auto-adjust the CubeEdge to the actual system size
  !cubeEdge = unitCellEdge / dble(2 ** lowestLevel)

  ! perform replication if nRepLayers > 0; copy cell to coor otherwise
  ! policy: convert the input unit-cell cuboid to a cube if replication is requested (replicate a cube)

  ! convert the dimensions of the cuboid to those of a cube before replication
  !diff(1) = (unitCellEdge - vertexMax(1) + vertexMin(1)) / 2.0d0
  !diff(2) = (unitCellEdge - vertexMax(2) + vertexMin(2)) / 2.0d0
  !diff(3) = (unitCellEdge - vertexMax(3) + vertexMin(3)) / 2.0d0
  !vertexMin(1) = vertexMin(1) - diff(1)
  !vertexMin(2) = vertexMin(2) - diff(2)
  !vertexMin(3) = vertexMin(3) - diff(3)
  !vertexMax(1) = vertexMax(1) + diff(1)
  !vertexMax(2) = vertexMax(2) + diff(2)
  !vertexMax(3) = vertexMax(3) + diff(3)

  ! replicate particle coordinates in positive direction of the x, y, z axes
  !do ijk = 0, nReplicas - 1
  !  k = ijk / nRepLayers2
  !  j = (ijk - k*nRepLayers2) / nRepLayers1
  !  i = mod(ijk, nRepLayers1)
  !  translation(1) = unitCellEdge * dble(i) ! replicate the cube
  !  translation(2) = unitCellEdge * dble(j) ! replicate the cube
  !  translation(3) = unitCellEdge * dble(k) ! replicate the cube
  !  atomStart = nAtoms * ijk
  !  do iAtom = 1, nAtoms
  !    coor(1:3,atomStart+iAtom) = cell(1:3,iAtom) + translation(1:3)
  !    q(atomStart+iAtom)        = w(iAtom)
  !  enddo
  !enddo

  !! write replicated coordinates to a file for visualization
  !open(unit=33,file="coor.xyz")
  !write(33,'(i0/)') nAtomsTotal
  !k = 1
  !do iAtom = 1, nAtomsTotal
  !  write(33,'(a2,3x,3f9.3)') atName(k), (coor(j,iAtom),j=1,3)
  !  k = k + 1
  !  if(k > nAtoms) k = 1
  !enddo
  !close(33)
  !memLocal = memLocal - dble( nAtoms*2 )
  !memMax = dmax1(memMax, (memUsed + memLocal))
  !deallocate (atName, stat=allocStat)
  !memLocal = memLocal - dble( nAtoms*2 )

  !! determine the cube edge for the system that has replicated particle coordinates
  !vertexMin = minval(coor, dim = 2) - 0.1d-8
  !vertexMax = maxval(coor, dim = 2) + 0.1d-8
  !unitCellEdge = dmax1(vertexMax(1)-vertexMin(1), vertexMax(2)-vertexMin(2), vertexMax(3)-vertexMin(3))

  !! apply scaling and save the scaling parameter
  !if(scaleToUnitSize) then
  !  reciprocal = 1.0d0 / unitCellEdge
  !  eConv = eConv * reciprocal   ! scale energy factor to account for scaled coordinates
  !
  !  ! scale the coordinates down
  !  do i = 1, nAtomsTotal
  !    coor(1,i) = coor(1,i) * reciprocal
  !    coor(2,i) = coor(2,i) * reciprocal
  !    coor(3,i) = coor(3,i) * reciprocal
  !  enddo
  !  vertexMin   = -0.5d0
  !  vertexMax   =  0.5d0
  !  cubeEdge = cubeEdge * reciprocal
  !endif

#ifdef OUTPUT_XYZQ_BINFORMAT
  ! output charges and coordinate of dataset 1 in binary format
  if(iam == leadRank) then
    ! compose the file name, fname_xyzq
    write(fname_xyzq,'(a,a)') fname_pdb(1 : max0(0, index(fname_pdb,".") - 1) ), ".xyzq"
    write(*,'(a,a)') "Creating binary file ", trim(fname_xyzq)

    ! open binary file
    open(unit=inpUnit, file=trim(fname_xyzq), form='unformatted', access='stream', status='unknown')

    ! write number of atoms from the dataset 1
    write(inpUnit) atomRange(2,1)

    ! write coordinates and charges from the dataset 1
    do iAtom = 1, atomRange(2,1)
      write(inpUnit) cell(1,iAtom), cell(2,iAtom), cell(3,iAtom), w(iAtom)
    enddo

    ! close binary file
    close(inpUnit)
  endif
#endif

#ifdef OUTPUT_VIS
  ! output coordinates of the entire (replicated) system in xyz format for the purpose of visualization
  if(iam == leadRank) then
    ! open xyz file
    open(unit=inpUnit, file="vis.xyz", form='formatted', access='stream', status='unknown')
    write(*,'(a)') "Creating xyz file  vis.xyz"

    ! set the skip parameter: write every nSkip-th atom from each dataset
    nSkip = 200

    ! compute the number of atoms to be written to the file
    nWrite = 0
    do i = 1, nDataSets
      nWrite = nWrite + ((atomRange(2,i) - atomRange(1,i) + 1) / nSkip ) * count(datasetid == i)
    enddo

    ! write the number of atoms
    write(inpUnit,'(i0/)') nWrite

    ! replicate the input data
    do ijk = 0, nReplicas - 1
      k = ijk / nRepLayers2
      j = (ijk - k*nRepLayers2) / nRepLayers1
      i = mod(ijk, nRepLayers1)
      !
      translation(1) = primitiveCellEdge * dble(i) - vertexMin(1)
      translation(2) = primitiveCellEdge * dble(j) - vertexMin(2)
      translation(3) = primitiveCellEdge * dble(k) - vertexMin(3)
      !
      if(datasetid(ijk) == 1) then
        atName = "H"
      else
        atName = "O"
      endif
      !
      do iAtom = atomRange(1,datasetid(ijk)), atomRange(2,datasetid(ijk))
        if(mod(iAtom - atomRange(1,datasetid(ijk)) + 1, nSkip) == 0) &
        write(inpUnit,'(a,3f15.5)') atName, cell(1:3,iAtom) + translation(1:3)
      enddo
    enddo

    ! close xyz file
    close(inpUnit)
  endif
#endif

  memUsed = memUsed + memLocal
  memMax = dmax1(memMax, memUsed)
  if(iam == prtRank) call MemPrint("ReadInp", memLocal, memUsed, memMax, nodeSize)

  if(iam == leadRank) then
    ! deallocate the storage for dataname since the latter is no longer needed
    memUsed = memUsed - dble(nDataSets * 512_8)
    deallocate (dataname, stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)
    memUsed = memUsed - dble(nDataSets * 512)
  endif

#ifdef TIMER
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
#endif
  if(iam == prtRank) then
    write(*,'(a,f9.1,a)') "Timer: Reading input data:", TimerCheck(2), "s"
  endif
#endif

end subroutine readinp


subroutine ReadNextArg(iArg, nArgs, type_data, out_string, out_int, out_dble)
  !
  ! Read command-line arguments
  !
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_constants, only: TYPE_FILE, TYPE_INT, TYPE_DBLE
  use common_data, only: prtRank, lxyzq
  implicit none
  !
  integer, intent(inout)         :: iArg       ! tells which argument to read; autoincremented
  integer, intent(in)            :: nArgs      ! number of arguments in the command line
  integer, intent(in)            :: type_data  ! type of data to read and return
  character (len=512), intent(out) :: out_string
  integer, intent(out)           :: out_int
  double precision, intent(out)  :: out_dble
  !
  ! local variables
  character (len=512) :: prgName
  logical (kind=1)    :: exists
  integer             :: iam
  integer             :: pos
#ifdef MPI
  integer :: mpiError, status(MPI_STATUS_SIZE)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
#endif
  !
  ! initialize output values to zero
  out_string = ""
  out_int    = 0
  out_dble   = 0.0d0
  !
  call getarg(0, prgName)  ! get program name
  !
  ! exit if no arguments provided
  if(nArgs == 0) then
    call Usage(prgName)
    call exitPrg(0)
  endif
  if(iArg > nArgs) return
  !
  ! read program argument
  call getarg(iArg, out_string)

  ! strip everything after ".xyzq"
  pos = index(out_string, ".xyzq")
  if(pos > 0) then
    out_string = out_string(:pos+4)
    lxyzq = .true.  ! the input data are in XYZQ format
  endif

#if defined (DEBUG) && defined (DETAILED0)
  if(iam == prtRank) write(*,'(a,i0,2x,a)') "iArg = ", iArg, trim(out_string)
#endif

  select case (type_data)
    case (TYPE_FILE)
      ! file argument
      inquire (file=out_string, exist=exists)
      if(.not.exists) then
        call NotExist(out_string)
        call Usage(prgName)
        call exitPrg(0)
      endif
      !
    case (TYPE_INT)
      ! integer argument
      read(out_string,*) out_int
      !
    case (TYPE_DBLE)
      ! double precision argument
      read(out_string,*) out_dble
      !
    case default
      ! wrong data type
      call Usage(prgName)
      call exitPrg(7)
  end select
  !
  ! increment argument counter
  iArg = iArg + 1
  !
end subroutine ReadNextArg
