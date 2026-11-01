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

subroutine energy_pairwise_test(fmmEnergy, fmmPotential, fmmForce)
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  !
  use common_data, only: nAtomsTotal, bC, bQ, stdEnergy
  use common_constants, only: eConv, pConv, fConv
  use timer
  implicit none
  double precision, intent(in) :: fmmEnergy
  double precision, intent(in) :: fmmPotential(nAtomsTotal), fmmForce(3,nAtomsTotal)
  !
#ifdef MPI
  integer :: mpiError
  double precision    :: esum
#endif
  integer (kind=4)    :: iam, nRanks, allocStat, pdigits
  integer (kind=8)    :: iAtom, jAtom
  double precision    :: pot, potA, potB,  r2, rij(3)
  double precision    :: forceA(3), forceB(3)
  double precision, allocatable :: refPotential(:), refForce(:,:)
  double precision    :: calcTime
  logical             :: hError

#ifdef MPI
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
  nRanks = 1
#endif

#ifdef TIMER
  call TimerInit(2)
#endif 

  allocate(refPotential(nAtomsTotal), &
           refForce(3,nAtomsTotal), &
           stat=allocStat)

  calcTime = 0.0d0       ! time spent in doing conventional calculation

  if(nAtomsTotal == 59430) then 
    stdEnergy = -4.146567836684987506E+05  !  0
  elseif(nAtomsTotal == 3042816) then      !  8
    stdEnergy = -2.123203119672707920D+08
  elseif(nAtomsTotal == 24342528) then     ! 16
    stdEnergy = -1.698573436778047080D+09
  elseif(nAtomsTotal == 104454168) then    ! 26
    stdEnergy = -7.288623214552740100d+09
  elseif(nAtomsTotal == 657248256) then    ! 48
    stdEnergy = 0.0d0
  elseif(nAtomsTotal == 1043685888) then   ! 56    Bcast error
    stdEnergy = 0.0d0

  elseif(nAtomsTotal ==   1749000) then   ! 1  aon
    stdEnergy = -1.182348110712629410d+08
  elseif(nAtomsTotal == 111936000) then   ! 4  aon
    stdEnergy = -7.567769299861265180d+09

  else
    stdEnergy = 0.0d0

    ! MPI+OpenMP calculation of energy and potential
    !    accomplishes double summation of energy
    ! Compute conventional electrostatic interaction of iAtom with jAtom
    !    to produce reference data for validation of FMM computation
    ! This is not FMM computation and does not need optimization
    do iAtom = iam+1, nAtomsTotal, nRanks
      ! compute electrostatic potential in the position of particle i
      potA = 0.0d0
      forceA = 0.0d0
      !$omp parallel do reduction(+:potA,forceA) private(rij, r2, pot)
      do jAtom = 1, iAtom-1
        rij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
        r2       = rij(1)**2 + rij(2)**2 + rij(3)**2 
        pot      = bQ(jAtom) / sqrt(r2)
        potA     = potA + pot
        forceA(1:3)= forceA(1:3) + pot / r2 * rij(1:3)
      enddo
      potB = 0.0d0
      forceB = 0.0d0
      !$omp parallel do reduction(+:potB,forceB) private(rij, r2, pot)
      do jAtom = iAtom+1, nAtomsTotal
        rij(1:3) = bC(1:3,iAtom) - bC(1:3,jAtom)
        r2       = rij(1)**2 + rij(2)**2 + rij(3)**2 
        pot      = bQ(jAtom) / sqrt(r2)
        potB     = potB + pot
        forceB(1:3)= forceB(1:3) + pot / r2 * rij(1:3)
      enddo
      refForce(1:3,iAtom) = bQ(iAtom) * (forceA(1:3) + forceB(1:3))
      refPotential(iAtom) = potA + potB
      stdEnergy = stdEnergy + bQ(iAtom) * (potA + potB)
      !write(*,'(a,i6,f20.16)') "force = ", iAtom, sqrt(refForce(1,iAtom)**2 + refForce(2,iAtom)**2 + refForce(3,iAtom)**2)
    enddo

#ifdef MPI
    call MPI_AllReduce(stdEnergy, esum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(9)
    stdEnergy = esum
#endif

    stdEnergy = stdEnergy * 0.5d0 * eConv
    calcTime = TimerCheck(2)

  endif

  if(iam == 0) then
    write(*,'(/a)') "Pairwise Box Energy"
#ifdef TIMER
    write(*,'(a,22x,es26.18,12x,f9.1,a)') "std = ", stdEnergy, calcTime, "s"
#else
    write(*,'(a,22x,es26.18)')            "std = ", stdEnergy
#endif
    call show_residual(fmmEnergy, stdEnergy, pdigits)
    refPotential = refPotential * pConv  ! in Volts
    refForce = refForce * fConv          ! in nano-Newtons
    call histogram(fmmPotential, fmmForce, refPotential, refForce, hError)
  endif

  deallocate(refPotential, &
             refForce, &
             stat=allocStat)

end subroutine energy_pairwise_test


subroutine Demo_energy_pairwise_test(fmmEnergy, fmmPotential, fmmForce, pdigits, hError)
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: nAtomsTotal, bC, bQ, stdEnergy
  use common_constants, only: eConv, pConv, fConv
  use timer
  implicit none
  double precision, intent(in) :: fmmEnergy
  double precision, intent(in) :: fmmPotential(nAtomsTotal), fmmForce(3,nAtomsTotal)
  integer, intent(out)         :: pdigits
  logical, intent(out)         :: hError
  integer             :: iam, nRanks, allocStat
  double precision    :: pot, potA, potB,  r2, rij(3)
  double precision, allocatable :: refPotential(:), refForce(:,:)
  double precision    :: calcTime
  logical             :: lSkipDirectSum
  integer (kind=8)    :: iAtom, jAtom, nAtomsRank
  integer :: i, k

#ifdef MPI
  integer :: mpiError
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
  nRanks = 1
#endif

  ! Split all particles evenly among nRanks
  ! nAtomsRank can be used in non-FMM computation only
  nAtomsRank = nAtomsTotal / nRanks
  if(nAtomsRank * nRanks < nAtomsTotal) nAtomsRank = nAtomsRank + 1 ! the +1 is needed for MPI_Gather to work

  lSkipDirectSum = .TRUE.
  if(nAtomsTotal ==      3042816) then           !  7  to enable the non-FMM step change 3042816 to 13042816
    stdEnergy = -2.123236267029418945D+08
  elseif(nAtomsTotal ==  24342528) then          ! 15
    stdEnergy = -1.698601821201063871D+09
  elseif(nAtomsTotal == 194740224) then          ! 31
    stdEnergy = -1.358886578240979576E+10
  else
    lSkipDirectSum = .FALSE.
  endif

  if(lSkipDirectSum .eqv. .TRUE.) then
    ! skip the direct sum computation to save time
    if(iam == 0) then
      write(*,'(a)') "Direct Sum"
      write(*,'(a,22x,es26.18)') "std = ", stdEnergy
      call show_residual(fmmEnergy, stdEnergy, pdigits)
    endif
    hError = .FALSE.  ! assume no errors in the histogramm
  else
    ! Perform the direct sum computation
    call TimerInit(2)
    allocate(refPotential(nAtomsRank * nRanks), &
             refForce(3,nAtomsRank * nRanks), &
             stat=allocStat)

    ! TODO DemoDirectSum does not use nAtomsRank
    call DemoDirectSum(bC, bQ, nAtomsTotal, nAtomsRank, stdEnergy, refPotential, refForce)

    stdEnergy = stdEnergy * eConv  ! corrected against double summation inside DemoDirectSum()
    calcTime = TimerCheck(2)

    if(iam == 0) then
      write(*,'(a)') "Direct Sum"
      write(*,'(a,22x,es26.18,12x,f9.1,a)') "std = ", stdEnergy, calcTime, "s"
      call show_residual(fmmEnergy, stdEnergy, pdigits)
      refPotential = refPotential * pConv  ! in Volts
      refForce = refForce * fConv          ! in nano-Newtons
      !k = 0
      !do i = 1, nAtomsTotal
      !  if(dabs(fmmPotential(i)-refPotential(i)) > 0.0001d0) then
      !    k = k + 1
      !    if(k<11) write(*,'(i8, 3f20.8)') i,fmmPotential(i),refPotential(i),dabs(fmmPotential(i)-refPotential(i))
      !  endif
      !enddo
      call histogram(fmmPotential, fmmForce, refPotential, refForce, hError)
    endif

    deallocate(refPotential, &
               refForce, &
               stat=allocStat)
  endif

end subroutine Demo_energy_pairwise_test


subroutine DemoDirectSum(coor, q, nAtomsTotal, nAtomsRank, energy, potential, force)
  use common_data, only: nDevicesNode, rankAtoms, dispAtoms
  !
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  implicit none
  integer (kind=8), intent(in)  :: nAtomsTotal          ! total number of particles
  integer (kind=8), intent(in)  :: nAtomsRank           ! number of particles per MPI rank for non-FMM computation
  double precision, intent(in)  :: coor(3,nAtomsTotal)  ! particle coordinates
  double precision, intent(in)  :: q(nAtomsTotal)       ! particle charge
  double precision, intent(out) :: energy               ! energy of the system
  double precision, intent(out) :: potential(nAtomsTotal), force(3,nAtomsTotal)

#ifdef OFFL_GPU
  integer (kind=4) :: iDevice
#endif
  integer          :: allocStat
  integer          :: iam, nRanks
  integer (kind=8) :: iAtomsBegin, iAtomsEnd, iAtom, jAtom, aAtom
  double precision :: r2, rij(3), pot, field(3), energy_sum
#if defined(FETCH_POTFORCE) && defined(BCAST)
  integer          :: i
#endif
#ifdef MPI
  double precision, allocatable :: bufPotential(:), bufForce(:,:)
  integer :: mpiError
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam    = 0
  nRanks = 1
#endif

#ifdef OFFL_GPU
  iDevice = mod(iam,nDevicesNode)
#endif

  energy    = 0.0d0
  potential = 0.0d0
  force     = 0.0d0

  ! Compute electrostatic potential and force in the position of particle iAtom
  ! Use double summation to avoid writing into other processor memory
  ! It requires more FLOPS but scales better and provides faster parallel execution

  !iAtomsBegin = nAtomsRank * iam + 1
  !iAtomsEnd   = min0(nAtomsRank * iam + nAtomsRank, nAtomsTotal)
  iAtomsBegin = dispAtoms(iam) + 1
  iAtomsEnd   = dispAtoms(iam) + rankAtoms(iam)
  if(nAtomsRank /= iAtomsEnd - iAtomsBegin + 1) stop "Error in DemoDirectSum"

#ifdef MPI
  allocate(bufPotential(rankAtoms(iam)), &
           bufForce(3,rankAtoms(iam)), &
           stat=allocStat)
  bufPotential = 0.0d0
  bufForce     = 0.0d0
#endif

  ! Split the external loop in the double summation into nRanks groups
#ifdef OFFL_CPU
  !$omp parallel do &
#elif  OFFL_GPU
  !$omp target teams distribute parallel do device(iDevice) &
  !$omp map(to: coor, q) &
#ifdef MPI
  !$omp map(tofrom: bufPotential, bufForce) &
#else
  !$omp map(tofrom: potential, force) &
#endif
#endif
  !$omp private(iAtom,field,aAtom,rij,r2,pot,jAtom) &
  !$omp reduction(+:energy)
  do iAtom = iAtomsBegin, iAtomsEnd

    field = 0.0d0
    aAtom = iAtom - iAtomsBegin + 1   ! aAtom = iAtom - nAtomsRank * iam   is a local index of the same particle

    do jAtom = 1, iAtom-1
      rij(1:3)         = coor(1:3,iAtom) - coor(1:3,jAtom)
      r2               = rij(1)**2 + rij(2)**2 + rij(3)**2 
      pot              = q(jAtom) / sqrt(r2)
      field(1:3)       = field(1:3) + pot / r2 * rij(1:3)
#ifdef MPI
      bufPotential(aAtom) = bufPotential(aAtom) + pot
#else
      potential(iAtom) = potential(iAtom) + pot
#endif
    enddo

    ! skip jAtom == iAtom since a charged particle does not interact with its own charge

    do jAtom = iAtom+1, nAtomsTotal
      rij(1:3)         = coor(1:3,iAtom) - coor(1:3,jAtom)
      r2               = rij(1)**2 + rij(2)**2 + rij(3)**2 
      pot              = q(jAtom) / sqrt(r2)
      field(1:3)       = field(1:3) + pot / r2 * rij(1:3)
#ifdef MPI
      bufPotential(aAtom) = bufPotential(aAtom) + pot
#else
      potential(iAtom) = potential(iAtom) + pot
#endif
    enddo

#ifdef MPI
    bufForce(1:3,aAtom) = q(iAtom) * field(1:3)
    energy = energy + q(iAtom) * bufPotential(aAtom)
#else
    force(1:3,aAtom) = q(iAtom) * field(1:3)
    energy = energy + q(iAtom) * potential(iAtom)
#endif
  enddo

  energy = energy * 0.5d0  ! correct for double summation in energy

#ifdef MPI
  call MPI_Reduce(energy, energy_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  energy = energy_sum

#if defined(FETCH_POTFORCE) && !defined(BCAST)
  !call MPI_Gather(bufPotential, nAtomsRank, MPI_DOUBLE_PRECISION, &
  !                potential, nAtomsRank, MPI_DOUBLE_PRECISION, &
  !                0, MPI_COMM_WORLD, mpiError)
  !if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  !call MPI_Gather(bufForce, 3*nAtomsRank, MPI_DOUBLE_PRECISION, &
  !                force, 3*nAtomsRank, MPI_DOUBLE_PRECISION, &
  !                0, MPI_COMM_WORLD, mpiError)
  !if(mpiError /= MPI_SUCCESS) call exitPrg(9)

  ! retrieve potential and force to the rank 0
  call MPI_Gatherv(bufPotential, rankAtoms(iam), MPI_DOUBLE_PRECISION, &
                  potential, rankAtoms, dispAtoms, MPI_DOUBLE_PRECISION, &
                  0, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  rankAtoms = rankAtoms * 3
  dispAtoms = dispAtoms * 3
  call MPI_Gatherv(bufForce, rankAtoms(iam), MPI_DOUBLE_PRECISION, &
                  force, rankAtoms, dispAtoms, MPI_DOUBLE_PRECISION, &
                  0, MPI_COMM_WORLD, mpiError)
  if(mpiError /= MPI_SUCCESS) call exitPrg(9)
  rankAtoms = rankAtoms / 3
  dispAtoms = dispAtoms / 3
#endif

#if defined(FETCH_POTFORCE) && defined(BCAST)
  !iAtomsBegin = nAtomsRank * iam + 1
  !potential(iAtomsBegin : nAtomsRank * (iam+1)) = bufPotential(1:nAtomsRank)
  !force(1:3,iAtomsBegin : nAtomsRank * (iam+1)) = bufForce(1:3,1:nAtomsRank)
  ! 
  !do i = 0, nRanks-1
  !  iAtomsBegin = nAtomsRank * i + 1
  !  call MPI_Bcast(potential(iAtomsBegin), nAtomsRank, MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
  !  call MPI_Bcast(force(1,iAtomsBegin), nAtomsRank*3, MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
  !  if(mpiError /= MPI_SUCCESS) call exitPrg(6)
  !enddo

  ! copy the data from the temporary buffer to the right place in the destination array
  potential(iAtomsBegin : iAtomsEnd) = bufPotential(1:rankAtoms(iam))
  force(1:3,iAtomsBegin : iAtomsEnd) = bufForce(1:3,1:rankAtoms(iam))

  ! broadcast the sections of the destination array; each rank will hold a copy of the potential and force
  iAtomsBegin = 1
  do i = 0, nRanks-1
    call MPI_Bcast(potential(iAtomsBegin), rankAtoms(i), MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
    call MPI_Bcast(force(1,iAtomsBegin), rankAtoms(i)*3, MPI_DOUBLE_PRECISION, i, MPI_COMM_WORLD, mpiError)
    if(mpiError /= MPI_SUCCESS) call exitPrg(6)
    iAtomsBegin = iAtomsBegin + rankAtoms(i)
  enddo
  
#endif

! Alternative to MPI_Gather would be point-to-point communications
! Following is a possible implementation
  !if(iam > 0) then
  !  call MPI_Send(bufPotential, nAtomsRank, MPI_DOUBLE_PRECISION, 0, iam, MPI_COMM_WORLD, mpiError)
  !else
  !  potential = bufPotential   ! directly copy the send buffer to the receiver array
  !  do iRank = 1, nRanks-1
  !    iAtomsBegin = iRank * nAtomsRank + 1
  !    call MPI_Recv(potential(iAtomsBegin), nAtomsRank, MPI_DOUBLE_PRECISION, &
  !                  iRank, iRank, MPI_COMM_WORLD, mpiError)
  !  enddo
  !endif
  !if(mpiError /= MPI_SUCCESS) call exitPrg(9)

  deallocate(bufPotential, &
             bufForce, &
             stat=allocStat)
#endif
end subroutine DemoDirectSum


subroutine histogram(fmmPotential, fmmForce, refPotential, refForce, hError)
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: nAtomsTotal, bC, bQ
  use common_constants, only: eConv
  implicit none
  double precision, intent(in) :: fmmPotential(nAtomsTotal), fmmForce(3,nAtomsTotal)  ! data to test
  double precision, intent(in) :: refPotential(nAtomsTotal), refForce(3,nAtomsTotal)  ! reference data
  logical, intent(out)         :: hError
  integer, parameter          :: nhist = 16   ! number of elements in the histogram
  double precision, parameter :: hstart = 0.1d0,  dh = 0.03d0   ! histogram starting value and delta
  double precision            :: e(nhist), hdelta ! see initialization bellow
  integer                     :: rPhistogram(nhist), rFhistogram(nhist)  ! relative % error
  integer                     :: aPhistogram(nhist), aFhistogram(nhist)  ! absolute error

  integer             :: iam, nRanks
  integer             :: i, j, iAtom, allocStat
  double precision    :: difPotMax, difPotRel, difForceMax, difForceRel
  double precision    :: rmsdP, rmsdF, averP, averF, stdeP, stdeF
  double precision, allocatable :: relPotential,  relForce
  double precision, allocatable :: difPotential,  difForce
  double precision, parameter   :: derr = 1.0d-1  ! acceptable error in potential and force
  integer                       :: nErrors

! integer, parameter  :: nworst = 30
! integer             :: iworst(nworst), imin
! double precision    :: qA, cA(3)

#ifdef MPI
  integer :: mpiError
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, mpiError)
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
  nRanks = 1
#endif

  ! test the correctness of potential
  nErrors = 0
  do iAtom = 1, nAtomsTotal
    if(dabs(fmmPotential(iAtom) - refPotential(iAtom)) > 1.0d-3) then
      if(nErrors == 0) write(*,'(/45x,a)') "particle      FMM potential          Ref potential"
      if(nErrors < 10) write(*,'(a,i8,2x,f20.10,2x,f20.10)') "Error: FMM and Ref potentials do not match ", &
          iAtom, fmmPotential(iAtom), refPotential(iAtom)
      nErrors = nErrors + 1
    endif
  enddo
  if(nErrors > 0) write(*,"(/3(a,i0))") "Error: ", nErrors, " values out of ", nAtomsTotal, " do not match" 

  rPhistogram = 0        ! zero initialize the histogram
  rFhistogram = 0        ! zero initialize the histogram
  aPhistogram = 0        ! zero initialize the histogram
  aFhistogram = 0        ! zero initialize the histogram
  e(1) = hstart          ! histogram error x-axis
  hdelta = dh            ! histogram increment
  do i = 2 , nhist       ! initialize error x-axis
    e(i) = e(i-1) * 0.1d0
    !e(i) = e(i-1) - hdelta
    !if(e(i) < hdelta) hdelta = hdelta * 0.1d0
    !write(*,'(a,i4,f20.16)') "e(i) = ", i, e(i)
  enddo

  ! histogram analysis of potential and electric force
  difPotMax  = 0.0d0
  difPotRel  = 0.0d0
  difForceMax  = 0.0d0
  difForceRel  = 0.0d0
  rmsdP = 0.0d0
  rmsdF = 0.0d0
  averP = 0.0d0
  averF = 0.0d0
  do iAtom = 1, nAtomsTotal

    ! determine absolute unsigned difference for each data point
    difPotential = dabs(fmmPotential(iAtom) - refPotential(iAtom))
    difForce = sqrt((fmmForce(1,iAtom)-refForce(1,iAtom))**2  &
                          +(fmmForce(2,iAtom)-refForce(2,iAtom))**2  &
                          +(fmmForce(3,iAtom)-refForce(3,iAtom))**2)
    rmsdP = rmsdP + difPotential * difPotential
    averP = averP + difPotential
    rmsdF = rmsdF + difForce * difForce
    averF = averF + difForce

    ! determine relative fractional (%) unsigned difference for each data point
    relPotential = dabs( (fmmPotential(iAtom) - refPotential(iAtom)) / refPotential(iAtom) ) 
    relForce = sqrt(((fmmForce(1,iAtom)-refForce(1,iAtom))**2  &
                    +(fmmForce(2,iAtom)-refForce(2,iAtom))**2  &
                    +(fmmForce(3,iAtom)-refForce(3,iAtom))**2) &
                  / ((refForce(1,iAtom)**2 + refForce(2,iAtom)**2 + refForce(3,iAtom)**2)))

    ! map the error on the histogram
    do j = 1 , nhist
      if(relPotential < e(j))  rPhistogram(j) = rPhistogram(j) + 1  ! relative
      if(relForce     < e(j))  rFhistogram(j) = rFhistogram(j) + 1
      !phistogram(j) = phistogram(j) + min0(1, max0(0, idnint(e(j)-relPotential)))  ! to avoid using if
      if(difPotential < e(j))  aPhistogram(j) = aPhistogram(j) + 1  ! absolute
      if(difForce     < e(j))  aFhistogram(j) = aFhistogram(j) + 1
    enddo

    !if(difForce < 1.0d-10) write(*,'(a,i6,f20.16,4x,3f20.16)') "diff force: ", iAtom, difForce, &
    !force(1,iAtom)-refForce(1,iAtom), force(2,iAtom)-refForce(2,iAtom), force(3,iAtom)-refForce(3,iAtom)

    difPotMax = dmax1(difPotMax, dabs(fmmPotential(iAtom) - refPotential(iAtom)))  ! unsigned max error in absolute value
    difPotRel = dmax1(difPotRel, relPotential)                           ! unsigned max error in relative value
    difForceMax = dmax1(difForceMax, sqrt((fmmForce(1,iAtom)-refForce(1,iAtom))**2 &
                                         +(fmmForce(2,iAtom)-refForce(2,iAtom))**2 &
                                         +(fmmForce(3,iAtom)-refForce(3,iAtom))**2) )
    difForceRel = dmax1(difForceRel, relForce)
  enddo

  ! compute RMSD
  rmsdP = dsqrt(rmsdP / dble(nAtomsTotal))
  rmsdF = dsqrt(rmsdF / dble(nAtomsTotal))

  ! compute average error
  averP = averP / dble(nAtomsTotal)
  averF = averF / dble(nAtomsTotal)

  ! compute standard error
  stdeP = 0.0d0
  stdeF = 0.0d0
  do iAtom = 1, nAtomsTotal
    difPotential = dabs(fmmPotential(iAtom) - refPotential(iAtom))
    difForce = sqrt((fmmForce(1,iAtom)-refForce(1,iAtom))**2  &
                          +(fmmForce(2,iAtom)-refForce(2,iAtom))**2  &
                          +(fmmForce(3,iAtom)-refForce(3,iAtom))**2)
    stdeP = stdeP + (averP - difPotential) * (averP - difPotential)
    stdeF = stdeF + (averF - difForce)     * (averF - difForce)
  enddo
  stdeP = dsqrt(stdeP / dble(nAtomsTotal - 1))
  stdeF = dsqrt(stdeF / dble(nAtomsTotal - 1))

  ! print the histogramm for relative error
  write(*,'(/a)') "Relative error distribution normalized for number of particles (% of particles in each bin)"
  do i = 1, nhist
    write(*,'(a,i2,f20.16,2f8.2,4x,2i20)') "hist: ", i, e(i), dble(rPhistogram(i)*100)/dble(nAtomsTotal), &
                                           dble(rFhistogram(i)*100)/dble(nAtomsTotal), rPhistogram(i),rFhistogram(i)
  enddo
  write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "potential maxAbsErr maxRelErr = ", &
       difPotMax, int(dabs(dlog10(difPotMax))), difPotRel, int(dabs(dlog10(difPotRel)))
  write(*,'(a,2x,2(f20.14," (",i0,")",1x))') "force     maxAbsErr maxRelErr = ", &
       difForceMax, int(dabs(dlog10(difForceMax))), difForceRel, int(dabs(dlog10(difForceRel)))

  ! print the histogramm for absolute error
  write(*,'(/a)') "Absolute unsigned error distribution normalized for number of particles (% of particles in each bin)"
  do i = 1, nhist
    write(*,'(a,i2,f20.16,2f8.2,4x,2i20)') "hist: ", i, e(i), dble(aPhistogram(i)*100)/dble(nAtomsTotal), &
                                           dble(aFhistogram(i)*100)/dble(nAtomsTotal), aPhistogram(i),aFhistogram(i)
  enddo
  write(*,'(a,2x,3(f20.14," (",i0,")",1x))') "potential averErr rmsd stde = ", &
       averP, int(dabs(dlog10(averP))), rmsdP, int(dabs(dlog10(rmsdP))), stdeP, int(dabs(dlog10(stdeP)))
  write(*,'(a,2x,3(f20.14," (",i0,")",1x))') "force     averErr rmsd stde = ", &
       averF, int(dabs(dlog10(averF))), rmsdF, int(dabs(dlog10(rmsdF))), stdeF, int(dabs(dlog10(stdeF)))

  hError = (difPotMax > derr) .or. (difPotRel > derr) .or. (difForceMax > derr) .or. (difForceRel > derr)

  !! find 'nworst' particles having the largest error
  !do i = 1, nworst
  !  iworst(i) = i
  !enddo
  !do i = nworst+1, nAtomsTotal
  !  ! find the particle with the smallest error in the worst-error list
  !  imin = 1
  !  do j = 2, nworst
  !    if(relPotential(iworst(j)) < relPotential(iworst(imin))) imin = j
  !  enddo
  !  ! overwrite the smallest term
  !  if(relPotential(i) > relPotential(iworst(imin))) iworst(imin) = i
  !enddo
  !do i = 1, nworst
  !  write(*,'(a,i8,3f10.5)') "worst: ",iworst(i), &
  !        potential(iworst(i)), refPotential(iworst(i)), relPotential(iworst(i))
  !enddo

end subroutine histogram


subroutine show_residual(energy1, energy2, pdigits)
  implicit none
  double precision, intent(in) :: energy1, energy2
  integer, intent(out) :: pdigits
  !
  integer              :: i
  character (len=26)   :: str1, str2
  !
  write(str1,'(es26.18)') energy1
  write(str2,'(es26.18)') energy2

  pdigits = 0  ! number of matching significant digits
  do i = 1, 26
    if(str1(i:i) /= str2(i:i)) exit
    if(str1(i:i) /= ' ' .and. str1(i:i) /= '-' .and. str1(i:i) /= '.') pdigits = pdigits + 1
  enddo

  write(*,'(a,1x,f41.10," (",i0,")")') "Residual energy =", energy1-energy2, pdigits
end subroutine show_residual


subroutine MemPrint(routine, memLocal, memUsed, memMax, nRanks)
  use common_constants, only: iMB
  implicit none
  character (len=*), intent(in) :: routine
  double precision, intent(in)  :: memLocal, memUsed, memMax
  integer, intent(in)           :: nRanks
  character (len=25)            :: str
  !
  str = adjustL(trim(routine))
  write(*,'(a,a25,a,i6,a,i6,a,i6,a,i8,a)') &
          "Memory usage in ", str, "  Local: ", int(memLocal * iMB)+1,     &
                                   " Mb   Used: ",  int(memUsed * iMb)+1,  &
                                   " Mb   Max: ",   int(memMax * iMB)+1,   &
                                   " Mb   Node: ", (int(memUsed * iMB)+1) * nRanks, " Mb"
end subroutine MemPrint


subroutine MemNodePrint(routine, memLocal, memUsed, memMax, memNode)
  use common_constants, only: iMB
  implicit none
  character (len=*), intent(in) :: routine
  double precision, intent(in)  :: memLocal, memUsed, memMax, memNode
  character (len=25)            :: str
  !
  str = adjustL(trim(routine))
  write(*,'(a,a25,a,i6,a,i6,a,i6,a,i8,a)') &
          "Memory usage in ", str, "  Local: ", int(memLocal * iMB)+1,     &
                                   " Mb   Used: ",  int(memUsed * iMb)+1,  &
                                   " Mb   Max: ",   int(memMax  * iMB)+1,  &
                                   " Mb   Node: ",  int(memNode * iMB)+1, " Mb"
end subroutine MemNodePrint


subroutine associated_legendre_derivatives()
  !
  ! Compute derivatives of associated Legendre functions
  !
  implicit none
  double precision :: theta
  double precision :: u, w, der
  double precision :: p20, p21, p22
  double precision :: p30, p31, p32, p33
  double precision :: p40, p41, p42, p43, p44
  double precision :: L, m
  !
  theta = 0.6d0
  u = cos(theta)
  w = sin(theta)
  !
  L = 3.0d0
  m = 2.0d0
  !
  ! canonical associated Legendre functions
  p20 =    0.5d0 * (3.0d0 * u*u - 1.0d0)
  p21 =   -3.0d0 * w * u
  p22 =    3.0d0 * w*w
  p30 =    0.5d0 * (5.0d0 * u*u*u - 3.0d0 * u)
  p31 =   -1.5d0 * w * (5.0d0 * u*u - 1.0d0)
  p32 =   15.0d0 * w*w * u
  p33 =  -15.0d0 * w*w*w
  p40 =    0.125d0*(35.0d0 * u*u*u*u - 30.0d0 * u*u + 3.0d0)
  p41 =   -2.5d0 * w * (7.0d0 * u*u*u - 3.0d0 * u)
  p42 =    7.5d0 * w*w * (7.0d0 * u*u - 1.0d0)
  p43 = -105.0d0 * w*w*w * u
  p44 =  105.0d0 * w*w*w*w
  !
  ! canonical derivatives
  der = -m * u * p32 / (w*w) - p33 / w
  write(*,'(a,2i3,f16.6)') "L m, p der ",nint(L), nint(m), der
  !
  der = ((L+1.0d0)*u*p32 - (L+1.0d0-m)*p42) / (w*w)
  write(*,'(a,2i3,f16.6)') "L m, p der ",nint(L), nint(m), der
  !
  der = ((L+m)*p22 - L*u*p32) / (w*w)
  write(*,'(a,2i3,f16.6)') "L m, p der ",nint(L), nint(m), der
  !
  error stop
end subroutine associated_legendre_derivatives


subroutine electric_field()
  use common_data, only: nL, tri, nSzTri
  implicit none
  integer          :: allocStat
  double precision :: c(3,3), q(3), origin(3)
  double precision :: ref_e, ref_p, ref_f(3), pot, force(3)
  !
  integer          :: i, L, m
  double precision :: sump, sumf(3), apot, energy
  double precision :: x, y, z, x2, y2, z2, r2, Dr2, x2y2, r2x2y2, xr2, yr2
  double precision :: z2Dr2x2y2, xz2Dr2x2y2, yz2Dr2x2y2, xDr2, yDr2, zDr2
  double complex   :: Dr2xiy, xzDr2xiy, yzDr2xiy, xyDr2xiy, cprod
  double precision :: dL, dL1, dLm1, dm, bprod
  !
  double complex, allocatable :: mul1(:), mul2(:), mul3(:), texp(:)
  double complex, allocatable :: pot1(:), pot2(:), pot3(:)

  ! allocate memory
  allocate(mul1(0:nSzTri), &
           mul2(0:nSzTri), &
           mul3(0:nSzTri), &
           pot1(0:nSzTri), &
           pot2(0:nSzTri), &
           pot3(0:nSzTri), &
           texp(0:nSzTri), &
           stat=allocStat)

  ! initialize input data
  origin = 0.0d0
  c(1,1) = 1.1d0;  c(2,1) = 3.1d0;  c(3,1) = 0.7d0
  c(1,2) = 10.1d0; c(2,2) = 7.1d0;  c(3,2) = 5.7d0
  c(1,3) = 12.1d0; c(2,3) = 12.4d0; c(3,3) = 11.7d0
  q(1) = 2.1d0;    q(2) = 1.1d0;    q(3) = -0.5d0
  r2 = (c(1,1)-c(1,2))**2 + (c(2,1)-c(2,2))**2 + (c(3,1)-c(3,2))**2

  ! reference interaction energy for particle 1 and particle 2
  ref_e = q(1) * q(2) / sqrt(r2)

  ! reference potential in the position of particle 1 due to particle 2
  ref_p = q(2) / sqrt(r2)

  ! reference force acting on particle 1 due to particle 2 (measurement point - source point)
  ref_f(1) = q(1) * q(2) * (c(1,1)-c(1,2)) / (r2 * sqrt(r2))
  ref_f(2) = q(1) * q(2) * (c(2,1)-c(2,2)) / (r2 * sqrt(r2))
  ref_f(3) = q(1) * q(2) * (c(3,1)-c(3,2)) / (r2 * sqrt(r2))

  ! compute multipole expansion
  call RegularHarmonicsOMP(c(1,1), 1.0d0, 1, origin, texp, mul1)           ! particle 1, chargeless multipole
  call IrregularHarmonics(c(1,2), q(2), 1, 0.0d0, 0.0d0, 0.0d0, pot2)      ! particle 2

  ! force induced by particle 2 on particle 1 in the coordinate system of particle 1
  x = c(1,1);  y = c(2,1);  z = c(3,1)

  x2 = x*x;  y2 = y*y;  z2 = z*z;   r2 = x2 + y2 + z2

  xr2 = x * r2;   yr2 = y * r2;   x2y2 = x2 + y2;   r2x2y2 = r2 * x2y2

  if(x2y2 < 1.0d-4) then
    if(r2 < 1.0d-4) then
      write(*,'(/a)') "Error: Particle is too close to the center of coordinate system"
    else
      write(*,'(/a)') "Error: Particle resides on Z-axis"
    endif
    write(*,'(a)')  "Coordinate transformation is necessary"
    error stop
  endif

  Dr2  = 1.0d0 / r2
  xDr2 = x * Dr2;  yDr2 = y * Dr2;  zDr2 = z * Dr2

   z2Dr2x2y2 = z2 / r2x2y2
  xz2Dr2x2y2 = x * z2Dr2x2y2
  yz2Dr2x2y2 = y * z2Dr2x2y2

  Dr2xiy   = 1.0d0 / dcmplx(xr2, -yr2)
  xzDr2xiy =  x*z * Dr2xiy
  yzDr2xiy =  y*z * Dr2xiy
  xyDr2xiy = x2y2 * Dr2xiy

  ! L = 0; m = 0
  sump = dble(mul1(0) * pot2(0)) * 0.5d0
  sumf = 0.0d0

  do L = 1, nL
    i = tri(L)
    dL   = dble(L)
    dL1  = dble(L + 1) * 0.5d0
    !
    ! m = 0
    apot = dble(mul1(i) * pot2(i)) * 0.5d0
    sump = sump + apot
    sumf(1) = sumf(1) + dL * x/x2y2 * apot - x*z/x2y2 * dble(mul1(tri(L-1)) * pot2(i)) * 0.5d0
    sumf(2) = sumf(2) + dL * y/x2y2 * apot - y*z/x2y2 * dble(mul1(tri(L-1)) * pot2(i)) * 0.5d0
    sumf(3) = sumf(3) + dble(mul1(tri(L-1)) * pot2(i)) * 0.5d0
    !sumf(1) = sumf(1) + dL * xDr2 * apot - dL1 * dble(xzDr2xiy * mul1(i+1) * pot2(i))
    !sumf(2) = sumf(2) + dL * yDr2 * apot - dL1 * dble(yzDr2xiy * mul1(i+1) * pot2(i))
    !sumf(3) = sumf(3) + dL * zDr2 * apot + dL1 * dble(xyDr2xiy * mul1(i+1) * pot2(i))
    do m = 1, L-1
      i = i + 1
      dm = dble(m)
      dLm1 = dble(L + m + 1)
      sump = sump + dble(mul1(i) * pot2(i))
      bprod = dble(mul1(tri(L-1)+m) * pot2(i))
      cprod = mul1(i) * pot2(i)
      sumf(1) = sumf(1) + dble(dcmplx(dL*x/x2y2,  dm*y/x2y2) * cprod) - x*z/x2y2 * bprod
      sumf(2) = sumf(2) + dble(dcmplx(dL*y/x2y2, -dm*x/x2y2) * cprod) - y*z/x2y2 * bprod
      sumf(3) = sumf(3) + bprod
      !sumf(1) = sumf(1) + dble(dcmplx(dL*xDr2 + dm*x*z2Dr2x2y2,  dm*y/x2y2) * mul1(i) * pot2(i)) &
      !                  - dLm1 * dble(xzDr2xiy * mul1(i+1) * pot2(i))
      !sumf(2) = sumf(2) + dble(dcmplx(dL*yDr2 + dm*y*z2Dr2x2y2, -dm*x/x2y2) * mul1(i) * pot2(i)) &
      !                  - dLm1 * dble(yzDr2xiy * mul1(i+1) * pot2(i))
      !sumf(3) = sumf(3) + dble(L-m) * zDr2 * dble(mul1(i) * pot2(i)) &
      !                  + dLm1 * dble(xyDr2xiy * mul1(i+1) * pot2(i))
    enddo
      ! m = L
      i = i + 1
      sump = sump + dble(mul1(i) * pot2(i))
      sumf(1) = sumf(1) + dble(dcmplx(dL*x/x2y2,  dL*y/x2y2) * mul1(i) * pot2(i))
      sumf(2) = sumf(2) + dble(dcmplx(dL*y/x2y2, -dL*x/x2y2) * mul1(i) * pot2(i))
      !sumf(1) = sumf(1) + dL/x2y2 * dble(dcmplx(x, y) * mul1(i) * pot2(i)) 
      !sumf(2) = sumf(2) + dL/x2y2 * dble(dcmplx(y,-x) * mul1(i) * pot2(i)) 
      !sumf(3) is zero for m = L
  enddo
  pot = sump * 2.0d0
  force(1:3) = -q(1) * sumf(1:3) * 2.0d0
  energy = q(1) * pot

  write(*,'(a,2f20.12)')       "Energy_ref,    error = ", ref_e, energy-ref_e
  write(*,'(a,2f20.12)')       "Potential_ref, error = ", ref_p, pot-ref_p
  write(*,'(a,3(3f16.12,3x))') "Force_ref,     error = ", (ref_f(i),i=1,3), (force(i)-ref_f(i),i=1,3)

  ! deallocate memory
  deallocate(mul1, &
             mul2, &
             mul3, &
             pot1, &
             pot2, &
             pot3, &
             texp, &
             stat=allocStat)
  error stop
end subroutine electric_field


integer function idBox(x, y, z, iLevel)
  ! convert grid index x, y,z into box id
  implicit none

  ! input data
  integer, intent(in) :: x, y, z, iLevel  

  ! local variables
  integer :: cx,cy,cz, px,py,pz, nb, aLevel

  ! copy grid indices into temporary variables
  cx = x
  cy = y
  cz = z

  ! get the number of cubes along the edge at hierarchy level iLevel
  nb = 2 ** iLevel      

  ! scan the hierarchy levels from top down to iLevel
  idBox = 0
  do aLevel = 1, iLevel  
    ! get child grid indices
    cx = mod(cx,nb) ! nb is number of boxes along the edge at child level
    cy = mod(cy,nb)
    cz = mod(cz,nb)
    nb = nb / 2     ! nb is number of boxes along the edge at parent level
    ! get parent grid indices
    px = cx / nb
    py = cy / nb
    pz = cz / nb
    idBox = idBox + (pz * 4 + py * 2 + px) * nb * nb * nb
    idBox = idBox + 2 ** (3*(aLevel - 1))  ! sum up boxes at upper levels
  enddo

  ! idBox is the computed box id
end function idBox

