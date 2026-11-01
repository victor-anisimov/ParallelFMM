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

subroutine fmm
!$ use omp_lib
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  use common_data, only: memUsed, memMax, lowestLevel, M2Lene, coefflkm, nSzTri, bC, bQ, ws, &
                         nAtomsTotal, boxRegExpansion, boxIrrExpansion, nDevicesNode, nThreads, &
                         prtRank, nNodes, nodeSize, rmSize, nRanks
  use common_constants, only: eConv
  use fmm_plan
  use timer
  implicit none
  double precision :: energy
  integer :: iam, iBox, jBox, allocStat, kx,ky,kz, iAtom, jAtom
  double precision :: bx, by, bz
  double complex, allocatable :: tmpExpansion(:)
  double precision, allocatable :: aPotential(:), aForce(:,:)  ! currently used only in PBC test
  double precision :: iOrigin(3), pot, potA, potB, forceA(3), forceB(3), r2, rij(3), esum
  double precision :: EnergyPairwise
  character (len=24) :: datetime

#ifdef MPI
  integer :: mpiError
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
#endif

  ! test calculation of derivatives of associated Legendre functions
  !call associated_legendre_derivatives()

!  call TimerInit(1)
!
!  ! reset memory usage counters; using double precision to process integer values
!  memUsed = 0.0d0
!  memMax  = 0.0d0
!
!  ! read input data
!  call readinp()

  ! initialize matrix of square root coefficients (l, k, m) for the purpose of rotation
  if(M2Lene == 2) call initsq(iam)

  nDevicesNode = omp_get_num_devices()
! call omp_set_num_threads(min0(omp_get_max_threads(), 8))  ! limit the max number of threads to 8
  nThreads     = omp_get_max_threads()

  if(iam == prtRank) then
    write(*,'(/a,1x,i15)') "Number of nodes:                           ", nNodes
    write(*,'( a,1x,i15)') "Total Number of MPI ranks:                 ", nRanks
    write(*,'( a,1x,i15)') "Number of MPI ranks per node:              ", nodeSize
    write(*,'( a,1x,i15, 3x,i0,a)') &
                           "Number of MPI ranks in RM sub-communicator:", rmSize, rmSize / nodeSize, " node(s)"
    write(*,'( a,1x,i15)') "Number of RMG sub-communicators:           ", nRanks / rmSize
    write(*,'( a,1x,i15)') "Number of OMP threads per rank:            ", nThreads
    write(*,'( a,1x,i15)') "Number of GPU devices visible to the rank: ", nDevicesNode
    write(*,'( a,1x,i15)') "prtRank:                                   ", prtRank
  endif

#ifdef OFFL_GPU
  if(nDevicesNode == 0) then
    if(iam == prtRank) write(*,'(/a)') "ERROR: No GPU devices found"
    return
  endif
#endif

  if(ws+2 > 2**lowestLevel) then
    if(iam == prtRank) then
      write(*,'(/a,i8)') "Well-separation:      ", ws
      write(*,'(a,8x,i0)') "Cubes along the edge:", 2**lowestLevel
      write(*,'(/a)') "Error: The system is to small: Need ws+2 equal to or less than Cubes along the edge"
    endif
    return
  endif

#if defined (DEBUG) && defined (DETAILED0)
  !call energy_pairwise()
#endif

  ! create cubes
  call CreatePlan()

  ! initialize multipole arrays
  call MultipoleInit

  if(lowestLevel >= 2) then
    !regular FMM computation requires system decomposition into 2 or more hierarchy levels (not counting 0th level)
    !call ElectricFieldTest()          ! uses its own particle coordinates but needs nL
    !error stop
!    do iBox = 1, treeLevel(loLevel)%nboxes
!      do jBox = 1, treeLevel(loLevel)%nboxes
!        if(iBox /= jBox) then
!          call MultipoleTest(iBox,jBox,loLevel) ! boxes 1 and 31 at the (same) lowest hierarchy level
!          !call MultipoleTest(iBox,31,loLevel) ! boxes 1 and 31 at the (same) lowest hierarchy level
!        endif
!      enddo
!    enddo
!    error stop
    !call GeneratorSolidHarmonicsTest
 
    ! test calculation of electric field
    !call electric_field()

    call ElectrostaticEnergy
  else
    write(*,'(/a)') "Test of Periodic Boundary Condition"

    !allocate(box(0)%chargeExpansion(0:nSzTri), &
    !         box(0)%potentialExpansion(0:nSzTri), &
    allocate(tmpExpansion(0:nSzTri), &
             aPotential(nAtomsTotal), &
             aForce(3,nAtomsTotal), &
             stat=allocStat)
    if(allocStat /= 0) call exitPrg(2)

    ! energy of the central unit cell: use direct sum method
    energy = 0.0d0

    ! OpenMP calculation of energy, potential, and force in the central unit cell
    ! Computation of potential and force leads to double summation of energy
    ! Computation of energy, potential, and force undergoes in Gauss units
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

      aForce(1:3,iAtom) = bQ(iAtom) * (forceA(1:3) + forceB(1:3))
      aPotential(iAtom) = potA + potB
      energy = energy + bQ(iAtom) * (potA + potB)  ! energy carries double summation
    enddo

    ! OpenMP calculation of energy, potential, and force in NF using direct sum method
    ! box-box (central unit cell - image cell) interaction makes energy carry double summation
    ! there should be no force contribution from the identical replicated particle under PBC
    ! direct summation bellow includes potential and force from the mirror particles, though. How does FMM handle that?
    ! this question may be answered by applying direct summation to FF and comparing the result to FMM
    esum = 0.0d0
    !$omp parallel do reduction (+:esum) &
    !$omp private(potA, forceA, bx,by,bz, rij, r2, pot) default(shared)
    do iAtom = 1, nAtomsTotal    ! loop over atoms in the central unit cell
      potA = 0.0d0
      forceA = 0.0d0

      do kz = -ws, ws
        bz = bC(3,iAtom) - treeLevel(0)%edge * dble(kz)  ! moved the cube edge to here out of jAtom loop
        do ky = -ws, ws
          by = bC(2,iAtom) - treeLevel(0)%edge * dble(ky)
          kxCell: do kx = -ws, ws
            if(abs(kx)+abs(ky)+abs(kz) == 0) cycle kxCell  ! skip the ineraction of central unit cell with itself
     
            bx = bC(1,iAtom) - treeLevel(0)%edge * dble(kx)
            ! loop over atoms in the replicated cell 
            ! replicated coordinate is (bC(1,jAtom) + treeLevel(0)%edge * dble(kx))
            ! rij is the distance between iAtom (from central unit cell ) and jAtom (from replicated cell)
            do jAtom = 1, nAtomsTotal
              rij(1) = bx - bC(1,jAtom)
              rij(2) = by - bC(2,jAtom)
              rij(3) = bz - bC(3,jAtom)
              r2       = rij(1)**2 + rij(2)**2 + rij(3)**2
              pot      = bQ(jAtom) / sqrt(r2)
              potA     = potA + pot
              forceA(1:3)= forceA(1:3) + pot / r2 * rij(1:3)
            enddo
     
          enddo kxCell
        enddo
      enddo

      aPotential(iAtom) = aPotential(iAtom) + potA          ! potential on particle iAtom
      aForce(1:3,iAtom) = aForce(1:3,iAtom) + forceA(1:3)   ! force on particle iAtom
      esum = esum + potA * bQ(iAtom)
    enddo
    energy = energy + esum  ! part of parallel implementation

    if(iam == prtRank) write(*,'(a,es26.18)') &
      "Energy of the central unit cell + energy contribution of the near-field = ", energy * eConv * 0.5d0

    ! TODO broken; needs fixing
    !box(0)%chargeExpansion    = dcmplx(0.0d0, 0.0d0)
    !box(0)%potentialExpansion = dcmplx(0.0d0, 0.0d0)
    boxRegExpansion(0:nSzTri) = dcmplx(0.0d0, 0.0d0)
    boxIrrExpansion(0:nSzTri) = dcmplx(0.0d0, 0.0d0)
    iOrigin = 0.0d0  ! set the origin to the center of coordinate system

    ! TODO broken; needs fixing
    !call RegularExpansion(bC(1,box(0)%atomStart), bQ(box(0)%atomStart), &
    !            box(0)%nAtoms, iOrigin, tmpExpansion, box(0)%chargeExpansion)
    call RegularExpansion(bC(1,1), bQ(1), &
                nAtomsTotal, iOrigin, tmpExpansion, boxRegExpansion(0))

    ! energy carries double summation, all in Gauss units
! temporarily disable PBC
!    call LatticeSumFF(energy, aPotential, aForce)  ! lattice sum in far-field

    !deallocate(box(0)%chargeExpansion, box(0)%potentialExpansion, tmpExpansion, aPotential, aForce)
  endif

  call FreeMultipoles
  call FreeCubes

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD, mpiError)
#endif
  if(M2Lene == 2) deallocate(coefflkm)
  if(iam == prtRank) then
    call fdate(datetime)
    write(*,'(/a,f9.1,a,4x,a)') "Calculation time = ", TimerCheck(1), " seconds", datetime
    write(*,'(a)') "All Done"
  endif

! compiler independent conditional check:
! if(i == j)
! min0(1,abs(i-j))           returns 0 for .true. and 1 for .false.
! mod(min0(1,abs(i-j))+1,2)  returns 1 for .true. and 0 for .false.
end subroutine fmm
