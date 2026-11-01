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

subroutine MultipoleInit()
  !
  ! Initialize data structure for Multipole expansion
  !
  use common_data, only: nL, nSzTri, nSzWset, rExpan, sExpan, d, tri, zonal
  implicit none
  !
  integer :: iam
  integer :: allocStat
  integer :: L
  !
  ! count the number of spherical harmionics using zero-based index (0:nSzTri)
  nSzTri = (nL+1)*(nL+2)/2 - 1   ! used in allocate(0:nSzTri), therefore it needs - 1

  ! count the number of elements in Wigner set of d-matrices
  nSzWset = (nL + 1)*(2*nL + 1)*(2*nL + 3) / 3

  ! allocate memory
  allocate(rExpan(0:nSzTri), sExpan(0:nSzTri), &
           d(1:nSzWset), tri(0:nL), zonal(0:nL), stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

  ! initialize triangular pointers to rExpan terms (L,0); m=0
  tri(0) =  0  ! (L,m)=(0,0)
  do L = 1, nL 
    tri(L) = tri(L-1) + L
  enddo

end subroutine MultipoleInit


subroutine FreeMultipoles
  !
  ! Deallocate multipole expansion arrays
  !
  use common_data, only: rExpan, sExpan, d, tri, zonal
  implicit none
  !
  integer :: allocStat
  !
  ! deallocate memory
  deallocate(rExpan, sExpan, d, tri, zonal, stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)

end subroutine FreeMultipoles


subroutine RegularHarmonics(c, q, n, ox, oy, oz, chargeExpansion)
  !
  ! Compute solid harmonics
  !
  use common_data, only: nL, rExpan, tri, nSzTri, zonal
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: q(n)     ! atomic charges
  double precision, intent(in)  :: ox,oy,oz ! origin for multipole expansion
  double complex, intent(inout)    :: chargeExpansion(0:nSzTri)  ! accumulated charge expansion

  integer :: iam
  integer :: L, m, iAtom, aL, aL1, aL2, LL
  double precision :: x, y, z, x2, y2, z2, r2
  double precision :: z2lm1

  chargeExpansion = 0.0d0

  ! compute multipole expansion of charge for each particle
  do iAtom = 1, n

    x = c(1,iAtom) - ox
    y = c(2,iAtom) - oy
    z = c(3,iAtom) - oz
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r2 = dmax1(r2, 1.0d-10)

    ! initialize first few regular harmonics
    ! (0,0) = (L,m)
    zonal(0) = 1.0d0
    rExpan(0) = dcmplx(1.0d0, 0.0d0)

    ! (1,0)
    zonal(1) = z
    rExpan(1) = dcmplx(z, 0.0d0)
    ! (1,1)
    rExpan(2) = dcmplx(0.5d0*x, -0.5d0*y)

    ! (2,0)
    zonal(2) = (2.0d0*z2 - x2 - y2) * 0.25d0
    rExpan(3) = dcmplx(zonal(2), 0.0d0)
    ! (2,1)
    rExpan(4) = z * rExpan(2)
    ! (2,2)
    rExpan(5) = rExpan(2) * rExpan(2) * 0.5d0

    ! obtain higher-order harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL  = L*L

      ! (L,0)
      zonal(L) = (z2lm1 * zonal(L-1) - r2 * zonal(L-2)) / dble(L*L)
      rExpan(aL) = dcmplx(zonal(L), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        rExpan(aL+m) = (z2lm1 * rExpan(aL1+m) - r2 * rExpan(aL2+m)) / dble(LL-m*m)
      enddo

      ! (L,L-1)
      rExpan(aL+L-1) = z * rExpan(aL1+L-1)

      ! (L,L)
      rExpan(aL+L) = rExpan(aL1+L-1) * rExpan(2) / dble(L)
    enddo

    ! accumulate multipole expansion
    chargeExpansion = chargeExpansion + rExpan * q(iAtom)
  enddo

end subroutine RegularHarmonics


subroutine IrregularHarmonics(c, q, n, ox, oy, oz, potentialExpansion)
  !
  ! Compute irregular solid harmonics
  !
  use common_data, only: nL, sExpan, tri, nSzTri, zonal
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: q(n)     ! atomic charges
  double precision, intent(in)  :: ox,oy,oz ! origin for multipole expansion
  double complex, intent(inout) :: potentialExpansion(0:nSzTri)  ! accumulated potential expansion
  !
  integer :: iam
  integer :: L, m, iAtom, aL, aL1, aL2, LL1
  double precision :: x, y, z, x2, y2, z2, r, r2, ir, ir2, ir3
  double precision :: z2lm1
  double complex      :: xpy

  potentialExpansion = 0.0d0

  ! compute multipole expansion of charge for each particle
  do iAtom = 1, n

    x = c(1,iAtom) - ox
    y = c(2,iAtom) - oy
    z = c(3,iAtom) - oz
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r2 = dmax1(r2, 1.0d-10)
    r = sqrt(r2)
    ir  = 1.0d0 / r
    ir2 = ir * ir
    ir3 = ir * ir2
    xpy = dcmplx(x * ir2, y * ir2)

    ! initialize first few irregular harmonics
    ! (0,0)
    zonal(0) = ir
    sExpan(0) = dcmplx(ir, 0.0d0)

    ! (1,0)
    zonal(1) = z * ir3
    sExpan(1) = dcmplx(zonal(1), 0.0d0)
    ! (1,1)
    sExpan(2) = dcmplx(x * ir3, y * ir3)

    ! (2,0)
    zonal(2) = (2.0d0 * z2 - x2 - y2) * ir2 * ir3
    sExpan(3) = dcmplx(zonal(2), 0.0d0)
    ! (2,1)
    sExpan(4) = 3.0d0 * z * ir2 * sExpan(2)
    ! (2,2)
    sExpan(5) = 3.0d0 * xpy * sExpan(2)

    ! obtain higher-order harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL1 = (L-1)*(L-1)

      ! (L,0)
      zonal(L) = (z2lm1 * zonal(L-1) - dble(LL1) * zonal(L-2)) * ir2
      sExpan(aL) = dcmplx(zonal(L), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        sExpan(aL+m) = (z2lm1 * sExpan(aL1+m) - dble(LL1-m*m) * sExpan(aL2+m)) * ir2
      enddo

      ! (L,L-1)
      sExpan(aL+L-1) = z2lm1 * sExpan(aL1+L-1) * ir2

      ! (L,L)
      sExpan(aL+L) = dble(2*L-1) * xpy * sExpan(aL1+L-1)
    enddo

    ! accumulate multipole expansion
    potentialExpansion = potentialExpansion + sExpan * q(iAtom)
  enddo

end subroutine IrregularHarmonics


subroutine IrregularHarmonicsOMPmake(c, q, n, oxyz, potentialExpansion)
  !
  ! Compute irregular solid harmonics
  !
  use common_data, only: nL, sExpan, tri, nSzTri, zonal
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: q(n)     ! atomic charges
  double precision, intent(in)  :: oxyz(3)  ! origin for multipole expansion
  double complex, intent(inout) :: potentialExpansion(0:nSzTri)  ! accumulated potential expansion
  !
  integer :: iam
  integer :: L, m, iAtom, aL, aL1, aL2, LL1
  double precision :: x, y, z, x2, y2, z2, r, r2, ir, ir2, ir3
  double precision :: z2lm1
  double complex      :: xpy

  potentialExpansion = 0.0d0

  ! compute multipole expansion of charge for each particle
  do iAtom = 1, n

    x = c(1,iAtom) - oxyz(1)
    y = c(2,iAtom) - oxyz(2)
    z = c(3,iAtom) - oxyz(3)
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r2 = dmax1(r2, 1.0d-10)
    r = sqrt(r2)
    ir  = 1.0d0 / r
    ir2 = ir * ir
    ir3 = ir * ir2
    xpy = dcmplx(x * ir2, y * ir2)

    ! initialize first few irregular harmonics
    ! (0,0)
    zonal(0) = ir
    sExpan(0) = dcmplx(ir, 0.0d0)

    ! (1,0)
    zonal(1) = z * ir3
    sExpan(1) = dcmplx(zonal(1), 0.0d0)
    ! (1,1)
    sExpan(2) = dcmplx(x * ir3, y * ir3)

    ! (2,0)
    zonal(2) = (2.0d0 * z2 - x2 - y2) * ir2 * ir3
    sExpan(3) = dcmplx(zonal(2), 0.0d0)
    ! (2,1)
    sExpan(4) = 3.0d0 * z * ir2 * sExpan(2)
    ! (2,2)
    sExpan(5) = 3.0d0 * xpy * sExpan(2)

    ! obtain higher-order harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL1 = (L-1)*(L-1)

      ! (L,0)
      zonal(L) = (z2lm1 * zonal(L-1) - dble(LL1) * zonal(L-2)) * ir2
      sExpan(aL) = dcmplx(zonal(L), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        sExpan(aL+m) = (z2lm1 * sExpan(aL1+m) - dble(LL1-m*m) * sExpan(aL2+m)) * ir2
      enddo

      ! (L,L-1)
      sExpan(aL+L-1) = z2lm1 * sExpan(aL1+L-1) * ir2

      ! (L,L)
      sExpan(aL+L) = dble(2*L-1) * xpy * sExpan(aL1+L-1)
    enddo

    ! accumulate multipole expansion
    potentialExpansion = potentialExpansion + sExpan * q(iAtom)
  enddo

end subroutine IrregularHarmonicsOMPmake


subroutine IrregularHarmonicsSingle(c, n, ox, oy, oz, potentialExpansion)
  !
  ! Compute irregular solid harmonics
  ! Note, factor (-1)^m is not present, for simplicity
  !
  use common_data, only: nL, sExpan, tri, nSzTri, zonal
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: ox,oy,oz ! origin for multipole expansion
  double complex, intent(inout)    :: potentialExpansion(0:nSzTri)  ! accumulated potential expansion
  !
  integer :: iam
  integer :: L, m, iAtom, aL, aL1, aL2, LL1
  double precision :: x, y, z, x2, y2, z2, r, r2, ir, ir2, ir3
  double precision :: z2lm1
  double complex      :: xpy

  potentialExpansion = 0.0d0

  ! compute multipole expansion of charge for each particle
  do iAtom = 1, n

    x = c(1,iAtom) - ox
    y = c(2,iAtom) - oy
    z = c(3,iAtom) - oz
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r2 = dmax1(r2, 1.0d-10)
    r = sqrt(r2)
    ir  = 1.0d0 / r
    ir2 = ir * ir
    ir3 = ir * ir2
    xpy = dcmplx(x * ir2, y * ir2)

    ! initialize first few irregular harmonics
    ! (0,0)
    zonal(0) = ir
    sExpan(0) = dcmplx(ir, 0.0d0)

    ! (1,0)
    zonal(1) = z * ir3
    sExpan(1) = dcmplx(zonal(1), 0.0d0)
    ! (1,1)
    sExpan(2) = dcmplx(x * ir3, y * ir3)

    ! (2,0)
    zonal(2) = (2.0d0 * z2 - x2 - y2) * ir2 * ir3
    sExpan(3) = dcmplx(zonal(2), 0.0d0)
    ! (2,1)
    sExpan(4) = 3.0d0 * z * ir2 * sExpan(2)
    ! (2,2)
    sExpan(5) = 3.0d0 * xpy * sExpan(2)

    ! obtain higher-order harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL1 = (L-1)*(L-1)

      ! (L,0)
      zonal(L) = (z2lm1 * zonal(L-1) - dble(LL1) * zonal(L-2)) * ir2
      sExpan(aL) = dcmplx(zonal(L), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        sExpan(aL+m) = (z2lm1 * sExpan(aL1+m) - dble(LL1-m*m) * sExpan(aL2+m)) * ir2
      enddo

      ! (L,L-1)
      sExpan(aL+L-1) = z2lm1 * sExpan(aL1+L-1) * ir2

      ! (L,L)
      sExpan(aL+L) = dble(2*L-1) * xpy * sExpan(aL1+L-1)
    enddo

    ! accumulate multipole expansion
    potentialExpansion = potentialExpansion + sExpan
  enddo

end subroutine IrregularHarmonicsSingle


subroutine IrregularHarmonicsSingleL(c, n, ox, oy, oz, potentialExpansion, nL, nSzTri)
  !
  ! Compute irregular solid harmonics
  ! Note, factor (-1)^m is not present, for simplicity
  !
  use common_data, only: sExpan, tri, zonal
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: ox,oy,oz ! origin for multipole expansion
  integer                       :: nL, nSzTri
  double complex, intent(inout)    :: potentialExpansion(0:nSzTri)  ! accumulated potential expansion
  !
  integer :: iam
  integer :: L, m, iAtom, aL, aL1, aL2, LL1
  double precision :: x, y, z, x2, y2, z2, r, r2, ir, ir2, ir3
  double precision :: z2lm1
  double complex      :: xpy

  potentialExpansion = 0.0d0

  ! compute multipole expansion of charge for each particle
  do iAtom = 1, n

    x = c(1,iAtom) - ox
    y = c(2,iAtom) - oy
    z = c(3,iAtom) - oz
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = x2 + y2 + z2
    r2 = dmax1(r2, 1.0d-10)
    r = sqrt(r2)
    ir  = 1.0d0 / r
    ir2 = ir * ir
    ir3 = ir * ir2
    xpy = dcmplx(x * ir2, y * ir2)

    ! initialize first few irregular harmonics
    ! (0,0)
    zonal(0) = ir
    sExpan(0) = dcmplx(ir, 0.0d0)

    ! (1,0)
    zonal(1) = z * ir3
    sExpan(1) = dcmplx(zonal(1), 0.0d0)
    ! (1,1)
    sExpan(2) = dcmplx(x * ir3, y * ir3)

    ! (2,0)
    zonal(2) = (2.0d0 * z2 - x2 - y2) * ir2 * ir3
    sExpan(3) = dcmplx(zonal(2), 0.0d0)
    ! (2,1)
    sExpan(4) = 3.0d0 * z * ir2 * sExpan(2)
    ! (2,2)
    sExpan(5) = 3.0d0 * xpy * sExpan(2)

    ! obtain higher-order harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL1 = (L-1)*(L-1)

      ! (L,0)
      zonal(L) = (z2lm1 * zonal(L-1) - dble(LL1) * zonal(L-2)) * ir2
      sExpan(aL) = dcmplx(zonal(L), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        sExpan(aL+m) = (z2lm1 * sExpan(aL1+m) - dble(LL1-m*m) * sExpan(aL2+m)) * ir2
      enddo

      ! (L,L-1)
      sExpan(aL+L-1) = z2lm1 * sExpan(aL1+L-1) * ir2

      ! (L,L)
      sExpan(aL+L) = dble(2*L-1) * xpy * sExpan(aL1+L-1)
    enddo

    ! accumulate multipole expansion
    potentialExpansion = potentialExpansion + sExpan
  enddo

end subroutine IrregularHarmonicsSingleL


subroutine ChargeExpShift(Ra, Rb, newExpansion)
  !
  ! M2M operation: Shift the origin of regular multipole expansion
  !
  ! Ra - expansion to be shifted
  ! Rb - helper expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion
  !
  ! newExpansion(L,m) = sum(j=0,l) sum(k=-j,+j) { Ra(j,k) * Rb(L-j,m-k) }
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign   ! mSign(0:1) = (1.0, -1.0) 
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: Ra(0:nSzTri)
  double complex, intent(in)    :: Rb(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex :: csum
  integer     :: L, m, j, k, aJ, aLJ
  !
  ! (0,0)
  newExpansion(0) = Ra(0) * Rb(0)

  do L = 1, nL
    do m = 0, L

      csum = dcmplx(0.0d0)

      ! j <= L-m
      do j = 0, L-m   !  0 <= m <= L-j
        aJ  = tri(j)
        aLJ = tri(L-j)

        ! k < 0,  i.e. m-k > 0;  m+|k| <= L-j
        do k = 1, min0(j,L-j-m)
          csum = csum + Ra(aLJ+m+k) * dconjg(Rb(aJ+k)) * mSign(mod(k,2))
        enddo

        ! k = 0
          csum = csum + Ra(aLJ+m) * Rb(aJ)

        ! 0 < k <= m;  m-k >= 0
        do k = 1, min0(j,m)
          csum = csum + Ra(aLJ+m-k) * Rb(aJ+k)
        enddo

        ! k > m;  m-k < 0
        do k = min0(j,m)+1, min0(j,m+L-j)
          csum = csum + dconjg(Ra(aLJ+k-m)) * Rb(aJ+k) * mSign(mod(k-m,2))
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
          csum = csum + Ra(aLJ+m-k) * Rb(aJ+k)
        enddo

        ! k = m+1 ... j
        do k = m+1, min0(j,L-j+m)
          csum = csum + dconjg(Ra(aLJ+k-m)) * Rb(aJ+k) * mSign(mod(k-m,2))
        enddo

      enddo

      newExpansion(tri(L)+m) = csum
    enddo
  enddo
  !do L = 0, nL
  !  write(*,'(a,i3,2x,f20.16)') "Harmonics: ",L, aimag(newExpansion(tri(L)))
  !enddo

end subroutine ChargeExpShift


subroutine PotentialExpShift(S, R, Snew)
  !
  ! Transformation involving charge R and potential S expansions
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion; accumulate the result
  !
  ! Snew(L,m) = sum(j)sum(k=-j,j) {R(j-L,k-m) * S(j,k)}
  ! note, this R is complex-conjugate (R*)
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: S(0:nSzTri)
  double complex, intent(in)    :: R(0:nSzTri)
  double complex, intent(inout) :: Snew(0:nSzTri)
  !
  double complex      :: csum
  double precision :: rsum, mSm
  integer          :: L, m, j, k, aL, aJ, aJL
  !
  do L = 0, nL
    aL  = tri(L)

    ! m=0; then both R and S have same azimuthal number k
    !      that means we can use symmetry in the product of R and S
    !      and that their product is a real number
    ! R(j-L,k) * S(j,k)
      rsum = 0.0d0
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        ! k = 0
          rsum = rsum + dble( R(aJL) * S(aJ) ) * 0.5d0

        ! k > 0;  k-m <= j-L otherwise R=0; max(k) = nL-L = j-L
        do k = 1, j-L
          rsum = rsum + dble( R(aJL+k) * S(aJ+k) )
        enddo
      enddo
      Snew(aL) = Snew(aL) + rsum * 2.0d0

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
          csum = csum + dconjg(R(aJL+m-k)) * S(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + R(aJL) * S(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + R(aJL+k-m) * S(aJ+k)
        enddo
      enddo

      do j = L+m, nL
        aJ  = tri(j)
        aJL = tri(j-L)

        ! k < 0;  m+k <= j-L
        do k = 1, j-L-m
          csum = csum + dconjg(R(aJL+m+k)) * dconjg(S(aJ+k)) * mSm
        enddo

        ! k = 0; j >= L+m
          csum = csum + dconjg(R(aJL+m)) * S(aJ) * mSm

        ! k > 0; k < m;  m-k <= j-L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(R(aJL+m-k)) * S(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + R(aJL) * S(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + R(aJL+k-m) * S(aJ+k)
        enddo
      enddo

      Snew(aL+m) = Snew(aL+m) + csum
    enddo

  enddo

end subroutine PotentialExpShift


subroutine PotentialExpShiftOverwrite(oldExpansion, oShift, newExpansion)
  !
  ! Shift the origin of potential expansion to new location
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion; old data will be overwritten
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: oldExpansion(0:nSzTri)
  double complex, intent(in)    :: oShift(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex      :: csum
  double precision :: rsum, mSm
  integer          :: L, m, j, k, aL, aJ, aJL
  !
  do L = 0, nL
      aL  = tri(L)

      ! m=0
      rsum = 0.0d0
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        ! k = 0
          rsum = rsum + dble( oShift(aJL) * oldExpansion(aJ) ) * 0.5d0

        ! k > 0;  m+k <= j-L
        do k = 1, j-L
          rsum = rsum + dble( oShift(aJL+k) * oldExpansion(aJ+k) )
        enddo
      enddo
      newExpansion(aL) = rsum * 2.0d0

    do m = 1, L
      mSm = mSign(mod(m,2))

      csum = dcmplx(0.0d0)
      do j = L, min0(nL,L+m-1)
        aJ  = tri(j)
        aJL = tri(j-L)

        ! k < 0

        ! k = 0

        ! k > 0; k < m;  m-k <= j-L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      do j = L+m, nL
        aJ  = tri(j)
        aJL = tri(j-L)

        ! k < 0;  m+k <= j-L
        do k = 1, j-L-m
          csum = csum + dconjg(oShift(aJL+m+k)) * dconjg(oldExpansion(aJ+k)) * mSm
        enddo

        ! k = 0; j >= L+m
          csum = csum + dconjg(oShift(aJL+m)) * oldExpansion(aJ) * mSm

        ! k > 0; k < m;  m-k <= j-L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      newExpansion(aL+m) = csum
    enddo

  enddo

end subroutine PotentialExpShiftOverwrite


subroutine InvertExpansion(oldExpansion, invertedExpansion)
  !
  ! Invert multipole expansion about its origin; expansion(L,m) = (-1)^L * expansion(L,m)
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)  :: oldExpansion(0:nSzTri)
  double complex, intent(out) :: invertedExpansion(0:nSzTri)
  !
  double precision :: mS
  integer          :: L, m, aL
  !
  invertedExpansion(0) = oldExpansion(0)
  do L = 1, nL
    mS = mSign(mod(L,2))
    aL = tri(L)
    do m = 0, L
      invertedExpansion(aL+m) = oldExpansion(aL+m) * mS
    enddo
  enddo
  !
end subroutine InvertExpansion


subroutine InvertExpansionSelf(expansion)
  !
  ! Invert multipole expansion about its origin; expansion(L,m) = (-1)^L * expansion(L,m)
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(inout) :: expansion(0:nSzTri)
  !
  double precision :: mS
  integer          :: L, m, aL
  !
  do L = 1, nL
    mS = mSign(mod(L,2))
    aL = tri(L)
    do m = 0, L
      expansion(aL+m) = expansion(aL+m) * mS
    enddo
  enddo
  !
end subroutine InvertExpansionSelf


subroutine PotentialExpShiftInverted(oldExpansion, oShift, newExpansion)
  !
  ! Shift the origin of potential expansion to new location
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion; accumulate the result
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: oldExpansion(0:nSzTri)
  double complex, intent(in)    :: oShift(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex   :: csum
  double precision :: rsum, mSm, mS
  integer          :: L, m, j, k, aL, aJ, aJL
  !
  do L = 0, nL
      aL  = tri(L)

      ! m=0
      rsum = 0.0d0
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))
        ! k = 0
          rsum = rsum + mS * dble( oShift(aJL) * oldExpansion(aJ) ) * 0.5d0

        ! k > 0;  m+k <= j-L
        do k = 1, j-L
          rsum = rsum + mS * dble( oShift(aJL+k) * oldExpansion(aJ+k) )
        enddo
      enddo
      newExpansion(aL) = newExpansion(aL) + rsum * 2.0d0

    do m = 1, L
      mSm = mSign(mod(m,2))

      csum = dcmplx(0.0d0)
      do j = L, min0(nL,L+m-1)
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k < 0

        ! k = 0

        ! k > 0; k < m;  m-k <= j-L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      do j = L+m, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k < 0;  m+k <= j-L
        do k = 1, j-L-m
          csum = csum + dconjg(mS * oShift(aJL+m+k)) * dconjg(oldExpansion(aJ+k)) * mSm
        enddo

        ! k = 0; j >= L+m
          csum = csum + dconjg(mS * oShift(aJL+m)) * oldExpansion(aJ) * mSm

        ! k > 0; k < m;  m-k <= j-L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      newExpansion(aL+m) = newExpansion(aL+m) + csum
    enddo

  enddo

end subroutine PotentialExpShiftInverted


subroutine PotentialExpShiftInvertedM2L(oldExpansion, oShift, newExpansion)
  !
  ! Shift the origin of potential expansion to new location
  !
  ! unscaled form:
  ! S(l,m) = sum(j=l,inf) sum(k=-j,j) (-1)^(j-l) (j-k)!/[(l-m)! (j-l+k-m)!] R*(j-l,k-m) S(j,k)
  ! implementation uses scaled form:
  ! S(l,m) = sum(j=l,inf) sum(k=-j,j) (-1)^(j-l) R*(j-l,k-m) S(j,k)
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion; accumulate the result
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: oldExpansion(0:nSzTri)
  double complex, intent(in)    :: oShift(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex   :: csum
  double precision :: rsum, mSm, mS
  integer          :: L, m, j, k, aL, aJ, aJL
  !
  do L = 0, nL
    aL  = tri(L)

    ! m=0

      rsum = 0.0d0
      do j = L, nL+L     ! for L2L: j = L, nL  ! for M2L: j = L, nL+L since max(j-l) = nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))
        ! k = 0
          rsum = rsum + mS * dble( oShift(aJL) * oldExpansion(aJ) ) * 0.5d0

        ! |k| > 0;  k <= j-L
        do k = 1, j-L
          rsum = rsum + mS * dble( oShift(aJL+k) * oldExpansion(aJ+k) )
        enddo
      enddo
      newExpansion(aL) = newExpansion(aL) + rsum * 2.0d0   ! includes k > 0 and k < 0 for m = 0

    do m = 1, L
      mSm = mSign(mod(m,2))
      csum = dcmplx(0.0d0)

      do j = L, min0(nL+L,L+m-1)  ! to enable k > 0 so that m-k <= j-L  ! for L2L: max(j) = nL  ! for M2L: max(j) = nL+L
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k <= 0
        ! no contribution because factorial (j-l + k-m)! in denominator goes to infinity
        ! j-l = l+m-1 - l = m-1; m-1 + k-m = k-1; therefore, k must be >= 1

        ! 0 < k < m; it makes k-m negative;  m-k <= j-L; therefore k >= m-j+L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L; therefore k <= j-L+m
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      do j = L+m, nL+L  ! for L2L: max(j) = nL  ! for M2L: max(j) = nL+L
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k < 0;  it changes |k-m| <= j-L to k+m <= j-L for positive k; therefore, k runs from 1 to j-L-m
        do k = 1, j-L-m
          csum = csum + dconjg(mS * oShift(aJL+m+k)) * dconjg(oldExpansion(aJ+k)) * mSm ! may not need two dconjg
        enddo

        ! k = 0; j >= L+m; k-m is negative
          csum = csum + dconjg(mS * oShift(aJL+m)) * oldExpansion(aJ) * mSm

        ! k > 0; k < m; it keeps k-m negative, and makes m-k <= j-L for positive k
        do k = 1, m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      newExpansion(aL+m) = newExpansion(aL+m) + csum
    enddo

  enddo

end subroutine PotentialExpShiftInvertedM2L


subroutine PotentialExpShiftInvertedL2L(oldExpansion, oShift, newExpansion)
  !
  ! Shift the origin of potential expansion to new location; works for M2L and L2L
  !
  ! S(l,m) = sum(j=l,inf) sum(k=-j,j) (-1)^(j-l) (j-k)!/[(l-m)! (j-l+k-m)!] R*(j-l,k-m) S(j,k)
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion; accumulate the result
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: oldExpansion(0:nSzTri)
  double complex, intent(in)    :: oShift(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex   :: csum
  double precision :: rsum, mSm, mS
  integer          :: L, m, j, k, aL, aJ, aJL
  !
  do L = 0, nL
    aL  = tri(L)

    ! m=0

      rsum = 0.0d0
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))
        ! k = 0
          rsum = rsum + mS * dble( oShift(aJL) * oldExpansion(aJ) ) * 0.5d0

        ! k > 0;  k <= j-L
        do k = 1, j-L
          rsum = rsum + mS * dble( oShift(aJL+k) * oldExpansion(aJ+k) )
        enddo
      enddo
      newExpansion(aL) = newExpansion(aL) + rsum * 2.0d0   ! includes k > 0 and k < 0 for m = 0

    do m = 1, L
      mSm = mSign(mod(m,2))
      csum = dcmplx(0.0d0)

      do j = L, min0(nL,L+m-1)  ! to enable k > 0 so that m-k <= j-L
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k <= 0
        ! no contribution because factorial (j-l + k-m)! in denominator goes to infinity
        ! j-l = l+m-1 - l = m-1; m-1 + k-m = k-1; therefore, k must be >= 1

        ! 0 < k < m; it makes k-m negative;  m-k <= j-L; therefore k >= m-j+L
        do k = max0(1,m-j+L), m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L; therefore k <= j-L+m
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      do j = L+m, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = mSign(mod(j-L,2))

        ! k < 0;  it changes |k-m| <= j-L to k+m <= j-L for positive k; therefore, k runs from 1 to j-L-m
        do k = 1, j-L-m
          csum = csum + dconjg(mS * oShift(aJL+m+k)) * dconjg(oldExpansion(aJ+k)) * mSm ! may not need two dconjg
        enddo

        ! k = 0; j >= L+m; k-m is negative
          csum = csum + dconjg(mS * oShift(aJL+m)) * oldExpansion(aJ) * mSm

        ! k > 0; k < m; it keeps k-m negative, and makes m-k <= j-L for positive k
        do k = 1, m-1
          csum = csum + dconjg(mS * oShift(aJL+m-k)) * oldExpansion(aJ+k) * mSign(mod(m-k,2))
        enddo

        ! k = m;  m <= j
          csum = csum + mS * oShift(aJL) * oldExpansion(aJ+m)

        ! k > m;  k-m <= j-L
        do k = m+1, j-L+m
          csum = csum + mS * oShift(aJL+k-m) * oldExpansion(aJ+k)
        enddo
      enddo

      newExpansion(aL+m) = newExpansion(aL+m) + csum
    enddo

  enddo

end subroutine PotentialExpShiftInvertedL2L


subroutine M2Mz(Mb, Tz, Ma)
  !
  ! M2M translation of the origin of bar-scaled regular multipole expansion along the z-axis
  !
  ! Mb - expansion to be shifted
  ! Tz - translation distance along the z-axis
  !      negative sign - in direction of z-axis; positive - in opposite direction
  !      Tz = zB - zA for translation from box B to box A
  !      box centers have Cartesian coordinates A(0,0,zA) and B(0,0,zB)
  ! Mc - the resulting expansion
  !
  ! Ma(L,m) = sum(j=m,l) { Mb(j,m) * (Tz)^(L-j) / (L-j)! }
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)      :: Mb(0:nSzTri)
  double precision, intent(in) :: Tz
  double complex, intent(inout)   :: Ma(0:nSzTri)
  !
  integer          :: L, m, j, LmJ, aJ, aL
  double precision :: factor

  ! Optimized implementation
  ! Each new expansion element (L,m) includes the old expansion element (L,m) when j = L
  ! j = L
  Ma = Mb

  ! (L,m) = (0,0),  j = L = 0
  ! this case has already been accounted above

  ! Loop over LmJ virtually represents the loop over j
  ! This loop structure avoids repetitive computation of variable factor.

  factor = 1.0d0
  do  LmJ = 1, nL    ! loop over the values of (L-j) where j < L
    ! the multiplication factor enters into Ma(L,m) for each L >= LmJ
    factor = factor * Tz / dble(LmJ)
    do L = LmJ, nL   ! the loop structure guarantees that j < L
      j = L - LmJ
      aL = tri(L)
      aJ = tri(j)
      do m = 0, j
        Ma(aL+m) = Ma(aL+m) + Mb(aJ+m) * factor
      enddo
    enddo
  enddo

  ! Naive implementation
  !do L = 0, nL
  !  aL = tri(L)
  !  do m = 0, L
  !    Ma(aL+m) = dcmplx(0.0d0, 0.0d0)
  !    do j = m, L
  !      factor = 1.0d0
  !      do LmJ = 1, L-j
  !        factor = factor * Tz / dble(LmJ)
  !      enddo
  !      aJ = tri(j)
  !      Ma(aL+m) = Ma(aL+m) + Mb(aJ+m) * factor
  !    enddo
  !  enddo
  !enddo
end subroutine M2Mz


subroutine M2Lz(Mb, Tz, La)
  !
  ! M2L translation from box B to box A along the z-axis. Naive code.
  !
  ! La(L,m) = 1/|Tz| * sum(j=m,nL) { (-1)^(j+m) * (L+j)! * (1/Tz)^(L+j) Mb*(j,m) }
  !
  ! Mb - regular multipole expansion to be shifted
  ! Tz - translation distance
  ! La - the resulting irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)      :: Mb(0:nSzTri)
  double precision, intent(in) :: Tz
  double complex, intent(inout)   :: La(0:nSzTri)
  !
  double precision :: factor, Tzinv
  integer          :: L, m, j, aL, aJ, LpJ

  Tzinv = 1.0d0 / Tz    ! negative Tz implies translation in the positive direction to the z-axis
  factor = dabs(Tzinv)  ! l+j=0
  !
  La = dcmplx(0.0d0, 0.0d0)
  !
  ! L+j=0
  La(0) = Mb(0) * factor

  do LpJ = 1, nL
    factor = factor * dble(LpJ) * Tzinv
    do L = 0, LpJ
      j  = LpJ - L
      aJ = tri(j)
      aL = tri(L)
      ! m = 0
      La(aL) = La(aL) + Mb(aJ) * factor * mSign(mod(j,2))
      do m = 1, min0(j,L)
        La(aL+m) = La(aL+m) + dconjg(Mb(aJ+m)) * factor * mSign(mod(j+m,2))
      enddo
    enddo
  enddo

  !! extra accuracy
  !do LpJ = 1, nL+nL
  !  factor = factor * dble(LpJ) * Tzinv
  !  do L = max0(0,LpJ-nL), min0(LpJ,nL)   ! max/min insure that 0 <= L <= nL
  !    j  = LpJ - L
  !    aJ = tri(j)
  !    aL = tri(L)
  !    ! m = 0
  !    La(aL) = La(aL) + Mb(aJ) * factor * mSign(mod(j,2))
  !    do m = 1, min0(j,L)
  !      La(aL+m) = La(aL+m) + dconjg(Mb(aJ+m)) * factor * mSign(mod(j+m,2))
  !    enddo
  !  enddo
  !enddo

  !! alternative form of loop indices
  !do LpJ = 1, nL+nL
  !  factor = factor * dble(LpJ) * Tzinv
  !  do j = max0(0,LpJ-nL), min0(LpJ,nL)   ! max/min insure that 0 <= j <= nL
  !    L  = LpJ - j
  !    aJ = tri(j)
  !    aL = tri(L)
  !    ! m = 0
  !    La(aL) = La(aL) + Mb(aJ) * factor * mSign(mod(j,2))
  !    do m = 1, min0(j,L)
  !      La(aL+m) = La(aL+m) + dconjg(Mb(aJ+m)) * factor * mSign(mod(j+m,2))
  !    enddo
  !  enddo
  !enddo

  !! naive implementation
  !Tzinv = 1.0d0 / Tz    ! negative Tz implies translation in the positive direction to the z-axis
  !do L = 0, nL
  !  aL = tri(L)
  !  do m = 0, L
  !    La(aL+m) = dcmplx(0.0d0, 0.0d0)
  !    do j = m, nL
  !      aJ = tri(j)
  !      factor = dabs(Tzinv)   ! l+j=0
  !      do LpJ = 1, L+j
  !        factor = factor * dble(LpJ) * Tzinv
  !      enddo
  !      La(aL+m) = La(aL+m) + dconjg(Mb(aJ+m)) * factor * (-1.0d0)**(j+m)
  !    enddo
  !  enddo
  !enddo
  
end subroutine M2Lz


subroutine L2Lz(Lb, Tz, La)
  !
  ! L2L translation of the origin of tilde-scaled irregular multipole expansion along the z-axis
  !
  ! Lb - irregular multipole expansion to be shifted
  ! Tz - translation distance along the z-axis
  !      negative sign - in direction of z-axis; positive - in opposite direction
  !      Tz = zB - zA for translation from box B to box A
  !      box centers have Cartesian coordinates A(0,0,zA) and B(0,0,zB)
  ! La - the resulting shifted irregular multipole expansion
  !
  ! La(L,m) = sum(j=L,nL) { Lb(j,m) * (-1)^(j-L) * (Tz)^(j-L) / (j-L)! }
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)      :: Lb(0:nSzTri)
  double precision, intent(in) :: Tz
  double complex, intent(inout)   :: La(0:nSzTri)
  !
  double precision :: factor, zT
  integer          :: L, m, j, aL, aJ, JmL

  zT = -Tz  ! invert the sign of translation distance to account for (-1)^(j-L)

  La = Lb   ! handle the case of j-L=0 or j=L

  factor = 1.0d0
  do JmL = 1, nL
    factor = factor * zT / dble(JmL)
    do j = JmL, nL   ! to ensure j > L
      L  = j - JmL
      aJ = tri(j)
      aL = tri(L)
      do m = 0, L
        La(aL+m) = La(aL+m) + Lb(aJ+m) * factor
      enddo
    enddo
  enddo

  ! alternative form of loop indices
  !factor = 1.0d0
  !do JmL = 1, nL
  !  factor = factor * zT / dble(JmL)
  !  do L = 0, nL-JmL   ! to ensure j > L
  !    j  = L + JmL
  !    aJ = tri(j)
  !    aL = tri(L)
  !    do m = 0, L
  !      La(aL+m) = La(aL+m) + Lb(aJ+m) * factor
  !    enddo
  !  enddo
  !enddo

  ! naive implementation
  !do L = 0, nL
  !  aL = tri(L)
  !  do m = 0, L
  !    La(aL+m) = dcmplx(0.0d0, 0.0d0)
  !    do j = L, nL
  !      factor = 1.0d0
  !      do JmL = 1, j-L
  !        factor = factor * zT / dble(JmL)
  !      enddo
  !      aJ = tri(j)
  !      La(aL+m) = La(aL+m) + Lb(aJ+m) * factor
  !    enddo
  !  enddo
  !enddo
end subroutine L2Lz


subroutine IrregularHarmonicsLocal(x,y,z, sExpan)
  !
  ! Compute irregular solid harmonics
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  double precision, intent(in) :: x,y,z ! origin for multipole expansion
  double complex, intent(out)  :: sExpan(0:nSzTri) ! potential expansion
  !
  integer :: L, m, aL, aL1, aL2, LL1
  double precision :: x2, y2, z2, r, r2, ir, ir2, ir3
  double precision :: z2lm1
  double complex   :: xpy

  ! compute multipole expansion for a vector x,y,z

  x2 = x * x
  y2 = y * y
  z2 = z * z
  r2 = x2 + y2 + z2
  r2 = dmax1(r2, 1.0d-10)
  r = sqrt(r2)
  ir  = 1.0d0 / r
  ir2 = ir * ir
  ir3 = ir * ir2
  xpy = dcmplx(x * ir2, y * ir2)

  ! initialize first few irregular harmonics
  ! (0,0)
  sExpan(0) = dcmplx(ir, 0.0d0)

  ! (1,0)
  sExpan(1) = dcmplx(z * ir3, 0.0d0)
  ! (1,1)
  sExpan(2) = dcmplx(x * ir3, y * ir3)

  ! (2,0)
  sExpan(3) = dcmplx((2.0d0 * z2 - x2 - y2) * ir2 * ir3, 0.0d0)
  ! (2,1)
  sExpan(4) = 3.0d0 * z * ir2 * sExpan(2)
  ! (2,2)
  sExpan(5) = 3.0d0 * xpy * sExpan(2)

  ! obtain higher-order harmonics from recurrence relation
  do L = 3, nL
    z2lm1 = z * dble(2*L-1)
    aL2 = tri(L-2)
    aL1 = tri(L-1)
    aL  = tri(L)
    LL1 = (L-1)*(L-1)

    ! (L,0)
    sExpan(aL) = dcmplx((z2lm1 * dble(sExpan(aL1)) - dble(LL1) * dble(sExpan(aL2))) * ir2, 0.0d0)

    ! (L,1) ... (L,L-2)
    do m = 1, L-2
      sExpan(aL+m) = (z2lm1 * sExpan(aL1+m) - dble(LL1-m*m) * sExpan(aL2+m)) * ir2
    enddo

    ! (L,L-1)
    sExpan(aL+L-1) = z2lm1 * sExpan(aL1+L-1) * ir2

    ! (L,L)
    sExpan(aL+L) = dble(2*L-1) * xpy * sExpan(aL1+L-1)
  enddo

end subroutine IrregularHarmonicsLocal


subroutine RegularHarmonicsLocal(x, y, z, rExpan)
  !
  ! Compute regular solid harmonics
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  double precision, intent(in) :: x,y,z             ! origin for multipole expansion
  double complex, intent(out)  :: rExpan(0:nSzTri)  ! output charge expansion

  integer :: L, m, aL, aL1, aL2, LL
  double precision :: x2, y2, z2, r2
  double precision :: z2lm1

  ! compute multipole expansion for a vector x,y,z

  x2 = x * x
  y2 = y * y
  z2 = z * z
  r2 = x2 + y2 + z2
  r2 = dmax1(r2, 1.0d-10)

  ! initialize first few regular harmonics
  ! (0,0) = (L,m)
  rExpan(0) = dcmplx(1.0d0, 0.0d0)

  ! (1,0)
  rExpan(1) = dcmplx(z, 0.0d0)
  ! (1,1)
  rExpan(2) = dcmplx(0.5d0*x, -0.5d0*y)

  ! (2,0)
  rExpan(3) = dcmplx((2.0d0*z2 - x2 - y2) * 0.25d0, 0.0d0)
  ! (2,1)
  rExpan(4) = z * rExpan(2)
  ! (2,2)
  rExpan(5) = rExpan(2) * rExpan(2) * 0.5d0

  ! obtain higher-order harmonics from recurrence relation
  do L = 3, nL
    z2lm1 = z * dble(2*L-1)
    aL2 = tri(L-2)
    aL1 = tri(L-1)
    aL  = tri(L)
    LL  = L*L

    ! (L,0)
    rExpan(aL) = dcmplx((z2lm1 * dble(rExpan(aL1)) - r2 * dble(rExpan(aL2))) / dble(L*L), 0.0d0)

    ! (L,1) ... (L,L-2)
    do m = 1, L-2
      rExpan(aL+m) = (z2lm1 * rExpan(aL1+m) - r2 * rExpan(aL2+m)) / dble(LL-m*m)
    enddo

    ! (L,L-1)
    rExpan(aL+L-1) = z * rExpan(aL1+L-1)

    ! (L,L)
    rExpan(aL+L) = rExpan(aL1+L-1) * rExpan(2) / dble(L)
  enddo

end subroutine RegularHarmonicsLocal


subroutine M2M(oldExpansion, oShift, newExpansion)
  !
  ! Shift the origin of charge expansion to new location
  !
  ! oldExpansion - expansion to be shifted
  ! oShift       - expansion defining the translation (origin shift)
  ! newExpansion - the resulting expansion
  !
  use common_data, only: nL, tri, nSzTri
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: oldExpansion(0:nSzTri)
  double complex, intent(in)    :: oShift(0:nSzTri)
  double complex, intent(inout) :: newExpansion(0:nSzTri)
  !
  double complex :: csum
  integer        :: L, m, j, k, aJ, aLJ
  !
  ! (0,0)
  newExpansion(0) = newExpansion(0) + oldExpansion(0) * oShift(0)

  !$omp parallel do     &
  !$omp default(shared) &
  !$omp private(csum, aJ, aLJ)
  do L = 1, nL

    do m = 0, L

      csum = dcmplx(0.0d0)

      ! 0 <= j <= L-m
      do j = 0, L-m   ! L-j >= m  .and.  m >= 0
        aJ  = tri(j)
        aLJ = tri(L-j)

        ! k < 0;  m+|k| <= L-j
        do k = 1, min0(j,L-j-m)
          csum = csum + oldExpansion(aLJ+m+k) * dconjg(oShift(aJ+k)) * mSign(mod(k,2))
        enddo

        ! 0 <= k <= m;  k <= j; whether j or m is smaller is uncertain
        ! |m-k| <= L-j is always satisfied for k <= m .and. j <= L-m
        do k = 0, min0(j,m)
          csum = csum + oldExpansion(aLJ+m-k) * oShift(aJ+k)
        enddo

        ! k > m;  m-k < 0
        ! if j <= m, no more job to do, skip this loop
        ! it means, j >= m, so k starts from m+1 
        do k = m+1, min0(j,L-j+m)
          csum = csum + dconjg(oldExpansion(aLJ+k-m)) * oShift(aJ+k) * mSign(mod(k-m,2))
        enddo

      enddo

      do j = L-m+1, L   ! L-j < m  .and.  m > 0
        aJ  = tri(j)
        aLJ = tri(L-j)

        ! k <= 0;  m+|k| > L-j   (L-j,m-k) = 0

        ! 0 < k < m-L+j          (L-j,m-k) = 0

        ! k = m-L+j ... m
        do k = m-L+j, min0(m,j)
          csum = csum + oldExpansion(aLJ+m-k) * oShift(aJ+k)
        enddo

        ! k = m+1 ... L-j+m
        do k = m+1, min0(j,L-j+m)
          csum = csum + dconjg(oldExpansion(aLJ+k-m)) * oShift(aJ+k) * mSign(mod(k-m,2))
        enddo

      enddo

      newExpansion(tri(L)+m) = newExpansion(tri(L)+m) + csum
    enddo

  enddo
end subroutine M2M


subroutine RegularHarmonicsOMP(c, q, n, oxyz, rExpan, multipoleExpansion)
  !
  ! Compute regular multipole expansion (solid harmonics) for a set of particles
  ! Note, factor (-1)^m is removed from associated Legendre functions
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: q(n)     ! atomic charges
  double precision, intent(in)  :: oxyz(3)  ! origin for multipole expansion
  double complex, intent(out) :: rExpan(0:nSzTri)              ! temporary data storage
  double complex, intent(out) :: multipoleExpansion(0:nSzTri)  ! accumulated charge expansion

  integer :: L, m, iAtom, aL, aL1, aL2, LL
  double precision :: x, y, z, x2, y2, z2, r2
  double precision :: z2lm1

  multipoleExpansion = 0.0d0

  !$omp parallel do reduction (+:multipoleExpansion) &
  !$omp default(none) shared(c, q, tri, n, nL, oxyz, nSzTri) &
  !$omp private(x,y,z, x2,y2,z2,r2, z2lm1,aL2,aL1,aL,LL, rExpan)
  do iAtom = 1, n

    x = c(1,iAtom) - oxyz(1)
    y = c(2,iAtom) - oxyz(2)
    z = c(3,iAtom) - oxyz(3)
    x2 = x * x
    y2 = y * y
    z2 = z * z
    r2 = dmax1(x2 + y2 + z2, 1.0d-10)

    ! initialize first few regular harmonics
    ! (0,0) = (L,m)
    rExpan(0) = dcmplx(1.0d0, 0.0d0)

    ! (1,0)
    rExpan(1) = dcmplx(z, 0.0d0)
    ! (1,1)
    rExpan(2) = dcmplx(0.5d0*x, -0.5d0*y)

    ! (2,0)
    rExpan(3) = dcmplx((2.0d0*z2 - x2 - y2) * 0.25d0, 0.0d0)
    ! (2,1)
    rExpan(4) = z * rExpan(2)
    ! (2,2)
    rExpan(5) = rExpan(2) * rExpan(2) * 0.5d0

    ! obtain higher-degree harmonics from recurrence relation
    do L = 3, nL
      z2lm1 = z * dble(2*L-1)
      aL2 = tri(L-2)
      aL1 = tri(L-1)
      aL  = tri(L)
      LL  = L*L

      ! (L,0)
      rExpan(aL) = dcmplx((z2lm1 * dble(rExpan(aL1)) - r2 * dble(rExpan(aL2))) / dble(LL), 0.0d0)

      ! (L,1) ... (L,L-2)
      do m = 1, L-2
        rExpan(aL+m) = (z2lm1 * rExpan(aL1+m) - r2 * rExpan(aL2+m)) / dble(LL-m*m)
      enddo

      ! (L,L-1)
      rExpan(aL+L-1) = z * rExpan(aL1+L-1)

      ! (L,L)
      rExpan(aL+L) = rExpan(aL1+L-1) * rExpan(2) / dble(L)
    enddo

    ! accumulate multipole expansion
    multipoleExpansion = multipoleExpansion + rExpan * q(iAtom)
  enddo

end subroutine RegularHarmonicsOMP


subroutine RegularExpansion(c, q, n, oxyz, regHarmonics, regExpansion)
  !
  ! Compute regular multipole expansion for a set of point charges
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  integer, intent(in)           :: n        ! number of atoms
  double precision, intent(in)  :: c(3,n)   ! atomic coordinates
  double precision, intent(in)  :: q(n)     ! atomic charges
  double precision, intent(in)  :: oxyz(3)  ! origin for multipole expansion
  double complex, intent(out) :: regHarmonics(0:nSzTri)  ! scratch array
  double complex, intent(out) :: regExpansion(0:nSzTri)  ! accumulated expansion

  integer :: L, m, iAtom, aL, aL1, aL2, LL
  double precision :: x, y, z, r2
  double precision :: z2lm1

  regExpansion = dcmplx(0.0d0, 0.0d0)

  ! compute multipole expansion for each particle
  !$omp parallel do reduction (+:regExpansion) &
  !$omp default(none) shared(c, q, tri, n, nL, oxyz, nSzTri) &
  !$omp private(x,y,z, r2, z2lm1,aL,aL1,aL2,LL, regHarmonics)
  do iAtom = 1, n

    ! shift particle coordinates to the new coordinate system
    x  = c(1,iAtom) - oxyz(1)
    y  = c(2,iAtom) - oxyz(2)
    z  = c(3,iAtom) - oxyz(3)
    r2 = x*x + y*y + z*z

    ! initialize first few regular harmonics
    ! (0,0) = (L,m)
    regHarmonics(0) = dcmplx(1.0d0, 0.0d0)

    ! (1,0)
    regHarmonics(1) = dcmplx(z, 0.0d0)
    ! (1,1)
    regHarmonics(2) = dcmplx(-0.5d0*x, 0.5d0*y)

    ! obtain higher-degree harmonics from recurrence relations
    do L = 2, nL
      z2lm1 = z * dble(2*L-1)
      aL  = tri(L)
      aL1 = tri(L-1)
      aL2 = tri(L-2)
      LL  = L*L

      ! (L,0) ... (L,L-2)   use first recurrece relation
      do m = 0, L-2
        regHarmonics(aL+m) = (z2lm1 * regHarmonics(aL1+m) - r2 * regHarmonics(aL2+m)) / dble(LL-m*m)
      enddo

      ! (L,L-1)             use second recurrence relation
      regHarmonics(aL+L-1) = z * regHarmonics(aL1+L-1)

      ! (L,L)               use third recurrence relation
      regHarmonics(aL+L) = regHarmonics(aL1+L-1) * regHarmonics(2) / dble(L)
    enddo

    ! accumulate multipole expansion
    regExpansion = regExpansion + regHarmonics * q(iAtom)
  enddo
end subroutine RegularExpansion


subroutine DemoRegularExpansion(c, q, n, oxyz, nL, nSzTri, tri, tmpExpansion, regMultipole)
!
! Compute regular multipole expansion for a set of point charges
!
implicit none
integer,          intent(in) :: n                      ! number of particles
double precision, intent(in) :: c(3,n)                 ! particle coordinates
double precision, intent(in) :: q(n)                   ! particle charges
double precision, intent(in) :: oxyz(3)                ! box origin
integer,          intent(in) :: nL                     ! expansion truncation
integer,          intent(in) :: nSzTri                 ! triangular size
integer,          intent(in) :: tri(0:nL)              ! triangular pointers
double complex,      intent(out):: tmpExpansion(0:nSzTri) ! scratch array
double complex,      intent(out):: regMultipole(0:nSzTri) ! output

! Local variables
integer          :: L, m, iAtom, aL, aL1, aL2, LL
double precision :: x, y, z, r2
double precision :: z2lm1

regMultipole = dcmplx(0.0d0, 0.0d0)

!$omp parallel do reduction (+:regMultipole) &
!$omp default(none) shared(c, q, tri, n, nL, oxyz, nSzTri) &
!$omp private(x,y,z, r2, z2lm1,aL,aL1,aL2,LL, tmpExpansion)
do iAtom = 1, n

  ! Shift particle coordinates to the new coordinate system
  x  = c(1,iAtom) - oxyz(1) 
  y  = c(2,iAtom) - oxyz(2) 
  z  = c(3,iAtom) - oxyz(3)
  r2 = x*x + y*y + z*z

  ! Initialize first few regular harmonics
  ! (0,0) = (L,m)
  tmpExpansion(0) = dcmplx(1.0d0, 0.0d0)

  ! (1,0)
  tmpExpansion(1) = dcmplx(z, 0.0d0)
  ! (1,1)
  tmpExpansion(2) = dcmplx(-0.5d0*x, 0.5d0*y)

  ! Obtain higher-order harmonics from recurrence relations
  do L = 2, nL
    ! Initialize frequently used temporary variables
    z2lm1 = z * dble(2*L-1)
    aL  = tri(L)
    aL1 = tri(L-1)
    aL2 = tri(L-2) 
    LL  = L*L

    ! (L,0) ... (L,L-2)     use first recurrence relation
    do m = 0, L-2
      tmpExpansion(aL+m) = (z2lm1 * tmpExpansion(aL1+m) - r2 * tmpExpansion(aL2+m)) / dble(LL-m*m)
    enddo

    ! (L,L-1)               use second recurrence relation
    tmpExpansion(aL+L-1) = z * tmpExpansion(aL1+L-1)

    ! (L,L)                 use third recurrence relation
    tmpExpansion(aL+L) = tmpExpansion(aL1+L-1) * tmpExpansion(2) / dble(L)
  enddo

  ! Accumulate multipole expansion
  regMultipole = regMultipole + tmpExpansion * q(iAtom)
enddo
end subroutine DemoRegularExpansion


subroutine IrregularExpansion(c, q, n, oxyz, tmpExpansion, irrExpansion)
  !
  ! Compute irregular multipole expansion for a set of point charges
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  integer, intent(in)          :: n        ! number of atoms
  double precision, intent(in) :: c(3,n)   ! atomic coordinates
  double precision, intent(in) :: q(n)     ! atomic charges
  double precision, intent(in) :: oxyz(3)  ! origin for multipole expansion
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri)  ! scratch array
  double complex, intent(inout)   :: irrExpansion(0:nSzTri)  ! accumulated expansion
  !
  integer :: L, m, iAtom, aL, aL1, aL2, LL1
  double precision :: x, y, z, ir, ir2, ir3
  double precision :: z2lm1
  double complex      :: xpy

  irrExpansion = dcmplx(0.0d0, 0.0d0)

  ! compute multipole expansion for each particle
  !$omp parallel do reduction (+:irrExpansion) &
  !$omp default(none) shared(c, q, tri, n, nL, oxyz, nSzTri) &
  !$omp private(x,y,z, ir,ir2,ir3,xpy, z2lm1,aL,aL1,aL2,LL1, tmpExpansion)
  do iAtom = 1, n

    x = c(1,iAtom) - oxyz(1)
    y = c(2,iAtom) - oxyz(2)
    z = c(3,iAtom) - oxyz(3)
    ir  = 1.0d0 / sqrt( dmax1(x*x + y*y + z*z, 1.0d-10) )
    ir2 = ir * ir
    ir3 = ir * ir2
    xpy = dcmplx(-x * ir2, -y * ir2)

    ! initialize first few irregular harmonics
    ! (0,0)
    tmpExpansion(0) = dcmplx(ir, 0.0d0)

    ! (1,0)
    tmpExpansion(1) = dcmplx(z * ir3, 0.0d0)
    ! (1,1)
    tmpExpansion(2) = dcmplx(-x * ir3, -y * ir3)

    ! obtain higher-order harmonics from recurrence relations
    do L = 2, nL
      z2lm1 = z * dble(2*L-1)
      aL  = tri(L)
      aL1 = tri(L-1)
      aL2 = tri(L-2)
      LL1 = (L-1)*(L-1)

      ! (L,0) ... (L,L-2)   use first recurrence relation
      do m = 0, L-2
        tmpExpansion(aL+m) = (z2lm1 * tmpExpansion(aL1+m) - dble(LL1-m*m) * tmpExpansion(aL2+m)) * ir2
      enddo

      ! (L,L-1)             use second recurrence relation
      tmpExpansion(aL+L-1) = z2lm1 * tmpExpansion(aL1+L-1) * ir2

      ! (L,L)               use third recurrence relation
      tmpExpansion(aL+L) = dble(2*L-1) * xpy * tmpExpansion(aL1+L-1)
    enddo

    ! accumulate multipole expansion
    irrExpansion = irrExpansion + tmpExpansion * q(iAtom)
  enddo
end subroutine IrregularExpansion


subroutine M2Mnaive(Ma, Mb, Mc)
  !
  ! M2M operation: Shift the origin of a scaled regular multipole expansion
  ! Mc(L,m) = sum(j=0,l) sum(k=-j,+j) { Ma(j,k) * Mb(L-j,m-k) }
  !
  ! Ma - bar-scaled regular multipole expansion to be shifted
  ! Mb - bar-scaled chargeless helper expansion defining the origin shift
  ! Mc - resulting bar-scaled regular multipole expansion
  ! Ma and Mb are interchangeable
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  ! Formal arguments
  double complex, intent(in)    :: Ma(0:nSzTri)
  double complex, intent(in)    :: Mb(0:nSzTri)
  double complex, intent(inout) :: Mc(0:nSzTri)
  !
  ! Local variables
  integer :: L, m, j, k, aJ, aLJ, aLm
  !
  Mc = dcmplx(0.0d0, 0.0d0)
  !
  do L = 0, nL
    do m = 0, L
      do j = 0, L
        aJ  = tri(j)
        aLJ = tri(L-j)
        aLm = tri(L)+m
        do k = -j, j
          if(iabs(m-k) <= L-j) then
            if(k < 0 .and. m-k < 0) then
              Mc(aLm) = Mc(aLm) + dconjg(Ma(aJ+iabs(k))) * dconjg(Mb(aLJ+iabs(m-k))) * (-1.0d0)**m
            elseif(k < 0 .and. m-k >= 0) then
              Mc(aLm) = Mc(aLm) + dconjg(Ma(aJ+iabs(k))) * Mb(aLJ+m-k) * (-1.0d0)**k
            elseif(k >= 0 .and. m-k < 0) then
              Mc(aLm) = Mc(aLm) + Ma(aJ+k) * dconjg(Mb(aLJ+iabs(m-k))) * (-1.0d0)**(m-k)
            elseif(k >= 0 .and. m-k >= 0) then
              Mc(aLm) = Mc(aLm) + Ma(aJ+k) * Mb(aLJ+m-k)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  !
end subroutine M2Mnaive


subroutine M2Momp(Ma, Mb, Mc)
  !
  ! M2M operation: Shift the origin of regular multipole expansion
  ! Mc(L,m) = sum(j=0,nL) sum(k=-j,+j) { Ma(j,k) * Mb(l-j,m-k) }
  !
  ! Ma - expansion to be shifted
  ! Mb - helper expansion defining the translation (origin shift)
  ! Mc - the resulting expansion
  !
  ! Mc(L,m) = sum(j=0,l) sum(k=-j,+j) { Ma(j,k) * Mb(L-j,m-k) }
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  double complex, intent(in)    :: Ma(0:nSzTri)
  double complex, intent(in)    :: Mb(0:nSzTri)
  double complex, intent(inout) :: Mc(0:nSzTri)
  !
  integer     :: L, m, j, k, aJ, aLJ, aLm
  !
  ! (0,0)
  Mc(0) = Ma(0) * Mb(0)

  !$omp parallel do default(none) &
  !$omp shared(nL,tri,nSzTri,Ma,Mb,Mc) &
  !$omp private(L,m,j,k,aJ,aLJ,aLm)
  do L = 1, nL
    do m = 0, L
      aLm = tri(L) + m

      Mc(aLm) = dcmplx(0.0d0, 0.0d0)

      ! j <= L-m
      do j = 0, L-m   !  0 <= m <= L-j
        aJ  = tri(j)
        aLJ = tri(L-j)

        ! k < 0;  m-k > 0;  m+|k| <= L-j
        do k = 1, min0(j,L-j-m)
          Mc(aLm) = Mc(aLm) + dconjg(Ma(aJ+k)) * Mb(aLJ+m+k) * dble((-1)**k)
        enddo

        ! k = 0
          Mc(aLm) = Mc(aLm) + Ma(aJ) * Mb(aLJ+m)

        ! 0 < k <= m;  m-k >= 0
        do k = 1, min0(j,m)
          Mc(aLm) = Mc(aLm) + Ma(aJ+k) * Mb(aLJ+m-k)
        enddo

        ! k > m;  m-k < 0
        do k = m+1, min0(j,L-j+m)
          Mc(aLm) = Mc(aLm) + Ma(aJ+k) * dconjg(Mb(aLJ+k-m)) * dble((-1)**(k-m))
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
          Mc(aLm) = Mc(aLm) + Ma(aJ+k) * Mb(aLJ+m-k)
        enddo

        ! k = m+1 ... L-j+m
        do k = m+1, min0(L-j+m,j)
          Mc(aLm) = Mc(aLm) + Ma(aJ+k) * dconjg(Mb(aLJ+k-m)) * dble((-1)**(k-m))
        enddo
      enddo

    enddo  ! index m
  enddo  ! index L
  !
end subroutine M2Momp


subroutine M2Lnaive(Ma, Sb, Lc)
  !
  ! M2L operation: Transforms regular multipole expansion into irregular expansion about different origin
  ! Lc(L,m) = sum(j=0,nL) sum(k=-j,+j) { (-1)^j * Ma(j,k) * Sb(L+j,m+k) }
  !
  ! Ma - bar-scaled regular multipole expansion to be shifted
  ! Sb - tilde-scaled chargeless irregular helper expansion, which defines the origin shift
  ! Lc - resulting tilde-scaled irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  ! Formal arguments
  double complex, intent(in)    :: Ma(0:nSzTri)
  double complex, intent(in)    :: Sb(0:nSzTri)
  double complex, intent(inout) :: Lc(0:nSzTri)
  !
  ! Local variables
  integer :: L, m, j, k, aJ, aLJ, aLm
  !
  ! Zero-initialize the output array
  Lc = dcmplx(0.0d0, 0.0d0)
  !
  do L = 0, nL
    do m = 0, L
      do j = 0, nL
        if(L+j <= nL) then
          aJ  = tri(j)
          aLJ = tri(L+j)
          aLm = tri(L)+m
          do k = -j, j
            if(iabs(m+k) <= L+j) then
              if(k < 0 .and. m+k < 0) then
                Lc(aLm) = Lc(aLm) + dconjg(Ma(aJ+iabs(k))) * dconjg(Sb(aLJ+iabs(m+k))) * (-1.0d0)**(m+j)
              elseif(k < 0 .and. m+k >= 0) then
                Lc(aLm) = Lc(aLm) + dconjg(Ma(aJ+iabs(k))) * Sb(aLJ+m+k) * (-1.0d0)**(k+j)
              elseif(k >= 0 .and. m+k < 0) then
                Lc(aLm) = Lc(aLm) + Ma(aJ+k) * dconjg(Sb(aLJ+iabs(m+k))) * (-1.0d0)**(m+k+j)
              elseif(k >= 0 .and. m+k >= 0) then
                Lc(aLm) = Lc(aLm) + Ma(aJ+k) * Sb(aLJ+m+k) * (-1.0d0)**j
              endif
            endif
          enddo  ! index k
        endif
      enddo  ! index j
    enddo  ! index m
  enddo  ! index L
  !
end subroutine M2Lnaive


subroutine L2Lnaive(Ra, Lb, Lc)
  !
  ! L2L operation: Transforms irregular multipole expansion from one coordinate system to another
  ! Lc(L,m) = sum(j=L,nL) sum(k=-j,+j) { (-1)^(j-L) * Ra(j-L,k-m) * Lb(j,k) }
  !
  ! Ra - bar-scaled chargeless regular helper expansion
  ! Lb - tilde-scaled irregular helper expansion to be shifted
  ! Lc - resulting tilde-scaled irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  ! Formal arguments
  double complex, intent(in)    :: Ra(0:nSzTri)
  double complex, intent(in)    :: Lb(0:nSzTri)
  double complex, intent(inout) :: Lc(0:nSzTri)
  !
  ! Local variables
  integer :: L, m, j, k, aJ, aJL, aLm
  !
  ! Zero-initialize the output array
  Lc = dcmplx(0.0d0, 0.0d0)
  !
  do L = 0, nL
    do m = 0, L
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        aLm = tri(L)+m
        !do k = L-j+m, j-L+m
        do k = -j, j
          if(iabs(k-m) <= j-L) then
            if(k-m < 0 .and. k < 0) then
              Lc(aLm) = Lc(aLm) + dconjg(Ra(aJL+iabs(k-m))) * dconjg(Lb(aJ+iabs(k))) * (-1.0d0)**(j-L+k-m+k)
            elseif(k-m <  0 .and. k >= 0) then
              Lc(aLm) = Lc(aLm) + dconjg(Ra(aJL+iabs(k-m))) * Lb(aJ+iabs(k)) * (-1.0d0)**(j-L+k-m)
            elseif(k-m >= 0 .and. k <  0) then
              Lc(aLm) = Lc(aLm) + Ra(aJL+iabs(k-m)) * dconjg(Lb(aJ+iabs(k))) * (-1.0d0)**(j-L+k)
            elseif(k-m >= 0 .and. k >= 0) then
              Lc(aLm) = Lc(aLm) + Ra(aJL+iabs(k-m)) * Lb(aJ+iabs(k)) * (-1.0d0)**(j-L)
            endif
          endif
        enddo  ! index k
      enddo  ! index j
    enddo  ! index m
  enddo  ! index L
  !
end subroutine L2Lnaive


subroutine L2Lopt1(Ra, Sb, Sc)
  !
  ! Sc(L,m) = sum(j=L,nL) sum(k=-j,+j) { (-1)^(j-L) * Ra(j-L,k-m) * Sb(j,k) }
  ! Optimized implementation (variant 1)
  !
  ! L2L operation: 
  ! Transforms irregular multipole expansion from one coordinate system to another
  ! Ra - bar-scaled chargeless regular helper expansion
  ! Sb - tilde-scaled irregular multipole expansion to be shifted
  ! Sc - resulting tilde-scaled irregular multipole expansion
  !
  ! M2L operation: 
  ! Transforms regular multipole expansion into irregular expansion about different origin
  ! Ra - bar-scaled regular multipole expansion to be transformed
  ! Sb - tilde-scaled chargeless irregular helper expansion
  ! Sc - resulting tilde-scaled irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  double complex, intent(in)    :: Ra(0:nSzTri)
  double complex, intent(in)    :: Sb(0:nSzTri)
  double complex, intent(inout) :: Sc(0:nSzTri)
  !
  integer          :: L, m, j, k, aL, aJ, aJL, aLm
  double precision :: rsum, mSm, mS
  !
  Sc = dcmplx(0.0d0, 0.0d0);
  !
  !$omp parallel do default(none) &
  !$omp shared(nL,tri,nSzTri,Ra,Sb,Sc) &
  !$omp private(L,m,j,k,aL,aJ,aJL,aLm,rsum,mSm,mS)
  do L = 0, nL
    aL  = tri(L)

    ! m=0

      rsum = 0.0d0
      do j = L, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = dble((-1)**(j-L))
        ! k = 0
          rsum = rsum + mS * dble( Ra(aJL) * Sb(aJ) ) * 0.5d0

        ! k > 0;  k <= j-L
        do k = 1, j-L
          rsum = rsum + mS * dble( Ra(aJL+k) * Sb(aJ+k) )
        enddo
      enddo
      Sc(aL) = Sc(aL) + dcmplx(rsum * 2.0d0, 0.0d0)  ! includes k > 0 and k < 0 for m = 0

    do m = 1, L
      aLm = aL + m
      mSm = dble((-1)**m)

      do j = L, min0(nL,L+m-1)  ! to enable k > 0 so that m-k <= j-L
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = dble((-1)**(j-L))

        ! k <= 0
        ! no contribution because factorial (j-l + k-m)! in denominator goes to infinity
        ! j-l = l+m-1 - l = m-1; m-1 + k-m = k-1; therefore, k must be >= 1

        ! 0 < k < m; it makes k-m negative;  m-k <= j-L; therefore k >= m-j+L
        do k = max0(1,m-j+L), m-1
          Sc(aLm) = Sc(aLm) + dconjg(mS * Ra(aJL+m-k)) * Sb(aJ+k) * dble((-1)**(m-k))
        enddo

        ! m <= k-m <= j-L; therefore k <= j-L+m
        do k = m, j-L+m
          Sc(aLm) = Sc(aLm) + mS * Ra(aJL+k-m) * Sb(aJ+k)
        enddo
      enddo

      do j = L+m, nL
        aJ  = tri(j)
        aJL = tri(j-L)
        mS = dble((-1)**(j-L))

        ! k < 0;  it changes |k-m| <= j-L to k+m <= j-L for positive k; therefore, k runs from 1 to j-L-m
        do k = 1, j-L-m
          Sc(aLm) = Sc(aLm) + dconjg(mS * Ra(aJL+m+k)) * dconjg(Sb(aJ+k)) * mSm ! may not need two dconjg
        enddo

        ! 0 <= k < m; it keeps k-m negative, and makes m-k <= j-L for positive k
        do k = 0, m-1
          Sc(aLm) = Sc(aLm) + dconjg(mS * Ra(aJL+m-k)) * Sb(aJ+k) * dble((-1)**(m-k))
        enddo

        ! m <= k-m <= j-L
        do k = m, j-L+m
          Sc(aLm) = Sc(aLm) + mS * Ra(aJL+k-m) * Sb(aJ+k)
        enddo
      enddo

    enddo  ! index m
  enddo  ! index L
  !
end subroutine L2Lopt1


subroutine L2Lopt2(Ra, Sb, Sc)
  !
  ! Sc(L,m) = sum(j=L,nL) sum(k=-j,+j) { (-1)^(j-L) * Ra(j-L,k-m) * Sb(j,k) }
  ! Optimized implementation (variant 2)
  !
  ! L2L operation: 
  ! Transforms irregular multipole expansion from one coordinate system to another
  ! Ra - bar-scaled chargeless regular helper expansion
  ! Sb - tilde-scaled irregular multipole expansion to be shifted
  ! Sc - resulting tilde-scaled irregular multipole expansion
  !
  ! M2L operation: 
  ! Transforms regular multipole expansion into irregular expansion about different origin
  ! Ra - bar-scaled regular multipole expansion to be transformed
  ! Sb - tilde-scaled chargeless irregular helper expansion
  ! Sc - resulting tilde-scaled irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri
  implicit none
  !
  ! Formal arguments
  double complex, intent(in)    :: Ra(0:nSzTri)
  double complex, intent(in)    :: Sb(0:nSzTri)
  double complex, intent(inout) :: Sc(0:nSzTri)
  !
  ! Local variables
  integer :: L, m, j, k, aL, aJ, aJL, aLm
  double precision :: rsum
  !
  ! Zero-initialize the output array
  Sc = dcmplx(0.0d0, 0.0d0)
  !
  !$omp parallel do default(none) &
  !$omp shared(nL,tri,nSzTri,Ra,Sb,Sc) &
  !$omp private(L,m,j,k,aL,aJ,aJL,aLm,rsum)
  do L = 0, nL
    aL = tri(L)
    do j = L, nL
      aJ  = tri(j)
      aJL = tri(j-L)

      ! m=0
      ! k = 0
      rsum = dble( Ra(aJL) * Sb(aJ) ) * 0.5d0

      ! k > 0;  k <= j-L
      do k = 1, j-L
        rsum = rsum + dble( Ra(aJL+k) * Sb(aJ+k) )
      enddo
      Sc(aL) = Sc(aL) + dcmplx(rsum * 2.0d0, 0.0d0) * dble((-1)**(j-L))   ! includes k > 0 and k < 0 for m = 0

      do m = 1, min0(j-L-1,L)   ! m < j - L
        aLm = aL+m

        do k = 1, j-L-m                ! k < 0
          Sc(aLm) = Sc(aLm) + dconjg(Ra(aJL+k+m)) * dconjg(Sb(aJ+k)) * dble((-1)**(j-L-m))
        enddo

        do k = 0, min0(m-1,j-L+m)      ! 0 <= k < m
          Sc(aLm) = Sc(aLm) + dconjg(Ra(aJL-k+m)) * Sb(aJ+k) * dble((-1)**(j-L+k-m))
        enddo

        do k = m, j-L+m                ! k >= m
          Sc(aLm) = Sc(aLm) + Ra(aJL+k-m) * Sb(aJ+k) * dble((-1)**(j-L))
        enddo  ! index k
      enddo

      do m = max0(1,j-L), L             ! m >= j - L
        aLm = aL+m

        do k = L-j+m, min0(m-1,j-L+m)  ! 0 <= k < m
          Sc(aLm) = Sc(aLm) + dconjg(Ra(aJL-k+m)) * Sb(aJ+k) * dble((-1)**(j-L+k-m))
        enddo

        do k = m, j-L+m                ! k >= m
          Sc(aLm) = Sc(aLm) + Ra(aJL+k-m) * Sb(aJ+k) * dble((-1)**(j-L))
        enddo
      enddo

    enddo  ! index j
  enddo  ! index L
  !
end subroutine L2Lopt2


subroutine CompareExpansions(a, b, zero)
  ! Compare two multipole expansions, a and b
  use common_data, only: nSzTri
  implicit none
  !
  ! Formal arguments
  double complex, intent(in)      :: a(0:nSzTri)
  double complex, intent(in)      :: b(0:nSzTri)
  double precision, intent(in) :: zero
  !
  ! Local variables
  double complex :: diff
  integer    :: i

  do i = 0, nSzTri
    diff = a(i) - b(i)
    if(abs(dble(diff)) + abs(dimag(diff)) > zero) then
      write(*,'(/a,i6)') "Error: expansions do not match for array element ", i
      write(*,'(a,2f20.12)') "Expansion 1: ", dble(a(i)), dimag(a(i))
      write(*,'(a,2f20.12)') "Expansion 2: ", dble(b(i)), dimag(b(i))
      error stop
    endif
  enddo
  !
end subroutine CompareExpansions
