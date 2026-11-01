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

subroutine ComputeBoxCenter(lowestVertex, edge, center, ix, iy, iz, prt)
  !
  ! Compute coordiantes of center of box ix,iy,iz
  !
  implicit none
  double precision, intent(in)  :: lowestVertex(3), edge
  double precision, intent(out) :: center(3)
  integer, intent(in)           :: ix, iy, iz
  logical, intent(in)           :: prt
  !
  center(1) = lowestVertex(1) + (dble(ix) + 0.5d0) * edge
  center(2) = lowestVertex(2) + (dble(iy) + 0.5d0) * edge
  center(3) = lowestVertex(3) + (dble(iz) + 0.5d0) * edge
  if(prt) write(*,'(a,3f10.4)') "Origin:", center(1), center(2), center(3)
  !
end subroutine ComputeBoxCenter


subroutine GetEulerAngles(b,a, jx,jy,jz, ix,iy,iz, alpha, beta)
  ! Optimized implementation
  !
  use common_constants, only: pi, machineZero
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: a(3), b(3)   ! origin of A and B boxes
  integer, intent(in) :: ix,iy,iz, jx,jy,jz    ! iBox and jBox
  double precision, intent(out) :: alpha, beta ! Euler angles
  !
  ! Local variables
  double precision :: t(3)                     ! translation axis. See Fig. 12.5
  !
  ! Determine translation axis
  t = a - b

  ! Check whether z-axis and t-vector are collinear
  if(abs(t(1)) + abs(t(2)) > machineZero) then
    ! t and z-axis are not collinear
    ! Determine angle alpha, which is subtended between the x-axis and the rotated x'-axis, xp
    ! cos(alp) = tx / sqrt(tx^2 + ty^2)
    ! if ty >= 0 then 0 <= alpha <= 180; therefore, alpha = acos(cosalp)
    ! if ty <  0 then 180 < alpha < 360; therefore, alpha = 2Pi - acos(cosalp)
    ! this leads to: alpha = [1 - sign(1,ty)] * pi + sign(1,ty) * acos(cosalp)
    alpha  = (1.0d0 - dsign(1.0d0,t(2))) * pi + dsign(1.0d0,t(2)) * dacos( t(1) / sqrt(t(1)*t(1) + t(2)*t(2)) )

    ! Determine angle beta, which is subtended between the z-axis and the translation axis
    ! cos(bet) = tz / sqrt(tx^2 + ty^2 + tz^2)
    beta = dacos( t(3) / sqrt(t(1)*t(1) + t(2)*t(2) +t(3)*t(3)) )
  else
    ! t and z-axis are collinear
    alpha = 0.0d0
    beta  = pi * (0.5d0 - dsign(0.5d0,t(3)))  ! beta = 0 for positive; beta = 180 for negative
  endif

  write(*,'(a,3i4,4x,3i4,f10.2,f10.2)') "Euler ", ix,iy,iz, jx,jy,jz, alpha*180.0d0/pi, beta*180.0d0/pi
end subroutine GetEulerAngles


subroutine GetEulerAngles2(ax,ay,az, bx,by,bz, edge, alpha, beta)
  ! Origin translation from box B (0,0,0) to box A (x,y,z)
  !
  use common_data, only: vertexMin
  use common_constants, only: pi, machineZero
  implicit none
  !
  ! Function arguments
  integer, intent(in) :: ax,ay,az, bx,by,bz    ! donor aBox and recepient bBox
  double precision, intent(in)  :: edge        ! box side length
  double precision, intent(out) :: alpha, beta ! Euler angles
  !
  ! Local variables
  double precision :: t(3)                     ! translation axis. See Fig. 12.5
  double precision :: a(3), b(3)               ! origin of boxes A and B
  !
  ! Determine translation axis
  call ComputeBoxCenter(vertexMin, edge, a, ax, ay, az, .false.)
  call ComputeBoxCenter(vertexMin, edge, b, bx, by, bz, .false.)
  t = a - b

  ! Check whether z-axis and t-vector are collinear
  if(abs(t(1)) + abs(t(2)) > machineZero) then
    ! t and z-axis are not collinear
    ! Determine angle alpha, which is subtended between the x-axis and the rotated x'-axis, xp
    ! cos(alp) = tx / sqrt(tx^2 + ty^2)
    ! if ty >= 0 then 0 <= alpha <= 180; therefore, alpha = acos(cosalp)
    ! if ty <  0 then 180 < alpha < 360; therefore, alpha = 2Pi - acos(cosalp)
    ! this leads to: alpha = [1 - sign(1,ty)] * pi + sign(1,ty) * acos(cosalp)
    alpha  = (1.0d0 - dsign(1.0d0,t(2))) * pi + dsign(1.0d0,t(2)) * dacos( t(1) / sqrt(t(1)*t(1) + t(2)*t(2)) )

    ! Determine angle beta, which is subtended between the z-axis and the translation axis
    ! cos(bet) = tz / sqrt(tx^2 + ty^2 + tz^2)
    beta = dacos( t(3) / sqrt(t(1)*t(1) + t(2)*t(2) +t(3)*t(3)) )
  else
    ! t and z-axis are collinear
    alpha = 0.0d0
    beta  = pi * (0.5d0 - dsign(0.5d0,t(3)))  ! beta = 0 for positive; beta = 180 for negative
  endif

  !write(*,'(22x,a,3i4,4x,3i4,4x,3f8.3,4x,2f8.2)') &
  ! "Euler ", ax,ay,az, bx,by,bz, t(1),t(2),t(3), alpha*180.0d0/pi, beta*180.0d0/pi
end subroutine GetEulerAngles2


subroutine GetEulerAngles3(t, alpha, beta)
  ! Origin translation from box B (0,0,0) to box A (tx,ty,tz)
  !
  use common_constants, only: pi, machineZero
  implicit none
  !
  ! Function arguments
  double precision, intent(in)  :: t(3)        ! translation vector
  double precision, intent(out) :: alpha, beta ! Euler angles

  ! Check whether z-axis and t-vector are collinear
  if(abs(t(1)) + abs(t(2)) > machineZero) then
    ! t and z-axis are not collinear
    ! Determine angle alpha, which is subtended between the x-axis and the rotated x'-axis, xp
    ! cos(alp) = tx / sqrt(tx^2 + ty^2)
    ! if ty >= 0 then 0 <= alpha <= 180; therefore, alpha = acos(cosalp)
    ! if ty <  0 then 180 < alpha < 360; therefore, alpha = 2Pi - acos(cosalp)
    ! this leads to: alpha = [1 - sign(1,ty)] * pi + sign(1,ty) * acos(cosalp)
    alpha  = (1.0d0 - dsign(1.0d0,t(2))) * pi + dsign(1.0d0,t(2)) * dacos( t(1) / sqrt(t(1)*t(1) + t(2)*t(2)) )

    ! Determine angle beta, which is subtended between the z-axis and the translation axis
    ! cos(bet) = tz / sqrt(tx^2 + ty^2 + tz^2)
    beta = dacos( t(3) / sqrt(t(1)*t(1) + t(2)*t(2) +t(3)*t(3)) )
  else
    ! t and z-axis are collinear
    alpha = 0.0d0
    beta  = pi * (0.5d0 - dsign(0.5d0,t(3)))  ! beta = 0 for positive; beta = 180 for negative
  endif

  !write(*,'(22x,a,4x,3f8.3,4x,2f8.2)') &
  ! "Euler ", t(1),t(2),t(3), alpha*180.0d0/pi, beta*180.0d0/pi
end subroutine GetEulerAngles3


subroutine M2Mrotate(tvector, inpExpansion, tmpExpansion, outExpansion, alpha, beta)
  !- Rotate coordinate system of regular multipole expansion to align it with z-axis
  !- Translate origin of expansion along z-axis
  !- Rotate back the coordinate system of resulting regular multipole expansion
  !
  use timer
  use common_data, only: nL, tri, nSzTri, d
  use common_constants, only: mSign
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(in)      :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  double precision, intent(in) :: alpha, beta             ! Euler angles
  !
  ! Local variables
  integer     :: L, m, k  ! indices of d-matrix d(L,k,m)
  !integer     :: i        ! index of element d(L,L,L)
  !integer     :: im       ! index of element d(L,L,m)
  integer     :: aL
  double complex :: csum
  double precision :: kalp, malp, distance
  !
  ! Function declarations
  integer          :: fptr
  !double precision :: Coeff
  !
  ! Compute Wigner d-matrix for angle beta
  call TimerInit(11)
  call dmatk(beta)
  write(*,'(a,f9.3,a)') "Timer: dmatk   : ", TimerCheck(11), " sec"

  call TimerInit(12)
  call dmatm(beta)
  write(*,'(a,f9.3,a)') "Timer: dmatm   : ", TimerCheck(12), " sec"

  call TimerInit(13)
  call dmatkbar(beta)
  write(*,'(a,f9.3,a)') "Timer: dmatkbar: ", TimerCheck(13), " sec"

  call TimerInit(14)
  call dmatmbar(beta)
  write(*,'(a,f9.3,a)') "Timer: dmatmbar: ", TimerCheck(14), " sec"

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )
 
  ! rotate multipole expansion
  !
  ! forward rotation of regular expansion, R*
  ! use conventional d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,0,m) * d(fptr(L,0,m)) * inpExpansion(aL)
  !    kalp = 0.0d0
  !    do k = 1, L
  !      kalp = kalp + alpha
  !      csum = csum + Coeff(L,k,m) * ( &
  !             dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
  !             dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
  !    enddo
  !    outExpansion(aL+m) = csum
  !  enddo
  !enddo

  ! d(L,k,m) index computation
  ! fptr = 1
  ! fptr = fptr + L * (2*L - 1) * (2*L + 1) / 3  ! size of the matrices up to L - 1
  ! fptr = fptr + (L - m) * (2*L + 1)            ! add shift due to index m
  ! fptr = fptr + (L - k)                        ! add shift due to index k

  ! use bar-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !i =  L * (2*L - 1) * (2*L + 1) / 3 + 1       ! index of first element d(L,L,L)
    !
    do m = 0, L
      !im = i + (L - m) * (2*L + 1)             ! add shift due to index m for d(L,L,m)
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      kalp = 0.0d0
      do k = 1, L
        kalp = kalp + alpha
        csum = csum + ( &
        dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
        dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
        !dcmplx(dcos(kalp), dsin(kalp)) * d(im+L-k) * inpExpansion(aL+k) + &                         ! d(fptr(L, k, m))
        !dcmplx(dcos(kalp),-dsin(kalp)) * d(im+L+k) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) ) ! d(fptr(L,-k, m))
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo

  ! shift along z-axis in positive direction
  call M2Mz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

  ! Compute Wigner d-matrix for angle -beta
  call dmatkbar(-beta)

  ! backward rotation of regular expansion, R* by using transposed d-matrix
  ! reuse previously computed conventional d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,0,m) * d(fptr(L,m,0)) * tmpExpansion(aL)
  !    do k = 1, L
  !      csum = csum + Coeff(L,k,m) * ( &
  !             d(fptr(L, m, k)) * tmpExpansion(aL+k) + &
  !             d(fptr(L, m,-k)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
  !    enddo
  !    malp = alpha * dble(m)
  !    outExpansion(aL+m) = csum * dcmplx(dcos(malp), -dsin(malp))
  !  enddo
  !enddo
  ! use bar-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + ( &
               d(fptr(L, k, m)) * tmpExpansion(aL+k) + &
               d(fptr(L,-k, m)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
      enddo
      malp = alpha * dble(m)
      outExpansion(aL+m) = csum * dcmplx(dcos(malp), -dsin(malp))
    enddo
  enddo
  ! reuse previously computed bar-scaled d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,0,m)**2 * d(fptr(L,m,0)) * tmpExpansion(aL)
  !    do k = 1, L
  !      csum = csum + Coeff(L,k,m)**2 * ( &
  !             d(fptr(L, m, k)) * tmpExpansion(aL+k) + &
  !             d(fptr(L, m,-k)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
  !    enddo
  !    malp = alpha * dble(m)
  !    outExpansion(aL+m) = csum * dcmplx(dcos(malp), -dsin(malp))
  !  enddo
  !enddo
end subroutine M2Mrotate


subroutine M2Mrotate2(tvector, expalp, inpExpansion, tmpExpansion, outExpansion)
  !- Rotate the coordinate system of regular multipole expansion to be translated
  !    to align its z-axis with translation axis
  !- Translate the regular multipole expansion along the z-axis
  !- Rotate the coordinate system of produced regular multipole expansion 
  !    to put its z-axis back to the original position
  !
#ifdef TIMER
  use timer
#endif
  !
  use common_data, only: nL, tri, nSzTri, d, coefflkm
  use common_constants, only: mSign
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)              ! translation vector
  double complex, intent(inout)   :: expalp(0:nL)            ! temporary array for complex exponents
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  !
  ! Local variables
  double precision :: alpha, beta      ! Euler angles
  integer          :: L, k, m          ! indices of d-matrix d(L,k,m)
  integer          :: aL, Lsum
  double complex      :: csum
  double precision :: kalp, distance
  !
  ! Function declarations
  integer          :: fptr
!#if defined (KSTD) || defined (MSTD)
!  double precision :: Coeff
!#endif

  ! Determine Euler angles
  call GetEulerAngles3(tvector, alpha, beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )

  ! Compute Wigner d-matrix for angle beta
#ifdef TIMER
  call TimerInit(50)
#endif
#ifdef KSTD
  call dmatk(beta)
#endif
#ifdef MSTD
  call dmatm(beta)
#endif
#ifdef KSCA
  call dmatkbar(beta)
#endif
#ifdef MSCA
  call dmatmbar(beta)
#endif
#ifdef MKSCA
  call dmatmbar(beta)
#endif
#ifdef TIMER
  timeM2Mdmat = timeM2Mdmat + TimerCheck(50)
#endif
 
  ! rotate multipole expansion
  ! forward rotation of regular multipole expansion, R*
#ifdef TIMER
  call TimerInit(50)
#endif

  ! compute complex exponent exp(i k alpha)
  kalp = 0.0d0
  expalp(0) = dcmplx(1.0d0, 0.0d0)
  do k = 1, nL
    kalp = kalp + alpha
    expalp(k) = dcmplx(dcos(kalp), dsin(kalp))
  enddo

#if defined (KSTD) || defined (MSTD)
  ! use conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,0,m) * d(fptr(L,0,m)) * inpExpansion(aL)
      csum = coefflkm(Lsum + m*(L+1)) * d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,k,m) * ( &
        csum = csum + coefflkm(Lsum + m*(L+1) + k) * ( &
                      expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif
 
#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use bar-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        csum = csum + expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Mdrot = timeM2Mdrot + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! shift along z-axis in positive direction
  call M2Mz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

#ifdef TIMER
  timeM2Mdztr = timeM2Mdztr + TimerCheck(50)
#endif

  ! Compute Wigner d-matrix for angle -beta
#ifdef TIMER
  call TimerInit(50)
#endif
#if defined (KSCA)
  call dmatkbar(-beta)
#endif
#if defined (MSCA) || defined (MKSCA)
  call dmatmbar(-beta)
#endif
#ifdef TIMER
  timeM2Mdmat = timeM2Mdmat + TimerCheck(50)
#endif

  ! backward rotation of regular expansion, R* by using transposed d-matrix
#ifdef TIMER
  call TimerInit(50)
#endif

#if defined (KSTD) || defined (MSTD)
  ! reuse previously computed conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,0,m) * d(fptr(L,m,0)) * tmpExpansion(aL)
      csum = coefflkm(Lsum + m*(L+1)) * d(fptr(L,m,0)) * tmpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,k,m) * ( &
        csum = csum + coefflkm(Lsum + m*(L+1) + k) * ( &
                      tmpExpansion(aL+k)  * d(fptr(L, m, k)) + &
               dconjg(tmpExpansion(aL+k)) * d(fptr(L, m,-k)) * mSign(mod(k,2)))
      enddo
      outExpansion(aL+m) = csum * dconjg(expalp(m))
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use bar-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + tmpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(tmpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum * dconjg(expalp(m))
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Mdrot = timeM2Mdrot + TimerCheck(50)
#endif
end subroutine M2Mrotate2


subroutine M2Lrotate(tvector, inpExpansion, tmpExpansion, outExpansion, alpha, beta)
  !- Rotate coordinate system of regular multipole expansion to align it with z-axis
  !- Transform regular multipole expansion into irregular one via translation along z-axis
  !- Rotate back the coordinate system of resulting irregular multipole expansion
  !
  use timer
  use common_data, only: nL, tri, nSzTri, d
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(in)      :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  double precision, intent(in) :: alpha, beta             ! Euler angles
  !
  ! Local variables
  integer     :: L, m, k
  integer     :: aL
  double complex :: csum
  double precision :: kalp, malp, distance
  !
  ! Function declarations
  integer          :: fptr
  double precision :: Coeff
  !
  ! Compute Wigner d-matrix for angle beta
  !call dmatk(beta)
  call dmatm(beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )
 
  ! rotate multipole expansion
  ! forward rotation of regular multipole expansion, R*

  ! use conventional d(beta)
  call TimerInit(15)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = Coeff(L,0,m) * d(fptr(L,0,m)) * inpExpansion(aL)
      kalp = 0.0d0
      do k = 1, L
        kalp = kalp + alpha
        csum = csum + Coeff(L,k,m) * ( &
               dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
               dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
  write(*,'(a,f9.3,a)') "Timer: rotate:   ", TimerCheck(15), " sec"

  call dmatkbar(beta)
  
  ! use bar-scaled d(beta)
  call TimerInit(15)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      kalp = 0.0d0
      do k = 1, L
        kalp = kalp + alpha
        csum = csum + ( &
               dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
               dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
  write(*,'(a,f9.3,a)') "Timer: rotatebar:", TimerCheck(15), " sec"

  ! shift along z-axis in positive direction
  call M2Lz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

  ! Compute Wigner d-matrix for angle -beta
  call TimerInit(11)
  call dmatktilde(-beta)
  write(*,'(a,f6.3,a)') "Timer: dmatktilde : ", TimerCheck(11), " sec"

  call TimerInit(12)
  call dmatmtilde(-beta)
  write(*,'(a,f6.3,a)') "Timer: dmatmtilde : ", TimerCheck(12), " sec"

  ! backward rotation of irregular multipole expansion, S
  ! reuse previously computed conventional d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,m,0) * d(fptr(L,m,0)) * tmpExpansion(aL)
  !    do k = 1, L
  !      csum = csum + Coeff(L,m,k) * ( &
  !             d(fptr(L, m, k)) * tmpExpansion(aL+k) + &
  !             d(fptr(L, m,-k)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
  !    enddo
  !    malp = alpha * dble(m)
  !    outExpansion(aL+m) = csum * dcmplx(dcos(malp), dsin(malp))
  !  enddo
  !enddo
  ! use tilde-scaled d(-beta)
  call TimerInit(15)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + ( &
               d(fptr(L, k, m)) * tmpExpansion(aL+k) + &
               d(fptr(L,-k, m)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
      enddo
      malp = alpha * dble(m)
      outExpansion(aL+m) = csum * dcmplx(dcos(malp), dsin(malp))
    enddo
  enddo
  write(*,'(a,f9.3,a)') "Timer: rotatetil:", TimerCheck(15), " sec"
end subroutine M2Lrotate


subroutine M2Lrotate2(tvector, expalp, inpExpansion, tmpExpansion, outExpansion)
  !- Rotate the coordinate system of regular multipole expansion to be translated
  !    to align its z-axis with translation axis
  !- Transform regular multipole expansion into irregular one via translation along the z-axis
  !- Rotate the coordinate system of produced irregular multipole expansion 
  !    to put its z-axis back to the original position
  !
#ifdef TIMER
  use timer
#endif
  use common_data, only: nL, tri, nSzTri, d, coefflkm
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)              ! translation vector
  double complex, intent(inout)   :: expalp(0:nL)            ! temporary array for complex exponents
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  !
  ! Local variables
  double precision :: alpha, beta      ! Euler angles
  integer          :: L, k, m          ! indices of d-matrix d(L,k,m)
  integer          :: aL, Lsum
  double complex      :: csum
  double precision :: kalp, distance
  !
  ! Function declarations
  integer          :: fptr
!#if defined (KSTD) || defined (MSTD)
!  double precision :: Coeff
!#endif

  ! Determine Euler angles
  call GetEulerAngles3(tvector, alpha, beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )

  ! Compute Wigner d-matrix for angle beta
#ifdef TIMER
  call TimerInit(50)
#endif
#ifdef KSTD
  call dmatk(beta)
#endif
#ifdef MSTD
  call dmatm(beta)
#endif
#ifdef KSCA
  call dmatkbar(beta)
#endif
#ifdef MSCA
  call dmatmbar(beta)
#endif
#ifdef MKSCA
  call dmatmbar(beta)
#endif
#ifdef TIMER
  timeM2Ldmat = timeM2Ldmat + TimerCheck(50)
#endif
 
  ! rotate multipole expansion
  ! forward rotation of regular multipole expansion, R*
#ifdef TIMER
  call TimerInit(50)
#endif

  ! compute complex exponent
  kalp = 0.0d0
  expalp(0) = dcmplx(1.0d0, 0.0d0)
  do k = 1, nL
    kalp = kalp + alpha
    expalp(k) = dcmplx(dcos(kalp), dsin(kalp))
  enddo

#if defined (KSTD) || defined (MSTD)
  ! use conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,0,m) * d(fptr(L,0,m)) * inpExpansion(aL)
      csum = coefflkm(Lsum + m*(L+1)) * d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,k,m) * ( &
        csum = csum + coefflkm(Lsum + m*(L+1) + k) * ( &
                      expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use bar-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        csum = csum + expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Ldrot = timeM2Ldrot + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! shift along z-axis in positive direction
  call M2Lz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

#ifdef TIMER
  timeM2Ldztr = timeM2Ldztr + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! Compute Wigner d-matrix for angle -beta
#if defined (KSCA) || defined (MKSCA)
  call dmatktilde(-beta)
#endif
#if defined (MSCA)
  call dmatmtilde(-beta)
#endif

#ifdef TIMER
  timeM2Ldmat = timeM2Ldmat + TimerCheck(50)
#endif

  ! backward rotation of irregular multipole expansion, L
#ifdef TIMER
  call TimerInit(50)
#endif

#if defined (KSTD) || defined (MSTD)
  ! reuse previously computed conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,m,0) * d(fptr(L,m,0)) * tmpExpansion(aL)
      csum = coefflkm(Lsum + m) * d(fptr(L,m,0)) * tmpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,m,k) * ( &
        csum = csum + coefflkm(Lsum + k*(L+1) + m) * ( &
               tmpExpansion(aL+k)  * d(fptr(L, m, k)) + &
        dconjg(tmpExpansion(aL+k)) * d(fptr(L, m,-k)) * mSign(mod(k,2)))
      enddo
      outExpansion(aL+m) = csum * expalp(m)
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use tilde-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + tmpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(tmpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum * expalp(m)
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Ldrot = timeM2Ldrot + TimerCheck(50)
#endif
end subroutine M2Lrotate2


subroutine M2Lrotate3(tvector, expalp, inpExpansion, tmpExpansion, scrExpansion, outExpansion)
  !- Rotate the coordinate system of regular multipole expansion to be translated
  !    to align its z-axis with translation axis
  !- Transform regular multipole expansion into irregular one via translation along the z-axis
  !- Rotate the coordinate system of produced irregular multipole expansion 
  !    to put its z-axis back to the original position
  !- Accumulate outExpansion
  !
#ifdef TIMER
  use timer
#endif
  use common_data, only: nL, tri, nSzTri, d, coefflkm
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)              ! translation vector
  double complex, intent(inout)   :: expalp(0:nL)            ! temporary array for complex exponents
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(inout)   :: scrExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  !
  ! Local variables
  double precision :: alpha, beta      ! Euler angles
  integer          :: L, k, m          ! indices of d-matrix d(L,k,m)
  integer          :: aL, Lsum
  double complex      :: csum
  double precision :: kalp, distance
  !
  ! Function declarations
  integer          :: fptr
!#if defined (KSTD) || defined (MSTD)
!  double precision :: Coeff
!#endif

  ! Determine Euler angles
  call GetEulerAngles3(tvector, alpha, beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )

  ! Compute Wigner d-matrix for angle beta
#ifdef TIMER
  call TimerInit(50)
#endif
#ifdef KSTD
  call dmatk(beta)
#endif
#ifdef MSTD
  call dmatm(beta)
#endif
#ifdef KSCA
  call dmatkbar(beta)
#endif
#ifdef MSCA
  call dmatmbar(beta)
#endif
#ifdef MKSCA
  call dmatmbar(beta)
#endif
#ifdef TIMER
  timeM2Ldmat = timeM2Ldmat + TimerCheck(50)
#endif
 
  ! rotate multipole expansion
  ! forward rotation of regular multipole expansion, R*
#ifdef TIMER
  call TimerInit(50)
#endif

  ! compute complex exponent
  kalp = 0.0d0
  expalp(0) = dcmplx(1.0d0, 0.0d0)
  do k = 1, nL
    kalp = kalp + alpha
    expalp(k) = dcmplx(dcos(kalp), dsin(kalp))
  enddo

#if defined (KSTD) || defined (MSTD)
  ! use conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,0,m) * d(fptr(L,0,m)) * inpExpansion(aL)
      csum = coefflkm(Lsum + m*(L+1)) * d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,k,m) * ( &
        csum = csum + coefflkm(Lsum + m*(L+1) + k) * ( &
                      expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2)) )
      enddo
      scrExpansion(aL+m) = csum
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use bar-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        csum = csum + expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      scrExpansion(aL+m) = csum
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Ldrot = timeM2Ldrot + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! shift along z-axis in positive direction
  call M2Lz(scrExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

#ifdef TIMER
  timeM2Ldztr = timeM2Ldztr + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! Compute Wigner d-matrix for angle -beta
#if defined (KSCA) || defined (MKSCA)
  call dmatktilde(-beta)
#endif
#if defined (MSCA)
  call dmatmtilde(-beta)
#endif

#ifdef TIMER
  timeM2Ldmat = timeM2Ldmat + TimerCheck(50)
#endif

  ! backward rotation of irregular multipole expansion, L
#ifdef TIMER
  call TimerInit(50)
#endif

#if defined (KSTD) || defined (MSTD)
  ! reuse previously computed conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,m,0) * d(fptr(L,m,0)) * tmpExpansion(aL)
      csum = coefflkm(Lsum + m) * d(fptr(L,m,0)) * tmpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,m,k) * ( &
        csum = csum + coefflkm(Lsum + k*(L+1) + m) * ( &
               tmpExpansion(aL+k)  * d(fptr(L, m, k)) + &
        dconjg(tmpExpansion(aL+k)) * d(fptr(L, m,-k)) * mSign(mod(k,2)))
      enddo
      outExpansion(aL+m) = outExpansion(aL+m) + csum * expalp(m)
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use tilde-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + tmpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(tmpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = outExpansion(aL+m) + csum * expalp(m)
    enddo
  enddo
#endif

#ifdef TIMER
  timeM2Ldrot = timeM2Ldrot + TimerCheck(50)
#endif
end subroutine M2Lrotate3


subroutine L2Lrotate(tvector, inpExpansion, tmpExpansion, outExpansion, alpha, beta)
  !- Rotate coordinate system of irregular multipole expansion to align it with z-axis
  !- Shift the origin of irregular multipole expansion along z-axis
  !- Rotate back the coordinate system of resulting irregular multipole expansion
  !
  use common_data, only: nL, tri, nSzTri, d
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(in)      :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  double precision, intent(in) :: alpha, beta             ! Euler angles
  !
  ! Local variables
  integer     :: L, m, k
  integer     :: aL
  double complex :: csum
  double precision :: kalp, malp, distance
  !
  ! Function declarations
  integer          :: fptr
  !double precision :: Coeff
  !
  ! Compute Wigner d-matrix for angle beta
  !call dmatm(beta)
  !call dmatk(beta)
  call dmatktilde(beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )
 
  ! rotate multipole expansion
  !
  ! forward rotation of irregular multipole expansion, S
  ! use conventional d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,m,0) * d(fptr(L,0,m)) * inpExpansion(aL)
  !    kalp = 0.0d0
  !    do k = 1, L
  !      kalp = kalp + alpha
  !      csum = csum + Coeff(L,m,k) * ( &
  !             dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
  !             dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
  !    enddo
  !    outExpansion(aL+m) = csum
  !  enddo
  !enddo
  ! use tilde-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      kalp = 0.0d0
      do k = 1, L
        kalp = kalp + alpha
        csum = csum + ( &
               dcmplx(dcos(kalp),-dsin(kalp)) * d(fptr(L, k, m)) * inpExpansion(aL+k) + &
               dcmplx(dcos(kalp), dsin(kalp)) * d(fptr(L,-k, m)) * dconjg(inpExpansion(aL+k)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo

  ! shift along z-axis in positive direction
  call L2Lz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

  ! Compute Wigner d-matrix for angle -beta
  call dmatktilde(-beta)

  ! backward rotation of irregular multipole expansion, S
  ! reuse previously computed conventional d(beta)
  !do L = 0, nL
  !  aL = tri(L)
  !  !
  !  do m = 0, L
  !    ! k = 0
  !    csum = Coeff(L,m,0) * d(fptr(L,m,0)) * tmpExpansion(aL)
  !    do k = 1, L
  !      csum = csum + Coeff(L,m,k) * ( &
  !             d(fptr(L, m, k)) * tmpExpansion(aL+k) + &
  !             d(fptr(L, m,-k)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
  !    enddo
  !    malp = alpha * dble(m)
  !    outExpansion(aL+m) = csum * dcmplx(dcos(malp), dsin(malp))
  !  enddo
  !enddo
  ! use tilde-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + ( &
               d(fptr(L, k, m)) * tmpExpansion(aL+k) + &
               d(fptr(L,-k, m)) * dconjg(tmpExpansion(aL+k)) * mSign(mod(k,2)))
      enddo
      malp = alpha * dble(m)
      outExpansion(aL+m) = csum * dcmplx(dcos(malp), dsin(malp))
    enddo
  enddo
end subroutine L2Lrotate


subroutine L2Lrotate2(tvector, expalp, inpExpansion, tmpExpansion, outExpansion)
  !- Rotate the coordinate system of irregular multipole expansion to be translated
  !    to align its z-axis with translation axis
  !- Translate the irregular multipole expansion along the z-axis
  !- Rotate the coordinate system of produced irregular multipole expansion 
  !    to put its z-axis back to the original position
  !
#ifdef TIMER
  use timer
#endif
  !
  use common_data, only: nL, tri, nSzTri, d, coefflkm
  use common_constants, only: mSign
  implicit none
  !
#ifdef MPI
  include 'mpif.h'
  integer :: mpiError
#endif
  !
  ! Function arguments
  double precision, intent(in) :: tvector(3)              ! translation vector
  double complex, intent(inout)   :: expalp(0:nL)            ! temporary array for complex exponents
  double complex, intent(in)      :: inpExpansion(0:nSzTri)  ! expansion to be rotated
  double complex, intent(inout)   :: tmpExpansion(0:nSzTri)  ! temporary expansion
  double complex, intent(out)     :: outExpansion(0:nSzTri)  ! rotated expansion
  !
  ! Local variables
  double precision :: alpha, beta      ! Euler angles
  integer          :: L, k, m          ! indices of d-matrix d(L,k,m)
  integer          :: aL, Lsum
  double complex   :: csum
  double precision :: kalp, distance
  !
  ! Function declarations
  integer          :: fptr
!#if defined (KSTD) || defined (MSTD)
!  double precision :: Coeff
!#endif

  ! Determine Euler angles
  call GetEulerAngles3(tvector, alpha, beta)

  ! Compute translation distance
  distance = dsqrt( tvector(1)*tvector(1) + tvector(2)*tvector(2) + tvector(3)*tvector(3) )

  ! Compute Wigner d-matrix for angle beta
#ifdef TIMER
  call TimerInit(50)
#endif
#ifdef KSTD
  call dmatk(beta)
#endif
#ifdef MSTD
  call dmatm(beta)
#endif
#if defined (KSCA) || defined (MKSCA)
  call dmatktilde(beta)  ! for L2L, KSCA == MKSCA will enforce k-tilde
#endif
#if defined (MSCA)
  call dmatmtilde(beta)  ! use MSCA to enforce m-tilde
#endif
#ifdef TIMER
  timeL2Ldmat = timeL2Ldmat + TimerCheck(50)
#endif

  ! rotate multipole expansion
  ! forward rotation of irregular multipole expansion
#ifdef TIMER
  call TimerInit(50)
#endif

  ! compute complex exponent
  kalp = 0.0d0
  expalp(0) = dcmplx(1.0d0, 0.0d0)
  do k = 1, nL
    kalp = kalp + alpha
    expalp(k) = dcmplx(dcos(kalp), -dsin(kalp))
  enddo

#if defined (KSTD) || defined (MSTD)
  ! use conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,m,0) * d(fptr(L,0,m)) * inpExpansion(aL)
      csum = coefflkm(Lsum + m) * d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,m,k) * ( &
        csum = csum + coefflkm(Lsum + k*(L+1) + m) * ( &
                      expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2)) )
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use tilde-scaled d(beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * inpExpansion(aL)
      do k = 1, L
        csum = csum + expalp(k) * inpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(expalp(k) * inpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum
    enddo
  enddo
#endif

#ifdef TIMER
  timeL2Ldrot = timeL2Ldrot + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! shift along z-axis in positive direction
  call L2Lz(outExpansion, -distance, tmpExpansion)  ! always moves in positive direction due to rotation

#ifdef TIMER
  timeL2Ldztr = timeL2Ldztr + TimerCheck(50)
#endif

#ifdef TIMER
  call TimerInit(50)
#endif

  ! Compute Wigner d-matrix for angle -beta
#if defined (KSCA) || defined (MKSCA)
  call dmatktilde(-beta)
#endif
#if defined (MSCA)
  call dmatmtilde(-beta)
#endif

#ifdef TIMER
  timeL2Ldmat = timeL2Ldmat + TimerCheck(50)
#endif

  ! backward rotation of irregular multipole expansion, L
#ifdef TIMER
  call TimerInit(50)
#endif

#if defined (KSTD) || defined (MSTD)
  ! reuse previously computed conventional d(beta)
  do L = 0, nL
    aL = tri(L)
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    !
    do m = 0, L
      ! k = 0
      !csum = Coeff(L,m,0) * d(fptr(L,m,0)) * tmpExpansion(aL)
      csum = coefflkm(Lsum + m) * d(fptr(L,m,0)) * tmpExpansion(aL)
      do k = 1, L
        !csum = csum + Coeff(L,m,k) * ( &
        csum = csum + coefflkm(Lsum + k*(L+1) + m) * ( &
               tmpExpansion(aL+k)  * d(fptr(L, m, k)) + &
        dconjg(tmpExpansion(aL+k)) * d(fptr(L, m,-k)) * mSign(mod(k,2)))
      enddo
      outExpansion(aL+m) = csum * dconjg(expalp(m))
    enddo
  enddo
#endif

#if defined (KSCA) || defined (MSCA) || defined (MKSCA)
  ! use tilde-scaled d(-beta)
  do L = 0, nL
    aL = tri(L)
    !
    do m = 0, L
      ! k = 0
      csum = d(fptr(L,0,m)) * tmpExpansion(aL)
      do k = 1, L
        csum = csum + tmpExpansion(aL+k)  * d(fptr(L, k, m)) + &
               dconjg(tmpExpansion(aL+k)) * d(fptr(L,-k, m)) * mSign(mod(k,2))
      enddo
      outExpansion(aL+m) = csum * dconjg(expalp(m))
    enddo
  enddo
#endif

#ifdef TIMER
  timeL2Ldrot = timeL2Ldrot + TimerCheck(50)
#endif
end subroutine L2Lrotate2


subroutine dmatm(beta)
  ! Compute Wigner d-matrix d(beta) with increment in index m
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k   !, lkm, lkm0, ikm, ikm0
  double precision :: cosb, sinb, sinb2, sinbsq2, opcosb2, omcosb2
  double precision :: Lk1
  !
  ! Function declarations
  integer          :: fptr
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.5d0
  sinbsq2 = sinb / dsqrt(2.0d0)
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0
  d(1) = 1.0d0

  ! L = 1
  d(2) = opcosb2         ! k= 1  m= 1
  d(3) = sinbsq2         ! k= 0  m= 1
  d(4) = omcosb2         ! k=-1  m= 1
  d(5) = -sinbsq2        ! k= 1  m= 0
  d(6) = cosb            ! k= 0  m= 0
  d(7) = sinbsq2         ! k=-1  m= 0
  d(8) = omcosb2         ! k= 1  m= 1
  d(9) = -sinbsq2        ! k= 0  m= 1
  d(10)= opcosb2         ! k=-1  m= 1

  ! recursively compute d-matrix
  do L = 2, nL
    !lkm0 = 1 + L * (2*L - 1) * (2*L + 1) / 3          ! index of first element of d-matrix of degree L
    !ikm0 = 1 + (L - 1) * (2*L - 3) * (2*L - 1) / 3    ! index of first element of d-matrix of degree L-1
    !
    do k = L, L-1, -1
      Lk1  = 1.0d0 / dble((L+k)*(L+k-1))
      !ikm = ikm0 + (L - k) * (2*L - 1)      ! add shift due to index k-1; ikm = fptr(L-1,k-1,L-1)
      !lkm = lkm0 + (L - k) * (2*L + 1)      ! add shift due to index k;   lkm = fptr(L,k,L)

      !fptr = 1 + L * (2*L - 1) * (2*L + 1) / 3  &  ! index of first element of matrix d(L,k,m)
      !     + (L - k) * (2*L + 1)                &  ! add shift due to index k
      !     + (L - m)                               ! add shift due to index m

      ! m=L
      !d(lkm) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * opcosb2 * d(ikm)                       !m-1  L-1,k-1,L-1
      d(fptr(L,k,L)) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * opcosb2 * d(fptr(L-1,k-1,L-1))  !m-1  L-1,k-1,L-1

      ! m=L-1
      !d(lkm+1) = dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * opcosb2 * d(ikm+1) &       !m-1  L-1,k-1,L-2
      !                   - dsqrt( dble(2*L-1) * Lk1 ) * sinb    * d(ikm)           !m    L-1,k-1,L-1
      d(fptr(L,k,L-1)) = dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * opcosb2 * d(fptr(L-1,k-1,L-2)) &  !m-1  L-1,k-1,L-2
                         - dsqrt( dble(2*L-1) * Lk1 ) * sinb    * d(fptr(L-1,k-1,L-1))            !m    L-1,k-1,L-1

      ! m= L-2 ... -L+2
      do m = L-2, -L+2, -1
        !d(lkm+L-m) = dsqrt( dble((L+m)*(L+m-1)) * Lk1 ) * opcosb2 * d(ikm+L-m) &   !m-1  L-1,k-1,m-1
        !           - dsqrt( dble((L+m)*(L-m))   * Lk1 ) * sinb    * d(ikm+L-1-m) & !m    L-1,k-1,m
        !           + dsqrt( dble((L-m)*(L-m-1)) * Lk1 ) * omcosb2 * d(ikm+L-2-m)   !m+1  L-1,k-1,m+1
        d(fptr(L,k,m)) = dsqrt( dble((L+m)*(L+m-1)) * Lk1 ) * opcosb2 * d(fptr(L-1,k-1,m-1)) &   !m-1  L-1,k-1,m-1
                   - dsqrt( dble((L+m)*(L-m))   * Lk1 ) * sinb    * d(fptr(L-1,k-1,m)) &         !m    L-1,k-1,m
                   + dsqrt( dble((L-m)*(L-m-1)) * Lk1 ) * omcosb2 * d(fptr(L-1,k-1,m+1))   !m+1  L-1,k-1,m+1
      enddo

      ! m= -L+1
      !d(lkm+L+L-1) = -dsqrt( dble(2*L-1)     * Lk1 ) * sinb    * d(ikm+L-1+L-1) &  !m    L-1,k-1,-L+1
      !        + dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * omcosb2 * d(ikm+L-1+L-2)    !m+1  L-1,k-1,-L+2
      d(fptr(L,k,-L+1)) = -dsqrt( dble(2*L-1)     * Lk1 ) * sinb    * d(fptr(L-1,k-1,-L+1)) &  !m    L-1,k-1,-L+1
              + dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * omcosb2 * d(fptr(L-1,k-1,-L+2))    !m+1  L-1,k-1,-L+2

      ! m= -L
      !d(lkm+L+L) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * omcosb2 * d(ikm+L-1+L-1)          !m+1  L-1,k-1,-L+1
      d(fptr(L,k,-L)) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * omcosb2 * d(fptr(L-1,k-1,-L+1))     !m+1  L-1,k-1,-L+1
    enddo

    do k = L-2, -L, -1
      Lk1  = 1.0d0 / dble((L-k)*(L-k-1))
      ! m=L
      !ikm = fptr(L-1,k+1,L-1)
      !lkm = fptr(L,k,L)
      !d(lkm) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * omcosb2 * d(ikm)                      !m    L-1,k+1,L-1
      d(fptr(L,k,L)) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * omcosb2 * d(fptr(L-1,k+1,L-1)) !m    L-1,k+1,L-1

      ! m=L-1
      !d(lkm+1) = dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * omcosb2 * d(ikm+1) &       !m-1  L-1,k+1,L-2
      !         + dsqrt( dble(2*L-1)           * Lk1 ) * sinb    * d(ikm)           !m    L-1,k+1,L-1
      d(fptr(L,k,L-1)) = dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * omcosb2 * d(fptr(L-1,k+1,L-2)) &   !m-1  L-1,k+1,L-2
               + dsqrt( dble(2*L-1)           * Lk1 ) * sinb    * d(fptr(L-1,k+1,L-1))             !m    L-1,k+1,L-1

      ! m= L-2 ... -L+2
      do m = L-2, -L+2, -1
        !d(lkm+L-m) = dsqrt( dble((L+m)*(L+m-1)) * Lk1 ) * omcosb2 * d(ikm+L-m) &   !m-1  L-1,k+1,m-1
        !           + dsqrt( dble((L+m)*(L-m))   * Lk1 ) * sinb    * d(ikm+L-1-m) & !m    L-1,k+1,m
        !           + dsqrt( dble((L-m)*(L-m-1)) * Lk1 ) * opcosb2 * d(ikm+L-2-m)   !m=1  L-1,k+1,m+1
        d(fptr(L,k,m)) = dsqrt( dble((L+m)*(L+m-1)) * Lk1 ) * omcosb2 * d(fptr(L-1,k+1,m-1)) &   !m-1  L-1,k+1,m-1
                   + dsqrt( dble((L+m)*(L-m))   * Lk1 ) * sinb    * d(fptr(L-1,k+1,m)) &   !m    L-1,k+1,m
                   + dsqrt( dble((L-m)*(L-m-1)) * Lk1 ) * opcosb2 * d(fptr(L-1,k+1,m+1))   !m=1  L-1,k+1,m+1
      enddo

      ! m= -L+1
      !d(lkm+L+L-1) = dsqrt( dble(2*L-1)           * Lk1 ) * sinb    * d(ikm+L-1+L-1) & !m    L-1,k+1,-L+1
      !             + dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * opcosb2 * d(ikm+L-1+L-2)   !m+1  L-1,k+1,-L+2
      d(fptr(L,k,-L+1)) = dsqrt( dble(2*L-1) * Lk1 ) * sinb * d(fptr(L-1,k+1,-L+1)) &   !m    L-1,k+1,-L+1
                   + dsqrt( dble((2*L-1)*(2*L-2)) * Lk1 ) * opcosb2 * d(fptr(L-1,k+1,-L+2))  !m+1  L-1,k+1,-L+2
      ! m= -L
      !d(lkm+L+L) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * opcosb2 * d(ikm+L-1+L-1)            !m    L-1,k+1,-L+1
      d(fptr(L,k,-L)) = dsqrt( dble(2*L*(2*L-1)) * Lk1 ) * opcosb2 * d(fptr(L-1,k+1,-L+1)) !m    L-1,k+1,-L+1
    enddo
    !
    ! Validate computed d-matrix
    !call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .FALSE., .TRUE.)
  enddo
  !
end subroutine dmatm


subroutine dmatmbar(beta)
  ! Compute bar-scaled Wigner d-matrix with increment in index m
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2
  !
  ! Function declarations
  integer          :: fptr
  !double precision :: Coeff
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.50d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0         Eq. 12.109
  d(1) = 1.0d0

  ! L = 1         Eq. 12.111
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb2                   ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb                   ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb                    ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb2                  ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! recursively compute d-matrix
  do L = 2, nL
    !
    do k = L, L-1, -1
      ! m=L              Eq. 12.131  m-set
      d(fptr(L,k,L)) = opcosb2 * d(fptr(L-1,k-1,L-1))

      ! m=L-1            Eq. 12.133  m-set
      d(fptr(L,k,L-1)) = opcosb2 * d(fptr(L-1,k-1,L-2)) - sinb * d(fptr(L-1,k-1,L-1))

      ! m= L-2 ... -L+2  Eq. 12.51   m-set
      do m = L-2, -L+2, -1
        d(fptr(L,k,m)) = opcosb2 * d(fptr(L-1,k-1,m-1)) - sinb * d(fptr(L-1,k-1,m)) + omcosb2 * d(fptr(L-1,k-1,m+1))
      enddo

      ! m= -L+1          Eq. 12.135  m-set
      d(fptr(L,k,-L+1)) = -sinb * d(fptr(L-1,k-1,-L+1)) + omcosb2 * d(fptr(L-1,k-1,-L+2))

      ! m= -L            Eq. 12.137  m-set
      d(fptr(L,k,-L)) = omcosb2 * d(fptr(L-1,k-1,-L+1))
    enddo

    do k = L-2, -L, -1
      ! m=L              Eq. 12.132  m-set
      d(fptr(L,k,L)) = omcosb2 * d(fptr(L-1,k+1,L-1))

      ! m=L-1            Eq. 12.134  m-set
      d(fptr(L,k,L-1)) = omcosb2 * d(fptr(L-1,k+1,L-2)) + sinb * d(fptr(L-1,k+1,L-1))

      ! m= L-2 ... -L+2  Eq. 12.53   m-set
      do m = L-2, -L+2, -1
        d(fptr(L,k,m)) = omcosb2 * d(fptr(L-1,k+1,m-1)) + sinb * d(fptr(L-1,k+1,m)) + opcosb2 * d(fptr(L-1,k+1,m+1))
      enddo

      ! m= -L+1          Eq. 12.136  m-set
      d(fptr(L,k,-L+1)) = sinb * d(fptr(L-1,k+1,-L+1)) + opcosb2 * d(fptr(L-1,k+1,-L+2))

      ! m= -L            Eq. 12.138
      d(fptr(L,k,-L)) = opcosb2 * d(fptr(L-1,k+1,-L+1))
    enddo
    !
  enddo

  ! Validation test
  !do L = 1, nL
  !  ! unscale the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,m,k)  ! unscale
  !    enddo
  !  enddo
  !
  !  ! run validation test
  !  call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .false., .true.)
  !
  !  ! scale back the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,k,m)  ! scale back
  !    enddo
  !  enddo
  !enddo
  
end subroutine dmatmbar


subroutine dmatmtilde(beta)
  ! Compute tilde-scaled Wigner d-matrix with increment in index m
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Subroutine arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k, L2, L2m1, L2m2
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2, LpKm1LpK, LmKm1LmK
  !
  ! Function declarations
  integer          :: fptr
  !double precision :: Coeff
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.50d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0         Eq. 12.109
  d(1) = 1.0d0

  ! L = 1         Eq. 12.111
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb                    ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb2                  ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb2                   ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb                   ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! recursively compute d-matrix
  do L = 2, nL
    L2   = L+L
    L2m1 = L+L - 1
    L2m2 = L+L - 2
    !
    do k = L, L-1, -1
      LpKm1LpK = 1.0d0 / dble((L+k-1)*(L+k))

      ! m=L              Eq. 12.155 m-set
      d(fptr(L,k,L)) = opcosb2 * d(fptr(L-1,k-1,L-1)) * dble(L2 * L2m1) * LpKm1LpK

      ! m=L-1            Eq. 12.157 m-set
      d(fptr(L,k,L-1)) = ( opcosb2 * d(fptr(L-1,k-1,L-2)) * dble(L2m1 * L2m2) &
                         - sinb * d(fptr(L-1,k-1,L-1)) * dble(L2m1) ) * LpKm1LpK

      ! m= L-2 ... -L+2  Eq. 12.69  m-set
      do m = L-2, -L+2, -1
        d(fptr(L,k,m)) = ( opcosb2 * d(fptr(L-1,k-1,m-1)) * dble((L+m-1)*(L+m)) &
                         - sinb * d(fptr(L-1,k-1,m)) * dble((L-m)*(L+m)) &
                         + omcosb2 * d(fptr(L-1,k-1,m+1)) * dble((L-m-1)*(L-m)) ) * LpKm1LpK
      enddo

      ! m= -L+1          Eq. 12.159 m-set
      d(fptr(L,k,-L+1)) = ( -sinb * d(fptr(L-1,k-1,-L+1)) * dble(L2m1) & 
                          + omcosb2 * d(fptr(L-1,k-1,-L+2)) * dble(L2m1 * L2m2) ) * LpKm1LpK

      ! m= -L            Eq. 12.161 m-set
      d(fptr(L,k,-L)) = omcosb2 * d(fptr(L-1,k-1,-L+1)) * dble(L2m1 * L2) * LpKm1LpK
    enddo

    do k = L-2, -L, -1
      LmKm1LmK = 1.0d0 / dble((L-k-1)*(L-k))
      ! m=L              Eq. 12.156 m-set
      d(fptr(L,k,L)) = omcosb2 * d(fptr(L-1,k+1,L-1)) * dble(L2 * L2m1) * LmKm1LmK

      ! m=L-1            Eq. 12.158 m-set
      d(fptr(L,k,L-1)) = ( omcosb2 * d(fptr(L-1,k+1,L-2)) * dble(L2m1 * L2m2) &
                         + sinb * d(fptr(L-1,k+1,L-1)) * dble(L2m1) ) * LmKm1LmK

      ! m= L-2 ... -L+2  Eq. 12.71  m-set
      do m = L-2, -L+2, -1
        d(fptr(L,k,m)) = ( omcosb2 * d(fptr(L-1,k+1,m-1)) * dble((L+m-1)*(L+m)) &
                         + sinb * d(fptr(L-1,k+1,m)) * dble((L-m)*(L+m)) &
                         + opcosb2 * d(fptr(L-1,k+1,m+1)) * dble((L-m-1)*(L-m)) ) * LmKm1LmK
      enddo

      ! m= -L+1          Eq. 12.160 m-set
      d(fptr(L,k,-L+1)) = ( sinb * d(fptr(L-1,k+1,-L+1)) * dble(L2m1) &
                          + opcosb2 * d(fptr(L-1,k+1,-L+2)) * dble(L2m1 * L2m2) ) * LmKm1LmK

      ! m= -L            Eq. 12.162 m-set
      d(fptr(L,k,-L)) = opcosb2 * d(fptr(L-1,k+1,-L+1)) * dble(L2 * L2m1) * LmKm1LmK
    enddo
    !
  enddo

  ! Validation test
  !do L = 1, nL
  !  ! unscale the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,k,m)  ! unscale
  !    enddo
  !  enddo
  !
  !  ! run validation test
  !  call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .false., .true.)
  !
  !  ! scale back the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,m,k)  ! scale back
  !    enddo
  !  enddo
  !enddo
  
end subroutine dmatmtilde


subroutine dmatk(beta)
  ! Compute Wigner d-matrix d(beta) with increment in index k (d_km)
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k
  double precision :: cosb, sinb, sinb2, sinbsq2, opcosb2, omcosb2
  !
  ! Function declarations
  integer          :: fptr
  !
  ! initialize d-matrix
  cosb    = dcos(beta)
  sinb    = dsin(beta)
  sinb2   = sinb * 0.5d0
  sinbsq2 = sinb / dsqrt(2.0d0)
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0
  d(1) = 1.0d0

  ! L = 1
  d(2) = opcosb2         ! k= 1  m= 1
  d(3) = sinbsq2         ! k= 0  m= 1
  d(4) = omcosb2         ! k=-1  m= 1
  d(5) = -sinbsq2        ! k= 1  m= 0
  d(6) = cosb            ! k= 0  m= 0
  d(7) = sinbsq2         ! k=-1  m= 0
  d(8) = omcosb2         ! k= 1  m= 1
  d(9) = -sinbsq2        ! k= 0  m= 1
  d(10)= opcosb2         ! k=-1  m= 1
  !
  ! recursively compute d-matrix
  do L = 2, nL

    ! m = L;  k = L
    d(fptr(L,L,L)) = opcosb2 * d(fptr(L-1,L-1,L-1))

    ! m = L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      d(fptr(L,k,L)) = dsqrt( dble(2*L*(2*L-1)) / dble((L+k)*(L-k)) ) * sinb2 * d(fptr(L-1,k,L-1)) 
    enddo

    ! m = L;  k = -L
    d(fptr(L,-L,L)) = omcosb2 * d(fptr(L-1,-L+1,L-1))

    ! L-1 >= m >= -L+1
    do m = L-1, -L+1, -1
      ! k = L
      d(fptr(L,L,m)) = -dsqrt( dble(2*L-1) / dble(2*L*(l-m)*(L+m)) ) * dble(L) * sinb * d(fptr(L-1,L-1,m))

      ! k = L-1
      d(fptr(L,L-1,m)) = dsqrt( dble(2*L-1) / dble((L-m)*(L+m)) ) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))

      ! L-2 >= k >= -L+2
      do k = L-2, -L+2, -1
        d(fptr(L,k,m)) = 1.0d0 / dsqrt( dble((L-k)*(L-m)*(L+m)) ) * &
                         ( dsqrt(dble(L-k-1)) * dble(L) * sinb * d(fptr(L-1,k+1,m)) + &
                           dsqrt(dble(L+k)) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) )
      enddo

      ! k = -L+1
      d(fptr(L,-L+1,m)) = dsqrt(dble(2*L-1) / dble((L-m)*(L+m))) * (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))

      ! k = -L
      d(fptr(L,-L,m)) = dsqrt(dble(2*L-1) / dble(2*L*(L-m)*(L+m))) * dble(L) * sinb * d(fptr(L-1,-L+1,m))
    enddo

    ! m = -L;  k = L
    d(fptr(L,L,-L)) = omcosb2 * d(fptr(L-1,L-1,-L+1))

    ! m = -L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      d(fptr(L,k,-L)) = -dsqrt(dble(2*L*(2*L-1)) / dble((L+k)*(L-k))) * sinb2 * d(fptr(L-1,k,-L+1))
    enddo

    ! m = -L;  k = -L
    d(fptr(L,-L,-L)) = opcosb2 * d(fptr(L-1,-L+1,-L+1))

  enddo
  !
  ! Validation test
  !do L = 1, nL
  !  call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .false., .true.)
  !enddo
  !
end subroutine dmatk


subroutine dmatkbar(beta)
  ! Compute bar-scaled Wigner d-matrix d(beta) with increment in index k (d_km)
  ! Naive implementation
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2
  !double precision :: dtmp
  !double precision :: Coeff
  !
  ! Function declarations
  integer          :: fptr
 
  ! initialize trigonometric functions
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.5d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0
  d(1) = 1.0d0

  ! L = 1
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb2                   ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb                   ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb                    ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb2                  ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! d(L,k,m) index computation
  ! fptr = 1
  ! fptr = fptr + L * (2*L - 1) * (2*L + 1) / 3  ! size of the matrices up to L - 1
  ! fptr = fptr + (L - m) * (2*L + 1)            ! add shift due to index m
  ! fptr = fptr + (L - k)                        ! add shift due to index k

  ! recursively compute d-matrix
  do L = 2, nL

    ! m = L;  k = L
    !dtmp = opcosb2 * d(fptr(L-1,L-1,L-1))
    !call dcmpbar(dtmp,L,L,L,"vam1")
    d(fptr(L,L,L)) = opcosb2 * d(fptr(L-1,L-1,L-1))

    ! m = L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      !dtmp = sinb2 * d(fptr(L-1,k,L-1)) 
      !call dcmpbar(dtmp,L,k,L,"vam2")
      d(fptr(L,k,L)) = sinb2 * d(fptr(L-1,k,L-1)) 
    enddo

    ! m = L;  k = -L
    !dtmp = omcosb2 * d(fptr(L-1,-L+1,L-1))
    !call dcmpbar(dtmp,L,-L,L,"vam3")
    d(fptr(L,-L,L)) = omcosb2 * d(fptr(L-1,-L+1,L-1))

    ! L-1 >= m >= -L+1
    do m = L-1, -L+1, -1
      ! k = L
      !dtmp = -dble(2*L-1) / dble((l-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,L-1,m))
      !call dcmpbar(dtmp,L,L,m,"vam4")
      d(fptr(L,L,m)) = -dble(2*L-1) / dble((l-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,L-1,m))

      ! k = L-1
      !dtmp = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))
      !call dcmpbar(dtmp,L,L-1,m,"vam4")
      d(fptr(L,L-1,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))

      ! L-2 >= k >= -L+2
      do k = L-2, -L+2, -1
        !dtmp = ( dble(L-k-1) * dble(L) * sinb * d(fptr(L-1,k+1,m)) + &
        !         dble(L+k) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble((L-m)*(L+m))
        !call dcmpbar(dtmp,L,k,m,"vam5")
        ! Eq. 12.67
        d(fptr(L,k,m)) = ( dble((L-k-1) * L) * sinb * d(fptr(L-1,k+1,m)) + &
                 dble(L+k) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble((L-m)*(L+m))
      enddo

      ! k = -L+1
      !dtmp = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))
      !call dcmpbar(dtmp,L,-L+1,m,"vam6")
      d(fptr(L,-L+1,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))

      ! k = -L
      !dtmp = dble(2*L-1) / dble((L-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,-L+1,m))
      !call dcmpbar(dtmp,L,-L,m,"vam7")
      d(fptr(L,-L,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,-L+1,m))
    enddo

    ! m = -L;  k = L
    !dtmp = omcosb2 * d(fptr(L-1,L-1,-L+1))
    !call dcmpbar(dtmp,L,L,-L,"vam8")
    d(fptr(L,L,-L)) = omcosb2 * d(fptr(L-1,L-1,-L+1))

    ! m = -L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      !dtmp = -sinb2 * d(fptr(L-1,k,-L+1))
      !call dcmpbar(dtmp,L,k,-L,"vam9")
      d(fptr(L,k,-L)) = -sinb2 * d(fptr(L-1,k,-L+1))
    enddo

    ! m = -L;  k = -L
    !dtmp = opcosb2 * d(fptr(L-1,-L+1,-L+1))
    !call dcmpbar(dtmp,L,-L,-L,"vam0")
    d(fptr(L,-L,-L)) = opcosb2 * d(fptr(L-1,-L+1,-L+1))

  enddo

  ! Validation test
  !do L = 1, nL
  !  ! unscale the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,m,k)  ! unscale
  !    enddo
  !  enddo
  !
  !  ! run validation test
  !  call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .false., .true.)
  !
  !  ! scale back the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,k,m)  ! scale back
  !    enddo
  !  enddo
  !enddo
 
end subroutine dmatkbar


subroutine dmatkbarOpt(beta)
  ! Compute bar-scaled Wigner d-matrix d(beta) with recursion in index k (d_km)
  ! Optimized implementation
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer          :: L, m, k  ! indices of d-matrix d(L,k,m)
  integer          :: i,  j    ! index of elements d(L,L,L) and d(L-1,L-1,L-1)
  integer          :: im, jm   ! index of elements d(L,L,m) and d(L-1,L-1,m)
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2
  !
  ! Function declarations
  !integer          :: fptr
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.5d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0
  d(1) = 1.0d0

  ! L = 1
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb2                   ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb                   ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb                    ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb2                  ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! d(L,k,m) index computation
  ! fptr = 1
  ! fptr = fptr + L * (2*L - 1) * (2*L + 1) / 3  ! size of the matrices up to L - 1
  ! fptr = fptr + (L - m) * (2*L + 1)            ! add shift due to index m
  ! fptr = fptr + (L - k)                        ! add shift due to index k
 
  ! recursively compute d-matrix
  do L = 2, nL
    ! m = L
    i =     L * (2*L - 1) * (2*L + 1) / 3 + 1   ! index of first element d(L,L,L)
    j = (L-1) * (2*L - 3) * (2*L - 1) / 3 + 1   ! index of first element d(L-1,L-1,L-1)

    ! m = L;  k = L
    d(i) = opcosb2 * d(j)             ! d(fptr(L,L,L)) = opcosb2 * d(fptr(L-1,L-1,L-1))

    ! m = L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      d(i+L-k) = sinb2 * d(j+L-k-1)   ! d(fptr(L,k,L)) = sinb2 * d(fptr(L-1,k,L-1)) 
    enddo

    ! m = L;  k = -L
    d(i+L+L) = omcosb2 * d(j+L+L-2)   ! d(fptr(L,-L,L)) = omcosb2 * d(fptr(L-1,-L+1,L-1))

    ! L-1 >= m >= -L+1
    do m = L-1, -L+1, -1
      im = i + (L   - m) * (2*L + 1)    ! add shift due to index m for d(L,L,m)
      jm = j + (L-1 - m) * (2*L - 1)    ! add shift due to index m for d(L-1,L-1,m)

      ! k = L
      !d(fptr(L,L,m)) = -dble(2*L-1) / dble((l-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,L-1,m))
      d(im) = -dble(2*L-1) / dble((l-m)*(L+m)) * dble(L) * sinb * d(jm)

      ! k = L-1
      !d(fptr(L,L-1,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))
      d(im+1) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb - dble(m)) * d(jm)

      ! L-2 >= k >= -L+2
      do k = L-2, -L+2, -1
        !d(fptr(L,k,m)) = ( dble(L-k-1) * dble(L) * sinb * d(fptr(L-1,k+1,m)) + &
        !                   dble(L+k) * (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble((L-m)*(L+m))
        d(im+L-k) = ( dble(L-k-1) * dble(L) * sinb * d(jm+L-k-2) + &
                      dble(L+k) * (dble(L)*cosb - dble(m)) * d(jm+L-k-1) ) / dble((L-m)*(L+m))
      enddo

      ! k = -L+1
      !d(fptr(L,-L+1,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))
      d(im+L+L-1) = dble(2*L-1) / dble((L-m)*(L+m)) * (dble(L)*cosb + dble(m)) * d(jm+L+L-2)

      ! k = -L
      !d(fptr(L,-L,m)) = dble(2*L-1) / dble((L-m)*(L+m)) * dble(L) * sinb * d(fptr(L-1,-L+1,m))
      d(im+L+L) = dble(2*L-1) / dble((L-m)*(L+m)) * dble(L) * sinb * d(jm+L+L-2)
    enddo

    ! m = -L
    im = i + (L   + L) * (2*L + 1)      ! add shift due to index m for d(L,L,-L)
    jm = j + (L+L - 2) * (2*L - 1)      ! add shift due to index m for d(L-1,L-1,-L+1)

    ! m = -L;  k = L
    d(im) = omcosb2 * d(jm)             ! d(fptr(L,L,-L)) = omcosb2 * d(fptr(L-1,L-1,-L+1))

    ! m = -L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      d(im+L-k) = -sinb2 * d(jm+L-k-1)  ! d(fptr(L,k,-L)) = -sinb2 * d(fptr(L-1,k,-L+1))
    enddo

    ! m = -L;  k = -L
    d(im+L+L) = opcosb2 * d(jm+L+L-2)   ! d(fptr(L,-L,-L)) = opcosb2 * d(fptr(L-1,-L+1,-L+1))

  enddo
  !
end subroutine dmatkbarOpt


subroutine dmatktilde(beta)
  ! Compute tilde-scaled Wigner d-matrix with increment in index k
  ! Naive implementation
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer     :: L, m, k
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2
  !double precision :: dtmp
  !
  ! Function declarations
  integer          :: fptr
  !double precision :: Coeff
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.5d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0      Eq. 12.109
  d(1) = 1.0d0

  ! L = 1      Eq. 12.111
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb                    ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb2                  ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb2                   ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb                   ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! recursively compute d-matrix
  do L = 2, nL

    ! m = L;  k = L               Eq. 12.72  m-set
    !dtmp = opcosb2 * d(fptr(L-1,L-1,L-1))
    !call dcmptilde(dtmp,L,L,L,"vam1")
    d(fptr(L,L,L)) = opcosb2 * d(fptr(L-1,L-1,L-1))

    ! m = L;  L-1 >= k >= -L+1    Eq. 12.73  m-set
    do k = L-1, -L+1, -1
      !dtmp = dble(L*(2*L-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,L-1)) 
      !call dcmptilde(dtmp,L,k,L,"vam2")
      d(fptr(L,k,L)) = dble(L*(2*L-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,L-1)) 
    enddo

    ! m = L;  k = -L              Eq. 12.74  m-set
    !dtmp = omcosb2 * d(fptr(L-1,-L+1,L-1))
    !call dcmptilde(dtmp,L,-L,L,"vam3")
    d(fptr(L,-L,L)) = omcosb2 * d(fptr(L-1,-L+1,L-1))

    ! L-1 >= m >= -L+1
    do m = L-1, -L+1, -1
      ! k = L                     Eq. 12.102 k-set
      !dtmp = -sinb2 * d(fptr(L-1,L-1,m))
      !call dcmptilde(dtmp,L,L,m,"vam4")
      d(fptr(L,L,m)) = -sinb2 * d(fptr(L-1,L-1,m))

      ! k = L-1                   Eq. 12.103 k-set
      !dtmp = (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))
      !call dcmptilde(dtmp,L,L-1,m,"vam4")
      d(fptr(L,L-1,m)) = (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))

      ! L-2 >= k >= -L+2          Eq. 12.101 k-set
      do k = L-2, -L+2, -1
        !dtmp = ( dble(L)*sinb * d(fptr(L-1,k+1,m)) + &
        !        (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble(L-k) 
        !call dcmptilde(dtmp,L,k,m,"vam5")
        d(fptr(L,k,m)) = ( dble(L)*sinb * d(fptr(L-1,k+1,m)) + &
                (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble(L-k) 
      enddo

      ! k = -L+1                  Eq. 12.104 k-set
      !dtmp = (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))
      !call dcmptilde(dtmp,L,-L+1,m,"vam6")
      d(fptr(L,-L+1,m)) = (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))

      ! k = -L                    Eq. 12.105 k-set
      !dtmp = sinb2 * d(fptr(L-1,-L+1,m))
      !call dcmptilde(dtmp,L,-L,m,"vam7")
      d(fptr(L,-L,m)) = sinb2 * d(fptr(L-1,-L+1,m))
    enddo

    ! m = -L;  k = L              Eq. 12.75  m-set
    !dtmp = omcosb2 * d(fptr(L-1,L-1,-L+1))
    !call dcmptilde(dtmp,L,L,-L,"vam8")
    d(fptr(L,L,-L)) = omcosb2 * d(fptr(L-1,L-1,-L+1))

    ! m = -L;  L-1 >= k >= -L+1   Eq. 12.76
    do k = L-1, -L+1, -1
      !dtmp = -dble(L*(2*L-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,-L+1))
      !call dcmptilde(dtmp,L,k,-L,"vam9")
      d(fptr(L,k,-L)) = -dble(L*(2*L-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,-L+1))
    enddo

    ! m = -L;  k = -L             Eq. 12.77
    !dtmp = opcosb2 * d(fptr(L-1,-L+1,-L+1))
    !call dcmptilde(dtmp,L,-L,-L,"vam0")
    d(fptr(L,-L,-L)) = opcosb2 * d(fptr(L-1,-L+1,-L+1))

  enddo
  !
  !! Validation test for tilde-scaled Wigner d-matrix
  !! loop over d-matrices
  !do L = 1, nL
  !  ! unscale the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,k,m)  ! unscale
  !    enddo
  !  enddo
  ! 
  !  ! run validation test
  !  call checkdmat(d(fptr(L,L,L)), L, 2*L+1, .false., .true.)
  ! 
  !  ! scale back the d-matrix
  !  do m = L, -L, -1
  !    do k = L, -L, -1
  !      d(fptr(L,k,m)) = d(fptr(L,k,m)) * Coeff(L,m,k)  ! scale back
  !    enddo
  !  enddo
  !enddo
  !
end subroutine dmatktilde


subroutine dmatktildeOpt(beta)
  ! Compute tilde-scaled Wigner d-matrix d(beta) with recursion in index k (d_km)
  ! Optimized implementation
  !
  use common_data, only: nL, d
  implicit none
  !
  ! Function arguments
  double precision, intent(in) :: beta    ! Euler angle
  !
  ! Local variables
  integer          :: L, m, k  ! indices of d-matrix d(L,k,m)
  integer          :: i,  j    ! index of elements d(L,L,L) and d(L-1,L-1,L-1)
  integer          :: im, jm   ! index of elements d(L,L,m) and d(L-1,L-1,m)
  double precision :: cosb, sinb, sinb2, opcosb2, omcosb2
  integer          :: L2
  double precision :: dL, dLL2m1, dm
  !
  ! Function declarations
  !integer          :: fptr
  !
  ! initialize d-matrix
  cosb  = dcos(beta)
  sinb  = dsin(beta)
  sinb2 = sinb * 0.5d0
  opcosb2 = (1.0d0 + cosb) * 0.5d0
  omcosb2 = (1.0d0 - cosb) * 0.5d0

  ! zero all elements
  d = 0.0d0

  ! L = 0
  d(1) = 1.0d0

  ! L = 1
  d(2) = opcosb2                 ! k= 1  m= 1
  d(3) = sinb                    ! k= 0  m= 1
  d(4) = omcosb2                 ! k=-1  m= 1
  d(5) = -sinb2                  ! k= 1  m= 0
  d(6) = cosb                    ! k= 0  m= 0
  d(7) = sinb2                   ! k=-1  m= 0
  d(8) = omcosb2                 ! k= 1  m=-1
  d(9) = -sinb                   ! k= 0  m=-1
  d(10)= opcosb2                 ! k=-1  m=-1

  ! d(L,k,m) index computation
  ! fptr = 1
  ! fptr = fptr + L * (2*L - 1) * (2*L + 1) / 3  ! size of the matrices up to L - 1
  ! fptr = fptr + (L - m) * (2*L + 1)            ! add shift due to index m
  ! fptr = fptr + (L - k)                        ! add shift due to index k

  ! recursively compute d-matrix
  do L = 2, nL
    L2 = L + L
    dLL2m1 = dble(L*(L2-1))
    dL = dble(L)

    ! m = L
    i =     L * (L2 - 1) * (L2 + 1) / 3 + 1   ! index of first element d(L,L,L)
    j = (L-1) * (L2 - 3) * (L2 - 1) / 3 + 1   ! index of first element d(L-1,L-1,L-1)

    ! m = L;  k = L
    d(i) = opcosb2 * d(j)             ! d(fptr(L,L,L)) = opcosb2 * d(fptr(L-1,L-1,L-1))

    ! m = L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      !d(fptr(L,k,L)) = dble(L*(L2-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,L-1)) 
      d(i+L-k) = dLL2m1 / dble((L-k)*(L+k)) * sinb * d(j+L-k-1) 
    enddo

    ! m = L;  k = -L
    d(i+L2) = omcosb2 * d(j+L2-2)     ! d(fptr(L,-L,L)) = omcosb2 * d(fptr(L-1,-L+1,L-1))

    ! L-1 >= m >= -L+1
    do m = L-1, -L+1, -1
      im = i + (L   - m) * (L2 + 1)   ! add shift due to index m for d(L,L,m)
      jm = j + (L-1 - m) * (L2 - 1)   ! add shift due to index m for d(L-1,L-1,m)
      dm = dble(m)

      ! k = L
      d(im) = -sinb2 * d(jm)          ! d(fptr(L,L,m)) = -sinb2 * d(fptr(L-1,L-1,m))

      ! k = L-1
      ! d(fptr(L,L-1,m)) = (dble(L)*cosb - dble(m)) * d(fptr(L-1,L-1,m))
      d(im+1) = (dL*cosb - dm) * d(jm)

      ! L-2 >= k >= -L+2
      do k = L-2, -L+2, -1
        !d(fptr(L,k,m)) = ( dble(L)*sinb * d(fptr(L-1,k+1,m)) + &
        !        (dble(L)*cosb - dble(m)) * d(fptr(L-1,k,m)) ) / dble(L-k) 
        d(im+L-k) = ( dL*sinb * d(jm+L-k-2) + &
                (dL*cosb - dm) * d(jm+L-k-1) ) / dble(L-k) 
      enddo

      ! k = -L+1
      !d(fptr(L,-L+1,m)) = (dble(L)*cosb + dble(m)) * d(fptr(L-1,-L+1,m))
      d(im+L2-1) = (dL*cosb + dm) * d(jm+L2-2)

      ! k = -L
      d(im+L2) = sinb2 * d(jm+L2-2)   ! d(fptr(L,-L,m)) = sinb2 * d(fptr(L-1,-L+1,m))
    enddo

    ! m = -L
    im = i + L2       * (L2 + 1)      ! add shift due to index m for d(L,L,-L)
    jm = j + (L2 - 2) * (L2 - 1)      ! add shift due to index m for d(L-1,L-1,-L+1)

    ! m = -L;  k = L
    d(im) = omcosb2 * d(jm)           ! d(fptr(L,L,-L)) = omcosb2 * d(fptr(L-1,L-1,-L+1))

    ! m = -L;  L-1 >= k >= -L+1
    do k = L-1, -L+1, -1
      !d(fptr(L,k,-L)) = -dble(L*(2*L-1)) / dble((L-k)*(L+k)) * sinb * d(fptr(L-1,k,-L+1))
      d(im+L-k) = -dLL2m1 / dble((L-k)*(L+k)) * sinb * d(jm+L-k-1)
    enddo

    ! m = -L;  k = -L
    d(im+L2) = opcosb2 * d(jm+L2-2)   ! d(fptr(L,-L,-L)) = opcosb2 * d(fptr(L-1,-L+1,-L+1))

  enddo
  !
end subroutine dmatktildeOpt


subroutine dcmp(dtmp,L,k,m,message)
  ! Compare Wigner d-matrix elements
  !
  use common_data, only: d
  implicit none
  double precision, intent(in)   :: dtmp
  integer, intent(in)            :: L, k, m
  character (len=4), intent(in) :: message
  !
  ! Function declarations
  integer :: fptr
  !
  if(abs(dtmp - d(fptr(L,k,m))) > 1.0d-10) then
    write(*,'(a,a,a,3i3)') "error in ", message, " L,k,m = ",L,k,m
    write(*,'(a,f16.10)')  "dtmp = ", dtmp
    write(*,'(a,f16.10)')  "d    = ", d(fptr(L,k,m))
    error stop
  endif
end subroutine dcmp


subroutine dcmpbar(dtmp,L,k,m,message)
  ! Compare bar-scaled Wigner d-matrix elements
  !
  use common_data, only: d
  implicit none
  double precision, intent(in)   :: dtmp
  integer, intent(in)            :: L, k, m
  character (len=4), intent(in) :: message
  !
  ! Function declarations
  integer          :: fptr
  double precision :: Coeff
  !
  if(abs(dtmp / Coeff(L,k,m) - d(fptr(L,k,m))) > 1.0d-10) then
    write(*,'(a,a,a,3i3)') "error in ", message, " L,k,m = ",L,k,m
    write(*,'(a,f16.10)')  "dtmp = ", dtmp
    write(*,'(a,f16.10)')  "d    = ", d(fptr(L,k,m))
    error stop
  endif
end subroutine dcmpbar


subroutine dcmptilde(dtmp,L,k,m,message)
  ! Compare tilde-scaled Wigner d-matrix elements
  !
  use common_data, only: d
  implicit none
  double precision, intent(in)   :: dtmp
  integer, intent(in)            :: L, k, m
  character (len=4), intent(in) :: message
  !
  ! Function declarations
  integer          :: fptr
  double precision :: Coeff
  !
  if(abs(dtmp / Coeff(L,m,k) - d(fptr(L,k,m))) > 1.0d-10) then
    write(*,'(a,a,a,3i3)') "error in ", message, " L,k,m = ",L,k,m
    write(*,'(a,f16.10)')  "dtmp = ", dtmp
    write(*,'(a,f16.10)')  "d    = ", d(fptr(L,k,m))
    error stop
  endif
end subroutine dcmptilde


subroutine checkdmat(d, L, n, doprint, dovalidate)
  !
  ! Validate Wigner d-matrix
  ! Map the linear array d() to two-dimensional matrix d(k,m)
  !
  implicit none
  integer,          intent(in) :: L,n
  double precision, intent(in) :: d(n,n)
  logical,          intent(in) :: doprint, dovalidate
  integer                      :: m, k, i, j, a, b
  double precision, parameter  :: zero = 1.0d-10
  double precision             :: sum
  !
  if(doprint) then
    ! Print the d-matrix
    do k = 1, n
      write(*,'(20f12.8)') (d(k,m), m = 1, n)
    enddo
  endif
  !
  if(dovalidate) then
    !
    ! row normalization
    do i = 1, n
      k = L - i + 1
      sum = 0.0d0
      do j = 1, n
        sum = sum + d(i,j)**2
      enddo
      if(dabs(sum - 1.0d0) > zero) then
        write(*,'(a,2i3,f16.8)') "Normalization test failed for L, k = ", L, k, sum
        error stop
      endif
    enddo
    !
    ! column normalization
    do j = 1, n
      m = L - j + 1
      sum = 0.0d0
      do i = 1, n
        sum = sum + d(i,j)**2
      enddo
      if(dabs(sum - 1.0d0) > zero) then
        write(*,'(a,2i3,f16.8)') "Normalization test failed for L, m = ", L, m, sum
        error stop
      endif
    enddo
    !
    ! row orthogonalization
    do i = 1, n-1
      k = L - i + 1
      sum = 0.0d0
      do a = i+1, n
        b = L - a + 1
        do j = 1, n
          sum = sum + d(i,j)*d(a,j)
        enddo
        if(dabs(sum) > zero) then
          write(*,'(a,3i3,f16.8)') "Orthogonalization test failed for L, k1, k2 = ", L, k, b, sum
          error stop
        endif
      enddo
    enddo
    !
    ! column orthogonalization
    do j = 1, n-1
      m = L - j + 1
      sum = 0.0d0
      do b = j+1, n
        a = L - b + 1
        do i = 1, n
          sum = sum + d(i,j)*d(i,b)
        enddo
        if(dabs(sum) > zero) then
          write(*,'(a,3i3,f16.8)') "Orthogonalization test failed for L, m1, m2 = ", L, m, a, sum
          error stop
        endif
      enddo
    enddo
    !
    ! symmetry across the main diagonal: d(k,m) = d(m,k) *(-1)^(m-k)
    do j = 2, n
      m = L - j + 1
      do i = 1, j-1
        k = L - i + 1
        if(dabs(d(i,j) - d(j,i) * (-1.0d0)**(j-i)) > zero) then
          write(*,'(a,2f8.2,4x,2i3)') &
          "Symmetry test failed for d(k,m) d(m,k) k m =", d(i,j), d(j,i), k, m
          error stop
        endif
      enddo
    enddo
    !
    ! symmetry in the index sign: d(k,m) = d(-m,-k)
    do m = L, 0, -1
      j = L - m + 1   ! index of +m in 2D array d(i,j)
      b = L + m + 1   ! index of -m in 2D array d(i,b); m self is always positive
      do k = L, -L, -1
        i = L - k + 1 ! index of  k in 2D array d(i,j)
        a = 2*L+2 - i ! index of -k in 2D array d(a,j); k may be positive or negative
        if(dabs(d(i,j) - d(b,a)) > zero) then
          write(*,'(a,2f8.2,4x,2i3,4x,2i3)') &
          "Symmetry test failed for d(k,m) d(-m,-k)  k m -m -k", d(i,j), d(b,a), k, m, -m, -k
          error stop
        endif
      enddo
    enddo
    !
    write(*,'(a,i3)') "Validation test passed for d(L), L = ", L
  endif
  
end subroutine checkdmat


pure integer function fptr(L,k,m)
  ! compute pointer to element d(L,k,m) of Wigner d-matrix(k,m) of degree L
  ! k is row index:    L, L-1, ..., 0, ..., -L+1, -L
  ! m is column index: L, L-1, ..., 0, ..., -L+1, -L
  ! d-matrix assumes column-major order, where index k runs faster than index m
  implicit none
  integer, intent(in) :: L, k, m

  ! sanity check
!  if(abs(k) > L .or. abs(m) > L) then
!    write(*,'(a,3i4)') "Error in indices L, k, m: ", L, k, m
!    error stop
!  endif

  ! determine the pointer to d(L,L,L) as a sum of preceding elements plus 1
  fptr = L * (2*L - 1) * (2*L + 1) / 3 + 1

  ! add shift due to index m
  fptr = fptr + (L - m) * (2*L + 1)

  ! add shift due to index k
  fptr = fptr + (L - k)
end function fptr


!integer function fptr(L,k,m)
!  ! works with conventional d-matrix
!  ! compute pointer to element (L,k,m) of rectangular Wigner d-matrix(k,m) of degree L
!  ! k is row index: L, L-1, ..., 0, ..., -L+1, -L
!  ! m is row element index: L, L-1, ..., 0, ..., -L+1, -L
!  ! d-matrix is stored in row-major order, that is index m in d(L,k,m) runs faster than index k
!  implicit none
!  integer, intent(in) :: L, k, m
!
!  ! sanity check
!  if(abs(k) > L .or. abs(m) > L) then
!    write(*,'(a,3i4)') "Error in indices L, k, m: ", L, k, m
!    error stop
!  endif
!
!  fptr = 1
!  if(L == 0) return
!
!  ! determine cumulative number of elements in preceeding matrices up to L - 1
!  fptr = fptr + L * (2*L - 1) * (2*L + 1) / 3
!
!  ! add shift due to index k
!  fptr = fptr + (L - k) * (2*L + 1)
!
!  ! add shift due to index m
!  fptr = fptr + (L - m)
!end function fptr


pure double precision function Coeff(L,k,m)
  !
  ! k = |k| and m = |m|
  ! For k > m:
  ! compute SQRT[ {(L-k)! (L+k)!} / {(L-m)! (L+m)!} ] = SQRT[ (L+m+1)...(L+k) / (L-k+1)...(L-m) ]
  ! For k < m:
  ! compute SQRT[ {(L-k)! (L+k)!} / {(L-m)! (L+m)!} ] = SQRT[ (L-m+1)...(L-k) / (L+k+1)...(L+m) ]
  !
  implicit none
  integer, intent(in) :: L, k, m
  !
  ! Local variables
  integer :: i, ka, ma, ks, ms
  !double precision :: a
  !
  ka = iabs(k)                  ! k =  |k|
  ma = iabs(m)                  ! m =  |m|
  ks = isign(1, ma-ka) * ka     ! k = -|k|   if |k| > |m|
  ms = isign(1, ka-ma) * ma     ! m = -|m|   if |m| > |k|
  !
  Coeff = 1.0d0
  do i = 1, max0(ka,ma) - min0(ka,ma)   ! if k == m, the do-loop skips
    Coeff = Coeff * dble(L+ms+i) / dble(L+ks+i)
  enddo
  !
  !a = 1.0d0
  !if(ka > ma) then
  !  do i = L+ma+1, L+ka
  !    a = a * dble(i)
  !  enddo
  !  do i = L-ka+1, L-ma
  !    a = a / dble(i)
  !  enddo
  !else
  !  do i = L-ma+1, L-ka
  !    a = a * dble(i)
  !  enddo
  !  do i = L+ka+1, L+ma
  !    a = a / dble(i)
  !  enddo
  !endif
  !
  !if(abs(a - Coeff) > 1.0d-12) then
  !  write(*,'(a,f18.12)') "a = ", a
  !  write(*,'(a,f18.12)') "c = ", Coeff
  !  write(*,'(a,3i4)') "L,k,m = ", L,k,m
  !  error stop
  !endif
  !
  Coeff = dsqrt(Coeff)
  !
end function Coeff


subroutine initsq(iam)
  !
  ! k = |k| and m = |m|
  ! For k > m:
  ! compute SQRT[ {(L-k)! (L+k)!} / {(L-m)! (L+m)!} ] = SQRT[ (L+m+1)...(L+k) / (L-k+1)...(L-m) ]
  ! For k < m:
  ! compute SQRT[ {(L-k)! (L+k)!} / {(L-m)! (L+m)!} ] = SQRT[ (L-m+1)...(L-k) / (L+k+1)...(L+m) ]
  !
  use timer
  use common_data, only: nL, memUsed, memMax, coefflkm, ncoeff
  implicit none
  ! Function arguments
  integer, intent(in) :: iam   ! MPI rank
  !
  ! Local variables
  integer :: i, L, ka, ma, ks, ms, allocStat, Lsum
  double precision :: memLocal, coeff
  !
#ifdef TIMER
  call TimerInit(2)
#endif
  !
  ! determine the number of coefficients
  ncoeff = (nL + 1) * ( 2 * nL*nL + 7 * nL + 6 ) / 6
  write(*,'(a,i0)') "Number of precomputed coeffs:  ", ncoeff

  ! allocate arrays
  allocate (coefflkm(0:ncoeff-1), stat=allocStat)
  if(allocStat /= 0) call exitPrg(2)
  memLocal = dble( ncoeff*8 )
  memMax = dmax1(memMax, (memUsed + memLocal))

  coefflkm(0) = 1.0d0
  do L = 1, nL
    Lsum = L * ( 2 * (L-1)*(L-1) + 7 * (L-1) + 6 ) / 6  ! sum of elements for 0 ... (L-1)
    do ma = 0, L
      do ka = 0, L
        !ka = iabs(k)                 ! k =  |k|
        !ma = iabs(m)                 ! m =  |m|
        ks = isign(1, ma-ka) * ka     ! k = -|k|   if |k| > |m|
        ms = isign(1, ka-ma) * ma     ! m = -|m|   if |m| > |k|
        !
        coeff = 1.0d0
        do i = 1, max0(ka,ma) - min0(ka,ma)   ! if k == m, the do-loop skips
          coeff = coeff * dble(L+ms+i) / dble(L+ks+i)
        enddo
        coefflkm(Lsum + ma * (L+1) + ka) = dsqrt(coeff)
        !
      enddo
    enddo
  enddo
  !
#ifdef TIMER
  if(iam == 0) then
    write(*,'(a,f9.1,a)') "Timer: Precomputing coeff:", TimerCheck(2), "s"
  endif
#endif
  !
end subroutine initsq
