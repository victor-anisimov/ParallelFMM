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

subroutine EstimateNumNearFarNeighbors(ws, loLevel, sumNearNeighbors, sumFarNeighbors)
  !
  ! return: the total number of near- and far-neighbors in a cubic cell
  ! algorithm: count the number of neighbors in 1 octant and multiply by 8 in 3D
  !
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif
  !
  use timer
  use common_data, only: maxNNperBox, maxFNperBox 
  implicit none
  !
  integer, intent(in)             :: ws, loLevel
  integer (kind = 8), intent(out) :: sumNearNeighbors, sumFarNeighbors
  !
  integer :: iam
  integer :: iLevel, i, j, k, nCubesEdge, nCubesScan
  integer :: nN, nF
  !
#ifdef MPI
  integer :: mpiError
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
#endif
  !
  sumFarNeighbors  = 0
  sumNearNeighbors = 0  ! there are zero neighbors in level 0 (master cube)
  nCubesEdge = 1
  nCubesScan = 1
  !
  ! loop over the hierarchy of cubes
  do iLevel = 1, loLevel
    nCubesEdge = nCubesEdge * 2         ! number of cubes along a side
    !
    ! loop over each cube i,j,k in the hierarchy
    ! determine parent box indices fx,fy,fz for each child box i,j,k
    !
    maxNNperBox = 0
    maxFNperBox = 0
    do k = 1, nCubesScan
      do j = 1, nCubesScan
        do i = 1, nCubesScan
          !
          ! add the number of boxes in the NF super-cube minus one (self)
          nN = (min0(k+ws,nCubesEdge) - max0(k-ws,1) + 1) * &
               (min0(j+ws,nCubesEdge) - max0(j-ws,1) + 1) * &
               (min0(i+ws,nCubesEdge) - max0(i-ws,1) + 1)
          !
          sumNearNeighbors = sumNearNeighbors + nN - 1
          maxNNperBox = max0(maxNNperBox, nN)  ! record the highest number of near neighbors for a box
          !
          nF = ( min0((int((k+1)/2)+ws)*2, nCubesEdge) - max0((int((k-1)/2)-ws)*2, 0) ) * &
               ( min0((int((j+1)/2)+ws)*2, nCubesEdge) - max0((int((j-1)/2)-ws)*2, 0) ) * &
               ( min0((int((i+1)/2)+ws)*2, nCubesEdge) - max0((int((i-1)/2)-ws)*2, 0) ) - nN
          sumFarNeighbors = sumFarNeighbors + nF
          maxFNperBox = max0(maxFNperBox, nF)  ! record the highest number of far neighbors for a box
          !
        enddo
      enddo
    enddo
    nCubesScan = nCubesScan * 2
    !
    ! subtract the reference box from the near-field super-box size
    maxNNperBox = max0(maxNNperBox - 1, 0)
    maxFNperBox = max0(maxFNperBox, 0)       ! does not include reference box; no need to subtract
    !
    ! print the stats for hierarchy level
    if(iam == 0) then
      write(*,'(a,i3,2(a,i16),2(a,i6))')             &
        "iLevel = ", iLevel,                         &
        "  sumNearNeighbors = ", sumNearNeighbors*8, &
        "  sumFarNeighbors = ", sumFarNeighbors*8,   &
        "  maxNNperBox = ", maxNNperBox, "  maxFNperBox = ", maxFNperBox
    endif
  enddo
  !
  ! multiply by 8 octants in 3D
  sumNearNeighbors = sumNearNeighbors * 8   ! 2 * 2 * 2 in 3D
  sumFarNeighbors  = sumFarNeighbors  * 8   ! 2 * 2 * 2 in 3D
  !
end subroutine EstimateNumNearFarNeighbors
