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

subroutine exitPrg(stat)

#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif

  implicit none
  integer, intent(in) :: stat
  integer             :: iam

#ifdef MPI
  integer             :: mpiError
  call MPI_Comm_rank(MPI_COMM_WORLD, iam, mpiError)
#else
  iam = 0
#endif

  select case (stat)
  case(0)
    write(*,'(/a)') "FMM Error. Input data missing"
  case(1)
    write(*,'(/a)') "FMM Error. Too many atoms"
  case(2)
    write(*,'(/a)') "FMM Error. Cannot allocate memory"
  case(3)
    write(*,'(/a)') "FMM Error. Insufficiently deep hierarchy level"
  case(4)
    write(*,'(/a)') "FMM Error. Atom assignement to cubes is faulty"
  case(5)
    write(*,'(/a)') "FMM Error. Non-zero total charge"
  case(6)
    write(*,'(/a)') "FMM Error. MPI Broadcast"
  case(7)
    write(*,'(/a)') "FMM Error. fmmCubeEdge is undefined"
  case(8)
    write(*,'(/a)') "FMM Error. MPI Barrier"
  case(9)
    write(*,'(/a)') "FMM Error. MPI AllReduce"
  case(10)
    write(*,'(/a)') "FMM Error. Internal Error in LatticeSumsPPinNF"
  case(11)
    write(*,'(/a)') "FMM Error. Wrong input parameter"
  case(99)
    write(*,'(/a)') "FMM Error. Stopped on request"
  case default
    write(*,'(/a)') "FMM Program Error"
  end select

#ifdef MPI
  write(*,'(a,i0)') "Reported by rank ", iam
  call MPI_Abort(MPI_COMM_WORLD, -1, mpiError)
  call MPI_Finalize(mpiError)
#endif
  error stop 1
end subroutine exitPrg


subroutine Usage(prgName)
  implicit none
  character (len=*), intent(in) :: prgName
  write(*,'(/a,a,a)') "Usage: ", trim(prgName), &
  " fmmCubeEdge polynomialOrder ws input.psf input.pdb scale=0|1 nRepLayers [[[[prtRank=0] szSubComm=1] nSplitRanks=12] ddLevel]"
end subroutine Usage


subroutine NotExist(fname)
  implicit none
  character (len=*), intent(in) :: fname
  write(*,'(a,a,a)') "File '", trim(fname), "' doesn't exist"
end subroutine NotExist
