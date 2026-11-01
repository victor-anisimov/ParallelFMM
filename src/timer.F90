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

! timer 1 accounts for total time to solution
! timer 2 accounts for time spent in a routine
! timer 3 accounts for time breakdown in the routine
module timer
  implicit none
  save
  public
  integer, parameter :: nTimers = 50
  integer(kind=8)    :: clock_count, clock_rate, clock_max
  integer(kind=8)    :: timerStart(nTimers)
  logical            :: first = .true.
  double precision   :: timeM2Mdmat, timeM2Mdztr, timeM2Mdrot
  double precision   :: timeM2Ldmat, timeM2Ldztr, timeM2Ldrot
  double precision   :: timeL2Ldmat, timeL2Ldztr, timeL2Ldrot
  double precision   :: timeConsumed, timeConsumedMin, timeConsumedMax, timeConsumedSum

contains


subroutine TimerInit( iTimer )
#ifdef MPI
#ifdef CPLR_NVFTN
  use mpi
#else
  use mpi_f08
#endif
#endif

  implicit none
  integer, intent(in) :: iTimer
  integer             :: i

  if( iTimer.gt.0 .and. iTimer.le.nTimers ) then
#ifdef MPI
    timerStart(iTimer) = MPI_WTime()
#else
    call system_clock(clock_count, clock_rate, clock_max)
    timerStart(iTimer) = clock_count
#endif
    if(first) then
      do i = 1, nTimers
        timerStart(i) = timerStart(iTimer)
      enddo
      first = .false.
    endif
  endif

end subroutine TimerInit


subroutine TimerRotReset()
  implicit none

  timeM2Mdmat = 0.0d0
  timeM2Mdztr = 0.0d0
  timeM2Mdrot = 0.0d0

  timeM2Ldmat = 0.0d0
  timeM2Ldztr = 0.0d0
  timeM2Ldrot = 0.0d0

  timeL2Ldmat = 0.0d0
  timeL2Ldztr = 0.0d0
  timeL2Ldrot = 0.0d0

end subroutine TimerRotReset


double precision function TimerCheck( iTimer )
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  integer, intent(in) :: iTimer

  if( iTimer.gt.0 .and. iTimer.le.nTimers ) then
#ifdef MPI
    TimerCheck = MPI_WTime() - timerStart(iTimer)
#else
    call system_clock(clock_count, clock_rate, clock_max)
    TimerCheck = dble(clock_count - timerStart(iTimer)) / dble(clock_rate)
#endif
  else
    ! it should not get here
    TimerCheck = 0.0d0
  endif

end function TimerCheck


end module timer
