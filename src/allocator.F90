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

module allocator
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int32_t, c_int64_t, c_f_pointer, c_associated
  implicit none

  interface
    function getcpuid() bind(c, name='getcpuid') result(int_cpuid)
      import c_int32_t
      integer (kind=c_int32_t)        :: int_cpuid
    end function
  end interface

  interface
    function allocate_hbm_memory(size) bind(c, name='allocate_hbm_memory')
      import c_int64_t, c_ptr
      integer (kind=c_int64_t), value :: size
      type (c_ptr) :: allocate_hbm_memory
    end function
  end interface

  interface
    subroutine free_hbm_memory(hbm_ptr, size) bind(c, name='free_hbm_memory')
      import c_int64_t, c_ptr
      integer (kind=c_int64_t), value :: size
      type (c_ptr) :: hbm_ptr
    end subroutine
  end interface

  interface
    subroutine stat_numa_memory() bind(c, name='stat_numa_memory')
    end subroutine
  end interface

  interface
    function stat_memory(numa_node) bind(c, name='stat_memory') result(avail_mem)
      import c_int32_t, c_int64_t
      integer (kind=c_int32_t), value :: numa_node
      integer (kind=c_int64_t)        :: avail_mem
    end function
  end interface

  interface
    function stat_allmemory() bind(c, name='stat_allmemory') result(avail_allmem)
      import c_int32_t, c_int64_t
      integer (kind=c_int64_t)        :: avail_allmem
    end function
  end interface

contains

end module allocator
