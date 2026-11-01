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

module common_constants
  implicit none

  ! Fundamental constants:
  ! vacuum permittivitty eps0 = 8.854187817 * 10^-12 (A^2 * s^4) / (kg * m^3)
  ! Coulomb's constant k = 1/(4*PI*eps0)  [Newton * meter^2 / Coulomb^2]
  ! Avogadro number Nav = 6.02214129 * 10^23 (mol^-1)
  ! 1 e = 1.602176565 * 10^-19 C [Coulomb]
  ! 1 Angstrom = 10^-8 cm = 10^-10 m
  ! 1 cal = 4.1868 J
  !
  ! V = J/C = (N*m)/(A*s) = (kg * m^2) / (A * s^3)
  ! J = N * m;  C = A * s;  N = (kg * m) / s^2  [A - ampere; s - second]
  !
  ! eConv:
  ! Electrostatic energy: E = 331.84164001 * (q*q)/r [kcal/mol]; q in electron units; r in Angstroms
  ! E = (q * q)/(4*PI*eps0 * r) [Joule]
  ! (10^-19 * 10^-19) / (10^-12 * 10^-10) * 10^23 * 10^-3 = 10000 (kJ/mol)
  ! 331.84164001 = 1.602176565*1.602176565/4/3.141592653589793/8.854187817*6.02214129/4.1868*10000
  !
  ! pConv:
  ! Electrostatic potential: V = 1/(4*PI*eps0) * (q/r)
  ! V = 14.399644851 * q/r [Volt]; q in electron units; r in Angstroms
  ! 10^-19 / (10^-12 * 10^-10) = 1000
  ! 14.399644851 = 1.602176565/4/3.141592653589793/8.854187817*1000
  !
  ! fConv:
  ! Electric force: F = 1/(4*PI*eps0) * (q^2/r^2) = [ (N * m^2 / C^2) * (C^2 / m^2) ] = [N]
  ! F= (q * q) / (4*PI*eps0 * r * r)
  ! F = 2.307077352532*10^-8 * q^2/r^2 [Newton]; q in electron units; r in Angstroms
  ! (10^-19 * 10^-19) / (10^-12 * 10^-10 * 10^-10) * 10^9 = 1000 nano-Newtons
  ! 23.07077352532266825 = 1.602176565*1.602176565/4/3.141592653589793/8.854187817*1000

  double precision            :: eConv = 331.84164001d0    ! energy in kcal/mol: q in electrons, r in Angstroms
  double precision            :: pConv = 14.399644851d0    ! potential in Volts: q in electrons, r in Angstroms
  double precision            :: fConv = 23.07077352532266825d0 ! force in nano-Newtons: q in electrons, r in Angstroms
  double precision, parameter :: pi = 3.14159265358979323846d0
  double precision, parameter :: mSign(0:1) = (/ 1.0d0, -1.0d0 /)
 !double precision, parameter :: machineZero = 1.0d-11
  double precision, parameter :: machineZero = 1.0d-5
  double precision, parameter :: iMB = 1.0d0 / 1024.0d0 / 1024.0d0            ! convert bytes to megabytes
  double precision, parameter :: zeroCoor = 1.0d-8, zero = 1.0d-14
  double precision, parameter :: iGB = 1.0d0 / 1024.0d0 / 1024.0d0 / 1024.0d0 ! convert bytes to gigabytes

  integer, parameter :: TYPE_FILE = 0, TYPE_INT = 1, TYPE_DBLE = 2
  integer (kind=8), parameter :: izero = 0, ione = 1
  double precision, parameter :: invLog8 = 1.0d0 / dlog(8.0d0)

end module common_constants

