!> Test suite for ThermoModule (only Nose–Hoover chain components)
!!
!! This suite is intentionally RNG-free: it tests
!!   - thermo_NHC_local
!!   - thermo_NHC_global

module test_thermo
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt
   use ParticleModule
   use ThermoModule, only: thermo_NHC_local, thermo_NHC_global
   implicit none
   private

   public :: collect_thermo_suite

contains
   subroutine collect_thermo_suite(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)

      tests = [ &
         new_unittest("NHC local: dt=0 leaves state unchanged",        test_nhc_local_dt0), &
         new_unittest("NHC local: forward/backward reversibility",     test_nhc_local_reversible), &
         new_unittest("NHC global: dt=0 leaves state unchanged",       test_nhc_global_dt0), &
         new_unittest("NHC global: forward/backward reversibility",    test_nhc_global_reversible), &
         new_unittest("NHC global: bypass when any particle mass=0",   test_nhc_global_mass0_bypass) &
      ]
   end subroutine collect_thermo_suite


   subroutine test_nhc_local_dt0(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: M, l
      real(8) :: dt, p, p0, mass, beta, thermE, thermE_expected, tol
      real(8), allocatable :: mNHC(:), rNHC(:), pNHC(:), r0(:), pN0(:)

      tol  = 1.0d-12
      M    = 4
      dt   = 0.0d0
      beta = 0.7d0
      mass = 2.3d0
      p    = 1.234d0
      p0   = p

      allocate(mNHC(M), rNHC(M), pNHC(M), r0(M), pN0(M))

      ! deterministic initial thermostat state (no RNG)
      mNHC = [ 1.0d0, 1.7d0, 2.1d0, 3.3d0 ]
      rNHC = [ 0.10d0, -0.05d0, 0.03d0, 0.07d0 ]
      pNHC = [ 0.20d0, -0.10d0, 0.40d0, -0.30d0 ]

      r0  = rNHC
      pN0 = pNHC

      call thermo_NHC_local(dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      call check(error, abs(p - p0) <= tol, "p changed with dt=0 in thermo_NHC_local")
      call check(error, maxval(abs(rNHC - r0)) <= tol, "rNHC changed with dt=0 in thermo_NHC_local")
      call check(error, maxval(abs(pNHC - pN0)) <= tol, "pNHC changed with dt=0 in thermo_NHC_local")

      thermE_expected = 0.0d0
      do l = 1, M
         thermE_expected = thermE_expected + 0.5d0*pN0(l)*pN0(l)/mNHC(l) + r0(l)/beta
      end do

      call check(error, abs(thermE - thermE_expected) <= 1.0d-12, "thermE formula mismatch (local, dt=0)")

      deallocate(mNHC, rNHC, pNHC, r0, pN0)
   end subroutine test_nhc_local_dt0


   subroutine test_nhc_local_reversible(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: M
      real(8) :: dt, p, p_init, mass, beta, thermE, tol
      real(8), allocatable :: mNHC(:), rNHC(:), pNHC(:), r_init(:), pN_init(:)

      tol  = 5.0d-11
      M    = 4
      dt   = 5.0d-3
      beta = 0.9d0
      mass = 1.8d0
      p    = 0.77d0

      allocate(mNHC(M), rNHC(M), pNHC(M), r_init(M), pN_init(M))

      mNHC = [ 0.8d0, 1.1d0, 1.4d0, 2.0d0 ]
      rNHC = [ 0.02d0, -0.03d0, 0.01d0, 0.00d0 ]
      pNHC = [ 0.05d0, 0.02d0, -0.04d0, 0.01d0 ]

      p_init  = p
      r_init  = rNHC
      pN_init = pNHC

      ! Forward step
      call thermo_NHC_local( dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)
      ! Backward step (reversibility test)
      call thermo_NHC_local(-dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      call check(error, abs(p - p_init) <= tol, "local NHC not reversible: p did not return")
      call check(error, maxval(abs(rNHC - r_init)) <= tol, "local NHC not reversible: rNHC did not return")
      call check(error, maxval(abs(pNHC - pN_init)) <= tol, "local NHC not reversible: pNHC did not return")

      deallocate(mNHC, rNHC, pNHC, r_init, pN_init)
   end subroutine test_nhc_local_reversible


   subroutine test_nhc_global_dt0(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: M, ndim, natom, l, NTotDim
      real(8) :: dt, beta, thermE, thermE_expected, tol
      real(8), allocatable :: mNHC(:), rNHC(:), pNHC(:), r0(:), pN0(:)
      real(8) :: p(3,2), p0(3,2), mass(2)

      tol   = 1.0d-12
      M     = 4
      ndim  = 3
      natom = 2
      NTotDim = ndim*natom

      dt   = 0.0d0
      beta = 0.8d0

      mass = [ 12.0d0, 1.0d0 ]
      p(:,1) = [ 0.10d0, -0.20d0, 0.30d0 ]
      p(:,2) = [ -0.05d0, 0.04d0, -0.02d0 ]
      p0 = p

      allocate(mNHC(M), rNHC(M), pNHC(M), r0(M), pN0(M))
      mNHC = [ 1.0d0, 1.3d0, 1.9d0, 2.7d0 ]
      rNHC = [ 0.03d0, 0.01d0, -0.02d0, 0.04d0 ]
      pNHC = [ -0.10d0, 0.07d0, 0.02d0, -0.05d0 ]
      r0  = rNHC
      pN0 = pNHC

      call thermo_NHC_global(dt, M, ndim, natom, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      call check(error, maxval(abs(p - p0)) <= tol, "p changed with dt=0 in thermo_NHC_global")
      call check(error, maxval(abs(rNHC - r0)) <= tol, "rNHC changed with dt=0 in thermo_NHC_global")
      call check(error, maxval(abs(pNHC - pN0)) <= tol, "pNHC changed with dt=0 in thermo_NHC_global")

      thermE_expected = 0.0d0
      thermE_expected = thermE_expected + 0.5d0*pN0(1)*pN0(1)/mNHC(1) + dble(NTotDim)*r0(1)/beta
      do l = 2, M
         thermE_expected = thermE_expected + 0.5d0*pN0(l)*pN0(l)/mNHC(l) + r0(l)/beta
      end do

      call check(error, abs(thermE - thermE_expected) <= 1.0d-12, "thermE formula mismatch (global, dt=0)")

      deallocate(mNHC, rNHC, pNHC, r0, pN0)
   end subroutine test_nhc_global_dt0


   subroutine test_nhc_global_reversible(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: M, ndim, natom
      real(8) :: dt, beta, thermE, tol
      real(8) :: p(3,2), p_init(3,2), mass(2)
      real(8), allocatable :: mNHC(:), rNHC(:), pNHC(:), r_init(:), pN_init(:)

      tol   = 1.0d-10
      M     = 4
      ndim  = 3
      natom = 2

      dt   = 3.0d-3
      beta = 1.1d0

      mass = [ 14.0d0, 16.0d0 ]
      p(:,1) = [ 0.12d0, -0.08d0, 0.05d0 ]
      p(:,2) = [ -0.03d0, 0.07d0, -0.09d0 ]

      allocate(mNHC(M), rNHC(M), pNHC(M), r_init(M), pN_init(M))
      mNHC = [ 0.9d0, 1.2d0, 1.5d0, 2.4d0 ]
      rNHC = [ 0.01d0, -0.01d0, 0.02d0, 0.00d0 ]
      pNHC = [ 0.03d0, -0.02d0, 0.01d0, 0.02d0 ]

      p_init  = p
      r_init  = rNHC
      pN_init = pNHC

      call thermo_NHC_global( dt, M, ndim, natom, p, mass, beta, rNHC, pNHC, mNHC, thermE)
      call thermo_NHC_global(-dt, M, ndim, natom, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      call check(error, maxval(abs(p - p_init)) <= tol, "global NHC not reversible: p did not return")
      call check(error, maxval(abs(rNHC - r_init)) <= tol, "global NHC not reversible: rNHC did not return")
      call check(error, maxval(abs(pNHC - pN_init)) <= tol, "global NHC not reversible: pNHC did not return")

      deallocate(mNHC, rNHC, pNHC, r_init, pN_init)
   end subroutine test_nhc_global_reversible


   subroutine test_nhc_global_mass0_bypass(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: M, ndim, natom
      real(8) :: dt, beta, thermE, thermE0, tol
      real(8) :: p(3,2), p0(3,2), mass(2)
      real(8), allocatable :: mNHC(:), rNHC(:), pNHC(:)

      tol   = 0.0d0   ! should be exactly unchanged due to early RETURN
      M     = 3
      ndim  = 3
      natom = 2

      dt   = 1.0d-2
      beta = 1.0d0

      mass = [ 12.0d0, 0.0d0 ]   ! triggers bypass
      p(:,1) = [ 0.40d0, 0.10d0, -0.20d0 ]
      p(:,2) = [ 0.30d0, -0.50d0, 0.60d0 ]
      p0 = p

      allocate(mNHC(M), rNHC(M), pNHC(M))
      mNHC = [ 1.0d0, 1.0d0, 1.0d0 ]
      rNHC = [ 0.0d0, 0.0d0, 0.0d0 ]
      pNHC = [ 0.0d0, 0.0d0, 0.0d0 ]

      thermE0 = 123.456d0
      thermE  = thermE0

      call thermo_NHC_global(dt, M, ndim, natom, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      call check(error, maxval(abs(p - p0)) <= tol, "p changed even though mass=0 should bypass thermostat")
      call check(error, thermE == thermE0, "thermE changed even though mass=0 bypass should return early")

      deallocate(mNHC, rNHC, pNHC)
   end subroutine test_nhc_global_mass0_bypass



end module test_thermo
