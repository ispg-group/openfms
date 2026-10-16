!> Test suite for ThermoModule (Langevin and Nose–Hoover chain)
!!
!! it tests the momentum and position distributions for a 1D harmonic oscillator at T=300 K, comparing to the expected canonical distribution.
!! It prints to the folder /test_thermo the sampled (p,x) values to files for visual inspection of the distributions.
!!   - Langevin thermostat with OBABO splitting (LangevinThermo_O)
!!   - NHC thermostat with VRORV splitting (thermo_NHC_local)

module test_thermo
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt, BoltzK
   use RandomModule, only: initialize_fortran_prng
   use ParticleModule
   use TrajectoryModule
   use ThermoModule, only: thermo_NHC_local, thermo_NHC_global, LangevinThermo_O, thermo_MBDist
   implicit none
   private

   public :: collect_thermo_suite

contains
   subroutine collect_thermo_suite(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)

      tests = [ &
         new_unittest("Langevin: 1D HO OBABO canonical phase-space distribution", test_langevin_momentum_dist), &
         new_unittest("NHC: 1D HO VRORV canonical phase-space distribution", test_nhc_momentum_dist), &
         new_unittest("NHC IC: chain masses, zero positions, MB momenta, thermE consistency", test_nhc_ic_init), &
         new_unittest("Langevin IC: MB momenta per atom at correct temperature", test_langevin_ic_momenta), &
         new_unittest("NHC NVT: quartic potential canonical distribution", test_nhc_quartic) &
         ]
   end subroutine collect_thermo_suite

   subroutine test_langevin_momentum_dist(error)
   !! OBABO Langevin integrator for a 1D harmonic oscillator: V = 0.5*mass*omega^2*x^2
   type(error_type), allocatable, intent(out) :: error
   type(T_Trajectory) :: traj
   real(kind=DefReal) :: gamma, beta, dt, g1_0, g2_0
   real(kind=DefReal) :: mass, omega2, x, p, force, T
   real(kind=DefReal) :: psum, p2sum, pmean, pvar, target_var
   real(kind=DefReal) :: xsum, x2sum, xmean, xvar, target_var_x
   real(kind=DefReal) :: Ekin_sum, Epot_sum, Ekin_mean, Epot_mean, target_Emean
   integer, parameter :: N_equil = 2000, N_sample = 50000
   integer :: k

   call traj%create(numparticles=1, numstates=1)
   traj%StateID = 1 ! refer to ground state
   traj%Particle(1)%Mass = 1.0d0
   traj%ElecStruc%ModPot = 0.0d0   ! not set by FMS_CreateElectronicStructure

   mass   = traj%Particle(1)%Mass
   omega2 = 1.0d0   ! omega^2; V = 0.5*mass*omega2*x^2,  F = -mass*omega2*x
   gamma  = 5.0d0
   T = 300.0d0
   beta = 1.0d0 / (BoltzK*T)
   dt     = 0.1d0
   target_var   = mass / beta                     ! <p²> = m/β
   target_var_x = 1.0d0 / (mass * omega2 * beta)  ! <x²> = 1/(mω²β)
   target_Emean = 0.5d0 / beta                     ! <Ekin> = <Epot> = ½kT (equipartition)

   ! Start at origin with zero momentum
   call traj%set_pos(1, 1, 0.0d0)
   call traj%set_mom(1, 1, 0.0d0)

   ! DerivMat stores the gradient dV/dx = mass*omega2*x; force = -DerivMat.
   traj%ElecStruc%DerivMat(1, 1, 1) = 0.0d0   ! dV/dx at x=0
   traj%ElecStruc%PotEn(1)          = 0.0d0   ! V(0) = 0
   traj%ESFlags%ZDerivCurrent       = .true.
   traj%ESFlags%ZPotEnCurrent       = .true.

   call initialize_fortran_prng(77777)

   psum  = 0.0d0
   p2sum = 0.0d0
   xsum  = 0.0d0
   x2sum = 0.0d0
   Ekin_sum = 0.0d0
   Epot_sum = 0.0d0

   open(unit=43, file="test_thermo/lang_momentum_T=300.dat", status="replace", action="write")

   do k = 1, N_equil + N_sample

      ! --- O (dt/2): first thermostat half-step
      call LangevinThermo_O(traj, gamma, beta, dt, g1_0, g2_0)

      ! --- B (dt/2): momentum half-kick from HO force
      x = traj%get_pos(1, 1)
      force = -mass * omega2 * x
      p = traj%get_mom(1, 1) + 0.5d0 * dt * force
      call traj%set_mom(1, 1, p)

      ! --- A (dt): full position step
      p = traj%get_mom(1, 1)
      x = traj%get_pos(1, 1) + dt * p / mass
      call traj%set_pos(1, 1, x)   ! triggers geom_changed() -> flags reset
      traj%ElecStruc%DerivMat(1, 1, 1) = mass * omega2 * x
      traj%ElecStruc%PotEn(1)          = 0.5d0 * mass * omega2 * x * x
      traj%ESFlags%ZDerivCurrent       = .true.
      traj%ESFlags%ZPotEnCurrent       = .true.

      ! --- B (dt/2): momentum half-kick from HO force
      force = -mass * omega2 * x
      p = traj%get_mom(1, 1) + 0.5d0 * dt * force
      call traj%set_mom(1, 1, p)

      ! --- O (dt/2): second thermostat half-step
      call LangevinThermo_O(traj, gamma, beta, dt, g1_0, g2_0)

      if (k > N_equil) then
         p = traj%get_mom(1, 1)
         write(43,'(G15.8,1X,G15.8)') p, x
         psum     = psum     + p
         p2sum    = p2sum    + p * p
         xsum     = xsum     + x
         x2sum    = x2sum    + x * x
         Ekin_sum = Ekin_sum + p * p / (2.0d0 * mass)
         Epot_sum = Epot_sum + 0.5d0 * mass * omega2 * x * x
      end if

   end do

   close(43)

   pmean = psum  / real(N_sample, kind=DefReal)
   pvar  = p2sum / real(N_sample, kind=DefReal) - pmean**2
   xmean = xsum  / real(N_sample, kind=DefReal)
   xvar  = x2sum / real(N_sample, kind=DefReal) - xmean**2

   call check(error, abs(pmean) < 0.05d0, &
              "Momentum mean not near zero: " // fmt_real(pmean))
   if (allocated(error)) return
   call check(error, abs(pvar - target_var) / target_var < 0.05d0, &
              "Momentum variance wrong: got " // fmt_real(pvar) // &
              " expected " // fmt_real(target_var))
   if (allocated(error)) return
   call check(error, abs(xmean) < 0.05d0, &
              "Position mean not near zero: " // fmt_real(xmean))
   if (allocated(error)) return
   call check(error, abs(xvar - target_var_x) / target_var_x < 0.05d0, &
              "Position variance wrong: got " // fmt_real(xvar) // &
              " expected " // fmt_real(target_var_x))
   if (allocated(error)) return
   Ekin_mean = Ekin_sum / real(N_sample, kind=DefReal)
   Epot_mean = Epot_sum / real(N_sample, kind=DefReal)
   call check(error, abs(Ekin_mean - target_Emean) / target_Emean < 0.05d0, &
              "Mean kinetic energy wrong: got " // fmt_real(Ekin_mean) // &
              " expected " // fmt_real(target_Emean))
   if (allocated(error)) return
   call check(error, abs(Epot_mean - target_Emean) / target_Emean < 0.05d0, &
              "Mean potential energy wrong: got " // fmt_real(Epot_mean) // &
              " expected " // fmt_real(target_Emean))

   call traj%destroy()
   end subroutine test_langevin_momentum_dist


   subroutine test_nhc_momentum_dist(error)
   !! VRORV (NHC half-step / Velocity Verlet) for a 1D harmonic oscillator:
   !!    V(x) = 0.5 * mass * omega^2 * x^2
   !! Parameters: mass=1, omega=1, T=300 K
   !!    target_var_p = m/β = mkT,  target_var_x = 1/(mω²β) = kT/(mω²)
   !! NHC chain length M=4, mNHC = 1/(β*ω²)
   type(error_type), allocatable, intent(out) :: error
   real(kind=DefReal) :: mass, omega2, beta, dt, T
   real(kind=DefReal) :: x, p, force, thermE
   real(kind=DefReal) :: psum, p2sum, pmean, pvar, target_var
   real(kind=DefReal) :: xsum, x2sum, xmean, xvar, target_var_x
   real(kind=DefReal) :: Ekin_sum, Epot_sum, Ekin_mean, Epot_mean, target_Emean
   integer, parameter :: M = 4
   integer, parameter :: N_equil = 10000, N_sample = 50000
   real(kind=DefReal), dimension(M) :: rNHC, pNHC, mNHC
   integer :: k

   mass   = 1.0d0
   omega2 = 1.0d0
   T = 300.0d0
   beta   = 1.0d0 / (BoltzK*T)
   dt     = 0.1d0
   target_var   = mass / beta                     ! <p²> = m/β
   target_var_x = 1.0d0 / (mass * omega2 * beta)  ! <x²> = 1/(mω²β)
   target_Emean = 0.5d0 / beta                     ! <Ekin> = <Epot> = ½kT (equipartition)

   ! Standard NHC masses: mNHC = NTotDim * kT / omega_ref^2, NTotDim=1
   mNHC(:) = 1.0d0 / (beta * omega2)   ! = 333.3
   rNHC(:) = 0.5d0
   pNHC(:) = 0.0d0

   ! Non-zero x so the first B-step gives a nonzero impulse (x=0,p=0 is degenerate)
   x = 0.0d0
   p = 0.5d0

   psum  = 0.0d0
   p2sum = 0.0d0
   xsum  = 0.0d0
   x2sum = 0.0d0
   Ekin_sum = 0.0d0
   Epot_sum = 0.0d0


   open(unit=42, file="test_thermo/nhc_momentum_T=300.dat", status="replace", action="write")

   do k = 1, N_equil + N_sample

      ! --- NHC (dt/2): thermostat half-step
      call thermo_NHC_local(0.5d0 * dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      ! --- B (dt/2): momentum half-kick from HO force
      force = -mass * omega2 * x
      p = p + 0.5d0 * dt * force

      ! --- A (dt): full position step
      x = x + dt * p / mass

      ! --- B (dt/2): momentum half-kick from HO force (at new x)
      force = -mass * omega2 * x
      p = p + 0.5d0 * dt * force

      ! --- NHC (dt/2): thermostat half-step
      call thermo_NHC_local(0.5d0 * dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      if (k > N_equil) then
         write(42,'(G15.8,1X,G15.8)') p, x
         psum     = psum     + p
         p2sum    = p2sum    + p * p
         xsum     = xsum     + x
         x2sum    = x2sum    + x * x
         Ekin_sum = Ekin_sum + p * p / (2.0d0 * mass)
         Epot_sum = Epot_sum + 0.5d0 * mass * omega2 * x * x
      end if

   end do

   close(42)

   pmean = psum  / real(N_sample, kind=DefReal)
   pvar  = p2sum / real(N_sample, kind=DefReal) - pmean**2
   xmean = xsum  / real(N_sample, kind=DefReal)
   xvar  = x2sum / real(N_sample, kind=DefReal) - xmean**2

   call check(error, abs(pmean) < 0.05d0, &
              "NHC momentum mean not near zero: " // fmt_real(pmean))
   if (allocated(error)) return
   call check(error, abs(pvar - target_var) / target_var < 0.05d0, &
              "NHC momentum variance wrong: got " // fmt_real(pvar) // &
              " expected " // fmt_real(target_var))
   if (allocated(error)) return
   call check(error, abs(xmean) < 0.05d0, &
              "NHC position mean not near zero: " // fmt_real(xmean))
   if (allocated(error)) return
   call check(error, abs(xvar - target_var_x) / target_var_x < 0.05d0, &
              "NHC position variance wrong: got " // fmt_real(xvar) // &
              " expected " // fmt_real(target_var_x))
   if (allocated(error)) return
   Ekin_mean = Ekin_sum / real(N_sample, kind=DefReal)
   Epot_mean = Epot_sum / real(N_sample, kind=DefReal)
   call check(error, abs(Ekin_mean - target_Emean) / target_Emean < 0.05d0, &
              "NHC mean kinetic energy wrong: got " // fmt_real(Ekin_mean) // &
              " expected " // fmt_real(target_Emean))
   if (allocated(error)) return
   call check(error, abs(Epot_mean - target_Emean) / target_Emean < 0.05d0, &
              "NHC mean potential energy wrong: got " // fmt_real(Epot_mean) // &
              " expected " // fmt_real(target_Emean))

   end subroutine test_nhc_momentum_dist


   subroutine test_nhc_ic_init(error)
   !! Tests FMS_InitializeNHC for correct NHC chain initial conditions:
   !!   1. Thermostat flags: zActive=.true., type='NHC', NumNHCChain=M
   !!   2. Chain masses: mNHC(1) = Ndim*kT*tau^2,  mNHC(i>1) = kT*tau^2
   !!   3. Chain positions initialized to zero
   !!   4. Chain momenta from Maxwell-Boltzmann: <pNHC(i)^2> = mNHC(i)/beta
   !!   5. thermE consistent with sampled pNHC and rNHC=0
   type(error_type), allocatable, intent(out) :: error
   type(T_Trajectory) :: traj
   real(kind=DefReal) :: T, beta, tau, tol
   integer, parameter :: M = 4
   integer, parameter :: N_sample = 20000
   real(kind=DefReal) :: p2sum(M), p2mean(M), target_p2(M), thermE_expected
   integer :: k, i

   T    = 300.0d0
   tau  = 1000.0d0       ! NHC relaxation time in au (~24 fs)
   beta = 1.0d0 / (BoltzK * T)
   tol  = 1.0d-10

   ! 3-atom trajectory => NumDimensions = 9
   call traj%create(numparticles=3, numstates=1)
   traj%StateID = 1

   call initialize_fortran_prng(54321)

   p2sum(:) = 0.0d0
   do k = 1, N_sample
      call FMS_InitializeNHC(traj, M, T, tau)
      do i = 1, M
         p2sum(i) = p2sum(i) + traj%ThermoInfo%GlobalNHC%pNHC(i)**2
      end do
   end do

   ! --- 1. Flags ---
   call check(error, traj%ThermoInfo%zActive, &
              "NHC IC: ThermoInfo%zActive should be .true.")
   if (allocated(error)) return
   call check(error, trim(traj%ThermoInfo%thermostat_type) == 'NHC', &
              "NHC IC: thermostat_type wrong: " // trim(traj%ThermoInfo%thermostat_type))
   if (allocated(error)) return
   call check(error, traj%ThermoInfo%NumNHCChain == M, &
              "NHC IC: NumNHCChain wrong: " // fmt_int(traj%ThermoInfo%NumNHCChain))
   if (allocated(error)) return

   ! --- 2. Chain masses ---
   call check(error, &
              abs(traj%ThermoInfo%GlobalNHC%mNHC(1) - real(traj%NumDimensions, DefReal) / beta * tau**2) < tol, &
              "NHC IC: mNHC(1) wrong: " // fmt_real(traj%ThermoInfo%GlobalNHC%mNHC(1)))
   if (allocated(error)) return
   do i = 2, M
      call check(error, abs(traj%ThermoInfo%GlobalNHC%mNHC(i) - tau**2 / beta) < tol, &
                 "NHC IC: mNHC(" // fmt_int(i) // ") wrong: " // &
                 fmt_real(traj%ThermoInfo%GlobalNHC%mNHC(i)))
      if (allocated(error)) return
   end do

   ! --- 3. Positions zero ---
   do i = 1, M
      call check(error, traj%ThermoInfo%GlobalNHC%rNHC(i) == 0.0d0, &
                 "NHC IC: rNHC(" // fmt_int(i) // ") should be zero")
      if (allocated(error)) return
   end do

   ! --- 4. pNHC Maxwell-Boltzmann: <pNHC(i)^2> = mNHC(i)/beta ---
   do i = 1, M
      target_p2(i) = traj%ThermoInfo%GlobalNHC%mNHC(i) / beta
      p2mean(i)    = p2sum(i) / real(N_sample, DefReal)
      call check(error, abs(p2mean(i) - target_p2(i)) / target_p2(i) < 0.05d0, &
                 "NHC IC: <pNHC(" // fmt_int(i) // ")^2> wrong: got " // &
                 fmt_real(p2mean(i)) // " expected " // fmt_real(target_p2(i)))
      if (allocated(error)) return
   end do

   ! --- 5. thermE = sum(0.5*p^2/m + r/beta); with r=0 => purely kinetic ---
   thermE_expected = 0.0d0
   do i = 1, M
      thermE_expected = thermE_expected + &
         0.5d0 * traj%ThermoInfo%GlobalNHC%pNHC(i)**2 / traj%ThermoInfo%GlobalNHC%mNHC(i)
   end do
   call check(error, abs(traj%ThermoInfo%thermE - thermE_expected) < tol, &
              "NHC IC: thermE wrong: got " // fmt_real(traj%ThermoInfo%thermE) // &
              " expected " // fmt_real(thermE_expected))
   if (allocated(error)) return

   call traj%destroy()
   end subroutine test_nhc_ic_init


   subroutine test_langevin_ic_momenta(error)
   !! Tests that thermo_MBDist correctly samples nuclear momenta for Langevin IC.
   !! Uses a 2-atom heterogeneous system (mass 1 and 12 au) to verify per-atom
   !! Maxwell-Boltzmann distribution: <p_atom^2> = ndim * m_atom / beta.
   !! This is the momentum-sampling step that precedes Langevin dynamics.
   type(error_type), allocatable, intent(out) :: error
   real(kind=DefReal) :: T, beta
   integer, parameter :: ndim = 3, natom = 2
   integer, parameter :: N_sample = 30000
   real(kind=DefReal) :: p(ndim, natom), mass(natom)
   real(kind=DefReal) :: p2sum(natom), p2mean(natom), target_p2(natom)
   integer :: k, iatom

   T    = 300.0d0
   beta = 1.0d0 / (BoltzK * T)
   mass(1) = 1.0d0    ! light atom
   mass(2) = 12.0d0   ! heavy atom (like carbon)

   call initialize_fortran_prng(99999)

   p2sum(:) = 0.0d0
   do k = 1, N_sample
      call thermo_MBDist(ndim, natom, p, mass, beta, .true.)
      do iatom = 1, natom
         p2sum(iatom) = p2sum(iatom) + sum(p(:, iatom)**2)
      end do
   end do

   ! <p_atom^2> = ndim * m_atom / beta  (equipartition: ndim/2 * kT per atom)
   do iatom = 1, natom
      target_p2(iatom) = real(ndim, DefReal) * mass(iatom) / beta
      p2mean(iatom)    = p2sum(iatom) / real(N_sample, DefReal)
      call check(error, abs(p2mean(iatom) - target_p2(iatom)) / target_p2(iatom) < 0.05d0, &
                 "Langevin IC: atom " // fmt_int(iatom) // &
                 " <p^2> wrong: got " // fmt_real(p2mean(iatom)) // &
                 " expected " // fmt_real(target_p2(iatom)))
      if (allocated(error)) return
   end do

   end subroutine test_langevin_ic_momenta

   subroutine test_nhc_quartic(error)
   !! VRORV NHC integrator for a 1D quartic oscillator: V(x) = lambda * x^4
   !! At equilibrium the following exact results hold:
   !!   <p>    = 0             (symmetry)
   !!   <x>    = 0             (symmetry)
   !!   <p^2>  = m/beta        (kinetic temperature)
   !!   <Ekin> = 1/(2*beta)    (equipartition)
   !!   <V>    = 1/(4*beta)    (virial: x*dV/dx = 4*V = kT)
   type(error_type), allocatable, intent(out) :: error
   real(kind=DefReal) :: mass, lambda, beta, dt, T
   real(kind=DefReal) :: x, p, force, thermE
   real(kind=DefReal) :: psum, p2sum, pmean, pvar, target_var_p
   real(kind=DefReal) :: xsum, xmean
   real(kind=DefReal) :: Ekin_sum, Epot_sum, Ekin_mean, Epot_mean
   real(kind=DefReal) :: target_Ekin, target_Epot
   real(kind=DefReal) :: sigma_p, sigma_x, u1, u2
   integer, parameter :: M = 4
   integer, parameter :: N_equil = 20000, N_sample = 80000
   real(kind=DefReal), dimension(M) :: rNHC, pNHC, mNHC
   integer :: k

   mass   = 1.0d0
   lambda = 0.25d0        ! V(x) = 0.25*x^4,  F(x) = -x^3
   T      = 300.0d0
   beta   = 1.0d0 / (BoltzK * T)
   dt     = 0.1d0

   ! NHC masses: mNHC = kT * tau_ref^2 with tau_ref = 1 au
   mNHC(:) = 1.0d0 / beta
   rNHC(:) = 0.0d0
   pNHC(:) = 0.0d0

   ! Sample initial (x, p) from approximate canonical distribution using Box-Muller.
   ! p ~ N(0, sqrt(m/beta));  x ~ N(0, sigma_x) with sigma_x from the quartic virial:
   !   <x^2>_quartic = Gamma(3/4)/Gamma(1/4) / sqrt(beta*lambda) ≈ 0.33813/sqrt(beta*lambda)
   call initialize_fortran_prng(42)
   sigma_p = sqrt(mass / beta)
   sigma_x = sqrt(0.33813d0 / sqrt(beta * lambda))
   call random_number(u1); call random_number(u2)
   p = sigma_p * sqrt(-2.0d0 * log(u1)) * cos(8.0d0 * atan(1.0d0) * u2)
   call random_number(u1); call random_number(u2)
   x = sigma_x * sqrt(-2.0d0 * log(u1)) * cos(8.0d0 * atan(1.0d0) * u2)

   ! Exact canonical-ensemble targets
   target_var_p = mass / beta    ! <p^2>  = m/beta
   target_Ekin  = 0.5d0 / beta   ! <Ekin> = kT/2
   target_Epot  = 0.25d0 / beta  ! <V>    = kT/4  (virial theorem)

   psum  = 0.0d0;  p2sum    = 0.0d0
   xsum  = 0.0d0
   Ekin_sum = 0.0d0;  Epot_sum = 0.0d0

   open(unit=45, file="test_thermo/nhc_quartic_T=300.dat", status="replace", action="write")
   write(45,'(A)') "# p                x              Ekin             V=lambda*x^4"

   do k = 1, N_equil + N_sample

      ! --- NHC (dt/2)
      call thermo_NHC_local(0.5d0 * dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      ! --- B (dt/2): impulse from quartic force  F = -dV/dx = -4*lambda*x^3
      force = -4.0d0 * lambda * x**3
      p = p + 0.5d0 * dt * force

      ! --- A (dt): free drift
      x = x + dt * p / mass

      ! --- B (dt/2): impulse at new x
      force = -4.0d0 * lambda * x**3
      p = p + 0.5d0 * dt * force

      ! --- NHC (dt/2)
      call thermo_NHC_local(0.5d0 * dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)

      if (k > N_equil) then
         write(45,'(4(G15.8,1X))') p, x, p*p/(2.0d0*mass), lambda*x**4
         psum     = psum     + p
         p2sum    = p2sum    + p * p
         xsum     = xsum     + x
         Ekin_sum = Ekin_sum + p * p / (2.0d0 * mass)
         Epot_sum = Epot_sum + lambda * x**4
      end if

   end do

   close(45)

   pmean = psum / real(N_sample, kind=DefReal)
   pvar  = p2sum / real(N_sample, kind=DefReal) - pmean**2
   xmean = xsum / real(N_sample, kind=DefReal)
   Ekin_mean = Ekin_sum / real(N_sample, kind=DefReal)
   Epot_mean = Epot_sum / real(N_sample, kind=DefReal)

   call check(error, abs(pmean) < 0.05d0, &
              "Quartic NHC: momentum mean not near zero: " // fmt_real(pmean))
   if (allocated(error)) return
   call check(error, abs(pvar - target_var_p) / target_var_p < 0.05d0, &
              "Quartic NHC: <p^2> wrong: got " // fmt_real(pvar) // &
              " expected " // fmt_real(target_var_p))
   if (allocated(error)) return
   call check(error, abs(xmean) < 0.05d0, &
              "Quartic NHC: position mean not near zero: " // fmt_real(xmean))
   if (allocated(error)) return
   call check(error, abs(Ekin_mean - target_Ekin) / target_Ekin < 0.05d0, &
              "Quartic NHC: <Ekin> wrong: got " // fmt_real(Ekin_mean) // &
              " expected " // fmt_real(target_Ekin))
   if (allocated(error)) return
   call check(error, abs(Epot_mean - target_Epot) / target_Epot < 0.05d0, &
              "Quartic NHC: <V> wrong: got " // fmt_real(Epot_mean) // &
              " expected " // fmt_real(target_Epot))

   end subroutine test_nhc_quartic


   function fmt_int(i) result(s)
      integer, intent(in) :: i
      character(len=20) :: s
      write(s, '(I0)') i
   end function fmt_int

   function fmt_real(x) result(s)
      real(kind=DefReal), intent(in) :: x
      character(len=32) :: s
      write(s, '(G15.6)') x
   end function fmt_real

end module test_thermo