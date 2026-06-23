!> Test suite for ThermoModule (Langevin and Nose–Hoover chain)
!!
!! it tests
!!   - LangevinThermo_O (momentum updates only, phase ignored)
!!   - thermo_NHC_local
!!   - thermo_NHC_global

module test_thermo
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt, BoltzK
   use RandomModule, only: initialize_fortran_prng
   use ParticleModule
   use TrajectoryModule
   use ThermoModule, only: thermo_NHC_local, thermo_NHC_global, LangevinThermo_O
   implicit none
   private

   public :: collect_thermo_suite

contains
   subroutine collect_thermo_suite(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)

      tests = [ &
         new_unittest("Langevin: 1D HO OBABO canonical phase-space distribution", test_langevin_momentum_dist), &
         new_unittest("NHC: 1D HO VRORV canonical phase-space distribution", test_nhc_momentum_dist) &
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