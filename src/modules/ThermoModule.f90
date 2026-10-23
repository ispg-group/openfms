! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module ThermoModule
   use GlobalModule
   implicit none
   private
   public :: thermo_init, thermo, thermo_bussi_global, thermo_bussi_local, thermo_NormDist, thermo_MBDist, &
             thermo_NHC_MBDist, thermo_NHC_local, thermo_NHC_global
   character(len=8) :: therm = '' !which thermostat to use in simple interface
   double precision :: beta_s = 0.0 !inverse temperature in simple interface
   double precision :: thermtime = 0.0 !Thermostat relaxation time (in au) for simple interface
   logical :: zcom_s = .true. !Do center of mass removal in simple interface?
   integer(4), parameter :: i4one = 1 !Ensure type is consistent with uses of integer(4)

contains

   subroutine thermo_init(iseed, therm_in, beta_in, thermtime_in, zcom_in)
      integer(4), intent(in) :: iseed
      character(len=8), intent(in), optional :: therm_in
      double precision, intent(in), optional :: beta_in, thermtime_in
      logical, intent(in), optional :: zcom_in

      double precision :: foo
      integer(4) :: idum

      idum = sign(iseed, -i4one) !ensure seed is negative
      call thermo_ran(foo, idum) !Initialize random number generator

      if (present(therm_in)) therm = therm_in !Set the thermostat for the simple interface
      if (present(beta_in)) beta_s = beta_in !Set the temperature for the simple interface
      if (present(zcom_in)) zcom_s = zcom_in !Set the center of mass removal flag for the simple interface
      if (present(thermtime_in)) thermtime = thermtime_in !Set the thermostat relaxation time for the simple interface

   end subroutine thermo_init

   subroutine thermo(ndim, natom, p, mass, dt, thermE)
      integer, intent(in) :: ndim, natom
      double precision, intent(inout) :: p(ndim, natom), thermE
      double precision, intent(in) :: mass(natom), dt
      double precision :: tau

      tau = thermtime / dt

      select case (therm)
      case ('bussi_g')
         call thermo_bussi_global(ndim, natom, p, mass, beta_s, tau, thermE, zcom_s)
      case ('bussi_l')
         call thermo_bussi_local(ndim, natom, p, mass, beta_s, tau, thermE, zcom_s)
      case ('mbdist')
         call thermo_MBDist(ndim, natom, p, mass, beta_s, zcom_s)
      case default !do nothing, and return

      end select
   end subroutine thermo

   function thermo_NormDist(sigma) result(x)
!This subroutine generates a normal distribution exp(-x^2/2sigma^2)
!This function is currently not working because gasdev() is yet to be implemented.
      double precision :: x
      double precision, intent(in) :: sigma

      x = gasdev() * sigma

   end function thermo_NormDist

   subroutine thermo_MBDist(ndim, natom, p, mass, beta, zcom)
      integer, intent(in) :: ndim, natom
      double precision, intent(inout) :: p(ndim, natom)
      double precision, intent(in) :: mass(natom)
      double precision, intent(in) :: beta
      logical, intent(in) :: zcom

      double precision :: mxw_unitf, sigma
      integer :: iatom, idm

!Sample momenta from Maxwell-Boltzmann distribution
      mxw_unitf = 1.0 / sqrt(beta)
      do iatom = 1, natom
         sigma = mxw_unitf * sqrt(mass(iatom))
         do idm = 1, ndim
            p(idm, iatom) = thermo_NormDist(sigma)
         end do
      end do

!Project out center of mass velocity if not a degree of freedom
      if (.not. zcom) then
         call vcom_project(ndim, natom, p, mass)
      end if
   end subroutine thermo_MBDist

   subroutine thermo_bussi_global(ndim, natom, p, mass, beta, tau, thermE, zcom)
      integer, intent(in) :: ndim, natom
      double precision, intent(inout) :: p(ndim, natom)
      double precision, intent(in) :: mass(natom)
      double precision, intent(in) :: beta, tau
      double precision, intent(inout) :: thermE
      logical, intent(in) :: zcom

      double precision :: ekin_old, ekin_new, ekinfac, sigma, signfac, vscale
      integer :: ndeg, iatom, idm

      if (.not. zcom) then
!Do not include center of mass as degree of freedom
         call vcom_project(ndim, natom, p, mass)
         ndeg = ndim * natom - ndim !Rotations are not frozen out, so keep them as degrees of freedom
      else
         ndeg = ndim * natom
      end if

      ekin_old = 0.0
      do iatom = 1, natom
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_old = ekin_old + ekinfac * p(idm, iatom) * p(idm, iatom)
         end do
      end do

      sigma = 0.5 * real(ndeg) / beta

      ekin_new = resamplekin(ekin_old, sigma, ndeg, tau)

      signfac = sign(1.0d0, ekin_new)
      ekin_new = abs(ekin_new)

      thermE = thermE + ekin_old - ekin_new

      vscale = signfac * sqrt(ekin_new / ekin_old)

      do iatom = 1, natom
         do idm = 1, ndim
            p(idm, iatom) = p(idm, iatom) * vscale
         end do
      end do

   end subroutine thermo_bussi_global

   subroutine thermo_bussi_local(ndim, natom, p, mass, beta, tau, thermE, zcom)
      integer, intent(in) :: ndim, natom
      double precision, intent(inout) :: p(ndim, natom)
      double precision, intent(in) :: mass(natom)
      double precision, intent(in) :: beta, tau
      double precision, intent(inout) :: thermE
      logical, intent(in) :: zcom

      double precision :: vtot(ndim)
      double precision :: c1, c2, c2fac, ekin_old, ekin_new, ekinfac
      integer :: iatom, idm

      c1 = exp(-1.0 / tau)
      c2 = sqrt(1.0 - c1 * c1)

      ekin_old = 0.0
      do iatom = 1, natom
         c2fac = c2 * sqrt(mass(iatom) / beta)
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_old = ekin_old + ekinfac * p(idm, iatom) * p(idm, iatom)
            p(idm, iatom) = c1 * p(idm, iatom) + c2fac * gasdev()
         end do
      end do

      if (.not. zcom) then
!Do not include center of mass as degree of freedom
         vtot(1:ndim) = 0.0
         do iatom = 1, natom
            do idm = 1, ndim
               vtot(idm) = vtot(idm) + p(idm, iatom) / mass(iatom)
            end do
         end do
         vtot(1:ndim) = vtot(1:ndim) / real(natom)
         do iatom = 1, natom
            do idm = 1, ndim
               p(idm, iatom) = p(idm, iatom) - mass(iatom) * vtot(idm)
            end do
         end do
      end if

      ekin_new = 0.0
      do iatom = 1, natom
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_new = ekin_new + ekinfac * p(idm, iatom) * p(idm, iatom)
         end do
      end do
!Calculate new thermostat energy
      thermE = thermE + ekin_old - ekin_new

   end subroutine thermo_bussi_local

   function resamplekin(kk, sigma, ndeg, taut)
      double precision :: resamplekin
      double precision, intent(in) :: kk ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
      double precision, intent(in) :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
      integer, intent(in) :: ndeg ! number of degrees of freedom of the atoms to be thermalized
      double precision, intent(in) :: taut ! relaxation time of the thermostat, in units of 'how often this routine is called'
      double precision :: factor, rr
      double precision :: dndeg
      dndeg = real(ndeg)
      if (taut > 0.1) then
         factor = exp(-1.0 / taut)
      else
         factor = 0.0
      end if
      rr = gasdev()
      resamplekin = kk + (1.0 - factor) * (sigma * (sumnoises(ndeg - 1) + rr**2) / dndeg - kk) &
                    + 2.0 * rr * sqrt(kk * sigma / dndeg * (1.0 - factor) * factor)

! Transfer sign to resamplekin
      resamplekin = sign(resamplekin, rr + sqrt(dndeg * factor / (sigma * (1.0 - factor))))

   end function resamplekin

   double precision function sumnoises(nn)
      integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
      if (nn == 0) then
         sumnoises = 0.0
      else if (nn == 1) then
         sumnoises = gasdev()**2
      else if (modulo(nn, 2) == 0) then
         sumnoises = 2.0 * gamdev(nn / 2)
      else
         sumnoises = 2.0 * gamdev((nn - 1) / 2) + gasdev()**2
      end if
   end function sumnoises

! gamma-distributed random number
! TODO: Reimplement this
   function gamdev(ia)
      double precision :: gamdev
      integer, intent(in) :: ia

      gamdev = 0.0d0
      call FMS_DieError('ERROR: gamdev not implemented')
   end function gamdev

! TODO: Reimplement this
   function gasdev()
      double precision :: gasdev

      gasdev = 0.0d0
      call FMS_DieError('ERROR: gasdev not implemented')
   end function gasdev

   subroutine thermo_ran(rnd, iseed)
! interface to random number generators
      integer(4), intent(in), optional :: iseed
      double precision, intent(in) :: rnd

      call FMS_DieError('ERROR: thermo_ran not implemented')
      if (present(iseed)) then
         ! TODO: Use FMS_ranb
         !rnd=ran1(iseed)
      end if
!rnd=ran1()
   end subroutine thermo_ran

   subroutine vcom_project(ndim, natom, p, mass)
      integer, intent(in) :: ndim, natom
      double precision, intent(inout) :: p(ndim, natom)
      double precision, intent(in) :: mass(natom)

      double precision :: vtot(ndim)
      integer :: iatom, idm

      vtot(1:ndim) = 0.0
      do iatom = 1, natom
         do idm = 1, ndim
            vtot(idm) = vtot(idm) + p(idm, iatom) / mass(iatom)
         end do
      end do
      vtot(1:ndim) = vtot(1:ndim) / real(natom)
      do iatom = 1, natom
         do idm = 1, ndim
            p(idm, iatom) = p(idm, iatom) - mass(iatom) * vtot(idm)
         end do
      end do

   end subroutine vcom_project

   !> Draw independent Maxwell-Boltzmann momenta for the Nose-Hoover chain.
   !> Each chain momentum is sampled with variance \(m_i/\beta\), and
   !> thermE is initialized to the kinetic energy of the sampled chain.
   !> This is the canonical momentum distribution used with the NHC
   !> extended-system dynamics of Martyna et al. (1996).
   !>
   !> Reference: G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and
   !> M. L. Klein, "Explicit reversible integrators for extended systems
   !> dynamics," Molecular Physics, 87(5), 1117-1157 (1996).
   !> https://doi.org/10.1080/00268979600100761
   subroutine thermo_NHC_MBDist(beta, M, pNHC, mNHC, thermE)
      !> Number of thermostats in the Nose-Hoover chain.
      integer, intent(in) :: M
      !> Inverse temperature, \(1/(k_B T)\).
      real(8), intent(in) :: beta
      !> Chain kinetic energy, initialized and accumulated by this routine.
      real(8), intent(out) :: thermE
      !> Thermostat masses, one for each chain degree of freedom.
      real(8), intent(in), dimension(M) :: mNHC
      !> Sampled momenta of the Nose-Hoover chain.
      real(8), intent(out), dimension(M) :: pNHC
      !> Standard deviation of the Maxwell-Boltzmann distribution for the
      !> current chain momentum.
      real(8) :: sigma
      !> Index of the chain thermostat being initialized.
      integer :: i

      thermE = 0.0d0
      do i = 1, M
         sigma = sqrt(mNHC(i) / beta) !Maxwell-Boltzmann distribution of Nose momenta
         pNHC(i) = thermo_NormDist(sigma)
         thermE = thermE + 0.5d0 * pNHC(i) * pNHC(i) / mNHC(i)
      end do
   end subroutine thermo_NHC_MBDist

   !> Propagate a single scalar momentum coupled to a Nose-Hoover chain.
   !> The scalar degree of freedom is treated as a local thermostat group,
   !> so the first chain thermostat is driven by \(p^2/m-1/\beta\)
   !> (twice the particle kinetic energy minus the target value).
   !> The chain and particle momenta are advanced with
   !> thermo_NHC_prop; afterward thermE is recomputed as the chain kinetic
   !> plus chain potential energy. The reversible NHC propagator follows
   !> Eq. (35) of Martyna et al. (1996); see thermo_NHC_prop for the
   !> update equations.

   !> Reference: G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and
   !> M. L. Klein, "Explicit reversible integrators for extended systems
   !> dynamics," Molecular Physics, 87(5), 1117-1157 (1996).
   !> https://doi.org/10.1080/00268979600100761
   subroutine thermo_NHC_local(dt, M, p, mass, beta, rNHC, pNHC, mNHC, thermE)
      implicit none
      !> Number of thermostats in the chain.
      integer, intent(in) :: M
      !> timestep for the thermostat propagation
      real(8), intent(in) :: dt
      !> Particle momentum; updated in place by the chain scaling factor.
      real(8), intent(inout) :: p
      !> Mass associated with the particle momentum.
      real(8), intent(in) :: mass
      !> Inverse temperature, \(1/(k_B T)\).
      real(8), intent(in) :: beta
      !> Extended thermostat energy, overwritten with the current chain energy.
      real(8), intent(out) :: thermE
      !> Masses of the Nose-Hoover chain thermostats.
      real(8), intent(in), dimension(M) :: mNHC
      !> Positions of the Nose-Hoover chain thermostats.
      real(8), intent(inout), dimension(M) :: rNHC
      !> Momenta of the Nose-Hoover chain thermostats.
      real(8), intent(inout), dimension(M) :: pNHC
      !> Twice the particle kinetic energy, \(p^2/m\), coupling to chain 1.
      real(8) :: akin
      !> Accumulated multiplicative scaling applied to the particle momentum.
      real(8) :: scale
      !> Index of the chain thermostat used to accumulate its energy.
      integer :: l

      akin = p * p / mass
      call thermo_NHC_prop(dt, akin, beta, 1, M, rNHC, pNHC, mNHC, scale)

      !Apply the chain's accumulated scaling to the local momentum.
      p = p * scale

      !Include each chain momentum and chain-position contribution to the
      !extended thermostat energy.
      thermE = 0.0d0
      do l = 1, M
         thermE = thermE + 0.5d0 * pNHC(l) * pNHC(l) / mNHC(l) + rNHC(l) / beta
      end do

   end subroutine thermo_NHC_local

   !> Propagate all particle momenta as one global Nose-Hoover thermostat
   !> group. The chain propagator and
   !> common momentum scaling follow Eq. (35) of Martyna et al. (1996);
   !> see thermo_NHC_prop for the update equations and full citation.
   !> If any particle mass is effectively zero, this routine returns
   !> without changing the inputs.

   !> Reference: G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and
   !> M. L. Klein, "Explicit reversible integrators for extended systems
   !> dynamics," Molecular Physics, 87(5), 1117-1157 (1996).
   !> https://doi.org/10.1080/00268979600100761
   subroutine thermo_NHC_global(dt, M, ndim, natom, p, mass, beta, rNHC, pNHC, mNHC, thermE)
      implicit none
      !> Number of thermostats in the chain.
      integer, intent(in) :: M
      !> Number of Cartesian components per particle.
      integer, intent(in) :: ndim
      !> Number of particles in the global thermostat group.
      integer, intent(in) :: natom
      !> Particle momenta, updated in place by the common chain scaling.
      real(8), intent(inout), dimension(ndim, natom) :: p
      !> Mass of each particle.
      real(8), intent(in), dimension(natom) :: mass
      !> timestep
      real(8), intent(in) :: dt
      !> Inverse temperature, \(1/(k_B T)\).
      real(8), intent(in) :: beta
      !> Extended thermostat energy, recomputed when particle masses are valid;
      !> preserved when this routine returns early for a zero mass.
      real(8), intent(inout) :: thermE
      !> Masses of the Nose-Hoover chain thermostats.
      real(8), intent(in), dimension(M) :: mNHC
      !> Positions of the Nose-Hoover chain thermostats.
      real(8), intent(inout), dimension(M) :: rNHC
      !> Momenta of the Nose-Hoover chain thermostats.
      real(8), intent(inout), dimension(M) :: pNHC
      !> Twice the total particle kinetic energy, \(\sum_{a,\mu}p_{a\mu}^2/m_a\).
      real(8) :: akin
      !> Accumulated multiplicative scaling applied to all particle momenta.
      real(8) :: scale
      !> Index of the particle in the momentum and kinetic-energy sums.
      integer :: iatom
      !> Index of a Cartesian momentum component.
      integer :: idm
      !> Index of a chain thermostat used to accumulate chain energy.
      integer :: l
      !> Number of particle momentum components thermostatted globally.
      integer :: NTotDim

      if (any(mass < 1.0d-15)) return !Do not thermostat if mass is zero
      akin = 0.0d0
      do iatom = 1, natom
         do idm = 1, ndim
            akin = akin + p(idm, iatom) * p(idm, iatom) / mass(iatom)
         end do
      end do
      NTotDim = ndim * natom
      call thermo_NHC_prop(dt, akin, beta, NTotDim, M, rNHC, pNHC, mNHC, scale)
      !Apply the same accumulated scaling to every particle momentum.
      p = p * scale
      !Include the N-fold first-chain potential and the remaining chain
      !potential and kinetic contributions in the extended thermostat energy.
      thermE = 0.0d0
      thermE = thermE + 0.5d0 * pNHC(1) * pNHC(1) / mNHC(1) + dble(NTotDim) * rNHC(1) / beta
      do l = 2, M
         thermE = thermE + 0.5d0 * pNHC(l) * pNHC(l) / mNHC(l) + rNHC(l) / beta
      end do
   end subroutine thermo_NHC_global

   !> Advance the Nose-Hoover chain over dt and return the accumulated
   !> factor used to scale the coupled particle momenta according to Eq (35)
   !> in the paper.
   !> Uses a reversible chain propagator with nc subcycles and a
   !> five-stage Yoshida-Suzuki factorization.
   !>
   !> Reference: G. J. Martyna, M. E. Tuckerman, D. J. Tobias, and
   !> M. L. Klein, "Explicit reversible integrators for extended systems
   !> dynamics," Molecular Physics, 87(5), 1117-1157 (1996).
   !> https://doi.org/10.1080/00268979600100761
   subroutine thermo_NHC_prop(dt, akin, beta, N, M, rNHC, pNHC, mNHC, scale)
      implicit none
      !> Number of degrees of freedom coupled to the first chain thermostat.
      integer, intent(in) :: N
      !> Number of thermostats in the chain.
      integer, intent(in) :: M
      !> Duration of this thermostat propagation (usually half a system step).
      real(8), intent(in) :: dt
      !> Twice the kinetic energy of the degrees of freedom coupled to chain 1.
      real(8), intent(in) :: akin
      !> Inverse temperature, \(1/(k_B T)\).
      real(8), intent(in) :: beta
      !> Accumulated multiplicative scaling returned for system momenta.
      real(8), intent(out) :: scale
      !> Masses of the Nose-Hoover chain thermostats.
      real(8), intent(in), dimension(M) :: mNHC
      !> Positions of the Nose-Hoover chain thermostats, advanced in place.
      real(8), intent(inout), dimension(M) :: rNHC
      !> Momenta of the Nose-Hoover chain thermostats, advanced in place.
      real(8), intent(inout), dimension(M) :: pNHC

      !> Number of equal subcycles used to split the requested propagation.
      integer, parameter :: nc = 10
      !> Maximum number of stored Yoshida-Suzuki weights.
      integer, parameter :: nysmax = 7
      !> Number of Yoshida-Suzuki stages; five stages give a fourth-order
      !> symmetric composition.
      integer, parameter :: nys = 5

      !> Instantaneous generalized forces on the chain momenta.
      real(8), dimension(M) :: G
      !> Reciprocal masses of the chain thermostats.
      real(8), dimension(M) :: mNHC_inv
      !> Yoshida-Suzuki composition weights.
      real(8), dimension(nysmax) :: w
      !> Temperature in energy units, \(1/\beta\).
      real(8) :: beta_inv
      !> Duration of one of the nc equal integration subcycles.
      real(8) :: dt_nc
      !> Duration of the current weighted Yoshida-Suzuki stage.
      real(8) :: dts
      !> Exponential factor for a symmetric chain-momentum update or scaling.
      real(8) :: AA
      !> Target force contribution \(N/\beta\) for the first chain site.
      real(8) :: dNbeta_inv
      !> Index of the Yoshida-Suzuki stage.
      integer :: j
      !> Index of an equal propagation subcycle.
      integer :: k
      !> Index of a chain site in the update sweeps.
      integer :: l

      scale = 1.0d0
      !If any NHC mass is zero, we disable NHC propagation
      if (any(mNHC < 1.0d-15)) then
         write (6, *) 'mNHC is zero'
         return
      end if

      !Set up the selected Yoshida-Suzuki composition weights.
      select case (nys)
      case (1)
         w(1) = 1.0d0
      case (3)
         w(1) = 1.0d0 / (2.0d0 - 2.0d0**(1.0d0 / 3.0d0))
         w(2) = 1.0d0 - 2.0d0 * w(1)
         w(3) = w(1)
      case (5)
         w(1) = 1.0d0 / (4.0d0 - 4.0d0**(1.0d0 / 3.0d0))
         w(2) = w(1)
         w(3) = 1.0d0 - 4.0d0 * w(1)
         w(4) = w(1)
         w(5) = w(1)
      case (7)
         w(1) = -0.117767998417887d1
         w(2) = 0.235573213359357d0
         w(3) = 0.784513610477560d0
         w(4) = 1.0d0 - 2.0d0 * (w(1) + w(2) + w(3))
         w(5) = w(3)
         w(6) = w(2)
         w(7) = w(1)
      case default
         write (*, *) 'Error, nys=', nys, ' not coded'
      end select

      beta_inv = 1.0d0 / beta
      dt_nc = dt / dble(nc)
      mNHC_inv(1:M) = 1.0d0 / mNHC(1:M)
      dNbeta_inv = beta_inv * dble(N)

      !Initialize the generalized forces from the NHC equations of motion.
      G(1) = akin - dNbeta_inv
      G(2:M) = mNHC_inv(1:M - 1) * pNHC(1:M - 1) * pNHC(1:M - 1) - beta_inv

      do k = 1, nc
         do j = 1, nys
            dts = w(j) * dt_nc
            !Update NHC momenta
            pNHC(M) = pNHC(M) + 0.5d0 * dts * G(M)
            do l = M - 1, 1, -1
               AA = exp(-0.25d0 * dts * pNHC(l + 1) * mNHC_inv(l + 1))
               pNHC(l) = pNHC(l) * AA * AA + 0.5d0 * dts * G(l) * AA
            end do
            !Update particle momenta
            AA = exp(-dts * pNHC(1) * mNHC_inv(1))
            scale = scale * AA
            !Update NHC forces
            G(1) = scale * scale * akin - dNbeta_inv
            !Update NHC positions
            do l = 1, M
               rNHC(l) = rNHC(l) + pNHC(l) * dts * mNHC_inv(l)
            end do
            !Update NHC momenta
            do l = 1, M - 1
               AA = exp(-0.25d0 * dts * pNHC(l + 1) * mNHC_inv(l + 1))
               pNHC(l) = pNHC(l) * AA * AA + 0.5d0 * dts * G(l) * AA
               G(l + 1) = mNHC_inv(l) * pNHC(l) * pNHC(l) - beta_inv
            end do
            pNHC(M) = pNHC(M) + G(M) * 0.5d0 * dts
         end do
      end do
      return

   end subroutine thermo_NHC_prop

end module ThermoModule
