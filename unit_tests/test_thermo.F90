module test_thermo
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use ThermoModule, only: thermo_NHC_local
   implicit none
   private

   public :: collect_thermo_suite

contains

   subroutine collect_thermo_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest('nhc_ho', test_nhc_ho)]
   end subroutine collect_thermo_suite

   subroutine test_nhc_ho(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: chain_length = 2, number_of_steps = 5
      real(8), parameter :: particle_mass = 1.0d0, omega = 1.0d0
      real(8), parameter :: thermostat_mass = 1.0d0, beta = 1.0d0
      real(8), parameter :: timestep = 0.005d0
      real(8), parameter :: reference_tolerance = 1.0d-11
      real(8) :: position, momentum, chain_position(chain_length)
      real(8) :: chain_momentum(chain_length), chain_mass(chain_length)
      real(8) :: thermostat_energy
      real(8) :: force, max_position_error, max_momentum_error, max_thermostat_energy_error
      integer :: chain_index, step
      logical :: finite_state

      real(8), parameter :: ref_position(number_of_steps + 1) = [1.0d0, &
                                                                 1.0049750312109134d0, 1.0099003120142305d0, &
                                                                 1.0147760913640445d0, 1.0196026174098716d0, &
                                                                 1.0243801374121491d0]

      real(8), parameter :: ref_momentum(number_of_steps + 1) = [1.0d0, &
                                                                 9.9002499910751862d-1, 9.8009986223049894d-1, &
                                                                 9.7022443718255169d-1, 9.6039855451514111d-1, &
                                                                 9.5062202823696751d-1]

      real(8), parameter :: ref_thermostat_energy(number_of_steps + 1) = [1.0d0, &
                                                                          1.0049378750636817d0, 1.0097528721079072d0, &
                                                                          1.0144471059777167d0, 1.0190226832452651d0, &
                                                                          1.0234817009863844d0]

      position = 1.0d0
      momentum = 1.0d0
      chain_position = 0.0d0
      chain_momentum = 1.0d0
      chain_mass = thermostat_mass
      thermostat_energy = 0.0d0

      do chain_index = 1, chain_length
         thermostat_energy = thermostat_energy + 0.5d0 * chain_momentum(chain_index)**2 / chain_mass(chain_index) &
                             + chain_position(chain_index) / beta
      end do

      max_position_error = abs(position - ref_position(1))
      max_momentum_error = abs(momentum - ref_momentum(1))
      max_thermostat_energy_error = abs(thermostat_energy - ref_thermostat_energy(1))
      finite_state = .true.

      do step = 1, number_of_steps

         ! --- O (dt/2): NHC half-kick
         call thermo_NHC_local(0.5d0 * timestep, chain_length, momentum, particle_mass, beta, &
                               chain_position, chain_momentum, chain_mass, thermostat_energy)

         ! --- B (dt/2): momentum half-kick
         force = -particle_mass * omega**2 * position
         momentum = momentum + 0.5d0 * timestep * force
         ! --- A (dt): position full-kick
         position = position + timestep * momentum / particle_mass
         force = -particle_mass * omega**2 * position
         ! --- B (dt/2): momentum half-kick
         momentum = momentum + 0.5d0 * timestep * force

         ! --- O (dt/2): NHC half-kick
         call thermo_NHC_local(0.5d0 * timestep, chain_length, momentum, particle_mass, beta, &
                               chain_position, chain_momentum, chain_mass, thermostat_energy)

         max_position_error = max(max_position_error, abs(position - ref_position(step + 1)))
         max_momentum_error = max(max_momentum_error, abs(momentum - ref_momentum(step + 1)))
         max_thermostat_energy_error = max(max_thermostat_energy_error, &
                                           abs(thermostat_energy - ref_thermostat_energy(step + 1)))

      end do

      call check(error, finite_state, 'NHC oscillator state became non-finite')
      if (allocated(error)) return
      call check(error, max_position_error < reference_tolerance, &
                 'NHC position does not match reference values')
      if (allocated(error)) return
      call check(error, max_momentum_error < reference_tolerance, &
                 'NHC momentum does not match reference values')
      if (allocated(error)) return
      call check(error, max_thermostat_energy_error < reference_tolerance, &
                 'NHC thermostat energy does not match reference values')
   end subroutine test_nhc_ho

end module test_thermo
