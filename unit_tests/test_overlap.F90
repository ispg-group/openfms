module test_overlap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt, fmiOut
   use OverlapModule, only: overlap, overlap_S_dot
   use TrajectoryModule
   implicit none
   private

   public :: collect_overlap_suite
   real(kind=DefReal), parameter :: atol = 1.d-10
   real(kind=DefReal), parameter :: rtol = 1.d-9

contains

!> Collect all exported unit tests
   subroutine collect_overlap_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest('overlap_S_dot', test_overlap_S_dot) &
                  ]

   end subroutine collect_overlap_suite

   function get_numerical_time_derivative(T1, T2, integral_type, nr_steps, step_size) result(num_dt)
      use VerletModule, only: FMS_PropVV
      type(T_Trajectory), intent(in) :: T1, T2
      character, intent(in) :: integral_type
      integer(kind=DefInt), intent(in) :: nr_steps
      real(kind=DefReal), intent(in) :: step_size
      real(kind=DefReal) :: step
      type(T_Trajectory) :: T1_tmp, T2_tmp
      complex(kind=DefReal) :: num_dt(nr_steps)
      complex(kind=DefReal) :: tmp_num_deriv, tmp
      integer(kind=DefInt) :: j_step, k

      num_dt = 0.d0
      step = step_size
      do j_step = 1, nr_steps
         tmp_num_deriv = 0.d0

         do k = 1, 2
            T1_tmp = T1
            T2_tmp = T2

            !call FMS_PropVV(T1_tmp, (-1.d0)**k * step_size)
            call FMS_PropVV(T2_tmp, (-1.d0)**k * step)

            select case (integral_type)
            case ('S')
               tmp = overlap(T1_tmp, T2_tmp)
            case default
               print*,'You really should not be seeing this error message!'
               stop 1
            end select
            tmp_num_deriv = tmp_num_deriv + (-1)**k * tmp
         end do

         num_dt = 0.5d0 * tmp_num_deriv / step
         step = step * step_size
      end do
   end function get_numerical_time_derivative

   subroutine test_overlap_S_dot(error)
      use GlobalModule, only: gliModel, gliMethod
      use TrajectoryCalcsModule, only: FMS_Forces
      use ToyModelModule, only: IzmaylovParams
      use ElecStrucModule
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: T1, T2
      integer(kind=DefInt) :: num_particles, num_states
      complex(kind=DefReal) :: SDot_analytical(7)
      complex(kind=DefReal) :: SDot_numerical(7)
      real(kind=DefReal) :: pos(6), mom(6), vel(6)
      real(kind=DefReal) :: mass
      real(kind=DefReal) :: min_err

      num_particles = 2
      num_states = 1
      gliModel = 0
      gliMethod = 3
      fmiOut = open_dev_null()
      call IzmaylovParams%initialize(W1=0.009557d0, W2=0.003515d0, XA=20.07d0, YA=0.0d0, &
                                     deltaE=0.01984d0, coupC=0.0006127d0)
      call FMS_ESInit(num_particles, num_states)
      call T1%create(numparticles=num_particles, numstates=num_states)
      call T2%create(numparticles=num_particles, numstates=num_states)

      mass = 0.0005485
      pos = [-9.8126, 0.0000, 0.0000, 14.2616, 0.0000, 0.0000]
      mom = [0.0938, 0.0000, 0.0000, 0.0035, 0.0000, 0.0000]
      vel = mom / mass
      call T1%set_pos(pos)
      call T1%set_mom(mom)
      call T1%set_vel(vel)
      T1%StateID = 1
      T1%Phase = 6.1910
      T1%Particle(:)%mass = mass
      T1%Particle(1)%width = 0.004775
      T1%Particle(2)%width = 0.001675

      pos = [-9.4015, 0.0000, 0.0000, 14.2493, 0.0000, 0.0000]
      mom = [0.2116, 0.0000, 0.0000, 0.0022, 0.0000, 0.0000]
      vel = mom / mass
      call T2%set_pos(pos)
      call T2%set_mom(mom)
      call T2%set_vel(vel)
      T2%StateID = 1
      T2%Phase = 0.0320
      T2%Particle(:)%mass = 0.0005485
      T2%Particle(1)%width = 0.004775
      T2%Particle(2)%width = 0.001675

      SDot_analytical = overlap_S_dot(T1, T2)
      ! calculates SDot numerically, (<chi_2(t)|chi_2(t+h)> - <chi_2(t)|chi_2(t-h)>)/2
      ! with seven different step sizes h = 1.d-1, 1.d-2, ..., 1.d-7
      SDot_numerical = get_numerical_time_derivative(T1, T2, 'S', 7, 1.d-1)
      ! since the numerical derivative suffers from catastrophic cancelations
      ! at small h, we look for the smallest error
      ! see, e.g., https://en.wikipedia.org/wiki/Numerical_differentiation#/media/File%3AAbsoluteErrorNumericalDifferentiationExample.png
      min_err = minval(abs(SDot_analytical - SDot_numerical))
      call check(error, min_err < atol + rtol * abs(SDot_analytical(1)), 'Analytical SDot ' &
                 //'significantly from numerical one')
      if (allocated(error)) return

      call T1%destroy()
      call T2%destroy()
   end subroutine test_overlap_S_dot

   function open_dev_null() result(output_unit)
      integer :: output_unit

      open (newunit=output_unit, file='/dev/null', action='write')
   end function open_dev_null

end module test_overlap
