module test_overlap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt
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

   function get_numerical_time_derivative(T1, T2, integral_type) result(num_dt)
      use VerletModule, only: FMS_PropVV
      type(T_Trajectory), intent(in) :: T1, T2
      character, intent(in) :: integral_type
      type(T_Trajectory) :: T1_tmp, T2_tmp
      complex(kind=DefReal) :: num_dt(7)
      real(kind=DefReal) :: step_size
      complex(kind=DefReal) :: tmp_num_deriv, tmp
      integer(kind=DefInt) :: i_direction, j_step, nr_steps, k, l

      step_size = 1.d-1
      num_dt = 0.d0
      nr_steps = 7
      do j_step = 1, nr_steps
         tmp_num_deriv = 0.d0

         do k = 1, 2
            T1_tmp = T1
            T2_tmp = T2

            !call FMS_PropVV(T1_tmp, (-1.d0)**k * step_size)
            call FMS_PropVV(T2_tmp, (-1.d0)**k * step_size)

            select case (integral_type)
            case ('S')
               tmp = overlap(T1_tmp, T2_tmp)
            case default
               return
            end select
            tmp_num_deriv = tmp_num_deriv + (-1)**k * tmp
         end do

         num_dt = 0.5d0 * tmp_num_deriv / step_size
         step_size = step_size * 1.d-1
      end do
   end function get_numerical_time_derivative

   subroutine test_overlap_S_dot(error)
      use GlobalModule, only: gliModel, gliMethod, glIzmOmegax, glIzmOmegay, &
                              glIzmXshift, glIzmYshift, glIzmDeltaE, glIzmCoupC
      use TrajectoryCalcsModule, only: FMS_Forces
      use ElecStrucModule
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: T1, T2
      integer(kind=DefInt) :: num_particles, num_states
      complex(kind=DefReal) :: SDot_analytical(7)
      complex(kind=DefReal) :: SDot_numerical(7)
      real(kind=DefReal) :: pos(6), mom(6), vel(6), force(6)
      real(kind=DefReal) :: mass
      real(kind=DefReal) :: min_err

      num_particles = 2
      num_states = 1
      gliModel = 0
      gliMethod = 3
      glIzmOmegax = 0.009557
      glIzmOmegay = 0.003515
      glIzmXshift = 20.07
      glIzmYshift = 0.0
      glIzmDeltaE = 0.01984
      glIzmCoupC = 0.0006127
      call FMS_ESInit(num_particles, num_states, suppress_write=.true.)
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
      SDot_numerical = get_numerical_time_derivative(T1, T2, 'S')
      ! since the numerical derivative suffers from cancelations at small h,
      ! we look for the smallest
      min_err = minval(abs(SDot_analytical - SDot_numerical))
      call check(error, min_err < atol + rtol * abs(SDot_analytical(1)), 'Analytical SDot ' &
                 //'significantly from numerical one')
      if (allocated(error)) return

      call T1%destroy()
      call T2%destroy()
   end subroutine test_overlap_S_dot

end module test_overlap
