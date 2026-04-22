program tester
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
                        init_color_output, get_final_message
   use testutils, only: unit_test_error_handler, check_dieerror_not_called
   use test_particle, only: collect_particle_suite
   use test_trajectory, only: collect_trajectory_suite
   use test_bundle, only: collect_bundle_suite
   use GlobalModule, only: set_error_handler
   implicit none
   integer :: num_failed_tests, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'
   character(len=255) :: force_color

   num_failed_tests = 0
   ! Decide whether the output should be colourful:
   ! 1. Respect the FORCE_COLOR envvar per https://force-color.org/
   ! 2. Color things only if the output is a TTY (terminal)
   ! The latter is currently only checked using the GNU-instrinsic isatty
   call get_environment_variable("FORCE_COLOR", force_color)
   if (trim(force_color) /= '') then
      call init_color_output(.true.)
   end if
#ifdef __GNUC__
   if (isatty(output_unit) .and. isatty(error_unit)) then
      call init_color_output(.true.)
   end if
#endif

   testsuites = [ &
                new_testsuite("ParticleModule", collect_particle_suite), &
                new_testsuite("TrajectoryModule", collect_trajectory_suite), &
                new_testsuite("BundleModule", collect_bundle_suite) &
                ]

   ! Swap the default FMS_DieError handler for a unit-test friendly
   ! error handler defined in testutils.F90.
   call set_error_handler(unit_test_error_handler)

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, num_failed_tests)
      ! Make sure FMS_DieError was not called unexpectedly.
      call check_dieerror_not_called()
      write (error_unit, *)
   end do

   write (error_unit, '(a)') trim(get_final_message(num_failed_tests))

   if (num_failed_tests > 0) then
      stop 1
   end if

end program tester
