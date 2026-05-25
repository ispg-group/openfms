program tester
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
                        init_color_output, get_final_message
   use testutils, only: unit_test_error_handler, check_dieerror_not_called
   use test_particle, only: collect_particle_suite
   use test_trajectory, only: collect_trajectory_suite
   use test_bundle, only: collect_bundle_suite
   use test_restart, only: collect_restart_suite
   use GlobalModule, only: set_error_handler
   implicit none
   integer :: num_failed_tests, is, num_args, selected_suite_index
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'
   character(len=255) :: force_color, selected_suite

   num_failed_tests = 0

   num_args = command_argument_count()
   if (num_args > 1) then
      write (error_unit, '(a)') 'Usage: openfms_unit_tests [suite]'
      stop 1
   end if
   if (num_args == 1) call get_command_argument(1, selected_suite)

   ! Decide whether the output should be colourful:
   ! 1. Respect the FORCE_COLOR envvar per https://force-color.org/
   ! 2. Color things only if the output is a TTY (terminal)
   ! The latter is currently only checked using the GNU-instrinsic isatty
   call get_environment_variable('FORCE_COLOR', force_color)
   if (trim(force_color) /= '') then
      call init_color_output(.true.)
   end if
#ifdef __GNUC__
   if (isatty(output_unit) .and. isatty(error_unit)) then
      call init_color_output(.true.)
   end if
#endif

   testsuites = [ &
                new_testsuite('ParticleModule', collect_particle_suite), &
                new_testsuite('TrajectoryModule', collect_trajectory_suite), &
                new_testsuite('BundleModule', collect_bundle_suite), &
                new_testsuite('RestartModule', collect_restart_suite) &
                ]

   ! Swap the default FMS_DieError handler for a unit-test friendly
   ! error handler defined in testutils.F90.
   call set_error_handler(unit_test_error_handler)

   if (num_args == 1) then
      selected_suite_index = find_suite(testsuites, selected_suite)
      if (selected_suite_index == 0) then
         call print_available_suites(testsuites, selected_suite)
         stop 1
      end if

      call run_suite(testsuites(selected_suite_index), num_failed_tests)
   else
      do is = 1, size(testsuites)
         call run_suite(testsuites(is), num_failed_tests)
      end do
   end if

   write (error_unit, '(a)') trim(get_final_message(num_failed_tests))

   if (num_failed_tests > 0) then
      stop 1
   end if

contains

   integer function find_suite(testsuites, suite_name)
      type(testsuite_type), intent(in) :: testsuites(:)
      character(len=*), intent(in) :: suite_name
      integer :: is

      do is = 1, size(testsuites)
         if (trim(testsuites(is)%name) == trim(suite_name)) then
            find_suite = is
            return
         end if
      end do

      find_suite = 0
   end function find_suite

   subroutine print_available_suites(testsuites, suite_name)
      type(testsuite_type), intent(in) :: testsuites(:)
      character(len=*), intent(in) :: suite_name
      integer :: is

      write (error_unit, '(a)') 'Unknown unit test suite: '//trim(suite_name)
      write (error_unit, '(a)') 'Available suites:'
      do is = 1, size(testsuites)
         write (error_unit, '(2x, a)') testsuites(is)%name
      end do
   end subroutine print_available_suites

   subroutine run_suite(testsuite, num_failed_tests)
      type(testsuite_type), intent(in) :: testsuite
      integer, intent(inout) :: num_failed_tests

      write (error_unit, fmt) 'Testing:', testsuite%name
      call run_testsuite(testsuite%collect, error_unit, num_failed_tests)

      ! Make sure FMS_DieError was not called unexpectedly.
      call check_dieerror_not_called()
      write (error_unit, *)
   end subroutine run_suite

end program tester
