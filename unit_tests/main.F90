program tester
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use testdrive, only: run_testsuite, run_selected, new_testsuite, testsuite_type, &
                        init_color_output, get_final_message
   use testutils, only: unit_test_error_handler, check_dieerror_not_called
   use test_particle, only: collect_particle_suite
   use test_trajectory, only: collect_trajectory_suite
   use test_bundle, only: collect_bundle_suite
   use GlobalModule, only: set_error_handler
#ifndef __GNUC__
! Needed for isatty intrinsic
   use ifport, only: isatty
#endif
   implicit none
   integer :: num_failed_tests, is, num_args, selected_suite_index
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'
   character(len=255) :: force_color, selected_suite, selected_test

   num_failed_tests = 0
   selected_suite = ''
   selected_test = ''

   num_args = command_argument_count()
   if (num_args > 2) then
      write (error_unit, '(a)') 'Usage: openfms_unit_tests [suite [test]]'
      stop 1
   end if
   if (num_args >= 1) call get_command_argument(1, selected_suite)
   if (num_args >= 2) call get_command_argument(2, selected_test)

   ! Decide whether the output should be colourful:
   ! 1. Respect the FORCE_COLOR envvar per https://force-color.org/
   ! 2. Don't color things if the output is not a TTY (e.g. redirected to a file)
   call get_environment_variable("FORCE_COLOR", force_color)
   if (trim(force_color) /= '' .or. (isatty(output_unit) .and. isatty(error_unit))) then
      call init_color_output(.true.)
   end if

   testsuites = [ &
                new_testsuite("ParticleModule", collect_particle_suite), &
                new_testsuite("TrajectoryModule", collect_trajectory_suite), &
                new_testsuite("BundleModule", collect_bundle_suite) &
                ]

   ! Swap the default FMS_DieError handler for a unit-test friendly
   ! error handler defined in testutils.F90.
   call set_error_handler(unit_test_error_handler)

   if (num_args >= 1) then
      selected_suite_index = find_suite(testsuites, selected_suite)
      if (selected_suite_index == 0) then
         call print_available_suites(testsuites, selected_suite)
         stop 1
      end if

      call run_suite(testsuites(selected_suite_index), selected_test, num_args >= 2, num_failed_tests)
   else
      do is = 1, size(testsuites)
         call run_suite(testsuites(is), selected_test, .false., num_failed_tests)
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

   subroutine run_suite(testsuite, selected_test, only_selected_test, num_failed_tests)
      type(testsuite_type), intent(in) :: testsuite
      character(len=*), intent(in) :: selected_test
      logical, intent(in) :: only_selected_test
      integer, intent(inout) :: num_failed_tests

      write (error_unit, fmt) "Testing:", testsuite%name
      if (only_selected_test) then
         call run_selected(testsuite%collect, trim(selected_test), error_unit, num_failed_tests)
         if (num_failed_tests < 0) stop 1
      else
         call run_testsuite(testsuite%collect, error_unit, num_failed_tests)
      end if

      ! Make sure FMS_DieError was not called unexpectedly.
      call check_dieerror_not_called()
      write (error_unit, *)
   end subroutine run_suite

end program tester
