module testutils
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: error_type, check
   implicit none
   private

   public :: check_dieerror_called, check_dieerror_not_called
   public :: unit_test_error_handler

   ! A global private flag that indicates whether FMS_DieError has been called.
   ! during unit tests. The error message is stored in 'dieerror_msg'.
   ! The flag and the error message are checked by the 'check_dieerror_called' subroutine.
   logical, private :: dieerror_called = .false.
   character(len=:), private, allocatable :: dieerror_msg
   character(len=*), parameter :: newline = new_line("a")

contains

!> Monkey patch for FMS_DieError so that it does not stop the test suite.
   subroutine unit_test_error_handler(message)
      character(len=*), intent(in) :: message

      ! If the 'dieerror_called' flag is already set, it means there was
      ! an unhandled DieError call, and we should stop immediately.
      call check_dieerror_not_called()

      dieerror_msg = trim(message)
      dieerror_called = .true.
   end subroutine unit_test_error_handler

!> Call this if you expect DieError to be have been called
   subroutine check_dieerror_called(error, expected_msg)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: expected_msg

      call check(error, dieerror_called, "DieError not called", "expected error: "//expected_msg)
      if (allocated(error)) then
         call reset_unit_test_error()
         return
      end if

      call check(error, expected_msg, dieerror_msg, more="Unexpected error message in FMS_DieError()")

      call reset_unit_test_error()
   end subroutine check_dieerror_called

!> Subroutine to verify that DieError was not called unexpectedly.
!> This is called automatically at the end of each test suite in main.F90.
!> Because the failure here likely means an error in the test suite itself,
!> we stop immediately and don't run any more tests.
   subroutine check_dieerror_not_called()
      if (dieerror_called) then
         write (error_unit, '(2x, a)') "Uncaught DieError call"
         if (allocated(dieerror_msg)) then
            write (error_unit, '(7x, a)') "error message: "//dieerror_msg
         else
            write (error_unit, '(7x, a)') "no error message allocated (weird!)"
         end if
         error stop 1
      end if
   end subroutine check_dieerror_not_called

   subroutine reset_unit_test_error()
      dieerror_called = .false.
      if (allocated(dieerror_msg)) deallocate (dieerror_msg)
   end subroutine reset_unit_test_error

end module testutils
