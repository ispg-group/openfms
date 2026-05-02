# OpenFMS unit tests

Run from top-level directory as:

```console
make unittest
```

Each test suite should typically correspond to a single OpenFMS module,
and should live in a separate `test_<suite>.F90` file.

> [!IMPORTANT]
> If you don't adhere to this naming style, the new files will not be picked up in Makefile.

All new test suites also need to be registered in `main.F90` and in `cmake/OpenFMSUnitTests.cmake`.

The rest of this README is copied over from the
[test-drive](https://github.com/fortran-lang/test-drive) framework that we utilize here,
and serves as a guide how to write new tests.

## test-drive: The simple testing framework


This project offers a lightweight, procedural unit testing framework based on nothing but standard Fortran.
Integration with [meson](https://mesonbuild.com), [cmake](https://cmake.org) and [Fortran package manager (fpm)](https://github.com/fortran-lang/fpm) is available.
Alternatively, the [``testdrive.F90``](src/testdrive.F90) source file can be redistributed in the project's testsuite as well.


### Usage

Testsuites are defined by a ``collect_interface`` returning a set of ``unittest_type`` objects.
To create a new test use the ``new_unittest`` constructor, which requires a test identifier and a procedure with a ``test_interface`` compatible signature.
The error status is communicated by the allocation status of an ``error_type``.

The necessary boilerplate code to setup the test entry point is just

```fortran
program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_suite1, only : collect_suite1
  use test_suite2, only : collect_suite2
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("suite1", collect_suite1), &
    new_testsuite("suite2", collect_suite2) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester
```

Every test is defined in a separate module using a ``collect`` function, which is exported and added to the ``testsuites`` array in the test runner.
All tests have a simple interface with just an allocatable ``error_type`` as output to provide the test results.

```fortran
module test_suite1
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_suite1

contains

!> Collect all exported unit tests
subroutine collect_suite1(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("valid", test_valid), &
    new_unittest("invalid", test_invalid, should_fail=.true.) &
    ]

end subroutine collect_suite1

subroutine test_valid(error)
  type(error_type), allocatable, intent(out) :: error
  ! ...
end subroutine test_valid

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  ! ...
end subroutine test_invalid

end module test_suite1
```


### Checking test conditions

The procedures defining the tests can contain any Fortran code required for checking the correctness of the project.
An easy way to do so is provided by the generic ``check`` function.

```f90
subroutine test_valid(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, 1 + 2 == 3)
  if (allocated(error)) return

  ! equivalent to the above
  call check(error, 1 + 2, 3)
  if (allocated(error)) return
end subroutine test_valid
```

After each check, the status of the error should be checked.
Uncaught errors will not be silently dropped, instead the error will be caught, its message displayed and the run aborted.
Possible ways to use check are listed below

| available checks     | arguments                                                      |
| -------------------- | -------------------------------------------------------------- |
| logical check        | *error*, *logical*, ...                                        |
| status check         | *error*, *integer*, ...                                        |
| logical comparison   | *error*, *logical*, *logical*, ...                             |
| integer comparison   | *error*, *integer*, ...                                        |
| character comparison | *error*, *character*, *character*, ...                         |
| real comparison      | *error*, *real*, *real*, ..., thr=*real*, rel=*logical*        |
| real NaN check       | *error*, *real*, ...                                           |
| complex comparison   | *error*, *complex*, *complex*, ..., thr=*real*, rel=*logical*  |
| complex NaN check    | *error*, *complex*, ...                                        |

Each check will generate a meaningful error message based on the available arguments, but can also be provided with a custom error message instead.

To generate custom checks the ``test_failed`` procedure is available to generate error messages

```f90
subroutine test_custom(error)
  type(error_type), allocatable, intent(out) :: error

  ! ...

  if (.not.cond) then
    call test_failed(error, "Custom check failed")
    return
  end if

  ! ...

  if (.not.cond) then
    call test_failed(error, "Custom check failed", "Additional context")
    return
  end if

end subroutine test_custom
```

To conditionally skip a test use the ``skip_test`` procedure.
It uses the same signature as ``test_failed``, but will mark the respective test as skipped, this is useful to disable tests based on conditional compilation, *e.g.* by using a preprocessor or a different submodule.
An uncaught skipped test will fail regardless, therefore make sure to not run any other checks afterwards.


### License

This project is free software: you can redistribute it and/or modify it under the terms of the [Apache License, Version 2.0](LICENSE-Apache) or [MIT license](LICENSE-MIT) at your opinion.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an _as is_ basis, without warranties or conditions of any kind, either express or implied. See the License for the specific language governing permissions and limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this project by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
