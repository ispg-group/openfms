module test_restart
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use GlobalModule
   use BundleModule
   use ElecStrucModule, only: esBlobSize, esLCIVec, esNBasis, esnElecPhase
   use RestartModule
   implicit none
   private

   public :: collect_restart_suite

contains

!> Collect all exported unit tests
   subroutine collect_restart_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest('read_header_ci_orbital_sizes', test_read_header_ci_orbital_sizes), &
                  new_unittest('use_known_ci_orbital_sizes', test_use_known_ci_orbital_sizes) &
                  ]

   end subroutine collect_restart_suite

   subroutine test_read_header_ci_orbital_sizes(error)
      type(error_type), allocatable, intent(out) :: error

      call exercise_restart_size_path(error, strip_dimensions=.false., use_known_sizes=.false.)
   end subroutine test_read_header_ci_orbital_sizes

   subroutine test_use_known_ci_orbital_sizes(error)
      type(error_type), allocatable, intent(out) :: error

      call exercise_restart_size_path(error, strip_dimensions=.true., use_known_sizes=.true.)
   end subroutine test_use_known_ci_orbital_sizes

   subroutine exercise_restart_size_path(error, strip_dimensions, use_known_sizes)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: strip_dimensions
      logical, intent(in) :: use_known_sizes
      type(T_TrajectoryBundle) :: written_bundle, restarted_bundle
      real(kind=DefReal), dimension(4, 4) :: orbitals
      real(kind=DefReal), dimension(2, 3) :: ci_vectors
      integer :: i

      orbitals = reshape([(real(i, kind=DefReal), i=1, 16)], shape(orbitals))
      ci_vectors = reshape([(real(100 + i, kind=DefReal), i=1, 6)], shape(ci_vectors))

      call configure_restart_test_globals()
      call FMS_DeleteFile('Last_Bundle.txt')
      call FMS_DeleteFile('tmp.Last_Bundle.txt')
      call FMS_DeleteFile('Checkpoint.txt')

      esNBasis = size(orbitals, dim=1)
      esLCIVec = size(ci_vectors, dim=2)
      esBlobSize = 1
      esnElecPhase = 2

      call written_bundle%create(numtraj=1, numdeadtraj=0, numparticles=1, numstates=2, ncbfs=1)
      written_bundle%CurrentTime = 1.0_DefReal
      written_bundle%Trajectory(1)%TrajID = 1
      written_bundle%Trajectory(1)%CBF = 1
      written_bundle%Trajectory(1)%StateID = 1
      written_bundle%Trajectory(1)%Ms = 1
      written_bundle%Trajectory(1)%ElecStruc%OldOrbitals = orbitals
      written_bundle%Trajectory(1)%ElecStruc%OldCIVecs = ci_vectors
      written_bundle%Trajectory(1)%ElecStruc%ElecPhase = 1.0_DefReal

      call PutRestart(written_bundle)

      call check_sized_restart_header(error, '# Orbitals', size(orbitals, dim=1), size(orbitals, dim=2))
      if (allocated(error)) return
      call check_sized_restart_header(error, '# CI vectors', size(ci_vectors, dim=1), size(ci_vectors, dim=2))
      if (allocated(error)) return

      if (strip_dimensions) call strip_restart_matrix_dimensions()

      if (use_known_sizes) then
         esNBasis = size(orbitals, dim=1)
         esLCIVec = size(ci_vectors, dim=2)
      else
         esNBasis = 0
         esLCIVec = 0
      end if
      esBlobSize = 1
      call GetRestart(restarted_bundle, -1.0_DefReal)

      call check(error, allocated(restarted_bundle%Trajectory(1)%ElecStruc%OldOrbitals), &
                 'Restart did not allocate old orbitals')
      if (allocated(error)) return
      call check(error, allocated(restarted_bundle%Trajectory(1)%ElecStruc%OldCIVecs), &
                 'Restart did not allocate CI vectors')
      if (allocated(error)) return

      call check(error, size(restarted_bundle%Trajectory(1)%ElecStruc%OldOrbitals, dim=1), size(orbitals, dim=1))
      if (allocated(error)) return
      call check(error, size(restarted_bundle%Trajectory(1)%ElecStruc%OldOrbitals, dim=2), size(orbitals, dim=2))
      if (allocated(error)) return
      call check(error, size(restarted_bundle%Trajectory(1)%ElecStruc%OldCIVecs, dim=1), size(ci_vectors, dim=1))
      if (allocated(error)) return
      call check(error, size(restarted_bundle%Trajectory(1)%ElecStruc%OldCIVecs, dim=2), size(ci_vectors, dim=2))
      if (allocated(error)) return

      call check(error, all(restarted_bundle%Trajectory(1)%ElecStruc%OldOrbitals == orbitals), &
                 'Restart changed old orbital values')
      if (allocated(error)) return
      call check(error, all(restarted_bundle%Trajectory(1)%ElecStruc%OldCIVecs == ci_vectors), &
                 'Restart changed CI vector values')
      if (allocated(error)) return

      call check(error, esNBasis, size(orbitals, dim=1))
      if (allocated(error)) return
      call check(error, esLCIVec, size(ci_vectors, dim=2))

      call written_bundle%destroy()
      call restarted_bundle%destroy()
      call FMS_DeleteFile('Last_Bundle.txt')
      call FMS_DeleteFile('tmp.Last_Bundle.txt')
      call FMS_DeleteFile('Checkpoint.txt')
   end subroutine exercise_restart_size_path

   subroutine configure_restart_test_globals()
      character(len=256) :: test_tmpdir

      call get_environment_variable('OPENFMS_TEST_TMPDIR', test_tmpdir)
      if (len_trim(test_tmpdir) == 0) test_tmpdir = '/private/tmp/openfms-restart-unit'
      call execute_command_line('mkdir -p '//trim(test_tmpdir))

      fmiOut = error_unit
      FMSWorkingDir = trim(test_tmpdir)//'/'
      glzCentroids = .false.
      glzStoSwiss = .false.
      gldTimeStep = 1.0_DefReal
      inIRestartTraj = 0
      NSing = 2
      RestartStep = 100000
      zRedoRestartES = .false.
   end subroutine configure_restart_test_globals

   subroutine check_sized_restart_header(error, label, expected_rows, expected_cols)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer, intent(in) :: expected_rows, expected_cols

      character(len=512) :: file_name
      character(len=500) :: line
      character(len=500) :: adjusted_line
      integer :: unit, ios, rows, cols
      logical :: found

      found = .false.
      file_name = trim(FMSWorkingDir)//'Last_Bundle.txt'
      open (newunit=unit, file=file_name, status='old', action='read')
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         adjusted_line = adjustl(line)
         if (index(adjusted_line, trim(label)) == 1) then
            read (adjusted_line(len_trim(label) + 1:), *, iostat=ios) rows, cols
            found = ios == 0 .and. rows == expected_rows .and. cols == expected_cols
            exit
         end if
      end do
      close (unit)

      call check(error, found, trim(label)//' restart header did not include expected dimensions')
   end subroutine check_sized_restart_header

   subroutine strip_restart_matrix_dimensions()
      character(len=500), allocatable :: lines(:)
      character(len=500) :: line
      character(len=512) :: file_name
      integer :: unit, ios, nlines, capacity, i

      capacity = 128
      nlines = 0
      allocate (lines(capacity))

      file_name = trim(FMSWorkingDir)//'Last_Bundle.txt'
      open (newunit=unit, file=file_name, status='old', action='read')
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (nlines == capacity) call grow_line_buffer(lines, capacity)
         nlines = nlines + 1
         if (index(adjustl(line), '# Orbitals') == 1) then
            lines(nlines) = ' # Orbitals'
         else if (index(adjustl(line), '# CI vectors') == 1) then
            lines(nlines) = ' # CI vectors'
         else
            lines(nlines) = line
         end if
      end do
      close (unit)

      open (newunit=unit, file=file_name, status='replace', action='write')
      do i = 1, nlines
         write (unit, '(A)') trim(lines(i))
      end do
      close (unit)
   end subroutine strip_restart_matrix_dimensions

   subroutine grow_line_buffer(lines, capacity)
      character(len=500), allocatable, intent(inout) :: lines(:)
      integer, intent(inout) :: capacity
      character(len=500), allocatable :: tmp(:)

      allocate (tmp(capacity * 2))
      tmp(1:capacity) = lines
      call move_alloc(tmp, lines)
      capacity = size(lines)
   end subroutine grow_line_buffer

end module test_restart
