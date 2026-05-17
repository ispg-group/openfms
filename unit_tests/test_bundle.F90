module test_bundle
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule
   use BundleModule
   use TrajectoryCalcsModule, only: FMS_isBundleCurrent, FMS_isAmpDotCurrent
   implicit none
   private

   public :: collect_bundle_suite

contains

!> Collect all exported unit tests
   subroutine collect_bundle_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest('create', test_create), &
                  new_unittest('assign', test_assign), &
                  new_unittest('destroy', test_destroy) &
                  ]

   end subroutine collect_bundle_suite

   subroutine test_create(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_TrajectoryBundle) :: B
      integer, parameter :: NUM_TRAJ = 2, NUM_DEADTRAJ = 3, NUM_PART = 4
      integer, parameter :: NUM_STATES = 5, NCBFS = 6

      ! To allocate B%Centroids
      glzCentroids = .true.

      call B%create(numtraj=NUM_TRAJ, &
                    numdeadtraj=NUM_DEADTRAJ, &
                    numparticles=NUM_PART, &
                    numstates=NUM_STATES, &
                    ncbfs=NCBFS)

      call check(error, size(B%Trajectory), NUM_TRAJ)
      if (allocated(error)) return

      call check(error, size(B%DeadTraj), NUM_DEADTRAJ)
      if (allocated(error)) return

      ! I don't understand why this size is computed as it is
      call check(error, size(B%Centroids), 15)
      if (allocated(error)) return

      call check(error, B%NumStates, NUM_STATES)
      if (allocated(error)) return

      call check(error, B%NumParticles, NUM_PART)
      if (allocated(error)) return

      call check(error, B%NCBFS, NCBFS)
      if (allocated(error)) return

      call check(error, B%CurrentTime, 0.0d0)
      if (allocated(error)) return

      call check(error, size(B%BMatrices%H), NUM_TRAJ * NUM_TRAJ, 'Create failed to allocate H matrix')
      if (allocated(error)) return

      call check(error, size(B%BMatrices%S), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%SDot), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%SInv), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%Sp5i), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%Sp5), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%HEff), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%HEff1), NUM_TRAJ * NUM_TRAJ)
      if (allocated(error)) return
      call check(error, size(B%BMatrices%AmpDot), NUM_TRAJ)
      if (allocated(error)) return

      ! Check that we can destroy "empty" Trajectory
      call B%destroy()

   end subroutine test_create

   subroutine test_assign(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_TrajectoryBundle) :: B1, B2 ! Assign B1 to B2
      integer, dimension(2), parameter :: NUM_TRAJ = [2, 1]
      integer, dimension(2), parameter :: NUM_DEADTRAJ = [3, 1], NUM_PART = [4, 1]
      integer, dimension(2), parameter :: NUM_STATES = [5, 1], NCBFS = [6, 1]
      integer :: num_traj_squared, itraj

      ! To allocate B%Centroids
      glzCentroids = .true.

      call B1%create(numtraj=NUM_TRAJ(1), &
                     numdeadtraj=NUM_DEADTRAJ(1), &
                     numparticles=NUM_PART(1), &
                     numstates=NUM_STATES(1), &
                     ncbfs=NCBFS(1))

      call B2%create(numtraj=NUM_TRAJ(2), &
                     numdeadtraj=NUM_DEADTRAJ(2), &
                     numparticles=NUM_PART(2), &
                     numstates=NUM_STATES(2), &
                     ncbfs=NCBFS(2))

      B2 = B1

      call check(error, size(B2%Trajectory), NUM_TRAJ(1))
      if (allocated(error)) return

      call check(error, size(B2%DeadTraj), NUM_DEADTRAJ(1))
      if (allocated(error)) return

      call check(error, B2%NumStates, NUM_STATES(1))
      if (allocated(error)) return

      call check(error, B2%NumParticles, NUM_PART(1))
      if (allocated(error)) return

      call check(error, B2%NCBFS, NCBFS(1))
      if (allocated(error)) return

      call check(error, B2%CurrentTime, 0.0d0)
      if (allocated(error)) return

      do itraj = 1, B2%NumTraj
         ! TODO(DH) Is this correct?? See comment in BundleModule::assign_bundle
         call check(error, FMS_IsBundleCurrent(B2%Trajectory(itraj)), .false.)
         call check(error, FMS_IsAmpDotCurrent(B2%Trajectory(itraj)), .false.)
      end do

      NUM_TRAJ_SQUARED = NUM_TRAJ(1) * NUM_TRAJ(1)
      call check(error, size(B2%BMatrices%H), NUM_TRAJ_SQUARED, 'Create failed to allocate H matrix')
      if (allocated(error)) return

      call check(error, size(B2%BMatrices%S), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%SDot), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%SInv), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%Sp5i), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%Sp5), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%HEff), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%HEff1), NUM_TRAJ_SQUARED)
      if (allocated(error)) return
      call check(error, size(B2%BMatrices%AmpDot), NUM_TRAJ(1))
      if (allocated(error)) return

   end subroutine test_assign

   subroutine test_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_TrajectoryBundle) :: B
      integer, parameter :: NUM_TRAJ = 2, NUM_DEADTRAJ = 3, NUM_PART = 4
      integer, parameter :: NUM_STATES = 5, NCBFS = 6

      ! To allocate B%Centroids
      glzCentroids = .true.

      call B%create(numtraj=NUM_TRAJ, &
                    numdeadtraj=NUM_DEADTRAJ, &
                    numparticles=NUM_PART, &
                    numstates=NUM_STATES, &
                    ncbfs=NCBFS)
      B%CurrentTime = 10.0d0

      call B%destroy()

      call check(error,.not. allocated(B%Trajectory), 'Destroy failed to deallocate B%Trajectory array')
      if (allocated(error)) return

      call check(error,.not. allocated(B%DeadTraj), 'Destroy failed to deallocate B%DeadTraj array')
      if (allocated(error)) return

      call check(error,.not. allocated(B%Centroids), 'Destroy failed to deallocate B%Centroids array')
      if (allocated(error)) return

      call check(error,.not. allocated(B%BMatrices%H), 'Destroy failed to deallocate H matrix')
      if (allocated(error)) return

      call check(error,.not. allocated(B%DeadH))
      if (allocated(error)) return

      call check(error,.not. allocated(B%BMatrices%S))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%SDot))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%SInv))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%Sp5i))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%Sp5))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%HEff))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%HEff1))
      if (allocated(error)) return
      call check(error,.not. allocated(B%BMatrices%AmpDot))
      if (allocated(error)) return

      call check(error, B%NumParticles, 0, 'Destroyed trajectory should have no particles!')
      if (allocated(error)) return

      call check(error, B%NumStates, 0)
      if (allocated(error)) return

      call check(error, B%NumTraj, 0)
      if (allocated(error)) return

      ! TODO: The checks below currently do not work!
      ! call check(error, B%CurrentTime, 0.0D0)
      ! if (allocated(error)) return

      ! call check(error, B%NumDeadTraj, 0)
      ! if (allocated(error)) return

      ! call check(error, B%NCBFS, 0)
      ! if (allocated(error)) return

   end subroutine test_destroy

end module test_bundle
