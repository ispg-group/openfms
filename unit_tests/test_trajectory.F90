module test_trajectory
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use testutils, only: check_dieerror_called
  use GlobalModule
  use TrajectoryModule
  implicit none
  private

  public :: collect_trajectory_suite

contains

!> Collect all exported unit tests
subroutine collect_trajectory_suite(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("create", test_create), &
    new_unittest("create_errors", test_create_errors), &
    new_unittest("setget_time", test_setget_time), &
    new_unittest("get_mass", test_get_mass), &
    new_unittest("get_width", test_get_width), &
    new_unittest("setget_pos_comp", test_setget_pos_comp), &
    new_unittest("setget_mom_comp", test_setget_mom_comp), &
    new_unittest("setget_vec", test_setget_vec), &
    new_unittest("setget_all", test_setget_all) &
  ]

end subroutine collect_trajectory_suite

subroutine test_create(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=100, num_states=5

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, T%NumStates == num_states, "Invalid number of states")
  if (allocated(error)) return

  call check(error, size(T%Particle) == num_particles, "Invalid number of particles")
  if (allocated(error)) return
  call check(error, T%NumParticles == num_particles, "Invalid number of particles")
  if (allocated(error)) return

  ! Currently the number of dimensions is harcoded to 3 in FMS_Create_Trajectory
  call check(error, T%Particle(1)%NumDimensions == 3, "Invalid particle dimensions")
  if (allocated(error)) return

  call check(error, T%NumDimensions == 3 * num_particles, "Invalid trajectory dimension")
  if (allocated(error)) return

  ! Check that we can destroy "empty" Trajectory
  call T%destroy()

  call check(error, .not. allocated(T%Particle), "Destroy failed to deallocate T%Particle array")
  if (allocated(error)) return

  call check(error, T%NumParticles == 0, "Destroyed trajectory should have no particles!")
  if (allocated(error)) return

end subroutine test_create


subroutine test_create_errors(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T

  call T%create(numparticles=0, numstates=1)

  call check_dieerror_called(error, "CreateTrajectory: NumParticles must be > 0")

  call T%create(numparticles=1, numstates=-1)

  call check_dieerror_called(error, "CreateTrajectory: NumStates must be > 0")
  if (allocated(error)) return
end subroutine test_create_errors


subroutine test_setget_time(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=1, num_states=3

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, T%get_time() == 0.0D0, "Initial time is not 0")
  if (allocated(error)) return

  call T%set_time(10.0D0)
  call check(error, T%get_time() == 10.D0, "Failed to set current time")
  if (allocated(error)) return

end subroutine test_setget_time


subroutine test_get_mass(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3
  integer, parameter :: traj_dim = num_particles * 3
  integer :: i

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, all(T%get_mass() == [(0.0D0, i = 1, traj_dim)]), "Initial masses should be 0")
  if (allocated(error)) return

  T%Particle(2)%mass = 2.0D0

  ! Note: T%get_mass() returns a npart*numdim array due convenience reasons
  call check(error, all(T%get_mass() == [0.D0, 0.D0, 0.D0, 2.D0, 2.D0, 2.D0]), "Some masses should be non-zero!")
  if (allocated(error)) return

  call check(error, T%get_mass(1) == 0.0D0, "mass of particle 1 should still be 0")
  if (allocated(error)) return

  call check(error, T%get_mass(2) == 2.0D0, "mass of particle 2 should be non-zero")
  if (allocated(error)) return
end subroutine test_get_mass

subroutine test_get_width(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3
  integer, parameter :: traj_dim = num_particles * 3
  integer :: i

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, all(T%get_width() == [(0.0D0, i = 1, traj_dim)]), "Initial widths should be 0")
  if (allocated(error)) return

  T%Particle(2)%width = 2.0D0
  call check(error, all(T%get_width() == [0.D0, 0.D0, 0.D0, 2.D0, 2.D0, 2.D0]), "Some widths should be non-zero!")
  if (allocated(error)) return

  call check(error, T%get_width(1) == 0.0D0, "width of particle 1 should still be 0")
  if (allocated(error)) return

  call check(error, T%get_width(2) == 2.0D0, "width of particle 2 should be non-zero")
  if (allocated(error)) return
end subroutine test_get_width


subroutine test_setget_pos_comp(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3
  integer, parameter :: traj_dim = num_particles * num_states
  integer :: ipart, idim, i, j

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, all(T%get_pos() == [(0.0D0, i=1,traj_dim)]), "Initial positions should be 0")
  if (allocated(error)) return

  ipart = 2
  idim = 3
  call T%set_pos(ipart, idim, 1.0D0)

  do i = 1, num_particles
    do j = 1, 3
      if (i == ipart .and. j == idim) then
        call check(error, T%get_pos(i, j) == 1.0D0, "Expected position 1.0")
      else
        call check(error, T%get_pos(i, j) == 0.0D0, "Expected position 0.0")
      end if
      if (allocated(error)) return
    end do
  end do

end subroutine test_setget_pos_comp

subroutine test_setget_mom_comp(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3
  integer, parameter :: traj_dim = num_particles * num_states
  integer :: ipart, idim, i, j

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, all(T%get_mom() == [(0.0D0, i=1,traj_dim)]), "Initial positions should be 0")
  if (allocated(error)) return

  ipart = 2
  idim = 3
  call T%set_mom(ipart, idim, 1.0D0)

  do i = 1, num_particles
    do j = 1, 3
      if (i == ipart .and. j == idim) then
        call check(error, T%get_mom(i, j) == 1.0D0, "Expected momentum 1.0")
      else
        call check(error, T%get_mom(i, j) == 0.0D0, "Expected momentum 0.0")
      end if
      if (allocated(error)) return
    end do
  end do

end subroutine test_setget_mom_comp

subroutine test_setget_vec(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3
  integer, parameter :: traj_dim = num_particles * num_states
  integer :: i, ipart

  call T%create(numparticles=num_particles, numstates=num_states)

  call check(error, all(T%get_pos() == [(0.0D0, i=1,traj_dim)]), "Initial positions should be 0")
  if (allocated(error)) return
  call check(error, all(T%get_mom() == [(0.0D0, i=1,traj_dim)]), "Initial momentum should be 0")
  if (allocated(error)) return
  call check(error, all(T%get_mom2() == [(0.0D0, i=1,traj_dim)]), "Initial momentum should be 0")
  if (allocated(error)) return

  ipart = 2
  call T%set_pos(ipart, [1.0D0, 2.0D0, 3.0D0])
  call T%set_mom(ipart, [4.0D0, 5.0D0, 6.0D0])
  call T%set_mom2(ipart, [7.0D0, 8.0D0, 9.0D0])

  call check(error, all(T%get_pos(1) == [0.0D0, 0.0D0, 0.0D0]), "Expected position 0.0")
  if (allocated(error)) return
  call check(error, all(T%get_pos(2) == [1.0D0, 2.0D0, 3.0D0]), "Expected positions 1 2 3")
  if (allocated(error)) return

  call check(error, all(T%get_mom(1) == [0.0D0, 0.0D0, 0.0D0]), "Expected momentum 0.0")
  if (allocated(error)) return
  call check(error, all(T%get_mom(2) == [4.0D0, 5.0D0, 6.0D0]), "Expected momentum 4 5 6")
  if (allocated(error)) return

  call check(error, all(T%get_mom2(1) == [0.0D0, 0.0D0, 0.0D0]), "Expected momentum2 0.0")
  if (allocated(error)) return
  call check(error, all(T%get_mom2(2) == [7.0D0, 8.0D0, 9.0D0]), "Expected momentum2 7 8 9")
  if (allocated(error)) return
end subroutine test_setget_vec


subroutine test_setget_all(error)
  type(error_type), allocatable, intent(out) :: error
  type(t_trajectory) :: T
  integer, parameter :: num_particles=2, num_states=3, numdim = 3
  integer, parameter :: traj_dim = num_particles * numdim
  real(DefReal), dimension(traj_dim) :: pos, mom, mom2
  real(DefReal), parameter :: mass = 2.0D0
  integer :: i

  call T%create(numparticles=num_particles, numstates=num_states)

  ! First three elements are for particle 1, last three elements for particle 2
  pos = [1.0D0, 2.0D0, 3.0D0, 4.0D0, 5.0D0, 6.0D0]
  mom = [10.0D0, 20.0D0, 30.0D0, 40.0D0, 50.0D0, 60.0D0]
  mom2 = [100.0D0, 200.0D0, 300.0D0, 400.0D0, 500.0D0, 600.0D0]

  call T%set_pos(pos)
  call check(error, all(T%get_pos() == pos), "Unexpected position")

  if (allocated(error)) return
  call T%set_mom(mom)
  call check(error, all(T%get_mom() == mom), "Unexpected momentum")
  if (allocated(error)) return

  call T%set_mom2(mom2)
  call check(error, all(T%get_mom2() == mom2), "Unexpected momentum2")
  if (allocated(error)) return

  do i = 1, num_particles
    T%Particle(i)%mass = mass
  end do

  call T%set_vel(mom)
  call check(error, all(T%get_vel() == mom), "Unexpected velocity")
  if (allocated(error)) return

  call T%set_vel2(mom2)
  call check(error, all(T%get_vel2() == mom2), "Unexpected velocity2")
  if (allocated(error)) return

  ! Momentum should be twice the velocity now!
  call check(error, all(T%get_mom() == mom * mass), "Unexpected momentum after setting velocity")
  if (allocated(error)) return
  call check(error, all(T%get_mom2() == mom2 * mass), "Unexpected momentum2 after setting velocity2")
  if (allocated(error)) return

end subroutine test_setget_all

end module test_trajectory
