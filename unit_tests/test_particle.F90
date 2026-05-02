module test_particle
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule, only: DefReal, DefInt
   use ParticleModule
   implicit none
   private

   public :: collect_particle_suite

contains

!> Collect all exported unit tests
   subroutine collect_particle_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("create", test_create), &
                  new_unittest("destroy", test_destroy), &
                  new_unittest("assign", test_assign), &
                  new_unittest("set_pos_vec", test_set_pos_vec), &
                  new_unittest("set_mom_vec", test_set_mom_vec), &
                  new_unittest("set_pos_comp", test_set_pos_comp), &
                  new_unittest("set_mom_comp", test_set_mom_comp), &
                  new_unittest("get_comp_wrong_index", test_get_comp_wrong_index), &
                  new_unittest("set_comp_wrong_index", test_set_comp_wrong_index), &
                  new_unittest("set_vec_wrong_size", test_set_vec_wrong_size), &
                  new_unittest("vec_norm", test_vec_norm), &
                  new_unittest("cross", test_cross), &
                  new_unittest("distance", test_distance_particle), &
                  new_unittest("angle", test_angle_particle), &
                  new_unittest("dihedral", test_dihedral_particle), &
                  new_unittest("distance_angle_wrong_dims", test_distance_angle_wrong_dims) &
                  ]

   end subroutine collect_particle_suite

   subroutine test_create(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: ID = 100, NUMDIM = 3

      call P%create(id=ID, numdim=NUMDIM)

      call check(error, P%ParticleID, ID)
      if (allocated(error)) return

      call check(error, P%NumDimensions == NUMDIM, "Invalid particle dimensions")
      if (allocated(error)) return

      ! The allocatable arrays should be allocated upon Particle creation
      ! and be initialized to 0.0
      call check(error, size(P%get_pos()) == NUMDIM, "Invalid particle position dimensions")
      if (allocated(error)) return
      call check(error, all(P%get_pos() == [0.0d0, 0.0d0, 0.0d0]), "Expected zero particle momenta")
      if (allocated(error)) return

      call check(error, size(P%get_mom()) == numdim, "Invalid particle momenta dimensions")
      if (allocated(error)) return
      call check(error, all(P%get_mom() == [0.0d0, 0.0d0, 0.0d0]), "Expected zero particle momenta")
      if (allocated(error)) return

   end subroutine test_create

   subroutine test_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: ID = 100, NUMDIM = 3

      call P%create(id=id, numdim=numdim)
      P%elmnt = "AB"
      P%mass = 10.0d0
      P%width = 10.0d0
      P%charge = 10.0d0
      P%atomicnum = 10.0d0

      call P%destroy()

      call check(error, P%particleID, 0)
      if (allocated(error)) return

      call check(error, P%NumDimensions, 0)
      if (allocated(error)) return

      call check(error, P%elmnt, "XX")
      if (allocated(error)) return

      call check(error, P%mass, 0.0d0)
      if (allocated(error)) return

      call check(error, P%width, 0.0d0)
      if (allocated(error)) return

      call check(error, P%charge, 0.0d0)
      if (allocated(error)) return

      call check(error, P%atomicNum, 0.0d0)
      if (allocated(error)) return

      ! TODO: Somehow check that the private arrays are deallocated
      !call check(error, size(P%get_mom()) == 1, "particle momentum array not deallocated")
      !if (allocated(error)) return
   end subroutine test_destroy

   subroutine test_assign(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P1, P2
      integer, parameter :: ID1 = 10, ID2 = 20
      integer, parameter :: NUMDIM1 = 2, NUMDIM2 = 1

      call P1%create(id=ID1, numdim=NUMDIM1)
      call P2%create(id=ID2, numdim=NUMDIM2)

      call P1%set_pos([1.0d0, 1.0d0])
      call P1%set_mom([2.0d0, 2.0d0])
      call P1%set_mom2([3.0d0, 3.0d0])
      P1%Elmnt = "YY"
      P1%mass = 10.0d0
      P1%Width = 11.0d0
      P1%AtomicNum = 12.0d0
      P1%Charge = -1.0d0

      P2 = P1

      call check(error, P1%ParticleID, ID1)
      if (allocated(error)) return
      call check(error, P2%ParticleID, ID1)
      if (allocated(error)) return

      call check(error, P1%NumDimensions, NUMDIM1)
      if (allocated(error)) return
      call check(error, P2%NumDimensions, NUMDIM1)
      if (allocated(error)) return

      call check(error, P2%Elmnt, "YY")
      if (allocated(error)) return
      call check(error, P2%mass, 10.d0)
      if (allocated(error)) return
      call check(error, P2%width, 11.d0)
      if (allocated(error)) return
      call check(error, P2%atomicnum, 12.d0)
      if (allocated(error)) return
      call check(error, P2%charge, -1.0d0)
      if (allocated(error)) return

      call check(error, all(P1%get_pos() == P2%get_pos()))
      if (allocated(error)) return
      call check(error, all(P1%get_mom() == P2%get_mom()))
      if (allocated(error)) return
      call check(error, all(P1%get_mom2() == P2%get_mom2()))
      if (allocated(error)) return
   end subroutine test_assign

   subroutine test_set_pos_vec(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: id = 1, numdim = 2
      real(DefReal), allocatable :: positions(:)

      call P%create(id=id, numdim=numdim)

      call P%set_pos([0.0d0, 1.0d0])
      positions = P%get_pos()

      call check(error, size(positions) == 2, "Invalid particle position dimensions")
      if (allocated(error)) return
      call check(error, all(positions == [0.0d0, 1.0d0]), "Invalid particle positions")
      if (allocated(error)) return

   end subroutine test_set_pos_vec

   subroutine test_set_mom_vec(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: id = 1, numdim = 2

      call P%create(id=id, numdim=numdim)

      ! The allocatable arrays should be allocated upon Particle creation
      ! and be initialized to 0.0
      call check(error, size(P%get_mom()) == 2, "Invalid particle momenta dimensions")
      if (allocated(error)) return
      call check(error, size(P%get_mom2()) == 2, "Invalid particle momenta2 dimensions")
      if (allocated(error)) return

      call check(error, all(P%get_mom() == [0.0d0, 0.0d0]), "Expected zero particle momenta")
      if (allocated(error)) return
      call check(error, all(P%get_mom2() == [0.0d0, 0.0d0]), "Expected zero particle momenta2")
      if (allocated(error)) return

      call P%set_mom([-2.0d0, 2.0d0])
      call P%set_mom2([-1.0d0, 1.0d0])

      call check(error, all(P%get_mom() == [-2.0d0, 2.0d0]), "Invalid particle momenta")
      call check(error, all(P%get_mom2() == [-1.0d0, 1.0d0]), "Invalid particle momenta")
      if (allocated(error)) return
   end subroutine test_set_mom_vec

   subroutine test_set_pos_comp(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: id = 1, numdim = 2

      call P%create(id=id, numdim=numdim)

      call P%set_pos(2, -1.0d0)

      call check(error, size(P%get_pos()) == numdim, "Invalid particle position dimensions")
      if (allocated(error)) return

      call check(error, all(P%get_pos() == [0.0d0, -1.0d0]), "Invalid particle positions")
      if (allocated(error)) return
   end subroutine test_set_pos_comp

   subroutine test_set_mom_comp(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: id = 1, numdim = 2

      call P%create(id=id, numdim=numdim)

      call P%set_mom(2, -2.0d0)

      call check(error, size(P%get_mom()) == numdim, "Invalid particle momentum dimensions")
      if (allocated(error)) return

      call check(error, all(P%get_mom() == [0.0d0, -2.0d0]), "Invalid particle momenta")
      if (allocated(error)) return
   end subroutine test_set_mom_comp

   subroutine test_set_vec_wrong_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: ID = 1, NUMDIM = 1

      call P%create(id=id, numdim=NUMDIM)

      ! These should fail, passing an array of wrong size!
      call P%set_pos([0.0d0, 0.0d0])
      call check_dieerror_called(error, "Particle%set_position_vec: wrong array size")
      if (allocated(error)) return

      call P%set_mom([0.0d0, 0.0d0])
      call check_dieerror_called(error, "Particle%set_momentum_vec: wrong array size")

      call P%set_mom2([0.0d0, 0.0d0])
      call check_dieerror_called(error, "Particle%set_momentum2_vec: wrong array size")
   end subroutine test_set_vec_wrong_size

   subroutine test_get_comp_wrong_index(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: ID = 1, NUMDIM = 3
      real(DefReal) :: pos

      call P%create(id=id, numdim=NUMDIM)

      pos = P%get_pos(0)
      call check_dieerror_called(error, "Particle%get_position_component: index out of range")
      if (allocated(error)) return

      pos = P%get_pos(NUMDIM + 1)
      call check_dieerror_called(error, "Particle%get_position_component: index out of range")
      if (allocated(error)) return

      pos = P%get_mom(0)
      call check_dieerror_called(error, "Particle%get_momentum_component: index out of range")
      if (allocated(error)) return

      pos = P%get_mom(NUMDIM + 1)
      call check_dieerror_called(error, "Particle%get_momentum_component: index out of range")
      if (allocated(error)) return
   end subroutine test_get_comp_wrong_index

   subroutine test_set_comp_wrong_index(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P
      integer, parameter :: ID = 1, NUMDIM = 2

      call P%create(id=id, numdim=NUMDIM)

      call P%set_pos(0, 1.0d0)
      call check_dieerror_called(error, "Particle%set_position_component: index out of range")
      if (allocated(error)) return

      call P%set_pos(NUMDIM + 1, 1.0d0)
      call check_dieerror_called(error, "Particle%set_position_component: index out of range")
      if (allocated(error)) return

      call P%set_mom(0, 1.0d0)
      call check_dieerror_called(error, "Particle%set_momentum_component: index out of range")
      if (allocated(error)) return

      call P%set_mom(NUMDIM + 1, 1.0d0)
      call check_dieerror_called(error, "Particle%set_momentum_component: index out of range")
      if (allocated(error)) return
   end subroutine test_set_comp_wrong_index

   subroutine test_vec_norm(error)
      type(error_type), allocatable, intent(out) :: error
      real(DefReal), parameter :: test_vec(3) = [1.0, 2.0, 2.0]
      real(DefReal) :: norm

      norm = vector_norm(test_vec)
      call check(error, norm == 3.0, "Value of vector norm is incorrect")
      if (allocated(error)) return
   end subroutine test_vec_norm

   subroutine test_cross(error)
      type(error_type), allocatable, intent(out) :: error
      real(DefReal), parameter :: test_vec1(3) = [1.0, 0.0, 0.0]
      real(DefReal), parameter :: test_vec2(3) = [0.0, 1.0, 0.0]
      integer, parameter :: NUMDIM = 3
      real(DefReal) :: cross_vec(3)

      cross_vec = cross(test_vec1, test_vec2)
      call check(error, size(cross_vec) == NUMDIM, "Invalid cross product vector dimension")
      if (allocated(error)) return

      call check(error, all(cross_vec == [0.0d0, 0.0d0, 1.0d0]), "Invalid cross product")
      if (allocated(error)) return
   end subroutine test_cross

   subroutine test_distance_particle(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P1, P2
      real(DefReal) :: distance

      call P1%create(id=1, numdim=3)
      call P2%create(id=2, numdim=3)

      call P1%set_pos([1.0d0, 3.0d0, 3.0d0])
      call P2%set_pos([0.0d0, 1.0d0, 1.0d0])

      distance = FMS_Distance(P1, P2)
      call check(error, distance == 3.0, "Distance between particles is incorrect")
      if (allocated(error)) return
   end subroutine test_distance_particle

   subroutine test_angle_particle(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P1, P2, P3
      real(DefReal) :: angle

      call P1%create(id=1, numdim=3)
      call P2%create(id=2, numdim=3)
      call P3%create(id=3, numdim=3)

      call P1%set_pos([1.0d0, 0.0d0, 0.0d0])
      call P2%set_pos([0.0d0, 0.0d0, 0.0d0])
      call P3%set_pos([-1.0d0, 0.0d0, 0.0d0])

      angle = FMS_Angle(P1, P2, P3)
      call check(error, angle == acos(-1.0d0), "Angle between particles is incorrect")
      if (allocated(error)) return
   end subroutine test_angle_particle

   subroutine test_dihedral_particle(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P1, P2, P3, P4
      real(DefReal) :: dihedral

      call P1%create(id=1, numdim=3)
      call P2%create(id=2, numdim=3)
      call P3%create(id=3, numdim=3)
      call P4%create(id=3, numdim=3)

      call P1%set_pos([-1.0d0, 0.0d0, 0.0d0])
      call P2%set_pos([0.0d0, 0.0d0, 0.0d0])
      call P3%set_pos([0.0d0, 1.0d0, 0.0d0])
      call P4%set_pos([1.0d0, 1.0d0, 0.0d0])

      dihedral = FMS_Dihedral(P1, P2, P3, P4)
      call check(error, dihedral, acos(-1.0d0))
      if (allocated(error)) return
   end subroutine test_dihedral_particle

   subroutine test_distance_angle_wrong_dims(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_Particle) :: P1, P2, P3
      real(DefReal) :: result

      call P1%create(id=1, numdim=1)
      call P2%create(id=2, numdim=2)
      call P3%create(id=3, numdim=3)

      result = FMS_Distance(P1, P2)
      call check_dieerror_called(error, "Distance: Particles must have same dim!")
      if (allocated(error)) return

      result = FMS_Angle(P1, P2, P3)
      call check_dieerror_called(error, "Angle: Particles must have same dim!")
      if (allocated(error)) return

      result = FMS_Angle(P1, P1, P3)
      call check_dieerror_called(error, "Angle: Particles must have same dim!")
      if (allocated(error)) return

   end subroutine test_distance_angle_wrong_dims

end module test_particle
