module test_overlap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule
   use OverlapModule, only: overlap, overlap_S_dot, overlap_V
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
                  new_unittest('overlap_S_dot', test_overlap_S_dot), &
                  new_unittest('overlap_V_case_A', test_overlap_V_case_A), &
                  new_unittest('overlap_V_case_B_no_centroids', test_overlap_V_case_B_no_centroids), &
                  new_unittest('overlap_V_case_B_centroids', test_overlap_V_case_B_centroids), &
                  new_unittest('overlap_V_case_B_small_overlap_no_centroids', test_overlap_V_case_B_small_overlap_nocent), &
                  new_unittest('overlap_V_case_B_small_overlap_centroids', test_overlap_V_case_B_small_overlap_cent), &
                  new_unittest('overlap_V_centroid_mismatch', test_overlap_V_centroid_mismatch), &
                  new_unittest('overlap_V_centroid_missing', test_overlap_V_centroid_missing) &
                  ]

   end subroutine collect_overlap_suite

   function get_numerical_time_derivative(T1, T2, integral_type, nr_steps, step_size) result(num_dt)
      use VerletModule, only: FMS_PropVV
      type(T_Trajectory), intent(in) :: T1, T2
      character, intent(in) :: integral_type
      integer(kind=DefInt), intent(in) :: nr_steps
      real(kind=DefReal), intent(in) :: step_size
      real(kind=DefReal) :: step
      type(T_Trajectory) :: T1_tmp, T2_tmp
      complex(kind=DefReal) :: num_dt(nr_steps)
      complex(kind=DefReal) :: tmp_num_deriv, tmp
      integer(kind=DefInt) :: j_step, k

      num_dt = 0.d0
      step = step_size
      do j_step = 1, nr_steps
         tmp_num_deriv = 0.d0

         do k = 1, 2
            T1_tmp = T1
            T2_tmp = T2

            !call FMS_PropVV(T1_tmp, (-1.d0)**k * step_size)
            call FMS_PropVV(T2_tmp, (-1.d0)**k * step)

            select case (integral_type)
            case ('S')
               tmp = overlap(T1_tmp, T2_tmp)
            case default
               print*,'You really should not be seeing this error message!'
               stop 1
            end select
            tmp_num_deriv = tmp_num_deriv + (-1)**k * tmp
         end do

         num_dt = 0.5d0 * tmp_num_deriv / step
         step = step * step_size
      end do
   end function get_numerical_time_derivative

   subroutine test_overlap_S_dot(error)
      use GlobalModule, only: gliModel, gliMethod
      use TrajectoryCalcsModule, only: FMS_Forces
      use ToyModelModule, only: izmaylov_params
      use ElecStrucModule
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: T1, T2
      integer(kind=DefInt) :: num_particles, num_states
      complex(kind=DefReal) :: SDot_analytical(7)
      complex(kind=DefReal) :: SDot_numerical(7)
      real(kind=DefReal) :: pos(6), mom(6), vel(6)
      real(kind=DefReal) :: mass
      real(kind=DefReal) :: min_err

      num_particles = 2
      num_states = 1
      gliModel = 0
      gliMethod = 3
      fmiOut = open_dev_null()
      call izmaylov_params%initialize(W1=0.009557d0, W2=0.003515d0, XA=20.07d0, YA=0.0d0, &
                                      deltaE=0.01984d0, coupC=0.0006127d0)
      call FMS_ESInit(num_particles, num_states)
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
      SDot_numerical = get_numerical_time_derivative(T1, T2, 'S', 7, 1.d-1)
      ! since the numerical derivative suffers from catastrophic cancelations
      ! at small h, we look for the smallest error
      ! see, e.g., https://en.wikipedia.org/wiki/Numerical_differentiation#/media/File%3AAbsoluteErrorNumericalDifferentiationExample.png
      min_err = minval(abs(SDot_analytical - SDot_numerical))
      call check(error, min_err < atol + rtol * abs(SDot_analytical(1)), 'Analytical SDot ' &
                 //'significantly from numerical one')
      if (allocated(error)) return

      call T1%destroy()
      call T2%destroy()
   end subroutine test_overlap_S_dot

   !> Case A: energy of a single trajectory, <X_i^I|V|X_i^I> = V(X_i)
   subroutine test_overlap_V_case_A(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: T
      complex(kind=DefReal) :: V
      real(kind=DefReal), parameter :: pot(2) = [0.123d0, 0.456d0]

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.

      call T%create(numparticles=2, numstates=2)
      T%StateID = 2
      ! Mark the potential as up to date, so no electronic structure call is made
      T%ElecStruc%PotEn = pot
      T%ElecStruc%ModPot = 0.d0
      T%ESFlags%ZPotEnCurrent = .true.

      V = overlap_V(T)

      call check(error, real(V), pot(2), thr=atol, message='V should equal potential energy of the trajectory state')
      if (allocated(error)) return
      call check(error, aimag(V), 0.d0, thr=atol, message='V should be purely real in case A')
      if (allocated(error)) return

      ! Modified potential (e.g. SMD) must be included
      T%ElecStruc%ModPot = 0.01d0
      V = overlap_V(T)
      call check(error, real(V), pot(2) + 0.01d0, thr=atol, message='V should include ModPot')
      if (allocated(error)) return

      call T%destroy()
   end subroutine test_overlap_V_case_A

   !> Creates a singlet trajectory with a preset (up-to-date) potential energy,
   !> so that no electronic structure call is needed.
   subroutine create_singlet(T, traj_id, state, x_shift, pot)
      type(t_trajectory), intent(out) :: T
      integer(kind=DefInt), intent(in) :: traj_id, state
      real(kind=DefReal), intent(in) :: x_shift
      real(kind=DefReal), intent(in) :: pot(:)

      call T%create(numparticles=2_defint, numstates=size(pot, kind=DefInt))
      T%TrajID = traj_id
      T%CBF = traj_id
      T%StateID = state
      T%triplet = .false.
      T%Ms = 0
      T%Phase = 0.1d0 * traj_id
      T%Particle(:)%mass = 1000.d0
      T%Particle(:)%width = 4.7d0
      call T%set_pos([x_shift, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0])
      call T%set_mom([0.1d0 * traj_id, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0])
      T%ElecStruc%PotEn = pot
      T%ElecStruc%ModPot = 0.d0
      T%ESFlags%ZPotEnCurrent = .true.
   end subroutine create_singlet

   !> Case B without centroids: V = S_ij * 0.5 * (V(X_i) + V(X_j))
   subroutine test_overlap_V_case_B_no_centroids(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj
      complex(kind=DefReal) :: V, S, S_pre, expected

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .false.
      gldOLapThresh = 1.d-8

      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])

      S = overlap(Ti, Tj)
      call check(error, abs(S) > 1.d-3, 'Test setup: overlap should be non-negligible')
      if (allocated(error)) return

      V = overlap_V(Ti, Tj)
      expected = S * 0.5d0 * (0.10d0 + 0.20d0)
      call check(error, abs(V - expected) < atol, 'V should be S * average potential')
      if (allocated(error)) return

      ! precalculated overlap must be used instead of recomputing
      S_pre = (0.3d0, -0.2d0)
      V = overlap_V(Ti, Tj, S_ij_precalc=S_pre)
      expected = S_pre * 0.5d0 * (0.10d0 + 0.20d0)
      call check(error, abs(V - expected) < atol, 'V should use precalculated overlap')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
   end subroutine test_overlap_V_case_B_no_centroids

   !> Case B with centroids: V = S_ij * V(X_c)
   subroutine test_overlap_V_case_B_centroids(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj, Tc
      complex(kind=DefReal) :: V, S, expected

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .true.
      gldOLapThresh = 1.d-8

      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])
      call create_singlet(Tc, 3_defint, 1_defint, 0.025d0, [0.33d0, 0.90d0])
      Tc%CentID = [1, 2]

      S = overlap(Ti, Tj)
      V = overlap_V(Ti, Tj, Tc)
      expected = S * 0.33d0
      call check(error, abs(V - expected) < atol, 'V should be S * potential at the centroid')
      if (allocated(error)) return

      ! centroid ordering should not matter
      Tc%CentID = [2, 1]
      V = overlap_V(Ti, Tj, Tc)
      call check(error, abs(V - expected) < atol, 'V should not depend on centroid ID ordering')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
      call Tc%destroy()
   end subroutine test_overlap_V_case_B_centroids

   !> Case B, negligible overlap, no centroids: V = 0
   subroutine test_overlap_V_case_B_small_overlap_nocent(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj
      complex(kind=DefReal) :: V

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .false.
      ! |S| <= 1, so every overlap is below this threshold
      gldOLapThresh = 2.d0

      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])

      V = overlap_V(Ti, Tj)
      call check(error, abs(V) < atol, 'V should vanish for negligible overlap')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
   end subroutine test_overlap_V_case_B_small_overlap_nocent

   !> Case B, negligible overlap, with centroids: V = 0 and centroid errors are ignored
   subroutine test_overlap_V_case_B_small_overlap_cent(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj, Tc
      complex(kind=DefReal) :: V

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .true.
      gldOLapThresh = 2.d0

      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])
      call create_singlet(Tc, 3_defint, 1_defint, 0.025d0, [0.33d0, 0.90d0])
      Tc%CentID = [1, 2]
      Tc%ESFlags%zIgnoreErrors = .false.

      V = overlap_V(Ti, Tj, Tc)
      call check(error, abs(V) < atol, 'V should vanish for negligible overlap')
      if (allocated(error)) return
      call check(error, Tc%ESFlags%zIgnoreErrors, 'Centroid should be flagged to ignore errors')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
      call Tc%destroy()
   end subroutine test_overlap_V_case_B_small_overlap_cent

   !> With centroids enabled, passing a centroid for different trajectories is an error
   subroutine test_overlap_V_centroid_mismatch(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj, Tc
      complex(kind=DefReal) :: V

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .true.
      gldOLapThresh = 1.d-8

      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])
      call create_singlet(Tc, 3_defint, 1_defint, 0.025d0, [0.33d0, 0.90d0])
      Tc%CentID = [1, 3]

      V = overlap_V(Ti, Tj, Tc)
      call check_dieerror_called(error, 'overlap_V_trajectory :: centroids dont match')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
      call Tc%destroy()
   end subroutine test_overlap_V_centroid_mismatch

   !> With centroids enabled, omitting the centroid is an error
   subroutine test_overlap_V_centroid_missing(error)
      type(error_type), allocatable, intent(out) :: error
      type(t_trajectory) :: Ti, Tj
      complex(kind=DefReal) :: V

      glzSPA1_SOC_model = .false.
      glzxfaims = .false.
      glzCentroids = .true.
      gldOLapThresh = 1.d-8

      ! Same CBF, so that the (absent) centroid is not dereferenced before the check
      call create_singlet(Ti, 1_defint, 1_defint, 0.d0, [0.10d0, 0.50d0])
      call create_singlet(Tj, 2_defint, 1_defint, 0.05d0, [0.20d0, 0.70d0])
      Tj%CBF = Ti%CBF

      V = overlap_V(Ti, Tj)
      call check_dieerror_called(error, 'overlap_V_trajectory :: no centroid passed')
      if (allocated(error)) return

      call Ti%destroy()
      call Tj%destroy()
   end subroutine test_overlap_V_centroid_missing

   function open_dev_null() result(output_unit)
      integer :: output_unit

      open (newunit=output_unit, file='/dev/null', action='write')
   end function open_dev_null

end module test_overlap
