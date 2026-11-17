!> Test suite for RingModule: PA-CMD ring-polymer dynamics tests.
   !! Deterministic PA-CMD NVE trajectory test.
   !! Quartic potential V(x) = x^4/4, n=4 beads, beta=8, dt=0.001, 5 steps.
   !! Checks bead positions and momenta at every timestep.
module test_ring
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use RingModule, only: T_Ring, ring_create, ring_destroy, ring_setconsts, ringprop, &
                       ring_T, ring_V, ring_centvec, nmtran_forward, nmtran_backward, &
                       FMS_rings, FMS_ring_setup, FMS_ring_compact, get_centroidforce
   use TrajectoryModule, only: T_Trajectory
   use ToyModelModule, only: quartic_force
   private
   public :: collect_ring_suite

   real(8), parameter :: RING_MASS   = 1.0d0
   real(8), parameter :: RING_OMEGA2 = 1.0d0

contains

   subroutine collect_ring_suite(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)
      tests = [ &
         new_unittest("PA-CMD NVE quartic n=4 beta=8: 5-step trajectory", &
                      test_pacmd_nve_quartic), &
         new_unittest("PA-CMD quartic: two trajectories, explicit beads", &
                      test_pacmd_nve_quartic_2trajs_set_beads), &
         new_unittest("PA-CMD quartic: two trajectories, centroid beads", &
                      test_pacmd_nve_quartic_2trajs_centroid_beads), &
         new_unittest("centroid force: bead mean and trajectory lookup", &
                      test_get_centroidforce) &
      ]
   end subroutine collect_ring_suite

   subroutine test_get_centroidforce(error)
      type(error_type), allocatable, intent(out) :: error
      type(T_Trajectory) :: trajectory
      real(8), allocatable :: force(:)
      logical :: available, success

      call trajectory%create(numparticles=1, numstates=1)
      trajectory%TrajID = 42

      call FMS_ring_setup(2)
      success = ring_create(FMS_rings(1), 3, 1, 2, 0)
      call check(error, success, "ring_create failed for first centroid-force ring")
      if (allocated(error)) return
      success = ring_create(FMS_rings(2), 3, 1, 2, 0)
      call check(error, success, "ring_create failed for second centroid-force ring")
      if (allocated(error)) return

      FMS_rings(1)%traj_id = 7
      FMS_rings(2)%traj_id = 42
      FMS_rings(2)%time = 1.d0
      FMS_rings(2)%force_state = trajectory%StateID
      FMS_rings(2)%forces_current = .true.
      FMS_rings(2)%F_ext(:, 1, 1) = [1.d0, 2.d0, 3.d0]
      FMS_rings(2)%F_ext(:, 1, 2) = [3.d0, 4.d0, 5.d0]

      force = get_centroidforce(trajectory, available)
      call check(error, available, "centroid force was not available")
      if (allocated(error)) return
      call check(error, maxval(abs(force - [2.d0, 3.d0, 4.d0])) < 1.d-12, &
                 "centroid force is not the arithmetic bead mean")
      if (allocated(error)) return

      call FMS_ring_compact([.false., .true.])
      force = get_centroidforce(trajectory, available)
      call check(error, available, "trajectory lookup failed after ring compaction")
      if (allocated(error)) return
      call check(error, maxval(abs(force - [2.d0, 3.d0, 4.d0])) < 1.d-12, &
                 "centroid force changed after ring compaction")
      if (allocated(error)) return

      FMS_rings(1)%forces_current = .false.
      force = get_centroidforce(trajectory, available)
      call check(error, .not. available, "stale ring force was reported as available")
      if (allocated(error)) return

      call FMS_ring_setup(0)
      call trajectory%destroy()
   end subroutine test_get_centroidforce

   subroutine test_pacmd_nve_quartic(error)
      implicit none
   type(error_type), allocatable, intent(out) :: error

   type(T_Ring) :: ring_nve

   integer,  parameter :: nbead   = 4, n_steps = 5, nNHC = 0
   real(8),  parameter :: beta    = 8.0d0, dt = 0.001d0
   real(8),  parameter :: ttol    = 1.0d-8   ! trajectory-regression tolerance

   ! Per-bead trajectory at each step
   real(8) :: r_traj(nbead, n_steps), p_traj(nbead, n_steps)

   ! Reference bead positions at each step from previous run
   real(8), parameter :: ref_r(nbead, n_steps) = reshape([ &
      -5.45951147593796482d-1, -6.63405979951143809d-1, &
      -3.61077063929434394d-1, -6.06882999094416675d-1, &
      -5.46406534220790885d-1, -6.63725907574599261d-1, &
      -3.61613272064908287d-1, -6.07166660049868523d-1, &
      -5.46861713309533259d-1, -6.64045487221426645d-1, &
      -3.62149535035332537d-1, -6.07450095038639337d-1, &
      -5.47316684379841933d-1, -6.64364718515562780d-1, &
      -3.62685852558091937d-1, -6.07733303844699257d-1, &
      -5.47771446951185403d-1, -6.64683601081009545d-1, &
      -3.63222224349862621d-1, -6.08016286251884974d-1 &
      ], [nbead, n_steps])

   ! Reference bead momenta at each step from previous run
   real(8), parameter :: ref_p(nbead, n_steps) = reshape([ &
      -5.20761232724217371d-1, -2.59406359377475360d-1, &
      -5.84805268910893128d-1, -2.30573126787987775d-1, &
      -5.20642859580084605d-1, -2.59009276153129742d-1, &
      -5.84895062991512282d-1, -2.30272819345356622d-1, &
      -5.20524001832093575d-1, -2.58611858346560419d-1, &
      -5.84984529418392074d-1, -2.29972304417410844d-1, &
      -5.20404658947296261d-1, -2.58214106147113209d-1, &
      -5.85073667374723327d-1, -2.29671582036331701d-1, &
      -5.20284830393634046d-1, -2.57816019745762737d-1, &
      -5.85162476042685786d-1, -2.29370652235020900d-1 &
      ], [nbead, n_steps])


   real(8) :: phys_mass(1), tau0
   integer :: k, j
   logical :: ZSuccess

   phys_mass(1) = RING_MASS

   ! ── ring creation ──────────────────────────────────────────────────────────
   ZSuccess = ring_create(ring_nve, 1, 1, nbead, nNHC)
   call check(error, ZSuccess, "ring_create failed"); if (allocated(error)) return
   tau0 = 0.0d0
   call ring_setconsts(ring_nve, phys_mass, beta, 'PA-CMD  ', 'NONE       ', tau0, .false., gamma_sq_in=-1.0d0)
   ! ── deterministic initial conditions (bead space) ─────────────────────────
   ring_nve%r(1,1, 1) = -0.54549555390838034d0
   ring_nve%r(1,1, 2) = -0.66308570472718709d0
   ring_nve%r(1,1, 3) = -0.36054091091081814d0
   ring_nve%r(1,1, 4) = -0.60659911238817998d0

   ring_nve%p(1,1, 1) = -0.52087912179832518d0
   ring_nve%p(1,1, 2) = -0.25980310783187927d0
   ring_nve%p(1,1, 3) = -0.58471514799233082d0
   ring_nve%p(1,1, 4) = -0.23087322671384289d0

   ! ── propagate n_steps NVE steps ────────────────────────────────────────────
   do k = 1, n_steps
      call ringprop(ring_nve, dt, quartic_force, skip0=.false., skipt0=.true.)
      r_traj(:, k) = ring_nve%r(1, 1, :)
      p_traj(:, k) = ring_nve%p(1, 1, :)
   end do

   ZSuccess = ring_destroy(ring_nve)

   ! ── regression: all bead positions and momenta at every step ───────────────
   do k = 1, n_steps
      do j = 1, nbead
         call check(error, abs(r_traj(j,k) - ref_r(j,k)) < ttol, &
                    "r bead " // fmt_int(j) // " step " // fmt_int(k) // &
                    " mismatch: got " // fmt_real(r_traj(j,k)) // &
                    " exp " // fmt_real(ref_r(j,k)))
         if (allocated(error)) return
         call check(error, abs(p_traj(j,k) - ref_p(j,k)) < ttol, &
                    "p bead " // fmt_int(j) // " step " // fmt_int(k) // &
                    " mismatch: got " // fmt_real(p_traj(j,k)) // &
                    " exp " // fmt_real(ref_p(j,k)))
         if (allocated(error)) return
      end do
   end do

   end subroutine test_pacmd_nve_quartic

   subroutine test_pacmd_nve_quartic_2trajs_set_beads(error)
      implicit none
      type(error_type), allocatable, intent(out) :: error

      call test_pacmd_nve_quartic_2trajs(error, .false.)
   end subroutine test_pacmd_nve_quartic_2trajs_set_beads

   subroutine test_pacmd_nve_quartic_2trajs_centroid_beads(error)
      implicit none
      type(error_type), allocatable, intent(out) :: error

      call test_pacmd_nve_quartic_2trajs(error, .true.)
   end subroutine test_pacmd_nve_quartic_2trajs_centroid_beads

   subroutine test_pacmd_nve_quartic_2trajs(error, centroid_beads)
      implicit none
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: centroid_beads

      integer, parameter :: ntraj = 2, nbead = 4, n_steps = 5, nNHC = 2
      real(8), parameter :: beta = 8.0d0, dt = 0.001d0, ttol = 1.0d-12
      real(8), parameter :: initial_r(nbead) = [ &
         -0.54549555390838034d0, -0.66308570472718709d0, &
         -0.36054091091081814d0, -0.60659911238817998d0]
      real(8), parameter :: initial_p(nbead) = [ &
         -0.52087912179832518d0, -0.25980310783187927d0, &
         -0.58471514799233082d0, -0.23087322671384289d0]

      type(T_Ring) :: ring_nve(ntraj)
      real(8) :: phys_mass(1), tau0
      real(8) :: start_r(nbead, ntraj), start_p(nbead, ntraj)
      real(8) :: r_traj(nbead, n_steps, ntraj)
      real(8) :: p_traj(nbead, n_steps, ntraj)
      integer :: itraj, k, j
      logical :: ZSuccess

      phys_mass(1) = RING_MASS
      tau0 = 0.0d0

      do itraj = 1, ntraj
         ZSuccess = ring_create(ring_nve(itraj), 1, 1, nbead, nNHC)
         call check(error, ZSuccess, "ring_create failed for trajectory " // fmt_int(itraj))
         if (allocated(error)) return
         call ring_setconsts(ring_nve(itraj), phys_mass, beta, 'PA-CMD  ', &
                             'NHC-L     ', tau0, .false., gamma_sq_in=-1.0d0)
      end do

      if (centroid_beads) then
         start_r(:, 1) = sum(initial_r)/dble(nbead)
         start_p(:, 1) = sum(initial_p)/dble(nbead)
      else
         start_r(:, 1) = initial_r
         start_p(:, 1) = initial_p
      end if
      start_r(:, 2) = -start_r(:, 1)
      start_p(:, 2) = -start_p(:, 1)

      do itraj = 1, ntraj
         ring_nve(itraj)%r(1, 1, :) = start_r(:, itraj)
         ring_nve(itraj)%p(1, 1, :) = start_p(:, itraj)
      end do

      do k = 1, n_steps
         do itraj = 1, ntraj
            call ringprop(ring_nve(itraj), dt, quartic_force, &
                          skip0=.false., skipt0=.true.)
            r_traj(:, k, itraj) = ring_nve(itraj)%r(1, 1, :)
            p_traj(:, k, itraj) = ring_nve(itraj)%p(1, 1, :)
         end do
      end do

      do itraj = 1, ntraj
         call check(error, abs(ring_nve(itraj)%time - n_steps*dt) < ttol, &
                    "incorrect time for trajectory " // fmt_int(itraj))
         if (allocated(error)) return
         ZSuccess = ring_destroy(ring_nve(itraj))
         call check(error, ZSuccess, "ring_destroy failed for trajectory " // fmt_int(itraj))
         if (allocated(error)) return
      end do

      call check(error, maxval(abs(r_traj(:, n_steps, 1) - start_r(:, 1))) > ttol, &
                 "trajectory positions did not change")
      if (allocated(error)) return

      do k = 1, n_steps
         do j = 1, nbead
            call check(error, abs(r_traj(j, k, 1) + r_traj(j, k, 2)) < ttol, &
                       "position symmetry mismatch at bead " // fmt_int(j) // &
                       " step " // fmt_int(k))
            if (allocated(error)) return
            call check(error, abs(p_traj(j, k, 1) + p_traj(j, k, 2)) < ttol, &
                       "momentum symmetry mismatch at bead " // fmt_int(j) // &
                       " step " // fmt_int(k))
            if (allocated(error)) return
         end do
      end do

   end subroutine test_pacmd_nve_quartic_2trajs

   function fmt_int(n) result(s)
      integer,          intent(in) :: n
      character(len=16)            :: s
      write(s, '(I0)') n
   end function fmt_int

   function fmt_real(x) result(s)
      real(8),          intent(in) :: x
      character(len=32)            :: s
      write(s, '(ES23.15E3)') x
   end function fmt_real

end module test_ring

