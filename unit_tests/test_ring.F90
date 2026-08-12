!> Test suite for RingModule: PA-CMD ring-polymer dynamics tests.
   !! Deterministic PA-CMD NVE trajectory test.
   !! Quartic potential V(x) = x^4/4, n=16 beads, beta=8, dt=0.001, 5 steps.
   !! Checks: beads positions and momenta at each teimestep.
   !! The results have been benchmarked against I-PI 3.3.0 with the same initial conditions and parameters to 10^(-12)
module test_ring
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use RingModule, only: T_Ring, ring_create, ring_destroy, ring_setconsts, ringprop, &
                       ring_T, ring_V, ring_centvec, nmtran_forward, nmtran_backward
   private
   public :: collect_ring_suite

   real(8), parameter :: RING_MASS   = 1.0d0
   real(8), parameter :: RING_OMEGA2 = 1.0d0

contains

   subroutine collect_ring_suite(tests)
      type(unittest_type), allocatable, intent(out) :: tests(:)
      tests = [ &
         new_unittest("PA-CMD NVE quartic n=16 beta=8: 5-step trajectory", &
                      test_pacmd_nve_quartic) &
      ]
   end subroutine collect_ring_suite

   function quartic_force(r_vec, ndim, natom, n)
   !! V(x) = x^4/4,  F = -dV/dx = -x^3
      implicit none
      integer,  intent(in)                          :: ndim, natom, n
      real(8),  dimension(ndim,natom,n), intent(in) :: r_vec
      real(8),  dimension(ndim,natom,n)             :: quartic_force
      quartic_force = -r_vec**3
   end function quartic_force

   subroutine test_pacmd_nve_quartic(error)
      implicit none
   type(error_type), allocatable, intent(out) :: error

   type(T_Ring) :: ring_nve

   integer,  parameter :: nbead   = 16, n_steps = 5, nNHC = 0
   real(8),  parameter :: beta    = 8.0d0, dt = 0.001d0
   real(8),  parameter :: ttol    = 1.0d-8   ! trajectory-regression tolerance

   ! Per-bead trajectory at each step
   real(8) :: r_traj(nbead, n_steps), p_traj(nbead, n_steps)

   ! Reference bead positions at each step from previous run
   real(8), parameter :: ref_r(nbead, n_steps) = reshape([ &
      -0.54578260206119955d0, -0.66341302139117286d0, -0.36102983493832730d0, &
      -0.60689776677756480d0,  0.30107636441234642d0, -4.8998488979279220d-2, &
       7.8070960388122557d-2,  0.15087701483857641d0,  0.35601138107006619d0, &
       0.19124829318225023d0, -0.10705589418390984d0, -0.34541549742985256d0, &
      -0.33337905066873164d0, -0.80311495079547557d0, -0.96111423592642742d0,  0.36063951253804905d0, &
      -0.54606619352729868d0, -0.66373663328130239d0, -0.36151747198971329d0, &
      -0.60719425804984561d0,  0.30079892099773764d0, -4.9409432858347702d-2, &
       7.7709143900017597d-2,  0.15046407393268210d0,  0.35567425416382159d0, &
       0.19071971930743337d0, -0.10753371687236654d0, -0.34602679330632746d0, &
      -0.33367329464581041d0, -0.80340016077561449d0, -0.96121326869515455d0,  0.36064686323568462d0, &
      -0.54634632741313782d0, -0.66405653928139419d0, -0.36200382038669832d0, &
      -0.60748858605682332d0,  0.30051770120562954d0, -4.9822581129007226d-2, &
       7.7344119968138059d-2,  0.15004742634264198d0,  0.35533234788727930d0, &
       0.19018768982340506d0, -0.10801270127419735d0, -0.34663711633964134d0, &
      -0.33396575493741204d0, -0.80367985077818349d0, -0.96130541295817495d0,  0.36065280813565093d0, &
      -0.54662300287363330d0, -0.66437273831348942d0, -0.36248887846026689d0, &
      -0.60778075065286541d0,  0.30023270445970962d0, -5.0237934123181560d-2, &
       7.6975887695272860d-2,  0.14962707141095852d0,  0.35498566128735631d0, &
       0.18965220506428165d0, -0.10849284702251522d0, -0.34724646484980293d0, &
      -0.33425643123166771d0, -0.80395401999765581d0, -0.96139066897509662d0,  0.36065734629847945d0, &
      -0.54689621911219810d0, -0.66468522933789143d0, -0.36297264455069167d0, &
      -0.60807075169490943d0,  0.29994393022706589d0, -5.0655492131390090d-2, &
       7.6604446235227663d-2,  0.14920300853464585d0,  0.35463419346896019d0, &
       0.18911326540630904d0, -0.10897415373305189d0, -0.34785483716675059d0, &
      -0.33454532325008501d0, -0.80422266770486217d0, -0.96146903709940712d0,  0.36066047675360546d0 &
      ], [nbead, n_steps])

   ! Reference bead momenta at each step from previous run
   real(8), parameter :: ref_p(nbead, n_steps) = reshape([ &
      -0.51756201577266225d0, -0.25783104665962808d0, -0.58686189519348220d0, &
      -0.22603415503193588d0, -0.23459289156822416d0, -0.46966934632928575d0, &
      -0.27314851938853013d0, -0.43049892498217912d0, -0.15709876210696089d0, &
      -0.50960494036750792d0, -0.21535634588011474d0, -0.65434189316605851d0, &
      -0.12875480616976751d0, -0.46580782234381546d0, -0.29079634626213546d0, -6.3512209678668330d-2, &
      -0.51424364186212823d0, -0.25585904645023366d0, -0.58900703736745186d0, &
      -0.22119542930270505d0, -0.23965264900397518d0, -0.46776028381113416d0, &
      -0.27336624613233051d0, -0.42997277808639467d0, -0.15862394547013656d0, &
      -0.51014558678733235d0, -0.21511570066190372d0, -0.65329808566873082d0, &
      -0.13064541114743067d0, -0.46404224093698720d0, -0.28398953340837618d0, -7.2472618600780137d-2, &
      -0.51092402135142001d0, -0.25388712234430594d0, -0.59115056133147292d0, &
      -0.21635707885595337d0, -0.24471292578921203d0, -0.46585049291330827d0, &
      -0.27358436629239014d0, -0.42944609641208498d0, -0.16015006423569750d0, &
      -0.51068520460710076d0, -0.21487577606387714d0, -0.65225225931303232d0, &
      -0.13253714436329117d0, -0.46227540877264317d0, -0.27718279410575719d0, -8.1434593104754691d-2, &
      -0.50760317554171863d0, -0.25191528949655517d0, -0.59329245391523855d0, &
      -0.21151913302534853d0, -0.24977369099396801d0, -0.46393998387048579d0, &
      -0.27380287785598734d0, -0.42891888215953150d0, -0.16167710725178686d0, &
      -0.51122378989622108d0, -0.21463657253199994d0, -0.65120441897551606d0, &
      -0.13442999453685242d0, -0.46050734564088314d0, -0.27037618605072744d0, -9.0398080006258638d-2, &
      -0.50428112575043149d0, -0.24994356307552223d0, -0.59543270196099496d0, &
      -0.20668162114821984d0, -0.25483491369000544d0, -0.46202876691931244d0, &
      -0.27402177880851919d0, -0.42839113753486491d0, -0.16320506337018739d0, &
      -0.51176133873457641d0, -0.21439809050538985d0, -0.65015456954640249d0, &
      -0.13632395038328707d0, -0.45873807136146699d0, -0.26356976694775264d0, -9.9363026110840871d-2 &
      ], [nbead, n_steps])


   real(8) :: phys_mass(1), tau0
   integer :: k, j
   logical :: ZSuccess

   phys_mass(1) = RING_MASS

   ! ── ring creation ──────────────────────────────────────────────────────────
   ZSuccess = ring_create(ring_nve, 1, 1, nbead, nNHC)
   call check(error, ZSuccess, "ring_create failed"); if (allocated(error)) return
   tau0 = 0.0d0
   call ring_setconsts(ring_nve, phys_mass, beta, 'PA-CMD  ', 'NONE        ', tau0, .false., gamma_sq_in=-1.0d0)
   ! ── deterministic initial conditions (bead space) ─────────────────────────
   ring_nve%r(1,1, 1) = -0.54549555390838034d0
   ring_nve%r(1,1, 2) = -0.66308570472718709d0
   ring_nve%r(1,1, 3) = -0.36054091091081814d0
   ring_nve%r(1,1, 4) = -0.60659911238817998d0
   ring_nve%r(1,1, 5) =  0.30135003206915273d0
   ring_nve%r(1,1, 6) = -4.8589749118503056d-2
   ring_nve%r(1,1, 7) =  7.8429570380646352d-2
   ring_nve%r(1,1, 8) =  0.15128624977230398d0
   ring_nve%r(1,1, 9) =  0.35634372961706251d0
   ring_nve%r(1,1,10) =  0.19177341115585303d0
   ring_nve%r(1,1,11) = -0.10657923355834381d0
   ring_nve%r(1,1,12) = -0.34480323040012650d0
   ring_nve%r(1,1,13) = -0.33308302335139989d0
   ring_nve%r(1,1,14) = -0.80282422171960288d0
   ring_nve%r(1,1,15) = -0.96100831448622448d0
   ring_nve%r(1,1,16) =  0.36063075695114488d0

   ring_nve%p(1,1, 1) = -0.52087912179832518d0
   ring_nve%p(1,1, 2) = -0.25980310783187927d0
   ring_nve%p(1,1, 3) = -0.58471514799233082d0
   ring_nve%p(1,1, 4) = -0.23087322671384289d0
   ring_nve%p(1,1, 5) = -0.22953368441338173d0
   ring_nve%p(1,1, 6) = -0.47157767023511671d0
   ring_nve%p(1,1, 7) = -0.27293118807180394d0
   ring_nve%p(1,1, 8) = -0.43102453490498716d0
   ring_nve%p(1,1, 9) = -0.15557452530143814d0
   ring_nve%p(1,1,10) = -0.50906326928864942d0
   ring_nve%p(1,1,11) = -0.21559771126569846d0
   ring_nve%p(1,1,12) = -0.65538367694213440d0
   ring_nve%p(1,1,13) = -0.12686534070641242d0
   ring_nve%p(1,1,14) = -0.46757213323335917d0
   ring_nve%p(1,1,15) = -0.29760317498000299d0
   ring_nve%p(1,1,16) = -5.4553419512313378d-2

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

