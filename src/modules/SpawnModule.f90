!     Copyright Todd J. Martinez and Raphael D. Levine, 1994

!!    @brief Spawning parameters are stored in this module
!!    \see ReadNameList for description of parameters
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
module SpawnModule
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   use GlobalModule
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_KineticClass, FMS_PotentialT, &
                                    FMS_Coupling, FMS_SOCoupling, FMS_ClassEnergy
   use TrajectoryIOModule, only: FMS_WriteFXYZ
   use BundleModule
   use BundleCalcsModule, only: FMS_UpdateMulliken
   use OverlapModule, only: overlap
   use VerletModule, only: FMS_PropVV
   use XFAIMSModule, only: xfaims_params
   implicit none

   !> OMax is a legacy threshold value that functioned
   !> as OMax_intra, OMax_inter and Omin_parent.
   !> To preserve the semantics of existing input files,
   !> OMax_intra, OMax_inter and OMin_parent are assigned the OMax value
   !> unless the user provides their value in Control.dat
   real(kind=DefReal), parameter :: OMAX_DEFAULT = 0.8

   !> Parameters controlling the spawning algorithm
   type :: t_spawn_params
      !> threshold for entry to  spawn regions
      real(kind=DefReal) :: CSThresh
      !> threshold for exit from spawn regions
      real(kind=DefReal) :: CFThresh
      !> min pop to bother with a spawn
      real(kind=DefReal) :: PopToSpawn
      !> OMax_inter is a maximum overlap between parent trajectory
      !> and any existing trajectory on the target state.
      !> If this overlap is bigger than OMax_intra, we don't enter
      !> the spawning procedure.
      !> NOTE: This parameter is critical to prevent spawning
      !> right after a spawning event.
      real(kind=DefReal) :: OMax_inter
      !> OMax_intra is a maximum overlap that we allow between the child TBF
      !> and any other pre-existing TBFs on the same electronic state.
      !> If the overlap is larger than this threshold, we abort the spawning procedure
      !> so that we don't produce essentially equivalent TBFs.
      !> NOTE: The overlap check happens after child backpropagation.
      real(kind=DefReal) :: OMax_intra
      !> OMin_parent is a minimum overlap between child and parent
      !> If the overlap between the parent and the child TBFs at t_spawn
      !> is below this threshold, we abort the spawn.
      !> NOTE: The overlap between parent and child at spawning time
      !> in the coordinate space is 1 by definition,
      !> but can be less than 1 in momentum space due to
      !> velocity rescaling along the NAC vector.
      real(kind=DefReal) :: OMin_parent
      !> Maximum number of trajectories
      integer(kind=DefInt) :: MaxTraj
      !> Flag determining if continuous spawning (0) or single spawning (1) is employed
      integer(kind=DefInt) :: MultiSpawn
      !> Flag deteriming wether to return either coupling .dot. velocity or norm of coupling
      logical :: SpawnCoupV
      !> threshhold to SOC Norm for spawning
      real(kind=DefReal) :: SOCThresh
   contains
      procedure, public :: initialize => initialize_spawn_params
   end type t_spawn_params

   type(t_spawn_params) :: spawn_params

! TODO: Make the module private by default
   public
   public :: spawn_params, FMS_Spawn, spawn_couple, FMS_SpawnDCouple, print_spawning_parameters
   private :: propagate_forward_overlap, propagate_backward
   private :: propagate_cont_spawn, propagate_recursive
   private :: FMS_AdjustEnergy2, adjust_child_energy, write_spawn_log

contains

   subroutine initialize_spawn_params(self, CSThresh, CFThresh, PopToSpawn, OMax_inter, OMax_intra, &
                                      OMin_parent, SOCThresh, MaxTraj, MultiSpawn, SpawnCoupV)
      class(t_spawn_params), intent(inout) :: self
      real(kind=DefReal), intent(in) :: CSThresh, CFThresh, PopToSpawn, OMax_inter, &
                                        OMax_intra, OMin_parent, SOCThresh
      integer(kind=DefInt), intent(in) :: MaxTraj, MultiSpawn
      logical, intent(in) :: SpawnCoupV

      self%CSThresh = CSThresh
      self%CFThresh = CFThresh
      self%PopToSpawn = PopToSpawn
      self%OMax_inter = OMax_inter
      self%OMax_intra = OMax_intra
      self%OMin_parent = OMin_parent
      self%SOCThresh = SOCThresh
      self%MaxTraj = MaxTraj
      self%MultiSpawn = MultiSpawn
      self%SpawnCoupV = SpawnCoupV
   end subroutine initialize_spawn_params

   subroutine print_spawning_parameters(unit)
      integer, intent(in) :: unit
      character(len=*), parameter :: divider = &
                                     ' -----------------------------------------------------------'

      write (fmiOut, *)
      write (fmiOut, '(a)') divider
      write (fmiOut, '(a)') ' SPAWNING parameters'
      write (fmiOut, '(a)') divider

      write (unit, '(a14,i7  ,a)') 'maxTraj:     ', spawn_params%MaxTraj, ' (Max number of Trajectories)'
      write (unit, '(a14,f7.4,a)') 'OMin_parent: ', spawn_params%OMin_parent, ' (Min parent-child overlap to spawn)'
      write (unit, '(a14,f7.4,a)') 'OMax_inter:  ', spawn_params%OMax_inter, ' (Max overlap between parent and existing TBFs)'
      write (unit, '(a14,f7.4,a)') 'OMax_intra:  ', spawn_params%Omax_intra, ' (Max overlap between child and existing TBFs)'
      write (unit, '(a14,f7.4,a)') 'CSThresh:    ', spawn_params%CSThresh, ' (Coupling threshold to enter spawning region)'
      write (unit, '(a14,f7.4,a)') 'CFThresh:    ', spawn_params%CFThresh, ' (Coupling threshold to exit spawning region)'
      write (unit, '(a14,f7.4,a)') 'PopToSpawn:  ', spawn_params%PopToSpawn, ' (Min population to Spawn)'

      write (fmiOut, '(a)') divider
   end subroutine print_spawning_parameters

!!    @brief Driver for spawning algorithm
!!
!!    Determines if a spawn is required.  If so, this routine will
!!    perform all the necessary logic to tile the phase space with
!!    virtual trajectories, which are added to the trajectory bundle.
!!    Odds are that this is by far the most complex piece of code in the
!!    FMS program and it is certainly the most crucial routine.
!!
!!    If a timestep is rejected during the spawn, we divide the timestep
!!    in half and retry the spawn.
!!    \param TimeStep Time step to use during spawning
!!    @ingroup spawning

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_Spawn(B1, TimeStep)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  parent_i ---> parent_s ----> parent_f
!                |
!               \/
!  child_i  <--- child_s -----> child_f
!
      use ElecStrucModule
      use SelectionModule, only: FMS_CalculateSelectionTime

      type(T_TrajectoryBundle), target, intent(inout) :: B1
      real(kind=DefReal), intent(in) :: TimeStep

      type(T_Trajectory) :: parent_i, parent_s, parent_f, &
                            child_i, child_s, child_f

      integer(kind=DefInt) :: ntraj, i, j, & ! counters for the trajectories
                              nstate, cs, ps, & ! state indices, parent state, child state are ps, cs
                              npart ! number of particles:

      real(kind=DefReal) :: entry_time, spawn_time, exit_time, &
                            back_time, coup

      logical :: l_overlap, l_spawn

      logical :: COUP_FIELD, COUP_CI ! xf added

      ntraj = B1%NumTraj
      nstate = B1%NumStates
      npart = B1%NumParticles

! Calculate the Mulliken population for the trajectories
! this stores the pop in each Trajectory
      call FMS_UpdateMulliken(B1)

      traj_loop: do i = 1, ntraj

         ps = B1%Trajectory(i)%StateID

         state_loop: do cs = 1, nstate

            ! Don't spawn to the same state
            if (cs == ps) cycle state_loop

            ! Only spawn for Ms=0
            ! DH: I think this condition can be moved to the outer traj_loop,
            ! since it is independent of cs
            if (B1%Trajectory(i)%Ms /= 2) cycle

!  Update coupling history list
!    Note: this spawning subroutine is only called if timestep was accepted, so CoupHist really does
!    contain previous timestep couplings

! xf changed
            if (glzxfaims) then
               COUP_FIELD = .false.
               COUP_CI = .false.

               ! check if we are in the field coupling region
               if (glzxfactive) then
                  coup = spawn_couple_field(B1%Trajectory(i), cs)
                  B1%Trajectory(i)%CoupHist(:, cs) = eoshift(B1%Trajectory(i)%CoupHist(:, cs), -1, coup)
                  COUP_FIELD = .true.
               end if

               ! check if we are in the nonadiabatic coupling region
               if (spawn_couple(B1%Trajectory(i), cs) > spawn_params%CSThresh) then
                  coup = spawn_couple(B1%Trajectory(i), cs)
                  B1%Trajectory(i)%CoupHist(:, cs) = eoshift(B1%Trajectory(i)%CoupHist(:, cs), -1, coup)
                  COUP_CI = .true.
               end if

               ! currently, we cannot handle both field and nonadiabatic couplings
               if (COUP_CI .and. COUP_FIELD) then
                  write (fmiOut, *) 'ERROR! Field and nonadiabatic couplings are present simultaneously. ', &
                     'XFAIMS algorithm not ready yet!'
                  !jj - probably kill the code here or turn off one of the couplings (maybe the field one)
               end if

            else
               coup = spawn_couple(B1%Trajectory(i), cs)
               B1%Trajectory(i)%CoupHist(:, cs) = eoshift(B1%Trajectory(i)%CoupHist(:, cs), -1, coup)
            end if
!xf changed end

            do j = 1, ntraj
               if (B1%Trajectory(j)%StateID == cs) then
                  if (abs(overlap(B1%Trajectory(i), B1%Trajectory(j))) > spawn_params%OMax_inter) then
                     write (fmiOut, *) 'Significant overlap with traj', j
                     cycle state_loop
                  end if
               end if
            end do

!spMultiSpawn :   0 -- Continuous Spawning
!                 1 -- Single Spawning
            if (spawn_params%MultiSpawn == 0) then
               if (spawn_trajectory(B1%Trajectory(i), cs, spawn_params%MultiSpawn, COUP_FIELD, COUP_CI)) then
                  call parent_i%create(npart, nstate)
                  call parent_s%create(npart, nstate)
                  call parent_f%create(npart, nstate)
                  call child_s%create(npart, nstate)
                  call child_i%create(npart, nstate)
                  call child_f%create(npart, nstate)

                  ! create a new trajectory
                  parent_i = B1%Trajectory(i)
                  write (fmiOut, '(20("*"))')
                  write (fmiOut, '(a,i0,a)') 'SPAWNING: Trajectory ', parent_i%TrajID, ' exceeded threshold'
                  write (fmiOut, '(2(a,i2))') 'Parent on state ', ps, ' spawning to state', cs
                  ! new way
                  l_spawn = .true.
                  l_overlap = .true.
                  call propagate_cont_spawn(parent_i, cs, TimeStep, &
                                            parent_s, parent_f, &
                                            child_s, child_f, &
                                            l_spawn, COUP_FIELD, COUP_CI)
                  ! prevent parent spawning to child  state until time=exit_time
                  exit_time = parent_f%get_time()
                  B1%Trajectory(i)%LastSpawn(cs) = exit_time

                  if (l_spawn) then
                     ! backwards propagation of child
                     entry_time = parent_i%get_time()
                     spawn_time = parent_s%get_time()
                     B1%Trajectory(i)%SpawnTime = parent_i%SpawnTime

                     back_time = spawn_time - entry_time
                     child_i = child_s
                     call propagate_backward(child_i, -TimeStep, back_time)

                     ! if there is no overlap with the Bundle, then add the spawn
                     l_overlap = no_overlap_with_bundle(child_i, B1)
                     if (l_overlap) then
                        !Erase Coupling history in child
                        child_i%CoupHist = 0.0d0
                        if (glzStoSwiss) then
                           child_i%SWISS%BirthDate = spawn_time
                           call FMS_CalculateSelectionTime(parent_s, child_s, child_i)
                        end if
                        call B1%add_traj(child_i)
!bfec
                        if (glzStochastic) gldLastSpawnSto = spawn_time

                     else
                        write (fmiOut, *) 'Overlap with exisiting trajectories in Bundle too large'
                        write (fmiOut, *) 'No trajectory will be spawned'
                     end if
                  end if

                  call write_spawn_log(parent_i, parent_s, parent_f, &
                                       child_i, child_s, l_spawn, l_overlap)

                  call parent_i%destroy()
                  call parent_s%destroy()
                  call parent_f%destroy()
                  call child_i%destroy()
                  call child_s%destroy()
                  call child_f%destroy()
               end if

            else if (spawn_params%MultiSpawn == 1) then
               if (spawn_trajectory(B1%Trajectory(i), cs, spawn_params%MultiSpawn, COUP_FIELD, COUP_CI)) then
                  call parent_i%create(npart, nstate)
                  call parent_s%create(npart, nstate)
                  call parent_f%create(npart, nstate)
                  call child_s%create(npart, nstate)
                  call child_i%create(npart, nstate)
                  call child_f%create(npart, nstate)

                  ! create a new trajectory
                  parent_i = B1%Trajectory(i)
                  write (fmiOut, '(20("*"))')
                  write (fmiOut, '(a,i0,a)') 'SPAWNING: Trajectory ', parent_i%TrajID, ' exceeded threshold'
                  write (fmiOut, '(2(a,i2))') 'Parent on state ', ps, ' spawning to state ', cs

                  ! new way
                  l_spawn = .true.
                  l_overlap = .true.
                  call propagate_forward_overlap(parent_i, cs, TimeStep, &
                                                 parent_s, parent_f, &
                                                 child_s, child_f, &
                                                 l_spawn, COUP_FIELD, COUP_CI)
                  ! prevent parent spawning to child  state until time=exit_time
                  exit_time = parent_f%get_time()
                  B1%Trajectory(i)%LastSpawn(cs) = exit_time

                  if (l_spawn) then
                     ! backwards propagation of child
                     entry_time = parent_i%get_time()
                     spawn_time = parent_s%get_time()
                     B1%Trajectory(i)%SpawnTime = parent_i%SpawnTime

                     back_time = spawn_time - entry_time
                     child_i = child_s
                     call propagate_backward(child_i, -TimeStep, back_time)

                     ! if there is no overlap with the Bundle, then add the spawn
                     l_overlap = no_overlap_with_bundle(child_i, B1)
                     if (l_overlap) then
                        !Erase Coupling history in child
                        child_i%CoupHist = 0.0d0
                        if (glzStoSwiss) then
                           child_i%SWISS%BirthDate = spawn_time
                           call FMS_CalculateSelectionTime(parent_s, child_s, child_i)
                        end if

                        !GAIMS change
                        if (child_i%StateID > NSing) then
                           call B1%add_traj_triplet(child_i)
                        else
                           call B1%add_traj(child_i)
                        end if
                        !GAIMS change end

                        !bfec
                        if (glzStochastic) gldLastSpawnSto = spawn_time

                     else
                        write (fmiOut, *) 'Overlap with exisiting trajectories in Bundle too large'
                        write (fmiOut, *) 'No trajectory will be spawned'
                     end if
                  end if

                  call write_spawn_log(parent_i, parent_s, parent_f, &
                                       child_i, child_s, l_spawn, l_overlap)

                  call parent_i%destroy()
                  call parent_s%destroy()
                  call parent_f%destroy()
                  call child_i%destroy()
                  call child_s%destroy()
                  call child_f%destroy()
               end if
            end if

         end do state_loop
      end do traj_loop

   end subroutine FMS_Spawn

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function spawn_trajectory(T1, is, ispawn, COUP_FIELD, COUP_CI) result(spwn)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! result is true if the conditions to spawn to state is are satisfied for the
! trajectory passed in
      type(T_Trajectory), intent(inout) :: T1 ! Trajectory coupling calculated => T1 changes
      integer(kind=DefInt), intent(in) :: is
      integer(kind=DefInt), intent(in) :: ispawn
      logical, intent(in) :: COUP_FIELD, COUP_CI
      logical :: spwn

      spwn = .true.

!DO NOT SPAWN IF:
! The population is below threshold,
      if (T1%Pop < spawn_params%PopToSpawn) spwn = .false.

! We have already spawned in this time period,
      if (T1%get_time() <= T1%LastSpawn(is)) spwn = .false.

! The magnitude of the coupling is too low,
!xf added
      if (glzxfaims) then
         if (COUP_FIELD) then
            if (T1%get_time() < xfaims_params%sp_spwn_i .or. T1%get_time() > xfaims_params%sp_spwn_f) spwn = .false.
         end if
         if (COUP_CI) then
            if (spawn_couple(T1, is) < spawn_params%CSThresh) spwn = .false.
         end if
      else
         if (spawn_couple(T1, is) < spawn_params%CSThresh) spwn = .false.
      end if
! xf added end

      if (ispawn == 1) then
         ! The magnitude of the coupling is decreasing.
         if (T1%CoupHist(1, is) < T1%CoupHist(2, is)) spwn = .false.
      end if

   end function spawn_trajectory

! TODO(DH): Split this function into normal and GAIMS versions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function spawn_couple(T1, is) result(DCouple)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: is

      integer :: js ! the state that the trajectory is currently on

!GAIMS added
      logical :: soc
      complex(kind=DefComp) :: SOM(9)
      real(kind=DefReal) :: DCouple, EDiff, SONorm, au2invcm
!GAIMS end added

      real(kind=DefReal), dimension(T1%NumDimensions) :: vec

      au2invcm = 0.00000455633528
      Dcouple = 0.d0

! No coupling between state and itself
      js = T1%StateID
      if (is == js) return

!GAIMS changed
      soc = .false.
      if (js <= NSing) then
         if (is > NSing) then
            soc = .true.
         end if
      else
         if (is <= NSing) then
            soc = .true.
         end if
      end if

      if (soc) then
         ! It is useless to take all SOMs into account, as most of them are 0 anyways.
         SOM(1) = FMS_SOCoupling(T1, js, is, 1, 1)
         SOM(2) = FMS_SOCoupling(T1, js, is, 1, 2)
         SOM(3) = FMS_SOCoupling(T1, js, is, 1, 3)
         SOM(4) = FMS_SOCoupling(T1, js, is, 2, 1)
         SOM(5) = FMS_SOCoupling(T1, js, is, 2, 2)
         SOM(6) = FMS_SOCoupling(T1, js, is, 2, 3)
         SOM(7) = FMS_SOCoupling(T1, js, is, 3, 1)
         SOM(8) = FMS_SOCoupling(T1, js, is, 3, 2)
         SOM(9) = FMS_SOCoupling(T1, js, is, 3, 3)

         Ediff = T1%ElecStruc%PotEn(is) - T1%ElecStruc%PotEn(js)
         SONorm = sqrt(dreal(SOM(1) * conjg(SOM(1)) + SOM(2) * conjg(SOM(2)) + SOM(3) * conjg(SOM(3)) + &
                             SOM(4) * conjg(SOM(4)) + SOM(5) * conjg(SOM(5)) + SOM(6) * conjg(SOM(6)) + &
                             SOM(7) * conjg(SOM(7)) + SOM(8) * conjg(SOM(8)) + SOM(9) * conjg(SOM(9))))

         !cutoff at 5cm-1
         if (SONorm > spawn_params%SOCThresh * au2invcm) then
            DCouple = abs(SONorm / Ediff)
         else
            DCouple = 0.d0
         end if

      else

         vec = FMS_Coupling(T1, js, is)

         ! Return either coupling .dot. velocity or norm of coupling
         if (spawn_params%SpawnCoupV) then
            DCouple = abs(dot_product(T1%get_vel(), vec))
         else
            DCouple = sqrt(dot_product(vec, vec))
         end if

      end if

!GAIMS end changed
   end function spawn_couple

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine propagate_cont_spawn(parent_i, cs, TimeStep, &
                                   parent_s, parent_f, &
                                   child_s, child_f, success, &
                                   COUP_FIELD, COUP_CI)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(inout) :: parent_i
      integer(kind=DefInt), intent(in) :: cs
      real(kind=DefReal), intent(in) :: TimeStep

      type(T_Trajectory), intent(inout) :: parent_s, parent_f, &
                                           child_s, child_f

      logical, intent(out) :: success
      logical, intent(in) :: COUP_FIELD, COUP_CI

      type(T_Trajectory) :: child_a ! attempt at a spawn

      complex(kind=DefComp) :: S
      real(kind=DefReal) :: coup_max, coup, curr_time, &
                            coup_prev(3), time_tmp1, time_tmp2
      integer(kind=DefInt) :: nstep, &
                              nstate, npart, &
                              ps

      logical :: child_adjusted, child_created
      character :: created_this_step

      success = .false.

      call cpu_time(time_tmp1)

      coup_max = 0.
      child_created = .false.
      nstep = 0

      coup_prev = 1.0d-08

      nstate = parent_i%NumStates
      npart = parent_i%NumParticles
      ps = parent_i%StateID
      curr_time = parent_i%get_time()
      parent_f = parent_i

      call child_a%create(npart, nstate)

      write (fmiOut, *) ' time    coup    overlap  new'
      created_this_step = ' '

      ! attempt to set up the child
      child_a = parent_f
      child_a%StateID = cs
      child_a%ESFlags%ZDerivCurrent = .false.
      !child_a%ESFlags%ZPotEnCurrent = .false.  ! as we have changed state the
      !force is not current
      ! This seems to handled by the
      ! mismatch of StateID
      ! and NForceState. (needs to be
      ! very careful)
      child_a%ParentID = parent_f%TrajID
      child_a%Amplitude = (0., 0.)
      call adjust_child_energy(parent_f, child_a, child_adjusted, COUP_FIELD, COUP_CI)

      if (child_adjusted) then
         ! check the overlap with the parent, it must be big enough
         S = overlap(parent_f, child_a)

         success = .true.
         child_created = .true.
!            success           = .true.
         created_this_step = '*'

         child_a%SpawnTime(ps) = curr_time
         child_s = child_a
         child_f = child_a
         coup_max = coup

         parent_i%SpawnTime(cs) = curr_time
         parent_f%SpawnTime(cs) = curr_time
         parent_s = parent_f
      end if

      ! write out the coupling and the overlap
      if (child_created) then

         S = overlap(parent_f, child_f)
         write (fmiOut, '(f8.2,1x,f8.4,1x,f5.3,1x,a)') parent_f%get_time(), coup, abs(S), created_this_step

      else
         write (fmiOut, '(f8.2,1x,f8.4,1x,a)') parent_f%get_time(), coup, '  NA '

      end if

      call parent_f%set_time(curr_time)
      parent_f%DeadTime = curr_time !This is necessary to prevent parent trajectory from being killed during long spawning event.
      nstep = nstep + 1

      child_s%LastSpawn(ps) = curr_time
      child_f%LastSpawn(ps) = curr_time

      parent_s%LastSpawn(cs) = curr_time
      parent_f%LastSpawn(cs) = curr_time

      call child_a%destroy()

      call cpu_time(time_tmp2)
      ftime = ftime + time_tmp2 - time_tmp1
   end subroutine propagate_cont_spawn

!-------------------------------------------------------------------------------

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine propagate_forward_overlap(parent_i, cs, TimeStep, &
                                        parent_s, parent_f, &
                                        child_s, child_f, success, &
                                        COUP_FIELD, COUP_CI)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(inout) :: parent_i
      integer(kind=DefInt), intent(in) :: cs
      real(kind=DefReal), intent(in) :: TimeStep

      type(T_Trajectory), intent(inout) :: parent_s, parent_f, &
                                           child_s, child_f

      logical, intent(out) :: success
      logical, intent(in) :: COUP_FIELD, COUP_CI

      type(T_Trajectory) :: child_a ! attempt at a spawn

      complex(kind=DefComp) :: S
      real(kind=DefReal) :: coup_max, coup, curr_time, &
                            coup_prev(3), time_tmp1, time_tmp2
      integer(kind=DefInt) :: nstep, &
                              nstate, npart, &
                              ps

      logical :: child_adjusted, child_created
      character :: created_this_step

      success = .false.

      call cpu_time(time_tmp1)

      coup_max = 0.
      child_created = .false.
      nstep = 0

      coup_prev = 1.0d-08

      nstate = parent_i%NumStates
      npart = parent_i%NumParticles
      ps = parent_i%StateID
      curr_time = parent_i%get_time()
      parent_f = parent_i

      call child_a%create(npart, nstate)

      write (fmiOut, *) ' time    coup    overlap  new'
      do
         created_this_step = ' '

! xf changed
         if (glzxfaims) then
            if (COUP_FIELD) then
               coup = spawn_couple_field(parent_f, cs)
            end if
            if (COUP_CI) then
               coup = spawn_couple(parent_f, cs)
            end if
         else
            coup = spawn_couple(parent_f, cs)
         end if
! end xf changed

         coup_prev = eoshift(coup_prev, -1, coup)

         if (all(coup_prev(1) < coup_prev(2:))) then
            write (fmiOut, '(f8.2,1x,f8.4,1x,A)') parent_f%get_time(), coup, '-'
            if (.not. child_created) then
               write (fmiOut, *) 'Spawn failed : coupling peaked but child could not be created'
            else
               write (fmiOut, '(a,g0.2)') 'Spawn successful: child created at t = ', parent_s%get_time()
            end if
            exit
         end if

         if (all(coup_prev(1) > coup_prev(2:))) then

            ! attempt to set up the child
            child_a = parent_f
            child_a%StateID = cs
            child_a%ESFlags%ZDerivCurrent = .false.
            !child_a%ESFlags%ZPotEnCurrent = .false.  ! as we have changed state the force is not current
            ! This seems to handled by the mismatch of StateID
            ! and NForceState. (needs to be very careful)
            child_a%ParentID = parent_f%TrajID
            child_a%Amplitude = (0., 0.)
            call adjust_child_energy(parent_f, child_a, child_adjusted, COUP_FIELD, COUP_CI)

            if (child_adjusted) then
               ! check the overlap with the parent, it must be big enough
               S = overlap(parent_f, child_a)

               if (abs(S) > spawn_params%OMin_parent) then
                  success = .true.
                  child_created = .true.
                  success = .true.
                  created_this_step = '*'

                  child_a%SpawnTime(ps) = curr_time
                  child_s = child_a
                  child_f = child_a
                  coup_max = coup

                  parent_i%SpawnTime(cs) = curr_time
                  parent_f%SpawnTime(cs) = curr_time
                  parent_s = parent_f
               end if
            end if
         end if
         if (glzxfaims) then
! xf changed
            if (COUP_FIELD) then
               if (.not. child_created) then
                  write (fmiOut, *) 'Coup FIELD Spawn fail : coupling dropped below threshold'
                  write (fmiOut, *) 'Coup FIELD              and no child could be created'
                  success = .false.
                  exit
               end if
            end if

            if (COUP_CI) then
               if (.not. child_created .and. coup < spawn_params%CFThresh) then
                  write (fmiOut, *) 'Spawn fail : coupling dropped below threshold'
                  write (fmiOut, *) '             and no child could be created'
                  success = .false.
                  exit
               end if
            end if

         else

            if (.not. child_created .and. coup < spawn_params%CFThresh) then
               write (fmiOut, *) 'Spawn fail : coupling dropped below threshold'
               write (fmiOut, *) '             and no child could be created'
               success = .false.
               exit
            end if

         end if
! xf changed end

         ! write out the coupling and the overlap
         S = overlap(parent_f, child_f)
         if (child_created) then

            write (fmiOut, '(f8.2,1x,f8.4,1x,f5.3,1x,a)') parent_f%get_time(), coup, abs(S), created_this_step

         else

            write (fmiOut, '(f8.2,1x,f8.4,1x,f5.3)') parent_f%get_time(), coup, abs(S)

         end if

         call propagate_recursive(parent_f, TimeStep)
         curr_time = curr_time + TimeStep
         call parent_f%set_time(curr_time)
         ! This is necessary to prevent parent trajectory from being killed during long spawning event.
         parent_f%DeadTime = curr_time
         nstep = nstep + 1
      end do

      child_s%LastSpawn(ps) = curr_time
      child_f%LastSpawn(ps) = curr_time

      parent_s%LastSpawn(cs) = curr_time
      parent_f%LastSpawn(cs) = curr_time

      call child_a%destroy()

      call cpu_time(time_tmp2)
      ftime = ftime + time_tmp2 - time_tmp1
   end subroutine propagate_forward_overlap

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine propagate_backward(T1, TimeStep, TotalTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: TimeStep, TotalTime
      real(kind=DefReal) :: time_tmp1, time_tmp2
      integer(kind=DefInt) :: nstep, n

      call cpu_time(time_tmp1)
! work out the number of time steps to step back
      nstep = int(abs(TotalTime / TimeStep))
      write (fmiOut, '(a,i0,a)') 'Propagating child backwards for ', nstep, ' steps'
      do n = 1, nstep
         write (fmiOut, '(a,g0.2)') ' --Time: ', T1%get_time()
         flush (fmiOut)
         call propagate_recursive(T1, TimeStep)
         call T1%set_time(T1%get_time() + TimeStep)
      end do
      write (fmiOut, *) 'Done backpropagating'

      call cpu_time(time_tmp2)
      btime = btime + time_tmp2 - time_tmp1
   end subroutine propagate_backward

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   recursive subroutine propagate_recursive(T1, TimeStep)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! recursive caller for PropVV.  If the propagation fails the time step is divided in
! two and the propagator recalled
      type(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: TimeStep

      type(T_Trajectory) :: T2

      real(kind=DefReal) :: E1, E2

      gldCurrentTStep = TimeStep

! propagate forward
      call T2%create(T1%NumParticles, T1%NumStates)
      T2 = T1
      call FMS_PropVV(T2, TimeStep)

      E1 = FMS_ClassEnergy(T1)
      E2 = FMS_ClassEnergy(T2)

      if (abs(E2 - E1) > gldEnergyStepCons) then
         write (fmiOut, '(A,g0.4,A,i0)') &
            'WARNING: Energy jumped by ', E2 - E1, &
            ' Hartrees for trajectory #', T2%TrajID
         call FMS_RejectStep(.true.)
      end if
!If current timestep is at min step and step is still rejected then die
      if (FMS_StepRejected() .and. abs(TimeStep + FPZero) < (2.0d0 * gldMinTimeStep)) then
         call FMS_DieError('Propagate_recursive : time step too small')
      end if

      if (FMS_StepRejected()) then
         call FMS_RejectStep(.false.)
         T2 = T1
         call propagate_recursive(T2, TimeStep / 2.d0)
         call propagate_recursive(T2, TimeStep / 2.d0)
      end if

! copy the Trajectory back into 1
      T1 = T2

      call T2%destroy()
   end subroutine propagate_recursive

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine adjust_child_energy(parent, child, success, COUP_FIELD, COUP_CI)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(inout) :: parent, child

      integer(kind=DefInt) :: is, js
      real(kind=DefReal) :: nacme(child%NumDimensions)
      logical, intent(out) :: success
      logical, intent(in) :: COUP_FIELD, COUP_CI

      logical :: adjust_nacme, adjust_mom

      is = parent%StateID
      js = child%StateID

!xf added
      if (glzxfaims) then

! xf changed
         if (COUP_FIELD) then
            nacme = 0.d0
            call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_nacme, nacme)
            if (adjust_nacme) then
               success = adjust_nacme
            else
               call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_mom)
               success = adjust_mom
            end if
         end if

         if (COUP_CI) then
            nacme = FMS_Coupling(child, js, is)
            call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_nacme, nacme)
            if (adjust_nacme) then
               success = adjust_nacme
            else
               call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_mom)
               success = adjust_mom
            end if
         end if
      else
         nacme = FMS_Coupling(child, js, is)
         call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_nacme, nacme)
         if (adjust_nacme) then
            success = adjust_nacme
         else
            call FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, adjust_mom)
            success = adjust_mom
         end if
!xf added end
      end if

   end subroutine adjust_child_energy

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_AdjustEnergy2(child, parent, COUP_FIELD, COUP_CI, success, scale_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    Adjust kinetic energy of child such that total energy of
!!    child and parent are identical.
!!    @ingroup spawning
!!    @param child  Trajectory to be adjusted
!!    @param parent Trajectory providing reference energy
!!    @param scale_in (optional) dictates the direction along which
!!    the momentum adjustment will be applied.  The default is to choose
!!    this direction as the momentum vector itself, i.e. scale all
!!    momentum components equally.
!!    @return Success or failure. Fails only if
!!    kinetic energy becomes negative (frustrated) or
!!    scale factor becomes complex
!!
!<
      type(T_Trajectory), intent(inout) :: child
      type(T_Trajectory), intent(in) :: parent
      logical, intent(inout), optional :: success
      real(kind=DefReal), intent(in), optional :: scale_in(child%NumDimensions)

      logical :: split_momentum ! if true, a component of the momentum
      ! is scaled
      logical, intent(in) :: COUP_FIELD, COUP_CI

      real(kind=DefReal), dimension(child%NumDimensions) :: &
         P, & ! momentum
         scale_vec, & ! component to scale the momentum along
         P_par, & ! Paralell      component of the momentum
         P_per, & ! Perpendicular component of the momentum
         P_new, & ! the new momentum
         mass ! vector of masses

      real(kind=DefReal) :: E_chld, E_prnt, &
                            KE_child, KE_goal, &
                            scale_norm, &
                            KE_par_par, & ! KE of child in terms the
                            KE_par_per, & ! the components of P
                            KE_per_per

      real(kind=DefReal) :: a, b, c, Delta, x ! for solving the KE quadratic

      success = .false.

! calculate classical energy of child and parent
      E_chld = FMS_ClassEnergy(child)
      E_prnt = FMS_ClassEnergy(parent)
!write(fmiOut,'(a,f0.6)') "Energy of child  trajectory: ", E_chld
!write(fmiOut,'(a,f0.6)') "Energy of parent trajectory: ", E_prnt

! check for a frustrated spawn
!xf added
      if (glzxfaims) then
         if (COUP_FIELD) then
! changed xf. KE_goal was the kinetic energy of the child TBF in order to have
! Eparent=Echild.
! Here we do not want Eparent=Echild. In fact the energy difference between the
! two states can be large.
! All we want is KE child = KE parent
            KE_goal = E_prnt - FMS_PotentialT(parent)
            KE_child = FMS_KineticClass(child)
         end if

         if (COUP_CI) then
! check for a frustrated spawn
            KE_goal = E_prnt - FMS_PotentialT(child)
            KE_child = FMS_KineticClass(child)
         end if
      else
         KE_goal = E_prnt - FMS_PotentialT(child)
         KE_child = FMS_KineticClass(child)
      end if
!xf added end

      if (KE_goal < 0) then
         !write(fmiOut,*) "Cannot adjust child kinetic energy to parent."
         return
      end if

! set up the scaling vector, decide if we are splitting the momentum
      if (present(scale_in)) then
         scale_vec = scale_in
         scale_norm = sqrt(dot_product(scale_vec, scale_vec))

         ! check the size of the scaling vector
         if (scale_norm > FPZero) then
            split_momentum = .true.
            scale_vec = scale_vec / scale_norm
         else
            split_momentum = .false.
            write (fmiOut, *) 'WARNING: AdjustEnergy:'
            write (fmiOut, *) 'Scaling vector is zero! Will scale momentum.'
         end if

      else
         split_momentum = .false.
      end if

      P = child%get_mom()

      if (.not. split_momentum) then
         ! simple case, solve: KE_goal = KE . x^2
         x = sqrt(KE_goal / KE_child)

         P_new = x * P

      else
         ! The momentum is split into perpendicular and parallel components
         !     P = P_par + P_per
         P_par = dot_product(P, scale_vec) * scale_vec
         P_per = P - P_par

         ! The kinetic energy in terms of the components is
         !     KE = P/2m . P
         !        = ( P_par + P_per ) /2m . ( P_par + P_per )
         !        = (P_par/2m . P_par) + (P_par/m . P_per) + (P_per/2m . P_per)
         !        = KE_par_par         + KE_par_per        + KE_per_per
         mass = child%get_mass()
         KE_par_par = dot_product(P_par / (2.0 * mass), P_par)
         KE_par_per = dot_product(P_par / mass, P_per)
         KE_per_per = dot_product(P_per / (2.0 * mass), P_per)

         ! We scale the P_par component by x to make the KE equal to KE_goal
         !     (KE_par_par) x^2 + (KE_par_per) x + (KE_per_per) = KE_goal
         !     (KE_par_par) x^2 + (KE_par_per) x + (KE_per_per - KE_goal ) = 0
         ! which is a simple quadratic equation
         a = KE_par_par
         b = KE_par_per
         c = KE_per_per - KE_goal

         Delta = b**2 - 4.*a * c
         if (Delta < 0.) then
            write (fmiOut, *) 'Unable to scale the parallel component'
            return
         end if

         x = (-b + sqrt(Delta)) / (2.*a)

         P_new = x * P_par + P_per
      end if

      call child%set_mom(P_new)
      success = .true.

   end subroutine FMS_AdjustEnergy2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function no_overlap_with_bundle(T, B) result(no_overlap)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! checks that there is no overlap between trajectory T with the trajectories in
! the bundle B (on the same state )
      type(T_Trajectory), intent(in) :: T
      type(T_TrajectoryBundle), intent(in) :: B
      logical :: no_overlap

      integer(kind=DefInt) :: n
      complex(kind=DefComp) :: S

      no_overlap = .true.

      do n = 1, B%NumTraj

         S = overlap(T, B%Trajectory(n), same_state=.true.)

         if (abs(S) > spawn_params%OMax_intra) then
            no_overlap = .false.
            return
         end if
      end do

   end function no_overlap_with_bundle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine write_spawn_log(parent_i, parent_s, parent_f, &
                              child_i, child_s, l_spawn, l_overlap)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: parent_i, child_i, parent_f, &
                                        parent_s, child_s
      logical, intent(in) :: l_spawn, l_overlap

      character(len=256) :: file_name, spawn_sum ! for the spawn XYZ file
! File units for Spawn.log and FailSpawn.log
      integer, save :: u_spawn_log = 0, u_fail_log = 0

      real(kind=DefReal) :: entry_time, spawn_time, exit_time
      integer(kind=DefInt) :: id_c, state_c, &
                              id_p, state_p

      character :: c_spawn, c_overlap

      logical :: file_exists

      entry_time = parent_i%get_time()
      exit_time = parent_f%get_time()
      id_p = parent_i%TrajID
      state_p = parent_i%StateID

      if (l_spawn .and. l_overlap) then
         spawn_time = child_s%get_time()
         id_c = child_i%TrajID
         state_c = child_i%StateID

         ! write the spawned geometry
         write (file_name, '(a,i0)') 'Spawn.', id_c
         write (spawn_sum, '(f10.2,3(a,i4))') spawn_time, ' to state ', state_c, &
            ' from traj ', id_p, &
            ' on state ', state_p
         call FMS_WriteFXYZ(child_s, file_name, spawn_sum)
         ! set up the file and write the header
         if (u_spawn_log == 0) then
            file_name = trim(FMSWorkingDir)//'Spawn.log'
            inquire (file=file_name, exist=file_exists)
            if (file_exists) then
               open (newunit=u_spawn_log, file=file_name, position='append', action='write', status='old')
            else
               open (newunit=u_spawn_log, file=file_name, action='write', status='new')
               write (u_spawn_log, '(A)') '#EntryTime SpawnTime  ExitTime  CID'// &
                  '  CSt  PID  PSt   ChildKE   ChildPE  '// &
                  'ParentKE  ParentPE '
            end if
         end if

         ! write the spawn log
2        format(3(1x, f9.2), 4(1x, i4), 4(1x, f9.5))
         write (u_spawn_log, 2) entry_time, spawn_time, exit_time, &
            id_c, state_c, &
            id_p, state_p, &
            FMS_KineticCLass(child_i), FMS_PotentialT(child_i), &
            FMS_KineticCLass(parent_i), FMS_PotentialT(parent_i)

      else
         ! write the failed spawn log
         if (u_fail_log == 0) then
            file_name = trim(FMSWorkingDir)//'FailSpawn.log'
            inquire (file=file_name, exist=file_exists)
            if (file_exists) then
               open (newunit=u_fail_log, file=file_name, action='write', position='append', status='old')
            else
               open (newunit=u_fail_log, file=file_name, action='write', status='new')
               write (u_fail_log, '(A)') '#Entry time  ExitTime  PID  PSt'// &
                  'ParentKE  ParentPE  spawn fail overlap fail'
            end if
         end if

         c_spawn = 'X'
         c_overlap = 'X'
         if (l_spawn) c_spawn = ' '
         if (l_overlap) c_overlap = ' '

4        format(2(1x, f9.2), 2(1x, i4), 2(1x, f9.5), 9x, a1, 9x, a1)
         write (u_fail_log, 4) entry_time, exit_time, &
            id_p, state_p, &
            FMS_KineticCLass(parent_i), FMS_PotentialT(parent_i), &
            c_spawn, c_overlap
      end if

   end subroutine write_spawn_log

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function spawn_couple_field(T1, is) result(DCouple)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use GlobalModule
      use TrajectoryModule
      type(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: is
      real(kind=DefReal) :: DCouple, t

      integer :: js ! the state that the trajectory is currently on

      Dcouple = 0.d0
!     No coupling between state and itself
      js = T1%StateID

      if (is == js) return

! Computation of the time-dependent electric field.
! "coup" is the electric field in this case.
      t = T1%get_time()

      if (xfaims_params%onespawnonly) then
         Dcouple = abs(exp(-(t - xfaims_params%t0)**2 / (2 * xfaims_params%sigma * xfaims_params%sigma)))
      else
         Dcouple = &
            exp(-(t - xfaims_params%t0)**2 / (2 * xfaims_params%sigma**2)) * ( &
            cos(xfaims_params%freq * t + xfaims_params%CEP) - &
            sin(xfaims_params%freq * t + xfaims_params%CEP) * (t - xfaims_params%t0) &
            / (xfaims_params%sigma**2 * xfaims_params%freq))
      end if

   end function spawn_couple_field

!     TODO(DH): Remove this function in favour of SpawnModule::spawn_couple
!     https://github.com/ispg-group/fms90-redux/issues/79
!>
!!    Computes the interstate coupling dotted with velocity between trajectory T1 and a
!!    trajectory at the same coordinates in state JState.  This calculation
!!    has been moved here in order to make the subroutine FMS_Spawn.f more
!!    readable.
!!
!!    The returned quantity depends on spzSpawnCoupV - if this is
!!    .true., then we return the NACV dotted into the velocity vector.
!!    If not, we just return the norm of the NACV.
!!    \param iState Index of electronic state to calculate coupling with current electronic state
!!    @ingroup spawning
!<
   function FMS_SpawnDCouple(T1, IState) result(DCouple)
      use TrajectoryCalcsModule, only: FMS_CoupDotVel, Kinetic, Potential
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: IState
      real(kind=DefReal) :: DCouple

      integer :: JState ! the state that the trajectory is currently on
      real(kind=DefReal), dimension(T1%NumDimensions) :: vec

      JState = T1%StateID
      Dcouple = 0.d0

!  No coupling between state and itself
      if (T1%StateID == IState) return

!  Return either coupling .dot. velocity or norm of coupling
      if (spawn_params%SpawnCoupV) then
         DCouple = abs(FMS_CoupDotVel(T1, IState))
      else
         vec = FMS_Coupling(T1, JState, IState)
         DCouple = sqrt(dot_product(vec, vec))
      end if

      ! I really don't think this does much good
      if (DCouple /= DCouple) call FMS_DieError('SpawnDCouple: coupling norm is nan')
   end function FMS_SpawnDCouple

!>
!!    Adjust kinetic energy of TChild such that energy of
!!    TChild and TParent are identical.
!!    @ingroup spawning
!!    @param TChild Trajectory to be adjusted
!!    @param TParent Trajectory providing reference energy
!!    @param ScaleVector (optional) dictates the direction along which
!!    the momentum adjustment will be applied.  The default is to choose
!!    this direction as the momentum vector itself, i.e. scale all
!!    momentum components equally.
!!    @return Success or failure. Fails only if
!!    kinetic energy becomes negative (frustrated) or
!!    scale factor becomes complex
!<
! TODO(DH): Figure out the differences between this and FMS_AdjustEnergy2
! We probably don't need both!
   function FMS_AdjustEnergy(TChild, TParent, ScaleVector) result(Success)
      use TrajectoryCalcsModule, only: Kinetic, Potential
      type(T_Trajectory), intent(inout) :: TChild
      type(T_Trajectory), intent(in) :: TParent
      real(kind=DefReal), intent(in), optional :: ScaleVector(TChild%NumDimensions)
      logical :: Success

      real(kind=DefReal), allocatable :: Scaling(:), PPar(:), PPerp(:)
      real(kind=DefReal), allocatable :: Mass(:)
      real(kind=DefReal), allocatable :: PTotal(:)
      real(kind=DefReal) :: Energy, GoalEnergy, ScaleFactor, CurrentKinetic
      real(kind=DefReal) :: DesiredKinetic
      real(kind=DefReal) :: ParPar, PerpPerp, ParPerp, Discriminant, DTemp2
      real(kind=DefReal) :: a, b, c
      integer(kind=DefInt) :: IParticle, IDim, Index
      save Scaling, PPar, PPerp, Mass, PTotal

2003  format('Cannot adjust child kinetic energy to parent.')

!     Dynamic allocation
      Success = .true.
      if (.not. allocated(Scaling)) then
         allocate (Scaling(TChild%NumDimensions))
         allocate (PPar(TChild%NumDimensions))
         allocate (PPerp(TChild%NumDimensions))
         allocate (PTotal(TChild%NumDimensions))
         allocate (Mass(TChild%NumDimensions))
      end if
!
!     Check parent/child energies and write output
!
      Energy = Kinetic(TChild) + Potential(TChild)
      GoalEnergy = Kinetic(TParent) + Potential(TParent)
2000  format('Energy of parent trajectory: ', d12.3)
      write (fmiOut, 2000) GoalEnergy
2002  format('Energy of child trajectory:  ', d12.3)
      write (fmiOut, 2002) Energy
      flush (fmiOut)

      DesiredKinetic = GoalEnergy - Potential(TChild)
      CurrentKinetic = Kinetic(TChild)
!     Frustrated spawn
      if (DesiredKinetic < 0) then
         ScaleFactor = 0
         write (fmiOut, 2003)
         Success = .false.
         return
      end if

      if (.not. present(ScaleVector)) then
         Index = 1
         do IParticle = 1, TChild%NumParticles
            do IDim = 1, TChild%Particle(IParticle)%NumDimensions
               Scaling(Index) = TChild%Particle(IParticle)%get_mom(IDim)
               Index = Index + 1
            end do
         end do
      else
         Scaling = ScaleVector
      end if

      if (sqrt(dot_product(Scaling, Scaling)) < FPZero) then
         write (fmiOut, *) 'WARNING: AdjustEnergy:'
         write (fmiOut, *) 'Scaling vector is zero! Will scale uniformly.'
         Scaling = 1
      end if
      Scaling = Scaling / sqrt(dot_product(Scaling, Scaling))

      Index = 1
      do IParticle = 1, TChild%NumParticles
         do IDim = 1, TChild%Particle(IParticle)%NumDimensions
            PTotal(Index) = TChild%Particle(IParticle)%get_mom(IDim)
            Mass(Index) = TChild%Particle(IParticle)%Mass
            Index = Index + 1
         end do
      end do

      PPar = dot_product(PTotal, Scaling) * Scaling
      PPerp = PTotal - PPar

      ParPar = 0
      ParPerp = 0
      PerpPerp = 0
      Index = 1
      do IParticle = 1, TChild%NumParticles
         do IDim = 1, TChild%Particle(IParticle)%NumDimensions
            ParPerp = ParPerp + PPar(Index) * PPerp(Index) / (2.0 * Mass(Index))
            ParPar = ParPar + PPar(Index) * PPar(Index) / (2.0 * Mass(Index))
            PerpPerp = PerpPerp + PPerp(Index) * PPerp(Index) / (2.0 * Mass(Index))
            Index = Index + 1
         end do
      end do

!     Choose root which leads to minimal change in momentum
      a = ParPar
      b = 2.0 * ParPerp
      c = PerpPerp - DesiredKinetic
      if (abs(a) < FPZero) then
         if (abs(b) > FPZero) then
            ScaleFactor = -c / b
         else
            ScaleFactor = 0
            success = .false.
            write (fmiOut, 2003)
            return
         end if
      else
         Discriminant = b * b - 4 * a * c
         if (Discriminant >= 0) then
            Discriminant = sqrt(Discriminant)
            ScaleFactor = (-b - Discriminant) / (2.0 * a)
            DTemp2 = (-b + Discriminant) / (2.0 * a)
            ScaleFactor = max(ScaleFactor, DTemp2)
         else
            ScaleFactor = 0
            Success = .false.
            write (fmiOut, 2003)
            return
         end if
      end if

      Index = 1
      do IParticle = 1, TChild%NumParticles
         do IDim = 1, TChild%Particle(IParticle)%NumDimensions
            call TChild%Particle(IParticle)%set_mom(iDim, PPerp(Index) + ScaleFactor * PPar(Index))
            Index = Index + 1
         end do
      end do
2001  format('Scaling child momentum by: ', g0.8)
      write (fmiOut, 2001) ScaleFactor

   end function FMS_AdjustEnergy

end module SpawnModule
