!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Parameters and data for generating initial conditions
!<
module InitialModule
   use GlobalModule, only: DefInt, DefReal, DefInt4, fmiOut
   implicit none
   public
   save

!> Initial state to create wavefunctions on
   integer(kind=DefInt) :: inInitState

!> Random seed
   integer(kind=DefInt4) :: iniRndSeed ! this is tied to fms_ranb

!> Create mirror basis?
   logical :: inzMirrorBasis
!> State on which to create mirror basis
   integer(kind=DefInt) :: inMirrorState

!> If true, place initial basis functions on the state with the higher
!! transition dipole moment, either inInitState or inIDark
   logical :: inzInitBright
!> If true, place initial basis functions on the state with the lower
!! transition dipole moment, either inInitState or inIDark
   logical :: inzInitDark

!> Other state to put basis functions on, for inzInitBright or inzInitDark
   integer(kind=DefInt) :: inIDark

!> If greater than zero, require initial trajectories to have a certain ground-excited energy gap
   real(kind=DefReal) :: inInitGap, inInitGapWidth

!> Select State ID based on energy gap
   logical :: inSelectState

   logical :: inzEnergyAdjust
   logical :: inzNormInitial
   logical :: inzSharpEnergy
   logical :: inzOSAmp

!     QM/MM
! TODO: Remove these!

!      Equilibration and system prep controls
   integer(kind=DefInt) :: innEquiStepPrint !< During equilibration, print every this many steps
   integer(kind=DefInt) :: inNCycles !< Number of velocity resamplings during equicon
   integer(kind=DefInt) :: inNSteps !< Equilibration timesteps
   logical :: inzGenSolvent !< Randomly generate solvent positions?
   integer(kind=DefInt) :: inNRelaxSteps !< Relax MM system for this many steps
   integer(kind=DefInt) :: inNFixSteps !< Relax MM system while keeping QM fixed for this many steps
   real(kind=DefReal) :: indEquilTStep !< Timestep for equilibration
   logical :: inzEquiRes !< is this a restart of an equilibration?
   integer(kind=DefInt) :: inIBrown !< The brownian force will be applied every inIBrown steps
   real(kind=DefReal) :: inDGamma !< Friction coefficient for brown
   real(kind=DefReal) :: indQuenchToler ! Convergence threshold for MM relaxation

!     TODO: Remove these!
!     Equilibration procedures
   logical :: inzBrown !< Brownian dynamics for equilibration
   logical :: inzBrownCon !< Constrained brownian dynamics?
   logical :: inzEquiCon !< RATTLE constrained QM
   logical :: inzEqui !< Boltzmann resampling equilibration?

contains

!>
!!    Top-level subroutine to apply initial conditions to QM region of a trajectory.
!!
!!    Calls the appropriate initial conditions routine
!!    depending on the value of inInitialCond
!!    @ingroup initial
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_SelectInitial(B1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use GlobalModule, only: DefReal, fmiOut, FmsWorkingDir, FMS_DieError
      use TrajectoryModule
      use TrajectoryCalcsModule, only: Potential
      use BundleModule
      use QM_MM_Module, only: qczQMMM
      use SamplingModule
      use FMSModule, only: FMS_ReadGeometry
      use VerletModule, only: FMS_PropVV
      implicit none
      type(T_TrajectoryBundle), intent(inout) :: B1
      type(T_Trajectory) :: TTemp

      integer :: ntraj, n

      logical :: zExist, redo
      real(kind=DefReal) :: Etmp

      ntraj = B1%NumTraj

      call TTemp%create(B1%Trajectory(1)%NumParticles, B1%NumStates)

      ! Select sampling type
      select case (inInitialCond)
         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! no sampling
      case (NOSAMPLE)
         write (fmiOut, *) 'No sampling for initial conditions'

         if (qczQMMM) then

            ! QM/MM methods generally use the MM code's input file, rather than
            ! Geometry.dat
            ! For InitialCond==NOSAMPLE, we'll grab the initial trajectory
            ! state from Geometry.dat (this will already have been done for QM jobs)
            inquire (file=trim(FMSWorkingDir)//'Geometry.dat', exist=zExist)
            if (zExist) then
               write (fmiOut, *) 'Reading initial conditions from Geometry.dat'
               call FMS_ReadGeometry(B1%Trajectory(1))
!              It is assumed that the user has a good geometry, so don't check
!              energy gap criterion here
            end if
         end if

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! Wigner sampling at 0 or finite temperature
      case (WIGNER)
         write (fmiOut, *) 'Wigner Sampling'

         do n = 1, ntraj
            redo = .true.
            !For now, a hacky fix
            Etmp = Potential(B1%Trajectory(n)) !do ES calc if needed
            do while (redo)
               TTemp = B1%Trajectory(n) !There was a problem here, B1
               !doesn't yet have ES, so in second
               !cycle of this loop, TC would crash
               call FMS_InitialWigner(TTemp, n)
               call FMS_CheckInitial(TTemp, redo)
            end do
            B1%Trajectory(n) = TTemp
         end do

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! Husimi sampling at 0 or finite temperature
         ! Desigened for thermal Gaussian approximation
      case (HUSIMI)
         write (fmiOut, *) 'Husimi sampling with optimum widths at finite temperature for initial conditions'

         do n = 1, ntraj
            redo = .true.
            do while (redo)
               TTemp = B1%Trajectory(n)
               call FMS_InitialHusimi(TTemp, n)
               call FMS_CheckInitial(TTemp, redo)
            end do
            B1%Trajectory(n) = TTemp
         end do

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! Quasi-classical (action/angle)
      case (QUASICLASS)
         write (fmiOut, *) 'Quasi-classical initial conditions'
         do n = 1, ntraj
            redo = .true.
            do while (redo)
               TTemp = B1%Trajectory(n)
               call FMS_InitialQuasi(TTemp)
               call FMS_CheckInitial(TTemp, redo)
            end do
            B1%Trajectory(n) = TTemp
         end do

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! Boltzmann at finite temperature
      case (BOLTZ)
         write (fmiOut, *) 'Boltzmann sampling at T =', IndTemperature
         do n = 1, ntraj
            redo = .true.
            do while (redo)
               TTemp = B1%Trajectory(n)
               call FMS_InitialConstT(TTemp)
               call FMS_CheckInitial(TTemp, redo)
            end do
            B1%Trajectory(n) = TTemp
         end do

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! cold Wigner
      case (COLDWIG)
         write (fmiOut, *) '"Cold Wigner" at T =', indTemperature
         do n = 1, ntraj
            redo = .true.
            do while (redo)
               TTemp = B1%Trajectory(n)
               call FMS_InitialColdWig(TTemp, n)
               call FMS_CheckInitial(TTemp, redo)
            end do
            B1%Trajectory(n) = TTemp
         end do

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! excited state trajectories
      case (SWARM)
         write (fmiOut, *) 'Generating a swarm from Geometry.dat'
         call FMS_InitialSwarm(B1)

         ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
         ! Error message
      case default
         call FMS_DieError('Initial condition not recognized.')

      end select
      write (fmiOut, *)

      call TTemp%destroy()

   contains

!>
!!    Checks a Trajectory's initial condition by applying constraints, if requested,
!!    selecting bright state and checking energy gaps.
!!
!!    @ingroup initial
!<
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine FMS_CheckInitial(T1, redo)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         use GlobalModule
         use TrajectoryModule
         use TrajectoryCalcsModule, only: Potential
         use VerletModule, only: FMS_PropVV
         implicit none
         type(T_Trajectory), intent(inout) :: T1
         logical, intent(out) :: redo
         real(kind=DefReal) :: EGap(T1%NumStates)
         logical :: statemask(T1%NumStates)
         integer :: startstate, endstate, i

!     Default behaviour is to accept initial condition
         redo = .false.

!     First, apply any constraints, if requested.
         if (glzConstrain) then
            call FMS_PropVV(T1, +gldTimeStep)
            call FMS_PropVV(T1, -gldTimeStep)
         end if

!     Now, occupy bright or dark state, if requested.
         if (inzInitBright .or. inzInitDark) then
            call FMS_CalcBrightDark(T1)
         end if

!     Now check energy gap, if requested
         if (inInitGap > FPZero) then

! Do we want to consider all states, or just the current (bright?) state?
            if (inSelectState) then
               startstate = 2
               endstate = T1%NumStates
            else
               startstate = T1%StateID
               endstate = T1%StateID
            end if

            EGap(:) = 9.9d9 !should avoid accidental energy condition matching
            do i = startstate, endstate
! Calculate actual energy gap (will invoke ES the first time)
               EGap(i) = Potential(T1, i) - Potential(T1, 1)
            end do

!WJG for now, just use square
!        SELECT CASE (inzInitGapShape)
!        case ("square")
            write (fmiOut, *) 'Desired Energy Gap:', inInitGap, '+/-', 0.5d0 * inInitGapWidth, ' Hartree.'
            write (fmiOut, *) 'Initial Energy Gaps:', EGap(startstate:endstate), ' Hartree.'

            statemask(:) = .false.
            statemask = (abs(EGap - inInitGap) < 0.5d0 * inInitGapWidth)
            if (all(.not. statemask)) then
               write (fmiOut, *) 'Resampling...'
               redo = .true.
               !Figure out how many states met the initial condition
            else
               if (count(statemask) > 1) then
                  write (fmiOut, *) 'ERROR: ', count(statemask), 'states met the energy condition'
                  write (fmiOut, *) 'Need to implement multiple initial trajectories'
                  stop 1
               end if
               write (fmiOut, *) 'Initial condition met by state:'
               do i = startstate, endstate
                  if (statemask(i)) then
                     write (fmiOut, *) i
                     write (fmiOut, *) 'Setting new StateID'
                     T1%StateID = i
                  end if
               end do
               !Figure out how many states met the initial condition
            end if
!        case ("gaussian")
!          write(fmiOut,*)'Desired Energy Gap:',inInitGap,' FWHM',
!     &                    inInitGapWidth,' Hartree.'
!! Sample a normal distributed energy gap from laser properties
!        if (samplegapfind) then
!1         continue
!            dx1=2.d0*FMS_ranb(i4zero)-1.d0
!            dx2=2.d0*FMS_ranb(i4zero)-1.d0
!            rsq=dx1*dx1+dx2*dx2
!            if(rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
!            fac=sqrt(-2.d0*log(rsq)/rsq)
!            x1=abs(dx1*fac)*InitGapSigma_inv_au
!            write(fmiOut,*)'Sampled experimental Gap:',inInitGap,'+/-',
!        endif
!          write(fmiOut,*)'Initial Energy Gap:',EGap,' Hartree.'
!          If (abs(EGap-inInitGap).gt.x1) then
!          write(fmiOut,*)'Resampling...'
!          redo=.true.
!        else
!          write(fmiOut,*)'Initial condition accepted!'
!        endif
         end if

      end subroutine FMS_CheckInitial

   end subroutine FMS_SelectInitial

!!    Assigns the initial electronic state ID for a trajectory
!!
!!    Assign the brighter or darker (based on the transition dipole moment
!!    as requested) of InitState or inIDark
!!
!!    @ingroup initial
!<
   subroutine FMS_CalcBrightDark(T1)
      use TrajectoryModule
      use TrajectoryCalcsModule, only: FMS_TransDipole
      type(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal) :: Tmp1, Tmp2

!     Measure transition dipole moments
      Tmp2 = sum(FMS_TransDipole(T1, T1%StateID)**2)
      Tmp1 = sum(FMS_TransDipole(T1, inIDark)**2)
      write (fmiOut, *)
1000  format(' Transition dipole magnitude for state ', i3, ':', d15.5)
      write (fmiOut, 1000) inIDark, sqrt(tmp1)
      write (fmiOut, 1000) T1%StateID, sqrt(tmp2)

!     If we want bright state, change state index if current choice is dark
      if (inzInitBright .and. Tmp1 > Tmp2) then
         T1%StateID = inIDark
      end if

!      If we want dark state, change state index if current choice is bright
      if (inzInitDark .and. Tmp1 < Tmp2) then
         T1%StateID = inIDark
      end if

2000  format(' --- InitBrightDark: Selected state ', i3, ' for trajectory #', i4)
      write (fmiOut, 2000) T1%StateID, T1%TrajID
      write (fmiOut, *)

   end subroutine FMS_CalcBrightDark
!>
!!
!!    Creates a mirror basis with zero population on the other electronic
!!    state
!!
!!    @ingroup initial
!<
   subroutine FMS_SpawnMirrorMP(Bundle)
      use TrajectoryModule
      use TrajectoryCalcsModule, only: FMS_Coupling
      use BundleModule
      use SpawnModule, only: FMS_AdjustEnergy
      implicit none
      type(T_TrajectoryBundle), intent(inout) :: Bundle
      type(T_Trajectory) :: SpawnTrajectory
      logical :: MErr
      real(kind=DefReal), allocatable :: CoupVec(:)
      integer(kind=DefInt) :: ITraj, NumTraj, I

      allocate (CoupVec(Bundle%Trajectory(1)%Numdimensions))
      call SpawnTrajectory%create(Bundle%Trajectory(1)%NumParticles, Bundle%NumStates)

      NumTraj = Bundle%NumTraj
      do ITraj = 1, NumTraj

         SpawnTrajectory = Bundle%Trajectory(ITraj)
         call SpawnTrajectory%set_time(Bundle%CurrentTime)

!        Change electronic state index, adjust momentum appropriately,
!        and add the trajectory to the current bundle.
         SpawnTrajectory%SpawnMode = 0
         SpawnTrajectory%StateID = inMirrorState
         SpawnTrajectory%SpawnMode(SpawnTrajectory%StateID) = -1
         SpawnTrajectory%Phase = 0
         SpawnTrajectory%Amplitude = 0

!        Need an estimate of the energy adjustment direction.
         if (inzEnergyAdjust) then
            do i = 1, Bundle%Trajectory(1)%Numdimensions
               CoupVec(i) = 0
            end do
            CoupVec = FMS_Coupling(SpawnTrajectory, SpawnTrajectory%StateID, Bundle%Trajectory(ITraj)%StateID)
            MErr = FMS_AdjustEnergy(SpawnTrajectory, Bundle%Trajectory(ITraj), ScaleVector=CoupVec)
         end if
         call Bundle%add_traj(SpawnTrajectory)
         Bundle%Trajectory(ITraj)%SpawnMode(inMirrorState) = -1
         Bundle%Trajectory(ITraj)%LastSpawn = Bundle%CurrentTime
      end do

      call SpawnTrajectory%destroy()
      deallocate (CoupVec)

   end subroutine FMS_SpawnMirrorMP
end module InitialModule
