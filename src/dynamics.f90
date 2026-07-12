!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!
!!    Propagates a trajectory bundle forward in time from
!!    Bundle%CurrentTime to GoalTime
!!
!!    This routine is recursive; if a timestep is rejected (for
!!    instance, because of the presence of a CI), it will call itself to
!!    integrate that region again with a smaller timestep. If
!!    integration fails through several recursion cycles, the program
!!    will die.
!!
!!    @param Timestep time step for dynamics
!!    @param GoalTime end time
!!    @ingroup propagation
!<
      recursive subroutine FMS_Dynamics(Bundle, Timestep, GoalTime)
         use FMSModule, only: FMS_Shutdown
         use GlobalModule, only: DefReal, DefInt, FPZero, fmiOut, &
        &  FMS_StepRejected, FMS_PrintMessg, FMSWorkingDir, FMS_RejectStep, &
        &  fmzWriteEveryStep, gldCoupTimeStep, &
        &  gldCurrentTStep, gldTimeStep, gldMinTimeStep, &
        &  glzStochastic, glnStepsRejected, FMS_DieError
         use BundleModule
         use BundleIOModule, only: FMS_Output
         use SelectionModule, only: FMS_StochasticCollapse, FMS_RemoveDead
         use RestartModule, only: PutRestart
         use PropagationModule, only: FMS_Monitor
         implicit none
         type(T_TrajectoryBundle), intent(inout) :: Bundle
         real(kind=DefReal), intent(in) :: GoalTime, Timestep

         integer(kind=DefInt), parameter :: RecursionMax = 4
         real(kind=DefReal) :: RedoGoalTime, RedoTimeStep

         logical, save :: FirstTime = .true.
         type(T_TrajectoryBundle), save :: BSave

!     Initialize backup bundle
         if (FirstTime) then
            write (fmiOut, *) 'Saving the first bundle'
            call BSave%create(numtraj=Bundle%NumTraj, &
                              numdeadtraj=Bundle%NumDeadTraj, &
                              numstates=Bundle%NumStates, &
                              numparticles=Bundle%NumParticles, &
                              ncbfs=Bundle%NCBFs)
            FirstTime = .false.
         end if

         do while (Bundle%CurrentTime < Goaltime)

            gldCurrentTStep = TimeStep

!     Save previous timestep in case adaptive propagation is necessary
            BSave = Bundle

!     Write start of timestep
            write (fmiOut, *)
2130        format(' --Time: ', f10.2, ' Timestep: ', f6.2, ' (', i0, ' trajectory)')
2131        format(' --Time: ', f10.2, ' Timestep: ', f6.2, ' (', i0, ' trajectories)')

            if (Bundle%NumTraj == 1) then
               write (fmiOut, 2130) Bundle%CurrentTime, TimeStep, Bundle%NumTraj
            else
               write (fmiOut, 2131) Bundle%CurrentTime, TimeStep, Bundle%NumTraj
            end if
            flush (fmiOut)

!     Propagate bundle through 1 timestep
            call PropagateBundle(Bundle, Timestep)

!     Check for energy and norm conservation, and make sure
!     timestep was appropriate
            call FMS_Monitor(Bundle, BSave, TimeStep)

!     Clean up bundle. This is done here to avoid updatePES being called on
!     dead trajectories
            if (.not. FMS_StepRejected()) then
               if (glzStochastic) then
                  !Stochastic collapse nuclear wavefunction
                  !(RemoveDead must be called immediately after)
                  call FMS_StochasticCollapse(Bundle)
               end if
               call FMS_RemoveDead(Bundle)
               if (Bundle%NumTraj == 0) return
            end if

!     Write output here only if user requested output at every step
            if (fmzWriteEveryStep .and. (.not. FMS_StepRejected())) then
               call FMS_Output(Bundle)
            end if

!     If timestep was rejected, redo it with a smaller timestep
!
            if (FMS_StepRejected()) then
               glnStepsRejected = glnStepsRejected + 1

               ! Set timestep and integration time for the retry
               ! Looking for equality between real number is BAD IDEA!!
               RedoGoalTime = Bundle%CurrentTime
               if (Timestep == gldTimeStep .and. gldTimestep /= gldCoupTimeStep) then
                  RedoTimeStep = gldCoupTimeStep
               else
                  RedoTimeStep = TimeStep / 2.0d0
               end if

!           If the timestep is too small, crash
               if (RedoTimeStep + FPZero < gldMinTimeStep) then
                  call RecursionCrash()
               end if

!           Write output
2000           format('Timestep rejected, will integrate until t =', f0.3)
2001           format('with timestep of', f7.3)
               write (fmiOut, 2000) RedoGoalTime
               write (fmiOut, 2001) RedoTimeStep

!           Reset bundle and re-run
               Bundle = BSave
               call FMS_RejectStep(.false.)
               call FMS_Dynamics(Bundle, RedoTimeStep, RedoGoalTime)

!           If we're here, redo was succesful
2141           format('Adaptive integration until ', f7.3, ' successful.')
               write (fmiOut, 2141) Bundle%CurrentTime
               flush (fmiOut)

            end if

            ! save Bundle
            call PutRestart(Bundle)
            call FMS_CheckStop(Bundle)
            flush (fmiOut)

            ! end of MD loop
         end do

      contains

!>
!!    Top-level propagation subroutine.
!!    Propagates simultaneously both classical and quantum degrees
!!    of freedom of a trajectory bundle.
!!    \note Both FullyCoupled (FullyCoupled==True) and Split-Operator propagation
!!    (FullyCoupled==False) are supported
!<
         subroutine PropagateBundle(Bundle, TimeStep)
            use GlobalModule, only: glcIntegType, GlzFullyCoupled
            use PropagationModule, only: FMS_PropSplitOperator, FMS_PropQCVV
            type(T_TrajectoryBundle), intent(inout) :: Bundle
            real(kind=DefReal), intent(in) :: TimeStep

!     Split-Operator Option
            if (.not. GlzFullyCoupled) then

               call FMS_PropSplitOperator(Bundle, TimeStep)

!     Fully-Coupled propagation
            else

               select case (glcIntegType)
               case ('VV') ! Velocity Verlet
                  call FMS_PropQCVV(Bundle, TimeStep)
               case default
                  call FMS_DieError('Invalid integrator: '//glcIntegType)
               end select

            end if

         end subroutine PropagateBundle

!>
!!    Check for the presence of "stop" file that indicates the user
!!    wants to stop the simulation.
!<
         subroutine FMS_CheckStop(Bundle)
            implicit none

            type(T_TrajectoryBundle), intent(inout) :: Bundle
            character(len=*), parameter :: STOP_FILE_NAME = 'stop'
            character(len=256) :: FilePath
            logical :: zExist

            FilePath = trim(FMSWorkingDir)//STOP_FILE_NAME

            inquire (file=FilePath, exist=zExist)

            if (zExist) then
               write (fmiout, *) 'Found a "stop" file. Shutting down now!'
               call FMS_Shutdown(Bundle)
            end if
         end subroutine FMS_CheckStop
!>
!!    Error handler for FMS_Dynamics
!<
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         subroutine RecursionCrash()
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            call FMS_PrintMessg('ERROR in FMS_Dynamics.')
            call FMS_PrintMessg('Timestep below minimum timestep.')
            call FMS_PrintMessg('Simulation will stop ...')
            call BSave%destroy()
            call FMS_Shutdown(B1=Bundle)

         end subroutine RecursionCrash

      end subroutine FMS_Dynamics
