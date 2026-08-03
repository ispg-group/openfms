! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module PropagationModule
   use FMSModule
   use GlobalModule, only: DefInt, DefReal, DefComp, DP5, FPZero, &
                           glzConstrain, glzCentroids, glcIntegType, glzFullyCoupled, &
                           FMS_DieError
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_GetForce, FMS_PhaseDot, FMS_PosDot, FMS_MomDot
   use BundleModule
   use BundleCalcsModule, only: FMS_UpdateCentroid
   use VerletModule, only: FMS_PropVV_a, FMS_PropVV_b
   use SpawnModule, only: FMS_Spawn, FMS_SpawnDCouple, spawn_params
   implicit none

   private
   public :: FMS_PropSplitOperator, FMS_PropQCVV
   public :: FMS_Monitor, FMS_SetTimeStep

contains

!-------------------------------------------------------------
!     Split-Operator propagation
!-------------------------------------------------------------
   subroutine FMS_PropSplitOperator(Bundle, TimeStep)
      use BundleCalcsModule, only: FMS_Norm, FMS_Branching
      type(T_TrajectoryBundle), intent(inout) :: Bundle
      real(kind=DefReal), intent(in) :: TimeStep
      real(kind=DefReal) :: HalfTimeStep
      real(kind=DefReal) :: OldNorms(Bundle%NumStates)
      integer(kind=DefInt) :: ITraj

      real(kind=DefReal), save :: dNorm0 = 0.0d0

      HalfTimeStep = dp5 * TimeStep

!-------------------------------------------------------------
!     Propagate each trajectory classically for 1/2 TimeStep
!-------------------------------------------------------------
      if (dNorm0 == 0.0d0) dNorm0 = FMS_Norm(Bundle)

      call FMS_Branching(Bundle, OldNorms)

      do ITraj = 1, Bundle%NumTraj
         call FMS_PropClassNew(Bundle%Trajectory(ITraj), HalfTimeStep)
      end do

      ! Renormalize if multi-state problem
      if (Bundle%Trajectory(1)%NumStates /= 1) then
         call FMS_Renormalize(Bundle, OldNorms, DNorm0)
      end if

!-------------------------------------------------------------
!     Quantum Coefficients are only propagated if this is a multi-state
!     problem
!-------------------------------------------------------------
      if (Bundle%Trajectory(1)%NumStates /= 1) then
         ! Propagate Trajectory Amplitudes for TimeStep
         call FMS_PropQuantum(Bundle, TimeStep)
      end if

!-------------------------------------------------------------
!     Propagate each trajectory classically for 1/2 TimeStep
!-------------------------------------------------------------

      call FMS_Branching(Bundle, OldNorms)

      do ITraj = 1, Bundle%NumTraj
         call FMS_PropClassNew(Bundle%Trajectory(ITraj), HalfTimeStep)
      end do

!     Renormalize if multi-state problem

      if (Bundle%Trajectory(1)%NumStates /= 1) then
         call FMS_Renormalize(Bundle, OldNorms, DNorm0)
      end if

   end subroutine FMS_PropSplitOperator

!>
!!    Propagates trajectory bundle forward in time, with
!!    classical parameters (including phase) propagated with velocity Verlet
!!    and quantum mechanical amplitudes propagated using adaptive 2nd-order Runge-Kutta.
!!    \param Bundle Trajectory bundle to propagate
!!    \param TimeStep Total time step by which to propagate system forward
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_PropQCVV(Bundle, TimeStep)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use GlobalModule, only: FMS_StepRejected
      use BundleCalcsModule, only: FMS_bHeff, FMS_Get_Amplitude, &
                                   FMS_bSp5i, FMS_Set_Amplitude
      type(T_TrajectoryBundle), intent(inout) :: Bundle
      real(kind=DefReal), intent(in) :: TimeStep

      type(T_TrajectoryBundle), save :: BSave

      type(T_Trajectory) :: LastTraj

      integer(kind=DefInt) :: ntraj0, i, j, ntraj1

      real(kind=DefReal), dimension(Bundle%NumTraj) :: g1_0, g2_0

! GAIMS added
!     Temporary structures to hold the information which needs to be
!     preserved when copying the Ms=0 trajectory to the other ones of
!     the same CBF
      integer(kind=DefInt) :: id, ms, iCBF, jCBF
      real(DefReal), dimension(3, Bundle%NumStates) :: cohist
      complex(DefComp) :: amp
! GAIMS end added

      ntraj0 = Bundle%NumTraj

      BSave = Bundle
      call FMS_Set_Amplitude(Bundle, &
                             Prop_H(FMS_bHeff(Bundle), FMS_Get_Amplitude(Bundle), TimeStep / 2.d0, 10) &
                             )

! GAIMS changed
      do i = 1, ntraj0

         ! only propagate once for each CBF (for Ms=0)
         if (Bundle%Trajectory(i)%Ms /= 2) cycle

         call FMS_PropVV_a(Bundle%Trajectory(i), LastTraj, TimeStep, g1_0(i), g2_0(i))

         ! copy the calculated trajectory information to the other BF of the
         ! contracted triplet
         ! only safe information specific to the sublevels.
         ! If only certain SOM elemens are read, this should be changed so
         ! that SOMat is again preserved, together with ZSOMCurrent!!!!
         ! TODO(DH): Might be good to separate this into a function, since this code
         ! is repeated below
         if (Bundle%Trajectory(i)%triplet) then
            do j = 1, ntraj0
               if (Bundle%Trajectory(i)%CBF == Bundle%Trajectory(j)%CBF .and. &
                   Bundle%Trajectory(i)%TrajID /= Bundle%Trajectory(j)%TrajID) then
                  id = Bundle%Trajectory(j)%TrajID
                  amp = Bundle%Trajectory(j)%Amplitude
                  ms = Bundle%Trajectory(j)%Ms
                  ! som=Bundle%Trajectory(j)%ElecStruc%SOMat
                  cohist = Bundle%Trajectory(j)%CoupHist

                  call Bundle%Trajectory(j)%copy_from(Bundle%Trajectory(i))
                  Bundle%Trajectory(j)%TrajID = id
                  Bundle%Trajectory(j)%Amplitude = amp
                  Bundle%Trajectory(j)%Ms = ms
                  ! Bundle%Trajectory(j)%ElecStruc%SOMat=som
                  Bundle%Trajectory(j)%CoupHist = cohist
               end if
            end do
         end if

         ! ntraj0 end do
      end do
! GAIMS end changed

!Positions have been updated, so now update centroids here.
!  Note: parents dont have updated ES at this point, but each centroid
!  should have previous ES information
!bfec
      if (glzCentroids) then
         ! Centroids only updated for traj with Ms=2, so only once per CBF
         do i = 1, ntraj0
            do j = 1, i - 1
               if (Bundle%Trajectory(i)%Ms == 2 .and. Bundle%Trajectory(j)%Ms == 2) then
                  iCBF = Bundle%Trajectory(i)%CBF
                  jCBF = Bundle%Trajectory(j)%CBF
                  call FMS_UpdateCentroid(Bundle%Trajectory(i), Bundle%Trajectory(j), &
                                          Bundle%Centroids(((iCBF - 2) * (iCBF - 1)) / 2 + jCBF))
               end if
            end do
         end do
      end if

      do i = 1, ntraj0
         ! only propagate once for each CBF (for Ms=0)
         if (Bundle%Trajectory(i)%Ms /= 2) cycle

         call FMS_PropVV_b(Bundle%Trajectory(i), LastTraj, TimeStep, g1_0(i), g2_0(i))

         if (Bundle%Trajectory(i)%triplet) then
            do j = 1, ntraj0
               if (Bundle%Trajectory(i)%CBF == Bundle%Trajectory(j)%CBF .and. &
                   Bundle%Trajectory(i)%TrajID /= Bundle%Trajectory(j)%TrajID) then
                  id = Bundle%Trajectory(j)%TrajID
                  amp = Bundle%Trajectory(j)%Amplitude
                  ms = Bundle%Trajectory(j)%Ms
!                 som=Bundle%Trajectory(j)%ElecStruc%SOMat
                  cohist = Bundle%Trajectory(j)%CoupHist

                  call Bundle%Trajectory(j)%copy_from(Bundle%Trajectory(i))

                  Bundle%Trajectory(j)%TrajID = id
                  Bundle%Trajectory(j)%Amplitude = amp
                  Bundle%Trajectory(j)%Ms = ms
!                 Bundle%Trajectory(j)%ElecStruc%SOMat=som
                  Bundle%Trajectory(j)%CoupHist = cohist
               end if
            end do
         end if
      end do

      ! Increment time
      Bundle%CurrentTime = Bundle%CurrentTime + TimeStep
      do i = 1, ntraj0
         call Bundle%Trajectory(i)%set_time(Bundle%CurrentTime)
      end do

      ! check for errors
      call FMS_Monitor(Bundle, BSave, TimeStep)

      if (FMS_StepRejected()) return

      ! TODO: In order to decouple spawning and propagation module, we
      ! could stop the subroutine here, and instead call spawn in the caller
      ! Spawn - this will be done in serial on the master for now as it is not readily parallelizable.
      ! amv: We spawn here in order to calculate the full Hamiltonian at t+dt
      call FMS_Spawn(Bundle, Timestep)
      ntraj1 = Bundle%NumTraj

!Update newly created centroids here.
!bfec
      if (glzCentroids) then
         do i = ntraj0 + 1, ntraj1
            do j = 1, i - 1
               if (Bundle%Trajectory(i)%Ms == 2 .and. Bundle%Trajectory(j)%Ms == 2) then
                  iCBF = Bundle%Trajectory(i)%CBF
                  jCBF = Bundle%Trajectory(j)%CBF
                  call FMS_UpdateCentroid(Bundle%Trajectory(i), Bundle%Trajectory(j), &
                                          Bundle%Centroids(((iCBF - 2) * (iCBF - 1)) / 2 + jCBF))
               end if
            end do
         end do
      end if
!     Update ES for new centroids.
!     Note - these centroids will have duplicated ES calculations
!     at the beginning of the next timestep when UpdateCentroid gets called for them again.
!     One could add some logic to prevent this, but that seems overly cumbersome for a small performance hit
!     (the new centroid calculations only get doubled on the timestep immediately after a spawn)

      call FMS_Set_Amplitude(Bundle, &
                             Prop_H(FMS_bHeff(Bundle), FMS_Get_Amplitude(Bundle), TimeStep / 2.d0, 10) &
                             )

   contains

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function Prop_H(H, C_t, dt, nslice) result(C_tdt)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Solve:
         !  d/dt C = -i H C
         !
         ! Solution:
         !  C(t+dt) = exp( -i H(t) dt ) C(t)
         !  C_t     = exp( B )          Ct_tdt
         !
         ! Basic property of expontial:
         !   exp(B) = exp( B/n ) ** n
         !          = exp( Bn  ) ** n
         ! This reduces the size of the exponential we must take.
         !
         ! The exponential is written Taylor series expansion to 4th order
         !    exp(Bn) = I + Bn + 1/2 B2**2 + 1/3! Bn**3 + 1/4! Bn**4
         !
         ! n is varied until C_tdt is stable
         complex(kind=DefComp), intent(in) :: C_t(:)
         complex(kind=DefComp), intent(in) :: H(size(C_t), size(C_t))
         real(kind=DefReal), intent(in) :: dt
         integer(kind=DefInt), intent(in) :: nslice

         complex(kind=DefComp) :: C_tdt(size(C_t))

         complex(kind=DefComp) :: C_tdt_prev(size(C_t))

         integer(kind=DefInt) :: n, ndim, i

         complex(kind=DefComp), dimension(size(C_t), size(C_t)) :: Id, B, Bn, Bn2, Bn3, Bn4, tayl, tayl_prod

         real(kind=DefReal) :: error

         ndim = size(C_t)
         Id = (0., 0.)
         do n = 1, ndim
            Id(n, n) = (1., 0.)
         end do

         B = -(0., 1.) * H * dt

!     Integrate with increasing numbers of intermediate steps
!     until the norm converges

         C_tdt_prev = (0.d0, 0.d0)
         do n = 1, nslice

            Bn = B / 2**n
            Bn2 = matmul(Bn, Bn)
            Bn3 = matmul(Bn2, Bn)
            Bn4 = matmul(Bn2, Bn2)

            tayl = Id + Bn + Bn2 / 2.d0 + Bn3 / 6.d0 + Bn4 / 24.d0

            tayl_prod = tayl
            do i = 1, n
               tayl_prod = matmul(tayl_prod, tayl_prod)
            end do

            C_tdt = matmul(tayl_prod, C_t)

            ! work out the difference to the prevois step
            error = sqrt(sum(abs(C_tdt - C_tdt_prev)**2))
            if (error < 1.d-10) exit

            C_tdt_prev = C_tdt
         end do

      end function Prop_H

   end subroutine FMS_PropQCVV

!>
!!    Propagates classical parameters, including phase
!!    @ingroup propagation
!<
   subroutine FMS_PropClassNew(T1, TimeStep)
      type(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: TimeStep
      real(kind=DefReal), allocatable, save :: MomDotSave(:)
      integer(kind=DefInt) :: iparticle, idim, jDim

      if (.not. allocated(MomDotSave)) allocate (MomDotSave(T1%NumDimensions))
      MomDotSave = FMS_MomDot(T1)
      jDim = 1

      T1%Phase = T1%Phase + TimeStep * FMS_PhaseDot(T1)
      do IParticle = 1, T1%NumParticles
         do IDim = 1, T1%Particle(IParticle)%NumDimensions

            ! TODO: Simplify this expression using T_i => T1%Particle(IParticle)
            call T1%set_pos(IParticle, IDim, &
                            T1%get_pos(IParticle, IDim) + &
                            TimeStep * FMS_PosDot(T1%Particle(IParticle), IDim) + &
                            dp5 * TimeStep * TimeStep * MomDotSave(jDim) / T1%Particle(IParticle)%Mass &
                            )

            call T1%set_mom(IParticle, IDim, &
                            T1%get_mom(IParticle, IDim) + TimeStep * MomDotSave(jDim) &
                            )

            jDim = jDim + 1
         end do
      end do

      call T1%rescale_phases()
   end subroutine FMS_PropClassNew
!>
!!    Propagates quantum parameters using a direct exponential propagator.
!!    @ingroup propagation
!<
   subroutine FMS_PropQuantum(Bundle, TimeStep)
      use GlobalModule, only: C1I
      use BundleCalcsModule, only: FMS_bH, FMS_bSp5i
      use EispackModule, only: FMS_CH
      type(T_TrajectoryBundle), intent(inout) :: Bundle
      real(kind=DefReal), intent(in) :: TimeStep

      integer(kind=DefInt) :: ITraj, JTraj, MatZ
      complex(kind=DefComp), allocatable :: CTmpMat(:, :)
      complex(kind=DefComp), allocatable :: CTmp2(:, :)
      complex(kind=DefComp), allocatable :: CVec(:, :), CValues(:)
      real(kind=DefReal), allocatable :: OReal(:, :), OImag(:, :)
      real(kind=DefReal), allocatable :: EVReal(:, :), EVImag(:, :)
      real(kind=DefReal), allocatable :: EValues(:)
      complex(kind=DefComp), allocatable :: AmpNew(:)
      real(kind=DefReal), allocatable :: FV1(:), FV2(:), FM1(:, :)
      complex(kind=DefComp), allocatable :: H(:, :), sp5i(:, :)
      integer(kind=DefInt) :: LastDimension, ierr
      save LastDimension
      save CTmp2, CTmpMat, CVec, CValues, OReal, OImag
      save EVReal, EVImag, EValues, FV1, FV2, FM1, AmpNew

      if (LastDimension /= Bundle%NumTraj) then
         if (LastDimension /= 0) then
            deallocate (CTmp2)
            deallocate (CTmpMat)
            deallocate (CVec)
            deallocate (CValues)
            deallocate (OReal)
            deallocate (OImag)
            deallocate (EVReal)
            deallocate (EVImag)
            deallocate (EValues)
            deallocate (FV1)
            deallocate (FV2)
            deallocate (FM1)
            deallocate (AmpNew)
            deallocate (H)
            deallocate (sp5i)
         end if
         LastDimension = Bundle%NumTraj
         allocate (CTmp2(Bundle%NumTraj, Bundle%NumTraj))
         allocate (CTmpMat(Bundle%NumTraj, Bundle%NumTraj))
         allocate (CVec(Bundle%NumTraj, Bundle%NumTraj))
         allocate (CValues(Bundle%NumTraj))
         allocate (OReal(Bundle%NumTraj, Bundle%NumTraj))
         allocate (OImag(Bundle%NumTraj, Bundle%NumTraj))
         allocate (EVReal(Bundle%NumTraj, Bundle%NumTraj))
         allocate (EVImag(Bundle%NumTraj, Bundle%NumTraj))
         allocate (EValues(Bundle%NumTraj))
         allocate (FV1(Bundle%NumTraj))
         allocate (FV2(Bundle%NumTraj))
         allocate (FM1(2, Bundle%NumTraj))
         allocate (AmpNew(Bundle%NumTraj))
         allocate (H(Bundle%NumTraj, Bundle%NumTraj))
         allocate (sp5i(Bundle%NumTraj, Bundle%NumTraj))
      end if

!     Build and Diagonalize Htwiddle = S^-1/2 H S^-1/2
      EVReal = 0.0d0
      EVImag = 0.0d0
      Evalues = 0.0d0
      H = FMS_bH(Bundle)
      sp5i = FMS_bSp5i(Bundle)
      CVec = matmul(H, Sp5i)
      CTmpMat = matmul(Sp5i, CVec)
      OReal = real(CTmpMat)
      OImag = real(-c1i * CTmpMat)
      MatZ = 1 ! We want BOTH Evals and Evecs

      call fms_ch(Bundle%NumTraj, Bundle%NumTraj, OReal, OImag, EValues, &
                  MatZ, EVReal, EVImag, FV1, FV2, FM1, IErr)

!     Build exp(-i Htwiddle dt)
      EValues = TimeStep * EValues
      CValues = cos(EValues) - c1i * sin(EValues)
      CVec = EVReal + c1i * EVImag

      do ITraj = 1, Bundle%NumTraj
         do JTraj = 1, Bundle%NumTraj
            CTmpMat(ITraj, JTraj) = CValues(ITraj) * conjg(CVec(JTraj, ITraj))
         end do
      end do
      CTmpMat = matmul(CVec, CTmpMat)

!     Amplitudes(t+dt)=exp(-i Htwiddle dt)*Amplitudes(t)
      do ITraj = 1, Bundle%NumTraj
         AmpNew(ITraj) = 0.0d0
         do JTraj = 1, Bundle%NumTraj
            AmpNew(ITraj) = AmpNew(ITraj) + CTmpMat(ITraj, JTraj) * Bundle%Trajectory(JTraj)%Amplitude
         end do
      end do

      do ITraj = 1, Bundle%NumTraj
         Bundle%Trajectory(ITraj)%Amplitude = AmpNew(ITraj)
         if (AmpNew(Itraj) /= AmpNew(Itraj)) then
            call FMS_DieError('Died in PropQuantum')
         end if
      end do

   end subroutine FMS_PropQuantum
!>
!!    Renormalizes a trajectory bundle such that Branching(Bundle)=OldNorms.
!!
!!    This is necessary because frozen Gaussian propagation does not
!!    conserve wavefunction normalization.
!!    \param OldNorms Norms of each trajectory to normalize to
!!    \param DNorm0 Overall normalization at time t=0; this is usually
!!    1 except when pruning trajectories on restart.
!!    @ingroup propagation
!<
   subroutine FMS_Renormalize(B1, OldNorms, DNorm0)
      use BundleCalcsModule, only: FMS_Norm, FMS_Branching
      type(T_TrajectoryBundle), intent(inout) :: B1
      real(kind=DefReal), intent(in) :: OldNorms(B1%NumStates)
      real(kind=DefReal), intent(in) :: DNorm0

      real(kind=DefReal), allocatable, save :: NewNorms(:)
      integer(kind=DefInt), save :: LastDimension
      real(kind=DefReal) :: DNorm
      integer(kind=DefInt) :: IState, ITraj

      if (B1%NumStates /= LastDimension) then
         if (LastDimension /= 0) deallocate (NewNorms)
         allocate (NewNorms(B1%NumStates))
         LastDimension = B1%NumStates
      end if

      call FMS_Branching(B1, NewNorms)
      do IState = 1, B1%NumStates
         if (abs(NewNorms(IState)) > FPZero) then
            NewNorms(IState) = sqrt(OldNorms(IState) / NewNorms(IState))
         else
            NewNorms(IState) = 0.d0
         end if
      end do

      do ITraj = 1, B1%NumTraj
         B1%Trajectory(ITraj)%Amplitude = NewNorms(B1%Trajectory(ITraj)%StateID) * B1%Trajectory(ITraj)%Amplitude
      end do

!     Now make sure overall normalization is unity...
      DNorm = FMS_Norm(B1)
      DNorm = sqrt(DNorm0 / DNorm)

      do ITraj = 1, B1%NumTraj
         B1%Trajectory(ITraj)%Amplitude = DNorm * B1%Trajectory(ITraj)%Amplitude
      end do

   end subroutine FMS_Renormalize

!>
!!    Monitors energy and norm conservation.
!!
!!    If any quantities differ from the previous timestep by the preset
!!    threshold, set flag to reject the timestep
!!
!!    Also, monitor energy or norm over the entire simulation.  If they
!!    have increased or decreased by the set thresholds over the entire
!!    simulation, shut down.
!!    \param Bundle   Trajectory bundle
!!    \param BSave    Bundle at previous timestep
!!    \param TimeStep Current time step. Used for checking only, not for
!!                    dynamics within this subroutine.
!!    @ingroup propagation
!<
   subroutine FMS_Monitor(Bundle, BSave, Timestep)
      use GlobalModule
      use TrajectoryCalcsModule, only: Potential, Kinetic
      use BundleModule
      use BundleCalcsModule, only: FMS_Norm, FMS_Branching, &
                                   Potential, FMS_PotentialB, Kinetic, FMS_KineticB
      use FMSModule, only: FMS_Shutdown
      use SMDModule, only: smSMD, SMD_Completed
      type(T_TrajectoryBundle), intent(inout) :: Bundle, BSave
      real(kind=DefReal), intent(in) :: TimeStep

      integer(kind=DefInt) :: iTraj, jTraj
      real(kind=DefReal) :: NewTimeStep, &
                            Norm, LastNorm, &
                            Energy, LastEnergy
      logical :: zDie

      logical, save :: FirstTime = .true.
      real(kind=DefReal), save :: FirstNorm, FirstEnergy
      real(kind=DefReal) :: Prob(Bundle%NumStates)

!     Output specifiers
1999  format('========================================================')
2000  format('Norm jumped by ', f7.4, '%')
2001  format('Energy jumped by ', d12.3, ' Hartrees for trajectory #', i0)
2002  format('SIMULATION ERROR - norm has changed by ', f7.4, '%')
2003  format('SIMULATION ERROR - energy has changed by ', d15.6, ' Hartree')

!     If the step is already marked for rejection, skip the work below
      if (FMS_StepRejected()) return

!     Always calculate Norm first for checking fatal errors.
      Norm = FMS_Norm(Bundle)

!     If large-timescale dynamics just jumped into the small-timescale
!     dynamics, redo this step with the small timestep
      call FMS_SetTimeStep(Bundle, NewTimeStep)
      if (TimeStep == gldTimeStep .and. &
          NewTimeStep == gldCoupTimeStep .and. &
          gldTimeStep /= gldCoupTimeStep) then
         write (fmiOut, *) 'Trajectory jumped into coupling region.'
         call FMS_RejectStep(.true.)
         if (FMS_StepRejected()) return
      end if

!
!     Check quantities relative to previous timestep
!
!     Check bundle's norm
      LastNorm = FMS_Norm(BSave)

      if (abs(LastNorm - Norm) / LastNorm > gldNormStepCons) then
         write (fmiOut, 2000) 100.0d0 * (LastNorm - Norm) / Norm
         flush (fmiOut)
         call FMS_RejectStep(.true.)
         if (FMS_StepRejected()) return
      end if

!     Check classical energy of each trajectory against the previous timestep
      jTraj = 0
      if (.not. smSMD) then
         do iTraj = 1, BSave%NumTraj
            jTraj = jTraj + 1

            if (BSave%Trajectory(jTraj)%TrajID < Bundle%Trajectory(iTraj)%TrajID) then
               jTraj = jTraj - 1
               cycle
            end if

            if (BSave%Trajectory(iTraj)%StateID == glIgnoreState) cycle

            LastEnergy = Potential(BSave%Trajectory(iTraj)) + Kinetic(BSave%Trajectory(iTraj))
            Energy = Potential(Bundle%Trajectory(jTraj)) + Kinetic(Bundle%Trajectory(jTraj))
            if (abs(LastEnergy - Energy) > gldEnergyStepCons) then
               write (fmiOut, 2001) Energy - LastEnergy, Bundle%Trajectory(jTraj)%TrajID
               flush (fmiOut)
               call FMS_RejectStep(.true.)
               if (FMS_StepRejected()) return
            end if
         end do
      end if

      if (FMS_StepRejected()) return !If step still rejected, leave

!     NOTE: The previous 'if(FMS_StepRejected()) return's were for efficiency reasons.
!     This if(FMS_StepRejected()) return, however, is absolutely necessary to prevent
!     FirstTime logic from getting messed up if the first timestep happened to be rejected

!     Beyond here, step was allowed, so check for fatal conditions

!     The first time through, get energy and norm for time 0
      if (FirstTime) then
         FirstEnergy = FMS_KineticB(Bundle, Incoherent=.true.) &
                       + FMS_PotentialB(Bundle, Incoherent=.true.)
         call FMS_Branching(Bundle, Prob, LOrth=.false.)
         FirstEnergy = FirstEnergy / sum(Prob)

         FirstNorm = FMS_Norm(Bundle)

         FirstTime = .false.
      end if

!
!     Check for SMD completion
!
      if (smSMD) then
         if (Bundle%NumTraj > 1) then
            call FMS_DieError('SMD should only have 1 trajectory')
         end if
         if (SMD_Completed()) then
            call FMS_Shutdown(B1=Bundle)
         end if
      end if

!
!     Check norm and classical energy conservation since
!     the first timestep
!
      zDie = .false.

      if (abs((FirstNorm - Norm) / FirstNorm) > gldNormCons) then
         write (fmiOut, 1999)
         write (fmiOut, 2002) 100.0d0 * (FirstNorm - Norm) / Norm

         zDie = .true.
      end if

      ! Check that the total energy of the bundle did not drift too much
      ! from the initial value at E0.
      ! This quantity is not conserved in Steered MD, XFAIMS, or
      ! if we started with more than 1 initial trajectory.
      if (.not. smSMD .and. .not. glzxfaims .and. NumInitBasis == 1) then
         Energy = Potential(Bundle, incoherent=.true.) + &
                  Kinetic(Bundle, incoherent=.true.)
         call FMS_Branching(Bundle, Prob, LOrth=.false.)
         Energy = Energy / sum(Prob)

         if (abs(FirstEnergy - Energy) > gldEnergyCons) then
            write (fmiOut, 1999)
            write (fmiOut, 2003) Energy - FirstEnergy
            ! TODO(DH): Why are we not dying here???
            ! zDie=.true.
         end if
      end if

      if (zDie) then
         call FMS_PrintMessg('FMS_Monitor: Simulation stopping.')
         call BSave%destroy()
         call FMS_ShutDown(B1=Bundle)
      end if

   end subroutine FMS_Monitor

   !>
   !!    Determines a suitable timestep for multiple timestep integration scheme
   !!    \param dt (Output) Time step to be generated
   !!    @ingroup propagation
   !<
   subroutine FMS_SetTimeStep(B1, dt)
      use GlobalModule, only: DefReal, DefInt, gldTimeStep, gldCoupTimeStep, &
                              gldTStepThresh, glIgnoreState, glzXFAIMS, glzXFActive
      use XFAIMSModule, only: xfaims_params
      use BundleCalcsModule, only: FMS_bH
      type(T_TrajectoryBundle), intent(in) :: B1
      real(kind=DefReal), intent(out) :: dt

      integer(kind=DefInt) :: iTraj, jTraj, iState
      logical :: zCoupling

      zCoupling = .false.

!     Check for large non-adiabatic coupling matrix elements
!     This also checks for overlap between trajectories on the same state
      do iTraj = 1, B1%NumTraj - 1
         if (zCoupling) exit
         do jTraj = iTraj + 1, B1%NumTraj
            if (abs(FMS_bH(B1, iTraj, jTraj)) > gldTStepThresh) then
               zCoupling = .true.
               exit
            end if
         end do
      end do

      ! Check if any of the trajectories have large nonadiabatic coupling vector
      do iTraj = 1, B1%NumTraj
         if (zCoupling) exit
         do iState = 1, B1%Trajectory(iTraj)%NumStates
            if (iState == B1%Trajectory(iTraj)%StateID) cycle
            ! TODO(DH): Switch to SpawnModule::spawn_couple
            ! https://github.com/ispg-group/fms90-redux/issues/79
            if (abs(FMS_SpawnDCouple(B1%Trajectory(iTraj), iState)) > spawn_params%CSThresh) then
               zCoupling = .true.
               exit
            end if
         end do
      end do

      ! Set timestep accordingly
      if (zCoupling) then
         dt = gldCoupTimeStep
      else
         dt = gldTimeStep
      end if

      ! xfaims
      if (glzxfaims) then
         if (glzxfactive) then
            dt = xfaims_params%CoupFieldTimeStep
         end if
         if (B1%CurrentTime > xfaims_params%sp_spwn_f) then
            if (xfaims_params%IgnoreStateAferField) then
               glIgnoreState = 1
            end if
         end if
      end if

   end subroutine FMS_SetTimeStep

end module PropagationModule
