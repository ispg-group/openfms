!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Main driver for sampling initial conditions.
!!    @ingroup initial
!!    \callgraph
!<
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine FMS_InitialCond(B1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use GlobalModule, only: DefInt, DefReal, fmiOut, FMS_RejectStep, &
          gldTimeStep, glzCentroids, glzConstrain, glzTimeStepRejection
      use InitialModule, only: inzInitBright, inzInitDark, &
         inzSharpEnergy, inzOSAmp, inzNormInitial, inzMirrorBasis, &
         FMS_CalcBrightDark, FMS_SpawnMirrorMP, FMS_SelectInitial
      use TrajectoryCalcsModule, only: FMS_TransDipole, Potential
      use BundleModule
      use BundleCalcsModule, only: FMS_Norm, FMS_UpdateCentroid
      use SpawnModule, only: FMS_AdjustEnergy
      use RestartModule, only: iniRestart
      use VerletModule, only: FMS_PropVV
      implicit none
      type(T_TrajectoryBundle),intent(inout) :: B1

      integer (kind=DefInt) :: ITraj
      real (kind=DefReal) :: dNorm0
      logical :: MErr
      logical :: RejectionSave
      integer (kind=DefInt) :: i,j
      integer (kind=DefInt) :: IState
      real (kind=DefReal) :: EGap,OS,TDip(3)
      integer (kind=DefInt) :: iCBF,jCBF

!     Suppress timestep rejection during system preparation
      RejectionSave=glzTimeStepRejection
      glzTimeStepRejection=.false.

      write(fmiOut,'(a)') 'Generating initial conditions:'

!     Sample initial conditions

      call FMS_SelectInitial( B1 )
      do iTraj= 1, B1%NumTraj
!         B1%Trajectory(iTraj)%TrajID = itraj

         if( glzConstrain )then
            call FMS_PropVV( B1%Trajectory(iTraj),+gldTimeStep )
            call FMS_PropVV( B1%Trajectory(iTraj),-gldTimeStep )
         endif

      enddo

!     FMS_CalcBright will determine the initial state based on strength of
!     transition dipole
      if(inzInitBright.or.inzInitDark) Then
         do ITraj=1,B1%NumTraj
            call FMS_CalcBrightDark(B1%Trajectory(ITraj))
         enddo
      endif

!     Generate a mirror basis if needed
      if(inzMirrorBasis) then
         call FMS_SpawnMirrorMP(B1)
      endif

!     Adjust momenta of pre-spawned virtuals if desired
      if (inzSharpEnergy) then
         do ITraj=2,B1%NumTraj
            MErr=FMS_AdjustEnergy(B1%Trajectory(ITraj),B1%Trajectory(1))
         enddo
      endif

      ! Set amplitude of initial trajectories based on oscillator strength?
      ! TODO: Separate this into a function
      if (inzOSAmp) then

        write(fmiOut,*) 'Scaling trajectory amplitudes based on oscillator strength'

!       This is currently only correct for a single trajectory
!       (Independent First Generation)
        if (B1%NumTraj>1) then
           write(fmiOut, *) 'ERROR: setting amplitude according to oscillator strength '
           write(fmiOut, *) 'is correct only for a single trajectory'
           stop 1
        endif

        if (inzNormInitial) then
          write(fmiOut,*)'ERROR: both inzOSAmp and inzNormInitial set'
          stop 1
        endif

        do ITraj=1,B1%NumTraj
          IState=B1%Trajectory(ITraj)%StateID
          if (IState==1) then
            write(fmiOut,'(a,i0)')'Warning: trying to set initial amplitude of trajectory ', iTraj
            write(fmiOut, *) 'based on oscillator strength, but it is in ground state, so no scaling performed'
            cycle
          end if
          Tdip=FMS_TransDipole(B1%Trajectory(ITraj),IState)
          Egap = Potential(B1%Trajectory(ITraj)) - Potential(B1%Trajectory(ITraj), 1)
          OS = (2.0d0 / 3.0d0) * Egap * dot_product(Tdip, Tdip)
          B1%Trajectory(ITraj)%Amplitude = B1%Trajectory(ITraj)%Amplitude * dcmplx(sqrt(OS), 0.0d0)
        enddo
      endif
!     check energy and normalize amplitude
      DNorm0=FMS_Norm(B1)

!     This makes no sense at all
      if (inzNormInitial .and. inIRestart == 0) then
!     if (inzNormInitial )then
         do ITraj=1,B1%NumTraj
            B1%Trajectory(ITraj)%Amplitude = B1%Trajectory(ITraj)%Amplitude / sqrt(DNorm0)
         enddo
      endif

!     Update centroids
!bfec
      if (glzCentroids) then
         do i =2, B1%NumTraj
            do j = 1, i-1
               if( B1%Trajectory(i)%Ms == 2 .and. B1%Trajectory(j)%Ms == 2) then
                  iCBF=B1%Trajectory(i)%CBF
                  jCBF=B1%Trajectory(j)%CBF

                  call FMS_UpdateCentroid(B1%Trajectory(i), B1%Trajectory(j), &
                     B1%Centroids(((iCBF-2)*(iCBF-1))/2+jCBF))
           endif
         enddo
       enddo
      endif

      write(fmiOut,*)
      flush(fmiOut)

!     Restore timestep rejection state
      glzTimeStepRejection=RejectionSave
      call FMS_RejectStep(.false.)

      end subroutine FMS_InitialCond
