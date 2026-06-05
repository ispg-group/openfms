!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    @brief Procedured for calculations on the Bundle
!<
module BundleCalcsModule
   !--------------------------------------------------------------------------------
   use GlobalModule
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_Weight, FMS_WeightC, FMS_KineticT, FMS_IsBundleCurrent, &
                                    FMS_IsAmpDotCurrent, FMS_AmpDotUpdated, FMS_BundleUpdated
   use BundleModule
   use EispackModule, only: FMS_CH
   use OverlapModule, only: overlap, overlap_KE, overlap_V, overlap_S_Dot

   implicit none
   private

   public :: FMS_Norm, FMS_Branching, FMS_UpdateCentroid
   public :: FMS_Set_Amplitude, FMS_Get_Amplitude
   public :: FMS_UpdateBundle, FMS_BuildHS
   public :: FMS_bSDotComp, FMS_bSDot
   public :: FMS_bH, FMS_bHeff, FMS_bSp5i, FMS_bSInvMat
   public :: Kinetic, Potential, FMS_KineticB, FMS_PotentialB
   public :: FMS_Mulliken, FMS_UpdateMulliken

!     Overloading functions:

   interface FMS_bSp5i
      module procedure FMS_bSp5iComp
      module procedure FMS_bSp5iMat
   end interface

   interface FMS_bH
      module procedure FMS_bHComp
      module procedure FMS_bHMat
   end interface

   interface FMS_bS
      module procedure FMS_bSComp
      module procedure FMS_bSMat
   end interface

   interface FMS_bSInv
      module procedure FMS_bSInvComp
      module procedure FMS_bSInvMat
   end interface

   interface FMS_bSDot
      module procedure FMS_bSDotComp
      module procedure FMS_bSDotMat
   end interface

   interface FMS_bHeff
      module procedure FMS_bHeffComp
      module procedure FMS_bHeffMat
   end interface

   interface FMS_bHeff1
      module procedure FMS_bHeff1Comp
      module procedure FMS_bheff1Mat
   end interface

   interface Kinetic
      module procedure FMS_KineticB
   end interface

   interface Potential
      module procedure FMS_PotentialB
   end interface

contains
!--------------- Module Procedures ---------------!

! Procedures to return bundle quantities
!>
!!    Recalculate bundle matrices if necessary.
!<
   subroutine FMS_UpdateBundle(Bundle)
      type(T_TrajectoryBundle) :: Bundle
      integer(kind=DefInt) :: ITraj

!     If we will reject the current timestep, don't bother
!     doing any calculations
      if (FMS_StepRejected()) return

      do iTraj = 1, Bundle%NumTraj
         if (.not. FMS_IsBundleCurrent(Bundle%Trajectory(iTraj))) then
            call FMS_BuildHS(Bundle)
            exit
         end if
      end do

   end subroutine FMS_UpdateBundle

   function FMS_Get_Amplitude(Bundle) result(amps)
      type(T_TrajectoryBundle), intent(in) :: Bundle
      complex(kind=DefComp) :: amps(Bundle%NumTraj)

      integer :: i

      do i = 1, Bundle%NumTraj
         amps(i) = Bundle%Trajectory(i)%Amplitude
      end do

   end function FMS_Get_Amplitude

   subroutine FMS_Set_Amplitude(Bundle, amps)
      type(T_TrajectoryBundle) :: Bundle
      complex(kind=DefComp), intent(in) :: amps(Bundle%NumTraj)
      integer :: i

      do i = 1, Bundle%NumTraj
         Bundle%Trajectory(i)%Amplitude = amps(i)
      end do

   end subroutine FMS_Set_Amplitude

!>
!!    Return AmpDot for trajectory iTraj.
!!
!!    Scope: Public
!<
   function FMS_AmpDot(Bundle, ITraj) result(AmpDot)
      type(T_TrajectoryBundle) :: Bundle
      integer(kind=DefInt) :: ITraj, i
      complex(kind=DefComp) :: AmpDot

      ! See if an update is needed
      do i = 1, Bundle%NumTraj
         if (.not. FMS_IsAmpDotCurrent(Bundle%Trajectory(iTraj))) then
            call FMS_CalcAmpDot(Bundle)
            exit
         end if
      end do

      AmpDot = Bundle%BMatrices%Ampdot(ITraj)

   end function FMS_AmpDot
!>
!!    Computes centroid between T1 and T2 and places this information in Centroid
!!    Unlike other functions in this module, this operates on
!!    Trajectory types, not Bundles.
!!    Make this into a Bundle bound procedure, only
!!    passing in the T1 and T2 indexes, Centroid index should be computed here!
!<
   subroutine FMS_UpdateCentroid(T1, T2, Centroid)
      use OverlapModule, only: overlap_dx_M_trajectory
      type(T_Trajectory), intent(in) :: T1, T2
      type(T_Trajectory), intent(inout) :: Centroid
      integer(kind=DefInt) :: iPart
      integer(kind=DefInt) :: iTrajSave, iCentSave(2)
      complex(kind=DefComp) :: v(T1%NumDimensions)
      real(kind=DefReal), dimension(T1%NumDimensions) :: v_re, v_imag

!     Check that this is the correct centroid
      if (.not. ((Centroid%CentID(1) == T1%CBF .and. Centroid%CentID(2) == T2%CBF) .or. &
                 (Centroid%CentID(2) == T1%CBF .and. Centroid%CentID(1) == T2%CBF))) then
         write (fmiOut, *) 'CBF 1: ', T1%CBF
         write (fmiOut, *) 'CBF 2: ', T2%CBF
         write (fmiOut, *) 'Centroid: ', Centroid%TrajID
         write (fmiOut, *) Centroid%CentID
         call FMS_DieError('ERROR in UpdateCentroid: incorrect centroid passed.')
      end if

!     If there's no electronic structure for the centroid, just set
!     its wavefunction to that of a parent, and tell electronic
!     structure package to ignore phasing or diabatization problems
      if (.not. Centroid%ESFlags%zESExists) then
         iTrajSave = Centroid%TrajID
         iCentSave = Centroid%CentID

         Centroid = T1

         Centroid%zCent = .true.
         Centroid%TrajID = iTrajSave
         Centroid%CentID = iCentSave

         call Centroid%geom_changed()
         ! The following will prevent a rejected step on the centroid's first step
         ! which could othewise happen if two of its states are flipped relative to
         ! parent - i.e. its on oposite side of a CI.
         ! The code should  make sure phase is continuous nevertheless.
         Centroid%ESFlags%zIgnoreErrors = .true.
      end if

!     Reset centroid time
      call Centroid%set_time(T1%get_time())

!     Set centroid state to the first trajectory's
      Centroid%StateID = T1%StateID

!
!     Calculate coordinates of centroid
!
      do iPart = 1, T1%NumParticles
         call Centroid%set_pos(iPart, (T1%Particle(iPart)%Width * T1%Particle(iPart)%get_pos() + &
                                       T2%Particle(iPart)%Width * T2%Particle(iPart)%get_pos()) / &
                               (T1%Particle(iPart)%Width + T2%Particle(iPart)%Width))
      end do

!     Set velocity to be off-diagonal components of (p/m) between parents
      v = overlap_dx_M_trajectory(T1, T2)
      v_re = real(v)
      v_imag = aimag(v)
      call Centroid%set_vel(v_re)
      call Centroid%set_vel2(v_imag)

   end subroutine FMS_UpdateCentroid

!---------------------------------------------------------
!     Routines to get a component of the matrices
!---------------------------------------------------------

!>
!!    Returns matrix element I, J of Hamiltonian.
!!
!!    Scope: Public
!<
   function FMS_BHComp(bundle, iTraj, jTraj) result(xH)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xH

      call FMS_UpdateBundle(bundle)

      xH = Bundle%BMatrices%H(ITraj, JTraj)

   end function FMS_BHComp

!>
!!    Returns matrix element I, J of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSComp(bundle, iTraj, jTraj) result(xS)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xS

      call FMS_UpdateBundle(bundle)

      xS = Bundle%BMatrices%S(ITraj, JTraj)

   end function FMS_bSComp

!>
!!    Returns matrix element I, J of the inverse of the overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSInvComp(bundle, iTraj, jTraj) result(xSInv)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xSInv

      call FMS_UpdateBundle(bundle)

      xSInv = Bundle%BMatrices%SInv(ITraj, JTraj)

   end function FMS_bSInvComp

!>
!!    Returns matrix element I, J of time derivative of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSDotComp(bundle, iTraj, jTraj) result(xSDot)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xSDot

      call FMS_UpdateBundle(bundle)

      xSDot = Bundle%BMatrices%SDot(ITraj, JTraj)

   end function FMS_bSDotComp

!>
!|    Returns matrix element I, J of effective hamiltonian
!|
!|    Scope: Public
!<
   function FMS_bHeffComp(bundle, iTraj, jTraj) result(xHeff)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xHeff

      call FMS_UpdateBundle(bundle)

      xHeff = Bundle%BMatrices%Heff(ITraj, JTraj)

   end function FMS_bHeffComp

!>
!!    Returns matrix element I, J of square root of inverse of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSp5iComp(bundle, iTraj, jTraj) result(xSp5i)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xSp5i

      call FMS_UpdateBundle(bundle)

      xSp5i = Bundle%BMatrices%Sp5i(ITraj, JTraj)

   end function FMS_bSp5iComp

!>
!!    Returns matrix element I, J of square root of inverse of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSp5Comp(bundle, iTraj, jTraj) result(xSp5)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xSp5

      call FMS_UpdateBundle(bundle)

      xSp5 = Bundle%BMatrices%Sp5(ITraj, JTraj)

   end function FMS_bSp5Comp

!>
!!    Returns matrix element I, J of the effective hamiltonian w/o S^-1
!!
!!    Scope: Public
!<
   function FMS_bHeff1Comp(bundle, iTraj, jTraj) result(xHeff1)
      type(T_TrajectoryBundle) :: bundle
      integer(kind=DefInt) :: iTraj, jTraj
      complex(kind=DefComp) :: xHeff1

      call FMS_UpdateBundle(bundle)

      xHeff1 = Bundle%BMatrices%Heff1(ITraj, JTraj)

   end function FMS_bHeff1Comp

!---------------------------------------------------------
!     Full matrix versions of the same functions
!---------------------------------------------------------
!>
!!    Returns Hamiltonian matrix
!!
!!    Scope: Public
!<
   function FMS_BHMat(bundle) result(xH)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xH(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xH = Bundle%BMatrices%H

   end function FMS_BHmat

!>
!!    Returns overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSMat(bundle) result(xS)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xS(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xS = Bundle%BMatrices%S

   end function FMS_bSMat

!>
!!    Returns inverse of the overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSInvMat(bundle) result(xSInv)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xSInv(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xSInv = Bundle%BMatrices%SInv

   end function FMS_bSInvMat

!>
!!    Returns time derivative of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSDotMat(bundle) result(xSDot)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xSDot(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xSDot = Bundle%BMatrices%SDot

   end function FMS_bSDotMat

!>
!!    Returns effective Hamiltonian matrix
!!
!!    Scope: Public
!<
   function FMS_bHeffMat(bundle) result(xHeff)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xHeff(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xHeff = Bundle%BMatrices%Heff

   end function FMS_bHeffMat

!>
!!    Returns square root of inverse of overlap matrix
!!
!!    Scope: Public
!<
   function FMS_bSp5iMat(bundle) result(xSp5i)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xSp5i(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xSp5i = Bundle%BMatrices%Sp5i

   end function FMS_bSp5iMat

!>
!!    Returns effective Hamiltonian w/o S^-1
!!
!!    Scope: Public
!<
   function FMS_bHeff1Mat(bundle) result(xHeff1)
      type(T_TrajectoryBundle) :: bundle

      complex(kind=DefComp) :: xHeff1(bundle%NumTraj, Bundle%NumTraj)

      call FMS_UpdateBundle(bundle)

      xHeff1 = Bundle%BMatrices%Heff1

   end function FMS_bHeff1Mat

!>
!!    Calculates norm of trajectory bundle as sum of branching ratios.
!!    @ingroup analysis
!<
   function FMS_Norm(B1) result(DNorm)
      type(T_TrajectoryBundle) :: B1
      real(kind=DefReal) :: DNorm
      real(kind=DefReal), allocatable :: Ratios(:)
      save Ratios

      if (.not. allocated(Ratios)) then
         allocate (Ratios(B1%NumStates))
      end if

      call FMS_Branching(B1, Ratios)
      DNorm = sum(Ratios)
   end function FMS_Norm

! Driver for calculating bundle quantities
!>
!!    Build Hamiltonian and overlap matrices for a bundle of trajectories.
!!
!!    SSOnly is true if we compute only intra-surface (diagonal in
!!       electronic state index) terms
!!    SSOnly is false if only inter-surface terms are to be computed
!!
!!    SSOnly is false for split propagator, where quantum propagation is
!!       only for nonadiabatic coupling of basis functions
!!    SSOnly should be true for fully coupled propagation -- This means
!!       that the variable is misnamed.  SSOnly true means both single
!!       surface and multi-surface are computed.  SSOnly false means
!!       only multi-surface terms are computed.
!<
   subroutine FMS_BuildHS(Bundle)
      type(T_TrajectoryBundle), target :: Bundle

      real(kind=DefReal) :: Threshold
      integer(kind=DefInt) :: IState, JState
      integer(kind=DefIntBlas) :: MatZ, Jerr, I4NumTraj
      integer(kind=DefInt) :: ISaddle

      real(kind=DefReal), dimension(Bundle%NumTraj) :: FV1, FV2, EValues

      real(kind=DefReal) :: DSTempSp5i, DSTemp

      real(kind=DefReal), dimension(Bundle%NumTraj, Bundle%NumTraj) :: EVReal, EVImag, OReal, OImag

      real(kind=DefReal), dimension(2, Bundle%NumTraj) :: FM1

      complex(kind=DefComp), dimension(Bundle%NumTraj, Bundle%NumTraj) :: CTmpMat, CEvec, CTmpMatSp5i

      integer(kind=DefInt) :: NTossed, ICent
      logical :: IterInv
      common / junky / NTossed
      common / CGCommon / IterInv

      integer(kind=DefInt) :: ntraj, i, j, CBFi, CBFj

      type(T_Trajectory), pointer :: T_i, T_j, T_ij ! these point to trajectories
      ! and centroids from the bundle
      ! makes the code readable!
      complex(kind=DefComp), pointer, dimension(:, :) :: S, H, S_dot

      real(kind=DefReal) :: time_tmp1, time_tmp2

      call cpu_time(time_tmp1)

      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "We are building HS, check NumTraj, NCBFs"
      !write (fmiout, *) "Number trajs:", Bundle%NumTraj
      !write (fmiout, *) "Number CBFs:", Bundle%NCBFs
      !write (fmiout, *) " ----------------------------------------------------"

      Threshold = gldRegThresh

      ntraj = Bundle%NumTraj

      S => Bundle%BMatrices%S
      H => Bundle%BMatrices%H
      S_dot => Bundle%BMatrices%SDot

      S = (0., 0.)
      H = (0., 0.)
      S_dot = (0., 0.)

      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Checking Centroid indices before calculating overlap: "
      !write (fmiout, *) " ----------------------------------------------------"
      !do i = 1, (((Bundle%NCBFs - 1) * Bundle%NCBFs) / 2)
      !   write (fmiout, *) "CentID", Bundle%Centroids(i)%TrajID
      !   write (fmiout, *) "Position", Bundle%Centroids(i)%Particle(1)%get_pos()
      !end do
      !write (fmiout, *) " ----------------------------------------------------"

!     Load H,S, and SDot  matrices.
      do i = 1, ntraj
         do j = 1, i

            T_i => Bundle%Trajectory(i)
            T_j => Bundle%Trajectory(j)
            T_ij => null()
!bfec
!        GAIMS changed
!        The idea is that within one CBF you once update T_ij for
!        the Ms=2 TBF and then use the same centroid for every
!        TBF within the CBF. When the next Ms=2 TBF appears, we
!        know that a new CBF is reached.
            if (glzCentroids .and. T_i%CBF /= T_j%CBF) then
               !write(fmiOut,*) "i, j: ", i,j
               !write(fmiOut,*) "IDi, IDj: ",T_i%StateID,T_j%StateID
               !write(fmiOut,*) "Msi, Msj: ", T_i%Ms,T_j%Ms
               !write(fmiOut,*) "CBFi, CBFj: ", T_i%CBF,T_j%CBF

!!          if( T_i%Ms.eq.2 .and. T_j%Ms.eq.2 ) then
               CBFi = T_i%CBF
               CBFj = T_j%CBF
               ICent = ((CBFi - 2) * (CBFi - 1)) / 2 + CBFj
               !write (fmiout, *) "This is what we are trying to use as CentID in BundleCalcs: ", ICent
               T_ij => Bundle%Centroids(ICent)
               !write (fmiout, *) "Position", Bundle%Centroids(ICent)%Particle(1)%get_pos()
!!          endif
            end if
!        GAIMS changed end

            IState = T_i%StateID
            JState = T_j%StateID

            ISaddle = gliSaddle
            if (gliSaddle == 2) ISaddle = 0
            if (gliSaddle == 3) then
               if (i == j) then
                  ISaddle = 0
               else
                  ISaddle = 1
               end if
            end if

            S(i, j) = overlap(T_i, T_j)
            !write (fmiout, *) "I get here!!!"
            !write (fmiout, *) "Are there Centroids?", glzCentroids

! Overlap will return elements of the individual overlap, which
! can be passed to kinetic.  The overlaps themselves will be
! passed to Potential and to CGSDoTT instead of being recalculated.
! In case of Split-Operator: no need to calculate diagonal elements of H

            if (glzFullyCoupled .or. IState /= JState) then
               ! If the overlap is beneath threshold, set matrix element to 0
               if (abs(S(i, j)) >= FPZero) then
                  !write (fmiout, *) "----------------------------------------------------------------"
                  if (i == j) then

                     !write (fmiout, *) "Pair of twice the same trajectory: i, j", i, j
                     H(i, j) = overlap_KE(T_i, T_i, S(i, j)) + overlap_V(T_i, S_ij_precalc=S(i, j))

                  else

                     if (glzCentroids) then
                        !write (fmiout, *) "Pair of different trajectories: i, j", i, j
                        !write (fmiout, *) "Calling overlap_V with existing centroids"
                        !write (fmiout, *) "This is the Centroid: ICent", ICent
                        H(i, j) = overlap_KE(T_i, T_j, S(i, j)) + overlap_V(T_i, T_j, T_ij, S(i, j))
                        !write (fmiout, *) "Got H(i, j)"
                     else
                        H(i, j) = overlap_KE(T_i, T_j, S(i, j)) + overlap_V(T_i, T_j, S_ij_precalc=S(i, j))
                     end if

                     !write (fmiout, *) "----------------------------------------------------------------"
                  end if
               end if
            end if

            ! fill in the offdiagonal
            H(j, i) = conjg(H(i, j))
            S(j, i) = conjg(S(i, j))

            ! zero the overlaps bettwen trajectories on different states
            ! GAIMS update
            if ((IState /= JState) .or. (T_i%Ms /= T_j%Ms)) then
               S(i, j) = 0.d0
               S(j, i) = 0.d0
            end if
            ! GAIMS update end

            ! For SSOnly=false, meaning no single-surface terms should survive,
            ! make sure the appropriate H elements are zero.
            if (.not. glzFullyCoupled .and. IState == JState) then
               H(i, j) = 0.d0
               H(j, i) = 0.d0
            end if

            ! Matrix elements of time derivative
            ! If we are propagating with fixed nuclei, matrix element of
            ! time-derivative should be zero.  I am not sure this is true for
            ! non-split integrators, but it is easily seen empirically to be true
            ! when the operator splitting is employed.
            if (glzFullyCoupled) then
               if (abs(S(i, j)) > FPZero) then
                  S_dot(i, j) = overlap_S_dot(T_i, T_j, S(i, j))
                  S_dot(j, i) = overlap_S_dot(T_j, T_i, S(j, i))
               end if
            end if

         end do
      end do

      call cpu_time(time_tmp2)
      bhstime = bhstime + time_tmp2 - time_tmp1

      if (.not. IterInv) then
!---------------------------------------------------------------
!     Invert S matrix using a regularized inversion code.
!---------------------------------------------------------------
         call cpu_time(time_tmp1)
         MatZ = 1 ! We want BOTH Evals and Evecs
         OReal = real(Bundle%BMatrices%S)
         OImag = real(-c1i * Bundle%BMatrices%S)

         I4NumTraj = Bundle%NumTraj
         call fms_ch(I4NumTraj, I4NumTraj, OReal, OImag, EValues, MatZ, EVReal, EVImag, FV1, FV2, FM1, Jerr)

         CEvec = EVReal + c1i * EVImag

         NTossed = 0
         CTmpMat = 0.0d0
         do i = 1, Bundle%NumTraj
            if (abs(EValues(i)) > Threshold) then
               DSTemp = d1 / EValues(i)
               DSTempSp5i = d1 / sqrt(EValues(i))
            else
               NTossed = NTossed + 1
               DSTemp = 0.0d0
               DSTempSp5i = 0
            end if
            do j = 1, Bundle%NumTraj
               CTmpMat(i, j) = DSTemp * conjg(CEvec(j, i))
               CTmpMatSp5i(i, j) = DSTempSp5i * conjg(CEvec(j, i))
            end do
         end do

         ! Assemble the effective Hamiltonian which enters the coupled equations
         ! for the trajectory amplitudes
         Bundle%BMatrices%SInv = matmul(CEvec, CTmpMat)
         ! mbn: on the HP platform it is not possible to compile with the -O2 flag
         ! if Bundle%HEff is used instead of an intermediate matrix (CTmpMat).
         ! This should not be the case but the use of CTmpMat solves the problem
         CTmpMat = Bundle%BMatrices%H - c1i * Bundle%BMatrices%SDot
         Bundle%BMatrices%Sp5i = matmul(CEvec, CTmpMatSp5i)
         Bundle%BMatrices%HEff1 = Bundle%BMatrices%H - c1i * Bundle%BMatrices%SDot
         Bundle%BMatrices%HEff = matmul(Bundle%BMatrices%SInv, CTmpMat)

         call cpu_time(time_tmp2)
         invtime = invtime + time_tmp2 - time_tmp1
      else
         Bundle%BMatrices%HEff1 = Bundle%BMatrices%H - c1i * Bundle%BMatrices%SDot
      end if

!     Flip flags
      do i = 1, Bundle%NumTraj
         call FMS_BundleUpdated(Bundle%Trajectory(i))
      end do

      if (glzCentroids) then
         nullify (T_i, T_j, T_ij, S, H, S_dot)
      else
         nullify (T_i, T_j, S, H, S_dot)
      end if

   end subroutine FMS_BuildHS

!>
!!    Loads a trajectory bundle structure with time derivatives.
!!    We do not compute H and S here -- this is assumed already done.
!!    The result is loaded into B1%BMatrices%AmpDot
!<
   subroutine FMS_CalcAmpDot(B1)
      type(T_TrajectoryBundle) :: B1 !< Calculate time-derivatives quantum amplitude for this bundle
      complex(kind=DefComp), allocatable :: Amp(:)
      complex(kind=DefComp), allocatable :: CX(:), CB(:)
!     complex (kind=DefComp) CX(B1%NumTraj),CB(B1%NumTraj)
!     complex (kind=DefComp) CA(B1%NumTraj,B1%NumTraj)
      complex(kind=DefComp), allocatable :: SMat(:, :)
      integer(kind=DefInt) :: ITraj, JTraj, NumTraj
      integer(kind=DefInt) :: LastDimension
      real(kind=DefReal) :: time_tmp1, time_tmp2

      save LastDimension
      save Amp, CX, CB, SMat

      if (LastDimension /= B1%NumTraj) then
         if (allocated(Amp)) then
            deallocate (Amp)
            deallocate (CX)
            deallocate (CB)
            deallocate (SMat)
         end if
         LastDimension = B1%NumTraj
         allocate (Amp(B1%NumTraj))
         allocate (CX(B1%NumTraj))
         allocate (CB(B1%NumTraj))
         allocate (SMat(B1%NumTraj, B1%NumTraj))
      end if

      do ITraj = 1, B1%NumTraj
         Amp(ITraj) = B1%Trajectory(ITraj)%Amplitude
      end do

!     dAmplitude/dt = -i H' Amplitude
      if (.not. GlzIterInv) then
         do ITraj = 1, B1%NumTraj
            B1%BMatrices%AmpDot(ITraj) = 0
            do JTraj = 1, B1%NumTraj
               B1%BMatrices%AmpDot(ITraj) = B1%BMatrices%AmpDot(ITraj) - c1i * FMS_bHEff(B1, ITraj, JTraj) * Amp(JTraj)
            end do
         end do
      else
!     tjm  Alternative procedure is to build b=HEff1*Amp
!     tjm  Then solve S x ampdot = b
!     DH: Deactivate this branch to get rid of dependency on BLAS/LAPACK
         call FMS_DieError('GlzIterInv is not supported right now')

         call cpu_time(time_tmp1)

         do ITraj = 1, B1%NumTraj
            CX(ITraj) = B1%BMatrices%AmpDot(ITraj)
!     The process doesn't seem to converge at all if seeded with zero
            if (CX(ITraj) == (0.0, 0.0)) CX(ITraj) = (1.0d-10, 1.0d-10)
            CB(ITraj) = d0
         end do

         do ITraj = 1, B1%NumTraj
            do JTraj = 1, B1%NumTraj
               CB(ITraj) = CB(ITraj) + FMS_bHEff1(B1, ITraj, JTraj) * Amp(JTraj)
            end do
         end do

         do ITraj = 1, B1%NumTraj
            CB(ITraj) = -c1i * CB(ITraj)
         end do
         NumTraj = B1%NumTraj

         SMat = FMS_bS(B1)
         ! DH: Commented out, see above
         ! call FMS_augCG(SMat,CX,CB,B1%NumTraj,B1%NumTraj)

         do ITraj = 1, B1%NumTraj
            B1%BMatrices%AmpDot(ITraj) = CX(ITraj)
         end do

         call cpu_time(time_tmp2)
         cgtime = cgtime + time_tmp2 - time_tmp1

      end if

      do iTraj = 1, B1%NumTraj
         call FMS_AmpDotUpdated(B1%Trajectory(ITraj))
      end do

   end subroutine FMS_CalcAmpDot

!>
!!    Computes Mulliken-like populations \f$ P_i \f$ for a trajectory bundle.
!!
!!    \f[
!!         P_i= \sum_j \Re\left( \rho_{ij} S_{ij}\right)
!!    \f]
!!
!!    \return Population calculated populations for each trajectory
!!    @ingroup analysis
!<
   function FMS_Mulliken(B1) result(Population)
      type(T_TrajectoryBundle), intent(inout) :: B1 ! will will update S for fun
      ! and store the population
      real(kind=DefReal) :: Population(B1%NumTraj)

      integer(kind=DefInt) :: i

      call FMS_UpdateMulliken(B1)

      do i = 1, B1%NumTraj
         Population(i) = B1%Trajectory(i)%Pop
      end do
   end function FMS_Mulliken

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_UpdateMulliken(B1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! calculates the Mulliken populations and stores them in the Trajectory structure
      type(T_TrajectoryBundle), target, intent(inout) :: B1 ! will will update S for fun
      ! and store the population
      integer(kind=DefInt) :: ntraj, i, j, is, js

      complex(kind=DefComp), pointer :: S(:, :)
      complex(kind=DefComp) :: Ctemp

      real(kind=DefReal) :: pop

      ntraj = B1%NumTraj

      S => B1%BMatrices%S

      do i = 1, ntraj
         pop = 0.d0
         do j = 1, ntraj
            is = B1%Trajectory(i)%StateID
            js = B1%Trajectory(j)%StateID

            S(i, j) = overlap(B1%Trajectory(i), B1%Trajectory(j), same_state=.true.)
            CTemp = conjg(B1%Trajectory(i)%Amplitude) * S(i, j) * B1%Trajectory(j)%Amplitude
            pop = pop + real(CTemp)
         end do
         ! store the populations
         B1%Trajectory(i)%Pop = pop
      end do

   end subroutine FMS_UpdateMulliken

!>
!!    Computes the kinetic energy of a trajectory bundle
!!    \param Incoherent (Optional) If .True., compute an incoherent
!!    sum.  The default is .False., i.e. calculate a coherent sum.
!!    \return Energy Kinetic energy
!!    @ingroup analysis
!<
   function FMS_KineticB(Bundle, Incoherent) result(Energy)
      type(T_TrajectoryBundle), intent(in) :: Bundle
      real(kind=DefReal) :: Energy, EnergyInc
      logical, optional :: Incoherent
      integer(kind=DefInt) :: ITraj, JTraj
      logical :: IncSum

      IncSum = .false.
      if (present(Incoherent)) IncSum = Incoherent

      Energy = 0.0d0
      do ITraj = 1, Bundle%NumTraj
         if (.not. IncSum) then
            do JTraj = 1, ITraj - 1
               EnergyInc = 2.0 * real(conjg(Bundle%Trajectory(ITraj)%Amplitude) * Bundle%Trajectory(JTraj)%Amplitude * &
                                      overlap_KE(Bundle%Trajectory(ITraj), Bundle%Trajectory(JTraj)))
               Energy = Energy + EnergyInc
            end do
            EnergyInc = FMS_Weight(Bundle%Trajectory(ITraj)) * &
                        real(overlap_KE(Bundle%Trajectory(ITraj), Bundle%Trajectory(ITraj)))
            Energy = Energy + EnergyInc
         else
            Energy = Energy + FMS_Weight(Bundle%Trajectory(ITraj)) * real(FMS_KineticT(Bundle%Trajectory(ITraj)))
         end if
      end do

!
!     Add energy from dead trajectories
!
      do ITraj = 1, Bundle%NumDeadTraj
         if (.not. IncSum) then
            do JTraj = 1, ITraj - 1
               Energy = Energy + 2.0 * real( &
                        conjg(Bundle%DeadTraj(ITraj)%Amplitude) * Bundle%DeadTraj(JTraj)%Amplitude * &
                        overlap_KE(Bundle%DeadTraj(ITraj), Bundle%DeadTraj(JTraj)))
            end do
            Energy = Energy + FMS_Weight(Bundle%DeadTraj(ITraj)) * &
                     real(overlap_KE(Bundle%DeadTraj(ITraj), Bundle%DeadTraj(ITraj)))
         else
            Energy = Energy + FMS_Weight(Bundle%DeadTraj(ITraj)) * real(FMS_KineticT(Bundle%DeadTraj(ITraj)))
         end if
      end do

   end function FMS_KineticB

!>
!!    Computes potential energy of a trajectory bundle.
!!    \param Incoherent (Optional) If .True., computes an incoherent
!!    sum.  The default is .False., i.e. computes a coherent sum.
!!    @ingroup analysis
!<
   function FMS_PotentialB(Bundle, Incoherent) result(Energy)
      type(T_TrajectoryBundle) :: Bundle
      real(kind=DefReal) :: Energy
      logical, optional :: Incoherent
!      intent (in) Bundle
      integer(kind=DefInt) :: ITraj, JTraj
      logical :: IncSum

      complex(kind=DefComp) :: ctemp
      integer(kind=DefInt) :: CBFi, CBFj

      IncSum = .false.
      if (present(Incoherent)) IncSum = Incoherent

      ctemp = d1
      Energy = 0
!tjm  The following only works as expected when analytic integrals are
!tjm  both available and equivalent to second-order saddle point
      if (GliSaddle >= 2) GliSaddle = 0
      if (GliSaddle >= 2) GliSaddle = 0
      do ITraj = 1, Bundle%NumTraj

!     For incoherent sums, skip these off-diagonal terms
         if (.not. IncSum) then
            do JTraj = 1, ITraj - 1
               CBFi = Bundle%Trajectory(ITraj)%CBF
               CBFj = Bundle%Trajectory(JTraj)%CBF
!bfec
               if (glzCentroids .and. (CBFi /= CBFj)) then
                  Energy = Energy + 2.0 * real( &
                           conjg(Bundle%Trajectory(ITraj)%Amplitude) * Bundle%Trajectory(JTraj)%Amplitude * &
                           overlap_V(Bundle%Trajectory(ITraj), Bundle%Trajectory(JTraj), &
                                     Bundle%Centroids(((CBFi - 2) * (CBFi - 1)) / 2 + CBFj)))
!     &                 Bundle%Centroids(((ITraj-2)*(ITraj-1))/2+JTraj)))
               else
                  Energy = Energy + 2.0 * real( &
                           conjg(Bundle%Trajectory(ITraj)%Amplitude) * Bundle%Trajectory(JTraj)%Amplitude * &
                           overlap_V(Bundle%Trajectory(ITraj), Bundle%Trajectory(JTraj)))
               end if
!    $              FMS_CGVTT(
!    &                 Bundle%Trajectory(ITraj),
!    $                 Bundle%Trajectory(JTraj),centroid=
!    $                 Bundle%Centroids(((ITraj-2)*(ITraj-1))/2+JTraj)))
            end do
         end if

!     Diagonal terms
         Energy = Energy + FMS_Weight(Bundle%Trajectory(ITraj)) * real(overlap_V(Bundle%Trajectory(ITraj)))
!    &        FMS_CGVTT(Bundle%Trajectory(ITraj))
      end do

!     Add energy from dead trajectories
      do ITraj = 1, Bundle%NumDeadTraj

!        For incoherent sums, skip these off-diagonal terms
         if (.not. IncSum) then
            do JTraj = 1, ITraj - 1
               Energy = Energy + 2.0 * real( &
                        conjg(Bundle%DeadTraj(ITraj)%Amplitude) * &
                        Bundle%DeadTraj(JTraj)%Amplitude * &
                        (Bundle%DeadH(iTraj, jTraj) - overlap_KE(Bundle%DeadTraj(iTraj), Bundle%DeadTraj(jTraj))) &
                        )
            end do
         end if

!        Diagonal terms
         Energy = Energy + FMS_Weight(Bundle%DeadTraj(ITraj)) * real(overlap_V(Bundle%DeadTraj(ITraj)))
      end do

   end function FMS_PotentialB

!>
!!    Computes Population on each electronic state.
!!
!!    @ingroup analysis
!!    @param Ratios (Output) calculated populations
!!    @param LOrth (Optional) allows the user to force an incoherent sum
!!    LOrth=.true.   uses the true overlap matrix
!!    Recommended if amplitudes are symmetrically orthogonalized
!!    LOrth=.false.  approximates overlap matrix with identity
!<
   subroutine FMS_Branching(B1, Ratios, LOrth)
      type(T_TrajectoryBundle) :: B1
      real(kind=DefReal) :: Ratios(B1%NumStates)
      logical, optional :: LOrth
      complex(kind=DefComp) :: CTemp
      integer(kind=DefInt) :: ITraj, JTraj
      logical :: LOrthT

      LOrthT = .true.
      if (present(LOrth)) then
         LOrthT = LOrth
      end if

      Ratios = 0
      if (LOrthT) then
         do ITraj = 1, B1%NumTraj
            do JTraj = 1, ITraj - 1
               Ctemp = overlap(B1%Trajectory(ITraj), B1%Trajectory(JTraj), same_state=.true.)
               CTemp = CTemp * conjg(B1%Trajectory(ITraj)%Amplitude) * B1%Trajectory(JTraj)%Amplitude
               Ratios(B1%Trajectory(ITraj)%StateID) = Ratios(B1%Trajectory(ITraj)%StateID) + d2 * real(CTemp)
            end do
            Ratios(B1%Trajectory(ITraj)%StateID) = Ratios(B1%Trajectory(ITraj)%StateID) + &
                                                   real(conjg(B1%Trajectory(ITraj)%Amplitude) * &
                                                        B1%Trajectory(ITraj)%Amplitude)
         end do
      else
         do ITraj = 1, B1%NumTraj
            Ratios(B1%Trajectory(ITraj)%StateID) = Ratios(B1%Trajectory(ITraj)%StateID) + &
                                                   FMS_WeightC(B1%Trajectory(ITraj)%Amplitude)
         end do
      end if

!     Add terms for dead trajectories
      if (LOrthT) then
         do ITraj = 1, B1%NumDeadTraj
            do JTraj = 1, ITraj - 1
               CTemp = overlap(B1%DeadTraj(ITraj), B1%DeadTraj(JTraj), same_state=.true.)
               CTemp = CTemp * conjg(B1%DeadTraj(ITraj)%Amplitude) * B1%DeadTraj(JTraj)%Amplitude
               Ratios(B1%DeadTraj(ITraj)%StateID) = Ratios(B1%DeadTraj(ITraj)%StateID) + d2 * real(CTemp)
            end do
            Ratios(B1%DeadTraj(ITraj)%StateID) = Ratios(B1%DeadTraj(ITraj)%StateID) + &
                                                 real(conjg(B1%DeadTraj(ITraj)%Amplitude) * &
                                                      B1%DeadTraj(ITraj)%Amplitude)
         end do
      else
         do ITraj = 1, B1%NumDeadTraj
            Ratios(B1%DeadTraj(ITraj)%StateID) = Ratios(B1%DeadTraj(ITraj)%StateID) + &
                                                 FMS_WeightC(B1%DeadTraj(ITraj)%Amplitude)
         end do
      end if

   end subroutine FMS_Branching

end module BundleCalcsModule
