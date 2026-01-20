!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Trajectory class methods that related to electronic structure
!!         and related calculations.
!!
!! Methods for the trajectory class involve communication with the
!! electronic structure package, and retrieval of electronic structure
!! quantities at the given molecular configuration. These quantities are
!! encapsulated within the scope of this class, but are accessible
!! through public access routines FMS_PotEn, FMS_Forces, etc.
!!
!! Communication with the ESP takes place through the UpdatePES routine;
!! if an electronic quantity is requested which has not yet been
!! calculated at the current molecular configuration, UpdatePES is
!! called. This, in turn, calls the appropriate interface to the ESP,
!! which triggers an electronic structure calculation.
!! \image latex "../sources/TrajModule.png"
!<
module TrajectoryCalcsModule

   use GlobalModule
   use ParticleModule
   use ToyModelModule, only: FMS_ToyModel

   use TrajectoryModule

   use QM_MM_Module, only: qczQMMM, qczPCharge

   implicit none

   private
   public :: FMS_Coupling, FMS_SOCoupling
   public :: FMS_CIVec, FMS_Orbitals, FMS_SetIgnoreError
   public :: FMS_PosDot, FMS_MomDot, FMS_CoupDotVel
   public :: Potential, FMS_PotentialT, FMS_MMPot
   public :: Kinetic, FMS_KineticT
   public :: FMS_GetForce, FMS_Forces, FMS_MMForces
   public :: FMS_KineticClass, FMS_ClassEnergy
   public :: FMS_Weight, FMS_WeightC
   public :: FMS_Dipole, FMS_TransDipole
   public :: FMS_PhaseDot, FMS_TransDipoleIJxf
   public :: FMS_isBundleCurrent, FMS_BundleUpdated
   public :: FMS_isAmpDotCurrent, FMS_AmpDotUpdated

!--------------------------------------------------------------------!
!                 DERIVED TYPES FOR TRAJECTORIES                     !
!--------------------------------------------------------------------!

!     Function overloads:

   interface FMS_GetForce
      module procedure FMS_ForceVecAvH
   end interface

   interface FMS_Coupling
!         module procedure FMS_Coupling_sub,
      module procedure FMS_Coupling_func
   end interface

   interface FMS_SOCoupling
      module procedure FMS_SOCoupling_func
   end interface

   interface FMS_ClassEnergy
      module procedure FMS_ClassEnergy_Trajectory
   end interface

   interface FMS_Forces
      module procedure FMS_ForceComp
      module procedure FMS_ForceVec
   end interface

   interface FMS_MMForces
      module procedure FMS_MMForceComp
      module procedure FMS_MMForceVec
   end interface

   interface FMS_MomDot
      module procedure FMS_MomDotComp
      module procedure FMS_MomDotVec
      module procedure FMS_MomDotAll
   end interface

   interface Potential
      module procedure FMS_PotentialAvH
   end interface

   interface FMS_Weight
      module procedure FMS_WeightT
   end interface

   interface Kinetic
      module procedure FMS_KineticT
   end interface

contains
!--------------- Module Procedures ---------------!

!     Triggers a recalculation - public
!>
!! Run an electronic structure calculation at the current geometry.
!!
!!     \note The optional flags reflect whether a given quantity is required.
!!     Those quantities may be calculated anyway, depending on the
!!     specifics of the electronic structure interface. The electronic
!!     potential is always calculated.
!<
   subroutine FMS_UpdatePES(T1, zForce, zCoup, zDip, zQuad, zTDip, iState, jState)
      use SMDModule, only: SMD_Mechano, smSMD
      use TeraChemModule, only: RunTerachem
      use QM_MM_Module
      use ToyModelModule, only: FMS_ToyModel

      type(T_Trajectory) :: T1
      logical, optional :: zForce, zCoup, zDip, zQuad, zTDip
      integer(kind=DefInt), optional :: iState, jState
      logical :: zCalcForce, zCalcCoup, zCalcDip, zCalcQuad, zCalcTDip
      real(kind=DefReal) :: time_tmp1, time_tmp2
      integer(kind=DefInt) :: iCalcState, jCalcState

!     If this trajectory is dead, we should not be here
      if (T1%is_dead()) then
         write (fmiOut, *) T1%TrajID
         write (fmiOut, *) T1%ParentID
         write (fmiOut, *) T1%DeadTime
         call T1%print_id()
         write (fmiOut, *) 'Warning in UpdatePES: called with dead trajectory.'
      end if

!     If this timestep will be rejected, do not bother
!     updating anything
      if (FMS_StepRejected()) return

!bfec
      if (.not. glzCentroids) then
         if (T1%zCent) return
      end if

      call cpu_time(time_tmp1)
!     For QM/MM, partial charge interactions must be turned off
      if (qczQMMM .and. qczPCharge) call FMS_DieError('In UpdatePES - partial charge interactions must be off.')

!
!     Figure out what needs to be calculated
!
      if (present(zForce)) then
         zCalcForce = zForce
      else
         zCalcForce = .false.
      end if

      if (present(zCoup)) then
         zCalcCoup = zCoup
      else
         zCalcCoup = .false.
      end if

      if (present(zDip)) then
         zCalcDip = zDip
      else
         zCalcDip = .false.
      end if

      if (present(zQuad)) then
         zCalcQuad = zQuad
      else
         zCalcQuad = .false.
      end if

      if (present(zTDip)) then
         zCalcTDip = zTDip
      else
         zCalcTDip = .false.
      end if

      if (present(iState)) then
         iCalcState = iState
      else
         iCalcState = T1%StateID
      end if

      if (present(jState)) then
         jCalcState = jState
      else
         jCalcState = T1%StateID
      end if

!     For calculations requiring a 2nd state, make sure it was passed
      if ((.not. present(iState)) .and. (.not. present(jState))) then
         if (zCalcCoup) call FMS_DieError('Coupling requested, but second state not provided.')
         if (zCalcTDip) call FMS_DieError('Transition Dipole requested, but second state not provided.')
      end if

!     Copy current potential energy to old potential energy (used in check CI vectors)
      T1%ElecStruc%OldPotEn = T1%ElecStruc%PotEn

!
!     Do the calculations
!

      select case (glIModel)

      case (TOY)
         call FMS_ToyModel(T1)

      case (TEMPLATE)
         call FMS_InterfaceTemplate(T1)

      case (TC)
         call RunTerachem(T1, iCalcState, jCalcState, zCalcCoup)

      case (QUANTICS)
         call RunQuantics(T1, iCalcState, jCalcState, zCalcCoup)
      case default
         call FMS_DieError("Invalid iModel variable")

      end select

!
!     Calculate force modifications, if necessary
!
      if (smSMD) then
         if (.not. T1%ESFlags%zModPotCurrent) then
            T1%ElecStruc%ModForce = 0.d0
            T1%ElecStruc%ModPot = 0.d0
            if (smSMD) call SMD_Mechano(T1)
            T1%ESFlags%zModPotCurrent = .true.
         end if
      else
         T1%ElecStruc%ModForce = 0.d0
         T1%ElecStruc%ModPot = 0.d0
      end if

!     Record timing statistics
      glNESCalls = glNESCalls + 1
      call cpu_time(time_tmp2)
      estime = estime + time_tmp2 - time_tmp1
   end subroutine FMS_UpdatePES

! Public routines to return electronic structure quantities
!> Wrapper part for Average Hamiltonian stuff. Normal FMS will just go
!straight through this and return usual Potential. AvH dynamics will
!return the averaged potential

   function FMS_PotentialAvH(T1, IState) result(PotEn)
      type(T_Trajectory) :: T1
      integer(kind=DefInt), optional :: IState
      real(kind=DefReal) :: PotEn

      !locals
      integer(kind=DefInt) :: ICalcState
      integer(kind=DefInt) :: navstate, i

      PotEn = 0.0d0
      navstate = 1

      if (present(iState)) then !FMS is asking for specific state potential,
         !don't do AvH stuff
         PotEn = FMS_PotentialT(T1, IState)
      else
         if (glzAvH) then !Average the Potential
            navstate = gliAvHNStates
            do i = 1, navstate
               ICalcState = gliAvHStates(i)
               PotEn = PotEn + FMS_PotentialT(T1, ICalcState)
            end do
            PotEn = PotEn / real(navstate, kind=DefReal)
         else !don't include state information, and pass through
            PotEn = FMS_PotentialT(T1)
         end if
      end if

   end function FMS_PotentialAvH

!>
!! Recalculate potential energies if necessary and
!! return one of them.
!!
!! Scope: public
!<
   function FMS_PotentialT(T1, IState) result(PotEn)
      type(T_Trajectory) :: T1
      integer(kind=DefInt), optional :: IState
      real(kind=DefReal) :: PotEn

      if (.not. T1%ESFlags%ZPotEnCurrent) call FMS_UpdatePES(T1)

      if (present(IState)) then
         PotEn = T1%ElecStruc%PotEn(IState)
      else
         PotEn = T1%ElecStruc%PotEn(T1%StateID)
      end if
      if (FMS_StepRejected()) return

!     Force modifications
      PotEn = PotEn + T1%ElecStruc%ModPot

      if (.not. T1%ESFlags%ZPotEnCurrent) call FMS_DieError('ERROR: Called Update PES, but Potential not updated.')

   end function FMS_PotentialT

!> Wrapper part for Average Hamiltonian stuff. Normal FMS will just go
!straight through this and return usual force vector. AvH dynamics will
!return the averaged force

   function FMS_ForceVecAvH(T1, IState) result(ForceVector)
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: ForceVector(T1%NumDimensions)
      integer(kind=DefInt), optional :: IState

      !locals
      integer(kind=DefInt) :: ICalcState
      integer(kind=DefInt) :: navstate, i

      ForceVector = 0.0d0
      navstate = 1

      if (present(iState)) then !FMS is asking for specific state force,
         !don't do AvH stuff
         ForceVector = FMS_ForceVec(T1, IState)
      else
         if (glzAvH) then !Average the force
            navstate = gliAvHNStates
            do i = 1, navstate
               ICalcState = gliAvHStates(i)
               ForceVector = ForceVector + FMS_ForceVec(T1, ICalcState)
            end do
            ForceVector = ForceVector / real(navstate, kind=DefReal)
         else !don't include state information, and pass through
            ForceVector = FMS_ForceVec(T1)
         end if
      end if

   end function FMS_ForceVecAvH

!>
!! Recalculate ForceVector if not current, then return it.
!!
!! Scope: public
!<
   function FMS_ForceVec(T1, IState) result(ForceVector)
      use ElecStrucModule, only: eszAnalyticForces
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: ForceVector(T1%NumDimensions)
      integer(kind=DefInt), optional :: IState
      integer(kind=DefInt) :: iCalcState

      !Should be temporary
      !MSS: necessary to avoid NaN if step rejected
      ForceVector = 0.d0

!     Default is to calculate forces for current state
      if (present(iState)) then
         iCalcState = iState
      else
         iCalcState = T1%StateID
      end if

!     Redo electronic structure if we don't have the right state
!       Also, check Potential is current - can be silly instances where it thinks
!       Force is current, but it is not really
!       WJG: we should track down where those instances are...
!      if(.not.T1%ESFlags%ZPotEnCurrent.or.
!     $           .not.T1%ESFlags%ZDerivCurrent(iCalcState,iCalcState)) then
      if (.not. T1%ESFlags%ZDerivCurrent(iCalcState, iCalcState)) then
         if (eszAnalyticForces) then
            call FMS_UpdatePES(T1, zForce=.true., iState=iCalcState)
         else
            call FMS_NumGradient(T1, iCalcState)
         end if
      end if
      if (FMS_StepRejected()) return

!     Return forces
      ForceVector = -T1%ElecStruc%DerivMat(iCalcState, iCalcState, :)

!     Add any force modifications
      ForceVector = ForceVector + T1%ElecStruc%ModForce

      if (.not. T1%ESFlags%ZDerivCurrent(iCalcState, iCalcState)) then
         call FMS_DieError('ERROR: Called Update PES, but forces not updated.')
      end if

   end function FMS_ForceVec

!>
!! Recalculate ForceVector if not current, then return a component
!!
!! Scope: public
!<
   function FMS_ForceComp(T1, IParticle, IDim, IState) result(ForceComp)
      use ElecStrucModule, only: eszAnalyticForces
      type(T_Trajectory), intent(inout) :: T1
      integer, intent(in) :: iDim, iParticle
      integer(kind=DefInt), intent(in), optional :: IState
      real(kind=DefReal) :: ForceComp
      integer(kind=DefInt) :: jDim, iCalcState

      !Should be temporary
      !MSS: necessary to avoid NaN if step rejected
      ForceComp = 0.d0

!     Default is to calculate forces for current state
      if (present(iState)) then
         iCalcState = iState
      else
         iCalcState = T1%StateID
      end if

!     Find dimension index
      if (iParticle > 1) then
         jDim = sum(T1%Particle(1:iParticle - 1)%NumDimensions) + iDim
      elseif (iParticle == 1) then
         jDim = iDim
      else
         call FMS_DieError("Invalid IParticle passed to FMS_MMForceComp")
         return
      end if

!     Redo electronic structure if we don't have the right state
      if (.not. T1%ESFlags%ZDerivCurrent(iCalcState, iCalcState)) then
         if (eszAnalyticForces) then
            call FMS_UpdatePES(T1, zForce=.true., iState=iCalcState)
         else
            call FMS_NumGradient(T1, iCalcState)
         end if
      end if
      if (FMS_StepRejected()) return

!     Return force component
      ForceComp = -T1%ElecStruc%DerivMat(iCalcState, iCalcState, jDim)

!     Add any force modifications
      ForceComp = ForceComp + T1%ElecStruc%ModForce(jDim)

      if (.not. T1%ESFlags%ZDerivCurrent(iCalcState, iCalcState)) then
         call FMS_DieError('ERROR: Called Update PES, but forces not updated.')
      end if

   end function FMS_ForceComp
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Coupling_func(T1, IState, JState) result(CoupVec)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState, JState
      real(kind=DefReal) :: CoupVec(T1%NumDimensions)

      call FMS_Coupling_sub(T1, IState, JState, CoupVec)
      return
   end function FMS_Coupling_func

   subroutine FMS_Coupling_sub(T1, IState, JState, CoupVec)
!----------------------------------------------------------
!     Recalculate given coupling if not current, then return it.
!
!     Scope: public
!---------------------------------------------------------
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState, JState
      real(kind=DefReal) :: CoupVec(T1%NumDimensions)
!     Is this a useful calculation?
      if (IState == JState) then
!        call FMS_DieError('IState=JState in FMS_Coupling.')
         CoupVec = 0.0
         return
      end if

!     Set NACs between singlet and triplet states to 0. gets taken
!     care of by FMS_COCoup
      if ((IState <= NSing .and. JState > NSing) .or. (JState <= NSing .and. IState > NSing)) then
         CoupVec = 0.d0
         return
      end if

!WJG we can probably remove the following error seeing as all couplings can be stored
      if (IState /= T1%StateID .and. JState /= T1%StateID) then
         write (fmiOut, *) Istate, JState, T1%StateID
         call FMS_DieError('FMS_Coupling: Cannot calculate coupling unless one state is current state of trajectory.')
      end if

!     If the coupling is not current, calculate it
      if (.not. T1%ESFlags%ZDerivCurrent(IState, JState)) then
         call FMS_UpdatePES(T1, zCoup=.true., iState=iState, jState=jState)
      end if

!     Return Coupling for the requested state
      CoupVec(1:T1%NumDimensions) = T1%ElecStruc%DerivMat(IState, JState, 1:T1%NumDimensions)

      return
   end subroutine FMS_Coupling_sub

!>
!!    Calculates coupling dot velocity.
!!    Coupling is from \e T1%StateID to \e iState.
!!    @ingroup spawning
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_CoupDotVel(T, j) result(dCouple)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory) :: T
      integer(kind=DefInt) :: j
      real(kind=DefReal) :: dCouple

      integer(kind=DefInt) :: i

      i = T%StateID

      !  No coupling between a state and itself
      DCouple = 0.0d0
      if (j == i) return

      ! Calculate Coupling . Velocity
      dCouple = dot_product(T%get_vel(), FMS_Coupling(T, i, j))

   end function FMS_CoupDotVel

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_SOCoupling_func(T1, IState, JState, Msi, Msj) result(SOME)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState, JState, Msi, Msj
      complex(kind=DefComp) :: SOME

      call FMS_SOCoupling_sub(T1, IState, JState, Msi, Msj, SOME)
      return
   end function FMS_SOCoupling_func

   subroutine FMS_SOCoupling_sub(T1, IState, JState, Msi, Msj, SOME)
!----------------------------------------------------------
!     Recalculate given coupling if not current, then return it.
!
!     Scope: public
!---------------------------------------------------------
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState, JState, Msi, Msj
      complex(kind=DefComp) :: SOME

!     Is this a useful calculation?
!      REMOVED FOR GAIMS NOW. ADD LATER!!!!!
!      if (IState==JState) then
!        call FMS_DieError('IState=JState in FMS_Coupling.')
!         CoupVec = 0.0
!         return
!      endif

!WJG we can probably remove the following error seeing as all couplings can be stored
      if (IState /= T1%StateID .and. JState /= T1%StateID) then
         write (fmiOut, *) Istate, JState, T1%StateID
         call FMS_DieError('FMS_Coupling: Cannot calculate coupling unless one state is current state of trajectory.')
      end if

!     If the coupling is not current, calculate it
      if (.not. T1%ESFlags%ZSOMCurrent(IState, JState, Msi, Msj)) then
!           call FMS_UpdatePES(T1,zCoup=.true.,iState=iState,
!     $                     jState=jState)
         call FMS_UpdatePES(T1, zCoup=.false., iState=iState, jState=jState)
         T1%ESFlags%ZSOMCurrent(IState, JState, Msi, Msj) = .true.
      end if

!     Return Coupling for the requested state
      SOME = T1%ElecStruc%SOMat(IState, JState, Msi, Msj)
!      write(*,'(4i5,f18.12,sp,f18.12,"i")') IState,JState,Msi,Msj,SOME

   end subroutine FMS_SOCoupling_sub

!>
!! Recalculate transition dipole between state i and the
!! ground state if necessary, then return it.
!!
!! Scope: public
!<
   function FMS_TransDipole(T1, IState) result(TDip)
      use ElecStrucModule, only: eszTransDipole
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState
      real(kind=DefReal) :: TDip(3)

      if (IState <= 1 .or. IState > T1%NumStates) then
         call FMS_DieError("No transition dipole for requested state")
      end if

      if (.not. eszTransDipole) then
         call FMS_DieError('ERROR: in FMS_TransDipole, but T-Dipoles not available')
      end if

! Recalculate if necessary
      if (.not. T1%ESFlags%ZTransDipsCurrent) then
         call FMS_UpdatePES(T1, zTDip=.true., iState=iState, jState=1)
      end if
      if (FMS_StepRejected()) return

      TDip = T1%ElecStruc%TransDipole(IState, 1:3)

      if (.not. T1%ESFlags%ZTransDipsCurrent) then
         call FMS_DieError('ERROR: Called Update PES, but TDips not updated.')
      end if

   end function FMS_TransDipole

   function FMS_TransDipoleIJxf(T1) result(TDipxf)
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState
      real(kind=DefReal) :: TDipxf(4)

      IState = T1%StateID
      if (.not. T1%ESFlags%ZTransDipsCurrentxf) then
         call FMS_UpdatePES(T1, zTDip=.true., iState=IState)
      end if
      TDipxf = T1%ElecStruc%TransDipolexf(1:4)

      if (.not. T1%ESFlags%ZTransDipsCurrentxf) then
         call FMS_DieError('ERROR: Called Update PES, but xf TDips not updated.')
      end if
   end function FMS_TransDipoleIJxf

!>
!! Recalculate dipole for state I and return it
!!
!! Scope: public
!<
   function FMS_Dipole(T1, IState) result(dipole)
      use ElecStrucModule, only: eszDipoleMoment
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState
      real(kind=DefReal) :: dipole(3)

      if (.not. eszDipoleMoment) then
         call FMS_DieError('ERROR: in FMS_Dipole, but dipole moments are not avialable.')
      end if

      if (.not. T1%ESFlags%ZDipolesCurrent) then
         call FMS_UpdatePES(T1, zDip=.true., iState=iState)
      end if

      if (FMS_StepRejected()) return
      dipole = T1%ElecStruc%Dipole(IState, 1:3)

      if (.not. T1%ESFlags%ZDipolesCurrent) then
         call FMS_DieError('ERROR: Called Update PES, but dipoles not updated.')
      end if
   end function FMS_Dipole
!>
!! Recalculate wavefucntion if necessary and
!! return CIVector
!!
!! Scope: public
!<
   function FMS_CIVec(T1, IState, iCiVec) result(CICoeff)
      use ElecStrucModule, only: eslCIVec
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IState, ICIVec
      real(kind=DefReal) :: CICoeff

      if (eslCIVec <= 0) then
         call FMS_DieError('ERROR: in FMS_Orbitals, but no orbitals are available.')
      end if
      if (.not. T1%ESFlags%ZPotEnCurrent) call FMS_UpdatePES(T1)

      CICoeff = T1%ElecStruc%OldCIVecs(IState, ICIVec)

   end function FMS_CIVec

!>
!! Recalculate wavefucntion if necessary and
!! return CIVector
!!
!! Scope: public
!<
   function FMS_Orbitals(T1, iOrb, jOrb) result(OrbCoeff)
      use ElecStrucModule, only: esnBasis
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: iOrb, jOrb
      real(kind=DefReal) :: OrbCoeff

      if (esnBasis <= 0) then
         call FMS_DieError('ERROR: in FMS_Orbitals, but no orbitals are available.')
      end if

      if (.not. T1%ESFlags%ZPotEnCurrent) call FMS_UpdatePES(T1)

      OrbCoeff = T1%ElecStruc%OldOrbitals(iOrb, jOrb)

   end function FMS_Orbitals

!>
!! Recalculate wavefucntion if necessary and
!! return OverlapMatrix
!!
!! Scope: public
!<
   function FMS_OverlapMatrix(T1, iS, jS) result(OverlapElement)
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: iS, jS
      real(kind=DefReal) :: OverlapElement

      if (.not. T1%ESFlags%ZPotEnCurrent) call FMS_UpdatePES(T1)

      OverlapElement = T1%ElecStruc%OverlapMatrix(iS, jS)

   end function FMS_OverlapMatrix

!>
!! \todo Document
!<
   function FMS_PhaseDot(T1) result(PhaseDot)
      type(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: PhaseDot

!     dGamma/dt = T - V - alpha / ( 2 M )
      PhaseDot = FMS_KineticClass(T1) - FMS_PotentialT(T1) - 0.5d0 * sum(T1%get_width() / T1%get_mass())

!     Zero-point contribution -- Seems to help stabilize alpha-dependence
!     do IParticle=1,T1%NumParticles
!        do IDim=1,T1%Particle(IParticle)%NumDimensions
!           PhaseDot=PhaseDot-
!    $           dp5*T1%Particle(IParticle)%Width/
!    $           T1%Particle(IParticle)%Mass
!        enddo
!     enddo

   end function FMS_PhaseDot

   function FMS_ClassEnergy_Trajectory(T1) result(E)
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: E

      E = FMS_KineticClass(T1) + FMS_PotentialT(T1)
   end function FMS_ClassEnergy_Trajectory
!>
!! \todo Document
!<
   function FMS_KineticClass(Trajectory) result(Energy)
      type(T_Trajectory) :: Trajectory
      real(kind=DefReal) :: Energy
      intent(in) Trajectory
      integer(kind=DefInt) :: IDim, IParticle

      Energy = 0
      do IParticle = 1, Trajectory%NumParticles
         do IDim = 1, Trajectory%Particle(IParticle)%NumDimensions
            Energy = Energy + dp5 * Trajectory%Particle(IParticle)%get_mom(IDim) * &
                     Trajectory%Particle(IParticle)%get_mom(IDim) / Trajectory%Particle(IParticle)%mass
         end do
      end do
   end function FMS_KineticClass

!>
!! Returns component IDim of particle's PosDot
!!
!! \note Right now, PosDot is just the momentum/mass, but this
!! could be changed later.
!!
!! Scope: Public
!<
   function FMS_PosDot(P1, IDim) result(dPosDot)
      type(T_Particle), intent(in) :: P1
      integer(kind=DefInt), intent(in) :: IDim
      real(kind=DefReal) :: dPosDot

      dPosDot = P1%get_mom(IDim) / P1%mass
   end function FMS_PosDot

!>
!! Returns component IDim of particle's MomDot
!! Right now, we just call Forces, but this could be
!! changed later.
!!
!! Scope: Public
!<
   function FMS_MomDotComp(T1, IPart, IDim) result(dMomDot)
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: IPart, IDim
      real(kind=DefReal) :: dMomDot

      dMomDot = FMS_Forces(T1, IPart, IDim)

   end function FMS_MomDotComp

!>
!! Returns 3-vector of particle's MomDot.
!! Right now, we just call Forces, but this could be
!! changed later.
!!
!! Scope: Public
!<
   function FMS_MomDotVec(T1, iPart) result(dMomDotVec)
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: iPart
      real(kind=DefReal) :: dMomDotVec(T1%Particle(iPart)%NumDimensions)

      integer(kind=DefInt) :: iDim

      do iDim = 1, T1%Particle(iPart)%NumDimensions
         dMomDotVec(iDim) = FMS_Forces(T1, iPart, iDim)
      end do

   end function FMS_MomDotVec

!>
!! Returns vector of particle's MomDot.
!! Right now, we just call Forces, but this could be
!! changed later.
!!
!! Scope: Public
!<
   function FMS_MomDotAll(T1) result(dMomDotVec)
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: dMomDotVec(T1%NumDimensions)

      dMomDotVec = FMS_Forces(T1)

   end function FMS_MomDotAll

!>
!! Recalculate potential energies if necessary and
!! return one of them.
!!
!! Scope: public
!<
   function FMS_MMForceComp(T1, iParticle, iDim) result(ForceComp)
      use QM_MM_Module
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: ForceComp
      integer(kind=DefInt) :: iParticle, iDim, jDim

      if (.not. qczQMMM) call FMS_DieError('MMForce called, but this is not a QM/MM run.')

      if (.not. T1%ESFlags%zMMForceCurrent) then
         select case (GlIModel)

         case default
            call FMS_DieError('MMForce has no method to calculate MM potential for current model')

         end select

         if (qczConfine) call FMS_ConfinePot(T1)
         if (qczPCharge) call FMS_QMPCharge(T1)
      end if

      jDim = 0
!     Find dimension index
      if (iParticle > 1) then
         jDim = sum(T1%Particle(1:iParticle - 1)%NumDimensions) + iDim
      elseif (iParticle == 1) then
         jDim = iDim
      else
         call FMS_DieError("Invalid IParticle passed to FMS_MMForceComp")
      end if

      ForceComp = T1%ElecStruc%MMForce(jDim)
   end function FMS_MMForceComp

!>
!! Recalculate potential energies if necessary and
!! return one of them.
!!
!! Scope: public
!<
   function FMS_MMForceVec(T1) result(ForceVec)
      use QM_MM_Module
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: ForceVec(T1%NumDimensions)

      if (.not. qczQMMM) call FMS_DieError('MMForce called, but this is not a QM/MM run.')

      if (.not. T1%ESFlags%zMMForceCurrent) then
         select case (GlIModel)

         case default
            call FMS_DieError('MMForce has no method to calculate MM potential for current model')

         end select

         if (qczConfine) call FMS_ConfinePot(T1)
         if (qczPCharge) call FMS_QMPCharge(T1)
      end if

      ForceVec = T1%ElecStruc%MMForce

   end function FMS_MMForceVec
!>
!! Recalculate if necessary, and return, MM potential only.
!!
!! Scope: public
!<
   function FMS_MMPot(T1) result(PotEn)
      use QM_MM_Module
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: PotEn

      if (.not. qczQMMM) call FMS_DieError('MMPot called, but this is not a QM/MM run.')

      if (.not. T1%ESFlags%zMMPotCurrent) then
         select case (GlIModel)

         case default
            call FMS_DieError('MMPot has no method to calculate MM potential for current model')

         end select

         if (qczConfine) call FMS_ConfinePot(T1)
         if (qczPCharge) call FMS_QMPCharge(T1)

      end if

      PotEn = T1%ElecStruc%MMPot

   end function FMS_MMPot

!>
!!    Check that CIVector is consistent since last timestep.
!!
!!    Record any states flipping sign.
!!    Reject the timestep if two states have flipped character.
!<
   subroutine FMS_CheckCIVector(T1, CIV, OCIV, nstate, ncilen)
      type(T_Trajectory) :: T1

      integer(kind=DefInt) :: nstate, ncilen
      integer(kind=DefInt) :: IState, JState, ioccstate

      real(kind=DefReal) :: CIV(nstate, ncilen), OCIV(nstate, ncilen)
      real(kind=DefReal) :: SCI_II, SCI_IJ, SCI_JI, SCI_JJ
      real(kind=defReal) :: E_occ, E_I, E_J, E_occ_old, E_I_old, E_J_old

      logical :: ImportantState

1999  format('========================================================')
2000  format("WARNING: Trajectory jumped an intersection.")

!     Don't bother correcting phase if there was no electronic structure previously stored
      if (.not. T1%ESFlags%zESExists) then
         T1%ElecStruc%ElecPhase = 1.0d0
         return
      end if

      ioccstate = T1%StateID
!
!     Make sure the CI Vector is consistent since the last timestep
!

      do IState = 2, T1%NumStates
         do JState = 1, IState - 1

!
!    Calculate dot products SCI_KL=CIVec(State K)(t) .dot. CIVec(State L)(t-dt)
!
            SCI_II = dot_product(CIV(IState, :), OCIV(IState, :))
            SCI_JJ = dot_product(CIV(JState, :), OCIV(JState, :))
            SCI_IJ = dot_product(CIV(IState, :), OCIV(JState, :))
            SCI_JI = dot_product(CIV(JState, :), OCIV(IState, :))

            write (*, *) SCI_II, SCI_JI
            write (*, *) SCI_IJ, SCI_JJ
            write (*, *)

!            write(*,*) "Overlap elements for states:", IState, JState
!            write(*,*) "SCI_II", SCI_II, "SCI_JJ", SCI_JJ
!            write(*,*) "SCI_IJ", SCI_IJ, "SCI_JI", SCI_JI

            if (abs(SCI_II) + abs(SCI_JJ) < abs(SCI_IJ) + abs(SCI_JI)) then

!
!     States I and J  have flipped character
!     Attempt to reject the timestep

!     Is at least one of these states important (i.e. is it the occupied state, or
!     coupled to the occupied state)?  This is determined from the energy gaps.
!     Note: in the case of a centroid, there are two relevant states, but only
!     one is indexed in the data structure.  This is fine, however, as both states
!     will be considered important based on the energy criterion (if the energy
!     criterion is not met, then there is no coupling between the states, so it is
!     unnecessary to keep track of them.

               ImportantState = .false.
!     Gather energies
               E_occ = T1%ElecStruc%PotEn(ioccstate)
               E_I = T1%ElecStruc%PotEn(Istate)
               E_J = T1%ElecStruc%PotEn(Jstate)

!     Consider the previous timestep energy gaps as well, because a state flip
!     may have artificially taken a coupled state out of the coupling region
               E_occ_old = T1%ElecStruc%OldPotEn(ioccstate)
               E_I_old = T1%ElecStruc%OldPotEn(Istate)
               E_J_old = T1%ElecStruc%OldPotEn(Jstate)

               if (abs(E_I - E_occ) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_J - E_occ) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_I_old - E_occ_old) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_J_old - E_occ_old) < gldMaxEDiff) ImportantState = .true.

               write (fmiOut, 1999)
               write (fmiOut, 2000)
               call T1%print_id()
1001           format('SCI_', 2i2, ':', f10.6, ' SCI_', 2i2, ':', f10.6)
               write (fmiOut, 1001) jState, jState, SCI_JJ, jState, iState, SCI_JI
               write (fmiOut, 1001) iState, jState, SCI_IJ, iState, iState, SCI_II
               write (fmiOut, 1999)
!     If user option RejectAllstateFlip=.false., only reject the step if one of
!      the states was an important state
               if ((glzRejectAllStateFlip .or. ImportantState) .and. .not. T1%ESFlags%zIgnoreErrors) then
                  call FMS_RejectStep(.true.)
                  return
               end if
            end if

         end do
      end do

!     Correct phase based on diagonal overlap.
!     Note: this is also how we correct phase when two states have flipped.
!     The assumption that allows us to consider only the diagonal overlap in determining
!     the phase is that in a single step the states are not rotated by more than +-pi/4 .
!     Note - if this assumption does not hold then we would be unable to distinguish a rotation of
!     |t| > 90degree from a rotation of t-180 with two sign flips.

      do IState = 1, T1%NumStates
         SCI_II = dot_product(CIV(IState, :), OCIV(IState, :))
         if (SCI_II < 0.0d0) then
            write (6, '("  Sign of CIVecs changed in state",i3)') iState
            T1%ElecStruc%ElecPhase(iState) = -1.0d0 * T1%ElecStruc%ElecPhase(iState)
         end if
      end do

   end subroutine FMS_CheckCIVector

!      function FMS_CorrectPhaseMolpro(T1) Result (Success)
!-----------------------------------------------------------------------
!
!     This function ensures continuity of the phase of the NACME vector
!     returned by MPCoupVec.  It is also resposible for loading
!     the trajectories with orbitals and CI coefficients.
!
!-----------------------------------------------------------------------
!      use ElecStrucModule
!      type (T_Trajectory) T1
!      logical :: Success

!      success=.true.

!?>??ifdef (MolPro) then
!?>      real (kind=DefReal),allocatable,save::CIVecs(:,:), Orbitals(:,:)
!?>      real (kind=DefRealMP),allocatable,save:: getStuff(:)
!?>      real (kind=DefReal),allocatable,save::SAO(:,:)
!?>      real (kind=DefReal),allocatable,save::TmpDotP(:,:)
!?>      real (kind=DefReal) Tmp
!?>      integer (kind=DefInt) iorbs,ibasis,icasorbs,icaselec,icivec
!?>      integer (kind=DefInt) jorbs,jbasis,ITmp,i,IState,j
!?>      integer (kind=DefIntMP),parameter::mpzero=0
!?>      integer (kind=DefIntMP),parameter::mpone=1
!?>      integer (kind=DefIntMP),parameter::mptwo=2
!?>      integer (kind=DefIntMP),parameter::mpthree=3
!?>      integer (kind=DefIntMP),parameter::mpfive=5
!?>      integer (kind=DefIntMP),parameter::mpsix=6
!?>      integer (kind=DefIntMP) IDiff
!?>      character*256 string
!?>      real (kind=DefReal),allocatable,save::MSCIV(:,:),OMSCIV(:,:)
!?>! >>> debug option
!?>      integer (kind=DefIntMP) ndimtmp
!?>      integer (kind=DefIntMP),parameter::mphund=100
!?>! <<<
!?>
!?>      if (GlIMethod == 0) return !No phasing necessary
!?>
!?>      Success=.true.
!?>
!?>      if(.not.allocated(Orbitals)) then
!?>         allocate(Orbitals(esNBasis,esNBasis))
!?>         allocate(CIVecs(T1%NumStates,esLCiVec))
!?>         allocate(SAO(esNBasis,esNBasis))
!?>         allocate(TmpDotP(esNCasOrbs,esNCasOrbs))
!?>         allocate(GetStuff(max(esLCiVec*T1%NumStates,
!?>     $        esNBasis*esNBasis)))
!?>      endif
!?>!
!?>!     Read CI Vectors
!?>!
!?>      call readm(getStuff,int(esLCiVec*T1%NumStates,DefIntMP),mptwo,
!?>     $     int(2501,DefIntMP),mpzero,string)
!?>
!?>      do IState=1, T1%NumStates
!?>         do icivec=1, esLCiVec
!?>            CIVecs(IState,icivec)=getStuff((IState-1)*esLCiVec+icivec)
!?>         enddo
!?>      enddo
!?>
!?>
!?>!
!?>!     read MOs
!?>!
!?>      call read_info(int(2101,DefIntMP),mptwo,mpzero,idiff,string)
!?>      call read_orb(getStuff,mptwo)
!?>      do iorbs=1,esNBasis
!?>         do ibasis=1,esNBasis
!?>            Orbitals(iorbs,ibasis)=getStuff((iorbs-1)*esNBasis+ibasis)
!?>         enddo
!?>      enddo
!?>
!?>
!?>!
!?>!     Read overlap matrix
!?>!
!?>      call readm(getStuff,int((esNBasis*(esNBasis+1))/2,DefIntMP),mpone
!?>     $     ,int(1100,DefIntMP),mpzero,string)
!?>      do ibasis=1, esNBasis
!?>         do jbasis=1, ibasis
!?>            SAO(ibasis,jbasis)=getStuff((ibasis*
!?>     $         (ibasis-1))/2+jbasis)
!?>            if (ibasis.ne.jbasis) SAO(jbasis,ibasis)=
!?>     $         SAO(ibasis,jbasis)
!?>         enddo
!?>      enddo
!?>
!?>
!?>
!?>c     check to see if Old Orbitals exist.  There is no need for
!?>c     any of this on the first iterations, so just skip to the end
!?>      if(.not. T1%ESFlags%zESExists)then
!?>         T1%ElecStruc%ElecPhase=1.0d0
!?>         Success=.true.
!?>         goto 2100
!?>      endif
!?>
!?>
!?>c     calculate overlap matrix between MOs of current and previous
!?>c     timestep
!?>
!?>      do iorbs=esNOrbs-esNCasOrbs+1, esNOrbs
!?>         do jorbs=esNOrbs-esNCasOrbs+1,esNOrbs
!?>            i=iorbs-esNOrbs+esNCasOrbs
!?>            j=jorbs-esNOrbs+esNCasOrbs
!?>            TmpDotP(i,j)=0.0d0
!?>            do ibasis=1,esNBasis
!?>               do jbasis=1,esNBasis
!?>                  TmpDotP(i,j)=TmpDotP(i,j)+
!?>     $                 Orbitals(iorbs,ibasis)*
!?>     $                 T1%ElecStruc%OldOrbitals(jorbs,
!?>     $                 jbasis)*SAO(ibasis,jbasis)
!?>               enddo
!?>            enddo
!?>         enddo
!?>      enddo
!?>
!?>
!?>c     determine which orbitals of this time step correspond to which
!?>c     orbitals of the previous timestep, and also which orbital have
!?>c     changed sign
!?> 1999 format('========================================================')
!?> 2000 format('ERROR in CorrectPhase: Active space not diabatized.')
!?> 2001 format('Orbital ',i4,' has changed sign. Overlaps:')
!?> 2002 format('Orbitals ',i4,' and ',i4,' have flipped. Overlaps:')
!?>      do i=1,esNCasOrbs
!?>         ITmp=0
!?>         Tmp=0.0d0
!?>         do j=1,esNCasOrbs
!?>            if (abs(TmpDotP(i,j)).gt.Tmp) then
!?>               ITmp=j
!?>               Tmp=abs(TmpDotP(i,j))
!?>            endif
!?>         enddo
!?>         if ((ITmp.ne.i).and.(Tmp.gt.0.7d0)) then
!?>            Success=.false.
!?>            write(fmiOut,1999)
!?>            write(fmiOut,2000)
!?>            call T1%print_id()
!?>            write(fmiOut,2002) i, iTmp
!?>            do iorbs=1,esNCasOrbs
!?>               write(fmiOut,*) (TmpDotP(iorbs,j),j=1,esNCasOrbs)
!?>            enddo
!?>            write(fmiOut,1999)
!?>            return
!?>         endif
!?>         if (TmpDotP(i,i).lt.0.0d0) then
!?>            Success=.false.
!?>            write(fmiOut,1999)
!?>            write(fmiOut,2000)
!?>            call T1%print_id()
!?>            write(fmiOut,2001) i
!?>            do iorbs=1,esNCasOrbs
!?>               write(fmiOut,*) (TmpDotP(iorbs,j),j=1,esNCasOrbs)
!?>            enddo
!?>            write(fmiOut,1999)
!?>            return
!?>         endif
!?>
!?>      enddo
!?>
!?>
!?>!
!?>!     Make sure the CI Vector is consistent since the last timestep
!?>!
!?>!       In case of MSPT2, we try to keep track of MSPT2*CAS CI coefficients.
!?>!
!?>      if (gliModel==10 .and. GlIMethod==2) then
!?>         if(.not.allocated(MSCIV))  then
!?>            allocate(MSCIV(T1%NumStates,esLCiVec))
!?>            allocate(OMSCIV(T1%NumStates,esLCiVec))
!?>         endif
!?>         do icivec=1,esLCiVec
!?>            do IState=1,T1%NumStates
!?>               MSCIV(IState,icivec)=
!?>     &           DOT_PRODUCT(T1%ElecStruc%MSPT2C(IState,1:T1%NumStates),
!?>     &              CIVecs(1:T1%NumStates,icivec))
!?>               OMSCIV(IState,icivec)=
!?>     &           DOT_PRODUCT(T1%ElecStruc%OldMSPT2C(IState,
!?>     &              1:T1%NumStates),
!?>     &              T1%ElecStruc%OldCIVecs(1:T1%NumStates,icivec))
!?>            enddo
!?>         enddo
!?>
!?>! >>> debug output
!?>!     Check for state flipping
!?>         write(6,*) '>> MSPT2_NAC'
!?>         if(esLCiVec.gt.mphund) then
!?>            write(6,'("Print only first 100 CIvecs.")')
!?>            ndimtmp=mphund
!?>         else
!?>            ndimtmp=esLCiVec
!?>         endif
!?>         write(6,'(" >>> Check vecs in CorrectPhase (new and old):")')
!?>         write(6,'(" MSPT2 coefficients")')
!?>         do IState=1,T1%NumStates
!?>            do j=1,T1%NumStates
!?>               write(6,'(2i5,2f20.9)') IState,j,
!?>     &              T1%ElecStruc%MSPT2C(IState,j),
!?>     &              T1%ElecStruc%OldMSPT2C(IState,j)
!?>            enddo
!?>         enddo
!?>         write(6,'(" CI coefficients (only CI part)")')
!?>         do IState=1,T1%NumStates
!?>            do icivec=1,ndimtmp
!?>               write(6,'(2i5,2f20.9)') IState,icivec,
!?>     &              CIVecs(IState,icivec),
!?>     &              T1%ElecStruc%OldCIVecs(IState,icivec)
!?>            enddo
!?>         enddo
!?>         write(6,'(" MSPT2*CI vector")')
!?>         do IState=1,T1%NumStates
!?>           do icivec=1,ndimtmp
!?>               write(6,'(2i5,2f20.9)') IState,icivec,
!?>     &              MSCIV(IState,icivec),OMSCIV(IState,icivec)
!?>            enddo
!?>         enddo
!?>         write(6,'(" <<<")')
!?>         write(6,*) '<< MSPT2_NAC'
!?>! <<<
!?>
!?>         !Check for state flips and correct any sign flips
!?>         call FMS_CheckCIVector(T1,MSCIV,OMSCIV,T1%NumStates,esLCiVec)
!?>         if(FMS_StepRejected()) return
!?>
!?>      else
!?>!
!?>!       else do usual CASSCF phase check
!?>!
!?>!        Check for state flips and correct any sign flips
!?>         call FMS_CheckCIVector(T1,CIVecs,T1%ElecStruc%OldCIVecs,
!?>     &                          T1%NumStates,esLCiVec)
!?>         if(FMS_StepRejected()) return
!?>
!?>      endif
!?>
!?>!
!?>!     Copy electronic structure into memory
!?>!
!?> 2100 continue
!?>      do IState=1, T1%NumStates
!?>         do icivec=1, esLCiVec
!?>            T1%ElecStruc%OldCIVecs(IState,icivec)=
!?>     $           CIVecs(IState,icivec)
!?>         enddo
!?>      enddo
!?>      do iorbs=1,esNBasis
!?>         do ibasis=1,esNBasis
!?>            T1%ElecStruc%OldOrbitals(iorbs,ibasis)=
!?>     $           Orbitals(iorbs,ibasis)
!?>         enddo
!?>      enddo
!?>! save MSPT2 coefficients as well
!?>      if(GlIMethod == 2) then
!?>         do IState=1,T1%NumStates
!?>            do j=1,T1%NumStates
!?>               T1%ElecStruc%OldMSPT2S(j,Istate) =
!?>     &              T1%ElecStruc%MSPT2S(j,Istate)
!?>               T1%ElecStruc%OldMSPT2C(j,Istate) =
!?>     &              T1%ElecStruc%MSPT2C(j,Istate)
!?>            enddo
!?>         enddo
!?>      endif
!?>
!?>      T1%ESFlags%zESExists=.true.
!?>      T1%ESFlags%zIgnoreErrors=.false.
!?>
!?>??endif
!!      end function FMS_CorrectPhaseMolpro
!>
!!     Finite difference to get gradient - either forward
!!     (not recommended for dynamics!) or central.
!!     For "forward" gradients, the sign of the displacement is reversed
!!     at each timestep
!!
!! Scope: Trajectory Module
!<
   subroutine FMS_NumGradient(T1, IState)
      type(T_Trajectory) :: T1
      type(T_Trajectory), save :: TTemp
      real(kind=DefReal), allocatable, save :: Gradient(:)
      real(kind=DefReal) :: PosTemp, EPlus, EMinus, E0
      integer(kind=DefInt) :: IState, Index, IPtcle, IDim

      real(kind=DefReal), save :: LastTime = -1.d0, direction = 1.d0

!     Dynamic allocation
      if (.not. allocated(Gradient)) then
         call TTemp%create(T1%NumParticles, T1%NumStates)
         allocate (Gradient(T1%NumDimensions))
      end if

!     For non-central gradient, determine
!     whether to go forward or backward
      if (.not. glzCentNGrad .and. T1%get_time() /= LastTime) then
         Direction = -1.d0 * Direction
         LastTime = T1%get_time()
      end if

      TTemp = T1
      Index = 1
      E0 = 0.0d0
      if (.not. glzCentNGrad) E0 = FMS_PotentialT(T1, iState)
2000  format(' Numerically differentiating energy over ', i4, ' coordinates for trajectory # ', i3)
      write (fmiOut, 2000) T1%NumDimensions, T1%TrajID

      do IPtcle = 1, T1%NumParticles
         do IDim = 1, T1%Particle(IPtcle)%NumDimensions

            if (Direction > 0) then
               write (fmiOut, *) index, '+'
            else
               write (fmiOut, *) index, '-'
            end if
            flush (fmiOut)
            PosTemp = T1%get_pos(IPtcle, IDim) + gldNGradStep * Direction
            call TTemp%set_pos(IPtcle, IDim, PosTemp)
            EPlus = FMS_PotentialT(TTemp, IState)

!     For central finite difference
            if (glzCentNGrad) then
               write (fmiOut, *) index, '-'
               flush (fmiOut)
               PosTemp = T1%get_pos(IPtcle, IDim) - gldNGradStep
               call TTemp%set_pos(IPtcle, IDim, PosTemp)
               EMinus = FMS_PotentialT(TTemp, IState)

               Gradient(Index) = (EPlus - EMinus) / (2.0 * gldNGradStep)

!     For forward finite difference
            else
               Gradient(Index) = Direction * (EPlus - E0) / gldNGradStep
            end if

!     Reset TTEmp position
            TTemp = T1
            Index = Index + 1
         end do
      end do

      T1%ElecStruc%DerivMat(IState, IState, :) = Gradient
      T1%ESFlags%zDerivCurrent(IState, IState) = .true.
   end subroutine FMS_NumGradient

!>
!!    For QM/MM runs, add an external confining potential to the MM
!!    system.
!!
!!    \note Currently, the only implemented form of this potential is
!!    radial with
!!
!!    V(x)={0                                          , r < ConfineD
!!         {AtomicMass(amu) * ConfineK * (r-ConfineD)^2, r > ConfineD
!!
!!    The potential acts ONLY on the MM system. We multiply by atomic
!!    mass (in amu, not au) to avoid shearing apart polyatomic molecules.
!!
!!    This routine adds the potential and gradient to Trajectory
!!    structure (T1%ElecStruc%{MMPot,MMForce})
!<
   subroutine FMS_ConfinePot(T1)
      use QM_MM_Module
      type(T_Trajectory) :: T1

      integer(kind=DefInt) :: iPtcle
      real(kind=DefReal) :: DistSqr, Dist, force(3)

      do IPtcle = 1, T1%NumParticles

         DistSqr = sum(T1%Particle(IPtcle)%get_pos()**2)

         if (DistSqr < qcdConfineD**2) cycle

         dist = sqrt(distsqr)

         T1%ElecStruc%MMPot = T1%ElecStruc%MMPot + &
                              (T1%Particle(iPtcle)%Mass / 1822.89) * qcdConfineK * ((dist - qcdConfineD)**2)

         force = -2.0d0 * T1%Particle(iPtcle)%Mass * &
                 qcdConfineK * (dist - qcdConfineD) * T1%Particle(IPtcle)%get_pos() / (1822.89 * dist)

         T1%ElecStruc%MMForce(3 * IPtcle - 2:3 * IPtcle) = T1%ElecStruc%MMForce(3 * IPtcle - 2:3 * IPtcle) + force

      end do

   end subroutine FMS_ConfinePot

!>
!!     Point charge interactions between partial charges of QM system and
!!     point charges of MM system. This should not be called unless QM
!!     calculations have been turned off, e.g. for constrained
!!     equilibration.
!<
   subroutine FMS_QMPCharge(T1)
      use QM_MM_Module
      use ElecStrucModule, only: eszPartialCharges
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: iPtcle, jPtcle
      real(kind=DefReal) :: pForce(3), diff(3), dist

      if (.not. eszPartialCharges) call FMS_DieError('Cannot compute partial charge interactions.')

      do jPtcle = qcNumQM + 1, T1%NumParticles
         if (T1%Particle(jPtcle)%Charge == 0.0d0) cycle
         do iPtcle = 1, qcNumQM
            if (T1%Particle(iPtcle)%Charge == 0.0d0) cycle

            diff = T1%Particle(iPtcle)%get_pos() - T1%Particle(jPtcle)%get_pos()
            dist = sqrt(sum(diff**2))

!     Add potential contribution
            T1%ElecStruc%MMPot = T1%ElecStruc%MMPot + &
                                 T1%Particle(IPtcle)%Charge * T1%Particle(JPtcle)%Charge / dist

!     Add force contribution
            pForce = T1%Particle(iPtcle)%Charge * T1%Particle(jPtcle)%Charge * diff / (dist**3)
            T1%ElecStruc%MMForce(3 * IPtcle - 2:3 * IPtcle) = T1%ElecStruc%MMForce(3 * IPtcle - 2:3 * IPtcle) + pForce
            T1%ElecStruc%MMForce(3 * jPtcle - 2:3 * jPtcle) = T1%ElecStruc%MMForce(3 * jPtcle - 2:3 * jPtcle) - pForce
         end do
      end do

   end subroutine FMS_QMPCharge
!>
!! Run calculation using provided templates
!<
   subroutine FMS_InterfaceTemplate(T1)
      use TemplateModule, only: FMS_RunCommand, tpiGradSign, FMS_WriteInput, FMS_ReadOutput, FMS_TemplateClean
      use ElecStrucModule, only: esnBasis, eslCIVec, eszAnalyticForces, eszNACoupVec
      type(T_Trajectory) :: T1
      integer(kind=DefInt) :: iState, jState, iUnit

      real(kind=DefReal), allocatable, save :: Geom(:), NACV(:)
      real(kind=DefReal), allocatable, save :: CIVecAll(:), CIVecs(:, :)
      real(kind=DefReal), allocatable, save :: Orbitals(:, :)
      integer(kind=DefInt) :: iDim, jDim, iPtcle

!     Copy geometry and write input
      if (.not. allocated(Geom)) allocate (Geom(T1%NumDimensions))
      jDim = 0
      do iPtcle = 1, T1%NumParticles
         do iDim = 1, T1%Particle(iPtcle)%NumDimensions
            jDim = jDim + 1
            Geom(jDim) = T1%Particle(iPtcle)%get_pos(iDim)
         end do
      end do
      if (esnBasis > 0 .and. .not. T1%ESFlags%zESExists) then
         call FMS_WriteInput(Geom, T1%StateID, 'template.genorb.write')
      elseif (esnBasis > 0) then
         call FMS_WriteInput(Geom, T1%StateID, 'template.write', &
                             Orb=T1%ElecStruc%OldOrbitals, CIVec=T1%ElecStruc%OldCIVecs)
      else
         call FMS_WriteInput(geom, T1%StateID, 'template.write')
      end if

!     Call electronic structure routine
      call FMS_RunCommand(T1%StateID)

!     Get energies
      call FMS_ReadOutput(T1%ElecStruc%PotEn, T1%StateID, 'template.poten.read')
      T1%ElecStruc%PotEn = T1%ElecStruc%PotEn + gldEShift
      T1%ESFlags%zPotEnCurrent = .true.

!     Get forces
      if (eszAnalyticForces) then
         call FMS_ReadOutput(T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, :), T1%StateID, 'template.grad.read')
         T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, :) = -tpigradsign * &
                                                            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, :) !minus sign because DerivMat storees -force
         T1%ESFlags%zDerivCurrent(T1%StateID, T1%StateID) = .true.
      end if

!     Get non-adiabatic coupling vectors
      if (eszNACoupVec) then
!     allocate to read in all ci vecs at once
         if (.not. allocated(NACV)) allocate (NACV(T1%NumDimensions * (T1%NumStates - 1)))
         call FMS_ReadOutput(NACV, T1%StateID, 'template.coup.read')

         jState = 0
         do iState = 1, T1%NumStates
            if (iState == T1%StateID) cycle
            jState = jState + 1
            T1%ElecStruc%DerivMat(T1%StateID, IState, :) = &
               NACV((jState - 1) * T1%NumDimensions + 1:jState * T1%NumParticles)
            !Anti-Hermitian
            T1%ElecStruc%DerivMat(IState, T1%StateID, :) = -T1%ElecStruc%DerivMat(T1%StateID, IState, :)
            T1%ESFlags%zDerivCurrent(T1%StateID, iState) = .true.
            T1%ESFlags%zDerivCurrent(iState, T1%StateID) = .true.
         end do
      end if

!    Get CI vector
      if (.not. allocated(CIVecAll)) allocate (CIVecAll(eslCIVec * T1%NumStates))
      if (.not. allocated(CIVecs)) allocate (CIVecs(T1%NumStates, eslCIVec))
      call FMS_ReadOutput(CIVecAll, T1%StateID, 'template.civec.read')
      do iState = 1, T1%NumStates
         CIVecs(iState, :) = CIVecAll((iState - 1) * eslCIVec + 1:iState * eslCIVec)
      end do

!     Read orbitals from dump file
      if (esnBasis > 0) then
         if (.not. allocated(Orbitals)) then
            allocate (Orbitals(esnBasis, esnBasis))
         end if
         open (file=trim(FMSWorkingDir)//'basis.dat', status='old', newunit=iUnit)
         if (esnBasis * esnBasis > 99999) then
            call FMS_DieError('Them are some big orbitals.')
         end if
         read (iUnit, *) T1%ElecStruc%OldOrbitals
         close (iUnit, status='delete')
      end if

!c$$$!     Read sao matrix from dump file
!c$$$      if(esnBasis>0 .and. esnCASOrbs>0 .and. nActiveOrbs>0) then
!c$$$         if(.not. allocated(SAO)) allocate(SAO(esnBasis,esnBasis))
!c$$$         open(file=trim(FMSWorkingDir)//'SAO.dat',status='old',
!c$$$     $        newunit=iUnit)
!c$$$         if(esnBasis*esnBasis > 99999) call FMS_DieError(
!c$$$     $        'Them are some big orbitals.')
!c$$$         read(iUnit,*) SAO
!c$$$         close(iUnit,status='delete')
!c$$$
!c$$$!     Re-order orbitals
!c$$$         do iorbs=esNOrbs-esNCasOrbs+1, esNOrbs
!c$$$            do jorbs=esNOrbs-esNCasOrbs+1,esNOrbs
!c$$$               i=iorbs-esNOrbs+esNCasOrbs
!c$$$               j=jorbs-esNOrbs+esNCasOrbs
!c$$$               TmpDotP(i,j)=0.0d0
!c$$$               do ibasis=1,esNBasis
!c$$$                  do jbasis=1,esNBasis
!c$$$                     TmpDotP(i,j)=TmpDotP(i,j)+
!c$$$     $                    Orbitals(iorbs,ibasis)*
!c$$$     $                    T1%ElecStruc%OldOrbitals(jorbs,
!c$$$     $                    jbasis)*SAO(ibasis,jbasis)
!c$$$                  enddo
!c$$$               enddo
!c$$$            enddo
!c$$$         enddo
!c$$$      do i=1,esNCasOrbs
!c$$$          ITmp=0
!c$$$         Tmp=0.0d0
!c$$$         do j=1,esNCasOrbs
!c$$$            if (abs(TmpDotP(i,j)).gt.Tmp) then
!c$$$               ITmp=j
!c$$$               Tmp=abs(TmpDotP(i,j))
!c$$$            endif
!c$$$         enddo
!c$$$         if ((ITmp.ne.i).and.(Tmp.gt.0.7d0)) then
!c$$$            Success=.false.
!c$$$            write(fmiOut,1999)
!c$$$            write(fmiOut,2000)
!c$$$            call T1%print_id()
!c$$$            write(fmiOut,2002) i, iTmp
!c$$$            do iorbs=1,esNCasOrbs
!c$$$               write(fmiOut,*) (TmpDotP(iorbs,j),j=1,esNCasOrbs)
!c$$$            enddo
!c$$$            write(fmiOut,1999)
!c$$$            return
!c$$$         endif
!c$$$         if (TmpDotP(i,i).lt.0.0d0) then
!c$$$            Success=.false.
!c$$$            write(fmiOut,1999)
!c$$$            write(fmiOut,2000)
!c$$$            call T1%print_id()
!c$$$            write(fmiOut,2001) i
!c$$$            do iorbs=1,esNCasOrbs
!c$$$               write(fmiOut,*) (TmpDotP(iorbs,j),j=1,esNCasOrbs)
!c$$$            enddo
!c$$$            write(fmiOut,1999)
!c$$$            return
!c$$$         endif
!c$$$      enddo
!c$$$
!c$$$      endif
      call FMS_CheckCIVector(T1, CIVecs, T1%ElecStruc%OldCIVecs, T1%NumStates, esLCiVec)
      T1%ElecStruc%OldCIVecs = CIVecs
      T1%ESFlags%zESExists = .true.
      T1%ESFlags%zIgnoreErrors = .false.

!     Clean up template files
      call FMS_TemplateClean()

   end subroutine FMS_InterfaceTemplate

   subroutine FMS_DummySetPosition(Points, NumQM, NumDummy, Weight)
!-----------------------------------------------------------------------
!
!     set positions of dummy atoms
!
!-----------------------------------------------------------------------
      integer(kind=DefInt) :: NumQM, NumDummy, idum, jqm
      real(kind=DefReal) :: Points(3, NumQM + NumDummy), Weight(NumQM, NumDummy)

      !
      ! set geometry of dummy atoms
      !
      do idum = 1, NumDummy
         Points(:, NumQM + idum) = 0.0d0
         do jqm = 1, NumQM
            Points(:, NumQM + idum) = Points(:, NumQM + idum) + Points(:, jqm) * Weight(jqm, idum)
         end do
      end do

   end subroutine FMS_DummySetPosition

   subroutine FMS_DummyCollectGrad(Points, NumQM, NumDummy, Weight)
!-----------------------------------------------------------------------
!
!     collect gradient on dummy atoms & map them into real atoms
!
!-----------------------------------------------------------------------
      integer(kind=DefInt) :: NumQM, NumDummy, idum, jqm, ixyz
      real(kind=DefReal) :: Points(3, NumQM + NumDummy), Weight(NumQM, NumDummy), com(3), com2(3)

      !
      ! retrieve gradient info and pack into real atoms
      !
      ! just for debugging
      do ixyz = 1, 3
         com(ixyz) = sum(Points(ixyz, :))
      end do

      do idum = 1, NumDummy
         do jqm = 1, NumQM
            Points(:, jqm) = Points(:, jqm) + Weight(jqm, idum) * Points(:, NumQM + idum)
         end do
         Points(:, NumQM + idum) = 0.0d0
      end do

      !> debug: check if sum of gradient conserves
      do ixyz = 1, 3
         com2(ixyz) = sum(Points(ixyz, :))
      end do

      do ixyz = 1, 3
         if (abs(com(ixyz) - com2(ixyz)) > 1.0d-10) then
            write (6, *) 'Weight for dummy atom set improperly?'
            write (6, '(A,3f15.10)') 'sum of original gradient:', com(:)
            write (6, '(A,3f15.10)') 'sum of modified gradient:', com2(:)
            call FMS_DieError('Error in AdjustDummy')
         end if
      end do

   end subroutine FMS_DummyCollectGrad

   subroutine RunQuantics(T, iCalcState, jCalcState, CalcCoup)
      use QuanticsModule, only: run_quantics
      type(T_Trajectory), intent(inout) :: T
      integer(kind=DefInt), intent(in) :: iCalcState, jCalcState
      logical, intent(in) :: CalcCoup

      real(kind=DefReal), allocatable :: xyz(:, :), gra(:, :)
      real(kind=DefReal), allocatable :: en(:)
      real(kind=DefReal), allocatable :: nadvec(:, :, :)

      integer :: natoms, nstates

      natoms = T%NumParticles
      nstates = T%NumStates

      allocate (En(T%NumStates))
      allocate (xyz(natoms, 3))
      allocate (gra(natoms, 3))
      allocate (nadvec(natoms, natoms, nstates))

      call run_quantics(step=1, xyz0=xyz, cstate=iCalcState, en=En, gra=gra, nadvec=nadvec)

   end subroutine RunQuantics

!     Bundle scope routines
   function FMS_IsBundleCurrent(T1) result(zCurrent)
      type(T_Trajectory), intent(IN) :: T1
      logical :: zCurrent

      zCurrent = T1%BFlags%zBundleCurrent

   end function FMS_IsBundleCurrent

   subroutine FMS_BundleUpdated(T1)
      type(T_Trajectory) :: T1

      T1%BFlags%zBundleCurrent = .true.

   end subroutine FMS_BundleUpdated

   function FMS_IsAmpDotCurrent(T1) result(zCurrent)
      type(T_Trajectory), intent(IN) :: T1
      logical :: zCurrent

      zCurrent = T1%BFlags%zAmpDotCurrent

   end function FMS_IsAmpDotCurrent

   subroutine FMS_AmpDotUpdated(T1)
      type(T_Trajectory) :: T1

      T1%BFlags%zAmpDotCurrent = .true.

   end subroutine FMS_AmpDotUpdated

   subroutine FMS_SetIgnoreError(T1, flag)
      type(T_Trajectory) :: T1
      logical :: flag

      T1%ESFlags%ZIgnoreErrors = flag
   end subroutine FMS_SetIgnoreError

!>
!!     Computes the kinetic Energy of a trajectory
!!    @ingroup analysis
!<
   function FMS_KineticT(Trajectory) result(Energy)
      type(T_Trajectory) :: Trajectory
      real(kind=DefReal) :: Energy
      intent(in) Trajectory
      integer(kind=DefInt) :: IDim, IParticle

      Energy = 0
      do IParticle = 1, Trajectory%NumParticles
         do IDim = 1, Trajectory%Particle(IParticle)%NumDimensions
            Energy = Energy + dp5 * Trajectory%Particle(IParticle)%get_mom(IDim) * &
                     Trajectory%Particle(IParticle)%get_mom(IDim) / Trajectory%Particle(IParticle)%Mass
         end do
      end do
   end function FMS_KineticT

!>
!!    Calculates probability weight of a trajectory
!!    @ingroup analysis
!<
   function FMS_WeightT(T1) result(Peso)
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: Peso

      Peso = FMS_WeightC(T1%Amplitude)
      return
   end function FMS_WeightT

   function FMS_WeightC(Amplitude) result(Peso)
      complex(kind=DefComp) :: Amplitude
      real(kind=DefReal) :: Peso

      Peso = real(conjg(Amplitude) * Amplitude)
   end function FMS_WeightC

end module TrajectoryCalcsModule
