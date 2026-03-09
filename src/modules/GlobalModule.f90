!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!> @brief General simulation-defining paramaters and general access function
!!
!! This module contains hard-coded global parameters for every
!! program run, and a run-time namelist set of parameters that
!! characterise an individual run.
!!
!! It also contains some general use utility functions
!<
module GlobalModule
   implicit none
   public
   save

!---------------------------------------------------------------------------
!                       COMPILE-TIME CONSTANTS
!---------------------------------------------------------------------------
! Parameters that are hard-coded in many places, and define dimensions
! of arrays that take too little memory to bother with dynamic allocation.
   integer :: idefault
   integer, parameter :: DefInt = kind(idefault) !<Default size of integer in bytes - allows for compiler flags to change this.
   integer, parameter :: DefInt4 = 4 !<4 byte integer

   double precision :: rdefault ! (allows for compiler flags to change this)
   integer(kind=DefInt), parameter :: DefReal = kind(rdefault) !<Default size of real in bytes
   ! TODO: Remove defintblas
   ! integer, parameter::DefIntBlas=4 !< Integer size for BLAS ! may need to match to system requirement (for fms_ch)
   integer(kind=DefInt), parameter :: DefIntBlas = DefInt !< Better to set blas integers to be default length.
   integer(kind=DefInt), parameter :: DefComp = 8 !< Default size of complex in bytes
   integer(kind=DefInt), parameter :: MaxTrajLimit = 1000 !< Maximum number of trajectories; 0 for no limit

! Physical Constants in AU and Unit Conversion Factors
   real(kind=DefReal), parameter :: kcalPMtoH = 0.0015952353d0 !< Conversion factor from kilocalorie per mole to Hartree
   real(kind=DefReal), parameter :: FsToAU = 41.34137221718d0 !< Conversion factor from femtosecond to atomic unit of time
   real(kind=DefReal), parameter :: BohrToAng = 0.529177249d0 !< Conversion factor from bohr (atomic unit of distance) to Angstrom
   real(kind=DefReal), parameter :: eVToH = 1.0 / 27.2113962d0 !< Conversion factor from electron volt to Hartree (a.u. of energy)
   real(kind=DefReal), parameter :: Pi = 3.141592653589793
   real(kind=DefReal), parameter :: d2 = 2.0 !< Floating point 2.0
   real(kind=DefReal), parameter :: dp5 = 0.5 !< Floating point 0.5
   real(kind=DefReal), parameter :: d1 = 1.0 !< Floating point one
   real(kind=DefReal), parameter :: d0 = 0.0 !< Floating point zero
   real(kind=DefReal), parameter :: FPZero = 1.0d-10 !< Floating point zero; an approximation to machine epsilon
   real(kind=DefReal), parameter :: CMToAu = 0.000004556 !< Conversion factor from wavenumber to Hartree (a.u. of energy)
   real(kind=DefReal), parameter :: MassToAu = 1822.887 !< Conversion factor from atomic mass unit to mass of electron (au of mass)
   real(kind=DefReal), parameter :: BoltzK = 3.16681520371153d-6 !< Conversion factor from Kelvin to Hartree (a.u.)
   real(kind=DefReal), parameter :: DegToAu = BoltzK !< Conversion factor from degrees Celsius to Hartree
   integer(kind=DefInt), parameter :: i0 = 0 !< Integer zero
   integer(kind=DefInt4), parameter :: i4zero = 0 !< 4 byte Integer zero (used to call fms_ranb)
   complex(kind=DefComp), parameter :: c1i = (0.0, 1.0) !< Complex imaginary unit
   integer(kind=DefInt), parameter :: MaxParticles = 1000 !< Maximum number of particles
! Model integers
   integer, parameter :: FMSZERO = 0, &
                         TOY = 0, &
                         TEMPLATE = 1, &
                         QUANTICS = 13, &
                         TC = 14
!------------------------------------------------------------
!     FILE IDENTIFIERS
!------------------------------------------------------------
   integer(kind=DefInt) :: fmiOut = -1
   character(len=256) :: FMSWorkingDir = './'

!----------------------------------------------------------------
!     TIMERS
!----------------------------------------------------------------
   real(kind=DefReal) :: glRestartTime
   real(kind=DefReal) :: ftime = 0.0d0
   real(kind=DefReal) :: btime = 0.0d0
   real(kind=DefReal) :: cgtime = 0.0d0
   real(kind=DefReal) :: invtime = 0.0d0
   real(kind=DEfReal) :: bhstime = 0.0d0
   real(kind=DefReal) :: olaptime = 0.0d0
   real(kind=DefReal) :: pottime = 0.0d0
   real(kind=DefReal) :: kintime = 0.0d0
   real(kind=DefReal) :: estime = 0.0d0
   real(kind=DefReal) :: timei = 0.0d0
   real(kind=DefReal) :: timef = 0.0d0
   real(kind=DefReal) :: sdottime = 0.0d0
   integer :: glNESCalls = 0
   integer(kind=DefInt) :: glnStepsRejected = 0

!---------------------------------------------------------------------------
!                           NAMELIST PARAMETERS
!---------------------------------------------------------------------------
!>
!! \var integer (kind=DefInt) :: GlIModel
!! Program for calculating electronic structure
!! -  0 - Fake electronic structure
!! - 14 - TeraChem
!<
   integer(kind=DefInt) :: GlIModel

!> \var integer (kind=DefInt) :: GlIMethod
!! Method for calulating electronic structure.
!<
   integer(kind=DefInt) :: GlIMethod = 1

   character(len=2) :: GlCIntegType !< Selects the type of integrator.
   real(kind=DefReal) :: GlDEShift !< Global potential energy shift
   real(kind=DefReal) :: GlDTripletShift !< Global (additional) triplet energy shift
   real(kind=DefReal) :: GlDRegThresh !< Regularization Threshold.
   logical :: GlZAssumeHermitian !< if TRUE the code does not check for Hermiticity of H or S.
   logical :: GlZAdaptive !< If false, overrides adaptive integrators
   integer(kind=DefInt) :: GlISaddle !< Type of saddle point approximation to use in integral evaluation
   logical :: glzMinSearch
   real(kind=DefReal) :: GlMaxCoup !Max couplings
   logical :: glzIterInv
   integer(DefInt) :: NumInitBasis
   real(kind=DefReal) :: gldCurrentTStep
   real(kind=DefReal) :: gldTimeStep
   real(kind=DefReal) :: gldCoupTimeStep
   real(kind=DefReal) :: gldMinTimeStep
   logical :: glzFullyCoupled
   logical :: glzStochastic !< Use Stochastic Collapse algorithm?
!bfec
   logical :: glzStoOlap !< Use Overlap criterion for SS
   logical :: glzRejectAllStateFlip
   logical :: glzCentroids
!bfec
!-----------------------------------------------
!
!     Global parameters controlling the Izmaylov 2-D 2-state
!     toy model. Used in FMS_ToyModel, where their purpose is
!     explained in detail.
!
   real(kind=DefReal) :: glIzmOmegax
   real(kind=DefReal) :: glIzmOmegay
   real(kind=DefReal) :: glIzmXshift
   real(kind=DefReal) :: glIzmYshift
   real(kind=DefReal) :: glIzmDeltaE
   real(kind=DefReal) :: glIzmCoupC
!     Global parameter for the GAIMS model in ToyModelModule
   real(kind=DefReal) :: glGrsigma
   logical :: glzSPA_1Dmodel
!
!-----------------------------------------------
   logical :: glzStoSwiss
   logical :: glzStoStateSpecific
   real(kind=DefReal) :: gldStochaThresh
   real(kind=DefReal) :: gldLastSpawnSto
   real(kind=DefReal) :: gldMaxEDiff
   real(kind=DefReal) :: gldOLapThresh
   logical :: glzConstrain !< Use constraints for propagation?
   logical :: glzAnalysisMode
   real(kind=DefReal) :: gldSimulationTime
   logical :: gldNoisyGuess !< Add noise to the initial guess?
   real(kind=DefReal) :: gldNoisyGuessFac !< Factor by which noise is multiplied before being added to orbital guess
   real(kind=DefReal) :: gldNormThresh !< Threshold for norm convergence for VV propagation
   real(kind=DefReal) :: gldNGradStep !< Step size for numerical gradient
   logical :: glzCentNGrad !< Central numerical gradient? (if false, just forward)
   integer(kind=DefInt) :: glIgnoreState !< Ignore energy conservation problems on this state
   logical :: glzAvH !< Propagate classical trajectories on averaged Hamiltonian?
   integer(kind=DefInt) :: gliAvHNStates !< Number of states to average over?
   integer(kind=DefInt) :: gliAvHStates(100) !< List of states to average Hamiltonian over (defaults to all states)
   real(kind=DefReal) :: gldDecoherenceTime !< Kill useless trajectories after this amount of time
   integer(kind=DefInt) :: gliForceKill(50) !< User option to force-kill trajectories

   ! A unit for an open 'Control.dat',
   ! This is here because we need to pass it from ReadNameList.f to ParticleTypes.f.
   ! TODO(DH): Make this better!
   integer(kind=DefInt) :: ICUnit

!bfec
   integer(kind=DefInt) :: glirestTC !< Restart variable for Terachem
! xf added
   logical :: glzxfaims
   logical :: gldIgnoreStateAferField
   real(kind=DefReal) :: f0_xf
   real(kind=DefReal) :: polx_xf
   real(kind=DefReal) :: poly_xf
   real(kind=DefReal) :: polz_xf
   real(kind=DefReal) :: freq_xf
   real(kind=DefReal) :: t0_xf
   real(kind=DefReal) :: sigma_xf
   real(kind=DefReal) :: CEP_xf
   logical :: onespawnonly_xf
   real(kind=DefReal) :: sp_spwn_i, sp_spwn_f, gldCoupFieldTimeStep
! end xf added
! GAIMS added
   integer(kind=DefInt) :: NSing
   integer(kind=DefInt) :: NTrip
! GAIMS added end

!---------------------------------------------------------------------------
!                      OUTPUT CONTROL
!---------------------------------------------------------------------------
   integer, private, parameter :: MaxGeometry = 100
   integer(kind=DefInt) :: fmNStepToPrint !print output every n steps
   logical :: fmzPrintDipole
   logical :: fmzTrajFile
   logical :: fmzPotEnFile
   logical :: fmzNDatFile
   logical :: fmzEdatFile
   logical :: fmzCIVecFile
   logical :: fmzBundMatFile
   logical :: fmzCoupFile
   logical :: fmzChargeFile
   logical :: fmzForce
   logical :: fmzCorrFile
   logical :: fmzTDipoleFile
   logical :: fmzDipoleFile
   logical :: fmzQMRRFile
   logical :: fmzMMFile
   logical :: fmzAmpFile
   logical :: fmzAllText
   logical :: fmzXYZ
   logical :: fmzDCD
   logical :: fmzWriteEveryStep
   logical :: fmzPCOlap ! print the overlap (abs, Re, Im) between parent and child
   logical :: fmzSOME ! print the SOME at the centroid position
   logical :: fmzSOCeff ! print the effective SOC
   integer(kind=DefInt) :: fmIBond(2, MaxGeometry), fmNBonds
   integer(kind=DefInt) :: fmIAngle(3, MaxGeometry), fmNAngles
   integer(kind=DefInt) :: fmIDihedral(4, MaxGeometry), fmNDihedrals
   integer(kind=DefInt) :: fmIPyram(4, MaxGeometry), fmNPyrams

!>
!!    Holds open unit files
!<
   type File_Units
      integer(kind=DefInt) :: ICUnit
      integer(kind=DefInt) :: Unit_E
      integer(kind=DefInt) :: Unit_N
      integer(kind=DefInt) :: Unit_Spawn
      integer(kind=DefInt) :: Unit_Traj
      integer(kind=DefInt) :: Unit_Restart
   end type File_Units
!>
!!    Holds data for a grid to plot the wavefunction
!<
   type GridData
      real(DefReal) :: InX !< initial value of X
      real(DefReal) :: InY !< initial value of Y
      real(DefReal) :: InZ !< initial value of Z
      real(DefReal) :: FnX !< final value of X
      real(DefReal) :: FnY !< final value of Y
      real(DefReal) :: FnZ !< final value of Z
      integer(kind=DefInt) :: NXpoints !< Number of grid points in the X direction
      integer(kind=DefInt) :: NYpoints !< Number of grid points in the Y direction
      integer(kind=DefInt) :: NZpoints !< Number of grid points in the Z direction
      real(DefReal) :: DelX !< increment in the X direction
      real(DefReal) :: DelY !< increment in the X direction
      real(DefReal) :: DelZ !< increment in the X direction
   end type GridData

!---------------------------------------------------------------------------
!     ERROR THRESHOLDS
!---------------------------------------------------------------------------
   logical, private :: glzRejectStep = .false.
   logical :: glzTimestepRejection !< Reject timesteps and repropagate?
   logical :: glzDieOnMinStep
   real(kind=DefReal) :: gldNormCons
   real(kind=DefReal) :: gldNormStepCons
   real(kind=DefReal) :: gldEnergyCons
   real(kind=DefReal) :: gldEnergyStepCons
   real(kind=DefReal) :: gldTStepThresh

!     This is needed so that FMS_DieError does not depend on TeraChem module
   abstract interface
      subroutine terachem_finalizer()
         implicit none
      end subroutine terachem_finalizer
   end interface
   procedure(terachem_finalizer), pointer :: tc_finalize

   ! The default FMS_DieError handler prints the error message and stops the program.
   ! In order to be able to test error conditions in unit tests,
   ! we must use a different error handler. Thus we use
   ! a procedure pointer which can be swapped using the
   ! 'set_error_handler' subroutine defined in this module
   ! and called in unit_tests/main.F90 before running the unit tests.
   ! The unit test error handler is defined in unit_tests/testutils.F90.
   !
   ! One wrinkle here is that we can no longer assume that
   ! FMS_DieError stops the program. Thus, each call to FMS_DieError
   ! should be followed by a return statement.
   abstract interface
      subroutine error_handler(message)
         implicit none
         character(len=*), intent(in) :: message
      end subroutine error_handler
   end interface
   procedure(error_handler), pointer :: FMS_DieError => default_error_handler

contains

!-------------------------------------------------------------
! Global use routines
!-------------------------------------------------------------

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function to_lowercase(string) result(string_out)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      character(len=*), intent(in) :: string
      character(len(string)) :: string_out

      integer :: iA, iZ, i, n, shift

      iA = ichar('A')
      iZ = ichar('Z')
      shift = ichar('a') - ichar('A')

      string_out = string

      do i = 1, len(string_out)
         n = ichar(string_out(i:i))
         if (iA <= n .and. n <= iZ) string_out(i:i) = achar(n + shift)
      end do

   end function to_lowercase

!>
!!    Print Message and return
!<
   subroutine FMS_PrintMessg(Message, INum)

      character(len=*), intent(in) :: Message
      integer(kind=DefInt), optional, intent(in) :: INum

      if (present(INum)) then
         write (fmiOut, *) trim(Message), INum
      else
         write (fmiOut, *) trim(Message)
      end if
      flush (fmiOut)
   end subroutine FMS_PrintMessg

!>
!!    Delete a file
!<
   subroutine FMS_DeleteFile(FileName)
      character(len=*), intent(in) :: FileName
      integer(kind=DefInt) :: IUnit
      character(len=256) :: FilePath
      logical :: zExist

!     See if the file even exists. If not, just return.
      FilePath = trim(FMSWorkingDir)//trim(FileName)

      inquire (file=FilePath, exist=zExist)
      if (.not. zExist) return

      open (newunit=IUnit, file=FilePath, status='old')
      close (IUnit, status='delete')
   end subroutine FMS_DeleteFile

!>
!!    Build a filename of the form cfname.idno
!<
   function FMS_NumberedFileName(cfname, idno, realno) result(FileName)
      character(len=:), allocatable :: FileName
      integer(kind=DefInt), intent(in) :: idno
      character(len=*), intent(in) :: cfname
      real(kind=DefReal), intent(in), optional :: RealNo
      character(len=30) :: ctemp
      integer :: iStart

      if (.not. present(RealNo)) then
         write (ctemp, *) idno
      else
1000     format(f8.2)
         write (ctemp, 1000) realno
      end if

      iStart = verify(cTemp, ' ')
      FileName = trim(cfname)//'.'//trim(ctemp(iStart:))

   end function FMS_NumberedFileName

!>
!!    Build a filename of the form cfnameidno
!<
   function FMS_NumberedName(cfname, idno, realno) result(FileName)
      integer(kind=DefInt), intent(in) :: idno
      character(len=*), intent(in) :: cfname
      real(kind=DefReal), optional, intent(in) :: RealNo
      character(len=:), allocatable :: FileName
      character(len=30) :: ctemp
      integer :: iStart

      if (.not. present(RealNo)) then
         write (ctemp, *) idno
      else
1000     format(f8.2)
         write (ctemp, 1000) realno
      end if

      iStart = verify(cTemp, ' ')
      FileName = trim(cfname)//trim(ctemp(iStart:))
   end function FMS_NumberedName

!>
!!    Flags a step as rejected or not depending on SetVal
!!
!!    If user option DieOnMinStep is false, this subroutine also checks
!!    to see if the current (global) timestep is at a minimum.  If it is
!!    then step is never rejected.
!!    \param SetVal   Logical to flag step as rejected or not
!<
   subroutine FMS_RejectStep(SetVal)
      logical, intent(IN) :: SetVal

      if (.not. SetVal) then
         glzRejectStep = SetVal
         return
      end if

!     SetVal=.true. below
      if (FMS_StepRejected()) then !step was already rejected
         return
      end if

!     We are trying to change reject status to true below
      if (.not. glzDieOnMinStep .and. (gldCurrentTStep + FPZero) < gldMinTimeStep * 2.0d0) then
         write (fmiOut, *) 'Warning: Step not rejected as timestep at min'
         glzRejectStep = .false.
      else
         glzRejectStep = SetVal
      end if
   end subroutine FMS_RejectStep

   function FMS_StepRejected()
      logical :: FMS_StepRejected
      if (.not. glzTimestepRejection) then
         FMS_StepRejected = .false.
      else
         FMS_StepRejected = glzRejectStep
      end if
   end function FMS_StepRejected

!>
!!    Prints a complex matrix.
!!    \param A matrix to print
!!    \param Unit Fortran output unit
!!    @ingroup output
!<
   subroutine FMS_PrintCMat(A, iUnit)
      complex(kind=DefComp), intent(in) :: A(:, :)
      integer(kind=DefInt), intent(in) :: iUnit
      integer(kind=DefInt) :: i, j

      do i = 1, size(A, dim=1)
         write (iunit, '(10ES13.4E3)') (A(i, j), j=1, size(A, dim=2))
      end do
   end subroutine FMS_PrintCMat

!>
!!    Open output file for writing
!<
   subroutine FMS_OpenFile(fileName, iunit, fileExists)

      character(len=*), intent(in) :: fileName
      integer(kind=DefInt), intent(out) :: iunit
      logical, intent(out) :: fileExists
      character(len=:), allocatable :: filePath

      filePath = trim(FMSWorkingDir)//trim(fileName)
      inquire (file=filePath, exist=fileExists)
      if (fileExists) then
         open (newunit=iunit, file=filePath, position="append")
      else
         open (newunit=iunit, file=filePath)
      end if

   end subroutine FMS_OpenFile

!>
!!    Convert a word to lower case
!!    Works only for ASCII chars.
!<
   function lower_case(word) result(l_word)
      character(len=*), intent(in) :: word
      character(len=len(word)) :: l_word
      integer :: i, ic

      l_word = word

      do i = 1, len(l_word)
         ic = ichar(l_word(i:i))
         if (65 <= ic .and. ic < 90) l_word(i:i) = char(ic + 32)
      end do
   end function lower_case

!>
!!    Print error message to FMS.out and stdout and halt execution
!!    (This is FMS_DieError in disguise)
!<
   subroutine default_error_handler(message)
      use, intrinsic :: iso_fortran_env, only: error_unit
      character(len=*), intent(in) :: message

      write (error_unit, '(a)') 'ERROR: '//trim(message)

      ! Don't write to FMS.out if it hasn't been opened yet!
      if (fmiOut /= -1) then
         write (fmiOut, '(a)') 'ERROR: '//trim(message)
         flush (fmiOut)
      end if
      if (gliModel == TC .and. associated(tc_finalize)) then
         call tc_finalize()
      end if
      stop 1
   end subroutine default_error_handler

   subroutine set_terachem_finalizer(method)
      procedure(terachem_finalizer) :: method
      tc_finalize => method
   end subroutine set_terachem_finalizer

   ! Swap the error handler, this is used in unit tests
   subroutine set_error_handler(method)
      procedure(error_handler) :: method
      FMS_DieError => method
   end subroutine set_error_handler

end module GlobalModule
