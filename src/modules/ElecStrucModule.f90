!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    @brief Stores parameters pertaining to the electronic wavefunction
!!
!!    ElecStrucModule holds variables describing the electronic wavefunction.
!!    Specifically, it contains a set of flags specifying the ESP's
!!    computational capabilities, as well as data specifying the array
!!    parameters for storing the electronic wavefunction.
!!
!!    The actual arrays storing electronic structure data are allocated
!!    within the T_Trajectory%ElecStruc data structure; these arrays
!!    will not be allocated the array lengths in this module are 0.
!!
!!    @note Module code: es
!<
module ElecStrucModule
   use GlobalModule, only: DefReal, DefInt, fmiOut
   implicit none
   public

!     What can this electronic structure method do?

!> True if electronic structure package has analytic forces available
   logical :: eszAnalyticForces

!> True if electronic structure package can return analytic non-adiabatic
!! coupling vector
   logical :: eszNACoupVec

!> True if electronic structure package can return atomic partial charges
   logical :: eszPartialCharges

!> True if electronic structure package has analytic forces available
   logical :: eszTransDipole

!> True if electronic structure package has analytic forces available
   logical :: eszDipoleMoment

!> True if electronic structure package has quadrupole moments available
   logical :: eszQuadpoleMoment

!> True if electronic structure package has analytic forces available
   logical :: eszMMForce

!
!     Wavefunction storage specifications
!

!> Number of orbitals in electronic wavefunction
   integer(kind=DefInt) :: esNOrbs

!> Number of basis AO primitives in electronic wavefunction
   integer(kind=DefInt) :: esNBasis

!> Number of active space orbitals
   integer(kind=DefInt) :: esNCasorbs

!> Number of active space electrons
   integer(kind=DefInt) :: esNCaselec

!> Length of CI vector
   integer(kind=DefInt) :: esLCivec

!> Size of the electronic phase vector
   integer(kind=DefInt) :: esnElecPhase

!> Size of the electronic structure blob
   integer(kind=DefInt) :: esBlobSize

!TM
!
! Flags for dummy atoms
!
!> Number of dummy atoms
   integer(kind=DefInt) :: esNDummy, esNPart, esNMM

!> Coefficients
   real(kind=DefReal), allocatable :: esDummyWeight(:, :)

contains

!>
!!    Initializes the current electronic structure package. Also,
!!    grabs esLCiVec and esnBasis (the size of the CI vector and
!!    the number of basis functions, respectively) from the electronic
!!    structure package.
!!
!!    This routine also lets FMS know what abilities the electronic
!!    structure package has by setting the following flags:
!!
!!     logical :: eszAnalyticForces  - analytic gradients
!!     logical :: eszNACoupVec       - non-adiabatic couplings
!!     logical :: eszPartialCharges  - partial charge analysis
!!     logical :: eszTransDipole     - transition dipole moments
!!     logical :: eszDipoleMoment    - regular dipole moments
!!    \param NumParticles Number of particles in system
!!    \param NumStates Number of electronic states
!!    @ingroup ESP
!<
   subroutine FMS_ESInit(NumParticles, NumStates)
      use GlobalModule, only: gliModel, gliMethod, &
                              FMSZERO, QUANTICS, TC, TEMPLATE, FMS_DieError

      integer(kind=DefInt), intent(in) :: NumParticles
      integer(kind=DefInt), intent(inout) :: NumStates
      character(len=50) :: errmsg

!
!     Set defaults
!
      eszAnalyticForces = .true.
      eszNACoupVec = .true.
      eszPartialCharges = .false.
      eszTransDipole = .false.
      eszDipoleMoment = .false.
      eszMMForce = .false.
      esnOrbs = 0
      esnBasis = 0
      esnCasOrbs = 0
      esnCasElec = 0
      eslCIVec = 0
      esnElecPhase = 0

      select case (gliModel)

!
!     FMSZero: fake electronic structure
!
      case (FMSZERO)
         write (fmiOut, *) 'FMSZero: Using toy potentials in adiabatic basis.'
         eszAnalyticForces = .true.
         eszNACoupVec = .true.
         eszPartialCharges = .false.
         eszTransDipole = .false.
         eszDipoleMoment = .false.
         select case (gliMethod)

         case (1)
            call FMS_DieError('Tully model 1 currently not implemented')

         case (2)
            write (fmiOut, *) 'gliMethod=2: Persico Model 2D 2-state CI'
            if (NumStates /= 2) then
               write (fmiOut, *) 'Number of states set to 2.'
               NumStates = 2
            end if
            if (NumParticles < 2) then
               call FMS_DieError('Persico requires >=2 particles.')
            end if

         case (3)
            write (fmiOut, *) 'gliMethod=3: Izmaylov Model 2D 2-state CI'
            if (NumStates /= 2) then
               write (fmiOut, *) 'Number of states set to 2.'
               NumStates = 2
            end if
            if (NumParticles < 2) then
               call FMS_DieError('Izmaylov model requires 2 particles.')
            end if

         case (4) !GAIMS_model
            write (fmiOut, *) 'gliMethod=4: GAIMS Model 1D 2-state SOC'
            if (NumStates /= 2) then
               write (fmiOut, *) 'Number of states set to 2.'
               NumStates = 2
            end if

         case (5) !XFAIMS CH3I_model
            write (fmiOut, *) 'gliMethod=5: XFAIMS CH3I Model 1D 2-state'
            if (NumStates /= 2) then
               write (fmiOut, *) 'Number of states set to 2.'
               write (fmiOut, *) 'Extension to 3 states not available.'
               NumStates = 2
            end if

            if (NumParticles /= 1) then
               call FMS_DieError('CH3I model requires 1 particle.')
            end if

            ! allowing dipoles and TDMs
            eszTransDipole = .true.
            eszDipoleMoment = .true.

         case default
            write (errmsg, '(a, i0)') 'Invalid iMethod value: ', gliMethod
            call FMS_DieError(errmsg)
         end select
!     end of fmszero initialization

!
!     FMSTemplate: write electronic structure input decks with templates
!
      case (TEMPLATE)
         write (fmiOut, *) 'FMSTemplate: system call ES interface.'
         if (eslCIVec > 0) esnElecPhase = NumStates

      case (TC)
         eszAnalyticForces = .true.
         eszNACoupVec = .true.
         esnOrbs = esNBasis
         esNElecPhase = NumStates
         eszTransDipole = .true.
         eszDipoleMoment = .true.
         eszPartialCharges = .true.

      case (QUANTICS)
         eszAnalyticForces = .true.
         eszNACoupVec = .true.
         esnOrbs = esNBasis
         esNElecPhase = NumStates
         eszTransDipole = .false.
         eszDipoleMoment = .false.
         eszPartialCharges = .false.

      case default
         call FMS_DieError('Invalid iModel')

      end select

   end subroutine FMS_ESInit
end module ElecStrucModule
