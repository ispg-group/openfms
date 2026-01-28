!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Trajectory class specification and related methods
!!
!! This module specifies the trajectory type, and all related
!! methods. A given instance of a trajectory object contains the
!! classical position and momenta for a given trajectory basis function
!! at a given instant in time. Also included in the trajectory data
!! structure are electronic structure quantities, in the "T_ElecStruc"
!! datatype, and flags describing those quantities in "T_ESFlags".
!!
!! This module also defines several other derived types, some of which
!! are stored in T_Trajectory: T_ElecStruc, T_ESFlags, T_BFlags
!!
!! Methods in this module pertain to basic manipulation of T_Trajectory
!! derived type. Methods that pertain to electronic structure are
!! in the TrajectoryCalcsModule, while the methods related to file IO
!! are in TrajectoryIOModule.
!!
!! \image latex "../sources/TrajModule.png"
!<
module TrajectoryModule

   use GlobalModule, only: D2, Pi, DefInt, DefReal, DefComp, fmiOut, FPZero, &
                           glzStoSwiss, gldDecoherenceTime, FMS_DieError
   use ParticleModule
   use QM_MM_Module, only: qczQMMM
   implicit none
   private
   public :: t_Trajectory, t_ElecStruc
   public :: assign_trajectory, assignment(=)
   public :: FMS_Distance, FMS_Angle, FMS_Dihedral, FMS_CalcPyrAngle
   public :: FMS_D_Distance, FMS_D_Angle, FMS_D_Dihedral

!--------------------------------------------------------------------!
!                 DERIVED TYPES FOR TRAJECTORIES                     !
!--------------------------------------------------------------------!

!--------------------------------------------------------------------!
! Flags to track electronic structure status
   type T_ESFlags
!        private
      logical :: zIgnoreErrors, & !< For the first timestep of centroids only - we can ignore phasing and diabatization errors.
                 zESExists, & !< Does a wavefunction exist to restart from?
                 ZPotEnCurrent, & !< Is the potential calculated for the current geometry
                 ZTransDipsCurrent, & !< If transition dipoles are calculated for this geometry
                 ZTransDipsCurrentxf, & !< xf added
                 ZDipolesCurrent, & !< If dipoles are calculated for this geometry
                 ZQuadpolesCurrent, & !< If quadrupoles are calculated for this geometry
                 zMMPotCurrent, & !< If MM potential is calculated for the current geometry
                 zMMForceCurrent, & !< If MM forces are calculated for the current geometry
                 zModPotCurrent !< Are external force modifications current?

      !< Is the coupling current with the geometry for this state?
      !< Are the derivatives current with the geometry for this state?
      logical, allocatable :: ZDerivCurrent(:, :)
      logical, allocatable :: ZSOMCurrent(:, :, :, :)
      !< This should be replaced by zForcesCurrent
   end type T_ESFlags

!---------------------------------------------------------------------!
!>
!!    Flags to track bundle status
!<
   type T_BFlags
!     private
      logical :: zBundleCurrent, & !< Has the trajectory changed since the last bundle calculation?
                 ZAmpDotCurrent !< Is AmpDot current?

   end type T_BFlags

!----------------------------------------------------------------------------
!>
!!     Electronic structure data for a trajectory
!!     TODO: Should this type live in ElecStrucModule?
!<
   type T_ElecStruc
!        private
      real(kind=DefReal) :: ModPot, & !< External potential, for SMD, etc.
                            MMPot !< Potential of MM system, for QM/MM runs

      real(kind=DefReal), allocatable :: ElecPhase(:), & !< Electronic phase
                                         PotEn(:), & !< Potential Energy for given nuclear configuration.
                                         OldPotEn(:), & !< Potential Energy from previous timestep.
                                         Dipole(:, :), & !< Dipole of each state
                                         QMRR(:), & !< Quadruple moment
                                         DerivMat(:, :, :), & !< holds derivative matrix (diagonals=-force, off-diagonals=coupling)
                                         ModForce(:), & !< External forces, for SMD, etc.
                                         MMForce(:), & !< Forces from MM sytem, for QM/MM runs
                                         OldOrbitals(:, :), & !< Used to ensure continuity of phase
                                         OldCIVecs(:, :), & !< The CI vectors obtained at the last time  step.  We will also use these to ensure continuity of
                                         OverlapMatrix(:, :), & !< Overlap matrix between old and new wavefunctions
                                         OldBlob(:), & !< Sent back and forth
                                         TransDipole(:, :), & !< Transition dipole between various states at current geom
                                         TransDipolexf(:), & !< xf added
                                         !           MSPT2S(:,:),  &      !< CASPT2 overlap matrix
                                         !           MSPT2C(:,:),  &      !< MSPT2 Mixing Coefficient
                                         !           OldMSPT2S(:,:),   &  !< CASPT2 overlap last timestep
                                         !           Need this one for restart file
                                         OldMSPT2C(:, :) !< MSPT2 Mixing Coeff last timestep

      complex(kind=DefComp), allocatable :: SOMat(:, :, :, :) !< Spin Orbit Matrix (i,j,Ms)

   end type T_ElecStruc

!-------------------------------------------------------------------------------
!>
!!    Spawning with informed stochastic selections (SWISS)
!!    data for a trajectory
!<
   type T_SWISS
      real(kind=DefReal) :: SelectionTime, & !< Time at which selection time it to be performed
                            ParentOverlap, & !< Overlap between trajecory and its parent
                            BirthDate !< Time at which trajectory was created
   end type T_SWISS

!---------------------------------------------------------------------------
!>
!!    This data structure holds an individual trajectory time step
!<
   type T_Trajectory

      integer(kind=DefInt) :: NumDimensions, & !< Total dimensionality
                              NumParticles, & !< Number of particles
                              NumStates, & !< Number of electronic states
                              StateID, & !< Electronic state index
                              TrajID, & !< Trajectory identifier
                              ParentID, & !< ID of parent
                              CentID(2) !< Trajectories that this is centroid for

! New stuff for triplet states
      integer(kind=DefInt) :: NumSing, & !< Number of electronic singlet states
                              NumTrip, & !< Number of electronic triplet states
                              CBF, & !< Contracted trajectory identifier
                              Ms !< Ms value of triplet trajectory (-1,0,1)

      real(kind=DefReal), private :: CurrentTime !< Time

      real(DefReal) :: DeadTime !< Time that the trajectory became uncoupled
      !< and is marked for death
      logical :: zCent !< .true. for centroids, false otherwise

! logical for triplet states
      logical :: triplet !< .true. when trajectory is in

      type(T_ElecStruc) :: ElecStruc !< electronic structure information
      type(T_SWISS) :: SWISS !< electronic structure information
      type(T_ESFlags) :: ESFlags !< Flags for electronic structure
      type(T_BFlags) :: BFlags !< Flags for bundle quantities

      !< Particle data
      type(T_Particle), allocatable :: Particle(:)

      ! Phase propagation
      complex(kind=DefComp) :: Amplitude !< Trajectory amplitude
      real(kind=DefReal) :: Phase !< Semiclassical phase

      ! Spawning quantities
      real(kind=DefReal) :: Pop !< The Mulliken population

      integer(kind=DefInt), allocatable :: SpawnMode(:) !< Spawning mode with respect to a given electronic state

      real(kind=DefReal), allocatable :: CoupHist(:, :), & !< History of Coupling amplitudes.  Used for spawning criterion.
                                         LastSpawn(:), & !< Point in time at which we exited the last nonadiabatic coupling region
                                         SpawnTime(:) !< Time of last spawn

      integer(kind=DefInt), allocatable :: SpawnCoupled(:) !< this keeps track of which Trajectory this is coupled to
      !< Size is nstates.
      !< If SpawnCoupled(n)==0 the traj is not involed in
      !< spawning with state n and may spawn to that state
      !< If SpawnCoupled(n)==K the trajectory is coupled with trajID K
      !< as a result of spawning and no further spawns are allowed to
      !< this state.  SpawnCoupled is set in SpawnModule as a result of
      !< a successful spawn.  It is reset to 0 is when the overlaps
      !< have dropped beow threshold
   contains
      private
      generic, public :: get_mass => get_mass_comp
      generic, public :: get_mass => get_mass_vec
      generic, public :: get_width => get_width_comp
      generic, public :: get_width => get_width_vec
      generic, public :: get_mom => get_mom_comp
      generic, public :: get_mom => get_mom_vec
      generic, public :: get_mom => get_mom_all
      generic, public :: set_mom => set_mom_comp
      generic, public :: set_mom => set_mom_vec
      generic, public :: set_mom => set_mom_all
      generic, public :: get_vel => get_vel_all
      generic, public :: get_vel2 => get_vel2_all
      generic, public :: set_vel => set_vel_all
      generic, public :: set_vel2 => set_vel2_all
      generic, public :: get_mom2 => get_mom2_vec
      generic, public :: set_mom2 => set_mom2_vec
      generic, public :: get_mom2 => get_mom2_all
      generic, public :: set_mom2 => set_mom2_all
      generic, public :: get_pos => get_pos_comp
      generic, public :: get_pos => get_pos_vec
      generic, public :: get_pos => get_pos_all
      generic, public :: set_pos => set_pos_comp
      generic, public :: set_pos => set_pos_vec
      generic, public :: set_pos => set_pos_all
      procedure, public :: create => create_trajectory
      procedure, public :: destroy => destroy_trajectory
      procedure, public :: copy_from
      procedure, public :: get_time, set_time
      procedure, public :: print_id
      procedure, public :: get_pop, rescale_phases
      procedure, public :: is_dead
      procedure, public :: center_of_mass
      procedure, public :: geom_changed
      procedure :: get_mass_comp, get_mass_vec
      procedure :: get_width_comp, get_width_vec
      procedure :: get_pos_vec, set_pos_vec
      procedure :: get_pos_all, set_pos_all
      procedure :: get_mom_comp, set_mom_comp
      procedure :: get_mom_vec, set_mom_vec
      procedure :: get_mom_all, set_mom_all
      procedure :: get_vel_all, set_vel_all
      procedure :: get_vel2_all, set_vel2_all
      procedure :: get_pos_comp, set_pos_comp
      procedure :: get_mom2_all, set_mom2_all
      procedure :: get_mom2_vec, set_mom2_vec

   end type T_Trajectory

   interface assignment(=)
      module procedure assign_trajectory
   end interface

!     Function overloads:
   interface FMS_Distance
      module procedure FMS_Distance_Trajectory
   end interface

   interface FMS_D_Distance
      module procedure FMS_D_Distance_Trajectory
   end interface

   interface FMS_Angle
      module procedure FMS_Angle_Trajectory
   end interface

   interface FMS_D_Angle
      module procedure FMS_D_Angle_Trajectory
   end interface

   interface FMS_Dihedral
      module procedure FMS_Dihedral_Trajectory
   end interface

   interface FMS_D_Dihedral
      module procedure FMS_D_Dihedral_Trajectory
   end interface

contains

!     Print the trajectory's identifier
!>
!! Identify this trajectory in FMS.out
!<
   subroutine print_id(T1)
      class(T_Trajectory) :: T1

2000  format('    Centroid of trajectories ', i3, ' and ', i3, '. (Cent #', i3, ')')
2001  format('    Trajectory #', i3)

      if (T1%zCent) then
         write (fmiOut, 2000) T1%CentID(1), T1%CentID(2), T1%TrajID
      else
         write (fmiOut, 2001) T1%TrajID
      end if
   end subroutine print_id

!     Constructors/Destructors/Assignment
!>
!!    Memory allocation for creating a trajectory structure
!<
   subroutine create_trajectory(T1, NumParticles, NumStates)
      class(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: NumParticles, NumStates
      integer(kind=DefInt) :: n

      if (allocated(T1%Particle) .and. allocated(T1%LastSpawn)) then
         if (NumParticles == size(T1%Particle) .and. NumStates == size(T1%LastSpawn)) then
            return
         end if
      end if

      call T1%destroy()

      if (NumStates <= 0) then
         call FMS_DieError('CreateTrajectory: NumStates must be > 0')
         return
      end if
      if (NumParticles <= 0) then
         call FMS_DieError('CreateTrajectory: NumParticles must be > 0')
         return
      end if

!  Trajectory information initialization
      T1%TrajID = 0
      T1%CentID = -1
      T1%ParentID = 0
      T1%zCent = .false.

      T1%NumStates = NumStates
      T1%NumDimensions = 3 * NumParticles
      T1%NumParticles = NumParticles
      T1%StateID = 0

      call T1%set_time(0.d0)
      T1%Amplitude = 0.d0
      T1%Phase = 0.d0
      T1%DeadTime = 0.d0

      call FMS_CreateElectronicStructure(T1%ElecStruc, NumStates, NumParticles)

!   Flags creation and initialization
      allocate (T1%ESFlags%ZDerivCurrent(NumStates, NumStates))
      allocate (T1%ESFlags%ZSOMCurrent(NumStates, NumStates, 3, 3))
      call T1%geom_changed()
      T1%ESFlags%ZDerivCurrent = .false.
      T1%ESFlags%ZSOMCurrent = .false.
      T1%ESFlags%zESExists = .false.
      T1%ESFlags%zIgnoreErrors = .false.

!   Allocate parameter arrays
      allocate (T1%SpawnMode(NumStates), T1%CoupHist(3, NumStates), T1%SpawnTime(NumStates), &
                T1%LastSpawn(NumStates), T1%SpawnCoupled(NumStates))
      T1%CoupHist = 0.0d0
      T1%SpawnTime = 0.0d0
      T1%LastSpawn = 0.0d0
      T1%SpawnMode = 0
      T1%SpawnCoupled = 0

      if (glzStoSwiss) then
         T1%SWISS%SelectionTime = 0.d0
         T1%SWISS%ParentOverlap = 0.d0
         T1%SWISS%BirthDate = 0.d0
      end if

!     Create particles:
      allocate (T1%Particle(NumParticles))
!     TODO(DH): Why is numdim hardcoded here?
      do n = 1, NumParticles
         call T1%Particle(n)%create(id=n, numdim=3)
      end do

   end subroutine create_trajectory
!>
!!    Memory deallocation to destroy a trajectory structure
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine destroy_trajectory(T1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(inout) :: T1

      integer(kind=DefInt) :: n

      T1%NumStates = 0
      T1%NumSing = 0
      T1%NumTrip = 0
      T1%NumParticles = 0
      if (glzStoSwiss) then
         T1%SWISS%SelectionTime = 0.d0
         T1%SWISS%ParentOverlap = 0.d0
         T1%SWISS%BirthDate = 0.d0
      end if
      if (allocated(T1%Particle)) then
         do n = 1, size(T1%Particle)
            call T1%Particle(n)%destroy()
         end do
         deallocate (T1%Particle)
      end if
      if (allocated(T1%SpawnMode)) deallocate (T1%SpawnMode)
      if (allocated(T1%SpawnCoupled)) deallocate (T1%SpawnCoupled)
      if (allocated(T1%CoupHist)) deallocate (T1%CoupHist)
      if (allocated(T1%LastSpawn)) deallocate (T1%LastSpawn)
      if (allocated(T1%SpawnTime)) deallocate (T1%SpawnTime)

      if (allocated(T1%ESFlags%ZDerivCurrent)) deallocate (T1%ESFlags%ZDerivCurrent)
      if (allocated(T1%ESFlags%ZSOMCurrent)) deallocate (T1%ESFlags%ZSOMCurrent)

      call FMS_DestroyElectronicStructure(T1%ElecStruc)
   end subroutine destroy_trajectory

   subroutine FMS_CreateElectronicStructure(ES, NumStates, NumParticles)
      use ElecStrucModule, only: esBlobSize, esnBasis, esnelecphase, eslcivec
      type(T_ElecStruc), intent(inout) :: ES
      integer, intent(in) :: NumStates, NumParticles

      allocate (ES%PotEn(NumStates))
      allocate (ES%OldPotEn(NumStates))
      allocate (ES%DerivMat(NumStates, NumStates, NumParticles * 3))
      allocate (ES%Dipole(NumStates, 4))
      allocate (ES%TransDipole(NumStates, 3))
      allocate (ES%QMRR(NumStates))
      ES%PotEn = 0.0d0
      ES%OldPotEn = 0.0d0
      ES%DerivMat = 0.0d0
      ES%Dipole = 0.0d0
      ES%TransDipole = 0.0d0
      ES%QMRR = 0.0d0

      allocate (ES%SOMat(NumStates, NumStates, 3, 3))
      ES%SOMat = 0.0d0

! xf added
      allocate (ES%TransDipolexf(4))
      ES%TransDipolexf = 0.0d0
! xf added end

      if (qczQMMM) then
         allocate (ES%MMForce(NumParticles * 3))
         ES%MMForce = 0.0d0
      end if
      allocate (ES%ModForce(NumParticles * 3))
      ES%ModForce = 0.0d0

!     Initialize electronic structure storage
      if (esLCiVec > 0) then
         allocate (ES%OldCIVecs(NumStates, esLCiVec))
         ES%OldCIVecs = 0.0d0
      end if
      if (esnBasis > 0) then
         allocate (ES%OldOrbitals(esNBasis, esNBasis))
         ES%OldOrbitals = 0.0d0
      end if

      allocate (ES%OverlapMatrix(NumStates, NumStates))
      ES%OverlapMatrix = 0.0d0

      if (esBlobSize > 0) then
         allocate (ES%OldBlob(esBlobSize))
         ES%OldBlob = 0.0d0
      end if
      if (esnElecPhase > 0) then
         allocate (ES%Elecphase(esnElecPhase))
         ES%Elecphase = 1.0d0
      end if
      !allocate(ES%MSPT2S   (NumStates,NumStates), &
      !         ES%MSPT2C   (NumStates,NumStates), &
      !         ES%OldMSPT2S(NumStates,NumStates))
      !ES%MSPT2S    = 0.0d0
      !ES%MSPT2C    = 0.0d0
      !ES%OldMSPT2S = 0.0d0
      allocate (ES%OldMSPT2C(NumStates, NumStates))
      ES%OldMSPT2C = 0.0d0
   end subroutine FMS_CreateElectronicStructure

   subroutine FMS_DestroyElectronicStructure(ES)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_ElecStruc), intent(inout) :: ES
      if (allocated(ES%PotEn)) deallocate (ES%PotEn)
      if (allocated(ES%OldPotEn)) deallocate (ES%OldPotEn)
      if (allocated(ES%DerivMat)) deallocate (ES%DerivMat)
      if (allocated(ES%MMForce)) deallocate (ES%MMForce)
      if (allocated(ES%ModForce)) deallocate (ES%ModForce)
      if (allocated(ES%SOMat)) deallocate (ES%SOMat)

      if (allocated(ES%TransDipole)) deallocate (ES%TransDipole)
      if (allocated(ES%Dipole)) deallocate (ES%Dipole)
      if (allocated(ES%QMRR)) deallocate (ES%QMRR)

      if (allocated(ES%oldCIVecs)) deallocate (ES%oldCIVecs)
      if (allocated(ES%oldOrbitals)) deallocate (ES%oldOrbitals)
      if (allocated(ES%OverlapMatrix)) deallocate (ES%OverlapMatrix)
      if (allocated(ES%OldBlob)) deallocate (ES%OldBlob)
      if (allocated(ES%ElecPhase)) deallocate (ES%ElecPhase)
      if (allocated(ES%TransDipoleXF)) deallocate (ES%TransDipoleXF)

!     if (allocated(ES%MSPT2S))         deallocate(ES%MSPT2S)
!     if (allocated(ES%MSPT2C))         deallocate(ES%MSPT2C)
!     if (allocated(ES%OldMSPT2S))      deallocate(ES%OldMSPT2S)
      if (allocated(ES%OldMSPT2C)) deallocate (ES%OldMSPT2C)
   end subroutine FMS_DestroyElectronicStructure
!>
!!    Memory allocation and book-keeping for assignment operation T1=T2
!!    \param T1 Trajectory to be assigned to
!!    \param T2 Trajectory being assigned from
!<
   subroutine assign_trajectory(T1, T2)
      type(T_Trajectory), intent(inout) :: T1
      type(T_Trajectory), intent(in) :: T2
      integer(kind=DefInt) :: IState, IParticle
      integer(kind=DefInt) :: i

      if ((T1%NumParticles /= T2%NumParticles) .or. (T1%NumStates /= T2%NumStates)) then
         call T1%destroy()
         call T1%create(T2%NumParticles, T2%NumStates)
      end if
      T1%TrajID = T2%TrajID
      T1%CentID = T2%CentID
      T1%ParentID = T2%ParentID
      T1%zCent = T2%zCent
      T1%NumDimensions = T2%NumDimensions
      T1%StateID = T2%StateID
      T1%Amplitude = T2%Amplitude
      T1%Phase = T2%Phase
      T1%Pop = T2%Pop
      T1%DeadTime = T2%DeadTime

      call T1%set_time(T2%get_time())

      ! GAIMS added
      T1%CBF = T2%CBF
      T1%triplet = T2%triplet
      T1%Ms = T2%Ms
      T1%ElecStruc%SOMat = T2%ElecStruc%SOMat
      T1%ESFlags%zSOMCurrent = T2%ESFlags%zSOMCurrent
      T1%NumSing = T2%NumSing
      T1%NumTrip = T2%NumTrip
      ! GAIMS end added

      do IState = 1, T1%NumStates
         do i = 1, 3
            T1%ElecStruc%TransDipole(IState, i) = T2%ElecStruc%TransDipole(IState, i)
            T1%ElecStruc%Dipole(IState, i) = T2%ElecStruc%Dipole(IState, i)
         end do
         T1%ElecStruc%Dipole(IState, 4) = T2%ElecStruc%Dipole(IState, 4)
         T1%ElecStruc%QMRR(IState) = T2%ElecStruc%QMRR(IState)
      end do

      T1%SpawnMode = T2%SpawnMode
      T1%CoupHist = T2%CoupHist
      T1%LastSpawn = T2%LastSpawn
      T1%SpawnTime = T2%SpawnTime
      if (glzStoSwiss) then
         T1%SWISS%SelectionTime = T2%SWISS%SelectionTime
         T1%SWISS%ParentOverlap = T2%SWISS%ParentOverlap
         T1%SWISS%BirthDate = T2%SWISS%BirthDate
      end if

      do IParticle = 1, T1%NumParticles
         T1%Particle(IParticle) = T2%Particle(IParticle)
      end do

      T1%ElecStruc%PotEn = T2%ElecStruc%PotEn
      T1%ElecStruc%OldPotEn = T2%ElecStruc%OldPotEn
      T1%ElecStruc%DerivMat = T2%ElecStruc%DerivMat

      if (qczQMMM .and. allocated(T2%ElecStruc%MMForce)) then
         T1%ElecStruc%MMForce = T2%ElecStruc%MMForce
         T1%ElecStruc%MMPot = T2%ElecStruc%MMPot
      end if

      if (allocated(T2%ElecStruc%ModForce)) then
         T1%ElecStruc%ModForce = T2%ElecStruc%ModForce
      end if

      T1%ElecStruc%ModPot = T2%ElecStruc%ModPot

      if (allocated(T2%ElecStruc%OldCIVecs)) T1%ElecStruc%OldCIVecs = T2%ElecStruc%OldCIVecs
      if (allocated(T2%ElecStruc%OldOrbitals)) T1%ElecStruc%OldOrbitals = T2%ElecStruc%OldOrbitals
      if (allocated(T2%ElecStruc%OverlapMatrix)) T1%ElecStruc%OverlapMatrix = T2%ElecStruc%OverlapMatrix
      if (allocated(T2%ElecStruc%OldBlob)) T1%ElecStruc%OldBlob = T2%ElecStruc%OldBlob
      if (allocated(T2%ElecStruc%ElecPhase)) T1%ElecStruc%ElecPhase = T2%ElecStruc%ElecPhase

!     if (allocated(T2%ElecStruc%MSPT2C))        T1%ElecStruc%MSPT2C        = T2%ElecStruc%MSPT2C
      if (allocated(T2%ElecStruc%OldMSPT2C)) T1%ElecStruc%OldMSPT2C = T2%ElecStruc%OldMSPT2C
!     DH: Would probably also need this?
!     if (allocated(T2%ElecStruc%MSPT2S))        T1%ElecStruc%MSPT2S        = T2%ElecStruc%MSPT2S
!     if (allocated(T2%ElecStruc%OldMSPT2S))     T1%ElecStruc%OldMSPT2S     = T2%ElecStruc%OldMSPT2S

! Flags
      T1%ESFlags%zESExists = T2%ESFlags%zESExists
      T1%ESFlags%zPotEnCurrent = T2%ESFlags%zPotEnCurrent
      T1%ESFlags%zTransDipsCurrent = T2%ESFlags%zTransDipsCurrent
      T1%ESFlags%zDipolesCurrent = T2%ESFlags%zDipolesCurrent
      T1%ESFlags%zMMPotCurrent = T2%ESFlags%zMMPotCurrent
      T1%ESFlags%zMMForceCurrent = T2%ESFlags%zMMForceCurrent
      T1%ESFlags%zIgnoreErrors = T2%ESFlags%zIgnoreErrors
      T1%ESFlags%zDerivCurrent = T2%ESFlags%zDerivCurrent

! These flags are dependent upon other trajectories in the bundle
! so are set to false in general. Assign Bundle will overwrite them with
! the appropriate value if we are copying an entire bundle.
      T1%BFlags%zBundleCurrent = .false.
      T1%BFlags%zAmpDotCurrent = .false.

   end subroutine assign_trajectory

!>
!!    Copy Trajectory T2 to T1
!!    \param T1 Trajectory to be assigned to
!!    \param T2 Trajectory being assigned from
!<
   subroutine copy_from(T1, T2)
      class(T_Trajectory), intent(inout) :: T1
      class(T_Trajectory), intent(in) :: T2

      call assign_trajectory(T1, T2)
   end subroutine copy_from

   function get_time(T1) result(time)
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: time

      time = T1%CurrentTime
   end function get_time

   subroutine set_time(T1, time)
      class(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: time

      T1%CurrentTime = time
   end subroutine set_time

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_pos_comp(T1, n, i, NewPosition)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !>
      !! Change a component of a particle's position. Set all electronic structure
      !! and bundle matrix flags to false.
      !!
      !! Scope: Public
      !<
      class(T_Trajectory) :: T1
      integer(kind=DefInt), intent(in) :: n, i
      real(kind=DefReal), intent(in) :: NewPosition

      call T1%Particle(n)%set_pos(i, NewPosition)
      call T1%geom_changed()
   end subroutine set_pos_comp

   ! Change a component of a particle's momentum.
   subroutine set_mom_comp(T1, n, i, new_mom)
      class(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: n, i
      real(kind=DefReal), intent(in) :: new_mom

      call T1%Particle(n)%set_mom(i, new_mom)
   end subroutine set_mom_comp

   function get_mom_comp(T1, n, i) result(mom)
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: n, i
      real(kind=DefReal) :: mom

      mom = T1%Particle(n)%get_mom(i)
   end function get_mom_comp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_pos_comp(T1, n, i) result(pos)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      integer(kind=DefInt), intent(in) :: n, i
      real(kind=DefReal) :: pos

      pos = T1%Particle(n)%get_pos(i)
   end function get_pos_comp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_pos_vec(T1, n, NewPosition)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !>
      !! Change a particle's position. Set all electronic structure
      !! and bundle matrix flags to false.
      !!
      !! Scope: Public
      !<
      class(T_Trajectory) :: T1
      integer(kind=DefInt) :: n
      real(kind=DefReal) :: NewPosition(:)

      call T1%Particle(n)%set_pos(NewPosition)
      call T1%geom_changed()

   end subroutine set_pos_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_pos_vec(T1, n) result(vec)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: n
      real(kind=DefReal) :: vec(T1%Particle(n)%NumDimensions)

      vec = T1%Particle(n)%get_pos()
   end function get_pos_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_pos_all(T1, pos)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: pos(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         call T1%Particle(n)%set_pos(pos(ndim * (n - 1) + 1:ndim * n))
      end do

      call T1%geom_changed()

   end subroutine set_pos_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_pos_all(T1) result(pos)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: pos(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         pos(ndim * (n - 1) + 1:ndim * n) = T1%Particle(n)%get_pos()
      end do
   end function get_pos_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_mass_comp(T1, n) result(mass)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: n
      real(kind=DefReal) :: mass

      mass = T1%Particle(n)%mass

   end function get_mass_comp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_mass_vec(T1) result(mass)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: mass(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         mass(ndim * (n - 1) + 1:n * ndim) = T1%Particle(n)%mass
      end do

   end function get_mass_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_width_comp(T1, n) result(width)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: n
      real(kind=DefReal) :: width

      width = T1%Particle(n)%width
   end function get_width_comp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_width_vec(T1) result(width)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: width(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         width(ndim * (n - 1) + 1:n * ndim) = T1%Particle(n)%width
      end do

   end function get_width_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_mom_vec(T1, i) result(mom)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt) :: i
      real(kind=DefReal) :: mom(T1%Particle(i)%NumDimensions)

      mom = T1%Particle(i)%get_mom()

   end function get_mom_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_mom_all(T1) result(mom)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: mom(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         mom(ndim * (n - 1) + 1:ndim * n) = T1%Particle(n)%get_mom()
      end do

   end function get_mom_all

   function get_mom2_vec(T1, i) result(mom2)
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt) :: i
      real(kind=DefReal) :: mom2(T1%Particle(i)%NumDimensions)

      mom2 = T1%Particle(i)%get_mom2()
   end function get_mom2_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_mom2_all(T1) result(mom2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: mom2(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         mom2(ndim * (n - 1) + 1:ndim * n) = T1%Particle(n)%get_mom2()
      end do
   end function get_mom2_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_mom_vec(T1, n, mom)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      integer(kind=DefInt), intent(in) :: n
      real(kind=DefReal), intent(in) :: mom(T1%Particle(1)%NumDimensions)

      call T1%Particle(n)%set_mom(mom)
   end subroutine set_mom_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_mom_all(T1, mom)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      real(kind=DefReal), intent(in) :: mom(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         call T1%Particle(n)%set_mom(mom(ndim * (n - 1) + 1:ndim * n))
      end do

   end subroutine set_mom_all

   subroutine set_mom2_vec(T1, n, mom)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      integer(kind=DefInt), intent(in) :: n
      real(kind=DefReal), intent(in) :: mom(T1%Particle(1)%NumDimensions)

      call T1%Particle(n)%set_mom2(mom)
   end subroutine set_mom2_vec

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_mom2_all(T1, mom2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      real(kind=DefReal), intent(in) :: mom2(T1%NumDimensions)

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%Particle(1)%NumDimensions

      do n = 1, npart
         call T1%Particle(n)%set_mom2(mom2(ndim * (n - 1) + 1:ndim * n))
      end do

   end subroutine set_mom2_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_vel_all(T1) result(vel)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: vel(T1%NumDimensions)

      vel = T1%get_mom() / T1%get_mass()
   end function get_vel_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! - - - - - - - -
   function get_vel2_all(T1) result(vel)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: vel(T1%NumDimensions)

      vel = T1%get_mom2() / T1%get_mass()
   end function get_vel2_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_vel_all(T1, vel)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      real(kind=DefReal), intent(in) :: vel(T1%NumDimensions)

      call T1%set_mom(T1%get_mass() * vel)
   end subroutine set_vel_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_vel2_all(T1, vel)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory) :: T1
      real(kind=DefReal), intent(in) :: vel(T1%NumDimensions)

      call T1%set_mom2(T1%get_mass() * vel)
   end subroutine set_vel2_all

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Distance_Trajectory(T1, i, j) result(R)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j
      real(kind=DefReal) :: R

      integer(kind=DefInt) :: npart

      npart = T1%NumParticles
      if (i < 1 .or. npart < i) then
         call FMS_DieError("Distance: index i out of range")
      end if
      if (j < 1 .or. npart < j) then
         call FMS_DieError("Distance: index j out of range")
      end if
      if (i == j) then
         call FMS_DieError("Distance: i and j the same ")
      end if

      R = FMS_Distance(T1%Particle(i), T1%Particle(j))

   end function FMS_Distance_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_D_Distance_Trajectory(T1, i, j) result(dRdX)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j
      real(kind=DefReal) :: dRdX(T1%NumDimensions)

      real(kind=DefReal) :: R, Rij(3)

      R = FMS_Distance(T1, i, j)
      Rij = T1%get_pos(i) - T1%get_pos(j)

      dRdX = 0.d0
      dRdX(3 * i - 2:3 * i) = +Rij / R
      dRdX(3 * j - 2:3 * j) = -Rij / R

   end function FMS_D_Distance_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Angle_Trajectory(T1, i, j, k) result(theta)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j, k
      real(kind=DefReal) :: theta

      integer(kind=DefInt) :: npart

      npart = T1%NumParticles
      if (i < 1 .or. npart < i) call FMS_DieError("FMS_Angle_Trajectory : index i out of range")
      if (j < 1 .or. npart < j) call FMS_DieError("FMS_Angle_Trajectory : index j out of range")
      if (k < 1 .or. npart < k) call FMS_DieError("FMS_Angle_Trajectory : index k out of range")

      theta = FMS_Angle(T1%Particle(i), T1%Particle(j), T1%Particle(k))

   end function FMS_Angle_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_D_Angle_Trajectory(T1, i, j, k) result(dthetadX)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j, k
      real(kind=DefReal) :: dthetadX(T1%NumDimensions)

      real(kind=DefReal) :: theta, & ! angle between atoms
                            Xij(3), Nij, & ! vector & norm from j to i
                            Xkj(3), Nkj ! vector & norm from j to k

      theta = FMS_Angle_Trajectory(T1, i, j, k)

      Xij = T1%get_pos(i) - T1%get_pos(j)

      Xkj = T1%get_pos(k) - T1%get_pos(j)

      Nij = sqrt(dot_product(Xij, Xij))
      Nkj = sqrt(dot_product(Xkj, Xkj))

      dthetadX = 0.d0

      dthetadX(3 * i - 2:3 * i) = -(Nij**2 * Xkj - dot_product(Xij, Xkj) * Xij) / (sin(theta) * Nij**3 * Nkj)

      dthetadX(3 * k - 2:3 * k) = -(Nkj**2 * Xij - dot_product(Xij, Xkj) * Xkj) / (sin(theta) * Nij * Nkj**3)

      dthetadX(3 * j - 2:3 * j) = -dthetadX(3 * i - 2:3 * i) - dthetadX(3 * k - 2:3 * k)

   end function FMS_D_Angle_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Dihedral_Trajectory(T1, i, j, k, l) result(phi)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j, k, l
      real(kind=DefReal) :: phi
      integer(kind=DefInt) :: npart

      npart = T1%NumParticles
      if (i < 1 .or. npart < i) then
         call FMS_DieError("Dihedral: index i out of range")
      end if
      if (j < 1 .or. npart < j) then
         call FMS_DieError("Dihedral: index j out of range")
      end if
      if (k < 1 .or. npart < k) then
         call FMS_DieError("Dihedral: index k out of range")
      end if
      if (l < 1 .or. npart < l) then
         call FMS_DieError("Dihedral: index l out of range")
      end if

      phi = FMS_Dihedral(T1%Particle(i), T1%Particle(j), T1%Particle(k), T1%Particle(l))

!     Incorrect definition of dihedral, consistent with
!     The analytic derivative below
!     Rji=FMS_GetPosition(T1,j)-FMS_GetPosition(T1,i)
!     Rkj=FMS_GetPosition(T1,k)-FMS_GetPosition(T1,j)
!     Rlk=FMS_GetPosition(T1,l)-FMS_GetPosition(T1,k)
!     v1 = unit_vector(cross(Rji, Rkj))
!     v2 = unit_vector(cross(Rlk,-Rkj))

!     cosd=dot_product(v1,v2)
!     v3 = cross(v1,v2)
!     sind=dot_product(v3,Rkj)/vector_norm(Rkj)
!     phi=atan2(sind,cosd)

   end function FMS_Dihedral_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_D_Dihedral_Trajectory(T1, i, j, k, l) result(dphidX)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! This code only works for an incorrect definition of the dihederal that is
      ! commented out in the above function
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), intent(in) :: i, j, k, l
      real(kind=DefReal) :: dphidX(T1%NumDimensions)

      real(kind=DefReal) :: Rji(3), Rkj(3), Rlk(3), Rik(3), Rjl(3), & ! displacement vectors
                            Dji, Dkj, Dlk, Dik, Djl, & ! norms
                            djikj, dlkkj ! dot products

      Rji = T1%get_pos(j) - T1%get_pos(i)
      Rkj = T1%get_pos(k) - T1%get_pos(j)
      Rlk = T1%get_pos(l) - T1%get_pos(k)

      Rik = cross(Rji, Rkj)
      Rjl = cross(Rlk, -Rkj)

      Dji = vector_norm(Rji)
      Dkj = vector_norm(Rkj)
      Dlk = vector_norm(Rlk)
      Dik = vector_norm(Rik)
      Djl = vector_norm(Rjl)

      djikj = dot_product(Rji, Rkj)
      dlkkj = dot_product(Rlk, Rkj)

      dphidX = 0.d0

!     Working with alternate def
      dphidX(3 * i - 2:3 * i) = -Dkj / Dik**2 * Rik
      dphidX(3 * l - 2:3 * l) = +Dkj / Djl**2 * Rjl

      dphidX(3 * j - 2:3 * j) = (-djikj / Dkj**2 - 1.d0) * dphidX(3 * i - 2:3 * i) + dlkkj / Dkj**2 * dphidX(3 * l - 2:3 * l)
      dphidX(3 * k - 2:3 * k) = (-dlkkj / Dkj**2 - 1.d0) * dphidX(3 * l - 2:3 * l) + djikj / Dkj**2 * dphidX(3 * i - 2:3 * i)

   end function FMS_D_Dihedral_Trajectory

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function center_of_mass(T1, i1, i2, mass_out) result(com)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), optional, intent(in) :: i1, i2 ! first and last atoms
      real(kind=DefReal), optional, intent(out) :: mass_out ! optional total mass passed out
      real(kind=DefReal) :: com(3)

      integer(kind=DefInt) :: n, n1, n2
      real(kind=DefReal) :: mass

      ! default range
      n1 = 1
      n2 = T1%NumParticles

      ! override default
      if (present(i1)) n1 = i1
      if (present(i2)) n2 = i2

      mass = 0.d0
      com = 0.d0

      do n = n1, n2
         mass = mass + T1%Particle(n)%mass
         com = com + T1%Particle(n)%mass * T1%Particle(n)%get_pos()
      end do

      com = com / mass

      if (present(mass_out)) mass_out = mass

   end function center_of_mass

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_pop(T1) result(pop)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: pop

      pop = T1%pop
   end function get_pop

!>
!!    Rescales phases into the interval \f$ [ 0, 2\pi ] \f$
!!    @ingroup propagation
!<
   subroutine rescale_phases(T1)
      class(T_Trajectory) :: T1

      T1%Phase = mod(T1%Phase, D2 * Pi)
      if (T1%Phase < 0) T1%Phase = T1%Phase + D2 * Pi
   end subroutine rescale_phases

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   logical function is_dead(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      class(T_Trajectory), intent(in) :: T

      real(DefReal) :: CurrentTime, DeadTime

      CurrentTime = T%get_time()
      DeadTime = T%DeadTime

      is_dead = .not. T%zCent .and. ((CurrentTime - DeadTime) > (gldDecoherenceTime + FPZero))
   end function is_dead

!     Private routine to flip flags
!>
!! Flip all flags after a position has changed.
!!
!! Scope: TrajectoryModule
!<
   subroutine geom_changed(T1)
      class(T_Trajectory) :: T1

!     Flags for electronic structure calculations
      T1%ESFlags%ZPotEnCurrent = .false.
      T1%ESFlags%ZDerivCurrent = .false.
      T1%ESFlags%ZSOMCurrent = .false.

      T1%ESFlags%ZTransDipsCurrent = .false.
      T1%ESFlags%ZTransDipsCurrentxf = .false.
      T1%ESFlags%ZDipolesCurrent = .false.
      T1%ESFlags%zMMPotCurrent = .false.
      T1%ESFlags%zMMForceCurrent = .false.
      T1%ESFlags%zModPotCurrent = .false.

!     Flags for bundle quantity calculations
      T1%BFlags%ZBundleCurrent = .false.
      T1%BFlags%ZAmpDotCurrent = .false.

   end subroutine geom_changed

!>
!!    Calculates the vector describing the pyramidalization plane.
!!
!!    This plane is defined by the three atoms involved in pyramidalization process, e.g.\ H-C-H in ethylene.
!!
!!    <pre>
!!              3rd
!!             /
!!      4th=1st
!!             |
!!              2nd
!!    </pre>
!!    \author Ben G. Levine
!!    @param nAtoms Number of atoms in system
!!    @param ccold Matrix of cartesian coordinates
!!    @param dAngle (Output) Dihedral angle in radians
!!    @param iAtom1 ID of first atom
!!    @param iAtom2 ID of second atom
!!    @param iAtom3 ID of third atom
!!    @param iAtom4 ID of fourth atom
!!    @ingroup analysis
!<
!!    TODO(DH): Turn this functions into FMS_PyrAngle_Trajectory with the
!!    same interface as FMS_Dihedral_Trajectory et al
   subroutine FMS_CalcPyrAngle(natoms, ccold, dangle, iatom1, iatom2, iatom3, iatom4)
      integer(kind=DefInt) :: iAtom1, iAtom2, iAtom3, iAtom4, i, j, natoms
      real(kind=DefReal) :: ccold(natoms, 3), vec(3, natoms, natoms), dangle

      do i = 1, nAtoms
         do j = 1, nAtoms
            vec(1, i, j) = ccold(j, 1) - ccold(i, 1)
            vec(2, i, j) = ccold(j, 2) - ccold(i, 2)
            vec(3, i, j) = ccold(j, 3) - ccold(i, 3)
         end do
      end do

      dAngle = dihed(vec(1:3, iAtom4, iAtom1), vec(1:3, iatom2, iAtom3), vec(1:3, iAtom1, iAtom2), vec(1:3, iAtom1, iAtom3))

   contains
!>
!! Calculates the dihedral angle \f$ \theta \f$ between the planes defined by vectors a and b, and c and d respectively.
!! \f[
!! \theta = \cos^{-1} \frac{(\vec{a} \times \vec{b}) \cdot (\vec{c} \times \vec{d})}{\Vert \vec{e} \Vert \Vert \vec{f} \Vert}
!! \f]
!! @param a first vector of first plane
!! @param b second vector of second plane
!! @param c first vector of first plane
!! @param d second vector of second plane
!<
      function dihed(a, b, c, d)
         real(kind=DefReal), intent(in) :: a(3), b(3), c(3), d(3)
         real(kind=DefReal) :: dihed
         real(kind=DefReal) :: e(3), f(3), g(3), h(3)
         real(kind=DefReal) :: tmp

         call crossprod(a, b, e)
         call crossprod(c, d, f)
         g(1:3) = e(1:3) / sqrt(dot_product(e, e))
         h(1:3) = f(1:3) / sqrt(dot_product(f, f))

         tmp = dot_product(g, h)
         dihed = acos(tmp)
      end function dihed
!>
!! Calculates the cross product between two three-dimensional vectors \f$ \vec{a} \times \vec{b} \f$
!! @param a left argument of cross product
!! @param b right argument of cross product
!! @param cross result of cross product
!<
      subroutine crossprod(a, b, cross)
         real(kind=DefReal), intent(in) :: a(3), b(3)
         real(kind=DefReal), intent(out) :: cross(3)

         cross(1) = a(2) * b(3) - a(3) * b(2)
         cross(2) = a(3) * b(1) - a(1) * b(3)
         cross(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine crossprod

   end subroutine FMS_CalcPyrAngle

end module TrajectoryModule
