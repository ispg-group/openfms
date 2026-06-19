!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Specification of the "Bundle" datatype and methods.
!!
!!    The T_TrajectoryBundle datatype contains a specification of the
!!    full nuclear/electronic wavefunction at a given instant in
!!    time. The data-type itself is composed of several trajectory basis
!!    functions, and also holds representations of quantum mechanical
!!    operators (i.e. the Hamiltonian and overlap matrix) in this
!!    trajectory basis.
!!
!!    These operators are held in the protected BundleMatrices datatype;
!!    analogously to the T_Trajectory type, they must be retreived using
!!    publicly scoped accessor routines.
!!
!!    Methods in this class are centered around building, storing, and
!!    retrieving these operators.
!<
module BundleModule
   use GlobalModule
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_isBundleCurrent, &
                                    FMS_isAmpDotCurrent, FMS_BundleUpdated, FMS_AmpDotUpdated

   implicit none

   private
   public :: t_TrajectoryBundle
   public :: assignment(=)

!>
!!    Physical operators over the full wavefunction, in the trajectory basis
!<
   type T_BundleMatrices

      complex(kind=DefComp), allocatable, dimension(:, :) :: H, & !< Hamiltonian Matrix
                                                             S, & !< Overlap Matrix
                                                             SDot, & !< Time Derivative Matrix
                                                             SInv, & !< Overlap Inverse
                                                             Sp5i, & !< \f$S^{-1/2}\f$
                                                             Sp5, & !< \f$S^{ 1/2}\f$
                                                             HEff, & !< Effective Hamiltonian
                                                             HEff1 !< Effective Hamiltonian before multiplying by \f$S^{-1}\f$

      complex(kind=DefComp), allocatable :: AmpDot(:) !< Time derivative of amplitude for each Trajectory

   end type T_BundleMatrices

!>
!!    Datatype holding a collection of trajectories.
!!
!!    A bundle represents the full
!!    wavefunction, represented as a linear combination of the trajectory
!!    basis functions. It also stores various observables in this basis.
!<
   type T_TrajectoryBundle

      integer(kind=DefInt) :: NumStates, & !< Number of Electronic States
                              NumSing, & !< Number of Electronic Singlet States
                              NumTrip, & !< Number of Electronic Triplet States
                              NumParticles, & !< Number of atoms in this system
                              NumTraj, & !< Number of LIVE trajectories in this bundle
                              NumDeadTraj !< Number of non-propagating trajectories

      real(kind=DefReal) :: CurrentTime !< Time for bundle

      ! TODO: This type actually private
      type(T_BundleMatrices) :: BMatrices !< Private nuclear basis info

      integer(kind=DefInt) :: NCBFs !< Number of contracted basis functions

      type(T_Trajectory), allocatable :: Trajectory(:), & !< Storage of all live trajectories
                                         Centroids(:), & !< \f$N_{traj}(N_{traj}-1)/2 \f$ centroids for live trajectories
                                         DeadTraj(:) !< Storage for dead (non-propagating) trajectories
      complex(kind=DefComp), allocatable :: DeadH(:, :) !< Hamiltonian for non-propagating trajectories

   contains
      ! Type-bound subroutines, call as `call Bundle%subroutine()`
      procedure :: create => create_bundle !< Allocates internal bundle structures
      procedure :: destroy => destroy_bundle !< Deallocates internal bundle structures
      procedure :: add_traj => add_traj_to_bundle !< Adds trajectory to bundle
      procedure :: add_traj_triplet => add_traj_triplet_to_bundle !< Adds triplet trajectory to bundle
      procedure, public :: copy_from

   end type T_TrajectoryBundle

!     Overloading functions:
   interface assignment(=)
      module procedure assign_bundle
   end interface

contains
!--------------- Module Procedures ---------------!

! Constructors/Destructors/Assignment
!>
!>    Memory allocation for creating a trajectory bundle
!>
   subroutine create_bundle(B1, NumTraj, NumDeadTraj, NumStates, NumParticles, NCBFs)
      class(T_TrajectoryBundle), intent(inout) :: B1
      integer(kind=DefInt), intent(in) :: NumTraj, NumDeadTraj, NumStates, NumParticles, NCBFs

      integer(kind=DefInt) :: ITraj, iCent, ICBF, JCBF

      if (allocated(B1%Trajectory)) then
         if (B1%NumTraj == NumTraj .and. B1%NumDeadTraj == NumDeadTraj .and. &
             B1%NumStates == NumStates .and. B1%NCBFs == NCBFs .and. &
             B1%NumParticles == NumParticles) return
      end if

      call B1%destroy()

      if (NumTraj < 0) call FMS_DieError('create_bundle: NumTraj must be >= 0')
      if (NumStates < 0) call FMS_DieError('create_bundle: NumStates must be >= 0')
      if (NumParticles < 0) call FMS_DieError('create_bundle: NumParticles must be >= 0')

      B1%NumTraj = NumTraj
      B1%NumDeadTraj = NumDeadTraj
      B1%NumStates = NumStates
      B1%NumParticles = NumParticles
      B1%NCBFs = NCBFs
      B1%CurrentTime = 0.0d0

      allocate (B1%Trajectory(NumTraj))
      do ITraj = 1, NumTraj
         call B1%Trajectory(ITraj)%create(NumParticles, NumStates)
         B1%Trajectory(ITraj)%TrajID = -1
      end do

!bfec
      if (glzCentroids) then
         if (NCBFs > 1) then
            allocate (B1%Centroids(((NCBFs - 1) * NCBFs) / 2))
            do ICBF = 2, B1%NCBFs
               do JCBF = 1, ICBF - 1
                  iCent = ((iCBF - 1) * (iCBF - 2)) / 2 + jCBF
                  call B1%Centroids(iCent)%create(NumParticles, NumStates)
                  B1%Centroids(ICent)%zCent = .true.
                  B1%Centroids(iCent)%TrajID = iCent
                  B1%Centroids(iCent)%CentID = [-1, -1]
               end do
            end do
         end if
      end if

      allocate (B1%BMatrices%H(NumTraj, NumTraj), &
                B1%BMatrices%S(NumTraj, NumTraj), &
                B1%BMatrices%SDot(NumTraj, NumTraj), &
                B1%BMatrices%SInv(NumTraj, NumTraj), &
                B1%BMatrices%Sp5i(NumTraj, NumTraj), &
                B1%BMatrices%Sp5(NumTraj, NumTraj), &
                B1%BMatrices%HEff(NumTraj, NumTraj), &
                B1%BMatrices%HEff1(NumTraj, NumTraj), &
                B1%BMatrices%AmpDot(NumTraj))

      if (NumDeadTraj > 0) then
         allocate (B1%DeadTraj(NumDeadTraj))
         do ITraj = 1, NumDeadTraj
            call B1%DeadTraj(ITraj)%create(NumParticles, NumStates)
            B1%DeadTraj(ITraj)%TrajID = -1
         end do

         allocate (B1%DeadH(NuMDeadTraj, NumDeadTraj))
      end if

   end subroutine create_bundle

!>
!!    Memory deallocation to destroy a trajectory bundle structure
!<
   subroutine destroy_bundle(B)
      class(T_TrajectoryBundle), intent(inout) :: B

      integer(kind=DefInt) :: n

      B%NumTraj = 0
      B%NumStates = 0
      B%NumParticles = 0

      ! TODO: Do we need to zero-out the Current Time?
      ! What about NCBFS? NumDeadTraj?
      ! B%CurrentTime = 0.0D0
      ! B%NCBFS = 0
      ! B%NumDeadTraj = 0

      if (allocated(B%Trajectory)) then
         do n = 1, size(B%Trajectory)
            call B%Trajectory(n)%destroy()
         end do
         deallocate (B%Trajectory)
      end if

!bfec
      if (glzCentroids) then
         if (allocated(B%Centroids)) then
            do n = 1, size(B%Centroids)
               call B%Centroids(n)%destroy()
            end do
            deallocate (B%Centroids)
         end if
      end if

      if (allocated(B%DeadTraj)) then
         do n = 1, size(B%DeadTraj)
            call B%DeadTraj(n)%destroy()
         end do
         deallocate (B%DeadTraj)
      end if

      if (allocated(B%DeadH)) deallocate (B%DeadH)

      if (allocated(B%BMatrices%H)) deallocate (B%BMatrices%H)
      if (allocated(B%BMatrices%S)) deallocate (B%BMatrices%S)
      if (allocated(B%BMatrices%SDot)) deallocate (B%BMatrices%SDot)
      if (allocated(B%BMatrices%SInv)) deallocate (B%BMatrices%SInv)
      if (allocated(B%BMatrices%Sp5i)) deallocate (B%BMatrices%Sp5i)
      if (allocated(B%BMatrices%Sp5)) deallocate (B%BMatrices%Sp5)
      if (allocated(B%BMatrices%HEff)) deallocate (B%BMatrices%HEff)
      if (allocated(B%BMatrices%HEff1)) deallocate (B%BMatrices%HEff1)
      if (allocated(B%BMatrices%AmpDot)) deallocate (B%BMatrices%AmpDot)

   end subroutine destroy_bundle

!>
!!    Memory allocation and book-keeping for B1 = B2 assignment,
!!    where B1 and B2 are trajectory bundles.
!<
   subroutine assign_bundle(B1, B2)
      !> Copy to B1
      type(T_TrajectoryBundle), intent(inout) :: B1
      !> Copy from B2
      type(T_TrajectoryBundle), intent(in) :: B2
      integer(kind=DefInt) :: ITraj

      if ((B1%NumStates /= B2%NumStates) .or. (B1%NumTraj /= B2%NumTraj) .or. B1%NumDeadTraj /= B2%NumDeadTraj) then
         call B1%destroy()
         call B1%create(numtraj=B2%NumTraj, &
                        numdeadtraj=B2%NumDeadTraj, &
                        numstates=B2%NumStates, &
                        numparticles=B2%NumParticles, &
                        ncbfs=B2%NCBFs)
      end if

      B1%CurrentTime = B2%CurrentTime

      !GAIMS added
      B1%NCBFs = B2%NCBFs

      if (B2%NumTraj > 0) then
         B1%BMatrices%H = B2%BMatrices%H
         B1%BMatrices%S = B2%BMatrices%S
         B1%BMatrices%SDot = B2%BMatrices%SDot
         B1%BMatrices%SInv = B2%BMatrices%SInv
         B1%BMatrices%Sp5i = B2%BMatrices%Sp5i
         B1%BMatrices%Sp5 = B2%BMatrices%Sp5
         B1%BMatrices%HEff = B2%BMatrices%HEff
         B1%BMatrices%HEff1 = B2%BMatrices%HEff1
         B1%BMatrices%AmpDot = B2%BMatrices%AmpDot

         do ITraj = 1, B2%NumTraj
            B1%Trajectory(ITraj) = B2%Trajectory(ITraj)
            ! TODO(DH) Is this correct?? Or should it be:
            ! if (FMS_IsBundleCurrent(B2%Trajectory(ITraj))) call FMS_BundleUpdated(B1%Trajectory(ITraj))
            ! if (FMS_IsAmpDotCurrent(B2%Trajectory(ITraj))) call FMS_AmpDotUpdated(B1%Trajectory(ITraj))
            ! See also comment at the end of FMS_AssignTrajectory
            if (FMS_IsBundleCurrent(B1%Trajectory(ITraj))) call FMS_BundleUpdated(B2%Trajectory(ITraj))
            if (FMS_IsAmpDotCurrent(B1%Trajectory(ITraj))) call FMS_AmpDotUpdated(B2%Trajectory(ITraj))
         end do

         if (glzCentroids) then
            do ITraj = 1, (B2%NCBFs * (B2%NCBFs - 1)) / 2
               B1%Centroids(ITraj) = B2%Centroids(ITraj)
            end do
         end if
      end if

      if (B2%NumDeadTraj > 0) then
         B1%NumDeadTraj = B2%NumDeadTraj
         do iTraj = 1, B2%NumDeadTraj
            B1%DeadTraj(iTraj) = B2%DeadTraj(iTraj)
         end do
         B1%DeadH = B2%DeadH
      end if
   end subroutine assign_bundle

   subroutine copy_from(B1, B2)
      !> Copy to B1
      class(T_TrajectoryBundle), intent(inout) :: B1
      !> Copy from B2
      class(T_TrajectoryBundle), intent(in) :: B2

      call assign_bundle(B1, B2)
   end subroutine copy_from

!>
!!    Add a new trajectory to a bundle
!!
!!    Add a new trajectory to the bundle; handle memory allocation
!!    and expansion of trajectory-basis operators.
!!    This is a type-bound procedure that should be called as:
!!    Bundle%add_traj(Trajectory)
!<
   subroutine add_traj_to_bundle(B1, T1)
      class(T_TrajectoryBundle), intent(inout) :: B1
      type(T_Trajectory), intent(inout) :: T1

      integer(kind=DefInt) :: ITraj, NumTrajP1, iCent
      integer(kind=DefInt) :: ICBF, NumCBFP1
      type(T_TrajectoryBundle) :: BTemp

      NumTrajP1 = B1%NumTraj + 1
      NumCBFP1 = B1%NCBFs + 1
      call BTemp%create(numtraj=NumTrajP1, &
                        numdeadtraj=B1%NumDeadTraj, &
                        numstates=B1%NumStates, &
                        numparticles=T1%NumParticles, &
                        ncbfs=NumCBFP1)
      T1%TrajID = BTemp%NumTraj + BTemp%NumDeadTraj

!     Copy trajectories
      do ITraj = 1, B1%NumTraj
         BTemp%Trajectory(ITraj) = B1%Trajectory(ITraj)
      end do

!     Copy dead trajectories
      do iTraj = 1, B1%NumDeadTraj
         BTemp%DeadTraj(iTraj) = B1%DeadTraj(iTraj)
      end do

!     Copy centroids
!bfec
      if (glzCentroids) then
!       do ITraj=1,((B1%NumTraj*(B1%NumTraj-1))/2)
         do ITraj = 1, ((B1%NCBFs * (B1%NCBFs - 1)) / 2)
            BTemp%Centroids(ITraj) = B1%Centroids(ITraj)
         end do

         BTemp%Trajectory(NumTrajP1) = T1
!     GAIMS added
         BTemp%NCBFs = B1%NCBFs + 1
         BTemp%Trajectory(NumTrajP1)%CBF = BTemp%NCBFs
!     GAIMS end added

!     Create centroids for new trajectory
!       do iTraj=1,B1%NumTraj
!          iCent=(B1%NumTraj*(B1%NumTraj-1))/2+iTraj
         do iCBF = 1, B1%NCBFs
            iCent = (B1%NCBFs * (B1%NCBFs - 1)) / 2 + iCBF
            BTemp%Centroids(iCent)%zCent = .true.
            BTemp%Centroids(iCent)%TrajID = -iCent
            BTemp%Centroids(iCent)%CentID = [BTemp%Trajectory(NumTrajP1)%CBF, iCBF]
!     $        (/ T1%CBF, B1%Trajectory(iTraj)%CBF /)
!          BTemp%Centroids(iCBF)%NumStates=
!     $    BTemp%Trajectory(1)%NumStates
         end do
      end if
      BTemp%CurrentTime = B1%CurrentTime

      BTemp%Trajectory(NumTrajP1)%triplet = .false.
      BTemp%Trajectory(NumTrajP1)%Ms = 2

!     GAIMS end added

      if (allocated(B1%DeadH)) BTemp%DeadH = B1%DeadH

      call B1%copy_from(BTemp)

      call BTemp%destroy()

   end subroutine add_traj_to_bundle

!>
!!    Add a new triplet trajectory to a bundle.
!!
!!    Unlike in add_traj_to_bundle, the triplet trajectory
!!    is copied into three sublevels.
!<
   subroutine add_traj_triplet_to_bundle(B1, T1)
      class(T_TrajectoryBundle), intent(inout) :: B1
      type(T_Trajectory), intent(inout) :: T1
      type(T_Trajectory) :: T2, T3

      integer(kind=DefInt) :: ITraj, NumTrajP1, iCent
      integer(kind=DefInt) :: i, j, NumCBFP1
      type(T_TrajectoryBundle) :: BTemp

      T2 = T1
      T3 = T1

      NumTrajP1 = B1%NumTraj + 3 ! Add the three sublevels of Ms=-1,0,1
      NumCBFP1 = B1%NCBFs + 1

      call BTemp%create(numtraj=NumTrajP1, &
                        numdeadtraj=B1%NumDeadTraj, &
                        numstates=B1%NumStates, &
                        numparticles=T1%NumParticles, &
                        ncbfs=NumCBFP1)
      T1%TrajID = BTemp%NumTraj + BTemp%NumDeadTraj - 2
      T2%TrajID = BTemp%NumTraj + BTemp%NumDeadTraj - 1
      T3%TrajID = BTemp%NumTraj + BTemp%NumDeadTraj

!     Copy trajectories
      do ITraj = 1, B1%NumTraj
         BTemp%Trajectory(ITraj) = B1%Trajectory(ITraj)
      end do

!     Copy dead trajectories
      do iTraj = 1, B1%NumDeadTraj
         BTemp%DeadTraj(iTraj) = B1%DeadTraj(iTraj)
      end do

!     Copy centroids
      if (glzCentroids) then

         do ITraj = 1, (B1%NCBFs * (B1%NCBFs - 1)) / 2
            BTemp%Centroids(ITraj) = B1%Centroids(ITraj)
         end do

!        Add the new trajectories to the bundle
         BTemp%Trajectory(NumTrajP1 - 2) = T1
         BTemp%Trajectory(NumTrajP1 - 1) = T2
         BTemp%Trajectory(NumTrajP1) = T3
         BTemp%NCBFs = B1%NCBFs + 1
         BTemp%Trajectory(NumTrajP1)%CBF = BTemp%NCBFs
         BTemp%Trajectory(NumTrajP1 - 1)%CBF = BTemp%NCBFs
         BTemp%Trajectory(NumTrajP1 - 2)%CBF = BTemp%NCBFs

         do i = 2, BTemp%NCBFs
            do j = 1, i - 1
               iCent = ((i - 2) * (i - 1)) / 2 + j
               BTemp%Centroids(iCent)%zCent = .true.
               BTemp%Centroids(iCent)%TrajID = -iCent
               BTemp%Centroids(iCent)%CentID = [i, j]
            end do
         end do

      end if

      BTemp%CurrentTime = B1%CurrentTime

      if (allocated(B1%DeadH)) BTemp%DeadH = B1%DeadH

      call B1%copy_from(BTemp)

      B1%Trajectory(NumTrajP1)%Ms = 3
      B1%Trajectory(NumTrajP1 - 1)%Ms = 1
      B1%Trajectory(NumTrajP1 - 2)%Ms = 2
      B1%Trajectory(NumTrajP1)%triplet = .true.
      B1%Trajectory(NumTrajP1 - 1)%triplet = .true.
      B1%Trajectory(NumTrajP1 - 2)%triplet = .true.

      call BTemp%destroy()

   end subroutine add_traj_triplet_to_bundle

end module BundleModule
