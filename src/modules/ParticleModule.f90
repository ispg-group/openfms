!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!     No portion may be distributed, modified, or used without express
!     written permission from one of the copyright holders.
!>
!!    Derived types for particles
!<
module ParticleModule
   use GlobalModule, only: DefReal, DefInt, fmiOut, FMS_DieError

   implicit none
   private
!>
!!    Public data for particle structure
!<
   public :: T_Particle, assignment(=)
   public :: FMS_Distance, FMS_Angle, FMS_Dihedral
   public :: cross, vector_norm

   type T_Particle
      integer(kind=DefInt) :: ParticleID, & !< ID for potentials, etc.
                              NumDimensions !< Number of spatial dimensions known to this particle

      real(kind=DefReal), allocatable, private :: &
         Position(:), & !< Position,  an array of size NumDimensions
         Momentum(:), & !< Momentum, an array of size NumDimensions
         Momentum2(:) !< Second momentum, for centroids to store complex component of off-diagonal momentum
      real(kind=DefReal) :: Mass !< Particle mass
      real(kind=DefReal) :: Width, & !< Position uncertainty
                            Charge, & !< Atomic partial charge, e.g. Mulliken charge read from Molpro
                            AtomicNum !< Nuclear charge
      character(2) :: Elmnt !< Particle chemical identity, for MOPAC
   contains
      private
      generic, public :: get_pos => get_position_vec
      generic, public :: get_pos => get_position_component
      generic, public :: get_mom => get_momentum_vec
      generic, public :: get_mom => get_momentum_component
      generic, public :: get_mom2 => get_momentum2_vec
      generic, public :: set_pos => set_position_vec
      generic, public :: set_pos => set_position_component
      generic, public :: set_mom => set_momentum_vec
      generic, public :: set_mom => set_momentum_component
      generic, public :: set_mom2 => set_momentum2_vec
      procedure, public :: create => create_particle
      procedure, public :: destroy => destroy_particle
      procedure, public :: print_info => print_particle_info
      procedure :: get_position_component
      procedure :: get_position_vec
      procedure :: get_momentum_component
      procedure :: get_momentum_vec
      procedure :: get_momentum2_vec
      procedure :: set_position_component
      procedure :: set_position_vec
      procedure :: set_momentum_component
      procedure :: set_momentum_vec
      procedure :: set_momentum2_vec
   end type T_Particle

!>
!!     Overloaded definitions
!<
   interface FMS_Distance
      module procedure FMS_Distance_Particle
   end interface

   interface FMS_Angle
      module procedure FMS_Angle_Particle
   end interface

   interface FMS_Dihedral
      module procedure FMS_Dihedral_Particle
   end interface

   interface assignment(=)
      module procedure assign_particle
   end interface

contains

!     Constructor/destructor/assignment
!>
!!    Memory allocation for creating a particle structure
!!    \param P1 particle to create
!!    \param ID particle ID
!!    \param NumDim dimensionality of particle (must be >= 0)
!<
   subroutine create_particle(P1, ID, NumDim)
      class(T_Particle), intent(inout) :: P1
      integer(kind=DefInt), intent(in) :: ID, NumDim

      call P1%destroy()

      if (NumDim < 0) then
         call FMS_DieError('create_particle: NumDim must be >= 0')
         return
      end if

      P1%ParticleID = ID
      P1%NumDimensions = NumDim
      P1%Width = 0.0d0
      P1%Mass = 0.0d0
      P1%Charge = 0.0d0
      P1%AtomicNum = 0.0d0
      P1%Elmnt = 'XX'

      allocate (P1%Position(NumDim))
      P1%Position = 0.d0

      allocate (P1%Momentum(NumDim))
      P1%Momentum = 0.d0

      allocate (P1%Momentum2(NumDim))
      P1%Momentum2 = 0.d0

   end subroutine create_particle

!>
!!    Memory allocation and book-keeping for carrying out assignment P1=P2
!!    \param P1 Particle to be assigned to
!!    \param P2 Particle being assigned from
!<
   subroutine assign_particle(P1, P2)
      type(T_Particle), intent(inout) :: P1
      type(T_Particle), intent(in) :: P2

      integer(kind=DefInt) :: ndim1, ndim2

      ndim1 = P1%NumDimensions
      ndim2 = P2%NumDimensions

      if (ndim1 /= ndim2) then
         call P1%destroy()
         call P1%create(id=P2%ParticleID, numdim=ndim2)
      end if

      P1%Width = P2%Width
      P1%Mass = P2%Mass
      P1%NumDimensions = P2%NumDimensions
      P1%Elmnt = P2%Elmnt
      P1%Charge = P2%Charge
      P1%AtomicNum = P2%AtomicNum
      P1%Position = P2%Position
      P1%Momentum = P2%Momentum
      P1%Momentum2 = P2%Momentum2
   end subroutine assign_particle

!>
!!    Memory deallocation to destroy a particle structure
!!    \param Particle Particle to annihilate
!<
   subroutine destroy_particle(P1)
      class(T_Particle), intent(inout) :: P1

      P1%ParticleID = 0
      P1%NumDimensions = 0
      P1%Mass = 0.0d0
      P1%Width = 0.0d0
      P1%Charge = 0.0d0
      P1%AtomicNum = 0.0d0
      P1%Elmnt = 'XX'

      if (allocated(P1%Position)) then
         deallocate (P1%Position)
      end if

      if (allocated(P1%Momentum)) then
         deallocate (P1%Momentum)
      end if

      if (allocated(P1%Momentum2)) then
         deallocate (P1%Momentum2)
      end if

   end subroutine destroy_particle

!>
!!    Prints information about a particle.
!!    \param Particle particle to query
!!    \param Unit Fortran output unit
!!    @ingroup output
!<
   subroutine print_particle_info(Particle, iunit)
!   use QM_MM_Module, only: qcZQMMM, qcNumQM
      implicit none
      class(T_Particle), intent(in) :: Particle
      integer(kind=DefInt), intent(in) :: iunit
      integer(kind=DefInt) :: i

      write (IUnit, 1005) 'ID:       ', Particle%ParticleID
!   if(qcZQMMM) then
!      if (Particle%ParticleID <= qcNumQM) then
!         write(Iunit,1020) 'Type:     ','QM'
!      else
!         write(Iunit,1020) 'Type:     ','MM'
!      endif
!   endif

      write (IUnit, 1010) 'Mass:     ', Particle%Mass
      write (IUnit, 1010) 'Width:    ', Particle%Width
      write (IUnit, 1020) 'Element:  ', Particle%Elmnt
      write (IUnit, 1015) 'At. Num:  ', Particle%AtomicNum
      write (IUnit, 1005) 'Dimens:   ', Particle%NumDimensions

      write (IUnit, 1010) 'Position: ', Particle%get_pos(1)
      do i = 2, Particle%NumDimensions
         write (IUnit, 1000) Particle%get_pos(i)
      end do

      write (IUnit, 1010) 'Momentum: ', Particle%get_mom(1)
      do i = 2, Particle%NumDimensions
         write (IUnit, 1000) Particle%get_mom(i)
      end do

1000  format(16x, f15.8)
1005  format(a16, i5)
1010  format(a16, f15.8)
1015  format(a16, f4.0)
1020  format(a16, a5)
   end subroutine print_particle_info

!     Functions for position types
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_position_component(P1, i) result(pos)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Returns component i of particle's position
!!    \param P1 Particle to query
!!    \param IDim Cartesian dimension to query
!!    \return Position of particle in the IDim direction
!!    Scope: Public
!<
      class(T_Particle), intent(in) :: P1
      integer(kind=DefInt), intent(in) :: i
      real(kind=DefReal) :: pos

      if (i < 1 .or. i > size(P1%Position)) then
         call FMS_DieError('Particle%get_position_component: index out of range')
         pos = 0.0d0; return
      end if

      pos = P1%Position(i)
   end function get_position_component

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_position_vec(P1) result(pos)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Returns vector of particle's position
!!    \param P1 Particle to query
!!    \return Position vector of particle
!!    Scope: Public
!<
      class(T_Particle), intent(in) :: P1
      real(kind=DefReal) :: pos(P1%NumDimensions)

      pos = P1%Position

   end function get_position_vec

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_position_component(P1, i, NewPosition)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Sets the position of a particle in a trajectory.
!!
!!    \note sets the geomchanged flag to .true.
!!
!!    \param P1 Particle to move
!!    \param IDim Cartesian dimension to set
!!    Scope: ParticleModule, TrajectoryModule
!<
      class(T_Particle), intent(inout) :: P1
      integer(kind=DefInt), intent(in) :: i
      real(kind=DefReal), intent(in) :: NewPosition

      if (i < 1 .or. i > size(P1%Position)) then
         call FMS_DieError('Particle%set_position_component: index out of range ')
         return
      end if

      P1%Position(i) = NewPosition

   end subroutine set_position_component

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_position_vec(P1, NewPosition)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Sets the position of a particle in a trajectory.
!!    \note sets the geomchanged flag to .true.
!!
!!    \param P1 Particle to move
!!    \param IDim Cartesian dimension to set
!!    Scope: ParticleModule, TrajectoryModule
!<
      class(T_Particle), intent(inout) :: P1
      real(kind=DefReal), intent(in) :: NewPosition(:)

      if (size(P1%Position) /= size(NewPosition)) then
         call FMS_DieError('Particle%set_position_vec: wrong array size')
         return
      end if

      P1%Position = NewPosition
   end subroutine set_position_vec

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_momentum_component(P1, i) result(mom)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Returns vector of particle's position
!!    \param P1 Particle to query
!!    \return Position vector of particle
!!    Scope: Public
!<
      class(T_Particle), intent(in) :: P1
      integer, intent(in) :: i
      real(kind=DefReal) :: mom

      if (i < 1 .or. i > size(P1%Momentum)) then
         call FMS_DieError('Particle%get_momentum_component: index out of range')
         mom = 0.0d0; return
      end if

      mom = P1%Momentum(i)

   end function get_momentum_component

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function get_momentum_vec(P1) result(mom)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Returns vector of particle's position
!!    \param P1 Particle to query
!!    \return Position vector of particle
!!    Scope: Public
!<
      class(T_Particle), intent(in) :: P1
      real(kind=DefReal) :: mom(P1%NumDimensions)

      mom = P1%Momentum

   end function get_momentum_vec

   function get_momentum2_vec(P1) result(mom2)
      class(T_Particle), intent(in) :: P1
      real(kind=DefReal) :: mom2(P1%NumDimensions)

      mom2 = P1%Momentum2
   end function get_momentum2_vec

   subroutine set_momentum_component(P1, i, NewMomentum)
      class(T_Particle), intent(inout) :: P1
      integer(kind=DefInt), intent(in) :: i
      real(kind=DefReal), intent(in) :: NewMomentum

      if (i < 1 .or. i > size(P1%Momentum)) then
         call FMS_DieError('Particle%set_momentum_component: index out of range ')
         return
      end if

      P1%Momentum(i) = NewMomentum
   end subroutine set_momentum_component

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine set_momentum_vec(P1, new_momentum)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Sets the momentum of a particle in a trajectory.
!<
      class(T_Particle), intent(inout) :: P1
      real(kind=DefReal), intent(in) :: new_momentum(:)

      if (size(P1%Momentum) /= size(new_momentum)) then
         call FMS_DieError('Particle%set_momentum_vec: wrong array size')
         return
      end if

      P1%Momentum = new_momentum
   end subroutine set_momentum_vec

   subroutine set_momentum2_vec(P, new_mom2)
      class(T_Particle), intent(inout) :: P
      real(kind=DefReal), intent(in) :: new_mom2(:)

      if (size(P%Momentum) /= size(new_mom2)) then
         call FMS_DieError('Particle%set_momentum2_vec: wrong array size')
         return
      end if
      P%Momentum2 = new_mom2
   end subroutine set_momentum2_vec

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function vector_norm(v1) result(norm_out)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(kind=DefReal), intent(in) :: v1(:)
      real(kind=DefReal) :: norm_out

      norm_out = sqrt(dot_product(v1, v1))
   end function vector_norm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function unit_vector(v1) result(v2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(kind=DefReal), intent(in) :: v1(:)
      real(kind=DefReal) :: v2(size(v1))

      v2 = v1 / sqrt(dot_product(v1, v1))
   end function unit_vector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function cross(v1, v2) result(v3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real(kind=DefReal), intent(in) :: v1(3), v2(3)
      real(kind=DefReal) :: v3(3)

      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
   end function cross

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Distance_Particle(P1, P2) result(R)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Particle), intent(in) :: P1, P2
      real(kind=DefReal) :: R

      if (P1%NumDimensions /= P2%NumDimensions) then
         call FMS_DieError('Distance: Particles must have same dim!')
         R = 0.0d0; return
      end if

      R = vector_norm(P1%get_position_vec() - P2%get_position_vec())

   end function FMS_Distance_Particle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Angle_Particle(Pi, Pj, Pk) result(theta)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Particle), intent(in) :: Pi, Pj, Pk
      real(kind=DefReal) :: theta

      real(kind=DefReal) :: Xij(3), Xkj(3) ! normalized bond vectors

      if (Pi%NumDimensions /= Pj%NumDimensions) then
         call FMS_DieError('Angle: Particles must have same dim!')
         theta = 0.0d0; return
      end if

      if (Pi%NumDimensions /= Pk%NumDimensions) then
         call FMS_DieError('Angle: Particles must have same dim!')
         theta = 0.0d0; return
      end if

      Xij = unit_vector(Pi%get_pos() - Pj%get_pos())
      Xkj = unit_vector(Pk%get_pos() - Pj%get_pos())

      theta = acos(dot_product(Xij, Xkj))

   end function FMS_Angle_Particle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Dihedral_Particle(Pi, Pj, Pk, Pl) result(phi)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Particle), intent(in) :: Pi, Pj, Pk, Pl
      real(kind=DefReal) :: phi

      real(kind=DefReal), parameter :: tol_lo = 0.0872, & ! error if <   5 deg
                                       tol_hi = 305.0543 ! error if > 175 deg

      real(kind=DefReal), dimension(3) :: &
         Xji, Xjk, Xkj, Xkl, & ! bond vectors from i to j, etc
         Xt_kji, Xt_lkj ! unit vector normal to plane formed by {k,j,i}, {l,k,j}

      real(kind=DefReal) :: theta, & ! angle between bond vectors
                            cos_phi, sin_phi ! cosine and sine of angle angle
      ! between the planes {k,j,i} and
      ! {l,k,j}

! if the angle betwenn i,j,k or j,k,l are less than 5 deg. print a warning
      theta = FMS_Angle_Particle(Pi, Pj, Pk)
      if (theta < tol_lo .or. tol_hi < theta) then
         write (fmiOut, *) 'FMS_Dihedral_Particle : small angle between i,j,k'
      end if

      theta = FMS_Angle_Particle(Pj, Pk, Pl)
      if (theta < tol_lo .or. tol_hi < theta) then
         write (fmiOut, *) 'FMS_Dihedral_Particle : small angle between j,k,l'
      end if

      Xji = Pi%get_pos() - Pj%get_pos()
      Xkj = Pj%get_pos() - Pk%get_pos()
      Xjk = -Xkj
      Xkl = Pl%get_pos() - Pk%get_pos()

      Xt_kji = unit_vector(cross(Xjk, Xji))
      Xt_lkj = unit_vector(cross(Xkl, Xkj))

      cos_phi = dot_product(Xt_kji, Xt_lkj)
      sin_phi = dot_product(cross(Xt_kji, Xt_lkj), unit_vector(Xjk))

      phi = atan2(sin_phi, cos_phi)

   end function FMS_Dihedral_Particle

end module ParticleModule
