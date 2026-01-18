! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module RattleModule
! Implements the Rattle algorithm as outlined in :
!  R. Kutteh "RATTLE Recipe for General Holonomic Constraints: Angle and Torsion
!             Constraints"
!  CCP5, Newsletter No. 46 page 9
use GlobalModule
use TrajectoryModule
implicit none
private

public :: Rattle_ReadConstraints,        &
          Rattle_SetConstraints,         &
          Rattle_Constrain_Position,     &
          Rattle_Constrain_Momentum
public :: constraint, D_constraint, all_position_constrained

save
logical :: constraints_set = .false.

integer(kind=DefInt), public :: nconstraint  ! Number of constraints

! Enumerate type for constraints
! Types of constraints, size nconstraint
integer(kind=DefInt), public, allocatable :: cn_type_list(:)
integer,parameter    :: cn_bond     = 1, &
                        cn_angle    = 2, &
                        cn_dihedral = 3, &
                        cn_com      = 4      ! center of mass constraints

real    (kind=DefReal) :: cn_tol     = 1.d-8 ! convergence for constraints
integer (kind=DefInt)  :: cn_max_its = 10    ! max number of iterations

integer(kind=DefInt )  :: cn_max

! Internal coordinate constraints, (at most 4 atoms)
! List of atoms involved
integer(kind=DefInt ), allocatable, public :: cn_atom_list(:,:)

real(kind=DefReal), allocatable :: cn_value(:) ! Value of of the constraint

! Center of Mass constraints
integer(kind=DefInt ), allocatable :: cn_com_comp(:) ! which component x,y,z = 1,2,3

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine Rattle_ReadConstraints( NumQMMMs_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer(DefInt), optional :: NumQMMMs_in

integer (kind=DefInt) :: iUnit
integer :: ios

character(200)       :: constraint_file

integer, parameter   ::   cnMaxConstr = 1000, cnMaxAtoms = 1000

integer(kind=DefInt) :: NumBonds     = 0,   &
                        NumAngles    = 0,   &
                        NumDihedrals = 0,   &
                        NumCOMs      = 0,   &
                        numQMMMs     = 0,   &
                        MaxIts       = 100, & ! max number of iterations
                        nc, n, m,           & ! counter of constraints, counter
                        i                     ! counter for Cartesian component

integer(kind=DefInt) :: iBondPtcles    (         2,cnMaxConstr)=0, &
                        iAnglePtcles   (         3,cnMaxConstr)=0, &
                        iDihedralPtcles(         4,cnMaxConstr)=0, &
                        iCOMPtcles     (cnMaxAtoms,cnMaxConstr)=0

real (kind=DefReal)  :: toler = 0.0001

! work out the max number of atoms needed
integer(kind=DefInt) :: max_val(1)


! For reading the file
namelist/Constraints/ NumBonds,     iBondPtcles,     &
                      NumAngles,    iAnglePtcles,    &
                      NumDihedrals, iDihedralPtcles, &
                      NumCOMs,      iCOMPtcles,      &
                      MaxIts,Toler

if( allocated(cn_value      )) deallocate( cn_value     )
if( allocated(cn_type_list  )) deallocate( cn_type_list )
if( allocated( cn_atom_list )) deallocate( cn_atom_list )
if( allocated( cn_com_comp  )) deallocate( cn_com_comp  )


! Open the constraints file
constraint_file = trim(FMSWorkingDir)//'Constraints.dat'

constraints_set = .false.

open(newunit=iUnit, file=constraint_file, status='old', iostat=ios)

if(ios == 0 )then
   ! Read the namelist
   read(iUnit,Constraints)
endif

close(iUnit)

if( present(NumQMMMs_in) ) NumQMMMs = NumQMMMs_in


! work out the number of constraints
nconstraint = NumBonds + NumAngles + NumDihedrals + 3*NumCOMs + NumQMMMs*(NumQMMMS-1)/2

allocate( cn_value    (    nconstraint ) )
allocate( cn_type_list(    nconstraint ) )

! work out the max number of atoms
max_val = maxval( [ (count(iCOMPtcles(:,n)>0), n=1,NumCOMs) ] )
cn_max  = max( max_val(1), 4 )

allocate( cn_atom_list( cn_max, nconstraint ) )
allocate( cn_com_comp( nconstraint) )

cn_atom_list = 0
cn_com_comp  = 0

cn_value     = 0.d0
cn_atom_list = 0

nc = 0

! Store the bond constraints
do n = 1, NumBonds
   nc = nc + 1
   cn_type_list(    nc) = cn_bond
   cn_atom_list(1:2,nc) = iBondPtcles(:,n)
enddo

! Store the angle constraints
do n = 1, NumAngles
   nc = nc + 1
   cn_type_list(    nc) = cn_angle
   cn_atom_list(1:3,nc) = iAnglePtcles(:,n)
enddo

! Store the dihderal constraints
do n = 1, NumDihedrals
   nc = nc + 1
   cn_type_list(    nc) = cn_dihedral
   cn_atom_list(1:4,nc) = iDihedralPtcles(:,n)
enddo

! Store the CoM constraint
do n = 1, NumCOMs
do i = 1, 3
   nc = nc + 1
   cn_type_list(  nc) = cn_com
   cn_atom_list(:,nc) = ICOMPtcles(:,n)
   cn_com_comp (  nc) = i
enddo
enddo

! Store the QMMM constraints
do n = 1,   NumQMMMs-1
do m = n+1, NumQMMMs
   nc = nc + 1
   cn_type_list(    nc) = cn_bond
   cn_atom_list(1:2,nc) = [n,m]
enddo
enddo

! Copy values into module
cn_max_its = MaxIts
cn_tol     = Toler

end subroutine Rattle_ReadConstraints

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine Rattle_SetConstraints( T1 )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!!    Using the current geometry, set the constraint values
!!    \param T1 Current trajectory
!<
type(T_Trajectory), intent(in) :: T1

integer(kind=DefInt) :: nc

do nc = 1, nconstraint
   cn_value(nc) = constraint( T1, nc )
enddo

constraints_set = .true.

end subroutine Rattle_SetConstraints

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine Rattle_Constrain_Position( T1, T0, dt )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This subroutine modifies the positions of T1 to maintain the constraints
! T1 is T0 propagated by one time step
type(T_Trajectory)             :: T1         ! Trajectory at current time  ( t+dt)
type(T_Trajectory), intent(in) :: T0         ! Trajectory at previous step (  t  )
real(kind=DefReal), intent(in) :: dt         ! Time step

integer(kind=DefInt) :: nc,  &               ! constraint counter
                        its                  ! iteration counter

real(kind=DefReal) :: &
                   gamma,              &     ! How far to step along constraint deriv
!  time t    time t+dt    update t+dt
                   S1                        ! value of constraint

real(kind=DefReal),dimension(T1%NumDimensions) :: &
   inv_2M,                             &     ! 1/(2M)
!  time t    time t+dt    update t+dt
                   R1,      R1_new,    &     ! position
                   V1,      V1_new,    &     ! momentum
   dS0dX,         dS1dX                      ! derivative of constraints


if( .not. constraints_set )then
   call Rattle_ReadConstraints()
   call Rattle_SetConstraints(T0)
   constraints_set = .true.
endif

inv_2M = 1.d0 / ( 2.d0 * T0%get_mass()  )

do its = 1, cn_max_its

   if( all_position_constrained(T1) ) return

   do nc = 1, nconstraint

       R1   = T1%get_pos()
       V1   = T1%get_vel()
       S1   =   constraint   ( T1, nc )
      dS1dX = D_constraint   ( T1, nc )
      dS0dX = D_constraint   ( T0, nc )

      gamma  = (1.d0/dt**2) * S1 / &
              ( dot_product(  inv_2M*dS0dX, dS1dX ) )

      R1_new = R1 - dt**2 * gamma * inv_2M * dS0dX
      V1_new = V1 - dt    * gamma * inv_2M * dS0dX

      call T1%set_pos(R1_new)
      call T1%set_vel(V1_new )
   enddo

enddo

end subroutine Rattle_Constrain_Position

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine Rattle_Constrain_Momentum( T1, T0, dt )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This subroutine modifies the momentum of T1 to maintain the constraints
! T1 is T0 propagated by one time step
type(T_Trajectory)             :: T1         ! Trajectory at current time  ( t+dt)
type(T_Trajectory), intent(in) :: T0         ! Trajectory at previous step (  t  )
real(kind=DefReal), intent(in) :: dt         ! Time step

integer(kind=DefInt) :: nc,  &               ! constraint counter
                        its                  ! iteration counter

real  (kind=DefReal) :: eta                  ! How far to step along constraint deriv

real(kind=DefReal),dimension(T1%NumDimensions) :: &
   inv_2M,                             &     ! 1/(2M)
!  time t    time t+dt    update t+dt
                  V1,                  &     ! velocity
                  P1,       P1_new,    &     ! momenta
   dS0dX,         dS1dX                      ! derivative of constraints

if( .not. constraints_set )then
   call Rattle_ReadConstraints()
   call Rattle_SetConstraints(T0)
   constraints_set = .true.
endif

inv_2M = 1.d0 /( 2.d0 * T0%get_mass() )

do its = 1, cn_max_its

   if( all_momentum_constrained(T1) ) return

   do nc = 1, nconstraint
       V1   = T1%get_vel()
       P1   = T1%get_mom()
      dS1dX = D_constraint   ( T1, nc )
      dS0dX = D_constraint   ( T0, nc )

      eta  = (1.d0/dt) * dot_product(       V1    , dS1dX ) / &
                         dot_product( inv_2M*dS1dX, dS1dX )

      P1_new = P1 - dt * inv_2M * eta * dS1dX

      call T1%set_mom(P1_new)
   enddo
enddo

end subroutine Rattle_Constrain_Momentum

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                        Private below here

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                        Calculate constrains & derivatives

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function constraint( T1, nc ) result( const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! calculates the constraint value
!! eg. if a bond, dist(i,j) - R_ij = const
!! if 0 is the constrain is satisfied
!<
type(T_Trajectory)   , intent(in) :: T1
integer(kind=DefInt ), intent(in) :: nc    ! constraint number
real   (kind=DefReal)             :: const

integer(kind=DefInt) :: i,j,k,l           ! atom counters

if( nc<1 .or. nconstraint<nc )then
   write(fmiOUT,*) "nc          ", nc
   write(fmiOUT,*) "nconstraint ", nconstraint
   call FMS_DieError( "constraint : nc out of range ")
endif

select case( cn_type_list(nc) )
case(cn_bond)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   const = FMS_distance( T1, i, j ) - cn_value(nc)

case(cn_angle)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   k = cn_atom_list(3,nc)
   const = FMS_angle( T1, i, j, k ) - cn_value(nc)

case(cn_dihedral)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   k = cn_atom_list(3,nc)
   l = cn_atom_list(4,nc)
   const = FMS_Dihedral( T1, i, j, k, l ) - cn_value(nc)

case(cn_com)
   const = COM_constraint( T1, nc ) - cn_value(nc)

case default
   const = 0.0D0
   write(fmiOut,*) "cn_type_list(nc)", cn_type_list(nc)
   call FMS_DieError("Invalid cn_type_list")

end select

end function constraint

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function D_constraint( T1, nc ) result( D_const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! calculates the constraint value with repect to cartesian
!! eg. if a bond, dist(i,j) - R_ij = const
!! if 0 is the constrain is satisfied
!<
type(T_Trajectory)  , intent(in) :: T1
integer(kind=DefInt), intent(in) :: nc    ! constraint number

integer(kind=DefInt) :: i,j,k,l           ! atom counters
real   (kind=DefReal)             :: D_const(T1%NumDimensions)

if (nc < 1 .or. nconstraint < nc) then
   call FMS_DieError("constraint: nc out of range")
end if

D_const = 0.d0

select case( cn_type_list(nc) )
case(cn_bond)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   D_const = FMS_D_distance( T1, i,j )

case(cn_angle)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   k = cn_atom_list(3,nc)
   D_const = FMS_D_angle( T1, i,j,k )

case(cn_dihedral)
   i = cn_atom_list(1,nc)
   j = cn_atom_list(2,nc)
   k = cn_atom_list(3,nc)
   l = cn_atom_list(4,nc)
   D_const = FMS_D_Dihedral( T1, i,j,k,l )

case(cn_com)
   D_const = D_COM_constraint( T1, nc )

case default
   call FMS_DieError("Invalid cn_type_list in D_constraint")

end select

end function D_constraint

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function COM_constraint( T1, nc ) result( const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory)  , intent(in) :: T1
integer(kind=DefInt), intent(in) :: nc
real   (kind=DefReal)            :: const

real   (kind=DefReal) :: sum_mass
integer(kind=DefInt)  :: n, np, i

sum_mass = 0.d0
const    = 0.d0

i = cn_com_comp(nc)

do n = 1, cn_max
   np = cn_atom_list(n,nc)
   if( np==0 ) exit

   sum_mass = sum_mass + T1%get_mass( np )
   const    = const    + T1%get_mass( np ) * T1%get_pos(np, i)
enddo

const = const  / sum_mass

end function COM_constraint

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function D_COM_constraint( T1, nc ) result( vec )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory)  , intent(in) :: T1
integer(kind=DefInt), intent(in) :: nc
real   (kind=DefReal)            :: vec(T1%NumDimensions)

real   (kind=DefReal) :: sum_mass
integer(kind=DefInt)  :: n, np, i

sum_mass = 0.d0
vec      = 0.d0

i = cn_com_comp(nc)

do n = 1, cn_max
   np = cn_atom_list(n,nc)
   if( np==0 ) exit

   sum_mass = sum_mass + T1%get_mass( np )
   vec(3*(np-1)+i) =     T1%get_mass( np )
enddo

vec = vec  / sum_mass

end function D_COM_constraint

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function D_constraint_finite( T1, nc ) result( D_const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculates the derivative of the constraint by finite difference
! Useful for checking analytic expression for very complicated constaints
! which are worh deriving the constraint for
type(T_Trajectory)  , intent(in) :: T1
integer(kind=DefInt), intent(in) :: nc
real   (kind=DefReal)            :: D_const(T1%NumDimensions)

type(T_Trajectory)  :: T2        ! temp trajectory that will get moved
integer :: natom, n, &           ! number of atoms and counter
           i,        &           ! counter for Cartesian components
           nx                    ! counter for Cartesian positions

real(kind=DefReal) :: dx=1.d-5,& ! set size for finite disp
                      xold       ! Cartesian position from T1

natom = T1%NumParticles
call T2%create(natom, 1)

call T2%set_pos(T1%get_pos())

nx = 0
do n = 1, natom
do i = 1, 3
   nx = nx + 1

   xold = T2%get_pos(n, i)

   D_const(nx) = 0.d0

   ! step up
   call T2%set_pos(n, i,xold + dx)
   D_const(nx) = D_const(nx) + constraint( T2, nc )

   ! step down
   call T2%set_pos(n, i, xold - dx)
   D_const(nx) = D_const(nx) - constraint( T2, nc )

   D_const(nx) = D_const(nx)/(2.d0*dx)

   call T2%set_pos(n, i, xold)
enddo
enddo

call T2%destroy()

end function D_constraint_finite


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!             Constraint checks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function position_constrained( T1, nc ) result( satisfied )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory)  ,intent(in) :: T1
integer(kind=Defint),intent(in) :: nc
logical :: satisfied

satisfied = .false.
if( abs( constraint( T1, nc) ) < cn_tol ) satisfied = .true.

end function position_constrained

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function all_position_constrained( T1 )      result( satisfied )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory),intent(in) :: T1
logical                       :: satisfied

integer(kind=DefInt) :: nc

satisfied = .true.
do nc = 1, nconstraint
   satisfied = satisfied .and. position_constrained( T1, nc )
enddo

end function all_position_constrained

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function momentum_constrained( T1, nc ) result( satisfied )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory)  ,intent(in) :: T1
integer(kind=Defint),intent(in) :: nc
logical :: satisfied

real(kind=DefReal),dimension(T1%NumDimensions) :: &
         dRdX,  & ! velocity
         dSdX     ! derivative of constraint

dRdX = T1%get_vel()
dSdX = D_constraint   ( T1, nc )

satisfied = .false.
if( abs( dot_product( dRdX, dSdX ) ) < cn_tol ) satisfied = .true.

end function momentum_constrained

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function all_momentum_constrained( T1 )      result( satisfied )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory),intent(in) :: T1
logical                       :: satisfied

integer(kind=DefInt) :: nc

satisfied = .true.
do nc = 1, nconstraint
   satisfied = satisfied .and. momentum_constrained( T1, nc )
enddo

end function all_momentum_constrained

end module RattleModule
