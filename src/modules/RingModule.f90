module RingModule
implicit none
real(8),parameter ::  pi=4.0d0*datan(1.0d0)
type T_Ring

integer  :: n !< Number of beads
integer  :: natom !< Number of atoms
integer  :: ndim !< Dimensions per atom
integer  :: nNHC !< Number of Nose-Hoover chain variables
real (8),allocatable :: r(:,:,:)      !< ndim,natom,n Positions (bohr)
real (8),allocatable :: p(:,:,:)      !< ndim,natom,n Momenta (a.u.)
real (8),allocatable :: rNHC(:,:,:,:) !< nNHC,ndim,natom,n (bohr)
real (8),allocatable :: pNHC(:,:,:,:) !< nNHC,ndim,natom,n (bohr)
real (8),allocatable :: mNHC(:,:)     !< nNHC,n NHC masses for ring normal modes
real (8),allocatable :: F_ext(:,:,:)  !< ndim,natom,n External force (au)
real (8),allocatable :: m(:,:)        !< natom,n Mass of ring normal modes
real (8),allocatable :: m_inv(:,:)    !< natom,n Inverse mass of ring normal modes
real (8),allocatable :: omega(:)      !< normal mode frequencies
real (8) :: beta !< System inverse temperature (a.u.)
real (8) :: beta_n !< Ring inverse temperature (a.u.)
real (8) :: thermE !< Thermostat internal energy
real (8) :: time !< Current time
character(len=8) :: method !< Which method to use (RPMD, PA-CMD)
character(len=10) :: thermostat !< Which thermostat to use
real (8) :: tau0 !< Thermostat relaxation time of centroid mode
real (8) :: gamma_sq !< PA-CMD adiabaticity γ²; negative → Manolopoulos Ω
logical :: nm !< Are we in normal modes or not
logical :: zcom !< Do we include center of mass as degrees of freedom
end type T_Ring

public
private :: pi

type(T_Ring), allocatable :: FMS_rings(:)

interface assignment (=)
  module procedure ring_assign
end interface


contains

subroutine ring_assign(ring1,ring2)
implicit none
type (T_Ring) :: ring1, ring2
intent (in) ring2
intent (inout) ring1
logical Zsuccess

if ((ring1%n /= ring2%n).or.(ring1%ndim /= ring2%ndim).or.(ring1%natom /= ring2%natom).or.(ring1%nNHC /= ring2%nNHC)) then
  Zsuccess=ring_destroy(ring1)
  Zsuccess=ring_create(ring1,ring2%ndim,ring2%natom,ring2%n,ring2%nNHC)
endif

ring1%r          = ring2%r
ring1%p          = ring2%p
ring1%F_ext      = ring2%F_ext
ring1%m          = ring2%m
ring1%m_inv      = ring2%m_inv
ring1%omega      = ring2%omega
ring1%beta       = ring2%beta
ring1%beta_n     = ring2%beta_n
ring1%thermE     = ring2%thermE
ring1%time       = ring2%time
ring1%method     = ring2%method
ring1%thermostat = ring2%thermostat
ring1%rNHC       = ring2%rNHC
ring1%pNHC       = ring2%pNHC
ring1%mNHC       = ring2%mNHC
ring1%tau0       = ring2%tau0
ring1%gamma_sq   = ring2%gamma_sq
ring1%nm         = ring2%nm
ring1%zcom       = ring2%zcom

end subroutine ring_assign


function ring_create(ring,ndim,natom,n,nNHC) result (ZSuccess)
implicit none
integer ndim,natom,n,nNHC
type(T_Ring) ring
logical ZSuccess

integer istat

ZSuccess=.true.
Zsuccess=ring_destroy(ring)
istat=0

allocate(ring%r(ndim,natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%p(ndim,natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%F_ext(ndim,natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%m(natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%m_inv(natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%omega(n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%rNHC(nNHC,ndim,natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%pNHC(nNHC,ndim,natom,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.
allocate(ring%mNHC(nNHC,n),stat=IStat) ; if(IStat /= 0) ZSuccess=.false.

ring%r=0.0d0
ring%p=0.0d0
ring%F_ext=0.0d0
ring%m=0.0d0
ring%m_inv=0.0d0
ring%omega=0.0d0
ring%rNHC=0.0d0
ring%pNHC=0.0d0
ring%mNHC=0.0d0

ring%ndim=ndim
ring%natom=natom
ring%n=n
ring%nNHC=nNHC

ring%beta=0.0d0
ring%beta_n=0.0d0
ring%thermE=0.0d0
ring%time=0.0d0
ring%method=''
ring%thermostat=''
ring%gamma_sq=0.005d0
ring%nm=.false.
ring%zcom=.false.

if (.not. ZSuccess) then
   write(*,*)' Memory allocation error in ring_create'
   stop
endif

end function  ring_create

function ring_destroy(ring) result (ZSuccess)
implicit none
type(T_Ring) ring
logical ZSuccess

integer istat

ZSuccess=.true.
istat=0

ring%n=0
ring%ndim=0
ring%natom=0

if(allocated(ring%r)) deallocate(ring%r,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%p))  deallocate(ring%p,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%F_ext))  deallocate(ring%F_ext,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%m))  deallocate(ring%m,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%m_inv))  deallocate(ring%m_inv,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%omega))  deallocate(ring%omega,stat=istat) ; if (istat /= 0) ZSuccess=.false.
if(allocated(ring%rNHC)) deallocate(ring%rNHC,stat=istat) ; if (istat /=0) ZSuccess=.false.
if(allocated(ring%pNHC)) deallocate(ring%pNHC,stat=istat) ; if (istat /=0) ZSuccess=.false.
if(allocated(ring%mNHC))  deallocate(ring%mNHC,stat=istat) ; if (istat /= 0) ZSuccess=.false.

if (.not. ZSuccess) then
   write(*,*)' Memory deallocation error in Ring_destroy'
   stop
endif

end function ring_destroy

subroutine FMS_ring_setup(nTraj)
  integer, intent(in) :: nTraj
  if (allocated(FMS_rings)) deallocate(FMS_rings)
  allocate(FMS_rings(nTraj))
end subroutine FMS_ring_setup

function ring_write(ring,iunit) result(ZSuccess)
implicit none
type(T_Ring) ring
integer iunit
logical ZSuccess

ZSuccess=.false.

if (ring%nm) then
  write(*,*)'Error: normal modes in ring_write'
  ZSuccess=.false.
  return
end if
!Note, this subroutine writes out only non-redundant information, i.e. m_inv is omitted. Therefore, ring_read must recreate these variables.

1    format(i11,      34x,a)    !4byte integer 
2    format(ES23.15E3, 22x,a)   !double precision real with label
3    format(3(1x,ES23.15E3))    !double precision real array
4    format(1x,l1,43x,a)        !logical format with label
5    format(a10,35x,a)          !String with label
write(iunit,2)ring%time,' /Time (a.u.)'
write(iunit,1)ring%n,' /Number of Beads'
write(iunit,1)ring%natom,' /Number of Atoms'
write(iunit,1)ring%ndim,' /Number of Dimensions'
write(iunit,2)ring%beta,' /Inverse temperature (a.u.)'
write(iunit,5)ring%method,' /Method'
write(iunit,5)ring%thermostat,' /Thermostat'
write(iunit,1)ring%nNHC,' /Number of Nose-Hoover Chain variables'
write(iunit,2)ring%tau0,' /Centroid thermostat relaxation constant (a.u.)'
write(iunit,2)ring%gamma_sq,' /PA-CMD adiabaticity gamma^2 (negative = Manolopoulos)'
write(iunit,2)ring%thermE,' /Thermostat energy (a.u.)'
write(iunit,4)ring%zcom,' /Treat COM as degrees of freedom?'
write(iunit,*)'Masses:'
write(iunit,3)ring%m
write(iunit,*)'Positions:'
write(iunit,3)ring%r
write(iunit,*)'Momenta:'
write(iunit,3)ring%p
write(iunit,*)'External forces:'
write(iunit,3)ring%F_ext
write(iunit,*)'NHC Positions:'
write(iunit,3)ring%rNHC
write(iunit,*)'NHC Momenta:'
write(iunit,3)ring%pNHC
write(iunit,*)'NHC Masses:'
write(iunit,3)ring%mNHC
write(iunit,*)"###"

ZSuccess=.true.
return
end function ring_write

function ring_read(ring,iunit) result(ZSuccess)
implicit none
type(T_Ring) ring
integer iunit
logical ZSuccess
ZSuccess=.false.

if (ring%nm) then
  write(*,*)'Error: normal modes in ring_read'
  ZSuccess=.false.
  return
end if

read(iunit,*)ring%time
read(iunit,*)ring%n
read(iunit,*)ring%natom
read(iunit,*)ring%ndim
read(iunit,*)ring%beta
read(iunit,*)ring%method
read(iunit,*)ring%thermostat
read(iunit,*)ring%nNHC
read(iunit,*)ring%tau0 ! NHC fictitious mass
read(iunit,*)ring%gamma_sq
read(iunit,*)ring%thermE
read(iunit,*)ring%zcom
read(iunit,*)
read(iunit,*)ring%m
read(iunit,*)
read(iunit,*)ring%r
read(iunit,*)
read(iunit,*)ring%p
read(iunit,*)
read(iunit,*)ring%F_ext
read(iunit,*)
read(iunit,*)ring%rNHC
read(iunit,*)
read(iunit,*)ring%pNHC
read(iunit,*)
read(iunit,*)ring%mNHC

!Fill in remaining variables
call ring_setconsts(ring,ring%m(:,1),ring%beta,ring%method,ring%thermostat,ring%tau0,ring%zcom,gamma_sq_in=ring%gamma_sq)

ZSuccess=.true.
return
end function ring_read

subroutine ring_Etot(ring,V_func,ETot)
use FFTModule
implicit none
type (T_Ring) ring
real(8) ETot
real(8), dimension(ring%n) :: V_vec
real(8) V_ext_tot
integer ibead

interface
  function V_func(r_vec,ndim,natom,n)
  implicit none
  integer,intent(in) ::  ndim,natom,n
  real(8), dimension(ndim,natom,n), intent(in) :: r_vec
  real(8), dimension(n) :: V_func
  end function V_func
end interface

V_vec=V_func(ring%r,ring%ndim,ring%natom,ring%n)
V_ext_tot=0.0d0
do ibead=1,ring%n
  V_ext_tot=V_ext_tot+V_vec(ibead)
enddo
if (ring%method(1:7)=='GAUSSAV') then
  V_ext_tot=V_ext_tot/dble(ring%n)
endif
Etot=ring_V(ring)+ring_T(ring)+V_ext_tot+ring%thermE

end subroutine ring_Etot

subroutine ring_Eout(ring,V_func,Force_func,iunit)
implicit none
type (T_Ring) ring
integer iunit
interface
  function V_func(r_vec,ndim,natom,n)
  implicit none
  integer,intent(in) ::  ndim,natom,n
  real(8), dimension(ndim,natom,n), intent(in) :: r_vec
  real(8), dimension(n) :: V_func
  end function V_func

  function Force_func(r_vec,ndim,natom,n)
  implicit none
  integer, intent(in) :: ndim,natom,n
  real(8), dimension(ndim,natom,n), intent(in) :: r_vec
  real(8), dimension(ndim,natom,n):: Force_func
  end function Force_func
end interface
real(8), dimension(ring%ndim,ring%natom,ring%n) :: F
real(8), dimension(ring%ndim,ring%natom) :: r_cent
real(8), dimension(ring%n) :: V_vec
real(8) V_ext_tot,V_ring,Tkin,T0_ring,E_ring_tot,V_ring_nm,T_ring_nm,V_system,T_system
real(8) dn_inv
integer ibead,iatom,idm,ndof
dn_inv = 1.0d0/dble(ring%n)

!Calculate total energy of ring (conserved quantity)
V_vec=V_func(ring%r,ring%ndim,ring%natom,ring%n)
V_ext_tot=0.0d0
do ibead=1,ring%n
  V_ext_tot=V_ext_tot+V_vec(ibead)
enddo
if (ring%method(1:7)=='GAUSSAV') then
  V_ext_tot=V_ext_tot/dble(ring%n)
endif

V_ring=ring_V(ring)
Tkin=ring_T(ring)
T0_ring=ring_T0(ring)
E_ring_tot=ring_V(ring)+ring_T(ring)+V_ext_tot+ring%thermE
!Calculate ring energies in normal mode representation
call nmtran_forward(ring)
V_ring_nm=ring_V(ring)
T_ring_nm=ring_T(ring)
call nmtran_backward(ring)

!Calculate system potential energy (Eq. 9)
V_system=V_ext_tot*dn_inv

!Calculate system kinetic energy estimator (Eq. 11)
if (ring%zcom) then
  ndof=ring%ndim*ring%natom
else
  ndof=ring%ndim*ring%natom - ring%ndim
endif

F=Force_func(ring%r,ring%ndim,ring%natom,ring%n)
T_system=0.0d0
r_cent=ring_centvec(ring%r,ring%ndim,ring%natom,ring%n)
do ibead=1,ring%n
  do iatom=1,ring%natom
    do idm=1,ring%ndim
      T_system=T_system-(ring%r(idm,iatom,ibead)-r_cent(idm,iatom))*F(idm,iatom,ibead)
    enddo
  enddo
enddo

T_system=0.5d0*T_system*dn_inv
T_system=T_system + dble(ndof)*0.5d0/ring%beta

!Print out
write(iunit,100)ring%time,V_ring,Tkin,T0_ring,E_ring_tot,V_system,T_system,ring%thermE

100 format(ES13.5E3,1x,8(ES23.15E3,1x))
end subroutine ring_Eout

subroutine ring_Rout(ring,iunit)
!Prints out ring positions.
! It is anticipated this would only be used for low dimensional systems
implicit none
type (T_Ring) :: ring
integer iunit,ibead

do ibead=1,ring%n
  write(iunit,'(ES23.15E3,1x,I6,1x,8(ES23.15E3,1x))')ring%time,ibead,ring%r(:,:,ibead)
enddo
end subroutine ring_Rout

subroutine ring_Pout(ring,iunit)
!Prints out ring momenta.
! It is anticipated this would only be used for low dimensional systems
implicit none
type (T_Ring) :: ring
integer iunit,ibead

do ibead=1,ring%n
  write(iunit,'(ES23.15E3,1x,I6,1x,8(ES23.15E3,1x))')ring%time,ibead,ring%p(:,:,ibead)
enddo
end subroutine ring_Pout

subroutine ring_centRout(ring,iunit)
!Prints out centroid position.
! It is anticipated this would only be used for low dimensional systems
implicit none
type (T_Ring) :: ring
integer iunit
real(8),dimension(ring%ndim,ring%natom) :: r_cent

r_cent=ring_centvec(ring%r,ring%ndim,ring%natom,ring%n)
write(iunit,'(9(ES23.15E3,1x))')ring%time,r_cent
end subroutine ring_centRout

subroutine ring_centPout(ring,iunit)
implicit none
type (T_Ring) :: ring
integer iunit
real(8),dimension(ring%ndim,ring%natom) :: p_cent

p_cent=ring_centvec(ring%p,ring%ndim,ring%natom,ring%n)
write(iunit,'(9(ES23.15E3,1x))')ring%time,p_cent
end subroutine ring_centPout


function ring_centvec(vec,ndim,natom,n) result(centvec)
!Calculates a centroid vector from a ring vector
implicit none
integer ndim,natom,n
real(8),dimension(ndim,natom,n) :: vec
real(8),dimension(ndim,natom) :: centvec

real(8) dn_inv
integer ibead

dn_inv=1.0d0/dble(n)

centvec(:,:)=0.0d0

do ibead=1,n
  centvec(:,:)=centvec(:,:)+vec(:,:,ibead)
enddo
centvec(:,:)=centvec(:,:)*dn_inv

end function ring_centvec

function ring_T(ring)
!Calculates ring kinetic energy in either bead or normal mode representation
implicit none
type (T_ring) :: ring
real(8) ring_T
integer k,ibead,iatom,idm
logical orignm

ring_T=0.0d0

orignm=ring%nm
if ((ring%method=='PA-CMD'.or.ring%method(1:7)=='GAUSSAV').and.(.not.orignm)) then
  call nmtran_forward(ring) !ensure we are in normal mode rep
endif

if (ring%nm) then
!Loop over normal modes (Eq. 19)
  do k=1,ring%n
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        ring_T=ring_T+ring%m_inv(iatom,k)*ring%p(idm,iatom,k)*ring%p(idm,iatom,k)
      enddo
    enddo
  enddo
else
!Loop over beads (Eq. 5)
  do ibead=1,ring%n
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        ! TL: but your ring%m_inv is the inverse of the normal mode mass, not the bead mass. 
        ! So this is not correct. It should be ring%m_inv(iatom,ibead) = 1.0d0/m(iatom) for all beads. 
        ! But I will leave it as is for now.
        ring_T=ring_T+ring%m_inv(iatom,ibead)*ring%p(idm,iatom,ibead)*ring%p(idm,iatom,ibead)
      enddo
    enddo
  enddo
endif
ring_T=0.5d0*ring_T

!Make sure we are back in the original representation
if (.not.orignm) then
  call nmtran_backward(ring)
endif

end function ring_T

function ring_T0(ring)
!Calculates kinetic energy of zero-frequency mode in normal mode representation
implicit none
type (T_ring) :: ring
real(8) ring_T0
integer iatom,idm
logical orignm

orignm=ring%nm
if (.not.orignm) then
  call nmtran_forward(ring) !ensure we are in normal mode rep
endif

ring_T0=0.0d0
do iatom=1,ring%natom
  do idm=1,ring%ndim
    ring_T0=ring_T0+ring%m_inv(iatom,1)*ring%p(idm,iatom,1)*ring%p(idm,iatom,1)
  enddo
enddo
ring_T0=0.5d0*ring_T0

!Make sure we are back in the original representation
if (.not.orignm) then
  call nmtran_backward(ring)
endif

end function ring_T0

function ring_V(ring)
!Calculates ring potential in either bead or normal mode representation
implicit none
type (T_ring) :: ring
real(8) ring_V

real(8) omegafac,m_omegafac
integer k,ibead,ibead_nei,iatom,idm
logical orignm

ring_V=0.0d0

!Fix for 1bead, zero temperature:
if (ring%n.eq.1) return

orignm=ring%nm
if ((ring%method=='PA-CMD'.or.ring%method(1:7)=='GAUSSAV').and.(.not.orignm)) then
  call nmtran_forward(ring) !ensure we are in normal mode rep
endif

if (ring%nm) then
!Run over normal modes (Eq. 19)
  do k=1,ring%n
    omegafac=0.5d0*ring%omega(k)*ring%omega(k)
    do iatom=1,ring%natom
      m_omegafac=ring%m(iatom,k)*omegafac 
      do idm=1,ring%ndim
        ring_V=ring_V+m_omegafac*ring%r(idm,iatom,k)*ring%r(idm,iatom,k)
      enddo
    enddo
  enddo
else
  if (ring%method=='PA-CMD'.or.ring%method(1:7)=='GAUSSAV') then
    write(6,*)'can only evaluate Vring in normal mode rep'
    stop
  endif
  omegafac=0.5d0/(ring%beta_n*ring%beta_n)
  !Run over bead links (Eq. 5)
  do ibead=1,ring%n
    ibead_nei=ring_nei(ibead,ring%n)
    do iatom=1,ring%natom
      m_omegafac=ring%m(iatom,ibead)*omegafac
      do idm=1,ring%ndim
        ring_V=ring_V + m_omegafac*(ring%r(idm,iatom,ibead) - ring%r(idm,iatom,ibead_nei))**2 
      enddo
    enddo
  enddo
endif

if (.not.orignm) then
  call nmtran_backward(ring)
endif

contains
function ring_nei(ibead,n)
!Finds lower index of neighbouring bead with periodicity
implicit none
integer ibead,n,ring_nei
if (ibead.lt.1.or.ibead.gt.n) then
  write(*,*)'Error in ring_nei. ibead out of bounds'
  stop
endif
ring_nei=ibead-1
if (ibead.eq.1) ring_nei=n
end function ring_nei

end function ring_V

subroutine ring_init(ring,RefGeom,m,beta,method,InitCond,thermostat,tau0,zcom,idum,omega_in,widths_in)
use ThermoModule
implicit none
type (T_Ring) ring
real(8) :: RefGeom(ring%ndim,ring%natom)
real(8) :: m(ring%natom)
character(len=8) method
character(len=10) InitCond
integer idum
real(8) :: beta,tau0,thermE
character(len=10) thermostat
logical zcom
real(8),intent(in),optional :: omega_in
real(8),intent(in),optional :: widths_in(ring%natom)

integer ibead,k,iatom,idm

!Initialize thermostat (random seed)
call thermo_init(idum)
!Set constants
call ring_setconsts(ring,m,beta,method,thermostat,tau0,zcom,omega_in,widths_in)

!Copy reference geometry into beads
do ibead=1,ring%n
  ring%r(:,:,ibead)=RefGeom(:,:)
enddo

select case(InitCond)
case ("MB")
!Sample momenta from the Maxwell-Boltzmann distribution in normal mode Rep
  call nmtran_forward(ring)
  call thermo_MBDist(ring%ndim,ring%natom,ring%p(:,:,1),ring%m(:,1),ring%beta_n,ring%zcom) !only centroid gets center of mass projected out
!In ring normal modes so dont project out vcom
  do k=2,ring%n
    call thermo_MBDist(ring%ndim,ring%natom,ring%p(:,:,k),ring%m(:,k),ring%beta_n,.true.) !In ring normal modes so dont project out vcom
  enddo

!Sample thermostat variables for all modes and degrees of freedom
  select case(thermostat)
  case("NHC-L")
    ring%thermE=0.0d0
    do k=1,ring%n
      do iatom=1,ring%natom
        do idm=1,ring%ndim
          call thermo_NHC_MBDist(ring%beta_n,ring%nNHC,ring%pNHC(:,idm,iatom,k),ring%mNHC(:,k),thermE)
          ring%thermE=ring%thermE+thermE
        enddo
      enddo
    enddo

  case("NHC-G")
    !centroid mode is globally thermostatted, so chain only coupled to 1 degree of
    !freedom
    ring%thermE=0.0d0
    call thermo_NHC_MBDist(ring%beta_n,ring%nNHC,ring%pNHC(:,1,1,1),ring%mNHC(:,1),thermE)
    ring%thermE=ring%thermE+thermE

    !Remaining modes (local)
    do k=2,ring%n
      do iatom=1,ring%natom
        do idm=1,ring%ndim
          call thermo_NHC_MBDist(ring%beta_n,ring%nNHC,ring%pNHC(:,idm,iatom,k),ring%mNHC(:,k),thermE)
          ring%thermE=ring%thermE+thermE
        enddo
      enddo
    enddo
  end select
  call nmtran_backward(ring)

end select

end subroutine ring_init

subroutine ring_setconsts(ring,m,beta,method,thermostat,tau0,zcom,omega_in,widths_in,gamma_sq_in)
implicit none
type (T_Ring) ring
real(8) :: m(ring%natom)
real(8) :: beta,tau0
character(len=8) method
character(len=10) thermostat
logical zcom
real(8),intent(in),optional :: omega_in !Desired frequency of ring modes for GAUSSAV
real(8),intent(in),optional :: widths_in(ring%natom) !Desired widths for GAUSSAV
real(8),intent(in),optional :: gamma_sq_in !< PA-CMD γ²; negative → Manolopoulos Ω; default 0.005
integer ibead,k
real(8) omega,sigma,two_omega_n,dn_inv
real(8) :: gamma_sq_local
logical :: use_manolopoulos
real(8) :: widths(ring%natom)

dn_inv=1.0d0/dble(ring%n)
! Resolve adiabaticity: negative gamma_sq_in selects Manolopoulos
if (present(gamma_sq_in)) then
  gamma_sq_local = gamma_sq_in
else
  gamma_sq_local = 0.005d0
end if
use_manolopoulos = (gamma_sq_local < 0.0d0)
ring%gamma_sq = gamma_sq_local
!Set constants
ring%beta=beta
ring%beta_n=ring%beta*dn_inv
ring%thermostat=thermostat
ring%tau0=tau0
ring%zcom=zcom
ring%method=method

!Calculate normal mode frequencies as if masses were physical in all modes
two_omega_n=2.0d0/ring%beta_n

do k=0,ring%n-1
  ring%omega(k+1)=two_omega_n*dsin(dble(k)*pi*dn_inv)
enddo

select case(method)
case('RPMD')
  !All beads/modes have physical mass
  do ibead=1,ring%n
    ring%m(:,ibead)=m(:)
    ring%m_inv(:,ibead)=1.0d0/m(:)
  enddo

case('PA-CMD')
  !centroid takes physical masses
  ring%m(:,1)=m(:)
  ring%m_inv(:,1)=1.0d0/m(:)

  !All other modes have a constant frequency, so scale masses:
  !  m_k = gamma^2 * m * lambda_k,  lambda_k = 4*sin^2(k*pi/n)
  !  All non-centroid modes then share the common frequency Omega.
  if (ring%n.gt.1) then
    if (use_manolopoulos) then
      ! Manolopoulos choice: Omega = n^(n/(n-1)) / beta
      omega = dble(ring%n)**(dble(ring%n)/dble(ring%n-1)) / beta

      do k=2,ring%n
        sigma=ring%omega(k)/omega
        ring%m(:,k)=m(:)*sigma*sigma
        ring%m_inv(:,k)=1.0d0/ring%m(:,k)
        ring%omega(k)=omega             ! store Omega for thermostat mass computation
        !write(*,*) '  sigma = ',sigma,'  omega_k = ',ring%omega(k), '  PA-CMD: Manolopoulos Omega = ',omega
        !write(*,*) '  m_k = ',ring%m(:,k)
      enddo
    else
      ! Voth default (gamma^2 = 0.005):
      !   Omega = omega_P / gamma,  omega_P = n/beta  →  Omega = n / (beta * sqrt(gamma^2))
      omega = dsqrt(dble(ring%n)) / (beta * dsqrt(gamma_sq_local))

      do k=2,ring%n
        sigma=beta*beta*ring%omega(k)*ring%omega(k)/(dble(ring%n)*dble(ring%n))
        ring%m(:,k)=gamma_sq_local*m(:)*sigma   ! m_k = m*(omega_k/Omega)^2 = gamma^2*m*lambda_k
        ring%m_inv(:,k)=1.0d0/ring%m(:,k)
        !ring%omega(k)=omega             ! store Omega for thermostat mass computation
      enddo
    end if
    !write(0,*)'PA-CMD: Omega = ',omega,'  gamma = ',dsqrt(abs(gamma_sq_local))
  end if

  !write(0,*)'max ring mode frequency: ',ring%omega(ring%n/2 + 1)

case('GAUSSAVE','GAUSSAVF')
  if (.not.present(omega_in)) then
    write(6,*)'Error in ring_setconsts. Need input ring frequency'
    stop
  endif
  if (.not.present(widths_in)) then
    write(6,*)'Error in ring_setconsts. Need input Gaussian widths'
    stop
  endif

  !centroid takes physical masses
  ring%m(:,1)=m(:)
  ring%m_inv(:,1)=1.0d0/m(:)
!  ring%beta_n=1.0d0            !Set effective temperature as constant so that width only determined from masses
!  ring%beta=ring%beta_n*dble(ring%n)

  !All other modes have a constant frequency and mass chosen to give correct
  !gaussian width
  !Scale widths according to number of beads, since centroid is frozen:
  widths(:)=widths_in(:)*sqrt(dble(ring%n)/dble(ring%n-1))

  do k=1,ring%n
    ring%m(:,k)=1.0d0/(ring%beta_n*(omega_in*widths(:))**2)
!    ring%m(:,k)=1.0d0/(dsqrt(dble(ring%n))*ring%beta_n*(omega_in*widths_in(:))**2)
    ring%m_inv(:,k)=1.0d0/ring%m(:,k)
    ring%omega(k)=omega_in
  enddo

case default
  write(6,*)'Error, unknown method:',method
  stop

end select


!Set thermostat masses
select case(thermostat)
case('NHC-L', 'NHC-G')
!For each mode, all chains have same mass
!Set centroid NHC masses
  ring%mNHC(:,1)=tau0*tau0/ring%beta
  ring%mNHC(1,1)=ring%natom*tau0*tau0/ring%beta
  !ring%mNHC(:,1)=4.0d0*tau0*tau0/ring%beta_n

!Set NHC masses for other modes
do k=2,ring%n
  ring%mNHC(:,k)=1.0d0/(ring%beta_n*ring%omega(k)*ring%omega(k))
  ring%mNHC(1,k)=ring%natom/(ring%beta_n*ring%omega(k)*ring%omega(k))
end do

end select
end subroutine ring_setconsts

subroutine ring_modesample(ring,omega0,skipr0,skipp0)
!This subrouine samples the thermal density of the ring normal modes within a harmonic approximation to the external potential. For now, we assume this is a 1D HO, so a single external potential frequency. We should make this more general by taking in a Hessian and then transforming to normal modes.
!Allow the user to set a custom position and momentum for ring centroid.
use ThermoModule
implicit none
type (T_Ring) :: ring
real(8) omega0(ring%ndim,ring%natom)
logical,optional ::skipr0, skipp0
logical skipr0here,skipp0here

real(8) omega,sigma
integer kstart,k,iatom,idm

if (present(skipr0)) then
  skipr0here = skipr0
else
  skipr0here = .false.
endif

if (present(skipp0)) then
  skipp0here = skipp0
else
  skipp0here = .false.
endif

call nmtran_forward(ring)
!Sample positions
kstart=1
if (skipr0here) then !don't resample mode 0
  kstart=2
endif
do k=kstart,ring%n
  do iatom=1,ring%natom
    do idm=1,ring%ndim
      omega=dsqrt(omega0(idm,iatom)*omega0(idm,iatom) + ring%omega(k)*ring%omega(k))
      sigma=1.0d0/(omega*dsqrt(ring%beta_n*ring%m(iatom,k)))
      ring%r(idm,iatom,k)=thermo_NormDist(sigma)
    enddo
  enddo
enddo

!Sample momenta
if (.not. skipp0here) then
  call thermo_MBDist(ring%ndim, ring%natom, ring%p(:,:,1), ring%m(:,1), ring%beta_n, ring%zcom)
endif
!sample other modes
do k=2,ring%n
  call thermo_MBDist(ring%ndim,ring%natom,ring%p(:,:,k),ring%m(:,k),ring%beta_n,.true.) !In ring normal modes so dont project out vcom
enddo

call nmtran_backward(ring)

end subroutine ring_modesample

subroutine ringprop(ring,dt,force_func,skip0,skipt0)
implicit none
type (T_Ring) :: ring
real(8) dt
logical,optional:: skip0,skipt0
logical skip0here,skipt0here
interface
  function Force_func(r_vec,ndim,natom,n)
  implicit none
  integer, intent(in) ::  ndim,natom,n
  real(8), dimension(ndim,natom,n), intent(in) :: r_vec
  real(8), dimension(ndim,natom,n):: Force_func
  end function Force_func
end interface

if (present(skip0)) then
  skip0here=skip0
else
  skip0here=.false.
endif

if (present(skipt0)) then
  skipt0here=skipt0
else
  skipt0here=.false.
endif

!Update Force (strictly not necessary except for 1st timestep)
call ring_updateforce(ring,Force_func)

!Thermostat velocities
call ring_thermostat(ring,0.5d0*dt,skipt0here)


!Propagate momenta 1/2 step using current external forces
call ring_momupdate(ring,0.5d0*dt,skip0here)

!Free ring propagation
call ring_freeprop(ring,dt,skip0here)


!Update Force
call ring_updateforce(ring,Force_func)

!Propagate momenta 1/2 step using new external forces
call ring_momupdate(ring,0.5d0*dt,skip0here)

!Thermostat velocities
call ring_thermostat(ring,0.5d0*dt,skipt0here)


ring%time=ring%time+dt
end subroutine ringprop

subroutine ring_updateforce(ring,Force_func)
implicit none
type (T_Ring) :: ring
interface
  function Force_func(r_vec,ndim,natom,n)
  implicit none
  integer, intent(in) ::  ndim,natom,n
  real(8), dimension(ndim,natom,n), intent(in) :: r_vec
  real(8), dimension(ndim,natom,n):: Force_func
  end function Force_func
end interface
ring%F_ext=Force_func(ring%r,ring%ndim,ring%natom,ring%n)
return
end subroutine ring_updateforce

subroutine ring_thermostat(ring,dt,skip0)
use ThermoModule
implicit none
type(T_Ring) ring
real(8) dt,thermE
integer iatom,idm,k,kstart
!real(8) p(ring%ndim,ring%natom), m(ring%natom)
real(8) tau
logical nm_orig,skip0

!Will want some logic to select which thermostat (or no thermostat if rpmd)
select case (ring%thermostat)
  case ("NONE")

  case ("NHC-L")
  nm_orig=ring%nm
  
  call nmtran_forward(ring)
  ring%thermE=0.0d0 !Reset thermostat energy because it is accumulated below
  kstart=1
  if (skip0) kstart=2
  do k=kstart,ring%n
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        call thermo_NHC_local(dt,ring%nNHC,ring%p(idm,iatom,k),ring%m(iatom,k),&
          ring%beta_n,ring%rNHC(:,idm,iatom,k),ring%pNHC(:,idm,iatom,k),ring%mNHC(:,k),thermE)
        ring%thermE=ring%thermE+thermE
      enddo
    enddo
  enddo

  if (.not.nm_orig) call nmtran_backward(ring)

  
  case ("NHC-G")
  nm_orig=ring%nm
  call nmtran_forward(ring)
  ring%thermE=0.0d0 !Reset thermostat energy because it is accumulated below

  if (.not.skip0) then
!Centroid mode thermostatted globally
    call thermo_NHC_global(dt,ring%nNHC,ring%ndim,ring%natom,ring%p(:,:,1),ring%m(:,1),&
      ring%beta_n,ring%rNHC(:,1,1,1),ring%pNHC(:,1,1,1),ring%mNHC(:,1),thermE)
    ring%thermE=ring%thermE+thermE
  endif
!Remaining modes
  do k=2,ring%n
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        call thermo_NHC_local(dt,ring%nNHC,ring%p(idm,iatom,k),ring%m(iatom,k),&
          ring%beta_n,ring%rNHC(:,idm,iatom,k),ring%pNHC(:,idm,iatom,k),ring%mNHC(:,k),thermE)
        ring%thermE=ring%thermE+thermE
      enddo
    enddo
  enddo
  if (.not.nm_orig) call nmtran_backward(ring)

  case ("PILE-L")
  nm_orig=ring%nm
  call nmtran_forward(ring)
!Zero-frequency mode:
!  m(:)=ring%m(:,1)!  p(:,:)=ring%p(:,:,1)
  tau=ring%tau0/dt !Note dt is whatever this subroutine was called with, which may be half a regular timestep
  kstart=1
  if (.not.skip0) then
    do iatom=1,ring%natom 
      do idm=1,ring%ndim
        call thermo_bussi_local(1, 1, ring%p(idm:idm,iatom:iatom,1), ring%m(iatom:iatom,1), &
          ring%beta_n, tau, ring%thermE, .true.)
      enddo
    enddo
  endif
!Other modes
  do k=2,ring%n
    tau=1.0d0/(dt*2.0d0*ring%omega(k))
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        call thermo_bussi_local(1, 1, ring%p(idm:idm,iatom:iatom,k), ring%m(iatom:iatom,k), &
          ring%beta_n, tau, ring%thermE, .true.)
      enddo
    enddo
  enddo
  if (.not.nm_orig) call nmtran_backward(ring)
  
  case ("PILE-G")
  nm_orig=ring%nm
  call nmtran_forward(ring)
  if (.not.skip0) then
!Zero-frequency mode:
    tau=ring%tau0/dt !Note dt is whatever this subroutine was called with, which may be half a regular timestep

    call thermo_bussi_global(ring%ndim,ring%natom,ring%p(:,:,1),ring%m(:,1),ring%beta_n,tau,ring%thermE,ring%zcom)
  endif
!Other modes
  do k=2,ring%n
    tau=1.0d0/(dt*2.0d0*ring%omega(k))
    do iatom=1,ring%natom
      do idm=1,ring%ndim
        call thermo_bussi_local(1, 1, ring%p(idm:idm,iatom:iatom,k), ring%m(iatom:iatom,k), &
          ring%beta_n, tau, ring%thermE, .true.)
      enddo
    enddo
  enddo
  if (.not.nm_orig) call nmtran_backward(ring)
end select
end subroutine ring_thermostat

subroutine ring_momupdate(ring,dt,skip0)
implicit none
type (T_Ring) :: ring
real(8) dt
logical :: skip0
integer kstart,kend


if (skip0.and..not.ring%nm) then
  write(6,*)'Error in ring_momupdate. Want to skip centroid, but not in NM'
  stop
endif

kstart=1
if (skip0) kstart=2
kend=ring%n
if (ring%method(1:7)=='GAUSSAV') kend=1

ring%p(:,:,kstart:kend)=ring%p(:,:,kstart:kend) + dt*ring%F_ext(:,:,kstart:kend)

end subroutine ring_momupdate

subroutine ring_freeprop(ring,dt,skip0)
implicit none
type(T_ring) :: ring
real(8) :: dt
logical :: skip0
!real(8) :: wterm,vterm,phi,A,phase
real(8) :: dtm_inv,mwterm,mwterm_inv,wtterm,coswt,sinwt,mwsinwt,mw_invsinwt,pnew,rnew
integer :: iatom,idm,k
logical nm_orig

nm_orig=ring%nm
!Transform to normal modes
call nmtran_forward(ring)

if (.not.skip0) then
!0 frequency mode (centroid) - velocity unchanged (no force). Trivial position update.
  do iatom=1,ring%natom
    dtm_inv=dt*ring%m_inv(iatom,1)
    do idm=1,ring%ndim
      ring%r(idm,iatom,1)=ring%r(idm,iatom,1)+dtm_inv*ring%p(idm,iatom,1)
   enddo
  enddo
endif

!Propagate remaining modes exactly. Eq. 23. Mass taken to be constant in each cartesian direction
do k=2,ring%n
  do iatom=1,ring%natom
    mwterm=ring%m(iatom,k)*ring%omega(k)
    mwterm_inv=1.0d0/mwterm
    wtterm=ring%omega(k)*dt
    coswt=dcos(wtterm)
    sinwt=dsin(wtterm)
    mwsinwt=-mwterm*sinwt
    mw_invsinwt=mwterm_inv*sinwt

    do idm=1,ring%ndim
      pnew=coswt*ring%p(idm,iatom,k)+mwsinwt*ring%r(idm,iatom,k)
      rnew=mw_invsinwt*ring%p(idm,iatom,k)+coswt*ring%r(idm,iatom,k)
      ring%p(idm,iatom,k)=pnew
      ring%r(idm,iatom,k)=rnew
    enddo
  enddo
enddo 

!do ibead=2,ring%n
!  do iatom=1,ring%natom
!    wterm=twopi*ring%omega(ibead)
!    do idm=1,ring%ndim
!      vterm=ring%v(idm,iatom,ibead)/wterm
!      phi=datan(-vterm/ring%r(idm,iatom,ibead)
!      A=dsqrt(ring%r(idm,iatom,ibead)*ring%r(idm,iatom,ibead)+vterm*vterm)
!      phase=wterm*dt+phi
!      ring%r(idm,iatom,ibead)=A*dcos(phase)
!      ring%v(idm,iatom,ibead)=-A*wterm*dsin(phase)
!    enddo
!  enddo
!enddo
!Transform back
if (.not.nm_orig) call nmtran_backward(ring)

return
end subroutine ring_freeprop

subroutine nmtran_forward(ring)
use FFTModule
implicit none
type(T_Ring) :: ring
real(8) :: x(ring%n)

integer iatom,idm

if (ring%nm) then
  write(6,*)'ERROR in nmtran_forward. Already transformed'
  stop
  return
end if

do iatom=1,ring%natom
  do idm=1,ring%ndim
    !Copy ring positions into fft array
    x(:)=ring%r(idm,iatom,:)
    !forward transform
    call fft_r2nm(x,ring%n)
    !copy back into ring structure
    ring%r(idm,iatom,:)=x(:)

    !Copy ring momenta into fft array
    x(:)=ring%p(idm,iatom,:)
    !forward transform
    call fft_r2nm(x,ring%n)
    !copy back into ring structure
    ring%p(idm,iatom,:)=x(:)

    !Copy forces into fft array
    x(:)=ring%F_ext(idm,iatom,:)
    !forward transform
    call fft_r2nm(x,ring%n)
    !copy back into ring structure
    ring%F_ext(idm,iatom,:)=x(:)
  enddo
enddo
ring%nm=.true.

return
end subroutine nmtran_forward

subroutine nmtran_backward(ring)
use FFTModule
type(T_Ring) :: ring
real(8) :: x(ring%n)

integer iatom,idm

if (.not.ring%nm) then
  write(6,*)'ERROR in nmtran_backward. Already transformed'
  stop
  return
end if

do iatom=1,ring%natom
  do idm=1,ring%ndim
    !Copy ring positions into fft array
    x(:)=ring%r(idm,iatom,:)
    !Back transform
    call fft_nm2r(x,ring%n)
    !copy back into ring structure
    ring%r(idm,iatom,:)=x(:)

    !Copy ring momenta into fft array
    x(:)=ring%p(idm,iatom,:)
    !Back transform
    call fft_nm2r(x,ring%n)
    !copy back into ring structure
    ring%p(idm,iatom,:)=x(:)

    !Copy forces into fft array
    x(:)=ring%F_ext(idm,iatom,:)
    !Back transform
    call fft_nm2r(x,ring%n)
    !copy back into ring structure
    ring%F_ext(idm,iatom,:)=x(:)
  enddo
enddo

ring%nm=.false.

return
end subroutine nmtran_backward

end module RingModule
