! Copyright Todd J. Martinez and Raphael D. Levine, 1994
Module ThermoModule
use GlobalModule
implicit none
private
public :: thermo_init, thermo, thermo_bussi_global, thermo_bussi_local, thermo_NormDist, thermo_MBDist
character(len=8),save :: therm="" !which thermostat to use in simple interface
double precision, save :: beta_s=0.0 !inverse temperature in simple interface
double precision,save :: thermtime=0.0 !Thermostat relaxation time (in au) for simple interface
logical, save :: zcom_s=.true. !Do center of mass removal in simple interface?
integer(4),parameter :: i4one=1 !Ensure type is consistent with uses of integer(4)

contains

subroutine thermo_init(iseed,therm_in,beta_in,thermtime_in,zcom_in)
integer(4) :: iseed
character(len=8),optional :: therm_in
double precision,optional :: beta_in, thermtime_in
logical,optional :: zcom_in

double precision :: foo
integer(4) :: idum

idum=sign(iseed,-i4one) !ensure seed is negative
call thermo_ran(foo,idum) !Initialize random number generator

if (present(therm_in)) therm=therm_in !Set the thermostat for the simple interface
if (present(beta_in)) beta_s=beta_in !Set the temperature for the simple interface
if (present(zcom_in)) zcom_s=zcom_in !Set the center of mass removal flag for the simple interface
if (present(thermtime_in)) thermtime=thermtime_in !Set the thermostat relaxation time for the simple interface

end subroutine thermo_init

subroutine thermo(ndim,natom,p,mass,dt,thermE)
integer :: ndim,natom
double precision :: p(ndim,natom),mass(natom),dt,tau,thermE

tau=thermtime/dt

select case (therm)
  case ("bussi_g")
    call thermo_bussi_global(ndim,natom,p,mass,beta_s,tau,thermE,zcom_s)
  case ("bussi_l")
    call thermo_bussi_local(ndim,natom,p,mass,beta_s,tau,thermE,zcom_s)
  case ("mbdist")
    call thermo_MBDist(ndim,natom,p,mass,beta_s,zcom_s)
  case default !do nothing, and return

end select
end subroutine thermo

function thermo_NormDist(sigma) result(x)
!This subroutine generates a normal distribution exp(-x^2/2sigma^2)
double precision :: x
double precision,intent(in) :: sigma

x=gasdev()*sigma
end function thermo_NormDist

subroutine thermo_MBDist(ndim,natom,p,mass,beta,zcom)
integer :: ndim,natom
double precision :: p(ndim,natom),mass(natom)
double precision :: beta
logical :: zcom

double precision :: mxw_unitf,sigma
integer :: iatom,idm

!Sample momenta from Maxwell-Boltzmann distribution
mxw_unitf=1.0/sqrt(beta)
do iatom=1,natom
  sigma=mxw_unitf*sqrt(mass(iatom))
  do idm=1,ndim
    p(idm,iatom)=thermo_NormDist(sigma)
  enddo
enddo


!Project out center of mass velocity if not a degree of freedom
if (.not.zcom) then
  call vcom_project(ndim,natom,p,mass)
endif
end subroutine thermo_MBDist

subroutine thermo_bussi_global(ndim,natom,p,mass,beta,tau,thermE,zcom)
integer :: ndim,natom
double precision :: p(ndim,natom),mass(natom)
double precision :: beta,tau,thermE
logical :: zcom

double precision :: ekin_old,ekin_new,ekinfac,sigma,signfac,vscale
integer :: ndeg,iatom,idm

if (.not.zcom) then
!Do not include center of mass as degree of freedom
  call vcom_project(ndim,natom,p,mass)
  ndeg=ndim*natom-ndim !Rotations are not frozen out, so keep them as degrees of freedom
else
  ndeg=ndim*natom
endif

ekin_old=0.0
do iatom=1,natom
  ekinfac=0.5/mass(iatom)
  do idm=1,ndim
    ekin_old=ekin_old+ekinfac*p(idm,iatom)*p(idm,iatom)
  enddo
enddo

sigma=0.5*real(ndeg)/beta

ekin_new=resamplekin(ekin_old,sigma,ndeg,tau)

signfac=sign(1.0d0,ekin_new)
ekin_new=abs(ekin_new)

thermE=thermE+ekin_old-ekin_new

vscale=signfac*sqrt(ekin_new/ekin_old)

do iatom=1,natom
  do idm=1,ndim
    p(idm,iatom)=p(idm,iatom)*vscale
  enddo
enddo

end subroutine thermo_bussi_global

subroutine thermo_bussi_local(ndim,natom,p,mass,beta,tau,thermE,zcom)
integer :: ndim,natom
double precision :: p(ndim,natom),mass(natom)
double precision :: beta,tau,thermE
logical :: zcom

double precision :: vtot(ndim)
double precision :: c1,c2, c2fac, ekin_old, ekin_new, ekinfac
integer :: iatom,idm

c1=exp(-1.0/tau)
c2=sqrt(1.0-c1*c1)

ekin_old=0.0
do iatom=1,natom
  c2fac=c2*sqrt(mass(iatom)/beta)
  ekinfac=0.5/mass(iatom)
  do idm=1,ndim
    ekin_old=ekin_old+ekinfac*p(idm,iatom)*p(idm,iatom)
    p(idm,iatom)=c1*p(idm,iatom) + c2fac*gasdev()
  enddo
enddo

if (.not.zcom) then
!Do not include center of mass as degree of freedom
  vtot(1:ndim)=0.0
  do iatom=1,natom
    do idm=1,ndim
      vtot(idm)=vtot(idm)+p(idm,iatom)/mass(iatom)
    enddo
  enddo
  vtot(1:ndim)=vtot(1:ndim)/real(natom)
  do iatom=1,natom
    do idm=1,ndim
      p(idm,iatom)=p(idm,iatom)-mass(iatom)*vtot(idm)
    enddo
  enddo
endif

ekin_new=0.0
do iatom=1,natom
  ekinfac=0.5/mass(iatom)
  do idm=1,ndim
    ekin_new=ekin_new+ekinfac*p(idm,iatom)*p(idm,iatom)
  enddo
enddo
!Calculate new thermostat energy
thermE=thermE+ekin_old-ekin_new

end subroutine thermo_bussi_local

function resamplekin(kk,sigma,ndeg,taut)
  double precision               :: resamplekin
  double precision,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  double precision,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  double precision,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  double precision :: factor,rr
  double precision ::dndeg
  dndeg=real(ndeg)
  if(taut>0.1) then
    factor=exp(-1.0/taut)
  else
    factor=0.0
  end if
  rr = gasdev()
  resamplekin = kk + (1.0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/dndeg-kk) &
               + 2.0*rr*sqrt(kk*sigma/dndeg*(1.0-factor)*factor)

! Transfer sign to resamplekin
  resamplekin=sign(resamplekin,rr+sqrt(dndeg*factor/(sigma*(1.0-factor))))

end function resamplekin

double precision function sumnoises(nn)
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0*gamdev(nn/2)
  else
    sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

! gamma-distributed random number
! TODO: Reimplement this
function gamdev(ia)
double precision :: gamdev
integer, intent(in) :: ia

gamdev = 0.0D0
call FMS_DieError("ERROR: gamdev not implemented")
end function gamdev

! TODO: Reimplement this
function gasdev()
double precision :: gasdev

gasdev = 0.0D0
call FMS_DieError("ERROR: gasdev not implemented")
end function gasdev

SUBROUTINE thermo_ran(rnd,iseed)
! interface to random number generators
integer(4), optional :: iseed
double precision :: rnd

call FMS_DieError("ERROR: thermo_ran not implemented")
if (present(iseed)) then
  ! TODO: Use FMS_ranb
  !rnd=ran1(iseed)
endif
!rnd=ran1()
end subroutine thermo_ran

subroutine vcom_project(ndim,natom,p,mass)
integer :: ndim,natom
double precision :: p(ndim,natom)
double precision :: mass(natom)

double precision :: vtot(ndim)
integer :: iatom,idm

vtot(1:ndim)=0.0
do iatom=1,natom
  do idm=1,ndim
    vtot(idm)=vtot(idm)+p(idm,iatom)/mass(iatom)
  enddo
enddo
vtot(1:ndim)=vtot(1:ndim)/real(natom)
do iatom=1,natom
  do idm=1,ndim
    p(idm,iatom)=p(idm,iatom)-mass(iatom)*vtot(idm)
  enddo
enddo

end subroutine vcom_project

end module ThermoModule
