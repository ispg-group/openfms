Module ThermoModule
implicit none
private
public :: thermo_init, thermo_NHC_MBDist, thermo_NHC_local, thermo_NHC_global, &
        thermo_bussi_global, thermo_bussi_local, thermo_NormDist, thermo_MBDist, vcom_project

contains

subroutine thermo_init(iseed)
implicit none
integer iseed

real(8) foo
integer idum

idum=sign(iseed,-1) !ensure seed is negative
call thermo_ran(foo,idum) !Initialize random number generator

end subroutine thermo_init

function thermo_NormDist(sigma) result(x)
!This subroutine generates a normal distribution exp(-x^2/2sigma^2)
implicit none
real(8) :: x
real(8),intent(in) :: sigma

x=gasdev()*sigma
end function thermo_NormDist

subroutine thermo_MBDist(ndim,natom,p,mass,beta,zcom)
implicit none
integer ndim,natom
real(8) p(ndim,natom),mass(natom)
real(8) beta
logical zcom

real(8) mxw_unitf,sigma
integer iatom,idm

!Sample momenta from Maxwell-Boltzmann distribution
mxw_unitf=1.0d0/dsqrt(beta)
do iatom=1,natom
  sigma=mxw_unitf*dsqrt(mass(iatom))
  do idm=1,ndim
    p(idm,iatom)=thermo_NormDist(sigma)
  enddo
enddo


!Project out center of mass velocity if not a degree of freedom
if (.not.zcom) then
  call vcom_project(ndim,natom,p,mass)
endif
end subroutine thermo_MBDist

subroutine thermo_NHC_MBDist(beta,M,pNHC,mNHC,thermE)
!Initialize Nose variables
implicit none
integer M
real(8) beta,thermE
real(8),dimension(M) :: mNHC,pNHC
real(8) sigma
integer i

thermE=0.0d0
do i=1,M
  sigma=dsqrt(mNHC(i)/beta) !Maxwell-Boltzmann distribution of Nose momenta
  pNHC(i)=thermo_NormDist(sigma)
  thermE=thermE+0.5d0*pNHC(i)*pNHC(i)/mNHC(i)
enddo
end subroutine thermo_NHC_MBDist

subroutine thermo_NHC_local(dt,M,p,mass,beta,rNHC,pNHC,mNHC,thermE)
implicit none
integer M
real(8) dt,p,mass,beta,thermE
real(8),dimension(M) :: mNHC,rNHC,pNHC

real(8) akin,scale
integer l

akin=p*p/mass !Note factor of 2. Grr Glenn....
call thermo_NHC_prop(dt,akin,beta,1,M,rNHC,pNHC,mNHC,scale)

!scale momenta
p=p*scale

!Calculate thermostat energy
thermE=0.0d0
do l=1,M
  thermE=thermE+0.5d0*pNHC(l)*pNHC(l)/mNHC(l) + rNHC(l)/beta
enddo

end subroutine thermo_NHC_local

subroutine thermo_NHC_global(dt,M,ndim,natom,p,mass,beta,rNHC,pNHC,mNHC,thermE)
implicit none
integer M,ndim,natom
real(8),dimension(ndim,natom) :: p
real(8),dimension(natom) :: mass
real(8) dt,beta,thermE
real(8),dimension(M) :: mNHC,rNHC,pNHC

real(8) akin,scale
integer iatom,idm,l,NTotDim

if(any(mass.lt.1.0d-15)) return !Do not thermostat if mass is zero
akin=0.0d0
do iatom=1,natom
  do idm=1,ndim
    akin=akin+p(idm,iatom)*p(idm,iatom)/mass(iatom)
  enddo
enddo
NTotDim=ndim*natom
call thermo_NHC_prop(dt,akin,beta,NTotDim,M,rNHC,pNHC,mNHC,scale)
!scale momenta
p=p*scale
!Calculate thermostat energy
thermE=0.0d0
!1st chain
thermE=thermE+0.5d0*pNHC(1)*pNHC(1)/mNHC(1) + dble(NTotDim)*rNHC(1)/beta
do l=2,M
  thermE=thermE+0.5d0*pNHC(l)*pNHC(l)/mNHC(l) + rNHC(l)/beta
enddo
end subroutine thermo_NHC_global

subroutine thermo_NHC_prop(dt,akin,beta,N,M,rNHC,pNHC,mNHC,scale)
!Propagates NHC variables from t to dt and scales momenta
implicit none
integer N     !Number of degrees of freedom
integer M     !Length of chain
real(8) dt    !Time step for this update (usually half of system timestep)
real(8) akin  !twice system kinetic energy to which chain 1 is coupled
real(8) beta  !Inverse temperature
real(8) scale !Scaling factor for system momenta
real(8),dimension(M) :: mNHC & ! NHC masses
                       ,rNHC &! NHC positions
                       ,pNHC ! NHC momenta

!hard coded parameters. These have been tested and appear to be optimal, atleast
!for 100 bead harmonic oscillator with 4 NHC variables per mode
integer,parameter :: nc=10 ! Number of Trotter Factorizations
integer,parameter :: nysmax=7
integer,parameter :: nys=5 !Order of Yoshida-Suzuki Integration

!local variables
real(8),dimension(M) :: G,mNHC_inv
real(8),dimension(nysmax) :: w
real(8) :: beta_inv,dt_nc,dts,AA,dNbeta_inv
integer j,k,l

scale=1.0d0
!If any NHC mass is zero, we disable NHC propagation
if (any(mNHC.lt.1.0d-15)) then
  write(6,*)'mNHC is zero'
  return
endif

!Set up higher-order integration weights
select case(nys)
case (1)
  w(1)=1.0d0
case (3)
  w(1)=1.0d0/(2.0d0-2.0d0**(1.0d0/3.0d0))
  w(2)=1.0d0-2.0d0*w(1)
  w(3)=w(1)
case (5)
  w(1)=1.0d0/(4.0d0-4.0d0**(1.0d0/3.0d0))
  w(2)=w(1)
  w(3)=1.0d0-4.0d0*w(1)
  w(4)=w(1)
  w(5)=w(1)
case (7)
  w(1)=-0.117767998417887D1
  w(2)=0.235573213359357D0
  w(3)=0.784513610477560D0
  w(4)=1.0D0 - 2.0D0*(w(1)+w(2)+w(3))
  w(5)=w(3)
  w(6)=w(2)
  w(7)=w(1)
case default
  write(*,*) 'Error, nys=',nys,' not coded'
end select

!Some useful variables
beta_inv=1.0d0/beta
dt_nc=dt/dble(nc)
mNHC_inv(1:M)=1.0d0/mNHC(1:M)
dNbeta_inv=beta_inv*dble(N)

!Update NHC forces
G(1)=akin-dNbeta_inv
G(2:M)=mNHC_inv(1:M-1)*pNHC(1:M-1)*pNHC(1:M-1)-beta_inv

do k=1,nc
  do j=1,nys
    dts=w(j)*dt_nc
    !Update NHC momenta
    pNHC(M)=pNHC(M)+0.5d0*dts*G(M)
    do l=M-1,1,-1
      AA=dexp(-0.25d0*dts*pNHC(l+1)*mNHC_inv(l+1))
      pNHC(l)=pNHC(l)*AA*AA+0.5d0*dts*G(l)*AA
    enddo
    !Update particle momenta
    AA=dexp(-dts*pNHC(1)*mNHC_inv(1))
    scale=scale*AA
    !Update NHC forces
    G(1)=scale*scale*akin-dNbeta_inv
    !Update NHC positions
    do l=1,M
      rNHC(l)=rNHC(l)+pNHC(l)*dts*mNHC_inv(l)
    enddo
    !Update NHC momenta
    do l=1,M-1
      AA=dexp(-0.25d0*dts*pNHC(l+1)*mNHC_inv(l+1))
      pNHC(l)=pNHC(l)*AA*AA+0.5d0*dts*G(l)*AA
      G(l+1)=mNHC_inv(l)*pNHC(l)*pNHC(l)-beta_inv
    enddo
    pNHC(M)=pNHC(M)+G(M)*0.5d0*dts
  enddo
enddo
return

end subroutine thermo_NHC_prop


subroutine thermo_bussi_global(ndim,natom,p,mass,beta,tau,thermE,zcom)
implicit none
integer ndim,natom
real(8) p(ndim,natom),mass(natom)
real(8) beta,tau,thermE
logical zcom

real(8) ekin_old,ekin_new,ekinfac,sigma,signfac,vscale
integer ndeg,iatom,idm

ekin_old=0.0d0
do iatom=1,natom
  ekinfac=0.5d0/mass(iatom)
  do idm=1,ndim
    ekin_old=ekin_old+ekinfac*p(idm,iatom)*p(idm,iatom)
  enddo
enddo

if (.not.zcom) then
!Do not include center of mass as degree of freedom
  call vcom_project(ndim,natom,p,mass)
  ndeg=ndim*natom-ndim !Rotations are not frozen out, so keep them as degrees of freedom
else
  ndeg=ndim*natom
endif

if (tau.gt.1.0d-15) then !only apply thermostat with nonzero relaxation time
  sigma=0.5d0*dble(ndeg)/beta

  ekin_new=resamplekin(ekin_old,sigma,ndeg,tau)

  signfac=sign(1.0d0,ekin_new)
  ekin_new=dabs(ekin_new)

  vscale=signfac*dsqrt(ekin_new/ekin_old)

  do iatom=1,natom
    do idm=1,ndim
      p(idm,iatom)=p(idm,iatom)*vscale
    enddo
  enddo
end if

ekin_new=0.0d0
do iatom=1,natom
  ekinfac=0.5d0/mass(iatom)
  do idm=1,ndim
    ekin_new=ekin_new+ekinfac*p(idm,iatom)*p(idm,iatom)
  enddo
enddo
!Calculate new thermostat energy
thermE=thermE+ekin_old-ekin_new

end subroutine thermo_bussi_global

subroutine thermo_bussi_local(p,mass,beta,tau,thermE)
implicit none
real(8) p,mass
real(8) beta,tau,thermE

real(8) c1,c2, c2fac, ekin_old, ekin_new, ekinfac

if (tau.lt.1.0d-15) return !only apply thermostat with nonzero relaxation time

ekinfac=0.5d0/mass
ekin_old=ekinfac*p*p

c1=dexp(-1.0d0/tau)
c2=dsqrt(1.0d0-c1*c1)

c2fac=c2*dsqrt(mass/beta)
p=c1*p+ c2fac*gasdev()

ekin_new=ekinfac*p*p

!Calculate new thermostat energy
thermE=thermE+ekin_old-ekin_new

end subroutine thermo_bussi_local

function resamplekin(kk,sigma,ndeg,taut)
  implicit none
  real*8               :: resamplekin
  real*8,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  real*8,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  real*8,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  real*8 :: factor,rr
  real*8 ::dndeg
  dndeg=dble(ndeg)
  if(taut>0.1d0) then
    factor=dexp(-1.0d0/taut)
  else
    factor=0.0d0
  end if
  rr = gasdev()
  resamplekin = kk + (1.0d0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/dndeg-kk) &
               + 2.0d0*rr*dsqrt(kk*sigma/dndeg*(1.0d0-factor)*factor)

! Transfer sign to resamplekin
  resamplekin=sign(resamplekin,rr+dsqrt(dndeg*factor/(sigma*(1.0d0-factor))))

end function resamplekin

double precision function sumnoises(nn)
  implicit none
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
  if(nn==0) then
    sumnoises=0.0d0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0d0*gamdev(nn/2)
  else
    sumnoises=2.0d0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

function gamdev(ia)
! gamma-distributed random number, implemented as described in numerical recipes

implicit none
real(8) :: gamdev
integer, intent(in) :: ia
integer j
real(8) am,e,s,v1,v2,x,y,z
if(ia.lt.1)pause 'bad argument in gamdev'
if(ia.lt.6)then
  x=1.
  do 11 j=1,ia
    call thermo_ran(z)
    x=x*z
11  continue
  x=-log(x)
else
1 call thermo_ran(z)
    v1=2.0d0*z-1.0d0
    call thermo_ran(z)
    v2=2.0d0*z-1.0d0
  if(v1**2+v2**2.gt.1.0d0)goto 1
    y=v2/v1
    am=ia-1
    s=dsqrt(2.0d0*am+1.0d0)
    x=s*y+am
  if(x.le.0.0d0)goto 1
    e=(1.0d0+y**2)*dexp(am*dlog(x/am)-s*y)
  call thermo_ran(z)
  if(z.gt.e)goto 1
endif
gamdev=x
return
end function gamdev

function gasdev()

implicit none
real(8) gasdev
real(8) fac,gset,rsq,v1,v2,z
INTEGER iset
SAVE iset,gset
DATA iset/0/

if (iset.eq.0) then
1 call thermo_ran(z) 
 v1=2.0d0*z-1.0d0
    call thermo_ran(z)
    v2=2.0d0*z-1.0d0
    rsq=v1**2+v2**2
    if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)goto 1
    fac=sqrt(-2.0d0*log(rsq)/rsq)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
else
    gasdev=gset
    iset=0
endif
return
END FUNCTION GASDEV

SUBROUTINE THERMO_RAN(rnd,iseed)
! interface to random number generators
implicit none
integer, optional :: iseed
real(8) rnd, foo

if (present(iseed)) then
!  call init_random_seed(iseed)

  foo=ran1(iseed)
endif

rnd=ran1()
!  call RANDOM_NUMBER(rnd)
return
end SUBROUTINE THERMO_RAN

SUBROUTINE init_random_seed(iseed)
IMPLICIT NONE
INTEGER :: i, n, iseed
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
          
seed = iseed * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
          
DEALLOCATE(seed)
END SUBROUTINE init_random_seed


FUNCTION ran1(iseed)
! random number generator
integer, optional :: iseed
INTEGER IA,IM,IQ,IR,NTAB,NDIV
REAL ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
INTEGER, SAVE :: idum=0 
if (present(iseed)) idum=iseed
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
  do 11 j=NTAB+8,1,-1
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11  continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
end function ran1

subroutine vcom_project(ndim,natom,p,mass)
implicit none
integer ndim,natom
real(8) p(ndim,natom)
real(8) mass(natom)

real(8) vtot(ndim)
real(8) Mtot_inv
integer iatom,idm

vtot(1:ndim)=0.0d0
Mtot_inv=1.0d0/sum(mass)

do iatom=1,natom
  do idm=1,ndim
    vtot(idm)=vtot(idm)+p(idm,iatom)*Mtot_inv
  enddo
enddo
do iatom=1,natom
  do idm=1,ndim
    p(idm,iatom)=p(idm,iatom)-mass(iatom)*vtot(idm)
  enddo
enddo

end subroutine vcom_project

end module ThermoModule
