! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module ThermoModule
   use GlobalModule
   implicit none
   private
   public :: thermo_init, thermo, thermo_bussi_global, thermo_bussi_local, thermo_NormDist, &
           thermo_MBDist, thermo_NHC_MBDist, thermo_NHC_local, thermo_NHC_global
   character(len=8), save :: therm = "" !which thermostat to use in simple interface
   double precision, save :: beta_s = 0.0 !inverse temperature in simple interface
   double precision, save :: thermtime = 0.0 !Thermostat relaxation time (in au) for simple interface
   logical, save :: zcom_s = .true. !Do center of mass removal in simple interface?
   integer(4), parameter :: i4one = 1 !Ensure type is consistent with uses of integer(4)

contains

   subroutine thermo_init(iseed, therm_in, beta_in, thermtime_in, zcom_in)
      integer(4) :: iseed
      character(len=8), optional :: therm_in
      double precision, optional :: beta_in, thermtime_in
      logical, optional :: zcom_in

      double precision :: foo
      integer(4) :: idum

      idum = sign(iseed, -i4one) !ensure seed is negative
      call thermo_ran(foo, idum) !Initialize random number generator

      if (present(therm_in)) therm = therm_in !Set the thermostat for the simple interface
      if (present(beta_in)) beta_s = beta_in !Set the temperature for the simple interface
      if (present(zcom_in)) zcom_s = zcom_in !Set the center of mass removal flag for the simple interface
      if (present(thermtime_in)) thermtime = thermtime_in !Set the thermostat relaxation time for the simple interface

   end subroutine thermo_init

   subroutine thermo(ndim, natom, p, mass, dt, thermE)
      integer :: ndim, natom
      double precision :: p(ndim, natom), mass(natom), dt, tau, thermE

      tau = thermtime / dt

      select case (therm)
      case ("bussi_g")
         call thermo_bussi_global(ndim, natom, p, mass, beta_s, tau, thermE, zcom_s)
      case ("bussi_l")
         call thermo_bussi_local(ndim, natom, p, mass, beta_s, tau, thermE, zcom_s)
      case ("mbdist")
         call thermo_MBDist(ndim, natom, p, mass, beta_s, zcom_s)
      case default !do nothing, and return

      end select
   end subroutine thermo

   function thermo_NormDist(sigma) result(x)
!This subroutine generates a normal distribution exp(-x^2/2sigma^2)
      double precision :: x
      double precision, intent(in) :: sigma

      x = gasdev() * sigma
   end function thermo_NormDist

   subroutine thermo_MBDist(ndim, natom, p, mass, beta, zcom)
      integer :: ndim, natom
      double precision :: p(ndim, natom), mass(natom)
      double precision :: beta
      logical :: zcom

      double precision :: mxw_unitf, sigma
      integer :: iatom, idm

!Sample momenta from Maxwell-Boltzmann distribution
      mxw_unitf = 1.0 / sqrt(beta)
      do iatom = 1, natom
         sigma = mxw_unitf * sqrt(mass(iatom))
         do idm = 1, ndim
            p(idm, iatom) = thermo_NormDist(sigma)
         end do
      end do

!Project out center of mass velocity if not a degree of freedom
      if (.not. zcom) then
         call vcom_project(ndim, natom, p, mass)
      end if
   end subroutine thermo_MBDist

   subroutine thermo_bussi_global(ndim, natom, p, mass, beta, tau, thermE, zcom)
      integer :: ndim, natom
      double precision :: p(ndim, natom), mass(natom)
      double precision :: beta, tau, thermE
      logical :: zcom

      double precision :: ekin_old, ekin_new, ekinfac, sigma, signfac, vscale
      integer :: ndeg, iatom, idm

      if (.not. zcom) then
!Do not include center of mass as degree of freedom
         call vcom_project(ndim, natom, p, mass)
         ndeg = ndim * natom - ndim !Rotations are not frozen out, so keep them as degrees of freedom
      else
         ndeg = ndim * natom
      end if

      ekin_old = 0.0
      do iatom = 1, natom
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_old = ekin_old + ekinfac * p(idm, iatom) * p(idm, iatom)
         end do
      end do

      sigma = 0.5 * real(ndeg) / beta

      ekin_new = resamplekin(ekin_old, sigma, ndeg, tau)

      signfac = sign(1.0d0, ekin_new)
      ekin_new = abs(ekin_new)

      thermE = thermE + ekin_old - ekin_new

      vscale = signfac * sqrt(ekin_new / ekin_old)

      do iatom = 1, natom
         do idm = 1, ndim
            p(idm, iatom) = p(idm, iatom) * vscale
         end do
      end do

   end subroutine thermo_bussi_global

   subroutine thermo_bussi_local(ndim, natom, p, mass, beta, tau, thermE, zcom)
      integer :: ndim, natom
      double precision :: p(ndim, natom), mass(natom)
      double precision :: beta, tau, thermE
      logical :: zcom

      double precision :: vtot(ndim)
      double precision :: c1, c2, c2fac, ekin_old, ekin_new, ekinfac
      integer :: iatom, idm

      c1 = exp(-1.0 / tau)
      c2 = sqrt(1.0 - c1 * c1)

      ekin_old = 0.0
      do iatom = 1, natom
         c2fac = c2 * sqrt(mass(iatom) / beta)
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_old = ekin_old + ekinfac * p(idm, iatom) * p(idm, iatom)
            p(idm, iatom) = c1 * p(idm, iatom) + c2fac * gasdev()
         end do
      end do

      if (.not. zcom) then
!Do not include center of mass as degree of freedom
         vtot(1:ndim) = 0.0
         do iatom = 1, natom
            do idm = 1, ndim
               vtot(idm) = vtot(idm) + p(idm, iatom) / mass(iatom)
            end do
         end do
         vtot(1:ndim) = vtot(1:ndim) / real(natom)
         do iatom = 1, natom
            do idm = 1, ndim
               p(idm, iatom) = p(idm, iatom) - mass(iatom) * vtot(idm)
            end do
         end do
      end if

      ekin_new = 0.0
      do iatom = 1, natom
         ekinfac = 0.5 / mass(iatom)
         do idm = 1, ndim
            ekin_new = ekin_new + ekinfac * p(idm, iatom) * p(idm, iatom)
         end do
      end do
!Calculate new thermostat energy
      thermE = thermE + ekin_old - ekin_new

   end subroutine thermo_bussi_local

   function resamplekin(kk, sigma, ndeg, taut)
      double precision :: resamplekin
      double precision, intent(in) :: kk ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
      double precision, intent(in) :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
      integer, intent(in) :: ndeg ! number of degrees of freedom of the atoms to be thermalized
      double precision, intent(in) :: taut ! relaxation time of the thermostat, in units of 'how often this routine is called'
      double precision :: factor, rr
      double precision :: dndeg
      dndeg = real(ndeg)
      if (taut > 0.1) then
         factor = exp(-1.0 / taut)
      else
         factor = 0.0
      end if
      rr = gasdev()
      resamplekin = kk + (1.0 - factor) * (sigma * (sumnoises(ndeg - 1) + rr**2) / dndeg - kk) &
                    + 2.0 * rr * sqrt(kk * sigma / dndeg * (1.0 - factor) * factor)

! Transfer sign to resamplekin
      resamplekin = sign(resamplekin, rr + sqrt(dndeg * factor / (sigma * (1.0 - factor))))

   end function resamplekin

   double precision function sumnoises(nn)
      integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
      if (nn == 0) then
         sumnoises = 0.0
      else if (nn == 1) then
         sumnoises = gasdev()**2
      else if (modulo(nn, 2) == 0) then
         sumnoises = 2.0 * gamdev(nn / 2)
      else
         sumnoises = 2.0 * gamdev((nn - 1) / 2) + gasdev()**2
      end if
   end function sumnoises

! gamma-distributed random number
! TODO: Reimplement this
   function gamdev(ia)
      double precision :: gamdev
      integer, intent(in) :: ia

      gamdev = 0.0d0
      call FMS_DieError("ERROR: gamdev not implemented")
   end function gamdev

! TODO: Reimplement this
   function gasdev()
      double precision :: gasdev

      gasdev = 0.0d0
      call FMS_DieError("ERROR: gasdev not implemented")
   end function gasdev

   subroutine thermo_ran(rnd, iseed)
! interface to random number generators
      integer(4), optional :: iseed
      double precision :: rnd

      call FMS_DieError("ERROR: thermo_ran not implemented")
      if (present(iseed)) then
         ! TODO: Use FMS_ranb
         !rnd=ran1(iseed)
      end if
!rnd=ran1()
   end subroutine thermo_ran

   subroutine vcom_project(ndim, natom, p, mass)
      integer :: ndim, natom
      double precision :: p(ndim, natom)
      double precision :: mass(natom)

      double precision :: vtot(ndim)
      integer :: iatom, idm

      vtot(1:ndim) = 0.0
      do iatom = 1, natom
         do idm = 1, ndim
            vtot(idm) = vtot(idm) + p(idm, iatom) / mass(iatom)
         end do
      end do
      vtot(1:ndim) = vtot(1:ndim) / real(natom)
      do iatom = 1, natom
         do idm = 1, ndim
            p(idm, iatom) = p(idm, iatom) - mass(iatom) * vtot(idm)
         end do
      end do

   end subroutine vcom_project

   subroutine thermo_NHC_MBDist(beta,M,pNHC,mNHC,thermE)
   !Initialise Nose cariables
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
   
   !Eqs.35 in Martyna et al. 1996 for the NHC propagator
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










end module ThermoModule
