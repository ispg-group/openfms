!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Ground state sampling routines, including sampling Wigner, Husimi and Boltzmann distributions.
!<
module SamplingModule
   use GlobalModule
   use TrajectoryModule
   use QM_MM_Module, only: qcNumQM
   use RandomModule, only: fms_ranb
   implicit none
   private

   public :: FMS_InitialWigner, FMS_InitialHusimi, FMS_InitialConstT
   public :: FMS_InitialSwarm, FMS_InitialQuasi, FMS_InitialColdWig

!     Integer parameter definitions:
   integer(kind=DefInt), public, parameter :: &
      NOSAMPLE = 0, &
      WIGNER = 1, &
      QUASICLASS = 2, &
      BOLTZ = 3, &
      COLDWIG = 4, &
      SWARM = 5, &
      HUSIMI = 6

!>    System temperature
   real(kind=DefReal), public :: indTemperature

!>    Determines ground state sampling method
   integer(kind=DefInt), public :: inInitialCond

   logical, public :: inzFirstGauss

!> For ConstT initial conditions, number of special modes
   integer(kind=DefInt), public :: inNTemp
!> For ConstT initial conditions, index of each special mode
   integer(kind=DefInt), public :: inIModeSharp(100)
!> For ConstT initial conditions, temperature for each of the special modes
   real(kind=DefReal), public :: indFModeSharp(100)

!> For quasiclassical initial conditions, number of modes with extra quanta
   integer(kind=DefInt), public :: inNAddQuanta
!> For quasiclassical initial conditions, modes with extra quanta
   integer(kind=DefInt), public :: inIAddQuanta(100)

   save
contains

!>
!!    @brief This routine samples initial conditions from a Wigner Distribution.
!!    For each Normal mode we sample from:
!!    \f[
!!    \rho(x,p)=\exp\left[ -2\alpha_x(x-x0)^2 -\frac{(p-p_0)^2}{2*\alpha_x} \right]
!!    \f]
!!    where x, and p are the normal mode position and momentum, respectively.
!!    The width \f$\alpha_x\f$ is Mass*Freq/2 [Heller JCP vol 80 p 5038 (1984)]
!!    But Normal modes have units of sqrt[mass]*Length and therefore
!!    \f$\alpha_x=Freq/2\f$.
!!    If indTemperature /=0 then the width is defined as
!!    AlphaX=Freq/2*Tanh(Beta*Freq/2) where Beta=1/(KbT)
!!    See Heller JCP, vol, 65, pp. 1289 (1976)
!!    This is the result for a Wigner at finite temperature
!!    @see FMS_InitialColdWig
!!    @ingroup initial
!<
   subroutine FMS_InitialWigner(T1, IDNum)
      type(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: IDNum

      real(DefReal), dimension(qcNumQM) :: Mass
      real(DefReal), dimension(3*qcNumQM) :: freq

      real(DefReal), dimension(3*qcNumQM, 3*qcNumQM) :: Umat

      ! this are all back to front!
      real(DefReal), dimension(qcNumQM, 3) :: atomcen, P, R, Vel

      real(DefReal) :: Beta
      real(DefReal) :: rsq, dx1, dx2, fac, x1, x2, alphax, sigp, sigx, ft

      real(kind=DefReal) :: PosVec(3 * qcNumQM)
      real(kind=DefReal) :: Momvec(3 * qcNumQM)

      real(kind=DefReal) :: TotMass, rhototal, RCM(3), VelCM(3)
      integer(kind=DefInt) :: I, J, K, L, natoms, nfreq, IDim, IParticle
      real(kind=DefReal) :: ZPE, Ekin

!     Beta is 1/KbT
      beta = 0.0d0
      if (indTemperature > d0) then
         write (fmiOut, *) 'Sampling at finite T =', indTemperature
         Beta = d1 / (indTemperature * BoltzK)
      end if

!     Read initial positions
      do i = 1, qcNumQM
         atomcen(i, 1:3) = T1%Particle(i)%get_pos()
      end do

!     now we can allocate memory
      natoms = qcNumQM

      do i = 1, natoms
         Mass(i) = T1%Particle(i)%Mass
      end do
      call read_frequencies(natoms, Mass, nfreq, freq, Umat)

!     Start sampling
      Ekin = 0.0d0
!     If this is the first trajectory and inzFirstGauss is true place the
!     basis function at the origin.  It will be shifted to atomcen/atvel below.
      if (inzFirstGauss .and. IDNum == 1) then
         PosVec = 0.0d0
         Momvec = 0.0d0
      else
!        Sample from a Gaussian distribution
         rhototal = 1.0d0
         do I = 1, nfreq
            Alphax = freq(i) / 2.0d0
            if (indTemperature <= d0) then
               ft = 1.0d0
            else
               ft = tanh(Beta * freq(i) * dp5)
            end if
            Sigx = sqrt(1.0d0 / (4.0d0 * Alphax * ft))
            Sigp = sqrt(1.0d0 * Alphax / ft)

1           continue
            dx1 = 2.d0 * FMS_ranb(i4zero) - 1.d0
            dx2 = 2.d0 * FMS_ranb(i4zero) - 1.d0
            rsq = dx1 * dx1 + dx2 * dx2
            if (rsq >= 1.d0 .or. rsq == 0.d0) goto 1
            fac = sqrt(-2.d0 * log(rsq) / rsq)
            x1 = dx1 * fac
            x2 = dx2 * fac
! Each Normal mode is in Equilibrium in both position and
! mometum, i.e., x0 and p0 are zero so we just need to "stretch"
! the random numbers using the proper Sigma.
! Note that p is the momentum but the normal modes have a unit mass
! so we do not need to divide by the mass when we go to velocity

            PosVec(I) = Sigx * x1
            Momvec(I) = Sigp * x2
! For debugging purpose calculate Kinetic energy
            Ekin = Ekin + Momvec(i) * Momvec(I)
         end do
      end if
! Kinetic energy in Normal-Modes
      Ekin = Ekin * dp5
!     multiply by U (transform to mass-weighted Cartesian--position/momentum)
      L = 0
      do I = 1, natoms
         do J = 1, 3
            L = L + 1
            R(I, J) = d0
            P(I, J) = d0
            do K = 1, nfreq
               R(I, J) = R(I, J) + UMat(L, K) * PosVec(K)
               P(I, J) = P(I, J) + UMat(L, K) * MomVec(K)
            end do
         end do
      end do

      ! un-mass-weight -- Position and Momentum
      do I = 1, natoms
         do J = 1, 3
            R(I, J) = R(I, J) / sqrt(T1%Particle(I)%Mass)
            P(I, J) = P(I, J) * sqrt(T1%Particle(I)%Mass)
         end do
      end do

      ! Adding the displacments to the Equilibrium structure
      ! determine velocity using the momentum
      do i = 1, natoms
         do j = 1, 3
            R(I, J) = atomcen(I, J) + R(I, J)
            Vel(I, J) = P(I, J) / T1%Particle(I)%Mass
         end do
      end do

      ! make sure that COM is at origin and remove COM velocity
      RCM = d0
      VelCM = d0
      TotMass = d0

      ! compute the c.o.m and its velocity
      do i = 1, natoms
         totMass = TotMass + T1%Particle(i)%Mass
         do j = 1, 3
            RCM(j) = RCM(j) + T1%Particle(i)%Mass * R(i, j)
            VelCM(j) = VelCM(j) + T1%Particle(i)%Mass * Vel(i, j)
         end do
      end do
      RCM = RCM / TotMass
      VelCM = VelCM / TotMass

!     Remove COM velocity
!     amv: Turned off COM removal - it can interfere with QM/MM
!     amv: momentum removal is still cool, though
!     Use Velocity to determine Momentum
      do i = 1, natoms
         do j = 1, 3
!            R(i,j)=R(i,j)-rcm(j)
            Vel(i, j) = Vel(i, j) - velcm(j)
            P(I, J) = Vel(I, J) * T1%Particle(I)%Mass
         end do
      end do

      call write_initial_condition(T1, natoms, R, P)

      ! Assign position and momentum to trajctory
      do IParticle = 1, qcNumQM
         do IDim = 1, 3
            call T1%set_pos(IParticle, IDim, R(IParticle, IDim))
            call T1%set_mom(IParticle, IDim, P(IParticle, IDim))
         end do
      end do

!     Calculate the kinetic energy and compare it to the ZPE.
!     Should be about 1/2 the ZPE.
!     We also compare the kinetic energy using Normal Modes and Cartesian Coord.
      write (fmiOut, *) 'Ekin Normal Modes:', real(Ekin)

      ! TODO: Use a common function to calculate Ekin for a Trajectory
      EKin = 0.0d0
      do I = 1, natoms
         do J = 1, 3
            EKin = EKin + P(i, j) * P(i, j) / T1%Particle(I)%Mass
         end do
      end do
      EKin = dp5 * EKin

      ZPE = get_zpe(freq, nfreq)

      write (fmiOut, *) 'ZPE is: ', real(ZPE), 'Half ZPE is: ', real(dp5 * ZPE)
      write (fmiOut, *) 'Kinetic Energy is: ', real(Ekin)
      flush (fmiOut)
   end subroutine FMS_InitialWigner

!>
!!    Samples initial conditions from a Boltzmann distribution at a constant Temperature
!!    @ingroup initial
!<
   subroutine FMS_InitialConstT(T1)
      type(T_Trajectory), intent(inout) :: T1

      real(kind=DefReal) :: atcen(qcNumQM, 3)
      real(kind=DefReal), dimension(3*qcNumQM, 3*qcNumQM) :: Umat
      real(kind=DefReal), dimension(3*qcNumQM) :: freq
      real(kind=DefReal), dimension(qcNumQM) :: Mass
      real(kind=DefReal) :: atmom(qcNumQM, 3)
      real(kind=DefReal) :: PosVec(3 * qcNumQM)
      real(kind=DefReal) :: VelVec(3 * qcNumQM)
      real(kind=DefReal) :: R(qcNumQM, 3), Vel(qcNumQM, 3)
      real(kind=DefReal) :: Ener, A, phi, TMode
      real(kind=DefReal) :: VelCM(3), RCM(3), TotMass
      real(kind=DefReal) :: ZPE, EKin, tSharp
      integer(kind=DefInt) :: I, J, K, L, natoms, nfreq, IDim, IParticle

!     Converting Temperature to AU
!     The default is that all the modes are at temperature indTemperature
!     But, the user can use the vector indFModeSharp to assign a different
!     Temprature to certain modes. (The logic here is the same as in
!     the quasi-classical action/angle variables.)

      TSharp = indTemperature * DegToAu
      indFModeSharp = indFModeSharp * DegToAu

!     Read initial positions
      do i = 1, qcNumQM
         atcen(i, 1:3) = T1%Particle(i)%get_pos()
      end do

!     now we can allocate memory
      natoms = qcNumQM

      do i = 1, natoms
         Mass(i) = T1%Particle(i)%Mass
      end do
      call read_frequencies(natoms, Mass, nfreq, freq, Umat)

!     Finally select initial conditions
!     selecting energy for each mode from a Boltzmann Distribution
!     at a temperature tSharp.  Note that we scan the vector IModeTemp
!     to see if we want to have a different temperature in certain modes
      do I = 1, nfreq
         Tmode = tSharp
         do J = 1, inNTemp
            if (inIModeSharp(J) == I) TMode = indFModeSharp(J)
         end do
         Ener = -TMode * log(1.d0 - FMS_ranb(i4zero))

!     add zpe
         Ener = Ener + freq(I) * 0.5
         phi = 2.d0 * Pi * FMS_ranb(i4zero)
         A = 2.d0 * Ener / (freq(I) * freq(I))
         PosVec(I) = -sqrt(A) * sin(phi) !Position
         VelVec(I) = -freq(I) * sqrt(A) * cos(phi) !Velocity

      end do

!     Multiply by U (transform to mass-weighted cartesians)
      L = 0
      do I = 1, natoms
         do J = 1, 3
            L = L + 1
            R(I, J) = d0
            Vel(I, J) = d0
            do K = 1, nfreq
               R(I, J) = R(I, J) + Umat(L, K) * PosVec(K)
               Vel(I, J) = Vel(I, J) + Umat(L, K) * VelVec(K)
            end do
         end do
      end do

!     un-mass-weight
      do i = 1, natoms
         do j = 1, 3
            R(I, J) = R(I, J) / sqrt(T1%Particle(I)%Mass)
            Vel(I, J) = Vel(I, J) / sqrt(T1%Particle(I)%Mass)
         end do
      end do

!     Adding the displacments, R(i,j), to the Equilibrium structure

      do i = 1, natoms
         do j = 1, 3
            R(i, j) = atcen(i, j) + R(i, j)
         end do
      end do

! make sure that COM is at origin and remove COM velocity

      RCM = d0
      VelCM = d0
      TotMass = d0

!     compute the c.o.m and its velocity
      do i = 1, natoms
         TotMass = TotMass + T1%Particle(i)%Mass
         do j = 1, 3
            RCM(j) = RCM(j) + T1%Particle(i)%Mass * R(i, j)
            VelCM(j) = VelCM(j) + T1%Particle(i)%Mass * vel(i, j)
         end do
      end do
      RCM = RCM / TotMass
      VelCM = VelCM / TotMass
!amv: turned off COM removal
      do i = 1, natoms
         do j = 1, 3
!            R(i,j)=R(i,j)-rcm(j)
            vel(i, j) = vel(i, j) - velcm(j)
         end do
      end do

!     Convert velocities to momenta
      do i = 1, natoms
         do j = 1, 3
            atmom(i, j) = Vel(i, j) * T1%Particle(i)%Mass
         end do
      end do

      call write_initial_condition(T1, natoms, R, atmom)

      do IParticle = 1, natoms
         do IDim = 1, 3
            call T1%set_pos(IParticle, IDim, R(IParticle, IDim))
            call T1%set_mom(IParticle, IDim, atmom(IParticle, IDim))
         end do
      end do

! Lets calculate the kinetic energy and compare it to the ZPE.  Should be
! about 1/2 the ZPE.

      ZPE = get_zpe(freq, nfreq)

      EKin = 0.0d0
      do I = 1, natoms
         do J = 1, 3
            EKin = EKin + vel(i, j) * vel(i, j) * T1%Particle(I)%Mass
         end do
      end do
      EKin = dp5 * EKin

      write (fmiOut, *) 'ZPE is : ', real(ZPE), 'Half ZPE is: ', real(dp5 * ZPE)
      write (fmiOut, *) 'Kinetic Energy is: ', real(Ekin)
      flush (fmiOut)

   end subroutine FMS_InitialConstT
!>
!!    The subrouinte is written as comparison to 'InitialWigner.f'
!!    @brief This routine samples initial conditions from a Husimi Distribution.
!!    For each Normal mode we sample from:
!!    rho(x,p)=exp[ {- Freq^2*(x-x0)^2 - (p-p_0)^2 } / {2*alpha_x}]
!!    where x, and p are the normal mode position and momentum, respectively.
!!    The width alpha_x = Freq  for zero temperature
!!    Note that Normal modes have units of sqrt[mass]*Length
!!    If indTemperature /=0 then the width is defined as
!!    Alpha_X=Freq/2*[1+1/Tanh(Beta*Freq/2)] where Beta=1/(Kb*T)
!!    See [ Liu and Miller, JCP vol 134, Article 104101 (2011) ]
!!    This is the result for the Husimi distribution at finite temperature
!!    For imaginary frequency problem occured in some molecules or more often
!!     in clusters and solution, the code needs to be modified with the
!!     local Gaussian approximation. See [Liu and Miller, JCP vol 131, 074113 (2009)]
!!    @ingroup initial
!<
   subroutine FMS_InitialHusimi(T1, IDNum)
      type(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: IDNum

      real(DefReal), dimension(qcNumQM) :: Mass
      real(DefReal), dimension(3*qcNumQM) :: freq

      real(DefReal), dimension(3*qcNumQM, 3*qcNumQM) :: UMat

      ! this are all back to front!
      real(DefReal), dimension(qcNumQM, 3) :: atomcen, P, R, Vel

      real(DefReal) :: Beta
      real(DefReal) :: rsq, dx1, dx2, fac, x1, x2, alphax, sigp, sigx

      real(kind=DefReal) :: PosVec(3 * qcNumQM)
      real(kind=DefReal) :: Momvec(3 * qcNumQM)

      real(kind=DefReal) :: TotMass, rhototal, RCM(3), VelCM(3)
      integer(kind=DefInt) :: I, J, K, L, natoms, nfreq, IDim, IParticle
      real(kind=DefReal) :: ZPE, Ekin

!     Read initial positions
      do i = 1, qcNumQM
         atomcen(i, 1:3) = T1%Particle(i)%get_pos()
      end do

!     now we can allocate memory
      natoms = qcNumQM

      do i = 1, natoms
         Mass(i) = T1%Particle(i)%Mass
      end do
      call read_frequencies(natoms, Mass, nfreq, freq, Umat)

!     Finally start sampling
      Ekin = 0.0d0
!     If this is the first trajectory and inzFirstGauss is true place the
!     basis function at the origin.  It will be shifted to atomcen/atvel below.
      if (inzFirstGauss .and. IDNum == 1) then
         PosVec = 0.0d0
         Momvec = 0.0d0
      else
!        Sample from a Gaussian distribution
         rhototal = 1.0d0
         do I = 1, nfreq
            Alphax = freq(i)
            if (indTemperature > 0.0d0) then
               write (fmiOut, *) 'Sampling at finite T =', indTemperature
               Beta = 1.0d0 / (indTemperature * BoltzK)
               Alphax = Alphax / 2.0d0 * (1.0d0 + 1.0d0 / tanh(Beta * freq(i) * dp5))
            end if
            Sigx = sqrt(Alphax) / freq(i)
            Sigp = sqrt(Alphax)

1           continue
            dx1 = 2.d0 * FMS_ranb(i4zero) - 1.d0
            dx2 = 2.d0 * FMS_ranb(i4zero) - 1.d0
            rsq = dx1 * dx1 + dx2 * dx2
            if (rsq >= 1.d0 .or. rsq == 0.d0) goto 1
            fac = sqrt(-2.d0 * log(rsq) / rsq)
            x1 = dx1 * fac
            x2 = dx2 * fac
! Each Normal mode is in Equilibrium in both position and
! mometum, i.e., x0 and p0 are zero so we just need to "stretch"
! the random numbers using the proper Sigma.
! Note that p is the momentum but the normal modes have a unit mass
! so we do not need to divide by the mass when we go to velocity

            PosVec(I) = Sigx * x1
            Momvec(I) = Sigp * x2
! For debugging purpose calculate Kinetic energy
            Ekin = Ekin + Momvec(i) * Momvec(I)
         end do
      end if
! Kinetic energy in Normal-Modes
      Ekin = Ekin * dp5
!     multiply by U (transform to mass-weighted Cartesian--position/momentum)
      L = 0
      do I = 1, natoms
         do J = 1, 3
            L = L + 1
            R(I, J) = d0
            P(I, J) = d0
            do K = 1, nfreq
               R(I, J) = R(I, J) + UMat(L, K) * PosVec(K)
               P(I, J) = P(I, J) + UMat(L, K) * MomVec(K)
            end do
         end do
      end do

!     un-mass-weight -- Position and Momentum

      do I = 1, natoms
         do J = 1, 3
            R(I, J) = R(I, J) / sqrt(T1%Particle(I)%Mass)
            P(I, J) = P(I, J) * sqrt(T1%Particle(I)%Mass)
         end do

      end do

!     Adding the displacments to the Equilibrium structure
!     determine velocity using the momentum
      do i = 1, natoms
         do j = 1, 3
            R(I, J) = atomcen(I, J) + R(I, J)
            Vel(I, J) = P(I, J) / T1%Particle(I)%Mass
         end do
      end do

! make sure that COM is at origin and remove COM velocity

      RCM = d0
      VelCM = d0
      TotMass = d0

! compute the c.o.m and its velocity
      do i = 1, natoms
         TotMass = TotMass + T1%Particle(i)%Mass
         do j = 1, 3
            RCM(j) = RCM(j) + T1%Particle(i)%Mass * R(i, j)
            VelCM(j) = VelCM(j) + T1%Particle(i)%Mass * Vel(i, j)
         end do
      end do
      RCM = RCM / TotMass
      VelCM = VelCM / TotMass
!     Remove COM velocity
!     amv: Turned off COM removal - it can interfere with QM/MM
!     amv: momentum removal is still cool, though
!     Use Velocity to determine Momentum
      do i = 1, natoms
         do j = 1, 3
!            R(i,j)=R(i,j)-rcm(j)
            Vel(i, j) = Vel(i, j) - velcm(j)
            P(I, J) = Vel(I, J) * T1%Particle(I)%Mass
         end do
      end do

!     Write initial conditions to a file
      call write_initial_condition(T1, natoms, R, P)

!     Assign position and momentum to trajecotry
      do IParticle = 1, qcNumQM
         do IDim = 1, 3
            call T1%set_pos(IParticle, IDim, R(IParticle, IDim))
            call T1%set_mom(IParticle, IDim, P(IParticle, IDim))
         end do
      end do
!     Lets calculate the kinetic energy and compare it to the ZPE.  Should be
!     about 1/2 the ZPE.
!     We also compare the kinetic energy using Normal Modes and Cartesian Coord.

      ZPE = get_zpe(freq, nfreq)

      write (fmiOut, *) 'Ekin Normal Modes:', real(Ekin)

      EKin = 0.0d0
      do I = 1, natoms
         do J = 1, 3
            EKin = EKin + P(i, j) * P(i, j) / T1%Particle(I)%Mass
         end do
      end do
      EKin = dp5 * EKin

      write (fmiOut, *) 'ZPE is : ', real(ZPE), 'Half ZPE is: ', real(dp5 * ZPE)
      write (fmiOut, *) 'Kinetic Energy is: ', real(Ekin)
      call flush (fmiOut)

   end subroutine FMS_InitialHusimi

!>
!!    Samples initial conditions from a Wigner distribution using a randomly chosen unit vector.
!!    @note The resulting  distribution seems to be too cold!
!!    This version is general and allows for sampling at finite indTemp.
!!    In this case the width of the Wigner distribution is given by:
!!    \f[ \alpha_x = Freq/2*\tanh(\beta*Freq/2) \f] where \f$ \beta=1/(K_bT) \f$
!!    See Heller JCP, vol, 65, pp. 1289 (1976)
!!    If indTemperature=0 then \f$ \alpha_x=Freq/2\f$ which is the known result (Note
!!    that the normal modes have unit mass)
!!    \see FMS_InitialWigner
!!    @ingroup initial
!<
   subroutine FMS_InitialColdWig(T1, IDNum)
      use EispackModule, only: FMS_ch
      type(t_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), intent(in) :: IDNum
      real(kind=DefReal) :: Beta
      integer(kind=DefInt) :: I, J, l, IParticle, IDim
      integer(kind=DefInt) :: natoms, natomsT6, NPSDim, MatZ
      real(kind=DefReal) :: atcen(qcNumQM, 3), vel(qcNumQM, 3)
      real(kind=DefReal), dimension(3*qcNumQM, 3*qcNumQM) :: Umat
      real(kind=DefReal), dimension(3*qcNumQM) :: freq
      real(kind=DefReal), dimension(qcNumQM) :: Mass
      real(kind=DefReal), dimension(6*qcNumQM) :: fv1, fv2
      real(kind=DefReal), dimension(6*qcNumQM, 6*qcNumQM) :: OReal, OImag
      real(kind=DefReal) :: big, ysq, dx1, dx2, x1, x2, alphax, sigx, sigp, prob
      real(kind=DefReal) :: ZPE, EKin
      integer(kind=DefInt) :: nfreq, k, ione, ierr
      real(kind=DefReal) :: PosVec(3 * qcNumQM), MomVec(3 * qcNumQM)
      real(kind=DefReal) :: R(qcNumQM, 3), P(qcNumQM, 3)
      real(kind=DefReal) :: RCM(3), VelCM(3), TotMass
      real(kind=DefReal) :: UVec(6 * qcNumQM)
      real(kind=DefReal) :: EValues(6 * qcNumQM)
      real(kind=DefReal) :: EVReal(6 * qcNumQM, 6 * qcNumQM)
      real(kind=DefReal) :: EVImag(6 * qcNumQM, 6 * qcNumQM)
      real(kind=DefReal) :: FM1(18 * qcNumQM)
      complex(kind=DefComp) :: CRMat(6 * qcNumQM, 6 * qcNumQM)
      complex(kind=DefComp) :: Lambda(6 * qcNumQM, 6 * qcNumQM)
      complex(kind=DefComp) :: CTmp(6 * qcNumQM, 6 * qcNumQM)
      complex(kind=DefComp) :: CEvec(6 * qcNumQM, 6 * qcNumQM)
      complex(kind=DefComp) :: CSTemp

!     Beta is 1/kBT
      Beta = d1 / (indTemperature * BoltzK)

!     Read initial positions
      do i = 1, qcNumQM
         atcen(i, 1:3) = T1%Particle(i)%get_pos()
      end do

      natoms = qcNumQM

      do i = 1, natoms
         Mass(i) = T1%Particle(i)%Mass
      end do
      call read_frequencies(natoms, Mass, nfreq, freq, Umat)

      NPSDim = nfreq * 2 !6D phase space
      natomsT6 = natoms * 6
!     Finally, select initial conditions
!     If this is the first trajectory and inzFirstGauss is true the basis fxn
!     is placed at the origin.  Below it will be shifted to atcen
      if (inzFirstGauss .and. IDNum == 1) then
         PosVec = d0
         MomVec = d0
      else
!     Need to choose a random vector on the 2N-dimensional unit hypersphere
!     Construct a random purely imaginary Hermitian matrix.
!     Then, exponentiate the matrix
!     which is guaranteed to give an orthogonal matrix, i.e. a rotation matrix.
         Big = 100000.0d0
         do i = 1, NPSDim
            do j = 1, NPSDim
               Lambda(i, j) = (0.d0, 0.d0)
            end do
         end do
         do i = 1, NPSDim
            do J = i + 1, NPSDim
               Lambda(i, j) = c1i * (Big - 2 * Big * FMS_ranb(i4zero))
               Lambda(j, i) = -Lambda(i, j)
            end do
         end do

         do i = 1, NPSDim
            do j = 1, NPSDim
               OReal(i, j) = real(Lambda(i, j))
               OImag(i, j) = real(-c1i * Lambda(i, j))
            end do
         end do
         matz = 1
         call FMS_ch(natomsT6, NPSDim, OReal, OImag, EValues, &
                     MatZ, EVReal, EVImag, FV1, FV2, FM1, IErr)

         do i = 1, NPSDim
            do j = 1, NPSDim
               CEvec(i, j) = EVReal(i, j) + c1i * EVImag(i, j)
            end do
         end do

         do i = 1, NPSDim
            CSTemp = exp(-c1i * EValues(i))
            do j = 1, NPSDim
               CTmp(i, j) = CSTemp * conjg(CEvec(j, i))
            end do
         end do

         do i = 1, NPSDim
            do j = 1, NPSDim
               CRMat(i, j) = d0
               do k = 1, NPSDim
                  CRMat(i, j) = CRMat(i, j) + CEvec(i, k) * CTmp(k, j)
               end do
            end do
         end do

         ! Now CRMat should be a random rotation matrix in the NPSDim dimensional
         ! space.  Furthermore, it should be purely real although we are
         ! storing it in a complex matrix (because it is built from complex
         ! matrix-matrix products)
         do i = 1, NPSDim
            do j = 1, NPSDim
               if (abs(aimag(CRMat(i, j))) > 1.0d-5) then
                  write (fmiOut, *) 'Ouch ', i, j, CRMat(i, j)
               end if
            end do
         end do

!    Choose a random unit vector and rotate by the random rotation matrix
!    make sure that IOne is in the correct range: 1-NPSDim
!    we add 1.5 and not one in order to avoid a situation where IOne=NPSDim
!    only when ranb=1.0d0.  Once we have chosen a unit vector randomly, we
!    rotate it with our random rotation matrix.
1        continue
         IOne = int(dble(NPSDim - 1) * FMS_ranb(i4zero) + 1.50d0)
         if (IOne > NPSDim .or. IOne == 0) go to 1
!    Matrix-vector multiplication CRMat*UnitVector
         do i = 1, NPSDim
            UVec(i) = real(CRMat(i, IOne))
         end do

!    Check that we actually have a unit vector on the hypersurface
         Ysq = d0
         do i = 1, NPSDim
            Ysq = Ysq + UVec(i) * UVec(i)
         end do
         if (abs(Ysq - d1) > 1.0d-04) then
            write (fmiOut, *) 'Warning: Vector has non-unit norm?'
            write (fmiOut, *) Ysq
         end if

!    Choose a random number from a Gaussian distribution we actually choose
!    two and use only one. Scale the randomly oriented vector by a random
!    number drawn from a standard normal deviate.  We can do this because
!    we will scale by frequency later.
         dx1 = FMS_ranb(i4zero)
         dx2 = FMS_ranb(i4zero)
         x1 = sqrt(-2.d0 * log(dx1)) * cos(2.d0 * pi * dx2)
         x2 = sqrt(-2.d0 * log(dx1)) * sin(2.d0 * pi * dx2)
         do I = 1, NPSDim
            UVec(i) = UVec(i) * x1
         end do

!    Wigner Distribution. For each Normal mode we sample from:
!    rho(x,p)=exp[-2*AlphaX(x-x0)^2 -(p-p0)^2/(2*AlphaX)]
!    note that p is momentum and therefore we divide by the mass
!    of the normal mode to get the velocity
!    The width AlphaX is Mass*Freq/2 [Heller JCP vol 80 p 5038 (1984)]
!    But Normal modes have units of sqrt[mass]*Length and therefore
!    alphax=Freq/2!
         do I = 1, nfreq
            Alphax = freq(i) / 2.d0
            if (indTemperature /= d0) then
               Alphax = Alphax * tanh(Beta * freq(i) * dp5)
            end if
            Sigx = sqrt(1.d0 / (4.d0 * Alphax))
            Sigp = sqrt(1.d0 * Alphax)
! Each Normal mode is in Equilibrium in both position and
!!    momentum, i.e., x0 and p0 are zero so we just need to "stretch"
! the random numbers using the proper Sigma.
! Note that p is the momentum but the normal modes have a unit mass
! so we do not need to divide by the mass when we go to velocity
            PosVec(I) = Sigx * UVec(2 * I - 1)
            MomVec(I) = Sigp * UVec(2 * I)
         end do
      end if

! Compute the probability of this Initial condition
      Prob = d1
      do I = 1, nfreq
         Alphax = freq(i) / 2.d0
         Prob = Prob * exp(-2.d0 * Alphax * PosVec(I) * PosVec(I) - &
                           MomVec(I) * MomVec(I) / (2.d0 * Alphax))
      end do
      write (fmiOut, *) 'Probability of These Initial Conditions: ', Prob
      flush (fmiOut)
!!    transform to mass weighted cartesian coordinates
      L = 0
      do I = 1, natoms
         do J = 1, 3
            L = L + 1
            R(I, J) = d0
            P(I, J) = d0
            do K = 1, nfreq
               R(I, J) = R(I, J) + UMat(L, K) * PosVec(K)
               P(I, J) = P(I, J) + UMat(L, K) * MomVec(K)
            end do
         end do
      end do

! un-mass weight, position and momentum

      do i = 1, natoms
         do j = 1, 3
            r(i, j) = r(i, j) / sqrt(T1%Particle(I)%Mass)
            p(i, j) = p(i, j) * sqrt(T1%Particle(I)%Mass)
         end do
      end do

!    Adding the displacments to the Equilibrium structure
!    Determine Velocity from Momentum

      do i = 1, natoms
         do j = 1, 3
            R(I, J) = atcen(I, J) + R(I, J)
            Vel(I, J) = P(I, J) / T1%Particle(I)%Mass
         end do
      end do

! make sure that COM is at origin and remove COM velocity

      RCM = d0
      VelCM = d0
      TotMass = d0

! compute the c.o.m and its velocity
      do i = 1, natoms
         TotMass = TotMass + T1%Particle(i)%Mass
         do j = 1, 3
            RCM(j) = RCM(j) + T1%Particle(i)%Mass * atcen(i, j)
            VelCM(j) = VelCM(j) + T1%Particle(i)%Mass * vel(i, j)
         end do
      end do
      RCM = RCM / TotMass
      VelCM = VelCM / TotMass
!     amv: Do not remove COM - it interferes with QM/MM
      do i = 1, natoms
         do j = 1, 3
!            R(i,j)=R(i,j)-rcm(j)
            vel(i, j) = vel(i, j) - velcm(j)
            P(i, j) = vel(i, j) * T1%Particle(i)%Mass
         end do
      end do

      ! Write initial conditions to a file
      call write_initial_condition(T1, natoms, R, P)

      do iParticle = 1, qcNumQM
         do iDim = 1, 3
            call T1%set_pos(iParticle, iDim, R(iParticle, IDim))
            call T1%set_mom(iParticle, iDim, P(iParticle, IDim))
         end do
      end do

      ! Lets calculate the kinetic energy and compare it to the ZPE.  Should be
      ! about 1/2 the ZPE.
      ZPE = get_zpe(freq, nfreq)

      EKin = 0.0d0
      do I = 1, natoms
         do J = 1, 3
            EKin = EKin + P(i, j) * P(i, j) / T1%Particle(I)%Mass
         end do
      end do
      EKin = dp5 * EKin

      write (fmiOut, *) 'ZPE is : ', real(ZPE), 'Half ZPE is: ', real(dp5 * ZPE)
      write (fmiOut, *) 'Kinetic Energy is: ', real(Ekin)
      flush (fmiOut)
   end subroutine FMS_InitialColdWig

!>
!!    Samples initial conditions from a quasiclassical
!!    distribution (I.e. using action-angle variables)
!!    @ingroup initial
!<
   subroutine FMS_InitialQuasi(T1)
      type(T_Trajectory), intent(inout) :: T1

      integer(kind=DefInt) :: natoms, nfreq, iGUnit
      real(kind=DefReal), dimension(3*qcNumQM) :: freq, xmass
      real(kind=DefReal) :: atcen(3 * qcNumQM)
      real(kind=DefReal), dimension(3*qcNumQM, 3*qcNumQM) :: Umat
      real(kind=DefReal), dimension(qcNumQM) :: Mass
      real(kind=DefReal) :: atmom(3 * qcNumQM)
      real(kind=DefReal) :: atvel(3 * qcNumQM)
      real(kind=DefReal) :: qprime(3 * qcNumQM)
      real(kind=DefReal) :: pprime(3 * qcNumQM)
      real(kind=DefReal) :: qsave(3 * qcNumQM)
      real(kind=DefReal) :: tmpq(3 * qcNumQM)
      real(kind=DefReal) :: tmpp(3 * qcNumQM), tmpe, zpe, fac, targete
      real(kind=DefReal) :: RCM(3), VelCM(3), TotMass
      integer(kind=DefInt) :: i, j, k
      real(kind=DefReal) :: x
      character(len=:), allocatable :: output_file

!     Read initial positions
      do i = 1, qcNumQM
         atcen((i - 1) * 3 + 1:i * 3) = T1%Particle(i)%get_pos()
      end do

!     now we can allocate memory
      natoms = qcNumQM

      do i = 1, natoms
         Mass(i) = T1%Particle(i)%Mass
      end do
      call read_frequencies(natoms, Mass, nfreq, freq, Umat)

!     Finally select initial conditions
!     select random angle for each mode
      do i = 1, nfreq
         fac = d1
         do j = 1, inNAddQuanta
            if (i == inIAddQuanta(j)) fac = fac + 2.0d0
         end do
         fac = sqrt(fac)
         x = FMS_ranb(i4zero) * d2 * pi
         qprime(i) = cos(x) * fac
         qsave(i) = qprime(i)
         pprime(i) = -1.0d0 * sin(x) * fac
         tmpq(i) = qprime(i) / sqrt(freq(i)) !Position
         tmpp(i) = sqrt(freq(i)) * pprime(i) !Velocity
      end do

!     multiply by U (transform to mass-weighted cartesians)
      do i = 1, natoms * 3
         qprime(i) = d0
         pprime(i) = d0
         do j = 1, nfreq
            qprime(i) = qprime(i) + UMat(i, j) * tmpq(j)
            pprime(i) = pprime(i) + UMat(i, j) * tmpp(j)
         end do
      end do

      xmass = 0.0d0
      k = 1
      do i = 1, natoms
         do j = 1, 3
            xmass(k) = sqrt(Mass(i))
            k = k + 1
         end do
      end do

!     un-mass-weight
      do i = 1, natoms * 3
         qprime(i) = qprime(i) / xmass(i)
         atcen(i) = atcen(i) + qprime(i) !add to equil. position
         atvel(i) = pprime(i) / xmass(i) !velocity
         atmom(i) = pprime(i) * xmass(i) !momentum
      end do
!     make sure that COM is at origin and remove COM velocity

      RCM = d0
      VelCM = d0
      TotMass = d0

!     compute COM and its velocity

      do i = 1, natoms
         TotMass = TotMass + T1%Particle(i)%Mass
         do j = 1, 3
            RCM(j) = RCM(j) + T1%Particle(i)%Mass * atcen((i - 1) * 3 + j)
            VelCM(j) = VelCM(j) + T1%Particle(i)%Mass * atvel((i - 1) * 3 + j)
         end do
      end do
      RCM = RCM / TotMass
      VelCM = VelCM / TotMass

!     remove COM and its velocity
!amv     turned off COM removal - it interferes with QM/MM
      do i = 1, natoms
         do j = 1, 3
!           atcen((i-1)*3+j)=atcen((i-1)*3+j)-RCM(j)
            atvel((i - 1) * 3 + j) = atvel((i - 1) * 3 + j) - VelCM(j)
            atmom((i - 1) * 3 + j) = atmom((i - 1) * 3 + j) - VelCM(j) * T1%Particle(I)%Mass
         end do
      end do

      ! Write initial conditions to file
      ! TODO: call 'write_initial_condition' subroutine
      output_file = trim(FMSWorkingDir)//'Out.xyz'
      open (newunit=IGUnit, file=output_file, action='write')
      write (IGUnit, '(A10)') "UNITS=BOHR"
      write (IGUnit, *) natoms
      do I = 1, natoms
         write (IGUnit, 1005) T1%Particle(I)%Elmnt, (atcen((I - 1) * 3 + J), J=1, 3)
      end do
      write (IGUnit, '(A9)') "# momenta"
      do I = 1, natoms
         write (IGUnit, 1006) (atmom((I - 1) * 3 + J), J=1, 3)
      end do

      close (IGUnit)
1005  format(A2, 3es24.15)
1006  format(3es24.15)

      zpe = d0
      tmpe = d0

      do i = 1, nfreq
         zpe = zpe + freq(i)
      end do

      do i = 1, natoms * 3
         tmpe = tmpe + atmom(i) * atmom(i) / (xmass(i) * xmass(i))
      end do

      zpe = zpe * dp5
      tmpe = tmpe * dp5

      write (fmiOut, *) 'ZPE is: ', real(zpe), 'Half the ZPE is: ', real(dp5 * ZPE)
      write (fmiOut, *) 'Kinetic Energy is:', real(tmpe)
      do j = 1, nfreq
         tmpq(j) = d0
         tmpq(j + 1) = d0
         do i = 1, natoms * 3
            tmpq(j) = tmpq(j) + (qprime(i) * xmass(i) * UMat(i, j))
            tmpq(j + 1) = tmpq(j + 1) + UMat(i, j) * UMat(i, j)
         end do
         tmpq(j) = tmpq(j) * sqrt(freq(j))
         tmpe = tmpe + dp5 * freq(j) * tmpq(j) * tmpq(j)
      end do
      targete = zpe
      do i = 1, inNAddQuanta
         targete = targete + freq(inIAddQuanta(i))
      end do

!      write(fmiOut,*) 'ZPE, target E, estimated E'
!      write(fmiOut,1004) zpe, targete, tmpe

      do i = 1, natoms
         do j = 1, 3
            call T1%set_pos(i, j, atcen((i - 1) * 3 + j))
            call T1%set_mom(i, j, atmom((i - 1) * 3 + j))
         end do
      end do

   end subroutine FMS_InitialQuasi

   !>
   !!    Reads in Frequencies from whatever file format is present
   !!    @ingroup initial
   !<
   subroutine read_frequencies(natoms, mass, nfreq, freq, Umat)
      use EispackModule, only: FMS_rs
      integer(kind=DefInt), intent(in) :: natoms
      real(DefReal), intent(in) :: mass(natoms)
      integer(kind=DefInt), intent(out) :: nfreq

      real(DefReal), intent(out) :: Umat(3 * natoms, 3 * natoms)
      real(DefReal), intent(out) :: freq(3 * natoms)

      real(DefReal), dimension(3*natoms) :: xmass, fv1, fv2, vfreqs
      integer(DefInt), dimension(3*natoms) :: jfreq
      real(DefReal), dimension(3*natoms, 3*natoms) :: vmat

      real(DefReal), allocatable :: OReal(:, :)

      real(kind=DefReal) :: dmass, fac
      integer(kind=DefInt) :: natom3, IFUnit, ios
      integer(kind=DefInt) :: natomsin, i, j, k, ist, lst, ix, idof, ierr, ifreq
      character(len=256) :: filein, fileout, cbuffer
      character(len=10) :: method
      logical :: NeedStretch, file_exists

      natom3 = natoms * 3
      method = ''
      nfreq = 0
      UMat = 0.0d0
      freq = 0.0d0

! First we'll try and find what sort of Hessian information we have
! based on the available files
! The order of priority is:
! Hessian.dat < FrequenciesMP.dat < Frequencies.dat

! Regular-formatted Hessian
      inquire (file=trim(FMSWorkingDir)//'Hessian.dat', exist=file_exists)
      if (file_exists) then
         method = 'HESS'
      end if

! Molpro-formatted Frequencies
      inquire (file=trim(FMSWorkingDir)//'FrequenciesMP.dat', exist=file_exists)
      if (file_exists) then
         method = 'FREQMP'
      end if

! Regular-formatted Frequencies
      inquire (file=trim(FMSWorkingDir)//'Frequencies.dat', exist=file_exists)
      if (file_exists) then
         method = 'FREQ'
      end if

      select case (method)

      case ("HESS")
         filein = trim(FMSWorkingDir)//'Hessian.dat'
         open (newunit=ifunit, file=filein, status='old', action='read', iostat=ios)
         write (fmiOut, *) 'Using Hessian.dat to find normal mode frequencies'
         read (IFUnit, *, end=999) natomsin
         if (natomsin /= natoms) then
            write (fmiOut, *) 'Wrong Number of atoms in Hessian.dat', natoms, natomsin
            call FMS_DieError('Error in ReadFreq')
         end if
         do i = 1, natom3
            read (IFUnit, *, end=999) (vmat(i, j), j=1, natom3)
         end do
         close (IFUnit)

         do i = 1, natoms
            dmass = sqrt(Mass(i))
            dmass = MassToAu / dmass !Invert mass and convert to amu
            vmat = vmat * dmass !Scale Hessian
         end do

! diagonalize mass-weighted hessian with RS subroutine
         call FMS_rs(natom3, natom3, vmat, vfreqs, 1, Umat, fv1, fv2, ierr)

!      fac=(sqrt(hartree/tokg))/a0/(2.0d0*pi*clight)
         fac = 5140.4548d0

         nfreq = 0
         do i = 1, natom3
            if (vfreqs(i) > 1.e-4) then
               nfreq = nfreq + 1
               jfreq(nfreq) = i
               freq(nfreq) = vfreqs(i)
               freq(nfreq) = sqrt(freq(nfreq)) * fac
            end if
         end do
! now we are going back to normal-modes without mass-weighting..
         do i = 1, natom3
            dmass = Mass(ceiling(dble(i) / 3.0d0)) !yuk
            dmass = dmass / MassToAu
            do j = 1, natom3
               Umat(i, j) = Umat(i, j) / sqrt(dmass)
            end do
         end do

         fileout = trim(FMSWorkingDir)//'Frequencies.dat'
         open (newunit=ifunit, file=fileout, status='new', action='write')
         write (IFUnit, *) natoms, nfreq
         write (IFUnit, *) (freq(i), i=1, nfreq)
         do i = 1, nfreq
            write (IFUnit, *) (Umat(j, jfreq(i)), j=1, natom3)
         end do
         close (IFUnit)

      case ("FREQMP")
         write (fmiOut, *) 'Found FrequenciesMP.dat'
         filein = trim(FMSWorkingDir)//'FrequenciesMP.dat'
         open (newunit=ifunit, file=filein, status='old', action='read')

         ist = 1
         read (IFUnit, *, end=999) nfreq !number of modes
         if (nfreq /= (natom3 - 6)) then
            call FMS_DieError('FrequenciesMP.dat: Number of frequencies =/= 3N-6')
         end if
         do
            read (IFUnit, '(a256)', end=999) cbuffer !Mode label
            read (IFUnit, '(a256)', end=999) cbuffer !Frequencies
            lst = min(ist + 4, nfreq)
            read (cbuffer(24:), '(5f12.2)', end=999) (freq(ix), ix=ist, lst)
            read (IFUnit, '(a256)', end=999) cbuffer !IR Intensities
            read (IFUnit, '(a256)', end=999) cbuffer !IR Intensities
            do idof = 1, natom3
               read (IFUnit, '(a256)', end=999) cbuffer
               read (cbuffer(24:), *, end=999) (Umat(idof, ix), ix=ist, lst) !Modes
            end do
            if (lst == nfreq) exit ! We're done
            ist = ist + 5
            read (IFUnit, '(a256)', end=999) cbuffer !blank line
         end do
         close (IFUnit)
         !Now write out to Frequencies.dat in standard format
         fileout = trim(FMSWorkingDir)//'Frequencies.dat'
         open (newunit=ifunit, file=fileout, status='new', action='write')
         write (IFUnit, '(2i10)') natoms, nfreq
         write (IFUnit, '(4f20.10)') (freq(i), i=1, nfreq)
         do ifreq = 1, nfreq
            write (IFUnit, '(4f20.10)') (Umat(idof, ifreq), idof=1, natom3)
         end do
         close (IFUnit)

      case ("FREQ")
         ! Handled below
         write (fmiOut, *) 'Found Frequencies.dat'

      case default

         call FMS_DieError("No frequencies/Hessian found")

      end select

      !Now read from Frequencies.dat - note, this is always done, even if
      !we started with Hessian.dat or some other format

      filein = trim(FMSWorkingDir)//'Frequencies.dat'
      open (newunit=IFUnit, file=filein, status='old', action='read')
      read (IFUnit, *, end=999) natomsin, nfreq
      if (natomsin /= natoms) then
         write (fmiOut, '(a)') 'Invalid number of atoms in Frequencies.dat'
         write (fmiOut, '(a,i0,a,i0)') 'Expected: ', natoms, ' got: ', natomsin
         call FMS_DieError('Invalid number of atoms in Frequencies.dat')
      end if
! allocate some temporary memory

      allocate (Oreal(nfreq, nfreq))
      freq = 0.0d0
      Umat = 0.0d0
      read (IFUnit, *, end=999) (freq(I), I=1, nfreq)
      do I = 1, nfreq
         freq(I) = freq(I) * CMToAu
      end do
      do I = 1, nfreq
         read (IFUnit, *, end=999) (Umat(J, I), J=1, natom3)
      end do

      close (IFUnit)

! Now for something new.
! We do not know if U should be stretched by the mass or not.
! Let's first check it without mass stretching.  If we fail -- mass stretch
! If we fail here -- code should stop

!     Build vector containing M^1/2
      xmass = 0.0d0
      k = 1
      do i = 1, natoms
         do j = 1, 3
            xmass(k) = sqrt(Mass(i))
            k = k + 1
         end do
      end do

      NeedStretch = .false.
!Calculate U'U without mass stretch
      do i = 1, nfreq
         do j = 1, nfreq
            OReal(i, j) = 0.0d0
            do k = 1, natom3
               OReal(i, j) = OReal(i, j) + UMat(k, i) * UMat(k, j)
            end do
         end do
      end do

! Check that U'U == I
      do i = 1, nfreq
         if (abs(OReal(i, i) - 1.0d0) > 1.0d-01) then
            NeedStretch = .true.
            write (fmiOut, *) 'Warning -- Mass scaling problem detected'
            write (fmiOut, *) 'Program will Mass Stretch U and try again'
            exit
         end if
      end do
      if (NeedStretch) then
!     The matrix U needs to be "stretched" by sqrt(mass)
         do i = 1, natom3
            do j = 1, nfreq
               UMat(i, j) = (xmass(i) * UMat(i, j)) / sqrt(MassToAu)
            end do
         end do

!     calculate U'U
         do i = 1, nfreq
            do j = 1, nfreq
               OReal(i, j) = 0.0d0
               do k = 1, natom3
                  OReal(i, j) = OReal(i, j) + UMat(k, i) * UMat(k, j)
               end do
            end do
         end do

!        Check that U'U == I

         do i = 1, nfreq
            if (abs(OReal(i, i) - 1.0d0) > 1.0d-01) then
               write (fmiOut, '(a)') 'ERROR: Mass scaling problem detected'
               write (fmiOut, '(i0,a,E23.16)') i, 'th element of UtU is not 1: ', OReal(i, i)
               call FMS_DieError('Could not fix mass-scaling')
            end if
         end do

         write (fmiOut, *) 'Mass Scaling Problem Solved'
      end if

      flush (fmiOut)

      deallocate (Oreal)

      return
999   call fms_dieerror("Reached end of file "//trim(filein))

   end subroutine read_frequencies

   ! Write initial positions and momenta to a fileing positions to a file
   ! TODO: Write this as a standard (extended) XYZ format
   subroutine write_initial_condition(T, natoms, R, P)
      use GlobalModule
      type(t_Trajectory), intent(in) :: T
      real(kind=DefReal), dimension(natoms, 3), intent(in) :: R, P
      integer, intent(in) :: natoms
      character(len=:), allocatable :: output_file
      integer :: i, j, iGUnit

      output_file = trim(FMSWorkingDir)//'Out.xyz'
      open (newunit=IGUnit, file=output_file, action='write')

      write (IGUnit, '(a)') "UNITS=BOHR"
      write (IGUnit, *) natoms
      do I = 1, natoms
         write (IGUnit, '(a,3es24.15)') T%Particle(I)%Elmnt, (R(I, J), J=1, 3)
      end do
      write (IGUnit, '(a)') "# momenta"
      do I = 1, natoms
         write (IGUnit, '(3ES24.15)') (P(I, J), J=1, 3)
      end do

      close (IGUnit)
   end subroutine write_initial_condition

   ! Calculate zero-point energy from a list of normal mode frequencies
   function get_zpe(freq, nfreq) result(ZPE)
      real(DefReal), intent(in) :: freq(nfreq)
      integer(defInt), intent(in) :: nfreq
      real(defReal) :: ZPE
      integer :: i

      ZPE = 0.0d0
      do i = 1, nfreq
         ZPE = ZPE + freq(i)
      end do

      ZPE = 0.5d0 * ZPE
   end function get_zpe

   subroutine FMS_InitialSwarm(B1)
      use GlobalModule, only: DefComp, fmiOut, glzCentroids
      use FMSModule, only: FMS_ReadGeometry
      use TrajectoryModule
      use BundleModule
      use BundleCalcsModule, only: FMS_Norm, FMS_bsinvmat, FMS_BuildHS, &
                                   FMS_Set_Amplitude, FMS_UpdateCentroid
      use OverlapModule, only: overlap
      type(T_TrajectoryBundle), intent(inout) :: B1

      type(T_Trajectory) :: T_init ! trajectory from Geometry.dat

      integer(DefInt) :: ntraj, n, m, nstate, natom

      real(DefReal), dimension(3*B1%NumParticles) :: X, dX, P, dP
      real(DefReal), parameter :: SIGMA = 0.1d0 ! size of random perturbation
      real(defReal) :: P_norm ! norm of momentum
      real(defReal) :: B1_norm
      real(defReal) :: dspl   ! displacement in x direction of X, for tests w/ SOC model
      complex(DefComp), dimension(B1%NumTraj) :: S_if ! overlap between initial and final

      write (fmiout, *) "Initializing a swarm of trajectories"

      ntraj = B1%NumTraj
      natom = B1%NumParticles
      nstate = B1%NumStates

      call T_init%create(natom, nstate)
      T_init = B1%Trajectory(1)

      ! 1. get the intial momenta and positions from Geometry.dat
      call FMS_ReadGeometry(T_init)
      X = T_init%get_pos()
      P = T_init%get_mom()

      P_norm = sqrt(sum(P**2))

      ! 2. set up the Bundle trajectories

      ! -- Added just for usage with SOC model and sepSS testing: ---------------------------------------
      if (gliModel == FMSZERO) then
         write (fmiout, *) "Initialising for ToyModel in 1D"
         dX = [0.5d0, 0.d0, 0.d0]
         dP = [0.d0, 0.d0, 0.d0]
         write (fmiout, *) "dX, dP: ", dX, dP

         if (ntraj > 2) then
            !dspl = 0.1
            dspl = -0.25
            do n = 2, ntraj
               B1%Trajectory(n)%TrajID = n
               call B1%Trajectory(n)%set_pos(X + dX)
               call B1%Trajectory(n)%set_mom(P + dP) ! useless here because dP is zero, but for completeness...
               dX = [10.0d0, 0.d0, 0.d0]             ! used for placing two Gaussians in front of, one after
                                                     ! the crossing point
               !dX = [-0.25d0, 0.d0, 0.d0]            ! for having three Gaussians in front of crossing
               ! ------------------------------------------------------------------------------------------
               !dX = dX + [dspl, 0.d0, 0.d0]          ! for having many Gaussians floating around
               !dspl = dspl + 0.1                     ! somewhere and somehow
               !dX = dX + [dspl, 0.d0, 0.d0]          ! Used for comp with full GAIMS
               !dspl = dspl - 0.25
            end do
         else
            call B1%Trajectory(ntraj)%set_pos(X + dX)
            call B1%Trajectory(ntraj)%set_mom(P + dP) ! useless here because dP is zero, but for completeness...
         end if

      ! -------------------------------------------------------------------------------------------------

      else
         do n = 1, ntraj

            call random_number(dX)
            dX = sigma * (2.0d0 * dX - 1.d0)

            call random_number(dP)
            dP = sigma / P_norm * (2.0d0 * dP - 1.d0)

            B1%Trajectory(n)%TrajID = n
            call B1%Trajectory(n)%set_pos(X + dX)
            call B1%Trajectory(n)%set_mom(P + dP)
         end do
      end if

      ! 3. update centroids
!bfec
      if (glzCentroids) then
         do n = 2, ntraj
            do m = 1, n - 1
               call FMS_UpdateCentroid(B1%Trajectory(n), B1%Trajectory(m), &
                                       B1%Centroids(((n - 2) * (n - 1)) / 2 + m))
            end do
         end do
      end if

      do n = 1, B1%NumTraj
         write (fmiout, *) "Before normalising Bundle: Amplitude of traj ", n, " :", B1%Trajectory(n)%Amplitude
      end do

      !! normalize the Bundle
      !B1_norm = FMS_Norm(B1)
      !! TODO: This should be a bundle method
      !do n = 1, ntraj
      !   B1%Trajectory(n)%Amplitude = B1%Trajectory(n)%Amplitude / sqrt(B1_norm)
      !end do

      !! 3. work out the overlaps between the Bundle trajectories
      !call FMS_BuildHS(B1)

      !! 4. overlap the intial Gaussian from Geometry.dat onto the Bundle
      !do n = 1, ntraj
      !   S_if(n) = overlap(T_init, B1%Trajectory(n))
      !end do

      !! careful here, matmul( A, B ) = A^* . B
      !S_if = matmul(conjg(FMS_bSInvMat(B1)), S_if)

      !! 5. set the amplitudes
      !call FMS_Set_Amplitude(B1, S_if)
      call T_init%destroy()

      B1_norm = FMS_Norm(B1)
      !write (fmiOut, *) "Overlap between Initial Trajectory and the Bundle"
      !write (fmiOut, *) sqrt(B1_norm)

      ! normalize the Bundle again
      do n = 1, ntraj
         B1%Trajectory(n)%Amplitude = B1%Trajectory(n)%Amplitude / sqrt(B1_norm)
      end do

      do n = 1, B1%NumTraj
         write (fmiout, *) "After normalising Bundle: Amplitude of traj ", n, " :", B1%Trajectory(n)%Amplitude
      end do
   end subroutine FMS_InitialSwarm

end module SamplingModule
