! Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Model potentials in the adiabatic basis. Includes:
!!
!!      (gliModel==1) Tully 1-D avoided crossing
!!            Tully, JCP 93, 1061 (1990).
!!
!!      (gliModel==2) Persico 2-D 2-state CI (in stretched coordinates)
!!            Ferretti et al, J CP 104, 5517 (1996).
!!
!!      (gliModel==3) Izmaylov 2-D 2-state LVC
!!            Ryabinkin et al, J CP 140, 214116 (2014).
!!
!!      (gliModel==4) GAIMS Model 1D 2-state SOC
!!            Granucci et al, J CP 137, 22A501 (2012).
!!
!!    All dimensions higher than required for the model potential
!!    are ignored, and their forces and derivative coupling set to 0
!!
!!    Surfaces here are coded in the *diabatic* basis; they are then
!!    diagonalized, and the derivative coupling and force evaluated
!!    via finite differences.
!<
module ToyModelModule
   use GlobalModule
   use TrajectoryModule
   use ElecStrucModule
   implicit none
   private

   real(kind=DefReal), parameter :: dx = 0.005 !finite difference step

   ! Parameters for Tully 1
   real(kind=DefReal), parameter :: A = 0.01
   real(kind=DefReal), parameter :: B = 1.6
   real(kind=DefReal), parameter :: C = 0.005
   real(kind=DefReal), parameter :: D = 1.0

   ! Parameters for Persico
   real(kind=DefReal), parameter :: alpha = 3.d0
   real(kind=DefReal), parameter :: beta = 1.5d0
   real(kind=DefReal), parameter :: gamma = 0.08
   real(kind=DefReal), parameter :: delta = 0.01
   real(kind=DefReal), parameter :: KX = 0.02
   real(kind=DefReal), parameter :: KY = 0.10
   real(kind=DefReal), parameter :: X1 = 4.d0
   real(kind=DefReal), parameter :: X2 = 3.d0
   real(kind=DefReal), parameter :: X3 = 3.d0

!-----------------------------------------------
!
!  Parameters for Izmaylov:
!     * W1 ....... frequency in x-direction (curvature of parabola)
!     * W2 ....... frequency in y-direction (curvature of parabola)
!     * XA ....... shift of parabola in x-direction
!     * YA ....... shift of parabola in y-direction
!     * deltaE ... energy difference between parabolas
!     * coupC .... linear coupling in y-direction
!
!-----------------------------------------------
   type :: t_izmaylov_params
      real(kind=DefReal) :: W1
      real(kind=DefReal) :: W2
      real(kind=DefReal) :: XA
      real(kind=DefReal) :: YA
      real(kind=DefReal) :: deltaE
      real(kind=DefReal) :: coupC
   contains
      procedure, public :: initialize => initialize_izmaylov_params
   end type t_izmaylov_params

!  Parameters for GAIMS_model
!  https://doi.org/10.1063/1.4707737
!-----------------------------------------------
   real(kind=DefReal), parameter :: a_1 = 0.03452d0
   real(kind=DefReal), parameter :: a_2 = 0.5d0
   real(kind=DefReal), parameter :: alpha_1 = 0.35d0
   real(kind=DefReal), parameter :: alpha_2 = 0.25d0
   real(kind=DefReal), parameter :: dE = 0.04d0
   real(kind=DefReal), parameter :: dr_sigma = 2.d0
   real(kind=DefReal), parameter :: c_0 = 0.001d0
   complex, parameter :: c_1 = (0.0005, 0.0005)

   type :: t_GAIMS_params
      !> controlls the position of SOC sign change
      real(kind=DefReal) :: r_sigma
   contains
      procedure, public :: initialize => initialize_gaims_params
   end type t_GAIMS_params
!-----------------------------------------------

   type(t_izmaylov_params) :: izmaylov_params
   type(t_GAIMS_params) :: GAIMS_params

   public :: FMS_ToyModel, izmaylov_params, GAIMS_params

contains

   subroutine initialize_izmaylov_params(self, W1, W2, XA, YA, deltaE, coupC)
      class(t_izmaylov_params), intent(inout) :: self
      real(kind=DefReal), intent(in) :: W1, W2, XA, YA, deltaE, coupC
      self%W1 = W1
      self%W2 = W2
      self%XA = XA
      self%YA = YA
      self%deltaE = deltaE
      self%coupC = coupC
   end subroutine initialize_izmaylov_params

   subroutine initialize_gaims_params(self, r_sigma)
      class(t_GAIMS_params), intent(inout) :: self
      real(kind=DefReal), intent(in) :: r_sigma

      self%r_sigma = r_sigma
   end subroutine initialize_gaims_params

   subroutine FMS_ToyModel(T1)
      type(T_Trajectory), intent(inout) :: T1

      real(kind=DefReal) :: Coupling(T1%NumStates - 1, T1%NumDimensions)
      real(kind=DefReal) :: x, y, H_Diab(2, 2)
      real(kind=DefReal), dimension(2) :: Evec1, Evec2, EDisplace1, &
                                          EDisplace2, EvecTemp
      real(kind=DefReal) :: v11, v12, v21, v22
      real(kind=DefReal) :: dv11dx, dv12dx, dv22dx
      real(kind=DefReal) :: dv11dy, dv12dy, dv22dy
      real(kind=DefReal) :: MASSY, MASSX
      real(kind=DefReal) :: sigma_G
      real(kind=DefReal) :: W1, W2, XA, YA, coupC

      select case (gliMethod)

      case (1) ! Tully 1
         x = T1%particle(1)%get_pos(1)
!         T1%ElecStruc%PotEn(2) = 0.5d0*exp(-0.25d0*x)
!         T1%ElecStruc%PotEn(1) = 0.03452d0*exp(-0.35d0*x) + 0.04d0
         T1%ElecStruc%PotEn(2) = &
            0.237236d0 * (1 - exp(-3.45596d0 * (x - 0.75268d0)))**2 + 0.0364684d0
         T1%ElecStruc%PotEn(1) = &
            0.261047d0 * (1 - exp(-3.33327d0 * (x - 0.721627d0)))**2 + 0.0182831d0

         T1%esflags%zpotencurrent = .true.
!         write(fmiOut,*) x,T1%ElecStruc%PotEn(1)

         ! Force
         T1%ElecStruc%DerivMat = 0.d0
         if (T1%StateID == 2) then
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = &
               2.d0 * 0.237236d0 * (1 - exp(-3.45596d0 * (x - 0.75268d0))) * &
               3.45596d0 * exp(-3.45596d0 * (x - 0.75268d0))
         else
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = &
               2.d0 * 0.261047d0 * (1 - exp(-3.33327d0 * (x - 0.721727d0))) * &
               3.33327d0 * exp(-3.33327d0 * (x - 0.721627d0))
         end if
         T1%ESFlags%ZDerivCurrent(T1%StateID, T1%StateID) = .true.
         Coupling = 0.d0
         T1%ESFlags%zDerivCurrent = .true.

!        SOC part

         T1%ElecStruc%SOMat = (0.d0, 0.d0)
!       call tully(x,sigma)
         T1%ElecStruc%SOMat(1, 2, 2, 2) = (0.d0, 1.d0) * 0.000322182 * 10
         T1%ElecStruc%SOMat(2, 1, 2, 2) = (0.d0, -1.d0) * 0.000322182 * 10
         T1%ElecStruc%SOMat(1, 2, 2, 1) = (0.5d0, -0.5d0) * 0.000322182 * 10
         T1%ElecStruc%SOMat(2, 1, 1, 2) = (0.5d0, 0.5d0) * 0.000322182 * 10
         T1%ElecStruc%SOMat(1, 2, 2, 3) = (0.5d0, 0.5d0) * 0.000322182 * 10
         T1%ElecStruc%SOMat(2, 1, 3, 2) = (0.5d0, -0.5d0) * 0.000322182 * 10

         T1%ESFlags%zSOMCurrent = .true.

!!!     Diabatic surfaces (v22=-v11, v21=v12)
!!         x = T1%Particle(1)%get_pos(1)
!!         call Tully(x,H_Diab)
!!         if(x <= 0) then
!!            dv11dx = A*B*exp(B*x)
!!         else
!!            dv11dx = A*B*exp(-B*x)
!!         endif
!!         dv12dx = -2*D*x*v12
!!         v11=H_Diab(1,1)
!!         v12=H_Diab(1,2)
!!
!!!     Adiabatic energy
!!         T1%ElecStruc%PotEn(2) = sqrt(v11*v11 + v12*v12)
!!         T1%ElecStruc%PotEn(1) = -T1%ElecStruc%PotEn(2)
!!         T1%ElecStruc%PotEn    = T1%ElecStruc%PotEn+gldEShift
!!         T1%ESFlags%ZPotEnCurrent=.true.
!!
!!!     Force
!!         T1%ElecStruc%DerivMat(T1%StateID,T1%StateID,:)=0.d0
!!         if(T1%StateID==1) then
!!            T1%ElecStruc%DerivMat(T1%StateID,T1%StateID,1)=-
!!     $      (dv11dx*v11+dv12dx*v12)/sqrt(v11*v11+v12*v12)
!!         else
!!            T1%ElecStruc%DerivMat(T1%StateID,T1%StateID,1)=
!!     $       (dv11dx*v11+dv12dx*v12)/sqrt(v11*v11+v12*v12)
!!         endif
!!         T1%ESFlags%ZDerivCurrent(T1%StateID,T1%StateID)=.true.
!!!     Non-adiabatic coupling (numerical for now)
!!         Coupling=0.d0
!!         call diagABBC(H_Diab,EVec1,EVec2)
!!         call Tully(x+dx,H_Diab)
!!         call diagABBC(H_Diab,EvecTemp,EDisplace1)
!!         call Tully(x-dx,H_Diab)
!!         call diagABBC(H_Diab,EvecTemp,EDisplace2)
!!         EVecTemp=(EDisplace2-EDisplace1)/(2.d0*dx)
!!         Coupling(1,1)=dot_product(Evec1,EVecTemp)
!!         if(T1%StateID == 2) then
!!            Coupling = - Coupling
!!         endif
!!         T1%ElecStruc%DerivMat(1,2,1)=Coupling(1,1)
!!         T1%ElecStruc%DerivMat(2,1,1)=-Coupling(1,1)
!!
!!         T1%ESFlags%zDerivCurrent(1,2)=.true.
!!         T1%ESFlags%zDerivCurrent(2,1)=.true.

      case (2) ! Persico

         ! Adiabatic energy
         x = T1%Particle(1)%get_pos(1)
         y = T1%Particle(2)%get_pos(1)
         call Persico(x, y, H_Diab)
         v11 = H_diab(1, 1)
         v12 = H_diab(1, 2)
         v22 = H_diab(2, 2)
         T1%ElecStruc%PotEn(1) = 0.5d0 * (v11 + v22 - &
                                          sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
         T1%ElecStruc%PotEn(2) = 0.5d0 * (v11 + v22 + &
                                          sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
         T1%ESFlags%ZPotEnCurrent = .true.

         ! Force
         T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, :) = 0.d0

         dv11dx = KX * (x - X1)
         dv11dy = KY * y

         dv12dx = -2 * alpha * (x - X3) * v12
         dv12dy = gamma * exp(-alpha * (x - X3)**2 - beta * y * y) - 2 * beta * y * v12

         dv22dx = KX * (x - X2)
         dv22dy = KY * Y

         if (T1%StateID == 1) then
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = 0.5d0 * ( &
                                                               dv11dx + dv22dx - &
                                                               (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                                + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 4) = 0.5d0 * ( &
                                                               dv11dy + dv22dy - &
                                                               (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                                + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
         else
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = 0.5d0 * ( &
                                                               dv11dx + dv22dx + &
                                                               (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                                + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 4) = 0.5d0 * ( &
                                                               dv11dy + dv22dy + &
                                                               (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                                + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
         end if
         T1%ESFlags%zDerivCurrent(T1%StateID, T1%StateID) = .true.

         ! Non-adiabatic coupling (numerical)
         Coupling = 0.d0
         call diagABBC(H_Diab, EVec1, EVec2)

         call Persico(x + dx, y, H_Diab)
         call diagABBC(H_Diab, EvecTemp, EDisplace1)
         call Persico(x - dx, y, H_Diab)
         call diagABBC(H_Diab, EvecTemp, EDisplace2)
         EVecTemp = (EDisplace2 - EDisplace1) / (2.d0 * dx)
         Coupling(1, 1) = dot_product(Evec1, EVecTemp)

         call Persico(x, y + dx, H_Diab)
         call diagABBC(H_Diab, EvecTemp, EDisplace1)
         call Persico(x, y - dx, H_Diab)
         call diagABBC(H_Diab, EvecTemp, EDisplace2)
         EVecTemp = (EDisplace2 - EDisplace1) / (2.d0 * dx)
         Coupling(1, 4) = dot_product(Evec1, EVecTemp)

         if (T1%StateID == 2) then
            Coupling = -Coupling
         end if
         T1%ElecStruc%DerivMat(1, 2, :) = Coupling(1, :)
         T1%ElecStruc%DerivMat(2, 1, :) = -Coupling(1, :)

         T1%ESFlags%zDerivCurrent(1, 2) = .true.
         T1%ESFlags%zDerivCurrent(2, 1) = .true.

      case (3) ! Izmaylov

         W1 = izmaylov_params%W1
         W2 = izmaylov_params%W2
         XA = izmaylov_params%XA
         YA = izmaylov_params%YA
         coupC = izmaylov_params%coupC
         MASSX = T1%Particle(1)%Mass
         MASSY = T1%Particle(2)%Mass

         ! Adiabatic energy
         x = T1%Particle(1)%get_pos(1)
         y = T1%Particle(2)%get_pos(1)

         call Izmaylov(x, y, H_Diab)
         v11 = H_diab(1, 1)
         v12 = H_diab(1, 2)
         v22 = H_diab(2, 2)
         T1%ElecStruc%PotEn(1) = 0.5d0 * (v11 + v22 - &
                                          sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
         T1%ElecStruc%PotEn(2) = 0.5d0 * (v11 + v22 + &
                                          sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
         T1%ESFlags%ZPotEnCurrent = .true.

         ! Force
         T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, :) = 0.d0

         dv11dx = W1 * W1 * (x + 0.5d0 * XA)
         dv11dy = W2 * W2 * (y + 0.5d0 * YA)

         dv12dx = 0.d0
         dv12dy = coupC

         dv22dx = W1 * W1 * (x - 0.5d0 * XA)
         dv22dy = W2 * W2 * (y - 0.5d0 * YA)

         if (T1%StateID == 1) then
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = 0.5d0 * ( &
                                                               dv11dx + dv22dx - &
                                                               (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                                + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 4) = 0.5d0 * ( &
                                                               dv11dy + dv22dy - &
                                                               (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                                + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
         else
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = 0.5d0 * ( &
                                                               dv11dx + dv22dx + &
                                                               (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                                + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 4) = 0.5d0 * ( &
                                                               dv11dy + dv22dy + &
                                                               (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                                + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                               / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
         end if
         T1%ESFlags%zDerivCurrent(T1%StateID, T1%StateID) = .true.

         ! Analytical non-adiabatic coupling
         Coupling = 0.d0
         call diagABBC(H_Diab, EVec1, EVec2)

         Coupling(1, 1) = ((v11 - v22) * dv12dx - v12 * (dv11dx - dv22dx)) / &
                          ((v11 - v22)**2 + 4.d0 * v12 * v12)
         Coupling(1, 4) = ((v11 - v22) * dv12dy - v12 * (dv11dy - dv22dy)) / &
                          ((v11 - v22)**2 + 4.d0 * v12 * v12)
         if (T1%StateID == 2) then
            Coupling = -Coupling
         end if

         T1%ElecStruc%DerivMat(1, 2, :) = Coupling(1, :)
         T1%ElecStruc%DerivMat(2, 1, :) = -Coupling(1, :)

         T1%ESFlags%zDerivCurrent(1, 2) = .true.
         T1%ESFlags%zDerivCurrent(2, 1) = .true.

      case (4) ! GAIMS_model

         ! Adiabatic energy
         x = T1%Particle(1)%get_pos(1)

         call GAIMS_model_ham(x, H_Diab, sigma_G)
         v11 = H_diab(1, 1)
         v22 = H_diab(2, 2)
         v12 = H_diab(1, 2)
         v21 = H_diab(2, 1)
         T1%ElecStruc%PotEn(1) = v11
         T1%ElecStruc%PotEn(2) = v22
         T1%ESFlags%ZPotEnCurrent = .true.

         ! Force
         T1%ElecStruc%DerivMat = 0.d0

         if (T1%StateID == 2) then
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = &
               -alpha_2 * a_2 * exp(-alpha_2 * x)
         else
            T1%ElecStruc%DerivMat(T1%StateID, T1%StateID, 1) = &
               -alpha_1 * a_1 * exp(-alpha_1 * x)
         end if
         T1%ESFlags%ZDerivCurrent(T1%StateID, T1%StateID) = .true.
         Coupling = 0.d0
         T1%ESFlags%zDerivCurrent = .true.

         ! SOC,  Eq (11) in https://doi.org/10.1063/1.4707737
         T1%ElecStruc%SOMat = (0.d0, 0.d0)
         T1%ElecStruc%SOMat(1, 2, 2, 1) = conjg(c_1 * sigma_G) ! z*  : S,T-1
         T1%ElecStruc%SOMat(1, 2, 2, 2) = c1i * c_0 * sigma_G ! ib  : S,T0
         T1%ElecStruc%SOMat(1, 2, 2, 3) = c_1 * sigma_G ! z   : S,T1
         T1%ElecStruc%SOMat(2, 1, 1, 2) = c_1 * sigma_G ! z   : T-1,S
         T1%ElecStruc%SOMat(2, 1, 2, 2) = -c1i * c_0 * sigma_G ! -ib : T0,S
         T1%ElecStruc%SOMat(2, 1, 3, 2) = conjg(c_1 * sigma_G) ! z*  : T1,S

         T1%ESFlags%zSOMCurrent = .true.

      case default
         write (fmiOut, *) 'iMethod', gliMethod
         call FMS_DieError( &
            'Error in ToyModel.f; iMethod not recognized.')
      end select

   end subroutine FMS_ToyModel
!?>??ifdef(Debug) then
!?>      call FMS_ElecStrucCalls(T1)
!?>??endif

!>
!! Returns Tully 1 diabatic Hamiltonian
!<
   subroutine Tully(x, sigma)
      real(kind=DefReal), intent(in) :: x
      real(kind=DefReal), intent(out) :: sigma
      real(kind=DefReal) :: rsigma, drsigma

      rsigma = 8.0d0
      drsigma = 2.d0
      if (x <= (rsigma - (drsigma / 2.d0))) then
         sigma = 1.d0
      else if ((x > (rsigma - (drsigma / 2.d0))) .and. &
               (x < (rsigma + (drsigma / 2.d0)))) then
         sigma = 4.d0 * ((x - rsigma) / (drsigma))**3 - &
                 3.d0 * ((x - rsigma) / (drsigma))
      else
         sigma = -1.d0
      end if
!!    subroutine Tully(x,H)
!!    real (kind=DefReal), intent(IN) :: x
!!    real (kind=DefReal), intent(OUT) :: H(2,2)
!!    if(x <= 0) then
!!       H(1,1)    = A*(1-exp(B*x))
!!    else
!!       H(1,1)    = A*(1-exp(-B*x))
!!    endif
!!    H(2,2)=-H(1,1)
!!    H(1,2)=C*exp(-D*x*x)
!!    H(2,1)=H(1,2)
   end subroutine tully

!>
!! Returns Persico diabatic Hamiltonian
!<
   subroutine Persico(x, y, H)
      real(kind=DefReal), intent(in) :: x, y
      real(kind=DefReal), intent(out) :: H(2, 2)
      H(1, 1) = 0.5d0 * KX * (x - X1)**2 + 0.5d0 * KY * y * y
      H(2, 2) = 0.5d0 * KX * (x - X2)**2 + 0.5d0 * KY * y * y + Delta
      H(1, 2) = gamma * y * exp(-alpha * (x - X3)**2) * exp(-beta * y * y)
      H(2, 1) = H(1, 2)
   end subroutine persico

!>
!! Returns Izmaylov diabatic Hamiltonian
!<

   subroutine Izmaylov(x, y, H)
      real(kind=DefReal), intent(in) :: x, y
      real(kind=DefReal), intent(out) :: H(2, 2)
      real(kind=DefReal) :: W1, W2, XA, YA, deltaE, coupC
      W1 = izmaylov_params%W1
      W2 = izmaylov_params%W2
      XA = izmaylov_params%XA
      YA = izmaylov_params%YA
      deltaE = izmaylov_params%deltaE
      coupC = izmaylov_params%coupC

      H(1, 1) = 0.5d0 * (W1**2) * (x + 0.5d0 * XA)**2 + 0.5d0 * (W2**2) * &
                (y + 0.5d0 * YA)**2 + 0.5d0 * deltaE
      H(2, 2) = 0.5d0 * (W1**2) * (x - 0.5d0 * XA)**2 + 0.5d0 * (W2**2) * &
                (y - 0.5d0 * YA)**2 - 0.5d0 * deltaE
      H(1, 2) = coupC * y
      H(2, 1) = H(1, 2)
   end subroutine izmaylov

!>
!! Returns GAIMS_model diabatic Hamiltonian
!<

   subroutine GAIMS_model_ham(x, H, sigma_G)
      real(kind=DefReal), intent(in) :: x
      real(kind=DefReal), intent(out) :: H(2, 2)
      real(kind=DefReal), intent(out) :: sigma_G
      real(kind=DefReal) :: theta_G, gamma_G, r_sigma

      r_sigma = GAIMS_params%r_sigma
      ! Step function sigma(r) controlling SOC sign change, Eq (21)
      if (x <= r_sigma - dr_sigma / 2) then
         sigma_G = 1.d0
      else if ((r_sigma - dr_sigma / 2 < x) .and. &
               (x < r_sigma + dr_sigma / 2)) then
         sigma_G = 4.d0 * ((x - r_sigma) / dr_sigma)**3 &
                   - 3.d0 * (x - r_sigma) / dr_sigma
      else if (x >= r_sigma + dr_sigma / 2) then
         sigma_G = -1.d0
      end if

      gamma_G = sigma_G * sqrt(2.d0 * (real(c_1)**2 + aimag(c_1)**2) + c_0**2)

      ! Heaviside function for theta
      if (x - r_sigma >= 0) then
         theta_G = pi
      else if (x - r_sigma < 0) then
         theta_G = 0
      end if

      ! Eq (15) in https://doi.org/10.1063/1.4707737
      H(1, 1) = a_1 * exp(-alpha_1 * x) + dE
      H(2, 2) = a_2 * exp(-alpha_2 * x)
      H(1, 2) = gamma_G * exp(c1i * theta_G)
      H(2, 1) = gamma_G * exp(-c1i * theta_G)
   end subroutine GAIMS_model_ham

!>
!! Helper routine for numerical derivative couplings; diagonalize a 2x2
!! matrix of the form ( (A,B), (B,C)); the eigenvectors are normalized
!! This is done analytically, since we can.
!<
   subroutine diagABBC(H, EVec1, EVec2)
      real(kind=DefReal), intent(in) :: H(2, 2)
      real(kind=DefReal), intent(out) :: EVec1(2), EVec2(2)
      real(kind=DefReal) :: A, B, C, norm

      if (H(2, 1) /= H(1, 2)) call FMS_DieError( &
         'ToyModel:diagABBC: H_21 /= H12')

      A = H(1, 1)
      B = H(2, 1)
      C = H(2, 2)

      ! Avoid numerical problems if it's already diagonal
      if (B < FPZero) then
         if (A < C) then
            EVec1 = [1.d0, 0.d0]
            EVec2 = [0.d0, 1.d0]
         else
            EVec2 = [1.d0, 0.d0]
            EVec1 = [0.d0, 1.d0]
         end if
         return
      end if

      ! 1st eigenvector
      EVec1(1) = -B
      EVec1(2) = (A - C - sqrt(C * C - 2 * A * C + A * A + 4 * B * B)) / 2.d0
      norm = sqrt(EVec1(1)**2 + EVec1(2)**2)
      EVec1(1) = EVec1(1) / norm
      EVec1(2) = EVec1(2) / norm

      ! 2nd eigenvector
      EVec2(1) = -B
      EVec2(2) = (A - C + sqrt(C * C - 2 * A * C + A * A + 4 * B * B)) / 2.d0
      norm = sqrt(EVec2(1)**2 + EVec2(2)**2)
      EVec2(1) = EVec2(1) / norm
      EVec2(2) = EVec2(2) / norm
   end subroutine diagABBC

end module ToyModelModule
