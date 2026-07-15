! Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Model potentials in the adiabatic basis. Includes:
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

   type :: t_GAIMS_model_params
      !> controls the position of SOC sign change
      real(kind=DefReal) :: r_sigma
   contains
      procedure, public :: initialize => initialize_gaims_model_params
   end type t_GAIMS_model_params
!-----------------------------------------------

   type(t_izmaylov_params) :: izmaylov_params
   type(t_GAIMS_model_params) :: GAIMS_model_params

   public :: FMS_ToyModel, izmaylov_params, GAIMS_model_params

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

   subroutine initialize_gaims_model_params(self, r_sigma)
      class(t_GAIMS_model_params), intent(inout) :: self
      real(kind=DefReal), intent(in) :: r_sigma

      self%r_sigma = r_sigma
   end subroutine initialize_gaims_model_params

   subroutine FMS_ToyModel(T)
      type(T_Trajectory), intent(inout) :: T
      character(len=50) :: errmsg

      select case (gliMethod)

      case (2)

         call persico_model(T)

      case (3)

         call izmaylov_model(T)

      case (4)

         call GAIMS_model(T)

      case default
         write (errmsg, '(a, i0)') 'Invalid iMethod value: ', gliMethod
         call FMS_DieError(errmsg)
      end select

   end subroutine FMS_ToyModel

   subroutine GAIMS_model(T)
      type(T_Trajectory), intent(inout) :: T
      real(kind=DefReal) :: x, H_Diab(2, 2)
      real(kind=DefReal) :: v11, v12, v21, v22
      real(kind=DefReal) :: sigma_G
      integer(kind=DefInt) :: state_id

      x = T%Particle(1)%get_pos(1)
      state_id = T%StateID

      call GAIMS_model_ham(x, H_diab, sigma_G)
      v11 = H_diab(1, 1)
      v22 = H_diab(2, 2)
      v12 = H_diab(1, 2)
      v21 = H_diab(2, 1)

      ! Adiabatic energy
      T%ElecStruc%PotEn(1) = v11
      T%ElecStruc%PotEn(2) = v22
      T%ESFlags%ZPotEnCurrent = .true.

      ! Force
      T%ElecStruc%DerivMat = 0.d0

      if (state_id == 2) then
         T%ElecStruc%DerivMat(state_id, state_id, 1) = -alpha_2 * a_2 * exp(-alpha_2 * x)
      else
         T%ElecStruc%DerivMat(state_id, state_id, 1) = -alpha_1 * a_1 * exp(-alpha_1 * x)
      end if
      T%ESFlags%ZDerivCurrent(state_id, state_id) = .true.

      T%ESFlags%zDerivCurrent = .true.

      ! SOC,  Eq (11) in https://doi.org/10.1063/1.4707737
      T%ElecStruc%SOMat = (0.d0, 0.d0)
      T%ElecStruc%SOMat(1, 2, 2, 1) = conjg(c_1 * sigma_G) ! z*  : S,T-1
      T%ElecStruc%SOMat(1, 2, 2, 2) = c1i * c_0 * sigma_G ! ib  : S,T0
      T%ElecStruc%SOMat(1, 2, 2, 3) = c_1 * sigma_G ! z   : S,T1
      T%ElecStruc%SOMat(2, 1, 1, 2) = c_1 * sigma_G ! z   : T-1,S
      T%ElecStruc%SOMat(2, 1, 2, 2) = -c1i * c_0 * sigma_G ! -ib : T0,S
      T%ElecStruc%SOMat(2, 1, 3, 2) = conjg(c_1 * sigma_G) ! z*  : T1,S

      T%ESFlags%zSOMCurrent = .true.
   end subroutine GAIMS_model

   subroutine persico_model(T)
      type(T_Trajectory), intent(inout) :: T

      real(kind=DefReal) :: Coupling(T%NumStates - 1, T%NumDimensions)
      real(kind=DefReal) :: x, y, H_Diab(2, 2)
      real(kind=DefReal), dimension(2) :: Evec1, Evec2, EDisplace1, EDisplace2, EvecTemp
      real(kind=DefReal) :: v11, v12, v22
      real(kind=DefReal) :: dv11dx, dv12dx, dv22dx
      real(kind=DefReal) :: dv11dy, dv12dy, dv22dy
      integer(kind=DefInt) :: state_id

      ! Parameters for Persico (TODO: These are hardcoded and cannot be specified in Control.dat)
      real(kind=DefReal), parameter :: alpha = 3.d0
      real(kind=DefReal), parameter :: beta = 1.5d0
      real(kind=DefReal), parameter :: gamma = 0.08
      real(kind=DefReal), parameter :: delta = 0.01
      real(kind=DefReal), parameter :: KX = 0.02
      real(kind=DefReal), parameter :: KY = 0.10
      real(kind=DefReal), parameter :: X1 = 4.d0
      real(kind=DefReal), parameter :: X2 = 3.d0
      real(kind=DefReal), parameter :: X3 = 3.d0

      real(kind=DefReal), parameter :: dx = 0.005 ! finite difference step

      state_id = T%StateID

      ! Adiabatic energy
      x = T%Particle(1)%get_pos(1)
      y = T%Particle(2)%get_pos(1)

      call persico_ham(x, y, H_diab)

      v11 = H_diab(1, 1)
      v12 = H_diab(1, 2)
      v22 = H_diab(2, 2)
      T%ElecStruc%PotEn(1) = 0.5d0 * (v11 + v22 - &
                                      sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
      T%ElecStruc%PotEn(2) = 0.5d0 * (v11 + v22 + &
                                      sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
      T%ESFlags%ZPotEnCurrent = .true.

      ! Force
      T%ElecStruc%DerivMat(state_id, state_id, :) = 0.d0

      dv11dx = KX * (x - X1)
      dv11dy = KY * y

      dv12dx = -2 * alpha * (x - X3) * v12
      dv12dy = gamma * exp(-alpha * (x - X3)**2 - beta * y * y) - 2 * beta * y * v12

      dv22dx = KX * (x - X2)
      dv22dy = KY * Y

      if (state_id == 1) then
         T%ElecStruc%DerivMat(state_id, state_id, 1) = 0.5d0 * ( &
                                                       dv11dx + dv22dx - &
                                                       (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                        + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

         T%ElecStruc%DerivMat(state_id, state_id, 4) = 0.5d0 * ( &
                                                       dv11dy + dv22dy - &
                                                       (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                        + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
      else
         T%ElecStruc%DerivMat(state_id, state_id, 1) = 0.5d0 * ( &
                                                       dv11dx + dv22dx + &
                                                       (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                        + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

         T%ElecStruc%DerivMat(state_id, state_id, 4) = 0.5d0 * ( &
                                                       dv11dy + dv22dy + &
                                                       (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                        + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
      end if
      T%ESFlags%zDerivCurrent(state_id, state_id) = .true.

      ! Non-adiabatic coupling (numerical)
      Coupling = 0.d0
      call diagABBC(H_Diab, EVec1, EVec2)

      call persico_ham(x + dx, y, H_Diab)
      call diagABBC(H_Diab, EvecTemp, EDisplace1)
      call persico_ham(x - dx, y, H_Diab)
      call diagABBC(H_Diab, EvecTemp, EDisplace2)
      EVecTemp = (EDisplace2 - EDisplace1) / (2.d0 * dx)
      Coupling(1, 1) = dot_product(Evec1, EVecTemp)

      call persico_ham(x, y + dx, H_Diab)
      call diagABBC(H_Diab, EvecTemp, EDisplace1)
      call persico_ham(x, y - dx, H_Diab)
      call diagABBC(H_Diab, EvecTemp, EDisplace2)
      EVecTemp = (EDisplace2 - EDisplace1) / (2.d0 * dx)
      Coupling(1, 4) = dot_product(Evec1, EVecTemp)

      if (state_id == 2) then
         Coupling = -Coupling
      end if
      T%ElecStruc%DerivMat(1, 2, :) = Coupling(1, :)
      T%ElecStruc%DerivMat(2, 1, :) = -Coupling(1, :)

      T%ESFlags%zDerivCurrent(1, 2) = .true.
      T%ESFlags%zDerivCurrent(2, 1) = .true.

   contains
      !>
      !! Returns Persico diabatic Hamiltonian
      !<
      subroutine persico_ham(x, y, H)
         real(kind=DefReal), intent(in) :: x, y
         real(kind=DefReal), intent(out) :: H(2, 2)

         H(1, 1) = 0.5d0 * KX * (x - X1)**2 + 0.5d0 * KY * y * y
         H(2, 2) = 0.5d0 * KX * (x - X2)**2 + 0.5d0 * KY * y * y + Delta
         H(1, 2) = gamma * y * exp(-alpha * (x - X3)**2) * exp(-beta * y * y)
         H(2, 1) = H(1, 2)
      end subroutine persico_ham

   end subroutine persico_model

   subroutine izmaylov_model(T)
      type(T_Trajectory), intent(inout) :: T

      real(kind=DefReal) :: Coupling(T%NumStates - 1, T%NumDimensions)
      real(kind=DefReal), dimension(2) :: Evec1, Evec2
      real(kind=DefReal) :: x, y, H_Diab(2, 2)
      real(kind=DefReal) :: v11, v12, v22
      real(kind=DefReal) :: dv11dx, dv12dx, dv22dx
      real(kind=DefReal) :: dv11dy, dv12dy, dv22dy
      real(kind=DefReal) :: massx, massy
      real(kind=DefReal) :: W1, W2, XA, YA, coupC
      integer(kind=DefInt) :: state_id

      W1 = izmaylov_params%W1
      W2 = izmaylov_params%W2
      XA = izmaylov_params%XA
      YA = izmaylov_params%YA
      coupC = izmaylov_params%coupC
      MASSX = T%Particle(1)%Mass
      MASSY = T%Particle(2)%Mass
      state_id = T%StateID

      ! Adiabatic energy
      x = T%Particle(1)%get_pos(1)
      y = T%Particle(2)%get_pos(1)

      call izmaylov_ham(x, y, H_Diab)
      v11 = H_diab(1, 1)
      v12 = H_diab(1, 2)
      v22 = H_diab(2, 2)
      T%ElecStruc%PotEn(1) = 0.5d0 * (v11 + v22 - &
                                      sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
      T%ElecStruc%PotEn(2) = 0.5d0 * (v11 + v22 + &
                                      sqrt(v11 * v11 - 2.d0 * v22 * v11 + v22 * v22 + 4.d0 * v12 * v12))
      T%ESFlags%ZPotEnCurrent = .true.

      ! Force
      T%ElecStruc%DerivMat(state_id, state_id, :) = 0.d0

      dv11dx = W1 * W1 * (x + 0.5d0 * XA)
      dv11dy = W2 * W2 * (y + 0.5d0 * YA)

      dv12dx = 0.d0
      dv12dy = coupC

      dv22dx = W1 * W1 * (x - 0.5d0 * XA)
      dv22dy = W2 * W2 * (y - 0.5d0 * YA)

      if (state_id == 1) then
         T%ElecStruc%DerivMat(state_id, state_id, 1) = 0.5d0 * ( &
                                                       dv11dx + dv22dx - &
                                                       (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                        + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

         T%ElecStruc%DerivMat(state_id, state_id, 4) = 0.5d0 * ( &
                                                       dv11dy + dv22dy - &
                                                       (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                        + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
      else
         T%ElecStruc%DerivMat(state_id, state_id, 1) = 0.5d0 * ( &
                                                       dv11dx + dv22dx + &
                                                       (v22 * dv22dx - v11 * dv22dx - v22 * dv11dx &
                                                        + v11 * dv11dx + 4.d0 * dv12dx * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))

         T%ElecStruc%DerivMat(state_id, state_id, 4) = 0.5d0 * ( &
                                                       dv11dy + dv22dy + &
                                                       (v22 * dv22dy - v11 * dv22dy - v22 * dv11dy &
                                                        + v11 * dv11dy + 4.d0 * dv12dy * v12) &
                                                       / sqrt(v22 * v22 - 2.d0 * v22 * v11 + v11 * v11 + 4.d0 * v12 * v12))
      end if
      T%ESFlags%zDerivCurrent(state_id, state_id) = .true.

      ! Analytical non-adiabatic coupling
      Coupling = 0.d0
      call diagABBC(H_Diab, EVec1, EVec2)

      Coupling(1, 1) = ((v11 - v22) * dv12dx - v12 * (dv11dx - dv22dx)) / &
                       ((v11 - v22)**2 + 4.d0 * v12 * v12)
      Coupling(1, 4) = ((v11 - v22) * dv12dy - v12 * (dv11dy - dv22dy)) / &
                       ((v11 - v22)**2 + 4.d0 * v12 * v12)
      if (state_id == 2) then
         Coupling = -Coupling
      end if

      T%ElecStruc%DerivMat(1, 2, :) = Coupling(1, :)
      T%ElecStruc%DerivMat(2, 1, :) = -Coupling(1, :)

      T%ESFlags%zDerivCurrent(1, 2) = .true.
      T%ESFlags%zDerivCurrent(2, 1) = .true.
   end subroutine izmaylov_model

!>
!! Returns Izmaylov diabatic Hamiltonian
!<
   subroutine izmaylov_ham(x, y, H)
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
   end subroutine izmaylov_ham

!>
!! Returns GAIMS_model diabatic Hamiltonian
!<
   subroutine GAIMS_model_ham(x, H, sigma_G)
      real(kind=DefReal), intent(in) :: x
      real(kind=DefReal), intent(out) :: H(2, 2)
      real(kind=DefReal), intent(out) :: sigma_G
      real(kind=DefReal) :: theta_G, gamma_G, r_sigma

      r_sigma = GAIMS_model_params%r_sigma
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
