! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module OverlapModule
   use GlobalModule
   use ParticleModule
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_Forces, FMS_PotentialT, FMS_KineticT, &
                                    FMS_Coupling, FMS_SOCoupling, FMS_Dipole, FMS_PhaseDot, FMS_TransDipoleIJxf, &
                                    FMS_SetIgnoreError
   implicit none

   private
   public :: overlap, overlap_KE, overlap_V, overlap_S_dot
   public :: overlap_dx_M_trajectory

   interface overlap
      module procedure overlap_trajectory
   end interface

   interface overlap_KE
      module procedure overlap_KE_trajectory
   end interface

   interface overlap_V
      module procedure overlap_V_trajectory
   end interface

   interface overlap_S_dot
      module procedure overlap_S_dot_trajectory
   end interface

contains

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!              Complex  Gaussians

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_CG(x_i, p_i, alpha_i, &
                       x_j, p_j, alpha_j) result(S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the overlap between two complex gaussian
      real(kind=DefReal), intent(in) :: x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j

      complex(kind=DefComp) :: S_ij

      real(kind=DefReal) :: delta_x, delta_p, & ! difference is position & momentum
                            x_cent, & ! the centroid position
                            prefactor, &
                            real_part, &
                            imag_part

      delta_x = x_i - x_j
      delta_p = p_i - p_j

      real_part = (alpha_i * alpha_j * delta_x**2 + 0.25d0 * delta_p**2) / &
         (alpha_i + alpha_j)

      if (real_part <= 10.0d0) then
         prefactor = sqrt(2.0 * sqrt(alpha_i * alpha_j) / &
                          (alpha_i + alpha_j))

         x_cent = (alpha_i * x_i + alpha_j * x_j) / (alpha_i + alpha_j)

         imag_part = (p_i * x_i - p_j * x_j) - x_cent * delta_p

         S_ij = prefactor * exp(-real_part + (0., 1.) * imag_part)
      else
         S_ij = (0., 0.)
      end if

   end function overlap_CG

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dp_CG(x_i, p_i, alpha_i, &
                          x_j, p_j, alpha_j, &
                          prefactor) result(dp_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expection of d/dp between two complex gaussians
! right acting
      real(kind=DefReal), intent(in) :: x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j
      logical, optional :: prefactor ! pre calculated overlap
      complex(kind=DefComp) :: dp_S_ij

      real(kind=DefReal) :: delta_x, delta_p

      delta_x = x_i - x_j
      delta_p = p_i - p_j

      dp_S_ij = (delta_p + 2.d0 * (0., 1.) * alpha_i * delta_x) &
                / (2.d0 * (alpha_i + alpha_j))

      if (present(prefactor)) then
         if (.not. prefactor) then
            dp_S_ij = dp_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                           x_j, p_j, alpha_j)
         end if
      else
         dp_S_ij = dp_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j)
      end if

   end function overlap_dp_CG

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dx_CG(x_i, p_i, alpha_i, &
                          x_j, p_j, alpha_j, &
                          prefactor) result(dx_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expection of d/dx between two complex gaussians
      real(kind=DefReal), intent(in) :: x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j
      logical, optional :: prefactor
      complex(kind=DefComp) :: dx_S_ij

      real(kind=DefReal) :: P_ij, delta_x

      delta_x = x_i - x_j
      P_ij = alpha_i * p_j + alpha_j * p_i

      dx_S_ij = ((0., 1.) * P_ij - 2.d0 * alpha_j * alpha_i * delta_x) &
                / (alpha_i + alpha_j)

      if (present(prefactor)) then
         if (.not. prefactor) then
            dx_S_ij = dx_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                           x_j, p_j, alpha_j)
         end if
      else
         dx_S_ij = dx_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j)
      end if

   end function overlap_dx_CG

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_d2x_CG(x_i, p_i, alpha_i, &
                           x_j, p_j, alpha_j, &
                           prefactor) result(d2x_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expectation of the d^2/dx^2  between two complex gaussians
      real(kind=DefReal), intent(in) :: x_i, p_i, alpha_i, &
                                        x_j, p_j, alpha_j
      logical, optional :: prefactor ! pre calculated overlap
      complex(kind=DefComp) :: d2x_S_ij
      real(kind=DefReal) :: P_ij, delta_x

      delta_x = x_i - x_j
      P_ij = alpha_i * p_j + alpha_j * p_i

      d2x_S_ij = -(+4.d0 * (0., 1) * alpha_i * alpha_j * delta_x * P_ij &
                   + 2.d0 * alpha_i * alpha_j * (alpha_i + alpha_j) &
                   - 4.d0 * delta_x**2 * alpha_i**2 * alpha_j**2 &
                   + P_ij**2) / (alpha_i + alpha_j)**2

      if (present(prefactor)) then
         if (.not. prefactor) then
            d2x_S_ij = d2x_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                             x_j, p_j, alpha_j)
         end if
      else
         d2x_S_ij = d2x_S_ij * overlap_CG(x_i, p_i, alpha_i, &
                                          x_j, p_j, alpha_j)
      end if

   end function overlap_d2x_CG

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!               Particles

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_particle(P1, P2) result(S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calculates the overlap between two particles
! optionally, a single component can be calculated.  This is a legacy feature
! will be phased out
      type(T_Particle), intent(in) :: P1, P2
      complex(kind=DefComp) :: S_ij

      integer(kind=DefInt) :: ndim, i

      real(kind=DefReal) :: x_1, p_1, a_1, &
                            x_2, p_2, a_2

      ndim = P1%NumDimensions

      a_1 = P1%width
      a_2 = P2%width

      S_ij = (1., 0.)

      do i = 1, ndim
         x_1 = P1%get_pos(i)
         x_2 = P2%get_pos(i)

         p_1 = P1%get_mom(i)
         p_2 = P2%get_mom(i)

         S_ij = S_ij * overlap_CG(x_1, p_1, a_1, &
                                  x_2, p_2, a_2)
      end do

   end function overlap_particle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dx_particle(P1, P2, prefactor) result(dx_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! this calculates the expectation of the derivative operator for a particle
! returns a ndim vector
      type(T_Particle), intent(in) :: P1, P2
      logical, optional :: prefactor
      complex(kind=DefComp) :: dx_S_ij(P1%NumDimensions)

      integer(kind=DefInt) :: ndim_1, ndim_2, n

      real(kind=DefReal) :: x_1, p_1, a_1, &
                            x_2, p_2, a_2

      ndim_1 = P1%NumDimensions
      ndim_2 = P2%NumDimensions

      if (ndim_1 /= ndim_2) &
         call FMS_DieError('overlap_dx_particle: dimension mismatch')

      a_1 = P1%width
      a_2 = P2%width

      do n = 1, ndim_1
         x_1 = P1%get_pos(n)
         x_2 = P2%get_pos(n)

         p_1 = P1%get_mom(n)
         p_2 = P2%get_mom(n)

         dx_S_ij(n) = overlap_dx_CG(x_1, p_1, a_1, &
                                    x_2, p_2, a_2, prefactor=.true.)
      end do

      if (present(prefactor)) then
         if (.not. prefactor) then
            dx_S_ij = dx_S_ij * overlap_particle(P1, P2)
         end if
      else
         dx_S_ij = dx_S_ij * overlap_particle(P1, P2)
      end if

   end function overlap_dx_particle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dp_particle(P1, P2, prefactor) result(dp_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calculates the expectation of d/dp.  Needed for calculating Sdot
! if prefactor is true the overlap part of the expression is not calculated
      type(T_Particle), intent(in) :: P1, P2
      logical, optional :: prefactor
      complex(kind=DefComp) :: dp_S_ij(P1%NumDimensions)

      integer(kind=DefInt) :: ndim_1, ndim_2, n

      real(kind=DefReal) :: x_1, p_1, a_1, &
                            x_2, p_2, a_2

      ndim_1 = P1%NumDimensions
      ndim_2 = P2%NumDimensions

      if (ndim_1 /= ndim_2) &
         call FMS_DieError('overlap_dp_particle: dimension mismatch')

      a_1 = P1%width
      a_2 = P2%width

      do n = 1, ndim_1
         x_1 = P1%get_pos(n)
         x_2 = P2%get_pos(n)

         p_1 = P1%get_mom(n)
         p_2 = P2%get_mom(n)

         dp_S_ij(n) = overlap_dp_CG(x_1, p_1, a_1, &
                                    x_2, p_2, a_2, prefactor=.true.)

      end do

      if (present(prefactor)) then
         if (.not. prefactor) then
            dp_S_ij = dp_S_ij * overlap_particle(P1, P2)
         end if
      else
         dp_S_ij = dp_S_ij * overlap_particle(P1, P2)
      end if

   end function overlap_dp_particle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_d2x_particle(P1, P2, prefactor) result(d2x_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calculates the Del^2 operator between two particles
      type(T_Particle), intent(in) :: P1, P2
      logical, optional :: prefactor
      complex(kind=DefComp) :: d2x_S_ij

      integer(kind=DefInt) :: ndim_1, ndim_2, i

      real(kind=DefReal) :: x_1, p_1, alpha_1, &
                            x_2, p_2, alpha_2

      ndim_1 = P1%NumDimensions
      ndim_2 = P2%NumDimensions

      if (ndim_1 /= ndim_2) &
         call FMS_DieError('overlap_d2x_particle: Particle Dimensionalities must be the same!')

      alpha_1 = P1%width
      alpha_2 = P2%width

      d2x_S_ij = (0., 0.)
      do i = 1, ndim_1
         x_1 = P1%get_pos(i)
         x_2 = P2%get_pos(i)

         p_1 = P1%get_mom(i)
         p_2 = P2%get_mom(i)

         d2x_S_ij = d2x_S_ij + overlap_d2x_CG(x_1, p_1, alpha_1, &
                                              x_2, p_2, alpha_2, &
                                              prefactor=.true.)
      end do

      if (present(prefactor)) then
         if (.not. prefactor) then
            d2x_S_ij = d2x_S_ij * overlap_particle(P1, P2)
         end if
      else
         d2x_S_ij = d2x_S_ij * overlap_particle(P1, P2)
      end if

   end function overlap_d2x_particle

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!                     Trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_trajectory(T1, T2, same_state) result(S)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      logical, optional :: same_state
      complex(kind=DefComp) :: S
      real(kind=DefReal) :: time_tmp1, time_tmp2
      integer :: n

      call cpu_time(time_tmp1)
! In some places of the code only the overlap between
! trajectories on the same state is neeeded.
      if (present(same_state)) then
         if (same_state .and. (T1%stateID /= T2%stateID)) then
            S = (0., 0.)
            return !note- this won't update olaptime, but that's ok
         end if
      end if

! GAIMS added
      if (present(same_state)) then
         if (same_state .and. (T1%Ms /= T2%Ms)) then
            S = (0., 0.)
            return !note- this won't update olaptime, but that's ok
         end if
      end if
! GAIMS added end

      S = exp((0., 1.) * (T2%Phase - T1%Phase))
      do n = 1, T1%NumParticles
         S = S * overlap_particle(T1%Particle(n), T2%Particle(n))
      end do

      call cpu_time(time_tmp2)
      olaptime = olaptime + time_tmp2 - time_tmp1
   end function overlap_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dp_trajectory(T1, T2, S_ij_precalc) result(dp_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: dp_S_ij(T1%NumDimensions)

      complex(kind=DefComp) :: dp_S_ij_tmp(T1%Particle(1)%NumDimensions, &
                                           T1%NumParticles)

      complex(kind=DefComp) :: S_ij

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%NumDimensions

      if (present(S_ij_precalc)) then
         S_ij = S_ij_precalc
      else
         S_ij = overlap_trajectory(T1, T2)
      end if

      do n = 1, npart
         dp_S_ij_tmp(:, n) = overlap_dp_particle(T1%Particle(n), T2%Particle(n), &
                                                 prefactor=.true.)
      end do

      dp_S_ij = reshape(dp_S_ij_tmp, [ndim]) * S_ij

   end function overlap_dp_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dx_trajectory(T1, T2, S_ij_precalc) result(dx_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: dx_S_ij(T1%NumDimensions)

      complex(kind=DefComp) :: dx_S_ij_tmp(T1%Particle(1)%NumDimensions, &
                                           T1%NumParticles)

      complex(kind=DefComp) :: S_ij

      integer(kind=DefInt) :: npart, ndim, n

      npart = T1%NumParticles
      ndim = T1%NumDimensions

      if (present(S_ij_precalc)) then
         S_ij = S_ij_precalc
      else
         S_ij = overlap_trajectory(T1, T2)
      end if

      do n = 1, npart
         dx_S_ij_tmp(:, n) = overlap_dx_particle(T1%Particle(n), T2%Particle(n), prefactor=.true.)
      end do

      dx_S_ij = reshape(dx_S_ij_tmp, [ndim]) * S_ij

   end function overlap_dx_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_dx_M_trajectory(T1, T2, S_ij_precalc) result(dx_M_S_ij)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! this is the expectation value of the momentum operator over mass
! this appears in the equations of motion on the off diagonal coupling
! different states together theough the NACME
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: dx_M_S_ij(T1%NumDimensions)

      complex(kind=DefComp) :: dx_M_S_ij_tmp(T1%Particle(1)%NumDimensions, &
                                             T1%NumParticles)

      complex(kind=DefComp) :: S_ij

      integer(kind=DefInt) :: npart, n, ndim

      npart = T1%NumParticles
      ndim = T1%NumDimensions

      if (present(S_ij_precalc)) then
         S_ij = S_ij_precalc
      else
         S_ij = overlap_trajectory(T1, T2)
      end if

      do n = 1, npart
         dx_M_S_ij_tmp(:, n) = overlap_dx_particle(T1%Particle(n), T2%Particle(n), prefactor=.true.) &
                               / T1%get_mass(n)
      end do

      dx_M_S_ij = reshape(dx_M_S_ij_tmp, [ndim]) * S_ij

   end function overlap_dx_M_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_V_trajectory(T_i, T_j, T_c, S_ij_precalc) result(V)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calculates expectation of the potential energy operator between two trajectories
! Updated for GAIMS
!
! This is the centroid:
!  X_c = 0.5 * ( X_i^I + X_j^J )
!
! In the spin-diabatic scheme there are 10 cases of possible interactions
! which can be grouped into singlet-singlet (1) singlet-triplet (2) and
! triplet-triplet (3)
! terms.  The indices used in the following scheme are:
! i,j       ... TrajectoryID
! I,J       ... StateID
! Msi,Msj   ... Spin-Multiplicity
! CBFi,CBFj ... Contracted basis function ID
!
! ###########################################################################
!
! 1. Singlet-Singlet
!
! A. energy of a trajectory (i=i, I=I)
!     < X_i^I | V | X_i^I > = V(X_i)
!
! B. energy between two distinct trajectories ( i!=j. I=J )
!  < X_i^I | V | X_j^I > = V(X_c) < X_i^I |X_j^I >
!
! C. coupling between trajectories on different states ( i!=j, I!=J )
!  < X_i^I | V | X_j^J > = F_x^{I,J}(X_c) .x. < X_i^I | P/M | X_j^I >
!
! ###########################################################################
!
! 2. Singlet-Triplet
!
! D. Coupling between trajectories for all three sublevels (Ms=-1,0,1)
!   < X_i^I | V | X_j^J > = SO_{I,J}^{Msi,Msj}(X_c) .x. < X_i^I |X_j^I >
!
! ###########################################################################
!
! 3. Triplet-Triplet
!
! E. energy of a trajectory (i=j, I=J, Msi=Msj) As this case is identical to
!    case(A), we use the same call
!   < X_i^I | V | X_i^I > = V(X_i)
!
! F. energy between two distinct trajectories of the same multiplicity
!  ( i!=j, I=J, CBFi!=CBFj, Msi=Msj)
!   < X_i^I | V | X_j^I > = V(X_c) < X_i^I |X_j^I >
!
! G. coupling between two distinct trajectories of different multiplicity
!   (i!=j, I=J, CBFi!=CBFj, Msi!=Msj)
!   < X_i^I | V | X_j^I > = SO_{I,I}^{Msi,Msj}(X_c) .x. < X_i^I |X_j^I >
!
! H. coupling between two trajectories on different states and same multiplicity
! (i!=j, I!=J, CBFi!=CBFj, Msi=Msj)
!  < X_i^I | V | X_j^J > = F_x^{I,J}(X_c) .x. < X_i^I | P/M | X_j^J >
!
! I. coupling between two trajectories on different states and different
! multiplicity
!  (i!=j, I!=J, CBFi!=CBFj, Msi!=Msj)
!  < X_i^I | V | X_j^J > = SO_{I,J}^{Msi,Msj}(X_c) .x. < X_i^I |  X_j^J >
!
!
!      Singlet-Singlet         Singlet-Triplet               Triplet-Triplet
!           i=j                                                      i=j
!          /   \                                                    /   \
!      yes/     \no                                             yes/     \no
!        /       \                                                /       \
!      case(A)  is=js                                         case(E)
!      __I=J_________
!               /   \               case(D)                           |
!               |
!           yes/     \no                                             yes
!           no
!             /       \                                               |
!             |
!       case(B)  case(C)                                           MsI=MsJ
!       MsI=MsJ
!                                                                   /  \
!                                                                   /  \
!                                                               yes/    \no
!                                                               yes/    \no
!                                                                 /      \
!                                                                 /      \
!                                                             case(F)  case(G)
!                                                             case(H)  case(I)
!
! ###########################################################################
!
! Open Questions for this procedure
!
! -) Symmetry of SOCs (SO{I,J} = SO{J,I})               ?????
! -) case(E and case (g)): Neglect diagonal SOC         ?????
! -) case(I): Is coupling NACV or NACV+SOC              ?????
!
! ###########################################################################

      type(T_Trajectory), intent(in) :: T_i
      type(T_Trajectory), intent(in), optional :: T_j, T_c
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: V
! enumerated type for cases
      integer :: mode
      integer, parameter :: &
         A = 1, &
         B = 2, &
         C = 3, &
         D = 4, &
         E = 5, &
         F = 6, &
         G = 7, &
         H = 8, &
         I = 9

      integer(kind=DefInt) :: it, jt, & ! label of trajectories
                              ic, jc, & ! centroid's labels
                              is, js ! state labels

      logical :: cent_match

!bfec
      real(kind=DefReal) :: CVec1(T_i%NumDimensions), CVec2(T_i%NumDimensions)

      real(kind=DefReal) :: F_ij(T_i%NumDimensions), time_tmp1, time_tmp2

      complex(kind=DefComp) :: S_ij ! overlap between trajectories

! xf added
      integer(kind=DefInt) :: bmis
      real(kind=DefReal) :: dipole(3), TDipxf(4), ft_xf(3)
      complex(kind=DefComp) :: dip_nuc_xf(3)
! xf end added

! GAIMS added
      integer(kind=DefInt) :: Msi, Msj, & ! multiplicities
                              CBFi, CBFj ! CBF label
      complex(kind=DefComp) :: SOC_ij ! SOCoupling
! GAIMS end added

      call cpu_time(time_tmp1)
! work out what we are doing
      if (.not. present(T_j)) then
         mode = A ! case(A) is equal to case(E), so only one of them is called
      else
! GAIMS changed
         it = T_i%TrajID; jt = T_j%TrajID
         Msi = T_i%Ms; Msj = T_j%Ms
         CBFi = T_i%CBF; CBFj = T_j%CBF
! GAIMS end changed
!bfec
         if (glzCentroids .and. (CBFi /= CBFj)) then
            ic = T_c%CentID(1); jc = T_c%CentID(2)
         end if
         is = T_i%StateID; js = T_j%StateID

         if ((T_i%triplet .eqv. .false.) .and. (T_j%triplet .eqv. .false.)) then ! (1)
            if (is == js) then
               mode = B
            else
               mode = C
            end if
         elseif (.not. (T_i%triplet .and. T_j%triplet)) then ! (2)
            mode = D
         else ! (3)
!    if( CBFi /= CBFj ) then
            if (is == js) then
               if (CBFi == CBFj) then
                  mode = F
               else
                  mode = G
               end if
            else
               if (Msi == Msj) then
                  mode = H
               else
                  mode = I
               end if
            end if
         end if
!  endif
      end if

!write(fmiOut,*) "Mode =", mode
!write(fmiOut,*) "Trajectory ",it, " on state ",is," coupled to"
!write(fmiOut,*) "Trajectory ",jt, " on state ",js," via mode",mode

! check the centroid
      select case (mode)
      case (B, C, D, G, H, I)
!bfec
         if (glzCentroids) then
            if (.not. present(T_c)) &
               call FMS_DieError('overlap_V_trajectory :: no centroid passed')

            ! check that the correct centroid was passed
!    cent_match = (it==ic .and. jt==jc) .or. (jt==ic .and. it==jc)
            cent_match = (CBFi == ic .and. CBFj == jc) .or. (CBFj == ic .and. CBFi == jc)
            write (fmiout, *) "(CBFij == ijc .and. CBFji == jic)", "CBFi, ic, CBFj, jc", CBFi, ic, CBFj, jc
            flush (fmiout)
!    cent_match = (i==ic .and. j==jc) .or. (j==ic .and. i==jc)
            if (.not. cent_match) then
               call FMS_DieError('overlap_V_trajectory :: centroids dont match')
            end if
         end if

!WJG - replaced below with UpdateCentroid within PropQCVV
!   call FMS_InitCentroid( T_i, T_j, T_c )

         ! See if the overlap passed
         if (present(S_ij_precalc)) then
            S_ij = S_ij_precalc
         else
            S_ij = overlap_trajectory(T_i, T_j)
         end if

         ! if the overlap is small, set V to 0
         if (abs(S_ij) < gldOLapThresh) then
            V = 0.d0
!bfec
            if (glzCentroids) then
               call FMS_SetIgnoreError(T_c, .true.)
            end if
            return
         end if

      case default
         continue
      end select

      select case (mode)
      case (A)
         V = FMS_PotentialT(T_i)
!   write(fmiOut,*) '#########AAAAAAAAAAA###########'
!   write(fmiOut,*) T_i%TrajID,T_i%StateID
!   write(fmiOut,*) FMS_PotentialT( T_i )
!   write(fmiOut,*) '###############################'

         if (glzxfaims) then
            bmis = T_i%StateID
            dipole(:) = FMS_Dipole(T_i, bmis)
            dip_nuc_xf(:) = dipole_nucl_xf(T_i, T_i)
            ft_xf(:) = field_xf(T_i)

            V = V &
                - (dipole(1) + dip_nuc_xf(1)) * ft_xf(1) &
                - (dipole(2) + dip_nuc_xf(2)) * ft_xf(2) &
                - (dipole(3) + dip_nuc_xf(3)) * ft_xf(3)
         end if

      case (B)
!bfec
         if (glzCentroids) then
            V = S_ij * FMS_PotentialT(T_c)
         else
            V = S_ij * 0.5d0 * (FMS_PotentialT(T_i) + FMS_PotentialT(T_j, is))
            write (fmiout, *) 'Using LI approximation for matrix elements.'
         end if

         if (glzxfaims) then
            bmis = T_i%StateID
            dipole(:) = FMS_Dipole(T_i, bmis)
            dip_nuc_xf(:) = dipole_nucl_xf(T_i, T_j)
            ft_xf(:) = field_xf(T_i)

            V = V &
                - (dipole(1) + dip_nuc_xf(1)) * ft_xf(1) * S_ij &
                - (dipole(2) + dip_nuc_xf(2)) * ft_xf(2) * S_ij &
                - (dipole(3) + dip_nuc_xf(3)) * ft_xf(3) * S_ij
         end if

!   write(fmiOut,*) '#########BBBBBBBBBBB###########'
!   write(fmiOut,*) T_i%TrajID,T_i%StateID
!   write(fmiOut,*) T_j%TrajID,T_j%StateID
!   write(fmiOut,*) FMS_PotentialT( T_C )
!   write(fmiOut,*) '###############################'

      case (C)
!bfec
         if (glzCentroids) then
            F_ij = FMS_Coupling(T_c, is, js)
            V = dot_product(F_ij, overlap_dx_M_trajectory(T_i, T_j, S_ij))
         else
            CVec1 = FMS_Coupling(T_i, is, js)
            CVec2 = FMS_Coupling(T_j, is, js)
! phase should be ok, need to be checked for FMS/TC with CASSCF.
            F_ij = 0.5d0 * (CVec1 + CVec2)
            write (fmiout, *) 'Using LI approximation for matrix elements.'
            V = dot_product(F_ij, overlap_dx_M_trajectory(T_i, T_j, S_ij))
         end if

         if (glzxfaims) then
            TDipxf(:) = FMS_TransDipoleIJxf(T_c)
            ft_xf(:) = field_xf(T_i)

            V = V &
                - TDipxf(1) * ft_xf(1) * S_ij &
                - TDipxf(2) * ft_xf(2) * S_ij &
                - TDipxf(3) * ft_xf(3) * S_ij

1987        format(a6, f7.3, a11, i2, i2, a11, i2, i2, a16, f7.3, f7.3, a12, f8.4, f8.4, f8.4)
            write (fmiOut, 1987) 'Time', T_i%get_time(), ' Traj i,j:', T_i%TrajID, T_j%TrajID &
               , ' State i,j:', T_i%StateID, T_j%StateID, ' Population i,j:', T_i%Pop, T_j%Pop, ' Field E(t):' &
               , ft_xf(1), ft_xf(2), ft_xf(3)

         end if

      case (D)
         SOC_ij = FMS_SOCoupling(T_c, is, js, Msi, Msj)
         V = S_ij * SOC_ij
!   write(fmiOut,*) '########DDDD#############'
!   write(fmiOut,*) 'd',T_i%TrajID,T_j%TrajID,T_i%StateID,T_j%StateID
!   write(fmiOut,*) Msi,Msj
!   write(fmiOut,*) 'SOCc',SOC_ij
!   write(fmiOut,*) S_ij,V
!   write(fmiOut,*) FMS_PotentialT( T_i ),FMS_PotentialT( T_j, is )
!   write(fmiOut,*) '#########################'

      case (F)
         V = (0.d0, 0.d0)

      case (G)
         if (Msi /= Msj) then
            V = (0.d0, 0.d0)
         else
            V = S_ij * FMS_PotentialT(T_c)
         end if
!   write(fmiOut,*) '#########GGGG############'
!   write(fmiOut,*) 'f',T_i%TrajID,T_j%TrajID,T_i%StateID,T_j%StateID
!   write(fmiOut,*) FMS_PotentialT( T_c )
!   write(fmiOut,*) FMS_PotentialT( T_i ),FMS_PotentialT( T_j, is )
!   write(fmiOut,*) S_ij,V
!   write(fmiOut,*) '#########################'

      case (H)
         if (CBFi == CBFj) then
            V = (0.d0, 0.d0)
         else
            F_ij = FMS_Coupling(T_c, is, js)
            ! Careful here, dot_product( A,B ) = sum( conjg(A)*B )
            V = dot_product(F_ij, overlap_dx_M_trajectory(T_i, T_j, S_ij))
         end if
!   write(fmiOut,*) '##########HHHH###########'
!   write(fmiOut,*) 'i',T_i%TrajID,T_j%TrajID,T_i%StateID,T_j%StateID
!   write(fmiOut,*) F_ij,V
!   write(fmiOut,*) '#########################'

      case (I)
         if (abs(Msi - Msj) > 1) then
            V = (0.d0, 0.d0)
         elseif (CBFi == CBFj) then
            V = (0.d0, 0.d0)
         else
            SOC_ij = FMS_SOCoupling(T_c, is, js, Msi, Msj)
            V = S_ij * SOC_ij
         end if
!   write(fmiOut,*) '#########IIII############'
!   write(fmiOut,*) 'j',T_i%TrajID,T_j%TrajID,T_i%StateID,T_j%StateID
!   write(fmiOut,*) Msi,Msj,CBFi,CBFj
!   write(fmiOut,*) 'SOCc',SOC_ij
!   write(fmiOut,*) S_ij,V
!   write(fmiOut,*) '#########################'

      case default
         call FMS_DieError("Invalid mode in overlap_V_trajectory")

      end select

      call cpu_time(time_tmp2)
      pottime = pottime + time_tmp2 - time_tmp1
   end function overlap_V_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_KE_trajectory(T1, T2, S_ij_precalc) result(KE)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: KE

      complex(kind=DefComp) :: S_ij

      integer :: npart, n

      KE = (0., 0.)
! the KE operator is diagonal in the electronic state basis
! GAIMS updated
      if ((T1%StateID == T2%StateID) .and. (T1%Ms == T2%Ms)) then
! GAIMS updated end
!else
         if (present(S_ij_precalc)) then
            S_ij = S_ij_precalc
         else
            S_ij = overlap_trajectory(T1, T2)
         end if

         npart = T1%NumParticles
         do n = 1, npart
            KE = KE - overlap_d2x_particle(T1%Particle(n), T2%Particle(n), prefactor=.true.) &
                 / (2.d0 * T1%get_mass(n))
         end do

         KE = KE * S_ij
      end if

   end function overlap_KE_trajectory

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function overlap_S_dot_trajectory(T1, T2, S_ij_precalc) result(S_dot)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory) :: T1, T2
      complex(kind=DefComp), optional :: S_ij_precalc
      complex(kind=DefComp) :: S_dot

      complex(kind=DefComp) :: S_ij

      if (present(S_ij_precalc)) then
         S_ij = S_ij_precalc
      else
         S_ij = overlap_trajectory(T1, T2)
      end if

! I don't understand the need for the minus
! but it is important, the norm is not conserved otherwise
!
! TODO(danielhollas): On HDF_free branch, there's this additional comment,
! but there's no minus sign anymore! We should make sure the formula below is correct.
!
! Tuesday, 20 December 2011 12:01:21
!  Worked out the problem:
!     overlap_dx   = < g1 | d/dx | g2 >
! This equation actually needs
!     overlap_dx_2 = < g1 | d/dx_2 | g2 >
! A little algebra shows
!     overlap_dx  = -overlap_dx_2
      S_dot = (-dot_product(T2%get_vel(), &
                            overlap_dx_trajectory(T1, T2, S_ij)) &
               + dot_product(FMS_Forces(T2), &
                             overlap_dp_trajectory(T1, T2, S_ij)) &
               + (0.d0, 1.d0) * FMS_PhaseDot(T2) * S_ij)

   end function overlap_S_dot_trajectory

   function field_xf(T1) result(ft_xf)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time-dependent electric field
      type(T_Trajectory), intent(in) :: T1
      real(kind=DefReal) :: ft_xf(3), t, pulse

      t = T1%get_time()
      pulse = exp(-(t - t0_xf)**2 / (2 * sigma_xf * sigma_xf)) * ( &
              cos(freq_xf * t + CEP_xf) - &
              sin(freq_xf * t + CEP_xf) * (t - t0_xf) / (sigma_xf * sigma_xf * freq_xf))

      ft_xf(1) = f0_xf * polx_xf * pulse
      ft_xf(2) = f0_xf * poly_xf * pulse
      ft_xf(3) = f0_xf * polz_xf * pulse

   end function field_xf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute the nuclear dipole moment between two TBF
   function dipole_nucl_xf(T1, T2) result(dip_nuc_xf)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      complex(kind=DefComp) :: dip_nuc_xf(3)
      integer :: npart, n

      npart = T1%NumParticles

      dip_nuc_xf = (0.0d0, 0.0d0)
      do n = 1, npart
         dip_nuc_xf(:) = dip_nuc_xf(:) + &
                         nuc_dip_particlexf(T1%Particle(n), &
                                            T2%Particle(n), T1%Particle(n)%AtomicNum)
      end do

   end function dipole_nucl_xf

   function nuc_dip_particlexf(P1, P2, charge) result(S_ij)
! Calculates the overlap between two particles
      type(T_Particle), intent(in) :: P1, P2
      complex(kind=DefComp) :: S_ij(P1%NumDimensions), Rterm

      integer(kind=DefInt) :: i

      real(kind=DefReal) :: x_1, p_1, a_1, &
                            x_2, p_2, a_2, &
                            delta_x, delta_p, x_cent, charge

      if (P1%NumDimensions /= 3) then
         call FMS_DieError("nuc_dip_particlexf only implemented for 3D systems")
      end if

      a_1 = P1%width
      a_2 = P1%width

      S_ij = (1.0d0, 0.0d0)

      do i = 1, P1%NumDimensions
         x_1 = P1%get_pos(i)
         x_2 = P2%get_pos(i)
         p_1 = P1%get_mom(i)
         p_2 = P2%get_mom(i)

         delta_x = x_1 - x_2
         delta_p = p_1 - p_2
         x_cent = (a_1 * x_1 + a_2 * x_2) / (a_1 + a_2) ! centroid position
         Rterm = cmplx(x_cent, (0.25d0 / a_1) * delta_p, kind=DefComp) * charge

         S_ij(i) = Rterm * overlap_CG(x_1, p_1, a_1, x_2, p_2, a_2)
      end do

   end function nuc_dip_particlexf

end module OverlapModule
