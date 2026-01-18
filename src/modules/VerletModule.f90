!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Subroutines for (velocity) Verlet integration
!!
!! NOTE: These subroutines cannot live in PropagationModule since they
!! would introduce cyclic dependency between SpawnModule and PropagationModule.
!!
!<
module VerletModule
      use GlobalModule, only: DefReal, glzConstrain, D2, Pi
      use TrajectoryModule
      use TrajectoryCalcsModule, only: FMS_GetForce, FMS_PhaseDot
      use RattleModule, only: rattle_constrain_position, rattle_constrain_momentum
      implicit none

      private
      public :: FMS_PropVV_a, FMS_PropVV_b, FMS_PropVV

      contains
!>
!!    Propagates classical positions
!!    using a velocity Verlet integrator
!!
!!    This routine has been modified by Aaron:
!!
!!    Unlike PropQCVV, which propagates
!!    \f{eqnarray*}
!!    R(t)      & \mapsto & R(t+\delta t) \\
!!    p(t-\delta t/2) & \mapsto & p(t+\delta t/2)
!!    \f}
!!    this routine propagates both R AND p from t to t+dt:
!!    \f{eqnarray*}
!!    R(t) & \mapsto & R(t+\delta t) \\
!!    p(t) & \mapsto & p(t+\delta t)
!!    \f}
!!
!!    This requires correction of the velocities - on the first step
!!    p(t-dt/2) is first moved up to p(t) before the standard VV
!!    propagation.
!<
      subroutine FMS_PropVV_a(T1,LastTraj,TimeStep,g1_0,g2_0)
      type (T_Trajectory),intent(inout) :: T1
      real (kind=DefReal),intent(in)    :: TimeStep

      type (T_Trajectory), intent(inout) :: LastTraj ! for Rattle

      logical, save :: FirstTime=.true.

      real (kind=DefReal),dimension(T1%NumDimensions) :: &
            M,           &! atomic masses
!        current    step
            R0,      R1, &! position
            V0,          &! velocity
            P0,      P1, &! momementum
            F0      ! forces

      real (kind=DefReal), intent(inout) :: g1_0, g2_0

      if(FirstTime) then
         if(glzConstrain) call LastTraj%create(T1%NumParticles, T1%NumStates)
         FirstTime = .false.
      endif

      if(glzConstrain) LastTraj=T1

      ! Get all the initial information
      M  = T1%get_mass()
      R0 = T1%get_pos()
      V0 = T1%get_vel()
      P0 = T1%get_mom()
      F0 = FMS_GetForce   (T1)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Nuclear phase  update

      !T1%Phase = T1%Phase + TimeStep/2.d0*FMS_PhaseDot(T1)
!     g1_0 = FMS_PhaseDot(T1)
!     g2_0 = 2.d0 * dot_product( FMS_GetForce   (T1),
!    &                           T1%get_vel() )
      g1_0 = FMS_PhaseDot(T1)
      g2_0 = 2.d0 * dot_product(F0,V0)

      !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Position update
      ! The momentum update is split into to steps for the rattle algorithm
      ! If not using rattle only one momentum update would be needed
      R1 = R0 + TimeStep*V0 + TimeStep**2/2.d0 * F0/M
      P1 = P0 + TimeStep/2.d0 * F0

      call T1%set_pos(R1)
      call T1%set_mom(P1)


      if (glzConstrain) then
         call Rattle_Constrain_Position(T1, LastTraj, TimeStep)
      end if

      end subroutine FMS_PropVV_a
!>
!!    Propagates classical positions
!!    using a velocity Verlet integrator
!!
!!    This routine has been modified by Aaron:
!!
!!    Unlike PropQCVV, which propagates
!!    \f{eqnarray*}
!!    R(t)      & \mapsto & R(t+\delta t) \\
!!    p(t-\delta t/2) & \mapsto & p(t+\delta t/2)
!!    \f}
!!    this routine propagates both R AND p from t to t+dt:
!!    \f{eqnarray*}
!!    R(t) & \mapsto & R(t+\delta t) \\
!!    p(t) & \mapsto & p(t+\delta t)
!!    \f}
!!
!!    This requires correction of the velocities - on the first step
!!    p(t-dt/2) is first moved up to p(t) before the standard VV
!!    propagation.
!<
      subroutine FMS_PropVV_b(T1,LastTraj,TimeStep,g1_0,g2_0)
      type (T_Trajectory),intent(inout) :: T1
      real (kind=DefReal),intent(in)    :: TimeStep

      type (T_Trajectory), intent(in) :: LastTraj ! for Rattle

      logical, save :: FirstTime=.true.

      real (kind=DefReal),dimension(T1%NumDimensions) :: &
            M,           &! atomic masses
!        current    step
            R0,          &! position
            P0,      P1, &! momementum
            F1            ! forces

      real (kind=DefReal), intent(inout) :: g1_0, g2_0
      real (kind=DefReal) :: g1_1, g2_1

      real(DefReal) :: alpha, beta, aa,bb,cc,dd, B(2)

      if(FirstTime) then
         FirstTime = .false.
      endif

      ! Get all the initial information
      M  = T1%get_mass()
      R0 = T1%get_pos()
      P0 = T1%get_mom()
      F1 = FMS_GetForce   (T1)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !                     Momentum update

      P1 = P0 + TimeStep/2.d0 * F1

      call T1%set_mom(P1)

      if (glzConstrain) then
         call Rattle_Constrain_Momentum( T1, LastTraj, TimeStep )
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Nuclear phase  update
      !T1%Phase = T1%Phase + TimeStep/2.d0*FMS_PhaseDot(T1)

      g1_1 = FMS_PhaseDot(T1)
      g2_1 = 2.d0 * dot_product(FMS_GetForce(T1), T1%get_vel())


      ! set up the linear system for the gammas
      aa = timestep**2/2.d0
      bb = timestep**3/6.d0
      cc = timestep
      dd = timestep**2/2.d0

      B = [ g1_1 - g1_0 - g2_0 * timestep, g2_1 - g2_1 ]

      alpha = ( dd*B(1)-bb*B(2) ) / (aa*dd-bb*cc)
      beta =  (-cc*B(1)-aa*B(2) ) / (aa*dd-bb*cc)

      T1%Phase = T1%Phase &
         + ( g1_0 + g1_1 ) / 2.d0 * TimeStep &
         - ( g2_0 - g2_1 ) / 8.d0 * TimeStep**2

      call T1%rescale_phases()
      end subroutine FMS_PropVV_b

!>
!!    Propagates classical parameters (and nuclear phase)
!!    using a velocity Verlet integrator
!!
!!    This routine has been modified by Aaron:
!!
!!    Unlike PropQCVV, which propagates
!!    \f{eqnarray*}
!!    R(t)      & \mapsto & R(t+\delta t) \\
!!    p(t-\delta t/2) & \mapsto & p(t+\delta t/2)
!!    \f}
!!    this routine propagates both R AND p from t to t+dt:
!!    \f{eqnarray*}
!!    R(t) & \mapsto & R(t+\delta t) \\
!!    p(t) & \mapsto & p(t+\delta t)
!!    \f}
!!
!!    This requires correction of the velocities - on the first step
!!    p(t-dt/2) is first moved up to p(t) before the standard VV
!!    propagation.
!<
subroutine FMS_PropVV(T1,TimeStep)
      type (T_Trajectory),intent(inout) :: T1
      real (kind=DefReal),intent(in)    :: TimeStep

      type (T_Trajectory)            :: LastTraj ! for Rattle

      logical, save :: FirstTime=.true.

      real (kind=DefReal), dimension(T1%NumDimensions) :: &
            M,      & ! atomic masses
            R0, R1, & ! position
            V0,     & ! velocity
            P0, P1, & ! momentum
            F0, F1    ! forces

      real (kind=DefReal) :: g1_0, g1_1, g2_0, g2_1

      real(DefReal) :: alpha, beta, aa, bb, cc, dd, B(2)

      ! TODO(DH): This is super weird, why is this needed???
      if(FirstTime) then
         if (glzConstrain) then
            call LastTraj%create(T1%NumParticles, T1%NumStates)
         end if
         FirstTime = .false.
      endif

      if(glzConstrain) LastTraj=T1

      ! Get all the initial information
      M  = T1%get_mass()
      R0 = T1%get_pos()
      V0 = T1%get_vel()
      P0 = T1%get_mom()
      F0 = FMS_GetForce(T1)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Nuclear phase  update

      !T1%Phase = T1%Phase + TimeStep/2.d0*FMS_PhaseDot(T1)
!     g1_0 = FMS_PhaseDot(T1)
!     g2_0 = 2.d0 * dot_product( FMS_GetForce   (T1),
!    &                           FMS_GetVelocity(T1) )
      g1_0 = FMS_PhaseDot(T1)
      g2_0 = 2.d0 * dot_product(F0,V0)

      !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Position update
      ! The momentum update is split into to steps for the rattle algorithm
      ! If not using rattle only one momentum update would be needed
      R1 = R0 + TimeStep*V0 + TimeStep**2/2.d0 * F0/M
      P1 = P0 + TimeStep/2.d0 * F0

      call T1%set_pos(R1)
      call T1%set_mom(P1)


      if (glzConstrain) then
         call Rattle_Constrain_Position( T1, LastTraj, TimeStep )
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !                     Momentum update
      ! Get the force at the new position (will trigger abinitio)
      P0 = T1%get_mom()
      F1 = FMS_GetForce   ( T1 )

      P1 = P0 + TimeStep/2.d0 * F1

      call T1%set_mom(P1)

      if (glzConstrain) then
         call Rattle_Constrain_Momentum( T1, LastTraj, TimeStep )
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !          Nuclear phase  update
      !T1%Phase = T1%Phase + TimeStep/2.d0*FMS_PhaseDot(T1)

      g1_1 = FMS_PhaseDot(T1)
      g2_1 = 2.d0 * dot_product( F1, P1/M )

      ! set up the linear system for the gammas
      aa = timestep**2/2.d0
      bb = timestep**3/6.d0
      cc = timestep
      dd = timestep**2/2.d0

      B = [ g1_1 - g1_0 - g2_0 * timestep, g2_1 - g2_1 ]

      alpha = ( dd*B(1)-bb*B(2) ) / (aa*dd-bb*cc)
      beta =  (-cc*B(1)-aa*B(2) ) / (aa*dd-bb*cc)

      T1%Phase = T1%Phase &
               + ( g1_0 + g1_1 ) / 2.d0 * TimeStep   &
               - ( g2_0 - g2_1 ) / 8.d0 * TimeStep**2

      call T1%rescale_phases()

   end subroutine FMS_PropVV

end module VerletModule
