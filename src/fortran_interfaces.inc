! Fortran interfaces for free-standing functions defined in src/.f90 files
! This file is included and exposed via the FMSModule
! IMPORTANT: This file needs to be kept in sync with the actual subroutines!

   recursive subroutine FMS_Dynamics(Bundle, Timestep, GoalTime)
      use GlobalModule, only: DefReal
      use BundleModule
      type(T_TrajectoryBundle), intent(inout) :: Bundle
      real(kind=DefReal), intent(in) :: GoalTime, Timestep
   end subroutine FMS_Dynamics

   subroutine FMS_InitialCond(B1)
      use BundleModule
      type (T_TrajectoryBundle), intent(inout) :: B1
   end subroutine FMS_InitialCond

   subroutine FMS_MinSearch(T1, iState)
      use GlobalModule, only: DefInt
      use TrajectoryModule
      type(T_Trajectory), intent(inout) :: T1
      integer (kind=DefInt), optional, intent(in) :: iState
   end subroutine FMS_MinSearch

   subroutine FMS_ParticleTypes(T1)
      use TrajectoryModule
      type(T_Trajectory), intent(inout) :: T1
   end subroutine FMS_ParticleTypes

   subroutine FMS_ReadGeometry(T1)
      use TrajectoryModule
      type (T_Trajectory), intent(inout) :: T1
   end subroutine FMS_ReadGeometry

   subroutine FMS_ReadNameList(NumParticles,NumStates, NumTraj, SimulationTime)
      use GlobalModule, only: DefInt, DefReal
      integer (kind=DefInt) :: NumParticles, NumStates, NumTraj
      real (kind=DefReal) :: SimulationTIme
   end subroutine FMS_ReadNameList

   subroutine FMS_Shutdown(B1, T1, terminate)
      use TrajectoryModule
      use BundleModule
      type(T_TrajectoryBundle), optional, intent(inout) :: B1
      type(T_Trajectory), optional, intent(inout) :: T1
      logical, optional, intent(in) :: terminate
   end subroutine FMS_Shutdown
