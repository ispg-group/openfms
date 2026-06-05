! Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!   @brief Stores parameters related to QM/MM simulations
!!   @note Module code: qc
!<
module QM_MM_Module
   use GlobalModule, only: defInt, defReal, MaxParticles
   implicit none
   public
   private :: defInt, defReal, MaxParticles

   logical :: qcZQMMM = .false. !< Is this a QM/MM simulation?
   integer(kind=DefInt) :: qcNumMM = 0 !< Number of MM atoms
   logical :: qcZFixAtom = .false. !< Fix certain atoms' positions
   integer(kind=DefInt) :: qciFixAtom(MaxParticles) !< Which atoms to fix
   integer(kind=DefInt) :: qcNumQM !< Number of QM particles
   integer(kind=DefInt) :: qcNumAS !< Atoms per solvent (for GenSolvent)
   integer(kind=DefInt) :: qcNumS !< Number of solvent molecules (for GenSolvent)

!< Add partial charge interactions between QM and MM systems? (only for
!! cases where the QM calculations are not used; for instance EquiCon and MM minimization)
   logical :: qczPCharge = .false.

!> Apply spherical confining potential of the form  \f$ V= k(R-R_0)^2 \f$?
   logical :: qczConfine
   real(kind=DefReal) :: qcdConfineD !< Radial confinement distance \f$ R_0 \f$
   real(kind=DefReal) :: qcdConfineK !< Radial confinement force constant \f$ k \f$
end module QM_MM_Module
