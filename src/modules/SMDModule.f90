!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!!    Steered Molecular Dynamics for Mechanochemistry Simulation
!!    \author Mitchell Ong and Todd J. Martinez
!>
!!    @brief Parameters for steered molecular dynamics
!<
module SMDModule
   use GlobalModule, only: DefReal, DefInt, FMS_DieError
   use TrajectoryModule
   implicit none

   private
   public :: SMD_Mechano, SMD_Print, SMD_Completed

   integer, parameter :: MaxAtoms = 100 !< Maximum number of atoms
   integer, parameter :: MaxSize = 100

   logical, public :: smSMD = .false. !< Use SMD?
   logical, public :: smVelPull = .false. !< Constant Velocity SMD?

   real(kind=DefReal), public :: smForceConst !< Force Constant (Const. Vel)
   real(kind=DefReal), public :: smPullRate !< Pulling Rate (Const. Vel)
   real(kind=DefReal), public :: smForce !< Force (Const. Force)
   integer(kind=DefInt), public :: smIAtoms(MaxAtoms), smNAtoms !< Label of SMD Atoms
   real(kind=DefReal), public :: smIDummy(3, MaxAtoms), smNDummy !< Position of dummy atoms
   real(kind=DefReal), public :: smElongCut !< Elongation factor cutoff
   logical, public :: smautodirect = .true. !< Direction calculated along SMD vector

   real(kind=DefReal) :: smGrad(3, MaxAtoms) !< External force applied to SMD atoms
   real(kind=DefReal) :: smElongFac(MaxSize) !< Elongation factor
   real(kind=DefReal) :: smNvec(3, MaxAtoms) !< Direction vector
   real(kind=DefReal) :: smExtEnergy !< SMD External Energy
   real(kind=DefReal) :: smIntEnergy !< Potential Energy w/o modification

contains

!     SMD public routine to calculate external force if specified
!>
!!    Calculate the elongation factor ((d-d0)/d0)*100
!<
   subroutine SMD_ElongFac(T1, CurrPos, IPos, NPairs)
      type(T_Trajectory) :: T1

!  Parameters
      integer, parameter :: maxsize = 100
      integer(kind=DefInt) :: NPairs
      real(kind=DefReal) :: CurrPos(3, smNAtoms)
      real(kind=DefReal) :: IPos(3, smNAtoms)
      real(kind=DefReal) :: sum(maxsize), isum(maxsize)
      real(kind=DefReal) :: length(maxsize), ilength(maxsize)
      integer :: i, j

!  Initialize parameters
      sum = 0.0d0
      isum = 0.0d0
      length = 0.0d0
      ilength = 0.0d0

!  Calculate the current and initial distance between 2 SMD Atoms
      do i = 1, NPairs
         do j = 1, 3
            sum(i) = sum(i) + (CurrPos(j, 2 * i - 1) - CurrPos(j, 2 * i))**2
            isum(i) = isum(i) + (IPos(j, 2 * i - 1) - IPos(j, 2 * i))**2
         end do
      end do

!  Calculate the elongation factor
      do i = 1, NPairs
         length(i) = sqrt(sum(i))
         ilength(i) = sqrt(isum(i))
         smElongFac(i) = ((length(i) - ilength(i)) / ilength(i)) * 100
      end do

   end subroutine SMD_ElongFac
!>
!!    Calculates time-dependent external force for contant-velocity pulling
!!
!!    Steered Molecular Dynamics for Mechanochemistry Simulation
!!
!!    This routine calculates a time dependent external force magnitude
!!    for two SMD atoms that will be added to the gradient calculated
!!    by MolPro.  Two atoms will be pulled in opposite directions by a
!!    harmonic spring with user-specified velocity and force constant.
!!
!!    \author Written by Mitchell Ong and Todd J. Martinez
!<
   function SMD_ConstVel(T1, CurrPos, IPos) result(FComp)
      use GlobalModule
      type(T_Trajectory) :: T1
      real(kind=DefReal) :: IPos(3, smNAtoms)
      real(kind=DefReal) :: CurrPos(3, smNAtoms)
      real(kind=DefReal) :: dotrn(smNAtoms)
      real(kind=DefReal) :: fc
      real(kind=DefReal) :: vel
      real(kind=DefReal) :: FComp(smNAtoms)
      integer :: i, j

!  Initialize variables
      dotrn = 0.0
      fc = 0.0
      vel = 0.0
      FComp = 0.0

!  Unit conversion of force constant (H/bohr^2) and pull rate (bohr/fs)
      fc = smForceConst * (BohrToAng**2) * kcalPMtoH
      vel = smPullRate / (BohrToAng * (1.0d3))

!  Calculate the dot product of the displacement and direction vectors
      do i = 1, smNAtoms
         do j = 1, 3
            dotrn(i) = dotrn(i) + (CurrPos(j, i) - IPos(j, i)) * smNvec(j, i)
         end do
      end do

!  Calculate the external force magnitude
      do i = 1, smNAtoms
         FComp(i) = fc * (vel * (T1%get_time() / FsToAU) - dotrn(i))
      end do

   end function SMD_ConstVel

!>
!!    Steered Molecular Dynamics for Mechanochemistry Simulation
!!    \author Mitchell Ong and Todd J. Martinez
!!    \todo Document
!<
   subroutine SMD_Direction(T1, CurrPos, NPairs)
      type(T_Trajectory) :: T1
      integer, parameter :: MaxSize = 100
      real(kind=DefReal) :: CurrPos(3, smNAtoms)
      real(kind=DefReal) :: sum(MaxSize)
      real(kind=DefReal) :: norm(MaxSize)
      integer :: i, j, NPairs

!  Initialize Parameters
      sum = 0.0d0
      norm = 0.0d0

      if (smautodirect) then

!  Calculate the distance between 2 SMD Atoms
         do i = 1, NPairs
            do j = 1, 3
               sum(i) = sum(i) + (CurrPos(j, 2 * i - 1) - CurrPos(j, 2 * i))**2
            end do
         end do

         do i = 1, NPairs
            norm(i) = sqrt(sum(i))
         end do

!  Calculate the direction vector of the force
         do i = 1, NPairs
            do j = 1, 3
               smNvec(j, 2 * i - 1) = (CurrPos(j, 2 * i - 1) - CurrPos(j, 2 * i)) / norm(i)
               smNvec(j, 2 * i) = (CurrPos(j, 2 * i) - CurrPos(j, 2 * i - 1)) / norm(i)
            end do
         end do

      else

!  Calculate the distance between SMD and Dummy atoms
         do i = 1, smNAtoms
            do j = 1, 3
               sum(i) = sum(i) + (smIDummy(j, i) - CurrPos(j, i))**2
            end do
         end do

         do i = 1, smNAtoms
            norm(i) = sqrt(sum(i))
         end do

!  Calculate the direction vector of the force
         do i = 1, smNAtoms
            do j = 1, 3
               smNvec(j, i) = (smIDummy(j, i) - CurrPos(j, i)) / norm(i)
            end do
         end do

      end if

   end subroutine SMD_Direction

!>
!!    Steered Molecular Dynamics for Mechanochemistry Simulation
!!    \author Mitchell Ong and Todd J. Martinez
!!    \todo Document
!<
   subroutine SMD_Mechano(T1)
      type(T_Trajectory) :: T1

!     Parameters
      integer, parameter :: maxsize = 100
      real(kind=DefReal) :: CurrPos(3, smNAtoms)
      real(kind=DefReal) :: IPos(3, smNAtoms)
      real(kind=DefReal) :: InitPos(T1%NumParticles, 3)
      real(kind=DefReal) :: ForceMag(smNAtoms)
      real(kind=DefReal) :: rdiff1, rdiff2
      integer :: i, j, NPairs

!     Initialize Parameters
      CurrPos = 0.0d0
      IPos = 0.0d0
      InitPos = 0.0d0
      ForceMag = 0.0d0
      smIntEnergy = 0.0d0
      smExtEnergy = 0.0d0

      do i = 1, smNAtoms
         IPos(1:3, i) = T1%Particle(smIAtoms(i))%get_pos()
      end do

!     Read in current position of SMD atoms
      do i = 1, smNAtoms
         CurrPos(1:3, i) = T1%Particle(smIAtoms(i))%get_pos()
      end do

!     Input error checking
      if (smNDummy /= smNAtoms .and. (.not. smautodirect)) then
         call FMS_DieError('Number of dummy and SMD atoms not equal')
      end if

      if (mod(smNAtoms, 2) /= 0) then
         call FMS_DieError('Number of SMD atoms must be even')
      end if

!     Determine if elongation factor is bigger than cutoff
!     If so, terminate the dynamics simulation
      NPairs = smNAtoms / 2
      call SMD_ElongFac(T1, CurrPos, IPos, NPairs)

!     Calculate the direction vector for each SMD atom
      call SMD_Direction(T1, CurrPos, NPairs)

!     Calculate the external force for constant velocity or constant force
      if (smVelPull) then
         ForceMag = SMD_ConstVel(T1, CurrPos, IPos)
         do i = 1, smNAtoms
            do j = 1, 3
               smGrad(j, i) = ForceMag(i) * smNvec(j, i)
            end do
         end do
      else
         do i = 1, smNAtoms
            do j = 1, 3
               smGrad(j, i) = smForce * smNvec(j, i)
            end do
         end do
      end if

!     Add external force to force vector
      do i = 1, smNAtoms
         do j = 1, 3
            T1%ElecStruc%ModForce(smIAtoms(i) * 3 + j - 3) = T1%ElecStruc%ModForce(smIAtoms(i) * 3 + j - 3) + smGrad(j, i)
         end do
      end do

!     Calculate external energy
      if (smautodirect) then
         do i = 1, NPairs
            rdiff1 = 0.0d0
            do j = 1, 3
               rdiff1 = rdiff1 + (CurrPos(j, 2 * i - 1) - CurrPos(j, 2 * i))**2
            end do
            smExtEnergy = smExtEnergy - smForce * sqrt(rdiff1)
         end do
         smIntEnergy = T1%ElecStruc%PotEn(T1%StateID)
         T1%ElecStruc%ModPot = smExtEnergy
      else
         do i = 1, smNAtoms
            rdiff1 = 0.0d0
            rdiff2 = 0.0d0
            do j = 1, 3
               rdiff1 = rdiff1 + (smIDummy(j, i) - CurrPos(j, i))**2
               rdiff2 = rdiff2 + (smIDummy(j, i) - IPos(j, i))**2
            end do
            smExtEnergy = smExtEnergy + smForce * (sqrt(rdiff1) - sqrt(rdiff2))
         end do
         smIntEnergy = T1%ElecStruc%PotEn(T1%StateID)
         T1%ElecStruc%ModPot = smExtEnergy
      end if

   end subroutine SMD_Mechano

   subroutine SMD_Print(current_time, write_header)
      use GlobalModule, only: FMSWorkingDir
      real(DefReal), intent(in) :: current_time
      logical, intent(in) :: write_header
      character(len=256) :: cfname
      integer :: iUnit, i, j

      cfname = trim(FMSWorkingDir)//'SMD.dat'
      open (newunit=IUnit, file=trim(cfname), position='APPEND', form='FORMATTED')
      if (write_header) then
         write (IUnit, '(3A)') '     #Time  ElongFac       dx(1)', &
            '       dy(1)       dz(1)       dx(2)       dy(2)', &
            '       dz(2)     ESPotEn    SMDExtEn    TotPotEn'
      end if
      write (IUnit, '(2f10.2,100f12.6)') current_time, (smElongFac(i), &
                                                        i=1, smNAtoms / 2), ((smGrad(j, i), j=1, 3), i=1, smNAtoms), &
         smIntEnergy, smExtEnergy, smIntEnergy + smExtEnergy

      close (IUnit)
   end subroutine SMD_Print

   function SMD_Completed() result(completed)
      use GlobalModule, only: FMS_PrintMessg
      logical :: completed
      integer(defInt) :: iPtcle

      completed = .false.

      do iPtcle = 1, smnAtoms / 2
         if (smElongFac(iPtcle) > smElongCut) then
            call FMS_PrintMessg('SMD: Elongation Factor cutoff reached')
            completed = .true.
            return
         end if
      end do

   end function SMD_Completed

end module SMDModule
