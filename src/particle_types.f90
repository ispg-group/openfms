!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Uses element names to assign atomic number, mass, and gaussian
!!    width to each atom
!!
!!
!!    If custom atom types are specified in Control.dat, they take
!!    precedence over the default values. The defaults are for the most
!!    common isotope, with gaussian width from Alexis Thompson et al
!!
!!    @ingroup initial
!<
subroutine FMS_ParticleTypes(T1)
   ! TODO(DH): Pass iCUnit as a parameter
   use GlobalModule, only: DefInt, DefReal, fmiOut, iCUnit, MassToAU, &
                           FMS_DieError, lower_case
   use TrajectoryModule
   implicit none
   type(T_Trajectory), intent(inout) :: T1

   integer(DefInt), allocatable :: nAtomicNum(:)
   real(DefReal), allocatable :: dMass(:), dWidth(:)
   character(len=2), allocatable :: cName(:)
   character(len=256) :: cjunk
   integer(DefInt) :: iType, nTypes, iPtcle, ios
   logical :: zFound

!     Get user-defined atom types from the bottom of Control.dat
!     NOTE: This is an ugly hack, we get an open file via a global
!     icUnit variable, and we suppose we are at the end of
!     the main 'control' namelist which is read in ReadNameList.f.
!     We should make this better somehow!
   nTypes = 0
   read (icUnit, *, iostat=ios) nTypes
   if (ios == 0) then
      read (icUnit, *) cjunk
      allocate (cName(nTypes))
      allocate (nAtomicNum(nTypes))
      allocate (dMass(nTypes), dWidth(nTypes))
      ! TODO: It would be nice to verify that cName does not have duplicates
      do iType = 1, nTypes
         read (icUnit, *) cName(iType), nAtomicNum(iType), dMass(iType), dWidth(iType)
      end do
   else
      write (fmiout, *) 'No user-defined atom types found in Control.dat.'
      write (fmiout, *) 'Defaults will be used.'
      nTypes = 0
   end if
   close (icUnit)

!
!     Assign Atomic Number, Gaussian Width and Atomic Mass to each
!     atom based on its name
!
   do iPtcle = 1, T1%NumParticles

!     Did user specify a width and mass?
      zFound = .false.
      do iType = 1, nTypes
         if (trim(T1%Particle(iPtcle)%Elmnt) == trim(cName(iType))) then
            zFound = .true.
            T1%Particle(iPtcle)%Mass = dMass(iType) * MassToAU
            T1%Particle(iPtcle)%AtomicNum = nAtomicNum(iType)
            T1%Particle(iPtcle)%Width = dWidth(iType)
         end if
      end do
      if (zFound) cycle

!        For no user-specified width and mass, use defaults
      select case (lower_case(trim(T1%Particle(iPtcle)%Elmnt)))
      case ('h')
         T1%Particle(iPtcle)%Width = 4.5d0
         T1%Particle(iPtcle)%Mass = 1.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 1.0d0
      case ('li')
         T1%Particle(iPtcle)%Width = 4.7d0
         T1%Particle(iPtcle)%Mass = 6.9d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 3.0d0
      case ('c')
         T1%Particle(iPtcle)%Width = 22.5d0
         T1%Particle(iPtcle)%Mass = 12.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 6.0d0
      case ('n')
         T1%Particle(iPtcle)%Width = 19.5d0
         T1%Particle(iPtcle)%Mass = 14.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 7.0d0
      case ('o')
         T1%Particle(iPtcle)%Width = 13.0d0
         T1%Particle(iPtcle)%Mass = 16.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 8.0d0
      case ('f')
         T1%Particle(iPtcle)%Width = 8.5d0
         T1%Particle(iPtcle)%Mass = 19.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 9.0d0
      case ('s')
         T1%Particle(iPtcle)%Width = 17.5d0
         T1%Particle(iPtcle)%Mass = 32.0d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 16.0d0
      case ('cl')
         T1%Particle(iPtcle)%Width = 8.5d0
         T1%Particle(iPtcle)%Mass = 35.45d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 17.0d0
      case ('cr')
         ! As reported in https://doi.org/10.1021/acs.jpclett.2c03295
         T1%Particle(iPtcle)%Width = 29.49d0
         T1%Particle(iPtcle)%Mass = 51.99d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 24.0d0
      case ('br')
         ! As reported in https://doi.org/10.1021/acs.jpca.9b00952
         T1%Particle(iPtcle)%Width = 36.69d0
         T1%Particle(iPtcle)%Mass = 79.90d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 35.0d0
      case ('cs')
         ! As reported in https://doi.org/10.1021/acs.jpca.9b00952
         T1%Particle(iPtcle)%Width = 20.96d0
         T1%Particle(iPtcle)%Mass = 132.91d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 55.0d0
      case ('pb')
         ! As reported in https://doi.org/10.1021/acs.jpca.9b00952
         T1%Particle(iPtcle)%Width = 24.63d0
         T1%Particle(iPtcle)%Mass = 207.00d0 * MassToAU
         T1%Particle(iPtcle)%AtomicNum = 82.0d0
      case default
         write (fmiOut, *) 'Particle number: ', iPtcle
         write (fmiOut, *) 'Atom type: ', T1%Particle(iPtcle)%Elmnt
         call FMS_DieError('Unrecognized atom type '//T1%Particle(iPtcle)%Elmnt)
      end select

   end do

end subroutine FMS_ParticleTypes
