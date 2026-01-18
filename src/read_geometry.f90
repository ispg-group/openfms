!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Reads in the Geometry.dat file, loads geometry into a trajectory.
!!    The file is formatted as follows:
!!    <pre>
!!    UNITS=[BOHR or ANG]
!!    [NumberOfAtoms]
!!    [Name1] [x1] [y1] [z1]
!!    [Name2] [x2] [y2] [z2]
!!    .       .    .    .
!!    .       .    .    .
!!    .       .    .    .
!!    # Comment line - momenta may be listed below
!!    [px1]  [py1]  [pz1]
!!    [px2]  [py2]  [pz2]
!!    .     .     .
!!    .     .     .
!!    .     .     .
!!    </pre>
!!     Momenta are optional unless inInitialCond=0
!!    @ingroup initial
!<
subroutine FMS_ReadGeometry(T1)
      use GlobalModule, only: DefInt, DefReal, fmiOut, FmsWorkingDir, &
         glzConstrain, BohrToAng, MassToAu, FMS_DieError
      use TrajectoryModule
      use RattleModule, only: Rattle_ReadConstraints, Rattle_SetConstraints
      implicit none
      type (T_Trajectory), intent(inout) :: T1
      character(len=120) :: atomName
      character(len=256) :: FilePath, unitname
      character(len=1) :: txt
      Logical :: ConUnits
      integer (kind=DefInt) :: I, J, IGUnit, natoms
      integer :: ios
      real (kind=DefReal) :: PosTemp(3)

!
!     Open file
!
      FilePath=trim(FMSWorkingDir)//'Geometry.dat'
      write(fmiOut,*) 'Reading Geometry from '//trim(FilePath)

      open(newunit=IGUnit, file=FilePath, status='old', action='read')
!
!     Make sure that unit system is specified
!
      ConUnits=.false.
      read(IGUnit,*)unitname
      txt=unitname(7:7)
      if (txt == "a" .or. txt == "A") Then
         ConUnits = .true.
      elseif (txt == "b" .or. txt == "B") Then
         ConUnits = .false.
      else ! If this line is absent code will stop
         write(fmiOut,*) 'First line in Geometry.dat must define the units'
         write(fmiOut,*)'It should be either: UNITS=ANG or UNITS=BOHR'
         call FMS_DieError('ERROR: no UNITS line in Geometry.dat')
      endif


!
!     Read and check the number of atoms
!
      read(IGUnit,*)natoms
      if (nAtoms /= T1%NumParticles) then
         call FMS_DieError('ERROR: Geometry.dat has different number of particles than Control.dat')
      end if

!
!     Read in the geometry
!
      Do I=1,natoms
         read(IGUnit,*)atomName,(PosTemp(J),J=1,3)
         if(ConUnits) PosTemp=PosTemp/BohrToAng
         call T1%Particle(I)%set_pos(PosTemp)
         T1%Particle(I)%Elmnt=AtomName(1:2)
      enddo

      if( glzConstrain )then
         call Rattle_ReadConstraints()
         call Rattle_SetConstraints(T1)
      endif

!
!     Read in momenta, if provided
!     Otherwise, set all momenta to 0
!
      read(igUnit,*,iostat=ios) unitname
      if(ios == 0) then
         write(fmiOut,*) 'Reading momenta from Geometry.dat...'
         Do I=1,natoms
            read(IGUnit,*,iostat=ios) (PosTemp(J),J=1,3)
            if(ConUnits) PosTemp=PosTemp/(BohrToAng*MassToAu)
            call T1%Particle(I)%set_mom(PosTemp)
            if(ios /= 0 .and. I/=nAtoms) then
               write(fmiOut,*) 'Error reading momenta from Geometry.dat'
               write(fmiOut,*) 'Setting all  momenta to 0.'
            endif
         enddo
      endif

      if(ios /= 0) then
         Do I=1,natoms
            call T1%Particle(I)%set_mom([0.0d0,0.0d0,0.0d0])
         enddo
      endif

      close(IGUnit)

end subroutine FMS_ReadGeometry
