!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! @brief Trajectory class IO methods.
!<
module TrajectoryIOModule

   use GlobalModule
   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_Weight, FMS_CIVec, FMS_Coupling, FMS_SOCoupling, FMS_CoupDotVel, &
                                    Kinetic, Potential, FMS_PotentialT, FMS_MMPot, FMS_Dipole, FMS_TransDipole
   use QM_MM_Module, only: qczQMMM

   implicit none

   private
   public :: FMS_WriteFXYZ
   public :: FMS_WriteTrajFiles

contains
!--------------- Module Procedures ---------------!
!>
!!    Writes XYZ file
!!    \param comment   Comment to write on comment line
!!    @ingroup output
!<
   subroutine FMS_WriteFXYZ(T1, filename, comment, firsttime)
      type(T_Trajectory), intent(in) :: T1
      character(len=*), intent(in) :: filename, comment
      logical, optional :: firsttime

      integer(kind=DefInt) :: IUnit, IPtcle
      character(len=256) :: file_name, pos_str

      file_name = trim(FMSWorkingDir)//filename

      pos_str = 'append'
      if (present(firsttime)) then
         if (firsttime) pos_str = 'rewind'
      end if

      open (newunit=IUnit, file=file_name, position=pos_str)

      write (IUnit, '(I8)') T1%NumParticles
      if (len(trim(comment)) == 0) then
         write (iUnit, *) 'Time ', T1%get_time()
      else
         write (IUnit, *) trim(Comment)
      end if

      do IPTcle = 1, T1%NumParticles
         write (IUnit, '(A3, 3F17.9)') T1%Particle(IPtcle)%Elmnt, T1%Particle(IPtcle)%get_pos() * BohrToAng
      end do
      flush (IUnit)

      close (IUnit)

   end subroutine FMS_WriteFXYZ

!>
!!    Driver for output for each trajectory
!!    Writes files for the classical trajectory, appended by '.suffix'.
!!
!!    \param   Suffix Append suffix to the end of each file to identify
!!          the trajectory
!!    @ingroup output
!<
   subroutine FMS_WriteTrajFiles(T1, suff, zFirst)
      use ElecStrucModule
      type(T_Trajectory), intent(in) :: T1
      character(len=*), intent(in) :: suff
      logical, optional, intent(in) :: zFirst
      character(len=:), allocatable :: suffix

      logical :: zHeader
      character(len=256) :: cfname, comment

      suffix = trim(adjustl(suff))

      zHeader = .false.
      if (present(zFirst)) then
         zHeader = zFirst
      end if

!     All traj files get Amp file and TrajDump, the rest only once per
!     CTBF
      if (T1%Ms /= 2) then
         if (fmzTrajFile .or. fmzAllText) call FMS_WriteFTrajDump(T1)
         if (fmzAmpFile .or. fmzAllText) call FMS_WriteFAmp(T1)

      else

         if (fmzTrajFile .or. fmzAllText) call FMS_WriteFTrajDump(T1)
         if (fmzAmpFile .or. fmzAllText) call FMS_WriteFAmp(T1)
         if (fmzPotEnFile .or. fmzAllText) call FMS_WriteFPotEn(T1)
         if (fmzQMRRFile .or. fmzAllText) call FMS_WriteFQMRR(T1)

!     GAIMS changed
         if ((fmzCoupFile .or. fmzAllText) .and. eszNACoupVec) then
            call FMS_WriteFCouple(T1)
            if (NTrip /= 0) call FMS_WriteFSOCouple(T1)
         end if
!     GAIMS end changed

!     Write MM energy files
         if ((fmzMMFile .or. fmzAllText) .and. qczQMMM .and. eszMMForce) then
            cfname = 'MM.'//trim(suffix)
            call FMS_WriteFMM(T1, cfname, zHeader)
         end if

!     Write partial charges
         if ((fmzChargeFile .or. fmzAllText) .and. eszPartialCharges) then
            cfname = 'Charge.'//trim(suffix)
            call FMS_WriteFCharge(T1, cfname, zHeader)
         end if

!     Write dipole files
         if ((fmzDipoleFile .or. fmzAllText) .and. eszDipoleMoment) then
            cfname = 'Dipole.'//trim(suffix)
            call FMS_WriteFDipole(T1, cfname, zHeader)
         end if

!     Write Transition dipole files
         if ((fmzTDipoleFile .or. fmzAllText) .and. eszTransDipole) call FMS_WriteFTDipole(T1)

!     Write Dihedral.x files
         if (fmNDihedrals /= 0) then
            cfname = 'Dihedrals.'//trim(suffix)
            call FMS_WriteFDihedral(T1, cfname, fmNDihedrals, fmIDihedral, zHeader)
         end if

!     Write Pyram.x files
         if (fmNPyrams /= 0) then
            cfname = 'Pyram.'//trim(suffix)
            call FMS_WriteFPyram(T1, cfname, fmNPyrams, fmIPyram, zHeader)
         end if

!     Write Angle.x files
         if (fmNAngles /= 0) then
            cfname = 'Angles.'//trim(suffix)
            call FMS_WriteFAngle(T1, cfname, fmNAngles, fmIAngle, zHeader)
         end if

!     Write Bond.x files
         if (fmNBonds /= 0) then
            cfname = 'Bonds.'//trim(suffix)
            call FMS_WriteFDistance(T1, cfname, fmNBonds, fmIBond, zHeader)
         end if

!     Write positions.x.xyz files
1030     format('Time: ', f8.2, ', Trajectory: ')
         if (fmzXYZ .or. fmzAllText) then
            cfname = 'positions.'//trim(suffix)//'.xyz'
            write (comment, 1030) T1%get_time()
            call FMS_WriteFXYZ(T1, cfname, trim(comment)//trim(suffix))
         end if

!     Write forces.x.xyz files (in XYZ format)
         if (fmzForce) then
            cfname = 'forces.'//trim(suffix)//'.xyz'
            write (comment, 1030) T1%get_time()
            call FMS_WriteFForces(T1, cfname, trim(comment)//trim(suffix))
         end if

      end if
   end subroutine FMS_WriteTrajFiles
!>
!!    Write formatted trajectory data.
!!    @ingroup output
!<
   subroutine FMS_WriteFTrajDump(T)
      type(T_Trajectory), intent(in) :: T

      character(len=256) :: file_name, tmp
      logical :: file_existed
      integer(DefInt) :: iunit, ntraj
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files

      integer(kind=DefInt) :: i, j

1     format(a10, 4000a10)
2     format(f10.2, 4000f10.4)

      ntraj = T%TrajID
      iunit = units(ntraj)

      ! First time setup
      if (iunit == 0) then
         write (tmp, *) ntraj
         file_name = trim(FMSWorkingDir)//'TrajDump.'//adjustl(trim(tmp))

         inquire (file=file_name, exist=file_existed)
         if (file_existed) then
            open (newunit=iunit, file=trim(file_name), position="append")
         else
            open (newunit=iunit, file=trim(file_name))
            write (iunit, 1) '# Time', ([ &
                                        FMS_NumberedFileName('pos', i)//'x', &
                                        FMS_NumberedFileName('pos', i)//'y', &
                                        FMS_NumberedFileName('pos', i)//'z'], i=1, T%NumParticles), &
               ([FMS_NumberedFileName('mom', i)//'x', &
                 FMS_NumberedFileName('mom', i)//'y', &
                 FMS_NumberedFileName('mom', i)//'z'], i=1, T%NumParticles), &
               'Phase', 'AmpReal', 'AmpImag', 'AmpNorm', 'StateID'
         end if
         units(ntraj) = iunit
      end if

      write (iunit, 2) T%get_time(), (T%Particle(i)%get_pos(), i=1, T%NumParticles), &
         ((T%Particle(i)%get_mom(j), j=1, T%Particle(i)%NumDimensions), i=1, T%NumParticles), &
         T%Phase, real(T%Amplitude), aimag(T%Amplitude), FMS_Weight(T), dble(T%StateID)
      flush (iunit)

   end subroutine FMS_WriteFTrajDump
!>
!!    Write amplitude for each trajectory
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFAmp(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      logical :: file_existed
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files
      ! First time setup
      iunit = units(T%TrajID)

      if (iunit == 0) then
         call FMS_OpenFile(FMS_NumberedFileName('Amp', T%TrajID), iunit, file_existed)
         if (.not. file_existed) then
            write (iunit, '(5(a10))') '#Time', 'Norm', 'Real', 'Imag'
         end if
         units(T%TrajID) = iunit
      end if

      write (iunit, '(f10.2,4(f10.4))') T%get_time(), FMS_Weight(T), T%Amplitude
      flush (iunit)

   end subroutine FMS_WriteFAmp
!>
!!    Writes out specified angles.
!!
!!    Columns are labeled with atom numbers.
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFAngle(T, file_name, NAngles, iAngle, first_time)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      character(len=256), intent(in) :: file_name
      integer(kind=DefInt), intent(in) :: NAngles, IAngle(:, :)
      logical, optional :: first_time

      integer(kind=DefInt), parameter :: MaxAngle = 20
      real(kind=DefReal), parameter :: rad2deg = 180.d0 / 3.141592654d0
      integer(kind=DefInt) :: i

      integer :: iunit
      character(len=256) :: long_file_name
      logical :: file_existed

      !         time             angle
1     format(a15, 20(4x, i3, ",", i3, ",", i3)) ! header format
2     format(3x, f12.2, 20(1x, f14.6)) ! angle format

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (present(first_time) .and. first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (iunit, 1) '# Time', (IAngle(1:3, i), i=1, NAngles)
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

      write (iunit, 2) T%get_time(), (rad2deg * FMS_Angle(T, iAngle(1, i), iAngle(2, i), iAngle(3, i)), i=1, NAngles)

      close (IUnit)
   end subroutine FMS_WriteFAngle
!>
!!   Writes CI vector for each trajectory (molpro only)
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFCIVec(T1, file_name, first_time)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use ElecStrucModule
      type(T_Trajectory) :: T1
      character(len=*), intent(in) :: file_name
      logical, intent(in) :: first_time

      character(len=256) :: long_file_name
      integer(kind=DefInt) :: iunit
      logical :: file_existed

      character(len=16) :: CDet

      integer(kind=DefInt) :: IDet(2 * esNCasOrbs), IPrintDet(1000, 2 * esNCasOrbs), ItoCPHelper(esNCasOrbs)
      real(kind=DefReal) :: PrintCoeff(1000, t1%NumStates)

      integer(kind=DefInt) :: iprint, i, icasorbs, icivec, istate, k

      logical :: zprint, MErr

      select case (gliModel)
      case (10)
         long_file_name = trim(FMSWorkingDir)//file_name

         file_existed = .true.
         if (first_time) then
            inquire (file=long_file_name, exist=file_existed)
         end if

         if (.not. file_existed) then
            open (newunit=iunit, file=trim(long_file_name))
         else
            open (newunit=iunit, file=trim(long_file_name), position='append')
         end if

         ! 4. Do whatever is happening here...
         iprint = 1
         do icasorbs = 1, esNCasElec / 2
            IDet(icasorbs) = 1
            IDet(icasorbs + esNCasOrbs) = 1
         end do
         do icasorbs = (esNCasElec / 2) + 1, esNCasOrbs
            IDet(icasorbs) = 0
            IDet(icasorbs + esNCasOrbs) = 0
         end do
         do icivec = 1, esLCiVec
            zprint = .false.
            do istate = 1, T1%NumStates
               PrintCoeff(iprint, istate) = FMS_CIVec(T1, istate, iCIVec)
               if (PrintCoeff(iprint, iState) > 0.025d0) then
                  do icasorbs = 1, esNCasOrbs * 2
                     IPrintDet(iprint, icasorbs) = IDet(icasorbs)
                  end do
                  zprint = .true.
               end if
            end do
            if (zprint) iprint = iprint + 1
            CDet = FMS_PrintDet(IDet, esNCasOrbs)
            do icasorbs = 1, esNCasOrbs
               ItoCPHelper(icasorbs) = IDet(icasorbs)
            end do
            call FMS_CPHelper(esNCasOrbs, esNCasElec / 2, ItoCPHelper, MErr)
            do icasorbs = 1, esNCasOrbs
               IDet(icasorbs) = ItoCPHelper(icasorbs)
            end do
            if (.not. MErr) then
               do icasorbs = 1, esNCasOrbs
                  ItoCPHelper(icasorbs) = IDet(icasorbs + esNCasOrbs)
               end do
               call FMS_CPHelper(esNCasOrbs, esNCasElec / 2, ItoCPHelper, MErr)
               do icasorbs = 1, esNCasOrbs
                  IDet(icasorbs + esNCasOrbs) = ItoCPHelper(icasorbs)
               end do
            end if
         end do

         write (IUnit, 1000) T1%get_time()
         do i = 1, iprint - 1
            do icasorbs = 1, esNCasOrbs * 2
               IDet(icasorbs) = IPrintDet(i, icasorbs)
            end do
            write (IUnit, 1100) FMS_PrintDet(IDet, esNCasOrbs), (PrintCoeff(i, k), k=1, T1%NumStates)
         end do
         close (IUnit)

      case default
         call FMS_DieError('WriteFCIVec only works with gliModel 10')
      end select

1000  format(f10.2)
1100  format(A14, 100f10.6)

   contains

!>
!!    Steps the IDet vector in CorrectPhase
!!    \param norbs Number of orbitals
!!    \param nelec Number of electrons
!!    \param Success Success or failure
!!    @todo Document iDet and explain when this subroutine fails
!<
      recursive subroutine FMS_CPHelper(norbs, nelec, IDet, Success)

         integer(kind=DefInt), intent(in) :: norbs, nelec
         integer(kind=DefInt) :: iorbs, itmp, i
         integer(kind=DefInt) :: IDet(norbs * 2)
         integer(kind=DefInt) :: ItoCPHelper(norbs - 1)
         logical(kind=DefInt) :: Success, MErr

         iorbs = norbs
100      continue
         if (IDet(iorbs) == 0) then
            iorbs = iorbs - 1
            if (iorbs == 0) then
               Success = .false.
               do iorbs = 1, norbs
                  if (iorbs <= nelec) then
                     IDet(iorbs) = 1
                  else
                     IDet(iorbs) = 0
                  end if
               end do
               return
            end if
            goto 100
         end if
         if (iorbs == norbs) then
            ItoCPHelper(1:norbs - 1) = IDet(1:norbs - 1)
            call FMS_CPHelper(norbs - 1, nelec - 1, ItoCPHelper, MErr)
            Success = MErr
            IDet(1:norbs - 1) = ItoCPHelper(1:norbs - 1)
            if (Success) then
               itmp = 0
               do i = 1, norbs
                  itmp = itmp + IDet(i)
                  if (itmp == nelec - 1) goto 300
               end do
300            continue
               IDet(norbs) = 0
               IDet(i + 1) = 1
               return
            end if
            do iorbs = 1, norbs
               if (iorbs <= nelec) then
                  IDet(iorbs) = 1
               else
                  IDet(iorbs) = 0
               end if
            end do
            return
         end if

         IDet(iorbs) = 0
         IDet(iorbs + 1) = 1
         Success = .true.
         return
      end subroutine FMS_CPHelper

!>
!!    @todo Needs documentation
!<
      function FMS_PrintDet(IDet, ncasorbs) result(Chars)
         integer(kind=DefInt) :: ncasorbs, IDet(2 * ncasorbs)
         character(len=12) :: Chars
         integer(kind=DefInt) :: i

         Chars = '            '
         do i = 1, ncasorbs
            if (IDet(i) == 1) then
               if (IDet(i + ncasorbs) == 1) then
                  Chars(i:) = '2'
               else
                  Chars(i:) = '-'
               end if
            else
               if (IDet(i + ncasorbs) == 1) then
                  Chars(i:) = '+'
               else
                  Chars(i:) = '0'
               end if
            end if
         end do
      end function FMS_PrintDet

   end subroutine FMS_WriteFCIVec

!>
!!    Writes formatted information about the Mulliken charges to file.
!!    @ingroup output
!<
   subroutine FMS_WriteFCharge(T, file_name, first_time)
      use QM_MM_Module
      type(T_Trajectory) :: T
      character(len=*), intent(in) :: file_name
      logical, intent(in) :: first_time

      character(len=256) :: long_file_name
      integer(kind=DefInt) :: iunit
      logical :: file_existed
      integer(kind=DefInt) :: k

1     format(a10, 100a10)
2     format(f10.2, 100(1x, f9.5))

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (iUnit, 1) '#Time', (FMS_NumberedFileName('P', k), k=1, qcNumQM)
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

      ! 4. Write partial charges
      write (iunit, 2) T%get_time(), (T%Particle(k)%Charge, k=1, qcNumQM)
      close (iunit)

   end subroutine FMS_WriteFCharge

!>
!!
!!    Writes formatted coupling norm, dot product of coupling and
!!    velocity to file.
!!
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFCouple(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      logical :: file_existed
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files

      integer :: nstate, i, j

1     format(a10, 50(a15))
2     format(f10.2, 50(1x, f14.7))

      nstate = T%NumStates
      i = T%StateID

      iunit = units(T%TrajID)

      if (iunit == 0) then
         call FMS_OpenFile(FMS_NumberedFileName('Coup', T%TrajID), iunit, file_existed)
         if (.not. file_existed) then
            write (iUnit, 1) '#Time', (FMS_NumberedFileName('Coup', j), j=1, nstate), &
               (FMS_NumberedFileName(' C*V', j), j=1, nstate)
         end if
         units(T%TrajID) = iunit
      end if

      ! 4. Write the couplings
      write (IUnit, 2) T%get_time(), (sqrt(sum(FMS_Coupling(T, i, j)**2)), j=1, nstate), &
         (FMS_CoupDotVel(T, j), j=1, nstate)
      flush (IUnit)

   end subroutine FMS_WriteFCouple
!>
!!
!!    Writes formatted coupling norm, dot product of coupling and
!!    velocity to file.
!!
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFSOCouple(T)
      type(T_Trajectory), intent(in) :: T

      logical :: file_existed
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files

      real(kind=DefReal) :: Coup(2 * NSing + 6 * (T%NumStates - NSing))
      integer :: nstate, i, j, k

1     format(a10, 50(a15))
2     format(f10.2, 50(1x, f14.7))

      nstate = T%NumStates
      i = T%StateID

      iunit = units(T%TrajID)

      if (iunit == 0) then
         call FMS_OpenFile(FMS_NumberedFileName('SOCoup', T%TrajID), iunit, file_existed)
         if (.not. file_existed) then
            write (iUnit, 1) '#Time' !,
         end if
         units(T%TrajID) = iunit
      end if

      k = 0
      Coup = 0.d0
      ! 4. Write the couplings
      do j = 1, NSing
         k = k + 1
         Coup(k) = real(FMS_SOCoupling(T, i, j, T%Ms, 2))
         k = k + 1
         Coup(k) = aimag(FMS_SOCoupling(T, i, j, T%Ms, 2))
      end do

      if (NTrip /= 0) then
         do j = nsing + 1, nstate
            k = k + 1
            Coup(k) = real(FMS_SOCoupling(T, i, j, T%Ms, 1))
            k = k + 1
            Coup(k) = aimag(FMS_SOCoupling(T, i, j, T%Ms, 1))
            k = k + 1
            Coup(k) = real(FMS_SOCoupling(T, i, j, T%Ms, 2))
            k = k + 1
            Coup(k) = aimag(FMS_SOCoupling(T, i, j, T%Ms, 2))
            k = k + 1
            Coup(k) = real(FMS_SOCoupling(T, i, j, T%Ms, 3))
            k = k + 1
            Coup(k) = aimag(FMS_SOCoupling(T, i, j, T%Ms, 3))
         end do
      end if

      write (IUnit, 2) T%get_time(), Coup(1:k)
      flush (IUnit)

   end subroutine FMS_WriteFSOCouple
!>
!!    Writes out specified angles.
!!
!!    Columns are labeled with atom numbers.
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFDihedral(T, file_name, NDihedrals, iDihedral, first_time)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      character(len=256), intent(in) :: file_name
      integer(kind=DefInt), intent(in) :: NDihedrals, IDihedral(:, :)
      logical, intent(in) :: first_time

      integer(kind=DefInt), parameter :: MaxAngle = 20
      real(kind=DefReal), parameter :: rad2deg = 180.d0 / 3.141592654d0
      integer(kind=DefInt) :: i

      integer :: iunit
      character(len=256) :: long_file_name
      logical :: file_existed

      !         time             angle
1     format(a15, 20(2x, i3, 3(",", i3))) ! header format
2     format(3x, f12.2, 20(1x, f16.6)) ! angle  format

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      open (newunit=iunit, file=trim(long_file_name), position='append')
      if (.not. file_existed) write (iunit, 1) '# Time', (IDihedral(1:4, i), i=1, NDihedrals)

      write (iunit, 2) T%get_time(), &
         (rad2deg * FMS_Dihedral(T, &
                                 iDihedral(1, i), iDihedral(2, i), iDihedral(3, i), iDihedral(4, i)), i=1, NDihedrals)

      close (IUnit)
   end subroutine FMS_WriteFDihedral
!>
!!    Writes out specified angles.
!!
!!    Columns are labeled with atom numbers.
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFDistance(T, file_name, NBonds, iBond, first_time)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      character(len=256), intent(in) :: file_name
      integer(kind=DefInt), intent(in) :: NBonds, iBond(:, :)
      logical, optional :: first_time

      integer(kind=DefInt) :: i

      integer :: iunit
      character(len=256) :: long_file_name
      logical :: file_existed

      !         time             angle
1     format(a15, 20(3x, i3, ",", i3)) ! header format
2     format(3x, f12.2, 20(1x, f9.4)) ! bond   format

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (present(first_time) .and. first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (iunit, 1) '# Time', (iBond(1:2, i), i=1, NBonds)
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

      ! 4.  write output
      write (iunit, 2) T%get_time(), (FMS_Distance(T, iBond(1, i), iBond(2, i)), i=1, NBonds)

      close (IUnit)

   end subroutine FMS_WriteFDistance
!>
!!    Writes dipole of each state to file
!!    @ingroup output
!<
   subroutine FMS_WriteFDipole(T, file_name, first_time)
      use ElecStrucModule
      type(T_Trajectory) :: T
      character(len=*), intent(in) :: file_name
      logical, intent(in) :: first_time

      character(len=256) :: long_file_name
      integer(kind=DefInt) :: iunit
      logical :: file_existed

      integer :: nstate
      integer(kind=DefInt) :: k

1     format(a10, 50(a10))
2     format(f10.2, 1000f11.6)

      nstate = T%NumStates

      long_file_name = trim(FMSWorkingDir)//file_name
      file_existed = .true.
      if (first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (iUnit, 1) '#Time', &
            (FMS_NumberedFileName('Mag', k), k=1, nstate), &
            ([FMS_NumberedFileName('D', k)//'x', &
              FMS_NumberedFileName('D', k)//'y', &
              FMS_NumberedFileName('D', k)//'z'], k=1, nstate)
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

      ! 4. Write dipoles,
      ! This was writing the square of the dipole, was that
      ! intentional?
      write (IUnit, 2) T%get_time(), &
         (sqrt(sum(FMS_Dipole(T, k)**2)), k=1, T%NumStates), &
         (FMS_Dipole(T, k), k=1, T%NumStates)

      close (IUnit)

   end subroutine FMS_WriteFDipole

!>
!!    Writes formatted QM/MM energy partitioning
!!
!!    @ingroup output
!<
   subroutine FMS_WriteFMM(T, file_name, first_time)
      use QM_MM_Module
      use ElecStrucModule
      type(T_Trajectory), intent(in) :: T
      character(len=256), intent(in) :: file_name
      logical, optional :: first_time

      integer :: iunit
      character(len=256) :: long_file_name
      logical :: file_existed

      integer(kind=DefInt) :: iPtcle

      real(kind=DefReal) :: MMKin, MMPot, MMTot
      real(kind=DefReal) :: QMKin, QMPot, QMTot

1     format(a10, 4000a12)
2     format(f10.2, 4000f12.6)

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (present(first_time) .and. first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (iunit, 1) '#Time', 'QMPot', 'QMKin', 'QMTotal', 'MMPot', 'MMKin', 'MMTotal'
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

      ! 4.  write output
      MMKin = 0.d0
      do iPtcle = qcNumQM + 1, qcNumMM
         MMKin = MMKin + sum(T%Particle(iPtcle)%get_mom()**2) / (2.d0 * T%Particle(iPtcle)%Mass)
      end do
      MMPot = FMS_MMPot(T)
      MMTot = MMKin + MMPot

      QMKin = 0.d0
      do iPtcle = 1, qcNumQM
         QMKin = QMKin + sum(T%Particle(iPtcle)%get_mom()**2) / (2.d0 * T%Particle(iPtcle)%Mass)
      end do
      QMPot = FMS_PotentialT(T) - MMPot
      QMTot = QMKin + QMPot

      write (IUnit, 2) T%get_time(), QMPot, QMKin, QMTot, MMPot, MMKin, MMTot
      close (IUnit)

   end subroutine FMS_WriteFMM
!>
!!    Writes out specified pyramidalization angles.
!!    Columns are labeled with atom numbers.
!!    @ingroup output
!<
   subroutine FMS_WriteFPyram(T, file_name, NPyrams, IPyram, first_time)
      type(T_Trajectory), intent(in) :: T
      character(len=256), intent(in) :: file_name
      logical, optional :: first_time

      integer :: iunit
      character(len=256) :: long_file_name
      logical :: file_existed

      integer(kind=DefInt), parameter :: MaxGeometry = 100
      real(kind=DefReal), parameter :: rad2deg = 180.d0 / 3.141592654d0
      integer(kind=DefInt) :: IPyram(4, MaxGeometry), NPyrams

      integer(kind=DefInt) :: IxPyram, IParticle
      real(kind=DefReal) :: DAngle(NPyrams)
      real(kind=DefReal) :: cctemp(4, 3)
      character(len=256), allocatable :: CHeader(:)
      logical :: ZErr

1     format(A12, 49a14)
2     format(f12.2, 49f14.6)

      long_file_name = trim(FMSWorkingDir)//file_name

      file_existed = .true.
      if (present(first_time) .and. first_time) then
         inquire (file=long_file_name, exist=file_existed)
      end if

!     create column labels
      allocate (CHeader(NPyrams))
      do ixpyram = 1, npyrams
         CHeader(ixpyram) = ''
         write (CHeader(ixpyram), '(I4)') IPyram(1, ixpyram)
         CHeader(ixpyram) = FMS_NumberedFileName(CHeader(ixpyram), IPyram(2, ixpyram))
         CHeader(ixpyram) = FMS_NumberedFileName(CHeader(ixpyram), IPyram(3, ixpyram))
         CHeader(ixpyram) = FMS_NumberedFileName(CHeader(ixpyram), IPyram(4, ixpyram))
      end do

      if (.not. file_existed) then
         open (newunit=iunit, file=trim(long_file_name))
         write (IUnit, 1) 'Time', (CHeader(ixpyram), ixpyram=1, NPyrams)
!        Make sure pyrams are reasonable
         ZErr = .false.
         do Ixpyram = 1, NPyrams
            do iparticle = 1, 4
               ZErr = ZErr .or. (IPyram(iparticle, ixpyram) <= 0) .or. (IPyram(iparticle, ixpyram) > T%NumParticles)
            end do
            if (ZErr) then
               write (fmiOut, *) 'Invalid pyram requested: ', ixpyram
               stop
            end if
         end do
      else
         open (newunit=iunit, file=trim(long_file_name), position='append')
      end if

!     calculate pyram angles
      do ixpyram = 1, npyrams
         do IParticle = 1, 4
            cctemp(IParticle, 1:3) = T%Particle(IPyram(IParticle, ixpyram))%get_pos()
         end do
         call FMS_CalcPyrAngle(4, cctemp, DAngle(ixpyram), 1, 2, 3, 4)
      end do

      write (IUnit, 2) T%get_time(), (rad2deg * DAngle(ixpyram), ixpyram=1, NPyrams)
      close (IUnit)

   end subroutine FMS_WriteFPyram
!>
!!    Writes transition dipole between state one and all other states
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFTDipole(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      logical :: file_existed
      integer :: nstate, j

1     format(a10, 50(a10))
2     format(f10.2, 50(1x, f9.5))

      nstate = T%NumStates
      iunit = units(T%TrajID)

      if (iunit == 0) then
         call FMS_OpenFile(FMS_NumberedFileName('TDip', T%TrajID), iunit, file_existed)
         if (.not. file_existed) then
            write (iUnit, 1) '#Time', &
               (FMS_NumberedFileName('Mag', j), j=2, nstate), &
               ([FMS_NumberedFileName('TD', j)//'x', &
                 FMS_NumberedFileName('TD', j)//'y', &
                 FMS_NumberedFileName('TD', j)//'z'], j=2, nstate)
         end if
         units(T%TrajID) = iunit
      end if

      ! 4. Write transition dipole
      write (iunit, 2) T%get_time(), (sqrt(sum(FMS_TransDipole(T, j)**2)), j=2, nstate), &
         (FMS_TransDipole(T, j), j=2, nstate)
      flush (iunit)

   end subroutine FMS_WriteFTDipole
!>
!!    Writes formatted potential energy data
!!
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFPotEn(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use QM_MM_Module
      use ElecStrucModule
      type(T_Trajectory), intent(in) :: T

      character(len=256) :: line
      logical :: file_existed
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files

      real(kind=DefReal) :: time
      integer(kind=DefInt) :: nstate, i, j

      time = T%get_time()
      nstate = T%NumStates
      i = T%StateID

1     format(21(a10))
2     format(f10.2, 20(1x, f12.6))

      iunit = units(T%TrajID)

      if (iunit == 0) then

         call FMS_OpenFile(FMS_NumberedFileName('PotEn', T%TrajID), iunit, file_existed)

         if (.not. file_existed) then
            write (line, 1) '#Time', (FMS_NumberedFileName('PotEn', j), j=1, nstate), 'EClass'
            if (qczQMMM .and. eszMMForce) line = trim(line)//'    EClass'
            write (iunit, '(a)') trim(line)
         end if
         units(T%TrajID) = iunit
      end if

      ! 3. write the energies
      if (qczQMMM .and. eszMMForce) then
         write (iunit, 2) time, (Potential(T, j), j=1, nstate), Kinetic(T) + Potential(T, i), FMS_MMPot(T)
      else
         write (iunit, 2) time, (Potential(T, j), j=1, nstate), Kinetic(T) + Potential(T, i)
      end if
      flush (iunit)

   end subroutine FMS_WriteFPotEn
!>
!!    Writes Force (in XYZ file format)
!!    \param comment   Comment to write on comment line
!!    @ingroup output
!<
   subroutine FMS_WriteFForces(T1, filename, comment, firsttime)
      type(T_Trajectory), intent(in) :: T1
      character(len=*), intent(in) :: filename, comment
      character(len=256) :: WriteComment
      logical, optional :: firsttime

      integer(kind=DefInt) :: IUnit, IPtcle, IStart, IState
      character(len=256) :: file_name, pos_str

      file_name = trim(FMSWorkingDir)//filename

      IState = T1%StateID
      if (len(trim(comment)) == 0) then
         write (WriteComment, *) 'Time ', T1%get_time()
      else
         WriteComment = trim(Comment)
      end if

      pos_str = 'append'
      if (present(firsttime)) then
         if (firsttime) pos_str = 'rewind'
      end if

      open (newunit=IUnit, file=file_name, position=pos_str)

      write (IUnit, *) T1%NumParticles
      write (IUnit, *) trim(WriteComment)

      do IPTcle = 1, T1%NumParticles
         IStart = (IPtcle - 1) * T1%Particle(IPtcle)%NumDimensions
         write (IUnit, '(A3, 3F17.9)') T1%Particle(IPtcle)%Elmnt, T1%ElecStruc%DerivMat(IState, IState, IStart + 1:Istart + 3)
      end do
      close (IUnit)

   end subroutine FMS_WriteFForces

!>
!!    Writes formatted potential energy data
!!
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFQMRR(T)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T

      character(len=256) :: line
      logical :: file_existed
      integer(DefInt) :: iunit
      integer(DefInt), save :: units(MaxTrajLimit) = 0 ! this keeps track of the
      ! unit number of the files

      real(kind=DefReal) :: time
      integer(kind=DefInt) :: nstate, i, j

      time = T%get_time()
      nstate = T%NumStates
      i = T%StateID

1     format(a10, 21(a13))
2     format(f10.2, 20(1x, f12.6))

      iunit = units(T%TrajID)
      if (iunit == 0) then

         call FMS_OpenFile(FMS_NumberedFileName('QMRR', T%TrajID), iunit, file_existed)

         if (.not. file_existed) then
            write (line, 1) '#Time', (FMS_NumberedFileName('QMRR', j), j=1, nstate) !,
!     &         'EClass'
!            if( qczQMMM .and. eszMMForce) line = line//'    EClass'
            write (iunit, '(a)') trim(line)
         end if
         units(T%TrajID) = iunit
      end if

      ! 3. write the energies
!      if( qczQMMM .and. eszMMForce )then
!         write(iunit,2) time,
!     &   ( Potential(T,j), j=1,nstate),
!     &   Kinetic(T) + Potential(T,i),
!     &   FMS_MMPot(T)
!      else
      write (iunit, 2) time, (T%ElecStruc%QMRR(j), j=1, nstate) !,
!     &   Kinetic(T) + Potential(T,i)
!      endif

      flush (iunit)
   end subroutine FMS_WriteFQMRR

end module TrajectoryIOModule
