!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    @brief IO Methods for the "Bundle" datatype.
!<
module BundleIOModule
   use GlobalModule, only: DefInt, DefReal, DefComp, FMSWorkingDir, FMS_DieError, &
                           fmzwriteeverystep, fmzndatfile, fmzbundmatfile, fmzcorrfile, fmzedatfile, fmzalltext, &
                           fmnsteptoprint, NTrip, fmzPCOlap, fmzSOME, fmzSOCeff, FMS_PrintCMat

   use TrajectoryModule
   use TrajectoryCalcsModule, only: FMS_BundleUpdated
   use TrajectoryIOModule, only: FMS_WriteTrajFiles
   use BundleModule
   use BundleCalcsModule, only: FMS_KineticB, FMS_PotentialB, FMS_bH, FMS_bSDot, FMS_Branching, FMS_Norm
   use OverlapModule, only: overlap, overlap_V
   use SMDModule, only: smSMD

   implicit none
   private

   public :: FMS_Output

contains
!>
!!    Writes Formatted Branching Ratios
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFBranching(B1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle) :: B1
      integer(kind=DefInt) :: nstate, i
      real(kind=DefReal) :: Prob(B1%NumStates)

      integer, save :: iunit = 0
      character(len=256) :: file_name, fmt_str
      logical :: file_exists

      nstate = B1%NumStates

!     write header the first time through
      if (iunit == 0) then
         ! 1. Make the file name
         file_name = trim(FMSWorkingDir)//'N.dat'

         ! 2. See if it exists
         inquire (file=file_name, exist=file_exists)
         if (file_exists) then
            open (newunit=iunit, file=trim(file_name), position='append')
         else
            open (newunit=iunit, file=trim(file_name))
            write (fmt_str, *) '(a12,', nstate, '(" State ",i3,2x ),a10)'
            write (iunit, fmt_str) '#Time', (i, i=1, nstate), 'Norm'
         end if
      end if

      call FMS_Branching(B1, Prob)

      write (IUnit, '(1x,f11.2,50(1x,f11.9))') B1%CurrentTime, (Prob(i), i=1, nstate), FMS_Norm(B1)

   end subroutine FMS_WriteFBranching

!>
!!    Writes formatted potential, kinetic, and total energies.
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_WriteFEnergy(B1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle) :: B1

      real(kind=DefReal) :: PotQM, KinQM, PotCl, KinCl

      integer, save :: iunit = 0
      character(len=256) :: file_name
      logical :: file_exists

!     write header the first time through
      if (iunit == 0) then
         ! 1. Make the file name
         file_name = trim(FMSWorkingDir)//'E.dat'

         ! 2. See if it exists
         inquire (file=file_name, exist=file_exists)
         if (file_exists) then
            open (newunit=iunit, file=trim(file_name), position='append')
         else
            open (newunit=iunit, file=trim(file_name))
            write (iunit, '(a10,6(a15))') '#Time', 'EPotentialQM', 'EKineticQM', 'ETotalQM', &
               'EPotentialCl', 'EKineticCl', 'ETotalCl'
         end if
      end if

      ! calculate the quantities
      PotQM = FMS_PotentialB(B1, Incoherent=.false.)
      KinQM = FMS_KineticB(B1, Incoherent=.false.)

      PotCl = FMS_PotentialB(B1, Incoherent=.true.)
      KinCl = FMS_KineticB(B1, Incoherent=.true.)

      ! write energies
      write (IUnit, '(f10.2,6(1x,f14.9))') B1%CurrentTime, PotQM, KinQM, PotQM + KinQM, PotCl, KinCl, PotCl + KinCl

   end subroutine FMS_WriteFEnergy

!>
!!    Write bundle matrices out to file
!!    \param FileName Name of file to output
!!          Valid values are 'H.dat', 'S.dat' and 'Sdot.dat'
!!    @ingroup output
!<
   subroutine FMS_WriteFMatrix(B1, FileName)
      type(T_TrajectoryBundle), intent(in) :: B1
      character(len=*), intent(in) :: FileName

      integer(kind=DefInt) :: IUnit, ITraj, JTraj
      character(len=:), allocatable :: FilePath
      complex(kind=DefComp) :: xMatTemp(B1%NumTraj, B1%NumTraj)

      FilePath = trim(FMSWorkingDir)//trim(FileName)
      open (newunit=IUnit, file=FilePath, position='append')

      write (IUnit, *) B1%CurrentTime

      do ITraj = 1, B1%NumTraj
         do JTraj = 1, B1%NumTraj
            select case (FileName)
            case ('H.dat')
               xMatTemp(ITraj, JTraj) = FMS_bH(B1, ITraj, JTraj)
            case ('S.dat')
!Want full overlap matrix including overlap between different states:
               xMatTemp(ITraj, JTraj) = overlap(B1%Trajectory(ITraj), B1%Trajectory(JTraj))
!               xMatTemp(ITraj,JTraj)=FMS_bS(B1,ITraj,JTraj)
            case ('SDot.dat')
               xMatTemp(ITraj, JTraj) = FMS_bSDot(B1, ITraj, JTraj)
            case default
               call FMS_DieError('Unknown matrix in WriteFMatrix')
            end select
         end do
      end do

      call FMS_PrintCMat(xMatTemp, IUnit)

      close (IUnit)

   end subroutine FMS_WriteFMatrix

!>
!!    Writes formatted autocorrelation function \f$ <\Psi(0)|\Psi(t)> \f$
!!    @ingroup output
!<
   subroutine FMS_WriteFCorr(B1)
      type(T_TrajectoryBundle), intent(in) :: B1
      type(T_TrajectoryBundle), save :: B_init

      complex(kind=DefComp) :: Corr
      logical, save :: initialize = .true.

      integer, save :: iunit = 0
      character(len=256) :: file_name
      logical :: file_exists

      if (initialize) then
         B_init = B1
         initialize = .false.
      end if

!     write header the first time through
      if (iunit == 0) then
         ! 1. Make the file name
         file_name = trim(FMSWorkingDir)//'CFxn.dat'

         ! 2. See if it exists
         inquire (file=file_name, exist=file_exists)
         if (file_exists) then
            open (newunit=iunit, file=trim(file_name), position='append')
         else
            open (newunit=iunit, file=trim(file_name))
            write (iunit, '(6(a10))') '#Time', 'Amp^2', 'CorrReal', 'CorrImag'
         end if
      end if

      Corr = overlap_bundle(B1, B_init)
      write (iUnit, '(f10.2,4(1x,f12.9))') B1%CurrentTime, abs(Corr)**2, Corr

   contains

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function overlap_bundle(B1, B2) result(S)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         type(T_TrajectoryBundle), intent(in) :: B1, B2
         complex(kind=DefComp) :: S
         integer(kind=defInt) :: i, j

         S = (0.d0, 0.d0)

         do i = 1, B1%NumTraj
            do j = 1, B2%NumTraj
               S = S + overlap(B1%Trajectory(i), B2%Trajectory(j))
            end do
         end do

      end function overlap_bundle

   end subroutine FMS_WriteFCorr
!>
!!    Writes formatted parent child overlap
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! - - - - - - - -
   subroutine FMS_WriteFParentChildOlap(T1, T2, filename, firsttime)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Trajectory), intent(in) :: T1, T2
      character(len=*), intent(in) :: filename
      logical, optional, intent(in) :: firsttime
      complex(defComp) :: olap
      logical :: file_exists, open_new_file

      integer :: IUnit
      character(len=256) :: file_name

      ! 1. get the file name and unit number
      file_name = trim(FMSWorkingDir)//'pcolap.'//filename

      ! 2. work out if we are opening a new file
      open_new_file = .false.
      if (present(firsttime)) then
         if (firsttime) then
            inquire (file=file_name, exist=file_exists)
            open_new_file = .not. file_exists
         end if
      end if

      ! 3. open the file
      if (open_new_file) then
         open (newunit=IUnit, file=trim(file_name))
         if (NTrip /= 0) then
            write (IUnit, '(a15,2(2x,i1), a15,2(2x,i1))') '# Child ID, Ms', T1%TrajID, T1%Ms, &
               'Parent ID, Ms', T2%TrajID, T2%Ms
         else
            write (IUnit, '(a10,i1,a9,i1)') '# Child ID', T1%TrajID, 'Parent ID', T2%TrajID
         end if
         write (IUnit, '(4(4x,a10))') '# Time', 'abs(S)', 'Re(S)', 'Im(S)'
      else
         open (newunit=iunit, file=trim(file_name), position='append')
      end if

      ! 4.  write output
      olap = overlap(T1, T2)
      write (IUnit, '(1x,f11.2,3(1x,f11.9))') T1%get_time(), abs(olap), real(olap), aimag(olap)

      close (IUnit)

   end subroutine FMS_WriteFParentChildOlap

!>
!!    Writes SOME for centroid position
!!    @ingroup output
!<
   subroutine FMS_WriteFSOME(T1, T2, Tc, filename, firsttime)
      type(T_Trajectory), intent(in) :: T1, T2, Tc
      character(len=*), intent(in) :: filename
      logical, optional, intent(in) :: firsttime
      logical :: file_exists

      integer :: IUnit
      character(len=256) :: file_name

      file_name = trim(FMSWorkingDir)//'SOME.'//filename

      ! 2. work out if we are opening a new file
      file_exists = .true.
      if (present(firsttime)) then
         if (firsttime) then
            inquire (file=file_name, exist=file_exists)
         end if
      end if

      ! 3. open the file
      if (.not. file_exists) then
         open (newunit=IUnit, file=trim(file_name))
         write (IUnit, '(a15,2(2x,i1), a15,2(2x,i1))') '# Child ID, Ms', T1%TrajID, T1%Ms, &
            'Parent ID, Ms', T2%TrajID, T2%Ms
         write (IUnit, '(4(4x,a10))') '# Time', 'Centroid pos.', 'SOME (Re, Im)', 'SOMENorm'
      else
         open (newunit=iunit, file=trim(file_name), position='append')
      end if

      ! 4.  write output
      write (IUnit, '(1x,f11.2,5(1x,f12.9))') T1%get_time(), Tc%Particle(1)%get_pos(1), &
         overlap_V(T1, T2, Tc), abs(overlap_V(T1, T2, Tc))

      close (IUnit)

   end subroutine FMS_WriteFSOME

!>
!!    Writes SOCeff
!!    @ingroup output
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! - - - - - - - -
   subroutine FMS_WriteFSOCeff(T1, T2, js, filename, firsttime)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! - - - - - - - -
      use SpawnModule, only: spawn_couple
      type(T_Trajectory), intent(in) :: T1, T2
      character(len=*), intent(in) :: filename
      logical, optional :: firsttime
      logical :: file_exists, open_new_file

      integer :: IUnit, js
      character(len=256) :: file_name

      ! 1. get the file name and unit number
      file_name = trim(FMSWorkingDir)//'SOCeff.'//filename

      ! 2. work out if we are opening a new file
      open_new_file = .false.
      if (present(firsttime)) then
         if (firsttime) then
            inquire (file=file_name, exist=file_exists)
            open_new_file = .not. file_exists
         end if
      end if

      ! 3. open the file
      if (open_new_file) then
         open (newunit=IUnit, file=trim(file_name))
         write (IUnit, '(a15,2(2x,i1), a15,2(2x,i1))') '# Child ID, Ms', T1%TrajID, T1%Ms, &
            'Parent ID, Ms', T2%TrajID, T2%Ms
      else
         open (newunit=iunit, file=trim(file_name), position='append')
      end if

      ! 4.  write output
      write (IUnit, '(1x,f11.2,2(1x,f15.9))') T1%get_time(), T1%Particle(1)%get_pos(1), spawn_couple(T1, js)

      close (IUnit)

   end subroutine FMS_WriteFSOCeff

!>
!!    Read H, S and SDot in from a file. Useful for analysis mode.
!!    NOTE: Currently not used anywhere
!<
   subroutine FMS_ReadFBundle(B1, IUnitH, iUnitS, iUnitSDot)
      type(T_TrajectoryBundle) :: B1
      integer(kind=DefInt) :: IUnitH, iUnitS, iUnitSDot, iTraj
      real(kind=DefReal) :: TimeTemp

!     Read in bundle matrices
      read (IUnitH, *) timetemp
      if (timetemp /= B1%CurrentTime) then
         call FMS_DieError('Time mismatch in H.dat')
      end if
      call FMS_ReadCMat(B1%BMatrices%H, IUnitH)

      read (IUnitS, *) timetemp
      if (timetemp /= B1%CurrentTime) then
         call FMS_DieError('Time mismatch in S.dat')
      end if
      call FMS_ReadCMat(B1%BMatrices%S, IUnitS)

      read (IUnitSDot, *) timetemp
      if (timetemp /= B1%CurrentTime) then
         call FMS_DieError('Time mismatch in SDot.dat')
      end if
      call FMS_ReadCMat(B1%BMatrices%SDot, IUnitSDot)

      do iTraj = 1, B1%NumTraj
         call FMS_BundleUpdated(B1%Trajectory(ITraj))
      end do

   contains

!>
!!    For reading complex matrix files
!<
      subroutine FMS_ReadCMat(A, Unit)
         complex(kind=DefComp), intent(inout) :: A(:, :)
         integer(kind=DefInt), intent(in) :: Unit
         integer(kind=DefInt) :: i, j

         do i = 1, size(A, dim=1)
            read (Unit, '(10f9.4)') (A(i, j), j=1, size(A, dim=2))
         end do
      end subroutine FMS_ReadCMat

   end subroutine FMS_ReadFBundle

!>
!!    Driver for all output during propagation
!!
!!    @ingroup output
!<
   subroutine FMS_Output(B1, FirstTime)
      use SMDModule, only: SMD_Print
      type(t_trajectoryBundle), intent(in) :: B1
      logical, optional :: FirstTime

      integer(kind=DefInt) :: ITraj, JTraj, CBFi, CBFj
      integer(kind=DefInt), save :: NPreviousTraj, nPrintStep
      character(len=32) :: suffix
      logical :: zWriteHeader

!
!     Initialization and output cleaning
!
      zWriteHeader = .false.
      if (present(FirstTime)) zWriteHeader = FirstTime

!
!     Leave if we are not writing output this timestep
!
      if (.not. fmzWriteEveryStep) then
         nPrintStep = nPrintStep + 1
         if (mod(nPrintStep, fmnStepToPrint) /= 0) return
      end if

      if (fmzCorrFile .or. fmzAllText) call FMS_WriteFCorr(B1) ! correlation function
      if (fmzEDatFile .or. fmzAllText) call FMS_WriteFEnergy(B1) ! bundle energies
      if (fmzNDatFile .or. fmzAllText) call FMS_WriteFBranching(B1) ! population

!     Write bundle matrices (alive trajectories only)
      if (fmzBundMatFile .or. fmzAllText) then
         call FMS_WriteFMatrix(B1, 'H.dat')
         call FMS_WriteFMatrix(B1, 'S.dat')
         call FMS_WriteFMatrix(B1, 'SDot.dat')
      end if

!
!     Write trajectory-level files
!
      do iTraj = 1, B1%NumTraj
         if (B1%Trajectory(iTraj)%TrajID <= nPreviousTraj) then
            zWriteHeader = .false.
         else
            nPreviousTraj = B1%Trajectory(iTraj)%TrajID
            zWriteHeader = .true.
         end if
         write (suffix, '(i32)') B1%Trajectory(iTraj)%TrajID
         call FMS_WriteTrajFiles(B1%Trajectory(iTraj), trim(adjustl(suffix)), zFirst=zWriteHeader)
![ydl
!     Write trajectory-pair-level files
         if ((fmzPCOlap) .and. (B1%Numtraj > 1)) then
            do jTraj = 1, B1%NumTraj
               if (iTraj == jTraj) cycle
               if (B1%Trajectory(iTraj)%ParentID == B1%Trajectory(jTraj)%TrajID) then
                  if (NTrip /= 0) then
                     write (suffix, '(i32)') ITraj
                  end if
                  call FMS_WriteFParentChildOlap(B1%Trajectory(iTraj), B1%Trajectory(jTraj), &
                                                 trim(adjustl(suffix)), zWriteHeader)
               end if
            end do
         end if
!ydl]
![vcb
         if ((fmzSOME) .and. (B1%Numtraj > 1)) then
            do jTraj = 1, B1%NumTraj
               if (iTraj == jTraj) cycle
               if (B1%Trajectory(iTraj)%ParentID == B1%Trajectory(jTraj)%TrajID) then
                  if (NTrip /= 0) then
                     write (suffix, '(i32)') ITraj
                  end if
                  CBFi = B1%Trajectory(iTraj)%CBF
                  CBFj = B1%Trajectory(jTraj)%CBF
                  call FMS_WriteFSOME(B1%Trajectory(iTraj), &
                                      B1%Trajectory(jTraj), &
                                      B1%Centroids(((CBFi - 2) * (CBFi - 1)) / 2 + CBFj), &
                                      trim(adjustl(suffix)), zWriteHeader)
               end if
            end do
         end if

         if ((fmzSOCeff) .and. (B1%Numtraj > 1)) then
            do jTraj = 1, B1%NumTraj
               if (iTraj == jTraj) cycle
               if (B1%Trajectory(iTraj)%ParentID == B1%Trajectory(jTraj)%TrajID .and. &
                   B1%Trajectory(iTraj)%Ms == 2) then
                  call FMS_WriteFSOCeff(B1%Trajectory(iTraj), B1%Trajectory(jTraj), &
                                        B1%Trajectory(jTraj)%StateID, &
                                        trim(adjustl(suffix)), zWriteHeader)
               end if
            end do
         end if
!vcb]
      end do

!     Write SMD output if requested
      if (smSMD) then
         call SMD_Print(B1%CurrentTime, zWriteHeader)
      end if

!     Save current number of trajectories
      NPreviousTraj = B1%NumTraj + B1%NumDeadTraj

   end subroutine FMS_Output

end module BundleIOModule
