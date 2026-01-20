! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module RestartModule
   use GlobalModule
   use ParticleModule
   use TrajectoryModule
   use BundleModule
   implicit none
   private

   public :: GetRestart, PutRestart
   integer(kind=DefInt), public :: inIRestart, & !> Determines whether this is a restart (1 = restart, 0 = not restart)
                                   RestartStep ! Number of steps between archived restarts
   integer(kind=DefReal), public :: inIRestartTraj(maxtrajlimit)
   real(kind=DefReal), public :: RestartTime
   logical, public :: zRedoRestartES !> Redo the electronic structure calculation for restarted step?
contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine GetRestart(B, time)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! if time < 0 then it opens Last_Bundle.txt
! if time > 0 then open Checkpoint.txt and try to find bundle with the
! time specfied.
      type(T_TrajectoryBundle), intent(inout) :: B
      real(kind=DefReal), intent(in) :: time

      integer(DefInt) :: nfile
      character(len=200) :: long_name, short_name, last_line
      logical :: file_exists

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Setup the file names
      if (time > 0) then
         short_name = "Checkpoint.txt"
         write (fmiOut, '(a,f10.2)') "Searching for time=", time, " in Checkpoint.txt"
      else
         short_name = "Last_Bundle.txt"
         write (fmiOut, *) "Taking restart from Last_Bundle.txt"
      end if
      long_name = trim(FMSWorkingDir)//trim(short_name)
!temp_name = trim(FMSWorkingDir)//".tmp."//trim(short_name)

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Check that the file exists
      inquire (file=long_name, exist=file_exists)
      if (.not. file_exists) then
         call FMS_DieError("Restart file '"//trim(long_name)//"' not found")
      end if

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Get the restart Bundle

      if (time > 0) then
         open (newunit=nfile, file=long_name, status='old')
         do
            call ReadBundle(B, nfile)

            if (B%numtraj == 0) then
               call FMS_DieError("GetRestart: end of Checkpoint.txt reached, requested time not found")
            end if

            if (abs(B%CurrentTime - time) < 0.25) exit

         end do
! attempt to avoid the use of unix commands
! this is done by reading & writing the last line of the current time step
         backspace (nfile)
         read (nfile, '(A200)') last_line
         backspace (nfile)
         write (nfile, '(A)') trim(last_line)
         close (unit=nfile)
      else
         open (newunit=nfile, file=long_name, action='read', status='old')
         write (fmiOut, '(a)') "Reading from: "//trim(long_name)
         call B%destroy()
         call ReadBundle(B, nfile)

         if (B%NumTraj == 0) then
            call FMS_DieError("could not read Last_Bundle.txt")
         end if

         close (unit=nfile)
      end if

!select which trajectories to restart from
      call prune_bundle(B)
   end subroutine GetRestart

   subroutine prune_bundle(B)
      type(T_TrajectoryBundle), intent(inout) :: B

      type(T_TrajectoryBundle) :: Btemp
      integer(kind=DefInt) :: ntraj, nlivetraj, ndeadtraj, itraj, jtraj, itrajnew, jtrajnew, iCent
      integer(kind=DefInt) :: CBFi, CBFj
      character(len=4*maxtrajlimit) :: ctrajrestart

      if (inirestarttraj(1) == 0) return !nothing to do

!figure out how many trajectories user requested
      ctrajrestart = ''
      ntraj = 0
      do itraj = 1, maxtrajlimit
         if (inirestarttraj(itraj) /= 0) then
            ntraj = ntraj + 1
            write (ctrajrestart((ntraj - 1) * 4 + 1:ntraj * 4), '(I3,1x)') inirestarttraj(itraj)
         end if
      end do

      write (fmiOut, '(A11,I0,A14)') 'Restarting ', ntraj, ' trajectories:'
      write (fmiOut, '(A)') trim(ctrajrestart)

!figure out how many live trajectories we want to restart
      ctrajrestart = ''
      nlivetraj = 0
      do itraj = 1, B%numtraj
         if (any(inirestarttraj(:) == B%trajectory(itraj)%trajid)) then
            nlivetraj = nlivetraj + 1
            write (ctrajrestart((nlivetraj - 1) * 4 + 1:nlivetraj * 4), '(I3,1x)') B%trajectory(itraj)%trajid
         end if
      end do

      if (nlivetraj > 0) then
         write (fmiOut, '(A6,I0,A19)') 'Found ', nlivetraj, ' live trajectories:'
         write (fmiOut, '(A)') trim(ctrajrestart)
      end if

!figure out if we are restarting a dead trajectory
      ctrajrestart = ''
      ndeadtraj = 0
      do itraj = 1, B%numDeadTraj
         if (any(inIRestartTraj(:) == B%DeadTraj(ITraj)%TrajID)) then
            NdeadTraj = NdeadTraj + 1
            write (ctrajrestart((ndeadtraj - 1) * 4 + 1:ndeadtraj * 4), '(I3,1x)') B%DeadTraj(itraj)%trajid
         end if
      end do

      if (ndeadtraj > 0) then
         write (fmiOut, '(A6,I0,A19)') 'Found ', ndeadtraj, ' dead trajectories:'
         write (fmiOut, '(A)') trim(ctrajrestart)
      end if

      if ((NDeadTraj == 1 .and. NTraj > 1) .or. (NDeadTraj > 1)) then
         call FMS_DieError('Error in PruneBundle. Can only restart one dead trajectory at a time')
      end if

      if (NDeadTraj == 0 .and. NLiveTraj /= NTraj) then
         call FMS_DieError('Error in PruneBundle. Did not find all requested live trajectories')
      end if

      select case (NDeadTraj)
      case (0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!For GAIMS it was assumed that all nlivetraj trajectories are restarted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call Btemp%create(numtraj=nlivetraj, &
                           numdeadtraj=B%NumDeadTraj, &
                           numstates=B%NumStates, &
                           numparticles=B%NumParticles, &
                           ncbfs=B%NCBFs)
         itrajnew = 0
         do itraj = 1, B%numTraj
!    TODO(danielhollas): It looks like this code is currently broken from non-GAIMS simulations?
!    The following if statement should be uncommented!
!    if (any(inirestarttraj(:).eq.B%trajectory(itraj)%trajid)) then
            itrajnew = itrajnew + 1
            Btemp%Trajectory(itrajnew) = B%Trajectory(itraj)

!bfec
            if (glzCentroids) then
               jtrajnew = 0
               CBFi = B%Trajectory(itraj)%CBF
               do jtraj = 1, B%numTraj - 1 !centroid
!         if (any(inirestarttraj(:).eq.B%trajectory(jtraj)%trajid)) then
                  CBFj = B%Trajectory(jtraj)%CBF !not ideal way, but for now let's see
                  jtrajnew = jtrajnew + 1
                  iCent = ((CBFi - 2) * (CBFi - 1)) / 2 + CBFj
!           jCent=((iTrajnew-2)*(iTrajnew-1))/2+jTrajnew
                  Btemp%Centroids(iCent) = B%Centroids(iCent)
!         endif
               end do
!       do jtraj=1,B%numTraj-1 !centroid
!         if (any(inirestarttraj(:).eq.B%trajectory(jtraj)%trajid)) then
!           jtrajnew=jtrajnew+1
!           iCent=((iTraj-2)*(iTraj-1))/2+jTraj
!           jCent=((iTrajnew-2)*(iTrajnew-1))/2+jTrajnew
!           Btemp%Centroids(jCent)=B%Centroids(iCent)
!         endif
!       enddo
            end if
!    endif
         end do
         if (B%NumDeadTraj > 0) then !copy dead trajectories
            do iTraj = 1, B%NumDeadTraj
               Btemp%DeadTraj(iTraj) = B%DeadTraj(iTraj)
            end do
            Btemp%DeadH = B%DeadH
         end if
         Btemp%CurrentTime = B%CurrentTime

      case (1)
         ! Create bundle with a single active trajectory
         call Btemp%create(numtraj=1, &
                           numdeadtraj=0, &
                           numstates=B%NumStates, &
                           numparticles=B%NumParticles, &
                           ncbfs=1)

!Find the dead trajectory to restart from
         do itraj = 1, B%numDeadTraj
            if (inIRestartTraj(1) == B%DeadTraj(ITraj)%TrajID) Btemp%Trajectory(1) = B%DeadTraj(ITraj)
         end do
!Set the current time to the trajectory's deadtime
         Btemp%CurrentTime = Btemp%Trajectory(1)%DeadTime
!Reset ES flags, if not using restarted ES
         if (zRedoRestartES) then
            Btemp%Trajectory(1)%ESFlags%zPotEnCurrent = .false.
            Btemp%Trajectory(1)%ESFlags%zDerivCurrent = .false.
            Btemp%Trajectory(1)%ESFlags%zTransDipsCurrent = .false.
            Btemp%Trajectory(1)%ESFlags%zDipolesCurrent = .false.
         end if

      case default
         call FMS_DieError("Whoops, we should not be here!")

      end select

      B = Btemp
      call Btemp%destroy()

   end subroutine prune_bundle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Made some changes here to make sure that a usable copy of the
! "Last_Bundle.txt" file exists at all times. This was becoming a problem for
! large calculations where writing this file takes non-negligible amounts
! of time. -- EGH
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine PutRestart(B)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle), intent(in) :: B

      character(len=200) :: temp_file, last_file, check_file
      integer(Defint) :: nfile

      temp_file = trim(FMSWorkingDir)//"tmp.Last_Bundle.txt"
      last_file = trim(FMSWorkingDir)//"Last_Bundle.txt"
      check_file = trim(FMSWorkingDir)//"Checkpoint.txt"

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Write the temp file
      open (newunit=nfile, file=temp_file)
      call WriteBundle(B, nfile)
      close (unit=nfile)

! Kill old bundle
      call FMS_DeleteFile('Last_Bundle.txt')

! Copy temp bundle to new bundle
      call system("cp "//trim(temp_file)//" "//trim(last_file))

! Kill temp bundle
      call FMS_DeleteFile('tmp.Last_Bundle.txt')

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Work out if we archive this one or not
      if (mod(int(B%CurrentTime / gldTimeStep), RestartStep) == 0) then
         open (newunit=nfile, file=check_file, position='append')
         call WriteBundle(B, nfile)
         close (unit=nfile)
      end if

   end subroutine PutRestart

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         BUNDLE
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine WriteBundle(B, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! writes the Bundle out to a flat, plain text file
! this is used for restarting the Trajectory
      type(T_TrajectoryBundle), intent(in) :: B
      integer(kind=DefInt), optional :: nfile_in

      integer(kind=DefInt) :: nf ! unit number to write to

      integer(kind=DefInt) :: ntraj, ndead, npart, n, ncbfs

! set which file we are writing to
      nf = 6
      if (present(nfile_in)) nf = nfile_in

1     format(i6, 34x, a) ! integer   format with label
2     format(f14.2, 26x, a) ! time type format with label

      ntraj = B%NumTraj
      ndead = B%NumDeadTraj
      npart = B%NumParticles
      nCBFs = B%NCBFs

      write (nf, *) "############### START BUNDLE #################"
      write (nf, 2) B%CurrentTime, " / Current time"
      write (nf, 1) ntraj, " / Number of Live Trajectories"
      write (nf, 1) ndead, " / Number of Dead Trajectories"
      write (nf, 1) B%NumStates, " / Number of States"
      write (nf, 1) NSing, " / Number of Singlet States"
      write (nf, 1) npart, " / Number of Particles"
      write (nf, 1) NCBFs, " / Number of Contracted Basis Functions"

      write (nf, *) "############ Common Particle info ###########"
      do n = 1, npart
         call WriteParticle(B%Trajectory(1)%Particle(n), nf)
      end do

      do n = 1, ntraj
         write (nf, '(a,i4,a)') "############### Live trajectory ", n, " ############"
         call WriteTraj(B%Trajectory(n), nf)
      end do

!bfec
      if (glzCentroids) then
! do n = 1, ntraj*(ntraj-1)/2
         do n = 1, nCBFs * (nCBFs - 1) / 2
            write (nf, '(a,i4,a)') "###############   Centroid      ", n, " ############"
            call WriteTraj(B%Centroids(n), nf)
         end do
      end if

      do n = 1, ndead
         write (nf, '(a,i4,a)') "############### Dead trajectory ", n, " ############"
         call WriteTraj(B%DeadTraj(n), nf)
      end do

   end subroutine WriteBundle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine ReadBundle(B, nf)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle), target, intent(inout) :: B
      integer(kind=DefInt), intent(in) :: nf ! unit number to read from
      integer(kind=DefInt) :: ntraj, ndead, n, &
                              nstate, npart, ncbfs

      type(T_Trajectory), pointer :: T

      real(kind=DefReal) :: time

      read (nf, *, end=20, err=20) ! ######## START BUNDLE ###########
      read (nf, *, end=20, err=20) time
      read (nf, *, end=20, err=20) ntraj
      read (nf, *, end=20, err=20) ndead
      read (nf, *, end=20, err=20) nstate
      read (nf, *, end=20, err=20) nsing
      read (nf, *, end=20, err=20) npart
      read (nf, *, end=20, err=20) ncbfs

      call B%create(numtraj=ntraj, &
                    numdeadtraj=ndead, &
                    numstates=nstate, &
                    numparticles=npart, &
                    ncbfs=ncbfs)

      B%CurrentTime = time
      B%NCBFs = ncbfs

      read (nf, *, end=20, err=20) ! # Common particle info

      T => B%Trajectory(1)

      do n = 1, npart
         call ReadParticle(T%Particle(n), nf)
      end do

! copy this info to the rest of the trajectories in the bundle
      do n = 2, ntraj
         B%Trajectory(n) = T
      end do
      do n = 1, ndead
         B%DeadTraj(n) = T
      end do
!bfec
      if (glzCentroids) then
         do n = 1, ncbfs * (ncbfs - 1) / 2
            B%Centroids(n) = T
         end do
      end if

! Now read in the Trajectories
      do n = 1, ntraj
         read (nf, *)
         call ReadTraj(B%Trajectory(n), nf)
         call B%Trajectory(n)%set_time(time)
      end do

!bfec
      if (glzCentroids) then
! do n = 1, ntraj*(ntraj-1)/2
         do n = 1, nCBFs * (nCBFs - 1) / 2
            read (nf, *)
            call ReadTraj(B%Centroids(n), nf)
            call B%Centroids(n)%set_time(time)
         end do
      end if

      do n = 1, ndead
         read (nf, *)
         call ReadTraj(B%DeadTraj(n), nf)
         call B%DeadTraj(n)%set_time(time)
         B%DeadTraj(n)%ESFlags%zPotEnCurrent = .true.
         B%DeadTraj(n)%ESFlags%zTransDipsCurrent = .true.
         B%DeadTraj(n)%ESFlags%zDipolesCurrent = .true.

      end do

      return

20    continue
      write (fmiout, *) "error reading Bundle"
      call B%destroy()
   end subroutine ReadBundle

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!           TRAJECTORY
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine WriteTraj(T1, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! writes out the trajectory structure to a specefied file
      type(T_Trajectory), intent(in) :: T1
      integer(kind=DefInt), optional :: nfile_in

      real(kind=DefReal) :: Coupling(T1%NumStates - 1, T1%NumDimensions), ForceVector(T1%NumDimensions)
      integer(kind=DefInt) :: nf, & ! unit number to write to
                              nstate, i, & !
                              npart, n, ns, ns2

      logical :: ZCouplingCurrent(T1%NumStates)

! set which file we are writing to
      nf = 6
      if (present(nfile_in)) nf = nfile_in

      npart = T1%NumParticles
      ns = T1%StateID
      nstate = T1%NumStates

!Translate new derivative matrix into old variables
      ForceVector(:) = -T1%ElecStruc%DerivMat(ns, ns, :)
      ZCouplingCurrent(:) = T1%ESFlags%ZDerivCurrent(ns, :)
      ZCouplingCurrent(ns) = .false.
      ns2 = 0 ! stupid index for coupling array, how I hate this!
      do i = 1, nstate
         if (i == ns) cycle
         ns2 = ns2 + 1
         Coupling(ns2, :) = T1%ElecStruc%DerivMat(ns, i, :)
      end do

1     format(i6, 34x, a) ! integer   format with label
2     format(f12.4, 28x, a) ! time type format with label
3     format(ES23.15, 22x, a) ! real      format with label
5     format(2(2x, i6), 24x, a) ! integer   format with label
6     format(1x, l1, 38x, a) ! logical   format with label

11    format(3(1x, ES23.15E3)) ! 3N vector format, full double precision
13    format(6(1x, f12.4))
14    format(10(1x, l1))

      write (nf, 1) T1%TrajID, ' / Traj ID'
      write (nf, 1) T1%CBF, ' / CBF ID'
      write (nf, 1) T1%StateID, ' / State ID'
      write (nf, 1) T1%Ms, ' / State Ms'
      write (nf, 1) T1%ParentID, ' / Parent ID'
      write (nf, 5) T1%CentID, ' / Centroid ID'
      write (nf, 2) T1%DeadTime, ' / Dead time'
      write (nf, 6) T1%zCent, ' / Centroid'
      write (nf, 6) T1%triplet, ' / Triplet'
      if (glzStoSwiss) then
         write (nf, 2) T1%SWISS%SelectionTime, ' / Selection time'
         write (nf, 2) T1%SWISS%BirthDate, ' / Birth date '
         write (nf, 3) T1%SWISS%ParentOverlap, ' / Overlap with Parent'
      end if
      write (nf, *) "# Last spawn"
      write (nf, 13) T1%LastSpawn
      write (nf, *) "# Spawn time"
      write (nf, 13) T1%SpawnTime

      write (nf, *) "# Positions"
      do n = 1, npart
         write (nf, 11) T1%Particle(n)%get_pos()
      end do
      write (nf, *) "# Momenta"
      do n = 1, npart
         write (nf, 11) T1%Particle(n)%get_mom()
      end do

      write (nf, *) "# Energies"
      write (nf, 11) T1%ElecStruc%PotEn

      write (nf, *) "# Dipoles Current"
      write (nf, 14) T1%ESFlags%ZDipolesCurrent
      if (T1%ESFlags%ZDipolesCurrent) then
         do i = 1, nstate
            write (nf, *) "# Dipole", i
            write (nf, 11) T1%ElecStruc%Dipole(i, :)
         end do
      end if

      write (nf, *) "# Transition Dipoles Current"
      write (nf, 14) T1%ESFlags%ZTransDipsCurrent
      if (T1%ESFlags%ZTransDipsCurrent) then
         do i = 2, nstate
            write (nf, *) "# Transition Dipole ", i
            write (nf, 11) T1%ElecStruc%TransDipole(i, :)
         end do
      end if

      write (nf, *) "# Force "
      write (nf, 11) ForceVector

      write (nf, *) "# Coupling status array"
      write (nf, 14) ZCouplingCurrent

      do i = 1, nstate
         if (i == ns) cycle
         write (nf, *) "# Coupling history ", ns, i
         write (nf, 11) T1%CoupHist(:, i)
      end do

      ns2 = 0 ! stupid index for coupling array, how I hate this!
      do i = 1, nstate
         if (i == ns) cycle
         ns2 = ns2 + 1
         if (ZCouplingCurrent(i)) then
            write (nf, *) "# coupling ", ns, i
            write (nf, 11) Coupling(ns2, :)
         end if
      end do

      write (nf, 3) T1%Phase, " / Phase"
      write (nf, *) T1%Amplitude, " / Amplitude"

      call WriteElecStruc(T1%ElecStruc, nf)

      return
   end subroutine WriteTraj

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine ReadTraj(T1, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! writes out the trajectory structure to a specefied file
      type(T_Trajectory), intent(inout) :: T1
      integer(kind=DefInt), optional :: nfile_in

      real(kind=DefReal) :: Coupling(T1%NumStates - 1, T1%NumDimensions), ForceVector(T1%NumDimensions)

      real(kind=DefReal), allocatable :: vec(:)
      integer(kind=DefInt) :: nf, & ! unit number to write to
                              nstate, i, & !
                              npart, ndim, n, ns, ns2

      logical :: ZCouplingCurrent(T1%NumStates)

! set which file we are reading from
      nf = 5
      if (present(nfile_in)) nf = nfile_in

      npart = T1%NumParticles
      nstate = T1%NumStates

      read (nf, *) T1%TrajID
      read (nf, *) T1%CBF
      read (nf, *) T1%StateID
      read (nf, *) T1%Ms
      read (nf, *) T1%ParentID
      read (nf, *) T1%CentID
      read (nf, *) T1%DeadTime
      read (nf, *) T1%zCent
      read (nf, *) T1%triplet
      if (glzStoSwiss) then
         read (nf, *) T1%SWISS%SelectionTime
         read (nf, *) T1%SWISS%BirthDate
         read (nf, *) T1%SWISS%ParentOverlap
      end if
      read (nf, *) ! Last spawn
      read (nf, *) T1%LastSpawn
      read (nf, *) ! Last time
      read (nf, *) T1%SpawnTime

      ndim = T1%Particle(1)%NumDimensions
      allocate (vec(ndim))
! read position
      read (nf, *) ! # Positions
      do n = 1, npart
         read (nf, *) (vec(i), i=1, ndim)
         call T1%Particle(n)%set_pos(vec)
      end do

! read momenta
      read (nf, *) ! # Momenta
      do n = 1, npart
         read (nf, *) (vec(i), i=1, ndim)
         call T1%Particle(n)%set_mom(vec)
      end do

      read (nf, *) ! # Energies
      read (nf, *) T1%ElecStruc%PotEn

      ns = T1%StateID

      read (nf, *) !# Dipoles Current
      read (nf, *) T1%ESFlags%ZDipolesCurrent
      if (T1%ESFlags%ZDipolesCurrent) then
         do i = 1, nstate
            read (nf, *) !# Dipole,i
            read (nf, *) T1%ElecStruc%Dipole(i, :)
         end do
      end if

      read (nf, *) !# Transition Dipoles Current
      read (nf, *) T1%ESFlags%ZTransDipsCurrent
      if (T1%ESFlags%ZTransDipsCurrent) then
         do i = 2, nstate
            read (nf, *) !# Transition Dipole , i
            read (nf, *) T1%ElecStruc%TransDipole(i, :)
         end do
      end if

      read (nf, *) ! # Force
      read (nf, *) ForceVector

      read (nf, *) ! # Coupling status array
      read (nf, *) ZCouplingCurrent

      do i = 1, nstate
         if (i == ns) cycle
         read (nf, *) !Coupling history, ns, i
         read (nf, *) T1%CoupHist(:, i)
      end do

      ns2 = 0 ! stupid index for coupling array, how I hate this!
      do i = 1, nstate
         if (i == ns) cycle
         ns2 = ns2 + 1
         if (ZCouplingCurrent(i)) then
            read (nf, *) ! # coupling ,ns, i
            read (nf, *) Coupling(ns2, :)
         end if
      end do

!translate old format into new derivmat
      T1%ElecStruc%DerivMat = 0.0d0
      T1%ElecStruc%DerivMat(ns, ns, :) = -ForceVector
      T1%ESFlags%ZDerivCurrent(ns, ns) = .true.
      ns2 = 0
      do i = 1, nstate
         if (i == ns) cycle
         ns2 = ns2 + 1
         T1%ElecStruc%DerivMat(ns, i, :) = Coupling(ns2, :)
         T1%ElecStruc%DerivMat(i, ns, :) = -Coupling(ns2, :)
         T1%ESFlags%ZDerivCurrent(i, ns) = ZCouplingCurrent(i)
         T1%ESFlags%ZDerivCurrent(ns, i) = ZCouplingCurrent(i)
      end do

      read (nf, *) T1%Phase
      read (nf, *) T1%Amplitude
      call ReadElecStruc(T1%ElecStruc, nf)
      T1%ESFlags%zESExists = .true.
      T1%ESFlags%zPotEnCurrent = .true.
      if (zRedoRestartES) then
         T1%ESFlags%zPotEnCurrent = .false.
         T1%ESFlags%zDerivCurrent = .false.
         T1%ESFlags%zTransDipsCurrent = .false.
         T1%ESFlags%zDipolesCurrent = .false.
      end if

   end subroutine ReadTraj

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!           ELECTRONIC STRUCTURE
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! This should be moved to a the electronic structure module
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine WriteElecStruc(ES, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_ElecStruc), intent(in) :: ES
      integer(kind=DefInt), optional :: nfile_in

      integer(kind=DefInt) :: nf ! unit number to write to
! WTF: We shouldn't just blidly ignore all errors!
      integer(kind=DefInt) :: ierr ! catch errors that may be thrown

! set which file we are writing to
      nf = 6
      if (present(nfile_in)) nf = nfile_in

3     format(10(1x, es15.8))

      write (nf, *) "# Orbitals"
      write (nf, 3, iostat=ierr) ES%OldOrbitals
      write (nf, *) "# CI vectors"
      write (nf, 3, iostat=ierr) ES%OldCIVecs
      write (nf, *) "# Overlap matrix"
      write (nf, 3, iostat=ierr) ES%OverlapMatrix
      write (nf, *) "# TC Blob"
      write (nf, 3, iostat=ierr) ES%OldBlob
      write (nf, *) "# MSPT2 Coeffs"
      write (nf, 3, iostat=ierr) ES%OldMSPT2C
      write (nf, *) "# Electronic Phases"
      write (nf, *, iostat=ierr) ES%ElecPhase

      return
   end subroutine WriteElecStruc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine ReadElecStruc(ES, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_ElecStruc), intent(inout) :: ES
      integer(kind=DefInt), optional :: nfile_in

      integer(kind=DefInt) :: nf ! unit number to write to
! WTF: We shouldn't just blidly ignore all errors!
      integer(kind=DefInt) :: ierr ! catch errors that may be thrown

! set which file we are writing to
      nf = 5
      if (present(nfile_in)) nf = nfile_in

3     format(10(1x, es15.8))

      read (nf, *)
      read (nf, 3, iostat=ierr) ES%OldOrbitals
      read (nf, *)
      read (nf, 3, iostat=ierr) ES%OldCIVecs
      read (nf, *)
      read (nf, 3, iostat=ierr) ES%OverlapMatrix
      read (nf, *)
      read (nf, 3, iostat=ierr) ES%OldBlob
      read (nf, *)
      read (nf, 3, iostat=ierr) ES%OldMSPT2C
      read (nf, *)
!DEBUG - changed v from read(nf,3,iostat=ierr)
      read (nf, *, iostat=ierr) ES%ElecPhase
      return
   end subroutine ReadElecStruc

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!           PARTICLE
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine WriteParticle(P, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Particle), intent(in) :: P
      integer(kind=DefInt), optional :: nfile_in

      integer(kind=DefInt) :: nf ! unit number to write to

! set which file we are writing to
      nf = 6
      if (present(nfile_in)) nf = nfile_in

1     format(i6, 34x, a) ! integer   format with label
3     format(e17.10, 24x, a) ! real      format with label

      write (nf, 1) P%ParticleID, " / Particle ID"
      write (nf, *) P%Elmnt, " / Element"
      write (nf, 3) P%Width, " / Width"
      write (nf, 3) P%Mass, " / Mass"
      write (nf, 3) P%Charge, " / Charge"
      write (nf, 3) P%AtomicNum, " / Atomic number"

   end subroutine WriteParticle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine ReadParticle(P, nfile_in)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_Particle), intent(inout) :: P
      integer(kind=DefInt), optional :: nfile_in

      integer(kind=DefInt) :: nf ! unit number to write to

! set which file we are writing to
      nf = 5
      if (present(nfile_in)) nf = nfile_in

      read (nf, *) P%ParticleID
      read (nf, *) P%Elmnt
      read (nf, *) P%Width
      read (nf, *) P%Mass
      read (nf, *) P%Charge
      read (nf, *) P%AtomicNum

   end subroutine ReadParticle

end module RestartModule
