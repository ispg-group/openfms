! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module TerachemModule
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use GlobalModule, only: DefReal, DefInt, DefComp, fmiOut, &
                           gldmaxediff, gldEShift, ntrip, nsing, glzxfaims, &
                           glirestTC, gldTripletShift, glzRejectAllStateFlip, &
                           FMS_RejectStep, FMS_DieError, set_terachem_finalizer
   use TrajectoryModule
#ifdef TeraChem
   use mpi, only: MPI_COMM_NULL, MPI_COMM_WORLD, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_SUCCESS, &
                  MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_CHARACTER, MPI_STATUS_SIZE
#endif
   implicit none

! EGH: There is no reason FMS needs to know every possible option in TeraChem.
!      All TeraChen options should be placed in a file named
!      "misc_options", which will be read by FMS and sent to
!      TeraChem. The number of states FMS wants and the threshold
!      to compute the coupling will be read from FMS's input and
!      sent to TeraChem. For all other variables, we will make it
!      TeraChem's job to square its internal variable with what FMS
!      thinks it will be getting.
   private
   public :: InitTerachem, tc_finalize, RunTerachem

   ! Sleep interval in miliseconds while waiting for TC calculation to finish.
   ! For now, this is hardcoded.
   integer, parameter :: MPI_MILISLEEP = 50

#ifdef TeraChem
   integer, public, save :: newcomm = MPI_COMM_NULL ! Initialized in connect_to_terachem subroutine

   ! Provide an explicit interface for MPI_Send/MPI_Recv with an assumed-type buffer
   ! to satisfy newer compiler versions with stricter argument checks across different buffer types.
   ! NOTE: The assumed-type (type(*)) concept has only been standardized in Fortran2018.
   ! If your compiler chokes on this code, you can simply comment out these interfaces.
   interface
      subroutine MPI_Send(buf, count, datatype, dest, tag, comm, ierr)
         implicit none
         type(*), dimension(*), intent(in) :: buf
         integer, intent(in) :: count, datatype, dest, tag, comm
         integer, intent(inout) :: ierr
      end subroutine MPI_Send

      subroutine MPI_Recv(buf, count, datatype, source, tag, comm, status, ierr)
         import MPI_STATUS_SIZE
         implicit none
         type(*), dimension(*), intent(inout) :: buf
         integer, intent(in) :: count, datatype, source, comm, tag
         integer, intent(inout) :: ierr
         integer, intent(inout) :: status(MPI_STATUS_SIZE)
      end subroutine MPI_Recv
   end interface
#endif

contains

#ifdef TeraChem
   subroutine InitTerachem(NumParticles, NumStates, tc_port_name)
      use GlobalModule, only: BohrToAng
      use ElecStrucModule
      integer(kind=DefInt), intent(in) :: NumParticles, NumStates
      character(len=:), allocatable, intent(inout) :: tc_port_name

      integer, parameter :: CLEN = 128
      integer(kind=DefInt) :: natoms, nqmmm
      real(kind=DefReal), allocatable :: atcoords(:, :)
      character(len=:), allocatable :: server_name
      character(len=2), allocatable :: atom_types(:)
      integer(kind=DefInt) :: bufints(20)
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=DefInt) :: i, ierr, FMSinit
      character(len=128) :: dbuffer(2, 128) ! 128 options, max
      character(len=128) :: tbuffer
      integer(kind=DefInt) :: noptions
      character(len=256) :: unitname
      character(len=1) :: txt
      logical :: GeomInAngs
      integer :: iunit
      character(len=:), allocatable :: tc_options_file
      logical :: file_exists
      character(len=500) :: errmsg

      ! Setup on first program call
      write (fmiOut, *) '    >>>> FMS / TC <<<<'
      server_name = get_server_name()

      natoms = NumParticles
      nqmmm = esNMM

      ! ---------------------------------------------------
      ! Initialization: Connect to "terachem_port", set
      ! newcomm (global), send relevant namelist variables.
      ! ---------------------------------------------------
      call connect_to_terachem(server_name, tc_port_name)

      ! -------------------
      ! Send job info to TC
      ! -------------------
      write (dbuffer(:, 1), '(a,/,E23.16)') 'coupthre', gldMaxEDiff
      write (dbuffer(:, 2), '(a,/,i0)') 'fmsnumstates', NumStates

      ! Read the "tc_options" file and send its contents to TeraChem
      ! If `tc_options` does not exist, try `misc_options` for backward compatibility
      tc_options_file = 'tc_options'
      inquire (file=tc_options_file, exist=file_exists)
      if (.not. file_exists) then
         tc_options_file = 'misc_options'
      end if
      open (newunit=iunit, file=tc_options_file, status='old', action='read', iostat=ierr, iomsg=errmsg)
      if (ierr /= 0) then
         call FMS_DieError(trim(errmsg))
      end if

      noptions = 2
      do
         read (iunit, '(A)', end=10) tbuffer
         noptions = noptions + 1
         write (dbuffer(:, noptions), '(A)') tbuffer, ''
      end do

10    close (iunit)

      noptions = noptions + 1
      write (dbuffer(:, noptions), '(a)') 'end', ''

      write (fmiOut, *) 'Initializing Electronic Structure...'

!   Start sending stuff to TeraChem. This is received by
!   Mpi::init in mpi_base.cpp on the TeraChem side.

!   Send input parameters to TC (the startfile)
      call MPI_Send(dbuffer, 2 * clen * size(dbuffer, 2), MPI_CHARACTER, 0, 2, newcomm, ierr)

      ! ---------------------------------------------
      ! Begin sending the inital geometry to terachem
      ! ---------------------------------------------
      if (.not. allocated(atcoords)) allocate (atcoords(3, natoms))
      if (.not. allocated(atom_types)) allocate (atom_types(natoms))

      ! Read initial coordinates from file
      ! TODO: Use the existing subroutine for this!
      open (newunit=iunit, file='Geometry.dat', action='read', status='old')

      GeomInAngs = .false.
      read (iunit, *) unitname
      txt = unitname(7:7)
      if (txt == 'a' .or. txt == 'A') then ! Angstrom
         GeomInAngs = .true.
      else if (txt == 'b' .or. txt == 'B') then ! Bohr
         GeomInAngs = .false.
      end if

!   TeraChem wants bohr, let's send it in bohr
      read (iunit, *)
      do i = 1, natoms
         read (iunit, *) atom_types(i), atcoords(:, i)
         if (GeomInAngs .eqv. .true.) then
            atcoords(:, i) = atcoords(:, i) / BohrtoAng
         end if
      end do
      close (iunit)

!   Just send the QM part here
!   Send natoms
      bufints(1) = natoms - nqmmm
      call MPI_Send(bufints, 1, MPI_INTEGER, 0, 2, newcomm, ierr)

!   Send atom types (needed to determine whether we have d-functions)
      call MPI_Send(atom_types, 2 * (natoms - nqmmm), MPI_CHARACTER, 0, 2, newcomm, ierr)

! ----------------------------------------------------------
!   TeraChem will catch from: Fms::fmsinit_recv in fms.cpp
! ----------------------------------------------------------

!   Send ESinit/natoms
      FMSinit = 1
      bufints(1) = FMSinit
      bufints(2) = natoms - nqmmm
      bufints(3) = nqmmm
      call MPI_Send(bufints, 3, MPI_INTEGER, 0, 2, newcomm, ierr)

!   Send atom types
      call MPI_Send(atom_types, 2 * size(atom_types), MPI_CHARACTER, 0, 2, newcomm, &
                    ierr)

!   Send coordinates
      call MPI_Send(atcoords, 3 * natoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
      write (fmiOut, '(a)') 'Sent initial data to TeraChem'

! ----------------------------------------------------------
!   TeraChem will pitch from: Fms::fmsinit_send in fms.cpp
! ----------------------------------------------------------

!   Receive CI length, basis functions and blob size
      call MPI_Recv(bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, newcomm, status, ierr)

      esLCiVec = bufints(1) ! Length of CI vector
      esNBasis = bufints(2) ! Number of basis functions
      esBlobSize = bufints(3) ! Blob size

      write (fmiOut, '(a)') 'TeraChem initialized!'

   end subroutine InitTerachem

   ! If the TeraChem server is named something other than "terachem_port",
   ! then read the server name from the file, "tc_input".
   function get_server_name() result(server_name)
      character(len=256) :: server_name
      integer :: iost, iunit
      namelist /tc/ server_name

      ! Default value
      server_name = 'terachem_port'

      ! Read Namelist
      open (newunit=iunit, file='tc_input', action='read', status='old', iostat=iost)
      if (iost /= 0) then
         write (fmiOut, *) '"tc_input" file does not exist, using default server name'
         return
      end if
      write (fmiOut, *) 'Reading server_name from file "tc_input"'
      read (iunit, nml=tc)
      close (iunit)
   end function get_server_name

   ! Connect to the TeraChem server.
   subroutine connect_to_terachem(server_name, port_name)
      use mpi, only: MPI_MAX_PORT_NAME, MPI_Init, &
                     MPI_Comm_size, MPI_COMM_CONNECT, MPI_INFO_NULL

      character(len=*), intent(in) :: server_name
      character(len=:), allocatable, intent(inout) :: port_name
      integer :: nproc
      integer :: ierr

      write (*, *) 'Terachem MPI Initialization'

      ! Initialize MPI.
      call MPI_Init(ierr)
      if (ierr /= 0) then
         call FMS_DieError('MPI_Init failed!')
      end if

      ! Check the number of processes, only 1 is allowed!
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      if (nproc /= 1) then
         write (*, '(A,I0)') 'Number of MPI processes: ', nproc
         call FMS_DieError('Only one MPI process should be running')
      end if

      ! -----------------------------------
      ! Look for server_name, get port name
      ! After 60 seconds, exit if not found
      ! -----------------------------------
      if (.not. allocated(port_name)) then
         call lookup_port_via_nameserver(server_name, port_name)
      end if

      ! Establish new communicator via port name
      call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, newcomm, ierr)

      ! Register tc_finalize to be called from FMS_DieError if needed
      call set_terachem_finalizer(tc_finalize)

      write (*, *) 'FMS <--> TeraChem communication established!'
      write (*, *)
   end subroutine connect_to_terachem

   ! Tell TeraChem to stop itself if FMS dynamics are finished.
   subroutine tc_finalize()
      use mpi, only: MPI_Finalize
      integer :: ierr
      logical :: initialized

      call MPI_Initialized(initialized, ierr)
      if (.not. initialized) then
         return
      end if

      if (newcomm /= MPI_COMM_NULL) then
         write (*, *) 'Finalizing communication with TC...'
         call MPI_Send([0.0d0], 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr)
      else
         write (*, *) 'Communication with TC not established, skipping finalization.'
      end if

      call MPI_Comm_free(newcomm, ierr)

      call MPI_Finalize(ierr)
   end subroutine tc_finalize

   ! Look for server_name via MPI nameserver, get port name
   subroutine lookup_port_via_nameserver(server_name, port_name)
      use mpi, only: MPI_wtime, MPI_lookup_name, MPI_comm_set_errhandler, &
                     MPI_MAX_PORT_NAME, MPI_ERRORS_RETURN, MPI_ERRORS_ARE_FATAL, MPI_INFO_NULL
      character(len=*), intent(in) :: server_name
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      ! Give up trying to find the port after this many seconds
      integer, parameter :: connection_timeout = 30
      real(DefReal) :: timer
      integer :: ierr

      port_name = ''

      write (fmiOut, '(a)') 'Looking up TeraChem server under name: '//trim(server_name)
      flush (fmiOut)

      timer = MPI_WTIME()

      ! Setting MPI_ERRORS_RETURN error handler allows us to retry
      ! failed MPI_LOOKUP_NAME() call.
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)

      do

         call MPI_LOOKUP_NAME(server_name, MPI_INFO_NULL, port_name, ierr)
         if (ierr == MPI_SUCCESS) then
            ! Workaround for a bug in hydra_nameserver for MPICH versions < 3.3
            if (len_trim(port_name) == 0) then
               write (fmiOut, '(A)') 'Found empty port, retrying...'
            else
               write (fmiOut, '(A)') 'Found port: '//trim(port_name)
               exit
            end if
         end if

         ! Timeout after connection_timeout seconds
         if ((MPI_WTIME() - timer) > connection_timeout) then
            call FMS_DieError('Server name "'//trim(server_name)//'" not found.')
         end if

         ! Let's wait a second since too many calls
         ! to MPI_LOOKUP_NAME() can crash the hydra_nameserver process
         write (fmiOut, '(A)') 'Waiting for TeraChem port...'
         flush (fmiOut)
         call milisleep(1000)

      end do

      ! Set the default error handler back
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)
      if (ierr /= 0) then
         write (fmiOut, *) 'WARNING: Could not set MPI error handler to MPI_ERRORS_ARE_FATAL!'
      end if

   end subroutine lookup_port_via_nameserver

   subroutine RunTerachem(T_FMS, iCalcState, jCalcState, CalcCoup)
      use ElecStrucModule
      type(T_Trajectory), intent(inout) :: T_FMS
      integer(kind=DefInt), intent(in) :: iCalcState, jCalcState
      logical, intent(in) :: CalcCoup

      integer(DefInt), save :: first_call = 1
      real(kind=DefReal), allocatable :: atcoords(:, :), tmpcoords(:)
      integer(kind=DefInt) :: bufints(20)
      real(kind=DefReal) :: buf(3)
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=DefInt) :: i, j, k, ierr
      integer(DefInt) :: FMSinit, natoms, doCoup, OldWfn
      real(kind=DefReal), allocatable :: VelRe(:, :), VelIm(:, :)
      real(kind=DefReal), allocatable :: Energ(:)
      real(kind=DefReal), allocatable :: CIvecs(:, :), SMat(:, :), TDip(:), Dip(:), Chg(:)
      real(kind=DefReal), allocatable :: dxyz(:), MO(:, :)
      real(kind=DefReal), allocatable :: tcBlob(:)

      !GAIMS added
      integer(kind=DefInt) :: l, shift
      real(kind=DefReal), allocatable :: SOMAT(:)
      !GAIMS end added

      ! If we're not at time zero, old wavefunction data exists
      OldWfn = 0
      if (first_call == 0 .or. glirestTC == 1) then
         T_FMS%ESFlags%zESExists = .true.
         OldWfn = 1
      end if

      write (*, '(a)') 'Sending data to TeraChem'

      natoms = T_FMS%NumParticles
      if (CalcCoup) then
         doCoup = 1
      else
         doCoup = 0
      end if

      ! ------------------------------------------------
      ! Begin sending data each step to terachem
      ! TeraChem catches from: Fms::receive in fms.cpp
      ! ------------------------------------------------

      FMSinit = 0
      bufints(1) = FMSinit ! Should already be initalized
      bufints(2) = natoms ! We already know this...
      bufints(3) = doCoup ! Does FMS want coupling vectors
      bufints(4) = T_FMS%TrajID ! Trajectory ID number
      bufints(5) = T_FMS%CentID(1) ! Centroid trajectory ID number
      bufints(6) = T_FMS%CentID(2) ! Centroid other trajectory ID number
      bufints(7) = T_FMS%StateID ! Electronic state of the trajectory
      ! FORTRAN indexing
      bufints(8) = OldWfn ! Do we have an old wavefunction?
      bufints(9) = iCalcState - 1 ! TC Target State, C indexing
      bufints(10) = jCalcState - 1 ! TC other state (for couplings), C indexing
      bufints(11) = first_call ! First time we're asking for wavefunctions?
      bufints(12) = glirestTC ! Are we restarting an old calculation?

      call MPI_Send(bufints, 12, MPI_INTEGER, 0, 2, newcomm, ierr)

      glirestTC = 0 ! Only restart TeraChem once

      ! Send TeraChem the current time
      buf(1) = T_FMS%get_time()
      call MPI_Send(buf, 1, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      ! Temp space for coordinates
      if (.not. allocated(tmpcoords)) allocate (tmpcoords(T_FMS%NumDimensions))
      ! Atomic coordinates
      if (.not. allocated(atcoords)) allocate (atcoords(3, natoms))
      ! Real velocity vector
      if (.not. allocated(VelRe)) allocate (VelRe(3, natoms))
      ! Imaginary velocity vector
      if (.not. allocated(VelIm)) allocate (VelIm(3, natoms))
      ! Temp space to recieve gradients/couplings
      if (.not. allocated(dxyz)) allocate (dxyz(3 * natoms))
      ! Temp space to recieve spin-orbit couplings
      if (.not. allocated(SOMAT)) allocate (SOMAT(18 * T_FMS%NumStates * T_FMS%NumStates))

      !--------------------
      !  State properties
      !--------------------
      ! Energies of each electronic state
      if (.not. allocated(Energ)) allocate (Energ(T_FMS%NumStates))
      ! Transition dipoles, always ground to excited states
      if (.not. allocated(TDip)) allocate (TDip((T_FMS%NumStates - 1) * 3))
      ! Dipoles for each electronic state
      if (.not. allocated(Dip)) allocate (Dip((T_FMS%NumStates) * 3))
      ! Atomic charges for the current electronic state
      if (.not. allocated(Chg)) allocate (Chg(natoms))

      !---------------------
      !  Wavefunction data
      !---------------------
      ! Diabatic molecular orbitals
      if (.not. allocated(MO)) allocate (MO(esNBasis, esNBasis))
      ! Diabatic CI vectors for each electronic state
      if (.not. allocated(CIvecs)) allocate (CIvecs(esLCiVec, T_FMS%NumStates))
      ! Overlap matrix of electronic wavefunctions
      if (.not. allocated(SMat)) allocate (SMat(T_FMS%NumStates, T_FMS%NumStates))
      ! The blob!!!
      if (.not. allocated(tcBlob)) allocate (tcBlob(esBlobSize))

      ! Get atomic coordinates from FMS
      tmpcoords = T_FMS%get_pos()

      ! Don't be a dipshit, work in bohr
      do i = 1, natoms
         atcoords(:, i) = tmpcoords(3 * (i - 1) + 1:3 * i)
      end do

      ! Get real momentum
      tmpcoords = T_FMS%get_mom()
      do i = 1, natoms
         do j = 1, 3
            VelRe(j, i) = tmpcoords((i - 1) * 3 + j)
         end do
      end do

      ! Get imaginary momentum
      tmpcoords = T_FMS%get_mom2()
      do i = 1, natoms
         do j = 1, 3
            VelIm(j, i) = tmpcoords((i - 1) * 3 + j)
         end do
      end do

      ! Old wavefunction data
      MO = T_FMS%ElecStruc%OldOrbitals
      CIVecs = transpose(T_FMS%ElecStruc%OldCIVecs)
      tcBlob = T_FMS%ElecStruc%OldBlob

      ! Zero properties
      Energ = 0.d0
      TDip = 0.d0
      Dip = 0.d0

      ! Send coordinates
      call MPI_Send(atcoords, 3 * natoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      ! Send previous MOs
      call MPI_Send(MO, esNBasis * esNBasis, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      ! Send previous CI vectors
      call MPI_Send(CIvecs, esLCiVec * T_FMS%NumStates, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      ! Send previous blob
      call MPI_Send(tcBlob, esBlobSize, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      ! Send momenta
      call MPI_Send(VelRe, 3 * natoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
      call MPI_Send(VelIm, 3 * natoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

      write (*, '(a)') 'Done sending data to TeraChem'

      ! ------------------------------------------------
      ! TeraChem pitches from: Fms::send in fms.cpp
      ! ------------------------------------------------

      call wait_for_terachem(newcomm)
      write (*, '(a)') 'Receiving data from TeraChem'

      !    Receive energies from TC
      call MPI_Recv(Energ, T_FMS%NumStates, &
     &    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    Store energies with the shift
      if (gldTripletShift /= 0d0) then
         write (*, '("Shifting triplets ",I2," to ",I2," by ",E15.7)') &
      &        NSing + 1, T_FMS%NumStates, gldTripletShift
         Energ(NSing + 1:T_FMS%NumStates) = Energ(NSing + 1:T_FMS%NumStates) + gldTripletShift
      end if
      T_FMS%ElecStruc%PotEn = Energ + gldEShift
      T_FMS%ESFlags%ZPotEnCurrent = .true.

      !    Receive transition dipoles from TC
      call MPI_Recv(TDip, (T_FMS%NumStates - 1) * 3, &
     &      MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG,  &
     &      newcomm, status, ierr)

      !    Store transition dipoles
      do i = 1, T_FMS%NumStates - 1
         T_FMS%ElecStruc%TransDipole(i + 1, :) = TDip(3 * (i - 1) + 1:3 * (i - 1) + 3)
      end do
      T_FMS%ESFlags%ZTransDipsCurrent = .true.

! xf added - preliminary
      if (glzxfaims) then
         ! Caution, two states for the moment only.
         T_FMS%ElecStruc%TransDipolexf(1) = T_FMS%ElecStruc%TransDipole(2, 1)
         T_FMS%ElecStruc%TransDipolexf(2) = T_FMS%ElecStruc%TransDipole(2, 2)
         T_FMS%ElecStruc%TransDipolexf(3) = T_FMS%ElecStruc%TransDipole(2, 3)
         T_FMS%ElecStruc%TransDipolexf(4) = T_FMS%get_time()
         T_FMS%ESFlags%ZTransDipsCurrentxf = .true.
      end if
! xf added end

      !    Receive dipole moments from TC
      call MPI_Recv(Dip, T_FMS%NumStates * 3,  &
     &      MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG,  &
     &      newcomm, status, ierr)

      !    Store dipole moments
      do i = 1, T_FMS%NumStates
         T_FMS%ElecStruc%Dipole(i, 1:3) = Dip(3 * (i - 1) + 1:3 * (i - 1) + 3)
      end do
      T_FMS%ESFlags%ZDipolesCurrent = .true.

      !    Receive atomic charges from TC
      call MPI_Recv(Chg, natoms, MPI_DOUBLE_PRECISION, &
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    Store atomic charges
      do i = 1, natoms
         T_FMS%Particle(i)%Charge = Chg(i)
      end do

      !    Receive diabatic MOs from TC
      call MPI_Recv(MO, esNBasis * esNBasis, MPI_DOUBLE_PRECISION, &
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    Store MOs
      T_FMS%ElecStruc%OldOrbitals = MO

      !    Receiving diabatic CI vectors from TC
      call MPI_Recv(CIvecs, T_FMS%NumStates * (esLCiVec),  &
     &    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    Receiving overlap matrix from TC
      call MPI_Recv(SMat, T_FMS%NumStates * T_FMS%NumStates, &
     &    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    FMS thinks it's phasing the CI vector, but we already
      !    phased it. Let FMS read the overlap matrix in order
      !    to detect state flipping.
      call FMS_CheckOverlap(T_FMS, SMat, T_FMS%NumStates)

      !    Store the CI vectors
      T_FMS%ElecStruc%OldCIVecs = transpose(CIvecs)

      !    Receiving the blob from TC
      call MPI_Recv(tcBlob, esBlobSize, MPI_DOUBLE_PRECISION,  &
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

      !    Store the blob!!!
      T_FMS%ElecStruc%OldBlob = tcBlob

      !    Receive derivative matrix from TeraChem
      !
      !    Let TeraChem send the entire derivative
      !    matrix, eventually we should fix this,
      !    but TeraChem is supposed to be fast, right?
      !
      !    TeraChem will maintain consistent phase, we
      !    can't really phase here, anyway. This only
      !    works if TeraChem calculates couplings between
      !    diabatic wavefunctions, which, in general,
      !    it doesn't and, in some cases, it can't.
      T_FMS%ElecStruc%DerivMat = 0.d0
      do i = 1, T_FMS%NumStates
         do j = i, T_FMS%NumStates
            call MPI_Recv(dxyz, 3 * natoms, MPI_DOUBLE_PRECISION,  &
       &        MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
            if (i == j) then
               T_FMS%ElecStruc%DerivMat(i, j, :) = dxyz
            else
               T_FMS%ElecStruc%DerivMat(i, j, :) = dxyz
               T_FMS%ElecStruc%DerivMat(j, i, :) = -1.d0 * dxyz
            end if
            T_FMS%ESFlags%ZDerivCurrent(i, j) = .true.
            T_FMS%ESFlags%ZDerivCurrent(j, i) = .true.
         end do
      end do

      !GAIMS added
      !     Recieve spin-orbit matrix from Terachem
      !
      !     Right now it would be easiest to get one array of
      !     3*nstates by 3*nstates. This is of course useless
      !     it this way, as we would only need half of the matrix,
      !     even less if we have many singlet states, but that's
      !     an improvement I can do later, once this works
      if (NTrip /= 0) then
         T_FMS%ElecStruc%SOMat = 0.d0
         SOMAT = 0.d0
         call MPI_Recv(SOMAT, 18 * T_FMS%NumStates * T_FMS%NumStates, &
                       MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
         shift = 1
         do i = 1, T_FMS%NumStates
            do j = 1, T_FMS%NumStates
               do k = 1, 3
                  do l = 1, 3
                     T_FMS%ElecStruc%SOMat(i, j, k, l) = cmplx(SOMAT(shift), SOMAT(shift + 1), kind=DefComp)
                     T_FMS%ESFlags%ZSOMCurrent(i, j, k, l) = .true.
!       write(*,'(i5,2f18.12)') shift, SOMAT(shift), SOMAT(shift+1)
!       write(*,'(4i5,f18.12,sp,f18.12,"i")') i,j,k,l,T_FMS%ElecStruc%SOMat(i,j,k,l)
                     shift = shift + 2
                  end do
               end do
            end do
         end do
      end if
      !GAIMS end added

      write (*, '(a)') 'Done receiving data from TeraChem'

      first_call = 0
   end subroutine RunTerachem
!>
!!    Check the overlap matrix for state flipping
!!
!!    Record any states flipping sign.
!!    Reject the timestep if two states have flipped character.
!<
   subroutine FMS_CheckOverlap(T1, SMat, nstate)
      use ElecStrucModule
      type(T_Trajectory), intent(inout) :: T1
      real(kind=DefReal), intent(in) :: SMat(nstate, nstate)
      integer(kind=DefInt), intent(in) :: nstate

      integer(kind=DefInt) :: IState, JState, ioccstate

      real(kind=DefReal) :: SCI_II, SCI_IJ, SCI_JI, SCI_JJ
      real(kind=defReal) :: E_occ, E_I, E_J, E_occ_old, E_I_old, E_J_old

      logical :: ImportantState

1999  format('========================================================')
2000  format('WARNING: Trajectory jumped an intersection.')

!     Don't bother correcting phase if there was no electronic structure previously stored
      if (.not. T1%ESFlags%zESExists) then
         T1%ElecStruc%ElecPhase = 1.0d0
         return
      end if

      ioccstate = T1%StateID
!
!     Make sure the CI Vector is consistent since the last timestep
!

      do IState = 2, T1%NumStates
         do JState = 1, IState - 1

!
!    Calculate dot products SCI_KL=CIVec(State K)(t) .dot. CIVec(State L)(t-dt)
!
            SCI_II = SMat(IState, IState)
            SCI_JJ = SMat(JState, JState)
            SCI_IJ = SMat(IState, JState)
            SCI_JI = SMat(JState, IState)

            write (*, *) 'Overlap matrix for states:', IState, JState
            write (*, *) SCI_II, SCI_JI
            write (*, *) SCI_IJ, SCI_JJ
            write (*, *)

            if (abs(SCI_II) + abs(SCI_JJ) < abs(SCI_IJ) + abs(SCI_JI)) then

!
!     States I and J  have flipped character
!     Attempt to reject the timestep

!     Is at least one of these states important (i.e. is it the occupied state, or
!     coupled to the occupied state)?  This is determined from the energy gaps.
!     Note: in the case of a centroid, there are two relevant states, but only
!     one is indexed in the data structure.  This is fine, however, as both states
!     will be considered important based on the energy criterion (if the energy
!     criterion is not met, then there is no coupling between the states, so it is
!     unnecessary to keep track of them.

               ImportantState = .false.
!     Gather energies
               E_occ = T1%ElecStruc%PotEn(ioccstate)
               E_I = T1%ElecStruc%PotEn(Istate)
               E_J = T1%ElecStruc%PotEn(Jstate)

!     Consider the previous timestep energy gaps as well, because a state flip
!     may have artificially taken a coupled state out of the coupling region
               E_occ_old = T1%ElecStruc%OldPotEn(ioccstate)
               E_I_old = T1%ElecStruc%OldPotEn(Istate)
               E_J_old = T1%ElecStruc%OldPotEn(Jstate)

               if (abs(E_I - E_occ) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_J - E_occ) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_I_old - E_occ_old) < gldMaxEDiff) ImportantState = .true.
               if (abs(E_J_old - E_occ_old) < gldMaxEDiff) ImportantState = .true.

               write (fmiOut, 1999)
               write (fmiOut, 2000)
               call T1%print_id()
1001           format('SCI_', 2i2, ':', f10.6, ' SCI_', 2i2, ':', f10.6)
               write (fmiOut, 1001) jState, jState, SCI_JJ, jState, iState, SCI_JI
               write (fmiOut, 1001) iState, jState, SCI_IJ, iState, iState, SCI_II
               write (fmiOut, 1999)
!     If user option RejectAllstateFlip=.false., only reject the step if one of
!      the states was an important state
               if ((glzRejectAllStateFlip .or. ImportantState) .and. .not.  &
           &        T1%ESFlags%zIgnoreErrors) then
                  call FMS_RejectStep(.true.)
                  return
               end if
            end if

         end do
      end do

!     Correct phase based on diagonal overlap.
!     Note: this is also how we correct phase when two states have flipped.
!     The assumption that allows us to consider only the diagonal overlap in determining
!     the phase is that in a single step the states are not rotated by more than +-pi/4 .
!     Note - if this assumption does not hold then we would be unable to distinguish a rotation of
!     |t| > 90degree from a rotation of t-180 with two sign flips.
      do IState = 1, T1%NumStates
         SCI_II = SMat(IState, IState)
         if (SCI_II < 0.0d0) then
            write (6, '("  Sign of CIVecs changed in state",i3)') iState
            T1%ElecStruc%ElecPhase(iState) = -1.0d0 * T1%ElecStruc%ElecPhase(iState)
         end if
      end do

   end subroutine FMS_CheckOverlap

   subroutine wait_for_terachem(tc_comm)
      integer, intent(in) :: tc_comm
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      logical :: ready

      if (MPI_MILISLEEP <= 0) return

      ! The idea here is to reduce the CPU usage of MPI_Recv() by taking a brief nap.
      ! In most MPI implementations, MPI_Recv() is actively polling the other end
      ! (in this case TeraChem) and consumes a whole CPU core. That's clearly wasteful,
      ! since we're typically waiting for a long time for the ab initio result.
      !
      ! Some implementation provide an option to change this behaviour,
      ! but I didn't figure out any for MPICH.
      ! Based according to an answer here:
      ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu

      ready = .false.
      ! TODO: we need to somehow make sure that
      ! we don't wait forever if TeraChem crashes.
      ! At this moment, this is ensured at the BASH script level.
      do while (.not. ready)
         call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, ready, status, ierr)
         call milisleep(MPI_MILISLEEP)
      end do
   end subroutine wait_for_terachem

#else
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Stub routines when FMS is compiled without TeraChem interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitTerachem(NumParticles, NumStates, tc_port_name)
      integer(kind=DefInt), intent(in) :: NumParticles, NumStates
      character(len=:), allocatable, intent(inout) :: tc_port_name
      call FMS_DieError('ERROR: not compiled for use with TeraChem.')
   end subroutine InitTerachem

   subroutine tc_finalize()
   end subroutine tc_finalize

   subroutine RunTerachem(T_FMS, iCalcState, jCalcState, CalcCoup)
      type(T_Trajectory), intent(inout) :: T_FMS
      integer(kind=DefInt), intent(in) :: iCalcState, jCalcState
      logical, intent(in) :: CalcCoup
      call FMS_DieError('ERROR: not compiled for use with TeraChem.')
   end subroutine RunTerachem

#endif

   !> Sleep for x miliseconds
   !> This uses a usleep syscall
   subroutine milisleep(milisec)
      use, intrinsic :: iso_c_binding, only: c_int, c_int32_t
      integer, intent(in) :: milisec
      integer(kind=c_int32_t) :: usec
      integer(kind=c_int) :: c_err

      ! Interface to usleep syscall
      interface
         ! https://cyber.dabamos.de/programming/modernfortran/sleep.html
         ! int usleep(useconds_t useconds)
         function usleep(useconds) bind(c, name='usleep')
            import :: c_int, c_int32_t
            implicit none
            integer(kind=c_int32_t), intent(in), value :: useconds
            integer(kind=c_int) :: usleep
         end function usleep
      end interface

      ! TODO: Based on usleep(2) manpage, we probably should not sleep more than a second
      usec = int(milisec * 1000, c_int32_t)
      c_err = usleep(usec)
      ! DH: If you ever see this warning, please report it!
      if (c_err /= 0) then
         write (error_unit, '(a,i0)') 'WARNING: usleep syscall returned an error: ', c_err
      end if
   end subroutine milisleep

end module TerachemModule
