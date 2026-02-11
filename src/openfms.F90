!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Ab Initio Multiple Spawning with Split-Operator Approximations
!!
!!    The main() routine (or entrance to FMS code from the electronic
!!    structure package). Drives the program from initialization to
!!    memory allocation to initial conditions selection to dynamics to
!!    shutdown.
!!    @ingroup memory
!<
program OpenFMS
   use FMSModule
   use GlobalModule, only: DefInt, DefReal, fmiOut, &
                           FMS_DeleteFile, FmsWorkingDir, timei, FMS_DieError, &
                           fmzWriteEveryStep, gldTimeStep, glzMinSearch
   use RandomModule, only: FMS_ranb, initialize_fortran_prng
   use BundleModule
   use BundleCalcsModule, only: FMS_UpdateBundle
   use BundleIOModule, only: FMS_Output
   use InitialModule, only: iniRndSeed, inInitState
   use ElecStrucModule, only: FMS_ESInit
   use RestartModule, only: inIRestart, RestartTime, getRestart
   use PropagationModule, only: FMS_SetTimeStep
   implicit none

   type(T_TrajectoryBundle) :: Bundle
   integer(kind=DefInt) :: NumTraj, l
   integer(kind=DefInt) :: numparticles, numstates
   real(kind=DefReal) :: simulationtime, dt
   character(len=256) :: tc_port_name
   real(8) :: RndNum

!-----------------------------------------------------------------------
!                          STARTUP
!-----------------------------------------------------------------------

   ! Get working directory
   call getenv('PWD', FMSWorkingDir)
   l = len(trim(FMSWorkingDir))
   if (l /= 0) then
      FMSWorkingDir(l + 1:) = '/'
   end if

   tc_port_name = ''
   call get_cmdline(tc_port_name)
   ! Start the timer
   call cpu_time(timei)

   ! Get parameters for this run and populate parameter modules
   call FMS_ReadNameList(NumParticles, NumStates, NumTraj, SimulationTime)

   if (inIRestart == 0) then
      call clean_output_files()
   end if

   call print_header()

!--------------------------------------------------------------------
!                        INITIALIZATION
!--------------------------------------------------------------------

   ! Initialize random number generator
   if (inIRndSeed == 0) then
      inIRndSeed = int(timei, 4)
      write (fmiOut, *) 'Random seed from clock: ', inIRndSeed
      if (inIRndSeed == 0) then
         call FMS_DieError('Failed to get non-zero random seed from clock.')
      end if
   end if
   call initialize_fortran_prng(iniRndSeed)
   RndNum = FMS_ranb(inIRndSeed)

   ! Initialize electronic structure package
   write (fmiOut, *) ' Initializing electronic structure ...'
   flush (fmiOut)
   ! This call initializes some internal variables for electronic structure,
   ! but does not actually call the electronic structure code yet...
   call FMS_ESInit(NumParticles, NumStates)

   ! Here we make the first call to TeraChem et al!
   ! Must be separate from FMS_ESInit to prevent circular module dependency
   call initialize_interface(NumParticles, NumStates, tc_port_name)

   ! Allocate memory structures, read initial geometry
   call initialize_bundle(Bundle, NumTraj, NumParticles, NumStates, inInitState)

!-----------------------------------------------------------------------
!     DYNAMICS SETUP
!-----------------------------------------------------------------------

   ! Construct initial system state
   if (inIRestart == 0) then
      call FMS_InitialCond(Bundle)
      call FMS_DeleteFile('Checkpoint.txt')
   else
      call GetRestart(Bundle, RestartTime)
      call FMS_UpdateBundle(Bundle)
   end if

!     Force delete any junk trajectories
   call FMS_RemoveDead(Bundle)

!     Write simulation information to FMS.out
   call print_simulation_info(Bundle)

!------------------------------------------------------------------
!     SPECIALTY CALCULATIONS
!------------------------------------------------------------------

!     Minimum search branches here
   if (glzMinSearch) call FMS_MinSearch(Bundle%Trajectory(1))

!-----------------------------------------------------------------------
!      BEGIN DYNAMICS
!-----------------------------------------------------------------------

!     Create output files and erase any old ones (if not restarting),
!     and write output for t=0
   if (inIRestart == 0) then
      call FMS_Output(Bundle, FirstTime=.true.)
   else
      call FMS_Output(Bundle, FirstTime=.false.)
   end if

   write (fmiOut, '(a,F0.3,a)') 'Propagate until t = ', SimulationTime, ' a.u.'
   flush (fmiOut)

   do while (Bundle%CurrentTime < SimulationTime)

      if (Bundle%NumTraj == 0) then
         write (fmiOut, *) 'No more live trajectories.'
         exit
      end if

!        Determine dt for this step
      call FMS_SetTimeStep(Bundle, dt)

!        Run dynamics from t to t + gldTimeStep
!        call cpu_time(timedi)
      call FMS_Dynamics(Bundle, dt, Bundle%CurrentTime + gldTimeStep)
!        call cpu_time(timede)
!        write(fmiOut,*)'Time for dynamics step:',timede-timedi

      ! Write output here by default
      if (.not. fmzWriteEveryStep) call FMS_Output(Bundle)

   end do

!-----------------------------------------------------------------------
!        EXIT
!-----------------------------------------------------------------------
   call FMS_Shutdown(Bundle, terminate=.false.)

contains

!>
!!    Erases old output files
!!    @ingroup output
!<
   subroutine clean_output_files
      use GlobalModule, only: DefInt, MaxTrajLimit, FMS_DeleteFile, FMS_NumberedFileName
      use RestartModule, only: iniRestart
      implicit none
      integer(kind=DefInt) :: iTraj

!     Do not erase old files if we're restarting
      if (iniRestart /= 0) return

!     Delete old files unless it's a restart
      call FMS_DeleteFile('CFxn.dat')
      call FMS_DeleteFile('Coup.dat') !written to in Spawn
      call FMS_DeleteFile('SOCoup.dat') !written to in Spawn
      call FMS_DeleteFile('E.dat')
      call FMS_DeleteFile('FailSpawn.log') !written to in Spawn
      call FMS_DeleteFile('H.dat')
      call FMS_DeleteFile('IterList') !written to in augCG
      call FMS_DeleteFile('N.dat')
      call FMS_DeleteFile('S.dat')
      call FMS_DeleteFile('SDot.dat')
      call FMS_DeleteFile('Spawn.log') !written to in Spawn
      call FMS_DeleteFile('TDip.dat')

!     Delete old numbered files
      do ITraj = 1, MaxTrajLimit

         call FMS_DeleteFile(FMS_NumberedFileName('Amp', ITraj)) !written to in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('Bonds', ITraj)) !written to in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('Angles', ITraj)) !written to in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('Dihedrals', ITraj)) !written in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('Pyram', ITraj)) !written to in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('Charge', ITraj))
         call FMS_DeleteFile(FMS_NumberedFileName('CIVec', ITraj))
         call FMS_DeleteFile(FMS_NumberedFileName('Coup', ITraj))
         call FMS_DeleteFile(FMS_NumberedFileName('SOCoup', ITraj))
         call FMS_DeleteFile(FMS_NumberedFileName('positions', ITraj)//'.xyz') !written in analysis
         call FMS_DeleteFile(FMS_NumberedFileName('positions', ITraj)//'.dcd') !written in analysis

         call FMS_DeleteFile(FMS_NumberedFileName('Spawn', ITraj)) !written to in Spawn
         call FMS_DeleteFile(FMS_NumberedFileName('TDip', ITraj))
         call FMS_DeleteFile(FMS_NumberedFileName('TrajDump', ITraj)) !written to in writefbundle
         call FMS_DeleteFile(FMS_NumberedFileName('PotEn', iTraj))
      end do

   end subroutine clean_output_files

!>
!!    Print program header to FMS.out
!!    @ingroup output
!<
   subroutine print_header()
      character(len=256) :: host
      character(len=*), parameter :: divider = &
                                     '--------------------------------------------------------------------------------'
      character(len=*), parameter :: title = 'OpenFMS'
      character(len=*), parameter :: subtitle = 'FMS for the Masses...'
      character(len=*), parameter :: credit = 'Courtesy TJM and MBN'

      call hostnm(host)

      write (fmiOut, '(a)') divider
      write (fmiOut, '(a)') repeat(' ', (len(divider) - len(title)) / 2)//title
      write (fmiOut, *)
      write (fmiOut, '(a)') repeat(' ', (len(divider) - len(subtitle)) / 2)//subtitle
      write (fmiOut, '(a)') repeat(' ', (len(divider) - len(credit)) / 2)//credit
      write (fmiOut, *)

      ! Defined in generated build_info.F90
      call print_build_info(fmiOut)

      write (fmiOut, '(1x, a, t25, a)') 'host:', trim(adjustl(host))
      write (fmiOut, *)
      write (fmiOut, '(a)') divider
      write (fmiOut, *)

!     NOTE(danielhollas): We do not print out the following list
!     of authors because it is incomplete, and also many features mentioned
!     below are not actually present in the open-source version of the code
!     so it would be misleading and confusing to mention them.
!     That said, we should think about how to credit all contributors properly,
!     both in the manual and in the source code.
!     A good start would be a list of unique commiters to the repo,
!     which one can get by running:
!
!     $ git shortlog -e -s -n
!
!     (this also prints the number of commits per author,
!     but that is a highly misleading number that should not be mistaken
!     for the "amount" of contribution of a given person)

!     write(fmiOut,*) 'Authors:                                    '
!     write(fmiOut,*) 'Michal Ben-Nun (FMS algorithm)              '
!     write(fmiOut,*) 'Christian R Evenhuis (Interpolated Dynamics)'
!     write(fmiOut,*) 'William J Glover (MPI)                      '
!     write(fmiOut,*) 'Benjamin J Levine (FMS-MolPro, MS-CASPT2)   '
!     write(fmiOut,*) 'Toshifumi Mori (MS-CASPT2)                  '
!     write(fmiOut,*) 'Mitchell T Ong (FMS-Columbus, SMD)          '
!     write(fmiOut,*) 'Michael S Schuurman (FMS-Columbus)          '
!     write(fmiOut,*) 'Hongli Tao (MS-CASPT2, laser fields)        '
!     write(fmiOut,*) 'Aaron M Virshup (QM/MM)                     '
!     write(fmiOut,*) 'Basile F. E. Curchod (Stochastic selection) '
!     write(fmiOut,*) 'Yorick Lassmann (AIMSWISS)                  '
!     write(fmiOut,*) 'and Todd J Martinez (FMS algorithm)         '
   end subroutine print_header

!>
!!    Print basic simulation information to FMS.out
!!    @ingroup output
!<
   subroutine print_simulation_info(Bundle)
      use GlobalModule
      use BundleModule
      use TrajectoryCalcsModule, only: FMS_Weight
      use SelectionModule, only: print_stochastic_selection_params
      use SpawnModule, only: print_spawning_parameters
      use RestartModule, only: iniRestart

      type(T_TrajectoryBundle), intent(in) :: Bundle

      integer(kind=DefInt) :: IParticle, ITraj

      select case (gliModel)
      case (FMSZERO)
         write (fmiOut, *) '     FMSZero: fake potential surface    '
      case (TEMPLATE)
         write (fmiOut, *) '     System call template interface     '
      case default
         write (fmiOut, *) '         Unknown potential model ', gliModel
      end select
      write (fmiOut, *) ' Number Electronic States: ', Bundle%NumStates
      write (fmiOut, *) ' Potential Model:          ', gliModel

      write (fmiOut, *)
      write (fmiOut, *)

      call print_spawning_parameters(fmiOut)

      if (glzStochastic) then
         call print_stochastic_selection_params()
      end if

      if (inIRestart /= 0) then
1000     format(' Restart at time: ', f10.2)
         write (fmiOut, 1000) Bundle%CurrentTime
         write (fmiOut, *) ' Live trajectories: ', Bundle%NumTraj
         write (fmiOut, *) ' Dead trajectories: ', Bundle%NumDeadTraj
      end if

      write (fmiOut, *)

      do ITraj = 1, Bundle%NumTraj
         write (fmiOut, '(a,i0  )') ' Trajectory #', Bundle%Trajectory(iTraj)%TrajID
         write (fmiOut, '(a,i0  )') ' Electronic State: ', Bundle%Trajectory(ITraj)%StateID
         write (fmiOut, '(a,f0.4)') ' Initial Weight:   ', FMS_Weight(Bundle%Trajectory(ITraj))

         if (inIRestart == 0) then
            write (fmiOut, *)
            do IParticle = 1, Bundle%Trajectory(iTraj)%NumParticles
               write (fmiOut, '(a,i0)') ' Particle #', IParticle
               call Bundle%Trajectory(iTraj)%Particle(IParticle)%print_info(fmiOut)
            end do
         end if

         write (fmiOut, *)

      end do

      write (fmiOut, *)
      flush (fmiOut)
   end subroutine print_simulation_info

   subroutine initialize_interface(NumParticles, NumStates, tc_port_name)
      use GlobalModule, only: gliModel, TEMPLATE, TC
      use TemplateModule, only: FMS_TemplateInit
      use TerachemModule, only: InitTerachem
      integer, intent(in) :: NumParticles, NumStates
      character(len=*), intent(in) :: tc_port_name

      select case (gliModel)
      case (TEMPLATE)
         call FMS_TemplateInit()
      case (TC)
         call InitTerachem(NumParticles, NumStates, tc_port_name)
      case default
         continue
      end select
   end subroutine initialize_interface

!>
!! Allocate and initialize all memory structures
!! \param B1 Trajectory bundle
!! \param NumTraj Number of trajectories to initialize
!! \param NumParticles Number of particles in system
!! \param NumStates Number of electronic states to calculate
!!    @ingroup memory
!<
   subroutine initialize_bundle(B1, NumTraj, NumParticles, NumStates, InitialState)
      use GlobalModule, only: glzAnalysisMode, FMS_DieError
      use ParticleModule
      use TrajectoryModule
      use BundleModule
      use FMSModule, only: FMS_ReadGeometry, FMS_ParticleTypes
      use QM_MM_Module, only: qcNumAS, qcNumQM, qcNumMM, qcNumS, qczQMMM
      implicit none

      type(T_TrajectoryBundle), intent(inout) :: B1
      integer(kind=DefInt), intent(in) :: NumTraj, NumParticles, NumStates, InitialState

      type(t_particle), pointer :: Particle(:)
      type(t_Trajectory), pointer :: T1(:)
      integer(kind=DefInt) :: ITraj, IParticle
      integer(kind=DefInt) :: TotalDim, NumDim
      real(kind=DefReal) :: DMTemp

!     Read basic particle data and create particle structures.
      TotalDim = 0
      NumDim = 3

      if (qczQMMM) call FMS_DieError('QM/MM is currently not available')

      if (.not. glzAnalysisMode) then
         ! Get particle data from Control.dat
         qcNumQM = NumParticles
         qcNumMM = 0
         qcNumAS = 0
         qcNumS = 0
         allocate (Particle(NumParticles))
         do IParticle = 1, NumParticles
            call Particle(IParticle)%create(id=IParticle, numdim=NumDim)
            TotalDim = TotalDim + NumDim
         end do
      end if

!     Allocate memory for trajectories that are initially running
      allocate (T1(NumTraj))
      call B1%create(numtraj=0, &
                     numdeadtraj=0, &
                     numstates=NumStates, &
                     numparticles=NumParticles, &
                     ncbfs=0)
      do ITraj = 1, NumTraj
         call T1(ITraj)%create(NumParticles, NumStates)
         T1%TrajID = iTraj
      end do

!     Copy particles into 1st trajectory
      do iParticle = 1, NumParticles
         T1(1)%Particle(iParticle) = Particle(iParticle)
      end do

!     Setup for starting trajectories
      if (.not. (qczQMMM .or. glzAnalysisMode)) then
         call FMS_ReadGeometry(T1(1))
      end if

      call FMS_ParticleTypes(T1(1))

      T1%StateID = InitialState
      DMTemp = 1.0
      T1(1)%Amplitude = DMTemp
      do ITraj = 2, NumTraj
         T1(ITraj) = T1(1)
         ! DH: This seems weird???
!        T1(ITraj)%Amplitude=0.0
         T1(ITraj)%Amplitude = DMTemp
      end do

!     Copy trajectories into bundle and clean up
      do ITraj = 1, NumTraj
         call B1%add_traj(T1(iTraj))
      end do

      do IParticle = 1, NumParticles
         call Particle(IParticle)%destroy()
      end do
      deallocate (Particle)

      do ITraj = 1, NumTraj
         call T1(ITraj)%destroy()
      end do
      deallocate (T1)

   end subroutine initialize_bundle

   subroutine get_cmdline(tc_port_name)
      use, intrinsic :: iso_fortran_env, only: output_unit
      character(len=*), intent(inout) :: tc_port_name
      character(len=256) :: arg
      integer :: i

      i = 0
      do while (i < command_argument_count())

         i = i + 1
         call get_command_argument(i, arg)

         select case (arg)
         case ('-h', '--help')
            call print_help()
            stop 0
         case ('-v', '--version')
            call print_version(output_unit)
            stop 0
         case ('--tc-port-name')
            i = i + 1
            call get_command_argument(i, tc_port_name)
            if (trim(tc_port_name) == '') then
               call FMS_DieError('Empty --tc-port-name argument. Provide port name of the TeraChem server')
            end if
         case default
            call print_help()
            call FMS_DieError('Invalid command line argument "'//trim(arg)//'"')
         end select

      end do
   end subroutine get_cmdline

   subroutine print_help()
      print '(a)', 'OpenFMS: Ab initio Multiple spawning for the masses'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -h, --help       print help and exit'
      print '(a)', '  -v, --version    print version'
   end subroutine print_help

end program OpenFMS
