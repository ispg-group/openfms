!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Driver for minimization routines
!!
!!    Bond and angle constrains may be set using the same input as for a
!!    dynamics run. They are enforced using quadratic energy penalties
!!    as follows:
!!       -# Penalty multipliers are initialized at AMBER TIP3P values
!!       -# A minimzation is run until convergence OR for MaxIter steps,
!!           whichever is sooner
!!       -# If the constraints are not satisfied to within cndToler,
!!            all penalty multipliers are doubled and the optimization
!!            is re-run for 20 steps
!!       -# If the constraints are still not satisfied after 20 doublings
!!            of the penalty multipliers, the optimization is aborted
!!
!!    In this way, we require that the constraints be satisfied EVEN if
!!    the optimization has not convereged. For sufficiently high
!!    multiplier values, the optimization will converge regardless of
!!    the rest of the system.
!!    \param iState State on which to find minimum
!!    @ingroup propagation
!<
subroutine FMS_MinSearch(T1, iState)
   use GlobalModule
   use MinSearchModule
   use TrajectoryModule
   use TrajectoryIOModule, only: FMS_WriteFXYZ
   use FMSModule, only: FMS_Shutdown
   use ElecStrucModule
   use QM_MM_Module
   implicit none
   type(T_Trajectory), intent(inout) :: T1
   integer(kind=DefInt), intent(in), optional :: iState
   integer(kind=DefInt) :: IUnit, iter, iProt, iStartProt
   integer(kind=DefInt) :: its
   integer(kind=DefInt) :: StepNum

   character(len=2) :: cType

   real(kind=DefReal) :: DefConfD, DefConfK

! Namelist variables for minimization protocol
   integer, parameter :: MaxMins = 50
   integer(kind=DefInt) :: NumMins
   integer(kind=DefInt) :: MinProtocol(MaxMins)
   integer(kind=DefInt) :: MinSteps(MaxMins)
   integer(kind=DefInt) :: QMCycleSteps(MaxMins)
   real(kind=DefReal) :: toler(MaxMins)
   real(kind=DefReal) :: ConfineD(MaxMins), ConfineK(MaxMins)
   logical :: zRestart
   integer(kind=DefInt) :: StepToPrint

   namelist /minsearch/ NumMins, MinProtocol, MinSteps, QMCycleSteps, toler, &
      zRestart, StepToPrint, ConfineD, ConfineK

   glzTimeStepRejection = .false.

!---------------------------------------------------------------
!     Read and check namelist input
!----------------------------------------------------------------

!     Set defaults
   StepToPrint = 5
   NumMins = 0
   toler = 0.0001
   MinSteps = -1
   MinProtocol = -1
   QMCycleSteps = 0
   zRestart = .false.
   ConfineD = -1.0d0
   ConfineK = -1.0d0

!     Read the namelist
   open (file=trim(FMSWorkingDir)//'MinSearch.dat', newunit=iUnit, status='old')
   read (iUnit, minsearch)
   close (iUnit)

!     Check input for errors
   if (NumMins <= 0) call FMS_DieError('Variable NumMins not set in MinSearch.dat.')
   if (NumMins > MaxMins) call FMS_DieError('NumMins too large - reset MaxMins in MinSearch.f.')
   do iProt = 1, NumMins
      if (MinSteps(iProt) <= 0) call FMS_DieError('MinSteps not set correctly in Minsearch.dat.')
      if (MinProtocol(iProt) <= 0) call FMS_DieError('MinProtocol not set correctly in Minsearch.dat.')
      if (MinProtocol(iProt) == Min_MM .and. (.not. qczQMMM)) then
         call FMS_DieError('ERROR in MinSearch: MM minimization requested, but this is not QM/MM.')
      end if
      if (MinProtocol(iProt) == Min_MM .and. (.not. eszMMForce)) then
         call FMS_DieError('ERROR in MinSearch: MM-only minimization not supported for this model.')
      end if
   end do

!     Set module parameters
   mnnStepToPrint = StepToPrint
   if (qczQMMM) then
      DefConfK = qcdConfineK
      DefConfD = qcdConfineD
   end if
   mncFileName = 'Optimization'

!     Allocation
   call mnT1%create(T1%NumParticles, T1%NumStates)
   mnT1 = T1

   if (zRestart) then
      call FMS_DieError('MinSearch: zRestart not implemented')
   else
      iStartProt = 1
   end if

!---------------------------------------------------------------
!     Run minimizations
!----------------------------------------------------------------

   if (present(iState)) T1%StateID = iState
4123 format('--------- Optimization on state ', i3, ' -----------')
   write (fmiOut, 4123) T1%StateID

   if (.not. zRestart) then
      call FMS_WriteFXYZ(mnT1, 'MinSearch.xyz', 'Initial Geometry', FirstTime=.true.)
   end if

   do iProt = iStartProt, NumMins

!     Write minimization info
      write (fmiOut, *) 'MinSearch step ', iProt

!     Set spherical confining potential
      if (ConfineD(iProt) > 0.0) then
         qcdConfineD = ConfineD(iProt)
         qcdConfineK = ConfineK(iProt)
      else
         qcdConfineD = DefConfD
         qcdConfineD = DefConfD
      end if

!     Set up MinSearchModule
      mnMinType = MinProtocol(iProt)

!     Restart if requested
      select case (mnMinType)
      case (Min_MM)
         cType = 'MM'
      case (Min_QM)
         cType = 'QM'
      case (Min_Both)
         ctype = 'mn'
      case (Min_FixQM)
         ctype = 'Fx'
      case default
         call FMS_DieError('ERROR in FMS_MinSearch: unknown minimization type')
      end select
!         MinName=FMS_NumberedName(trim(cType),iProt)
      if (zRestart .and. iStartProt == iProt) then

      else
         iter = MinSteps(iProt)
         StepNum = 0
      end if

!     Run the minimization
      its = FMS_Minimizer(iter, toler(iProt), nChargeCycle=QMCycleSteps(iProt), StepNum=StepNum)

   end do

!     We are done - exit FMS
   T1 = mnT1
   call mnT1%destroy()
   write (fmiOut, *) 'Minimization complete.'
   write (fmiOut, *) 'Result written to MinSearch.xyz'
   call FMS_WriteFXYZ(T1, 'MinSearch.xyz', 'Optimization result')
   call FMS_Shutdown(T1=T1)

end subroutine FMS_MinSearch
