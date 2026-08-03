!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!! Reads input from Control.dat
!<
subroutine FMS_ReadNameList(NumParticles, NumStates, NumTraj, SimulationTime)
   use, intrinsic :: iso_fortran_env, only: error_unit
   use GlobalModule
   use QM_MM_Module
   use ToyModelModule, only: izmaylov_params, GAIMS_model_params
   use SpawnModule, only: spawn_params, OMAX_DEFAULT
   use InitialModule
   use SamplingModule
   use SMDModule
   use RestartModule
   use XFAIMSModule, only: xfaims_params
!TM
! note: this can be eventually moved to another directory,
!       passing the NumParticles and Dummy-related variables.
   use ElecStrucModule
!

   implicit none
   ! Newline character
   character(len=*), parameter :: NL = new_line('A')
   character(len=*), parameter :: AVAILABLE_MODELS = 'TOY, TC, QUANTICS'

   integer(kind=DefInt), intent(inout) :: NumParticles, NumStates, NumTraj
   real(kind=DefReal), intent(inout) :: SimulationTime
   integer(kind=DefInt) :: NumMM

   logical :: zFatal
   integer :: ios
   character(len=500) :: errmsg

   integer, parameter :: MaxArray = 100 !WJG Note: this is rather small- we have systems with more atoms than this...

   character(len=2) :: IntegType

   logical :: TimestepRejection
   logical :: DieOnMinStep
   logical :: RejectAllStateFlip
   logical :: WriteEveryStep
   real(kind=DefReal) :: TimeStep
   real(kind=DefReal) :: NormCons, EnergyCons, CoupTimeStep
   real(kind=DefReal) :: MinTimeStep, OLapThresh
   real(kind=DefReal) :: NormStepCons, EnergyStepCons, tStepThresh
   real(kind=DefReal) :: CFThresh
   real(kind=DefReal) :: CSThresh
   real(kind=DefReal) :: DGamma
   real(kind=DefReal) :: NormThresh
   real(kind=DefReal) :: EquilTStep
   real(kind=DefReal) :: ExShift
   real(kind=DefReal) :: TripletShift
   real(kind=DefReal) :: FModeSharp(MaxArray)
   real(kind=DefReal) :: InitGap
   real(kind=DefReal) :: InitGapWidth
   logical :: SelectState
   real(kind=DefReal) :: MaxEDiff
   real(kind=DefReal) :: OMax
   real(kind=DefReal) :: Omin_parent, OMax_intra, OMax_inter
   real(kind=DefReal) :: NumGradStep
   real(kind=DefReal) :: PopToSpawn
   real(kind=DefReal) :: RegThresh
   real(kind=DefReal) :: DecoherenceTime
!bfec
   real(kind=DefReal) :: StochasticThresh
   ! Switch on AIMSWISS?
   logical :: StochasticSwiss
   logical :: StochasticStateSpecific

! xf added
   real(kind=DefReal) :: CoupFieldTimeStep
   real(kind=DefReal) :: tspwnixf
   real(kind=DefReal) :: tspwnfxf
   real(kind=DefReal) :: f0_xfr
   real(kind=DefReal) :: polx_xfr
   real(kind=DefReal) :: poly_xfr
   real(kind=DefReal) :: polz_xfr
   real(kind=DefReal) :: freq_xfr
   real(kind=DefReal) :: t0_xfr
   real(kind=DefReal) :: sigma_xfr
   real(kind=DefReal) :: CEP_xfr
   logical :: IgnoreStateAferField
   logical :: onespawnonly_xfr, XFAIMS
! xf added end

! GAIMS added
   integer(kind=DefInt) :: NumSinglets
   integer(kind=DefInt) :: NumTriplets
   real(kind=DefReal) :: SOCThresh
   real(kind=DefReal) :: ShiftTrip
   logical :: SPA1_SOC_model
! GAIMS added end

   logical :: zMMFile

   real(kind=DefReal) :: QuenchToler

   real(kind=DefReal) :: ConfineD
   real(kind=DefReal) :: ConfineK

   integer(kind=DefInt) :: ibrown
   integer(kind=DefInt) :: iAddQuanta(MaxArray)
   integer(kind=DefInt) :: IDark
   integer(kind=DefInt) :: IMethod
   integer(kind=DefInt) :: iModeSharp(MaxArray)
   character(len=32) :: InitialCond
   integer(kind=DefInt) :: InitState

   integer(kind=DefInt) :: IRestart

   integer(kind=DefInt) :: IRestartTraj(MaxTrajLimit)
   integer(kind=DefInt) :: ISaddle
   integer(kind=DefInt) :: MaxTraj
   integer(kind=DefInt) :: MirrorState
   character(len=32) :: Model
   integer(kind=DefInt) :: MoldenStep
   integer(kind=DefInt) :: MultiSpawn
   integer(kind=DefInt) :: naddquanta
   integer(kind=DefInt) :: NCycles
   integer(kind=DefInt) :: nRelaxSteps
   integer(kind=DefInt) :: nFixSteps
   integer(kind=DefInt) :: NSteps
   integer(kind=DefInt) :: NStepToPrint
   integer(kind=DefInt) :: NEquiStepPrint
   integer(kind=DefInt) :: NTemp
   integer(kind=DefInt) :: numas
   integer(kind=DefInt) :: NumIseed
   integer(kind=DefInt) :: nums
   real(kind=DefReal) :: Temperature
   integer :: ForceKill(50)

   real(kind=DefReal) :: IzmOmegax
   real(kind=DefReal) :: IzmOmegay
   real(kind=DefReal) :: IzmXshift
   real(kind=DefReal) :: IzmYshift
   real(kind=DefReal) :: IzmDeltaE
   real(kind=DefReal) :: IzmCoupC
   real(kind=DefReal) :: Grsigma

   integer*4 :: IRndSeed

   logical :: AnalysisMode
   logical :: Brown
   logical :: BrownCon
   logical :: StochasticOlap
   logical :: CentroidApprox
   logical :: Constrain
   logical :: EnergyAdjust
   logical :: Equi
   logical :: EquiRes
   logical :: EquiCon
   logical :: FirstGauss
   logical :: FullyCoupled
   logical :: Stochastic
   logical :: GenSolvent
   logical :: InitBright
   logical :: InitDark
   logical :: IterInv
   integer(kind=DefInt) :: IgnoreState
   logical :: AvH
   integer(kind=DefInt) :: AvHNStates
   integer(kind=DefInt) :: AvHStates(MaxArray)
   logical :: MinSearch
   logical :: MirrorBasis
   logical :: NormInitial
   logical :: OSAmp
   logical :: SharpEnergy
   logical :: SpawnCoupV
   logical :: WriteMolden
   logical :: ZQMMM
   logical :: ZTurnPoint
   logical :: CentNGrad

!     Output control
   logical :: zTrajFile
   logical :: zPotEnFile
   logical :: zNDatFile
   logical :: zEdatFile
   logical :: zAmpFile
   logical :: zCIVecFile
   logical :: zBundMatFile
   logical :: zCoupFile
   logical :: zChargeFile
   logical :: zCorrFile
   logical :: zTDipoleFile
   logical :: zDipoleFile
   logical :: zQMRRFile
   logical :: zAllText
   logical :: zXYZ
   logical :: zDCD
   logical :: zForce
   logical :: zPCOlap
   logical :: zSOME
   logical :: zSOCeff

! Constraint info
   integer(kind=DefInt) :: IBond(2, MaxArray), NBonds
   integer(kind=DefInt) :: IAngle(3, MaxArray), NAngles
   integer(kind=DefInt) :: IDihedral(4, MaxArray), NDihedrals
   integer(kind=DefInt) :: IPyram(4, MaxArray), NPyrams

!  List of SMD Input
   logical :: SMD
   logical :: VelPull
   real(kind=DefReal) :: ForceConst
   real(kind=DefReal) :: PullRate
   real(kind=DefReal) :: Force
   integer(kind=DefInt) :: IAtoms(MaxArray), SMDNAtoms
   real(kind=DefReal) :: IDummy(3, MaxArray), NDummy
   real(kind=DefReal) :: ElongCut
   logical :: autodirect

!TM
! Electronic structure info
   integer(kind=DefInt) :: NDummyParticles
   real(kind=DefReal), allocatable :: DummyCoeff(:, :)
!

   integer(kind=DefInt) :: i
   ! DH: These variables are just for backwards compatibility for old Control.dat files,
   ! and have no actual effect on the simulation.
   ! TODO: Print a warning when they are used in Control.dat
   logical :: tunnel, EShellOnly, SpawnAdaptive
   integer(kind=DefInt) :: SteepestDescent
   integer(kind=DefInt) :: nCubeOrbs
   integer(kind=DefInt) :: nCubeOrbIndex(MaxArray)
   integer(kind=DefInt) :: NCubeStep

   namelist /control/ AnalysisMode, Brown &
      , BrownCon, CFThresh, CSThresh, DGamma, EnergyAdjust, Equi &
      , EquilTStep, EquiRes, ExShift, FirstGauss &
      , FModeSharp, ForceKill, FullyCoupled, Stochastic &
      , IAddQuanta, IBrown, IDark, StochasticThresh, StochasticOlap &
      , CentroidApprox, TripletShift &
      , IMethod, IModeSharp, InitBright, InitDark &
      , InitialCond, InitState, InitGap, InitGapWidth, SelectState, IntegType &
      , IRestart, IRestartTraj, RestartTime, RestartStep, zRedoRestartES &
      , IRndSeed, ISaddle &
      , IterInv, MaxEDiff, MaxTraj, MirrorBasis, MirrorState &
      , Model, MoldenStep, MultiSpawn, NAddQuanta, NCycles, OSAmp &
      , NormInitial, NSteps, NStepToPrint, NTemp, numas &
      , NumIseed, NumParticles, NumMM, NumInitBasis &
      , nums, NumStates, OMax, OMin_parent, OMax_intra, OMax_inter &
      , PopToSpawn, RegThresh, SharpEnergy, SimulationTime &
      , SpawnCoupV, Temperature, TimeStep &
      , WriteMolden, ZQMMM, ZTurnpoint &
      , nCubeOrbs, nCubeOrbIndex, nCubeStep, MinSearch, constrain &
      , NAngles, NBonds, NPyrams, NDihedrals, IDihedral, IBond, IPyram &
      , IAngle, zForce &
      , zTrajFile, zPotEnFile, zNDatFile, zEdatFile &
      , zCIVecFile, zBundMatFile, zCoupFile, zChargeFile, zCorrFile &
      , zTDipoleFile, zDipoleFile, zQMRRFile, zAllText, zXYZ, zDCD &
      , zPCOlap, zSOME, zSOCeff &
      , EquiCon, GenSolvent, nRelaxSteps, nFixSteps, ConfineD, ConfineK &
      , SMD, VelPull, ForceConst, PullRate, Force, SMDNAtoms, IAtoms &
      , NDummy, IDummy, ElongCut, autodirect &
      , NormThresh, NormCons, NormStepCons &
      , EnergyCons, EnergyStepCons, CoupTimeStep, Timesteprejection &
      , DieOnMinStep, RejectAllStateFlip, tStepThresh, WriteEveryStep &
      , MinTimeStep, OLapThresh, nEquiStepPrint, NumGradStep, CentNGrad &
      , IgnoreState, QuenchToler, zAmpFile, zMMFile, DecoherenceTime &
      , AvH, AvHNStates, AvHStates &
      ! xf added
      , freq_xfr, t0_xfr, CEP_xfr, sigma_xfr, onespawnonly_xfr &
      , f0_xfr, polx_xfr, poly_xfr, polz_xfr &
      , XFAIMS, CoupFieldTimeStep, IgnoreStateAferField &
      , tspwnixf, tspwnfxf &
      ! xf added end
      ! GAIMS added
      , NumSinglets, NumTriplets, SOCThresh, ShiftTrip, SPA1_SOC_model &
      ! GAIMS added end
      , StochasticSwiss, StochasticStateSpecific &
      ! Toy models
      , IzmOmegax, IzmOmegay, IzmXshift, IzmYshift, IzmDeltaE, IzmCoupC &
      , Grsigma &
      ! Backwards compat
      , Tunnel, EShellOnly, SpawnAdaptive, SteepestDescent &
      !TM
      , NDummyParticles, DummyCoeff

   IRestart = 0
   IRestartTraj(:) = 0
   RestartTime = -100.d0
   RestartStep = 2147483647
   DecoherenceTime = 200.d0
   QuenchToler = 1.d-5
   NumGradStep = 0.05
   CentNGrad = .true.
   NormCons = 0.25
   NormStepCons = 0.05
   EnergyStepCons = 0.005
   EnergyCons = 0.02
   tStepThresh = 0.005
   WriteEveryStep = .true.
   OLapThresh = 1.d-3
   TimeStepRejection = .true.
   RejectAllStateFlip = .true.
   DieOnMinStep = .true.
   CoupTimeStep = -1.d0
   MinTimeStep = -1.d0
   NormThresh = 0.000005d0
   nCubeOrbs = 0
   nCubeOrbIndex = 0
   nCubeStep = 0
   AnalysisMode = .false.
   Brown = .false.
   BrownCon = .false.
   CFThresh = 0.005
   CSThresh = 0.1
   Constrain = .false.
   DGamma = 0.0d0
   EnergyAdjust = .false.
   Equi = .false.
   EquilTStep = -1.0
   EquiRes = .false.
   ExShift = 0.0
   TripletShift = 0.0
   FirstGauss = .false.
   IDark = 2
   IMethod = 1
   InitBright = .false.
   InitDark = .false.
   InitialCond = 'NOSAMPLE'
   ForceKill(:) = 0
   IzmOmegax = 0.009557d0
   IzmOmegay = 0.0033515d0
   IzmXshift = 20.07d0
   IzmYshift = 0.d0
   IzmDeltaE = 0.01984d0
   IzmCoupC = 0.0006127d0
   Grsigma = 10.0d0
   FullyCoupled = .true.
   Stochastic = .false.
!bfec
   StochasticThresh = 1.0d-10
   StochasticOlap = .false.
   StochasticStateSpecific = .false.
   StochasticSwiss = .false.
   CentroidApprox = .true.
   IBrown = 1
   InitState = -1
   InitGap = 0.0d0 !Disable initial energy gap selection
   InitGapWidth = 0.0d0
   SelectState = .false. !Disable state selection based on energy gap
   IntegType = 'VV'
   zRedoRestartES = .true.
   IRndSeed = 0
   IgnoreState = 0
   AvH = .false.
   AvHNStates = MaxArray !Will need to check this and reset it to NumStates later
   do i = 1, MaxArray
      AvHStates(i) = i
   end do
   ISaddle = 0
   IterInv = .false.
   MaxEDiff = 0.03d0
   MaxTraj = 100
   MinSearch = .false.
   MirrorBasis = .false.
   MirrorState = 1
   Model = 'UNDEF'
   MoldenStep = 200
   MultiSpawn = 1
   NCycles = 0
   NormInitial = .false.
   OSAmp = .false.
   NSteps = 0
   NStepToPrint = 1
   nEquiStepPrint = 10
   NTemp = 0
   numas = 0
   NumParticles = 0
   NumMM = 0
   NumInitBasis = 1
   nums = 0
   NumStates = -1
!     OMAX_DEFAULT is defined in SpawnModule
   OMax = OMAX_DEFAULT
   OMax_intra = -1
   OMax_inter = -1
   OMin_parent = -1
   PopToSpawn = 0.1
!      QIniKvec=24*0.0
!      QIniPos=24*0.0
   RegThresh = 0.0001
!      Reject=0.7
   SimulationTime = -1.0d0
   SharpEnergy = .false.
   SpawnCoupV = .false.
!      TDB=.false.
!      TDBOMax=0.5
   Temperature = 0.0d0
   TimeStep = -1.0d0
   WriteMolden = .false.
   ZQMMM = .false.
   ZTurnPoint = .true.
   NDihedrals = 0
   NBonds = 0
   NAngles = 0
   NPyrams = 0
   zXYZ = .false.
   zDCD = .false.
   zTrajFile = .false.
   zPotEnFile = .false.
   zNDatFile = .true.
   zEdatFile = .true.
   zAmpFile = .true.
   zCIVecFile = .false.
   zBundMatFile = .false.
   zCoupFile = .false.
   zMMFile = .false.
   zChargeFile = .false.
   zForce = .false.
   zCorrFile = .false.
   zTDipoleFile = .false.
   zDipoleFile = .false.
   zQMRRFile = .false.
   zAllText = .false.
   zPCOlap = .false.
   zSOME = .false.
   zSOCeff = .false.
   ConfineD = 100.0
   ConfineK = 0.0
   GenSolvent = .false.
   EquiCon = .false.
   nRelaxSteps = 0
   nFixSteps = 0
   SMD = .false.
   VelPull = .false.
   PullRate = 0.0
   ForceConst = 0.0
   Force = 0.0
   SMDNAtoms = 0
   NDummy = 0
   ElongCut = 100.00
   autodirect = .true.
   NAddQuanta = 0
   IAddQuanta = 0
!TM
   NDummyParticles = 0
   allocate (DummyCoeff(MaxParticles, MaxParticles))
   DummyCoeff = 0.0

! xf added
   CoupFieldTimeStep = 1.0d0
   tspwnixf = 0.0d0
   tspwnfxf = 0.0d0
   f0_xfr = 0.0d0
   polx_xfr = 0.0d0
   poly_xfr = 0.0d0
   polz_xfr = 0.0d0
   freq_xfr = 0.0d0
   t0_xfr = 0.0d0
   sigma_xfr = 0.0d0
   CEP_xfr = 0.0d0
   onespawnonly_xfr = .false.
   XFAIMS = .false.
   IgnoreStateAferField = .false.
! xf added end

! GAIMS added
   NumSinglets = -1
   NumTriplets = 0
   SOCThresh = 0.d0
   ShiftTrip = 0.d0
   SPA1_SOC_model = .false.
! GAIMS added end

!     Flag for input errors
   zFatal = .false.

!     Read namelist information
!     Note: 'icunit' is defined in GlobalModule and passed to ParticleTypes.f
!     where this file is closed.
   open (newunit=ICUnit, file=trim(FMSWorkingDir)//'Control.dat', &
         status='old', action='read', form='formatted', iostat=ios, &
         iomsg=errmsg)
   if (ios /= 0) then
      call FMS_DieError(trim(errmsg))
   end if
   read (ICUnit, control)

   fmiOut = open_fms_out(iniRestart)

!     Copy global parameters into Global_Module
   select case (trim(lower_case(Model)))
   case ('undef')
      write (error_unit, *) '"Model" must be defined in Control.dat'// &
         NL//'Available models:'//NL//'  '//AVAILABLE_MODELS
      zFatal = .true.
   case ('template')
      gliModel = TEMPLATE
   case ('zero', 'toy')
      gliModel = FMSZERO
   case ('quantics')
      gliModel = 13
   case ('tc')
      gliModel = 14
   case default
      write (error_unit, *) 'ERROR in Control.dat: Unknown model "' &
         //trim(Model)//'"'//NL//'Available models:'//NL &
         //'    '//AVAILABLE_MODELS
      zFatal = .true.
   end select

   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   ! Copy Initial Condition parameters into InitialModule
   select case (trim(lower_case(InitialCond)))
   case ('nosample'); inInitialCond = NOSAMPLE
   case ('wigner'); inInitialCond = WIGNER
   case ('husimi'); inInitialCond = HUSIMI
   case ('quasiclassical'); inInitialCond = QUASICLASS
   case ('boltzmann'); inInitialCond = BOLTZ
   case ('coldwigner'); inInitialCond = COLDWIG
   case ('swarm'); inInitialCond = SWARM
   case default
      write (error_unit, *) 'ERROR in Control.dat:'
      write (error_unit, *) '  Unrecognized initial condition "'//trim(InitialCond)//'"'
      write (error_unit, *) 'Available initial conditions methods:'
      write (error_unit, *) '  NOSAMPLE, WIGNER, QUASICLASSICAL, BOLTZMANN, COLDWIGNER'
      zFatal = .true.
   end select

   NumTraj = NumInitBasis
   gldDecoherenceTime = DecoherenceTime
   glzTimeStepRejection = TimeStepRejection
   glzDieOnMinStep = DieOnMinStep
   glzMinSearch = MinSearch
   glimethod = IMethod
   gldEShift = ExShift
   gldTripletShift = TripletShift
   glcIntegType = IntegType
   gldRegThresh = RegThresh
   gliSaddle = ISaddle
   glzIterInv = IterInv

   glRestartTime = RestartTime
   gldTimeStep = TimeStep
   gldCoupTimeStep = CoupTimeStep
   gldMinTimeStep = MinTimeStep
   gldTStepThresh = tStepThresh

   gliForceKill = ForceKill

   call izmaylov_params%initialize(W1=IzmOmegax, W2=IzmOmegay, XA=IzmXshift, YA=IzmYshift, &
                                   deltaE=IzmDeltaE, coupC=IzmCoupC)
   call GAIMS_model_params%initialize(r_sigma=Grsigma)

   glzFullyCoupled = FullyCoupled
   glzStochastic = Stochastic
   glzSPA1_SOC_model = SPA1_SOC_model
!bfec
   gldStochaThresh = StochasticThresh
   glzCentroids = CentroidApprox
   glzStoOlap = StochasticOlap
   glzStoSwiss = StochasticSwiss
   glzStoStateSpecific = StochasticStateSpecific
   if (glzStoSwiss) then
      gldDecoherenceTime = 10000.d0
      IgnoreState = 0
   end if
   gldMaxEDiff = MaxEDiff
   glzRejectAllStateFlip = RejectAllStateFlip
   glzConstrain = constrain
   glzAnalysisMode = AnalysisMode
   gldSimulationTime = SimulationTime
   gldNormThresh = NormThresh
   gldNormCons = NormCons
   gldNormStepCons = NormStepCons
   gldEnergyCons = EnergyCons
   gldEnergyStepCons = EnergyStepCons
   gldOLapThresh = OLapThresh
   gldNGradStep = NumGradStep
   glzCentNGrad = CentNGrad
   glIgnoreState = IgnoreState
   glzAvH = AvH
   gliAvHNStates = AvHNStates
   gliAvHStates = AvHStates
   fmzWriteEveryStep = WriteEveryStep
   fmNStepToPrint = NStepToPrint
   fmzAmpFile = zAmpFile
   fmiBond = iBond
   fmiDihedral = iDihedral
   fmiAngle = iAngle
   fmiPyram = iPyram
   fmNDihedrals = NDihedrals
   fmNBonds = NBonds
   fmNAngles = NAngles
   fmNPyrams = NPyrams
   fmzXYZ = zXYZ
   fmzTrajFile = zTrajFile
   fmzPotEnFile = zPotEnFile
   fmzMMFile = zMMFile
   fmzNDatFile = zNDatFile
   fmzEdatFile = zEdatFile
   fmzCIVecFile = zCIVecFile
   fmzBundMatFile = zBundMatFile
   fmzCoupFile = zCoupFile
   fmzForce = zForce
   fmzChargeFile = zChargeFile
   fmzCorrFile = zCorrFile
   fmzTDipoleFile = zTDipoleFile
   fmzDipoleFile = zDipoleFile
   fmzQMRRFile = zQMRRFile
   fmzAllText = zAllText
   fmzPCOlap = zPCOlap
   fmzSOME = zSOME
   fmzSOCeff = zSOCeff

!     Copy SMD parameters into SMDModule
   smSMD = SMD
   smVelPull = VelPull
   smPullRate = PullRate
   smForceConst = ForceConst
   smForce = Force
   smNAtoms = SMDNAtoms
   smIAtoms = IAtoms
   smNDummy = NDummy
   smIDummy = IDummy
   smElongCut = ElongCut
   smautodirect = autodirect

   indQuenchToler = QuenchToler
   inzSharpEnergy = SharpEnergy
   iniRndSeed = iRndSeed
   inMirrorState = MirrorState
!      inzTDB=TDB
   inzMirrorBasis = MirrorBasis
   inIDark = IDark
   inInitState = InitState
   inzEnergyAdjust = EnergyAdjust
   indTemperature = Temperature
   innCycles = nCycles
   innSteps = nSteps
   inzEqui = Equi
   inzBrown = Brown
   inzBrownCon = BrownCon
   inzEquiRes = EquiRes
   inIBrown = IBrown
   inDGamma = DGamma
   inzFirstGauss = FirstGauss
   if (initBright .and. initDark) then
      write (error_unit, *) 'ERROR in Control.dat: '//NL// &
         '    initBright and initDark cannot both be true'
      zFatal = .true.
   end if
   inzInitBright = initBright
   inzInitDark = initDark
   inInitGap = InitGap
   inInitGapWidth = InitGapWidth
   inSelectState = SelectState
   inzNormInitial = NormInitial
   inzOSAmp = OSAmp
   inIAddQuanta = IAddQuanta
   inNAddQuanta = NAddQuanta
   inIModeSharp = IModeSharp
   inNTemp = NTemp
   indFModeSharp = FModeSharp
   indEquilTStep = EquilTStep
   inIRestart = IRestart
!bfec
   if (IRestart /= 0) then
      glirestTC = 1
   else
      glirestTC = 0
   end if
   inIRestartTraj = IRestartTraj
   inzEquiCon = EquiCon
   innEquiStepPrint = nEquiStepPrint
   inzGenSolvent = GenSolvent
   inNRelaxSteps = nRelaxSteps
   inNFixSteps = nFixSteps
! XFAIMS parameters
   glzxfaims = XFAIMS
   xfaims_params%sp_spwn_i = tspwnixf
   xfaims_params%sp_spwn_f = tspwnfxf
   xfaims_params%f0 = f0_xfr
   xfaims_params%polx = polx_xfr
   xfaims_params%poly = poly_xfr
   xfaims_params%polz = polz_xfr
   xfaims_params%freq = freq_xfr
   xfaims_params%t0 = t0_xfr
   xfaims_params%sigma = sigma_xfr
   xfaims_params%CEP = CEP_xfr
   xfaims_params%onespawnonly = onespawnonly_xfr
   xfaims_params%CoupFieldTimeStep = CoupFieldTimeStep
   xfaims_params%IgnoreStateAferField = IgnoreStateAferField
! end of XFAIMS parameters

! Check XFASIM paremeters
   if (glzxfaims) then
      if (xfaims_params%sigma <= 0.0d0) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    sigma_xfr must be greater than zero'
         zFatal = .true.
      end if

      if (xfaims_params%freq <= 0.0d0) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    freq_xfr must be greater than zero'
         zFatal = .true.
      end if

      if (xfaims_params%f0 <= 0.0d0) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    f0_xfr must be greater than zero'
         zFatal = .true.
      end if

      if (xfaims_params%CoupFieldTimeStep <= 0.0d0) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    CoupFieldTimeStep must be greater than zero'
         zFatal = .true.
      end if

      if (xfaims_params%sp_spwn_i < 0.0d0) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    tspwnixf must be greater than or equal to zero'
         zFatal = .true.
      end if

      if (xfaims_params%sp_spwn_f <= xfaims_params%sp_spwn_i) then
         write (error_unit, *) 'ERROR in Control.dat: '//NL// &
            '    tspwnfxf must be greater than tspwnixf'
         zFatal = .true.
      end if
   end if
! end Check XFASIM paremeters

! GAIMS added
   !If NumSinglets and NumTriplets not specified, assume all states are singlet
   if (NumTriplets == 0 .and. NumSinglets == -1) NumSinglets = NumStates
   NSing = NumSinglets
   NTrip = NumTriplets

!  OMax parameter had originally multiple meanings that are now covered
!  by separate parameters: OMin_parent, OMax_inter and OMax_intra.
!  To preserve the meaning of existing input files, we only overwrite
!  them with OMax if they are not explicitly specified in Control.dat
   if (OMax_intra < 0) then
      OMax_intra = OMax
   end if
   if (OMax_inter < 0) then
      OMax_inter = OMax
   end if
   if (OMin_parent < 0) then
      OMin_parent = OMax
   end if
   if (MaxTraj == 0) then
      MaxTraj = MaxTrajLimit
   end if

   call spawn_params%initialize(CSThresh=CSThresh, CFThresh=CFThresh, PopToSpawn=PopToSpawn, OMax_inter=OMax_inter, &
                                OMax_intra=OMax_intra, OMin_parent=OMin_parent, SOCThresh=SOCThresh, MaxTraj=MaxTraj, &
                                MultiSpawn=MultiSpawn, SpawnCoupV=SpawnCoupV)

!     Copy QM/MM variables into QM_MM_Module
   qcZQMMM = ZQMMM
   qcNumAS = NumAS
   qcNumS = NumS
   qcdConfineD = ConfineD
   qcdConfineK = ConfineK

!TM
   ! setup for dummy atoms in electronic structure code
   esNDummy = NDummyParticles
   esNPart = NumParticles
   esNMM = NumMM
   if (esNDummy >= 1) then
      allocate (esDummyWeight(esNPart, esNDummy))
      esDummyWeight(:, :) = DummyCoeff(1:esNPart, 1:esNDummy)
      ! also do some parameter check & scaling here
      call FMS_CheckDummy(esNPart, esNDummy, esDummyWeight)
   end if
!

   if (EquilTStep < 0.0) then
      indEquilTStep = TimeStep
   else
      indEquilTStep = EquilTStep
   end if
   if (gliModel == 9) qcNumMM = numas * nums
   if (ConfineK /= 0) qczConfine = .true.
   indGamma = indGamma * (1.d-15 / FsToAu) !convert from 1/s to 1/au
   if (gldCoupTimestep <= 0.) gldCoupTimeStep = gldTimeStep / 4.0d0
   if (gldMinTimeStep <= 0.) gldMinTimeStep = gldCoupTimeStep / 4.0d0

!     Check AvH variables
   if (glzAvH) then
      if (gliAvHNStates > NumStates) gliAvHNStates = NumStates
      do i = 1, gliAvHNStates
         if (gliAvHStates(i) <= 0) then
            call NameListFail('AvHStates', int(gliAvHStates(i), DefInt))
         end if
      end do
      write (fmiOut, *) 'Performing Average Hamiltonian Dynamics over states:', gliAvHStates(1:gliAvHNStates)
   end if

!
!     Make sure required variables have been set
!
   if (NumParticles <= 0) call NameListFail('NumParticles', NumParticles)
   if (NumStates <= 0) call NameListFail('NumStates', NumStates)
   if (initState <= 0 .or. InitState > NumStates) call NameListFail('InitState', InitState)
   if (.not. MinSearch) then
      if (SimulationTime < 0) call NameListFail('SimulationTime', int(SimulationTime, DefInt))
      if (TimeStep <= 0) call NameListFail('TimeStep', int(TimeStep, DefInt))
   end if

   if (zFatal) call FMS_DieError('Namelist input failed.')

   if (zDCD) then
      call FMS_DieError('DCD output has been removed. If you need it, please comment on this this GitHub issue:'//NL// &
                        'https://github.com/ispg-group/openfms/issues/33')
   end if

!
!     Make sure output options make sense
!

   if (fmzSOME .eqv. .true. .and. NumTriplets == 0) then
      write (fmiOut, *) 'SOME cannot be saved because number of triplet states is set to zero.'
      fmzSOME = .false.
      write (fmiOut, *) 'zSOME now set to false.'
   end if
   if (fmzSOCeff .eqv. .true. .and. NumTriplets == 0) then
      write (fmiOut, *) 'SOCeff cannot be saved because number of triplet states is set to zero.'
      fmzSOCeff = .false.
      write (fmiOut, *) 'zSOCeff now set to false.'
   end if

contains

   ! Open FMS.out for writing
   ! If we're restarting, we're append to an existing file,
   ! otherwise we overwrite it (if it somehow already exists).
   function open_fms_out(restart) result(output_unit)
      integer, intent(in) :: restart
      integer :: output_unit
      character(len=256) :: filepath

      filepath = trim(FMSWorkingDir)//'FMS.out'
      if (restart == 0) then
         open (newunit=output_unit, file=filepath, status='replace')
      else
         open (newunit=output_unit, file=filepath, position='append')
      end if
   end function open_fms_out

!
!     Return information about the bad input
!
   subroutine NameListFail(name, value)
      character(len=*), intent(in) :: name
      integer(kind=DefInt), intent(in) :: value
      write (error_unit, *) 'ERROR in Control.dat: namelist variable '//trim(name)
      write (error_unit, *) 'has invalid value: ', value
      write (error_unit, *) 'or has not been set.'
      zFatal = .true.
   end subroutine NameListFail

!
!     check the coefficients of dummy atom, and also scale them
!
   subroutine FMS_CheckDummy(NPart, NDummy, Weight)
      integer(kind=DefInt), intent(in) :: NPart, NDummy
      real(kind=DefReal), intent(inout) :: Weight(NPart, NDummy)
      integer(kind=DefInt) :: i, j
      real(kind=DefReal) :: Tot, toler

      toler = 1.0d-08

      do i = 1, NDummy
         Tot = sum(Weight(:, i))
         if (Tot < toler) then
            write (*, '(A,200f10.8)') 'Weight for defining dummy atom:', Weight(:, i)
            call FMS_DieError('STOP: Weight is set incorrectly for use in nDummyAtom.')
         end if
         do j = 1, NPart
            Weight(j, i) = Weight(j, i) / Tot
         end do
      end do

   end subroutine FMS_CheckDummy

end subroutine FMS_ReadNameList
