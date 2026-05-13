!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    @brief Parameters and methods for energy minimizations
!!
!!    For computational efficiency, all evaluations during the minimum
!!    search should be done using a single trajectory. We'll pass MSTraj
!!    around using MinSearchModule, which is only useful to the
!!    functions below
!<
module MinSearchModule
   use GlobalModule
   use TrajectoryModule
   use TrajectoryCalcsModule
   use TrajectoryIOModule, only: FMS_WriteFXYZ
   use QM_MM_Module, only: qczPCharge, qczQMMM, qcNumQM
   implicit none
   private
   public :: FMS_Minimizer
   save

   ! enumerate tyoe for the minimization
   integer(DefInt), public :: mnMinType
   integer(DefInt), parameter, public :: &
  &   Min_MM = 1, &
  &   Min_QM = 2, &
  &   Min_Both = 3, &
  &   Min_FixQM = 4 ! freeze the Cartesian of the QM part, mniStartPtcle = qcnumQM+1

   ! the trajectory sets up the calculations of the energy
   type(T_Trajectory), public :: mnT1

   integer(DefInt), public :: mnnStepToPrint

   integer(DefInt), private ::  &
  &   mnNumPartSearch, & !
  &   mniStartPtcle,   & ! this allows the QM part to be completely frozen
  &   mniStepNum

   real(DefReal) :: mndBondPenalty, mndAnglePenalty

   real(DefReal) :: penalty(4) ! array of penalty for the constaints used
   ! 1=bond, 2=angle, 3=dihedal, 4=CoM

   character(len=256), public :: mncFileName

contains

!
!-----References to the subroutines below get passed to dfpmin-----!
!

!>
!!    The function (in this case, the energy) to be minimized
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_MSfunc(P) result(pot)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use RattleModule
      use TrajectoryModule
      real(kind=DefReal) :: P(:)
      real(kind=DefReal) :: pot
      integer(kind=DefInt) :: i, j, nc, nt

!     Set trajectory to current geometry
      do i = 1, mnNumPartSearch
         j = i + mniStartPtcle
         call mnt1%set_pos(j, P(3 * i - 2:3 * i))
      end do

!     Calculate energy
      if (mnMinType == Min_MM .or. mnMinType == Min_FixQM) then
         pot = FMS_MMPot(mnT1)
      else
         pot = FMS_PotentialT(mnT1)
      end if

!     Add penalty for constrained degrees of freedom
      if (glzConstrain .or. mnMinType == Min_MM) then
         do nc = 1, nconstraint
            nt = cn_type_list(nc)
            pot = pot + 0.5d0 * penalty(nt) * constraint(mnt1, nc)**2
         end do
      end if
   end function FMS_MSFunc

!>
!!    Returns the gradient of the function to be minimized
!!    (i.e. -1*the forces)
!<
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_MSdfunc(P, grad)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use TrajectoryModule
      use RattleModule
      real(DefReal), intent(in) :: P(:)
      real(DefReal), intent(out) :: grad(size(P))

      real(DefReal) :: ForceTemp(3 * mnT1%NumParticles)
      integer(kind=DefInt) :: i, j, nc, nt

!     Calculate forces
      do i = 1, mnNumPartSearch
         j = i + mniStartPtcle ! apply QM offset
         call mnT1%set_pos(j, P(3 * i - 2:3 * i))
      end do

      if (mnMinType == Min_MM .or. mnMinType == Min_FixQM) then
         ForceTemp = FMS_MMForces(mnT1)
      else
         ForceTemp = FMS_Forces(mnT1)
      end if
!
!     Add penalty for constrained degrees of freedom
!
      if (glzConstrain .or. mnMinType == MIN_MM) then
! Sum forces on constrained bonds
         do nc = 1, nconstraint
            nt = cn_type_list(nc)
            ForceTemp = ForceTemp - penalty(nt) * constraint(mnt1, nc) * D_constraint(mnt1, nc)
            ! note here that that derivative of the constraint is SUBTRACTED
            ! from the derivative and it was ADDED in the function value
         end do
      end if

!     Return the force on the minimized subsystem
      grad = -1.0d0 * ForceTemp(3 * mniStartPtcle + 1:3 * (mniStartPtcle + mnNumPartSearch))

   end subroutine FMS_MSdFunc

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   function FMS_Minimizer(iter, toler, nChargeCycle, StepNum)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !>
      !!    Setup and call an energy  minimization for mnT1.
      !!    Parameters for minimization are read from
      !!    read from MinSearchModule
      !!    \param iter Maximum number of iterations
      !!    \param toler Tolerance for convergence
      !!    \param type System to minimize
      !<
      use TrajectoryModule
      use ElecStrucModule, only: eszPartialCharges
      use RattleModule
      use QM_MM_Module
      integer(kind=DefInt) :: iter, FMS_Minimizer
      real(kind=DefReal) :: toler
      integer(kind=DefInt), optional :: nChargeCycle, StepNum

      integer(kind=DefInt) :: QMCycleSteps, iCycle, iTry
      integer(kind=DefInt) :: i, j, nc
      integer(kind=DefInt) :: iterconst, stepsave
      real(kind=DefReal) :: FirstEnergy, fmin !, PotEn(mnT1%NumStates)
      logical :: MErr

      character(len=2) :: ctype ! code for the type of minimization performed

      real(kind=DefReal), allocatable :: posvec(:)

!     Cycle for recalculating QM point charges
      if (present(nChargeCycle)) then
         QMCycleSteps = nChargeCycle
      else
         QMCycleSteps = 0
      end if

!     For restarts, this is the step to start on
      if (present(StepNum)) then
         mniStepNum = StepNum
      else
         mniStepNum = 0
      end if

!     Set up constraints, if necessary
      penalty = [1.5, 0.3, 0.3, 1.0]
      mndBondPenalty = 1.5d0
      mndAnglePenalty = 0.3d0

!     if(glzConstrain .and. .not. cnzFileRead)
      if (glzConstrain) call Rattle_ReadConstraints(qcNumQM)
      do nc = 1, nconstraint
         write (fmiOut, '(i1,":",4(1x,i3))') cn_type_list(nc), cn_atom_list(:, nc)
      end do

      !if(mnMinType==Min_MM) call FMS_FreezeQM(mnT1)

      if (glzConstrain .or. mnMinType == Min_MM) call Rattle_SetConstraints(mnT1)

!     Write minimization info

      select case (mnMinType)
      case (Min_MM)
         mniStartPtcle = 0
         mnNumPartSearch = mnT1%NumParticles
         ctype = 'MM'
         qczPCharge = .true.

         if (eszPartialCharges) then
            qczPCharge = .true.
            call mnT1%geom_changed()
         end if

      case (Min_QM)
         mniStartPtcle = 0
         mnNumPartSearch = qcNumQM
         ctype = 'QM'
         qczPCharge = .false.

      case (Min_Both)
         mniStartPtcle = 0
         mnNumPartSearch = mnT1%NumParticles
         ctype = 'mn'
         qczPCharge = .false.

      case (Min_FixQM)

         mniStartPtcle = qcNumQM
         mnNumPartSearch = qcNumMM
         ctype = "Fx"
         qczPCharge = .true.

      case default
         call FMS_DieError('ERROR in FMS_Minimizer: unknown minimization type.')

      end select

4343  format('Minimizing ', a2, ' system for ', i5, ' steps.')
      write (fmiOut, 4343) ctype, iter
      flush (fmiOut)

!    Copy positions into an array and setup output
      allocate (PosVec(3 * mnNumPartSearch))
      do i = 1, mnNumPartSearch
         j = i + mniStartPtcle
         PosVec(3 * i - 2:3 * i) = mnT1%get_pos(j)
      end do

      if (mnMinType == Min_Both) then
         glzMinSearch = .false. !calculate everything, like real dynamics
      else
         glzMinSearch = .true. !do not calculate electronic coordinate
         !forces for MM system
      end if

      FirstEnergy = FMS_MSFunc(PosVec)
!      FirstEnergy=FMS_PotentialT(mnT1)
      write (fmiOut, *) 'Initial energy (Hartrees): ', FirstEnergy
      flush (fmiOut)

      if (qczPCharge .and. QMCycleSteps > 0) then

!     Run minimization, recalculating QM partial charges every N steps
         do iCycle = 1, iter / QMCycleSteps

            qczPCharge = .false.
            FirstEnergy = FMS_PotentialT(mnT1)
            qczPCharge = .true.
            call mnT1%geom_changed()

            call FMS_dfpmin(posvec, 3 * mnNumPartSearch, toler, QMCycleSteps, fmin, &
     &            FMS_MSFunc, FMS_MSDFunc, OutFunc=FMS_MSWrite)
            if (iter < iter / QMCycleSteps .and. iCycle > 1) exit
         end do

      else

!     Run minimization
         call FMS_dfpmin(posvec, 3 * mnNumPartSearch, toler, iter, fmin, &
     &         FMS_MSFunc, FMS_MSDFunc, OutFunc=FMS_MSWrite)
      end if

!     Make sure constraints are satisfied
      if (glzConstrain .or. mnMinType == Min_MM) then
         iTry = 0
         do while (.not. all_position_constrained(mnT1))
            iTry = iTry + 1
            if (iTry == 21) then
               call FMS_DieError('ERROR in Minimizer: Constraints not satisfied.')
            end if
            IterConst = 100
2021        format('Constraints not satisfied. Tightening penalty functions for ', i4, ' iterations.')
            write (fmiOut, 2021) IterConst

            mndBondPenalty = mndBondPenalty * 100.0d0
            mndAnglePenalty = mndAnglePenalty * 100.0d0

            penalty = 100.d0 * penalty

            call mnT1%geom_changed()
            call FMS_dfpmin(posvec, 3 * mnNumPartSearch, toler &
                            , IterConst, fmin, FMS_MSFunc, FMS_MSDFunc, OutFunc=FMS_MSWrite)
         end do
      end if

!     Report the results
      FMS_Minimizer = iter
      write (fmiOut, *) 'Minimization ran for ', iter, ' steps'
1011  format('Energy reduced from ', f16.10, ' to ', f16.10)
      if (mnMinType == Min_MM .or. mnMinType == Min_FixQM) then
         !  call FMS_H5DReadStep('MMPot',1,mnT1%h5Group, FirstEnergy)
      else
         !  call FMS_H5DReadStep('PotEn',1,mnT1%h5Group,PotEn)
         FirstEnergy = mnT1%ElecStruc%PotEn(mnT1%StateID)
      end if
      write (fmiOut, 1011) FirstEnergy, FMin
      flush (fmiOut)

!     Make sure to write last step
      if (mod(mniStepNum - 1, mnnStepToPrint) /= 0) then
         StepSave = mnnStepToPrint
         mnnStepToPrint = 1
         MErr = FMS_MSWrite(PosVec, iter, fmin)
         mnnStepToPrint = StepSave
      end if

!     Switch off point charge calculations, if necessary
      if (qczPCharge) then
         qczPCharge = .false.
         call mnT1%geom_changed()
      end if

      if (mnMinType == MIN_MM) call Rattle_ReadConstraints()
      deallocate (PosVec)
      MErr = .true.

   end function FMS_Minimizer

!>
!!    Writes XYZ file for current geometry in minsearch.
!<
   function FMS_MSWrite(p, its, fret) result(MErr)
      use TrajectoryModule
      use ElecStrucModule
      real(kind=DefReal) :: p(:), fret
      character(len=64) :: comment
      logical :: MErr
      integer(kind=DefInt) :: jPtcle, its
      integer(kind=DefInt) :: iPtcle, iState, iOrb, jOrb
      real(kind=DefReal), allocatable, save :: DimArray(:, :), ForceVec(:)
      real(kind=DefReal), allocatable, save :: PotEn(:), d(:)
      real(kind=Defreal), allocatable, save :: Orbitals(:, :), CIVec(:, :)
      real(kind=DefReal) :: rms
      character(len=256) :: cTemp

      integer(kind=DefInt), save :: LastPSize = 0

      MErr = .true.

!     Allocate data buffers
      if (.not. allocated(DimArray)) allocate (DimArray(mnT1%NumParticles, 3))
      if (.not. allocated(PotEn)) allocate (PotEn(mnT1%NumStates))
      if (.not. allocated(Orbitals)) allocate (Orbitals(esnBasis, esnBasis))
      if (.not. allocated(CIVec)) allocate (CIVec(mnT1%NumStates, eslCIVec))
      if (.not. allocated(ForceVec)) allocate (ForceVEc(3 * mnT1%NumParticles))
      if (size(P) /= LastPSize) then
         if (allocated(d)) deallocate (d)
         allocate (d(size(P)))
         LastPSize = size(P)
      end if

!     Write to FMS.out
      call FMS_MSDFunc(p, d)

      rms = sqrt(sum(d**2) / mnT1%NumParticles)
1000  format('Iteration: ', i8, ' energy: ', e14.6, ' gradient rms:', e14.6)
      write (fmiOut, 1000) its, fret, rms
      flush (fmiOut)

!
!     Only write output every mnnStepToPrint steps
!
      mniStepNum = mniStepNum + 1
      if (mod(mniStepNum - 1, mnnStepToPrint) /= 0) return

!     Set position
      do IPtcle = 1, mnNumPartSearch
         jPtcle = iPtcle + mniStartPtcle
         call mnT1%set_pos(jPtcle, P(3 * IPtcle - 2:3 * IPtcle))
      end do

!     Write positions
      write (comment, *) 'Optimization step ', its
      call FMS_WriteFXYZ(mnT1, 'MinSearch.xyz', comment)
      do iPtcle = 1, mnT1%NumParticles
         DimArray(iPtcle, 1:3) = mnT1%Particle(iPtcle)%get_pos()
      end do

   end function FMS_MSWrite

!     TODO: We need to re-implement this function(originally taken from NR)
!!    The BFGS minimization algorithm,
!!    plus all the other subroutines needed to run it. Note
!!    that "func" and "dfunc" here are function references, not
!!    variables.
!!
!!    Everythere here is from NR 2.10. I modified it to use FMS
!!    DefTypes, and also made the arrays automatic, since there's no
!!    need for static allocation - this thing will only get called once
!!    or twice.
!!
!!    Will write parameters after each step, if you pass it a reference
!!    to a routine to write the output
!!
!!    \param p       Coordinate vector
!!    \param n       Number of components of p
!!    \param gtol    tolerance
!!    \param iter    Pass in the max number of iterations.
!!                Returns the actual number of iterations
!!    \param fret    Returns final value of function
!!    \param func    This is a reference to an external function which
!!                calculates the quantity to be minimized
!!    \param dfunc    This is a reference to an external function which
!!                returns the gradient of func
!!    \param Outfunc  A routine to write output for each minimization step.
!!                It's optional, and there will be no output if it's not provided.
!!
!!    The interfaces of fun,dfunc and OutFunc must conform to those
!!    in the interface block below
!<
   subroutine FMS_dfpmin(p, n, gtol, iter, fret, func, dfunc, OutFunc)
      integer(kind=defInt) :: iter, n
      real(kind=defReal) :: p(n), gtol, fret

      interface
         function func(p) result(FuncVal)
            use GlobalModule, only: DefReal
            real(kind=DefReal) :: p(:), FuncVal
         end function func

         subroutine dfunc(p, dFuncVec)
            use GlobalModule, only: DefReal
            real(kind=DefReal), intent(in) :: P(:)
            real(kind=DefReal), intent(out) :: dFuncVec(size(p))
         end subroutine dfunc

         function OutFunc(p, its, fret) result(Merr)
            use GlobalModule, only: defInt, defReal
            real(kind=DefReal) :: p(:), fret
            integer(kind=DefInt) :: ITS
            logical :: Merr
         end function OutFunc
      end interface

      call FMS_DieError('ERROR: dfpmin function not implemented')
   end subroutine FMS_dfpmin

end module MinSearchModule
