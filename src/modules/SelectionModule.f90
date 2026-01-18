! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module SelectionModule
use GlobalModule
use BundleModule
use TrajectoryModule
use BundleCalcsModule, only: FMS_bH, FMS_Norm
use OverlapModule, only: overlap
use RandomModule, only: fms_ranb
implicit none

private
public :: print_stochastic_selection_params
public :: FMS_StochasticCollapse, FMS_CalculateSelectionTime

! The standard overlap threshold is 1/e, and allows us to differentiate
! betweeen normal and premature selections, where the latter are defined
! as selections for which the parent-child TBF pair overlap is above
! 1/e of its initial value at spawning.
real(kind=DefReal), parameter :: swissThresh=3.678794411714d-1

contains

subroutine print_stochastic_selection_params()

if (glzStoSwiss) then

    write(fmiOut, *) 'Stochastic selection: Running in AIMSWISS mode!'

    write(fmiOut, *) 'SWISS mode sets DecoherenceTime to ', gldDecoherenceTime
    write(fmiOut, '(a,i0)') ' and IgnoreState to ', glIgnoreState
    write(fmiOut, *) 'so that Trajectories are only killed when they cannot spawn anymore.'

else

    if (glzStoOlap) then
        write (fmiOut, *) 'Stochastic selection: Running in OSSAIMS mode!'
        write (fmiOut, *) 'OSSAIMS threshold (overlap between TBFs):', gldStochaThresh
    else
        write (fmiOut, *) 'Stochastic selection: Running in ESSAIMS mode!'
        write (fmiOut, *) 'ESSAIMS threshold (coupling between TBFs):', gldStochaThresh
    end if

end if

end subroutine print_stochastic_selection_params

!!    @brief Main driver for stochastic selection algorithm
!!
!!    This subroutine is the public facing interface to the stochastic selection
!!    algorithm, independent of which variant of stochastic selection is used.
!!    It is called within the FMS_Dynamics() subroutine
!!
!!    @ingroup propagation
!!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_StochasticCollapse(B1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! FMS Bundle (molecular wavefunction) that will be modified
type(T_TrajectoryBundle), intent(inout) :: B1
! An array of temporary Bundles for single state stochastic selection
type(T_TrajectoryBundle), allocatable  :: BundleSS(:)
! AIMSWISS: Shall we perform stochastic selection
logical :: performSelection
! AIMSWISS: Current selection time
real (kind=DefReal) :: selectionTime
integer (kind=DefInt) :: iState


if (B1%NumTraj == 1) return

! Check that we are not just before a Spawning event
if ((gldLastSpawnSto>=B1%CurrentTime) &  ! When AIMSWISS is invoked don't
      .and. (.not. glzStoSwiss))then       ! wait for next Spawning event
    write(fmiout,*) 'Aborting Stochastic Collapse - Spawn is coming...'
    return
endif

if (glzStoSwiss) then
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!   Check if a selection time of any pair TBF is reached
    call FMS_CheckSelectionTime(B1, performSelection, selectionTime)
!   If not then there is nothing to do
    if (.not. performSelection) then
        return
    endif
endif

! State-Specific Stochastic Selection (aka 4S :-)
! 1. Create temporary Bundles that hold a set of trajectories for a given state
!    We essentically split the Bundle into NumStates smaller bundles
! 2. Perform selections within those bundles
! 3. Convert the results into the original bundle.
if (glzStoStateSpecific) then

    allocate(BundleSS(B1%NumStates))
    call fill_state_bundles(B1, BundleSS)

    do iState = 1, B1%NumStates
        if (BundleSS(iState)%NumTraj == 1) cycle
        call perform_stochastic_selection(BundleSS(iState), selectionTime)
    enddo

    call copy_state_bundles_to_original_bundle(B1, BundleSS)

    do iState = 1, B1%NumStates
        call BundleSS(iState)%destroy()
    enddo
    deallocate(BundleSS)

else

    call perform_stochastic_selection(B1, selectionTime)

endif

end subroutine FMS_StochasticCollapse

!!    @brief Implements stochastic selection algorithm
!!
!!    Partition nuclear trajectory basis functions into uncoupled blocks
!!    then stochastically collapse the wavefunction onto a single block.
!!    This is done by renormalizing population in collapsed block to one
!!    and then marking all other trajectories to be killed.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine perform_stochastic_selection(B1, selectionTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(inout) :: B1
! Coupling matrix for TBF basis
integer (kind=DefInt), dimension(B1%NumTraj,B1%NumTraj) :: Coupled
! Number of blocks
integer (kind=DefInt) :: nblock
! Matrix containing IDs of TBFs belonging to a block.
integer (kind=DefInt), dimension(B1%NumTraj+1,B1%NumTraj+1) :: blocktrajid
! Array containing number of TBFs belonging to a block.
integer (kind=DefInt), dimension(B1%NumTraj+1) :: ntrajblock
! Population of all blocks combined (not necessarily 1)
real (kind=DefReal) :: totBlockpop
! Array containing the populations of each block
real (kind=DefReal),dimension(B1%NumTraj+1) :: Blockpop
! Index of selected block
integer (kind=DefInt) :: iblockslct
! AIMSWISS: Current selection time
real (kind=DefReal), intent(in) :: selectionTime


! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! First, decompose the hamiltonian into a block diagonal representation
! work out the coupling matrix
Coupled = 0
call FMS_BuildCoupled(B1, Coupled, selectionTime)


! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Next, group trajectories into coupled blocks
blocktrajid =0
ntrajblock =0
nblock = 0
call FMS_GroupIntoBlocks(B1, Coupled, blocktrajid, ntrajblock, nblock)

! There should be at least one block of trajectories
if (nblock==0) then
    write(fmiOut,*)'Number of trajectory basis blocks is zero'
    call FMS_DieError("ERROR in FMS_StochasticCollapse")
endif
! If just a single block then no more work to be done
if (nblock==1) return
! At this point, more than one uncoupled block of trajectories,
! so we want to stochastically collapse to one of them

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Compute populations of each block
totBlockpop=0.0d0
Blockpop = 0
call FMS_ComputeBlockPopulations(B1, blocktrajid, ntrajblock, nblock, &
                                 totBlockpop, Blockpop)

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Select block randomly based on block populations via Monte Carlo procedure
call FMS_SelectBlock(B1%NumTraj, nblock, ntrajblock, totBlockpop, Blockpop, &
                     iblockslct)


! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Renormalize TBF amplitudes of selected block
call FMS_RenormalizeBlockAmplitudes(B1, Blockpop(iblockslct), &
                                    ntrajblock(iblockslct),   &
                                    blocktrajid(:,iblockslct))

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Remove all TBFs that do not belong selected block
call FMS_KillOtherBlocks(B1, nblock, ntrajblock, iblockslct, blocktrajid)

end subroutine perform_stochastic_selection


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_CalculateSelectionTime(parent_s, child_s, child_i)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Specific to AIMSWISS.
! Called by SpawnModule:FMS_Spawn subroutine
!
!     parent_s and child_s are the parent and child TBFs at the spawning time,
!     whose forces are used to calculate the selection time. Their state is not
!     changed in this subroutine!
!
!     child_i is the child TBF at the entry time, which will be added to the
!     Bundle object after this subroutine is finished. Its state is modified
!     in this subroutine, by setting its T1%SWISS%SelectionTime variable to
!     its current time + the decoherence time with respect to its parent.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type(T_Trajectory), intent(in) :: parent_s, child_s
type(T_Trajectory), intent(inout) :: child_i
real(kind=DefReal) :: decoherenceTime

decoherenceTime = FMS_CalculateDecoherenceTime(parent_s, child_s, parent_s%NumParticles)
child_i%SWISS%SelectionTime = child_i%SWISS%BirthDate + decoherenceTime
child_i%SWISS%ParentOverlap = abs(overlap(parent_s, child_s)) ** 2

write(fmiOut,'(a,i0,a,f0.2)') 'SWISS: Trajectory ', parent_s%TrajID, ' and '// &
                'its child will decohere at t = ', child_i%SWISS%SelectionTime
write(fmiOut,'(a,f5.3)') 'Their current absolute overlap is ', child_i%SWISS%ParentOverlap

end subroutine FMS_CalculateSelectionTime

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function FMS_CalculateDecoherenceTime(parent, child, npart) result(decoherenceTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Specific to AIMSWISS.
! Called by FMS_CalculateSelectionTime:
!
!   Caclculate the decoherence time (tau_D) in the following way:
!
!       First, we calculate the decoherence rate (Gamma_D) rate via
!
!          Gamma_D = (F_P - F_C)^T * alpha^-1 * (F_P - F_C) / 4,
!
!       where F_P and F_C are the gradients acting on the parent and child TBF,
!       respectively, and alpha is a matrix containing the widths of the TBFs.
!
!       Second, we determine tau_D by taking the square root of Gamma_D and
!       its reciprocal
!
!           tau_D = 1. / sqrt(Gamma_D)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory), intent(in)   :: parent, child
integer(kind=DefInt), intent(in) :: npart
real(kind=DefReal) :: decoherenceTime
integer(kind=DefInt) :: ipart, jdim, ishifted, istart, ndim, &
                        parentSt, childSt
real(kind=DefReal) ::  forceElementDiff2, decoherenceRate2, &
                       parentForceElement, childForceElement, &
                       summand2, width

! Calculation of the squared decoherence rate
parentSt = parent%StateID
childSt = child%StateID
! stochastic selection is currently only implemented for
! multi-state spawning (may change in the future)
if (childSt == parentSt) then
    write(fmiOut,*) "SWISS: States don't differ, decoherence not possible"
    call FMS_DieError("ERROR in FMS_CalculateDecoherenceTime")
endif

decoherenceRate2 = 0.d0
forceElementDiff2 = 0.d0
do ipart = 1, npart
    width = parent%Particle(ipart)%Width
    ndim = parent%Particle(ipart)%NumDimensions

    istart = (ipart - 1) * ndim

    do jdim = 1, ndim
        ishifted = jdim + istart
        parentForceElement = parent%ElecStruc%DerivMat(parentSt, &
                             parentSt, ishifted)
        childForceElement = child%ElecStruc%DerivMat(childSt, &
                            childSt, ishifted)
        forceElementDiff2 = (parentForceElement - &
                             childForceElement) ** 2
        summand2 = forceElementDiff2 / (4.d0*width)
        decoherenceRate2 = decoherenceRate2 + summand2
    enddo
enddo

decoherenceTime = 1.d0 / sqrt(decoherenceRate2)

return

end function FMS_CalculateDecoherenceTime

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_CheckSelectionTime(B1, performSelection, selectionTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Specific to AIMSWISS.
! Called by FMS_StochasticCollapse:
!
!    Check if parent-child TBF pair has decohered, and if so
!    change the module variable performSelection to true.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(in) :: B1
logical, intent(out) :: performSelection
real(kind=DefReal), intent(out) :: selectionTime
integer(kind=DefInt) :: i, ntraj

performSelection = .false.
selectionTime = 0.d0
ntraj = B1%NumTraj
do i = 2, ntraj
    if ( B1%Trajectory(i)%SWISS%SelectionTime <= B1%CurrentTime ) then
        performSelection = .true.
        selectionTime = B1%Trajectory(i)%SWISS%SelectionTime
        write(fmiout,*)'SWISS: Time to select! The current time is ', &
                        B1%CurrentTime, ' and the selection time is ', &
                        selectionTime
        exit
    endif
enddo

end subroutine FMS_CheckSelectionTime

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_BuildCoupled(B1, Coupled, selectionTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Called by FMS_StochasticCollapse:
!
!    Procedure handles the setup of the coupling matrix. The idea
!    behind the use of a coupling matrix is that the TBF basis can
!    be thought of as a graph, where the vertices are the TBFs
!    and the edges are the overlap or matrix elements between TBFs.
!    To determine the adjacency matrix of this graph, a.k.a. the
!    coupling matrix, we either use user defined thresholds
!    (as in O/ESSAIMS) or an adaptive procedure (a in AIMSWISS).
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(inout) :: B1
integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
                 intent(inout) :: Coupled
real(kind=DefReal), intent(in) :: selectionTime
integer(kind=DefInt) :: ntraj, i, j

ntraj   = B1%NumTraj

do i = 1, ntraj
    Coupled(i, i) = 1
    do j = i, ntraj
       !  X.Z. Keep different trajectories in the same CBF together.
       if(B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf)then
           Coupled(i,j)=1
           Coupled(j,i)=1
           cycle
       endif
    enddo
enddo


if (glzStoSwiss) then

    call FMS_BuildCoupled_SWISS(B1, Coupled, selectionTime)

else if (glzStoOlap) then

    call FMS_BuildCoupled_OSS(B1, Coupled)

else

    call FMS_BuildCoupled_ESS(B1, Coupled)

end if
end subroutine FMS_BuildCoupled


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_BuildCoupled_ESS(B1, Coupled)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(in) :: B1
integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
                 intent(inout) :: Coupled
integer(kind=DefInt) :: ntraj, i, j

ntraj   = B1%NumTraj

do i = 2, ntraj
    do j = 1 ,i
       if(B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf) cycle
       if(abs(FMS_bH(B1,i,j)) > gldStochaThresh) then
           Coupled(i,j) = 1
           Coupled(j,i) = 1
       endif
    enddo
enddo

! Uncover all indirect connections between TBFs
call FMS_ConvergeCoupled(ntraj, Coupled)

end subroutine FMS_BuildCoupled_ESS


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_BuildCoupled_OSS(B1, Coupled)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(in) :: B1
integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
                 intent(inout) :: Coupled

! Overlap between TBF i and j
complex (kind=DefComp) :: S_ij
integer(kind=DefInt) :: ntraj, i, j

ntraj   = B1%NumTraj

do i = 2, ntraj
    do j = 1 ,i
       if(B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf) cycle
       S_ij = overlap(B1%Trajectory(i), &
                                 B1%Trajectory(j))
       if(abs(S_ij) > gldStochaThresh) then
           Coupled(i,j) = 1
           Coupled(j,i) = 1
       endif
    enddo
enddo

! Uncover all indirect connections between TBFs
call FMS_ConvergeCoupled(ntraj, Coupled)

end subroutine FMS_BuildCoupled_OSS


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_BuildCoupled_SWISS(B1, Coupled, selectionTime)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(inout) :: B1
integer(kind=DefInt), intent(inout) :: Coupled(B1%NumTraj, B1%NumTraj)
real(kind=DefReal), intent(in) :: selectionTime
integer(kind=DefInt) :: ntraj, i, j
real(kind=DefReal) :: S_ik, S_ik_max
logical       :: success

ntraj   = B1%NumTraj
success = .true.
do i = 2, ntraj
    do j = 1, i
        if (B1%Trajectory(i)%ParentID == B1%Trajectory(j)%TrajID) then
            ! Add connection for all parent-child TBF pairs that haven't
            ! decohered!
            if (.not.(abs(B1%Trajectory(i)%SWISS%SelectionTime &
                          - selectionTime) <= FPZero)) then
                Coupled(i, j) = 1
                Coupled(j, i) = 1
            else
                S_ik = abs(overlap(B1%Trajectory(i), B1%Trajectory(j))) ** 2
                S_ik_max = swissThresh * B1%Trajectory(i)%SWISS%ParentOverlap
                write(fmiout,'(a,i0,a,i0,a)') 'Check if overlap of the pair ', &
                                 B1%Trajectory(i)%TrajID, ' and ', B1%Trajectory(j)%TrajID, &
                                 ' is not too large'
                if(S_ik > S_ik_max) then
                    success = .false.
                endif
                call FMS_WriteSelectionLog(B1%Trajectory(j), &
                                           B1%Trajectory(i), &
                                           S_ik, S_ik_max,   &
                                           success)
                success = .true.
            endif
        endif
    enddo
enddo

! Uncover all indirect connections between TBFs
call FMS_ConvergeCoupled(ntraj, Coupled)

end subroutine FMS_BuildCoupled_SWISS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_ConvergeCoupled(ntraj, Coupled)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     This matrix has the properties that it is non-zero is
!        Coup^1   : directy connected
!        Coup^2   : connected by 1 or less common trajectory
!        Coup^3   : connected by 2 or less common trajectories
!     We will iterate the matrix multiplication to convergence
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer(kind=DefInt), intent(in) :: ntraj
integer(kind=DefInt), dimension(ntraj,ntraj), &
                 intent(inout) :: Coupled
! Coupling matrix for TBF basis of the previous iteration
integer (kind=DefInt), dimension(ntraj,ntraj) :: Coupled_prev
do
    Coupled_prev = Coupled

    Coupled = min( matmul(Coupled,Coupled), 1 )

    if( all(Coupled==Coupled_prev) ) exit
enddo

end subroutine  FMS_ConvergeCoupled


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_GroupIntoBlocks(B1, Coupled, blocktrajid, &
                               ntrajblock, nblock)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_TrajectoryBundle), intent(in) :: B1
integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
                 intent(in) :: Coupled
integer (kind=DefInt), dimension(B1%NumTraj+1,B1%NumTraj+1),  &
                       intent(inout) :: blocktrajid
integer (kind=DefInt), dimension(B1%NumTraj+1), &
                       intent(inout) :: ntrajblock
integer (kind=DefInt), intent(inout) :: nblock

! Boolean array containing IDs of already considered TBFs
logical, dimension(B1%NumTraj+1) :: ztrajdone

integer(kind=DefInt) :: ntraj, i, j, l, itraj

ntraj = B1%NumTraj

ztrajdone(:)=.false.

! Loop over ntraj possible blocks
do i=1,ntraj+1

!   Add the first unsorted trajectory to this current block (for first block,
!   will always be first trajectory)
    do j=1,ntraj
        if (ztrajdone(j)) cycle
        blocktrajid(1,i)=j
        ztrajdone(j)=.true.
        exit
    enddo

    itraj=blocktrajid(1,i) !first trajectory in block

!   If no trajectories added, then exit loop
    if (itraj==0) then
        nblock=i-1
        exit
    endif


!   Add other coupled trajectories to this current block
    l=1
    do j=1,ntraj
        if (ztrajdone(j)) cycle !if this trajectory has already been sorted skip
        if (Coupled(j,itraj)==1) then
            l=l+1
            blocktrajid(l,i)=j
            ztrajdone(j)=.true.
        endif
    enddo
    ntrajblock(i)=l

enddo

end subroutine FMS_GroupIntoBlocks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_ComputeBlockPopulations(B1, blocktrajid, ntrajblock, &
                                       nblock, totBlockpop, Blockpop)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type(T_TrajectoryBundle), intent(in) :: B1
integer (kind=DefInt), dimension(B1%NumTraj+1,B1%NumTraj+1),  &
                       intent(in) :: blocktrajid
integer (kind=DefInt), dimension(B1%NumTraj+1), &
                       intent(in) :: ntrajblock
integer (kind=DefInt), intent(in) :: nblock
real (kind=DefReal), intent(inout) :: totBlockpop
real (kind=DefReal), dimension(B1%NumTraj+1), &
                     intent(inout) :: Blockpop
! Temporaty FMS Bundle needed to store blocks
type(T_TrajectoryBundle) :: BTemp

integer(DefInt) :: i, j, jtraj, nstate, npart, ncbf

nstate = B1%NumStates
npart = B1%Trajectory(1)%NumParticles
ncbf = B1%NCBFs

! Compute populations of each block:
!    Do this by building a temporary bundle of each block then use existing subroutine norm(bundle)
do i=1,nblock
    call BTemp%create(numtraj=ntrajblock(i), &
                      numdeadtraj=0, &
                      numstates=nstate, &
                      numparticles=npart, &
                      ncbfs=ncbf)

    do j=1,ntrajblock(i)
        jtraj=blocktrajid(j,i)
        BTemp%Trajectory(j)=B1%Trajectory(jtraj)
    enddo

    Blockpop(i)=FMS_Norm(BTemp)

    call BTemp%destroy()

    totBlockpop=totBlockpop+Blockpop(i)
enddo

end subroutine FMS_ComputeBlockPopulations

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_SelectBlock(ntraj,  nblock, ntrajblock, totBlockpop, &
                           Blockpop, iblockslct)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer (kind=DefInt), intent(in) :: ntraj
integer (kind=DefInt), intent(in) :: nblock
real (kind=DefReal), intent(in) :: totBlockpop
real (kind=DefReal), dimension(ntraj+1), &
                     intent(in) :: Blockpop
integer (kind=DefInt), dimension(ntraj+1), &
                       intent(in) :: ntrajblock
integer (kind=DefInt), intent(inout) :: iblockslct
! Character denoting selected block
character(len=1) :: cselect
! Random number used for Monte Carlo procedure
real (kind=DefReal) :: drand
! Dummy variable for storing relative population of block
real (kind=DefReal) :: Blockpopnorm
! Cumulative population of blocks used in Monte Carlo procedure
real (kind=DefReal) :: Blockpopnormsum

integer (kind=DefInt) :: i

!     Use random number to select which block to collapse to
drand=fms_ranb(i4zero)
!       find which block was randomly selected
Blockpopnormsum=0.0d0
do i=1,nblock
    Blockpopnorm=Blockpop(i)/totBlockpop
    Blockpopnormsum=Blockpopnormsum+Blockpopnorm
    if (drand<=Blockpopnormsum) then
        iblockslct=i
        exit
    endif
enddo

write(fmiOut,*)'Stochastically collapsing'
write(fmiOut,*)'Block   NTraj   Pop     Selected'
do i=1,nblock
    cselect=''
    if (i==iblockslct) cselect='*'
    write(fmiOut,'(2(I3,5X),f8.4,X,A1)')i,ntrajblock(i),Blockpop(i), &
          cselect
enddo

end subroutine FMS_SelectBlock


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_RenormalizeBlockAmplitudes(B1, slctBlockpop, slctNtrajblock,   &
                                          slctBlocktrajid)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type(T_TrajectoryBundle), intent(inout) :: B1

real (kind=DefReal), intent(in) :: slctBlockpop
integer (kind=DefInt), intent(in) :: slctNtrajblock
integer (kind=DefInt), dimension(B1%NumTraj+1), &
                       intent(in) :: slctBlocktrajid
! Norm of molecular wavefunction before selection
real (kind=DefReal) :: NormInit

! Renormalization factor
complex (kind=DefComp) ::cnorm

integer (kind=DefInt) :: j, jtraj

!     Renormalize selected block to total initial Norm (need not be 1)
NormInit=FMS_Norm(B1)
cnorm=dcmplx(sqrt(NormInit/slctBlockpop))
do j=1,slctNtrajblock
    jtraj=slctBlocktrajid(j)
    B1%Trajectory(jtraj)%Amplitude=B1%Trajectory(jtraj)%Amplitude*cnorm
enddo

end subroutine FMS_RenormalizeBlockAmplitudes

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_KillOtherBlocks(B1, nblock, ntrajblock, iblockslct, blocktrajid)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type(T_TrajectoryBundle), intent(inout) :: B1
integer (kind=DefInt), intent(in) :: nblock
integer (kind=DefInt), dimension(B1%NumTraj+1), &
                       intent(in) :: ntrajblock
integer (kind=DefInt), intent(in) :: iblockslct
integer (kind=DefInt), dimension(B1%NumTraj+1,B1%NumTraj+1),  &
                       intent(in) :: blocktrajid

integer (kind=DefInt) :: i, j, l, jtraj
! Mark all other trajectories for death
l=0
do i=1,nblock
    if (i==iblockslct) cycle
    do j=1,ntrajblock(i)
        l=l+1
        jtraj=blocktrajid(j,i)
        gliForceKill(l)=B1%Trajectory(jtraj)%TrajID
        B1%Trajectory(jtraj)%Amplitude=dcmplx(0.0d0,0.0d0)
    enddo
enddo


! Remove any dead trajectory amplitudes
do jtraj=1,B1%NumDeadTraj
    B1%DeadTraj(jtraj)%Amplitude=dcmplx(0.0d0,0.0d0)
enddo

end subroutine FMS_KillOtherBlocks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine FMS_WriteSelectionLog(parent, child, currentS, predictedS, l_select)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Specific to AIMSWISS.
! Called by FMS_BuildCoupled_SWISS:
!    Writer subroutine for selections, printing information about the selection
!    event. This is done in such a way that when the actual overlap is larger
!    than the predicted one, then the just mentioned information is printed
!    to FailSelect.log, while when this is not the case it is written to
!    Select.log. These files may then be used to assess the accuracy of AIMSWISS,
!    as done in Lassmann et al. J. Phys. Chem. Lett. 2022, 13, 51, 12011–12018,
!    https://doi.org/10.1021/acs.jpclett.2c03295.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(T_Trajectory), intent(in) :: parent, child
logical, intent(in) :: l_select
real(kind=DefReal), intent(in) :: currentS, predictedS
character(len=256) :: file_name
real(kind=DefReal) :: currentTime, locSelectionTime
integer(kind=DefInt ) :: id_c, state_c, &
                         id_p, state_p
integer(DefInt) :: iunit
logical :: file_exists

currentTime   = parent%get_time()
locSelectionTime = child%SWISS%SelectionTime

id_p       = parent%TrajID
state_p    = parent%StateID
id_c    = child%TrajID
state_c = child%StateID

if(l_select)then
    file_name = 'Select.log'
else
    file_name = 'FailSelect.log'
endif

call FMS_OpenFile(file_name, iunit, file_exists)

if(.not. file_exists) then
    write(iunit,'(A)') '#PredictedSTime ActualSTime'// &
                       '  CID  CSt  PID  PSt  PredictedS  ActualS'
endif

! write the selection log
101 format(2(1x,f9.2),4(1x,i4),2(1x,f9.5))
write(iunit, 101) locSelectionTime, currentTime, &
                     id_c, state_c, id_p, state_p,  &
                     predictedS, currentS
close(unit=iunit)

end subroutine FMS_WriteSelectionLog


subroutine fill_state_bundles(B1, BundleSS)
type(T_TrajectoryBundle), intent(in) :: B1
type(T_TrajectoryBundle), intent(inout) :: BundleSS(B1%NumStates)
integer (kind=DefInt) :: iState, iTraj, nStateTraj, iSSTraj, npart

npart = B1%Trajectory(1)%NumParticles

do iState =1, B1%NumStates

    ! Count the number of trajectories for a given state
    nStateTraj = 0
    do iTraj=1, B1%NumTraj
        if (B1%Trajectory(iTraj)%StateID == iState) then
            nStateTraj = nStateTraj + 1
        end if
    enddo

    call BundleSS(iState)%create(numtraj=nStateTraj, &
                                 numdeadtraj=0, &
                                 numstates=B1%NumStates, &
                                 numparticles=npart, &
                                 ncbfs=B1%NCBFs)

    ! Fill the the new bundle with trajs currently in iState
    iSSTraj = 1
    do iTraj = 1, B1%NumTraj
        if (B1%Trajectory(iTraj)%StateID == iState) then
            BundleSS(iState)%Trajectory(iSSTraj) = B1%Trajectory(iTraj)
            iSSTraj = iSSTraj + 1
        end if
    enddo
enddo

end subroutine fill_state_bundles

! Copy over trajectories after the selection from BundleSS
! to the original bundle. Note that we must copy over
! also the trajectories marked for death.
subroutine copy_state_bundles_to_original_bundle(B1, BundleSS)
type(T_TrajectoryBundle), intent(inout), target :: B1
type(T_TrajectoryBundle), intent(in), target :: BundleSS(B1%NumStates)
type(T_TrajectoryBundle), pointer :: BSS_i
type(T_Trajectory), pointer :: T_i
integer (kind=DefInt) :: iState, iTrj, iSSTraj

! Loop over all trajectories in the original bundle,
! find the matching trajectory in BundleSS
! and copy it over.
do iTrj = 1, B1%NumTraj

   T_i => B1%Trajectory(iTrj)

   do iState = 1, B1%NumStates

      if (T_i%StateID /= iState) cycle

      BSS_i => BundleSS(iState)

      do iSSTraj = 1, BSS_i%NumTraj
          if (T_i%TrajID == BSS_i%Trajectory(iSSTraj)%TrajID) then
              B1%Trajectory(iTrj) = BSS_i%Trajectory(iSSTraj)
              exit
          end if
      enddo

   enddo

enddo

end subroutine copy_state_bundles_to_original_bundle

end module SelectionModule
