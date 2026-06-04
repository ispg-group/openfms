! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module SelectionModule
   use GlobalModule
   use BundleModule
   use TrajectoryModule
   use BundleCalcsModule, only: FMS_bH, FMS_Norm, FMS_Branching
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
   real(kind=DefReal), parameter :: swissThresh = 3.678794411714d-1

contains

   subroutine print_stochastic_selection_params()

      if (glzStoSwiss) then

         write (fmiOut, *) 'Stochastic selection: Running in AIMSWISS mode!'

         write (fmiOut, *) 'SWISS mode sets DecoherenceTime to ', gldDecoherenceTime
         write (fmiOut, '(a,i0)') ' and IgnoreState to ', glIgnoreState
         write (fmiOut, *) 'so that Trajectories are only killed when they cannot spawn anymore.'

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
      type(T_TrajectoryBundle), allocatable :: BundleSS(:)
! sepSS: array of temporary Bundles for singlet and triplet separated stochastic selection
!        and the dimension of sep bundles array
      type(T_TrajectoryBundle), allocatable :: sepB(:)
! --- Added 4 sep SS
! sepSS: needed to determine total Coupled matrix before stochastic selection is performed
!        (same as in FMS_BuildCoupled and FMS_GroupIntoBlocks which are otherwise
!         only called from perform_stochastic_selection)
      integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj) :: Coupled
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1) :: blocktrajid
      integer(kind=DefInt), dimension(B1%NumTraj + 1) :: ntrajblock
      integer(kind=DefInt), dimension(B1%NumTraj + 1, 2) :: MergedBlockTrajID
      integer(kind=DefInt), dimension(B1%NCBFs + 1, 2) :: BlockCBFID
      integer(kind=DefInt), dimension(2) :: NumTrajSepBlock
      integer(kind=DefInt), dimension(2) :: NumCBFsBlock
      real(kind=DefReal), allocatable :: SaveNorms(:)
      real(kind=DefReal) :: RNorm
! sepSS: Shall we perform multiplet separated stochastic selection
      integer(kind=DefInt), allocatable :: StoSelMode(:)
      logical, dimension(2) :: IsSelection
      integer(kind=DefInt) :: iBlocks, j, NumBlockOuter
! --- Added 4 sep SS end ---
! AIMSWISS: Shall we perform stochastic selection
      logical :: performSelection
! AIMSWISS: Current selection time
      real(kind=DefReal) :: selectionTime
      integer(kind=DefInt) :: iState

      if (B1%NumTraj == 1) return

! Check that we are not just before a Spawning event
      if ((gldLastSpawnSto >= B1%CurrentTime) & ! When AIMSWISS is invoked don't
          .and. (.not. glzStoSwiss)) then ! wait for next Spawning event
         write (fmiout, *) 'Aborting Stochastic Collapse - Spawn is coming...'
         return
      end if

      if (glzStoSwiss) then
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!   Check if a selection time of any pair TBF is reached
         call FMS_CheckSelectionTime(B1, performSelection, selectionTime)
!   If not then there is nothing to do
         if (.not. performSelection) then
            return
         end if
      end if

! State-Specific Stochastic Selection (aka 4S :-)
! 1. Create temporary Bundles that hold a set of trajectories for a given state
!    We essentically split the Bundle into NumStates smaller bundles
! 2. Perform selections within those bundles
! 3. Convert the results into the original bundle.
      if (glzStoStateSpecific) then

         allocate (BundleSS(B1%NumStates))
         call fill_state_bundles(B1, BundleSS)

         do iState = 1, B1%NumStates
            if (BundleSS(iState)%NumTraj == 1) cycle
            call perform_stochastic_selection(BundleSS(iState), selectionTime)
         end do

         call copy_state_bundles_to_original_bundle(B1, BundleSS)

         do iState = 1, B1%NumStates
            call BundleSS(iState)%destroy()
         end do
         deallocate (BundleSS)

      end if

! Multiplicity-Specific Stochastic Selection
! 1. Decide whether we have Sing only, Trip only, or Sing and Trip
! 2. Create temporary, separate Bundles that only hold the trajectories within a multiplicity
! 3. Perform selections within those bundles
! 4. Reunite to reform the original bundle, get rid of temporary block Bundles

      if (NTrip /= 0) then

         ! (1) PREPARE B1 -----------------------------------------------------------------
         ! getting Coupled matrix to identify separated S only, or T only blocks

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "These trajs are in B1 (right at the start of FMS_StochasticCollapse)"
         !do i = 1, B1%NumTraj
         !   write (fmiout, *) "TrajID: ", B1%Trajectory(i)%TrajID, "CBF: ", B1%Trajectory(i)%CBF
         !   write (fmiout, *) "Ms: ", B1%Trajectory(i)%Ms
         !   write (fmiout, *) "Amp: ", B1%Trajectory(i)%Amplitude
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !write (fmiout, *) "-----------------------------------------------------------"
         !write (fmiout, *) "Check that Centroids of B1 remain unaffected"
         !write (fmiout, *) "Part 1: B1 Centroids before"
         !do iTraj = 1, (((B1%NCBFs - 1) * B1%NCBFs) / 2)
         !   write (fmiout, *) "Centroid number", iTraj
         !   write (fmiout, *) "Is centroid to trajectories", B1%Centroids(iTraj)%CentID
         !   write (fmiout, *) "And has position: ", B1%Centroids(iTraj)%Particle(1)%get_pos()
         !end do
         !write (fmiout, *) "-----------------------------------------------------------"

         ! - get the total converged coupled matrix, and the number of blocks
         !write (fmiout, *) "-----------------------------------------------------------"
         !write (fmiout, *) "Building coupled matrix for total Bundle (B1, unseparated)"

         ! Why setting gliForceKill to zero?
         ! Needed for if selection happens in Singlet and Triplet block in the same step
         ! (see line 420-480)
         gliForceKill = 0
         Coupled = 0
         call FMS_BuildCoupled(B1, Coupled, selectionTime)
         !write (fmiout, *) "This is the total, converged coupled matrix of size ", size(Coupled)
         !do iTraj = 1, B1%NumTraj
         !   write (fmiout, *) Coupled(iTraj,:)
         !end do

         !write (fmiout, *) "Get number of blocks for total Bundle (B1, unseparated)"
         blocktrajid = 0
         ntrajblock = 0
         NumBlockOuter = 0
         call FMS_GroupIntoBlocks(B1, Coupled, blocktrajid, ntrajblock, NumBlockOuter)
         !write (fmiout, *) "Number of Blocks: ", NumBlockOuter


         ! Get isSIng, isTrip, isSingTrip for each block in the original, total Coupled matrix
         allocate (StoSelMode(NumBlockOuter))
         StoSelMode = 0
         call FMS_GroupIntoMultiplicityBlocks(B1, StoSelMode, NumBlockOuter, blocktrajid)

         ! Merge blocks depending on isSing, isTrip
         ! (No action for isSingTrip because not to be touched in StoSel)
         !write (fmiout, *) "blocktrajid before merging: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) blocktrajid(i,:)
         !end do

         !! Set StoSelMode to zero if column of blocktrajid is fully zero
         !! i.e. if the block was copied over (the j column)
         !do i = 1, NumBlockOuter
         !   if (all(blocktrajid(:,i) == 0)) then
         !      StoSelMode(i) = 0
         !   end if
         !end do

         !write (fmiout, *) "blocktrajid after merging: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) blocktrajid(i,:)
         !end do

         MergedBlockTrajID = 0
         call FMS_getMergedBlockTrajID(B1, NumBlockOuter, StoSelMode, ntrajblock, blocktrajid, &
                                   MergedBlockTrajID)

         !write (fmiout, *) "MergedBlockTrajID after merging: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) MergedBlockTrajID(i, :)
         !end do

         !MergedBlockTrajID = blocktrajid
         !nblockMerged = NumBlockOuter        ! TODO: For now, nblockMerged is not used, but do we need to account for change
                                            ! in number of blocks after merging (so, nblockMerged < NumBlockOuter depending
                                            ! on how many blocks we merged???
                                            ! This is actually done two steps later I think

         ! Generate 2nd version of MergedBlockTrajID w/ indices ref to what will become the sep bundles
         ! This is the MergedSepBlockTrajID matrx (sorry complicated.. :( )
         ! In MergedSepBlockTrajID every column restarts from 1,
         ! while in MergedBlockTrajID indices continue
         ! ------------------------------------------------------------------------------------------------
         ! TODO: this doesn't seem robust at all, go back later and improve
         ! edit: seems a little more robust now, at least (still problematic: 
         !                                                 - is it fine to only adjust indices after the
         !                                                   first column = from the second block onwards??)

         !MergedSepBlockTrajID = MergedBlockTrajID
         !do i = 1, B1%NumTraj + 1
         !   do j = 2, B1%NumTraj + 1
         !      if (MergedBlockTrajID(i,j) == 0) cycle
         !      MergedSepBlockTrajID(i,j) = MergedBlockTrajID(i,j) - MergedBlockTrajID(1,j) + 1
         !   end do
         !end do

         ! Make sure that indices in MergedBlockTrajID leave out indices of trajs that died before
         ! This step is necessary because FMS_GroupIntoBlocks loops over B1%NumTraj
         !totDeadCount = 0
         !do i = 1, B1%NumTraj + 1
         !   do j = 2, B1%NumTraj + 1
         !      if (MergedBlockTrajID(i,j) == 0) cycle
         !      do k = 1, B1%NumDeadTraj
         !         if (MergedBlockTrajID(i,j) == B1%DeadTraj(k)%TrajID) then
         !            !write (fmiout, *) "This index belongs to a dead traj: ", MergedBlockTrajID(i,j), B1%DeadTraj(k)%TrajID
         !            !write (fmiout, *) "totDeadCount increased from ", totDeadCount
         !            totDeadCount = totDeadCount + 1
         !            !write (fmiout, *) "to... ", totDeadCount
         !         end if
         !      end do
         !      MergedBlockTrajID(i,j) = MergedBlockTrajID(i,j) + totDeadCount
         !   end do
         !end do


         ! ----------------------------- !!! --------------------------------------------------
         ! 21/04/26: Dead Traj IDs need to be incorporated to MergedBlockTrajID
         !           
         ! 22/04/26: Check with unsep SS: does BlockTrajID account for Dead Trajs?
         !           If so, why does it not in my case?
         ! ----------------------------- !!! --------------------------------------------------

         !write (fmiout, *) "This is the S, T, S/T StoSel mode merged blocktrajid: "
         !do j = 1, B1%NumTraj + 1
         !   write (fmiout, *) blocktrajid(j,:)
         !end do
         !write (fmiout, *) "And these are our StoSel (updated) modes: ", StoSelMode

         !write (fmiout, *) "This is the MergedBlockTrajID: "
         !do j = 1, B1%NumTraj + 1
         !   write (fmiout, *) MergedBlockTrajID(j,:)
         !end do

         !write (fmiout, *) "This is the MergedSepBlockTrajID with what will be sep bundle IDs: "
         !do j = 1, B1%NumTraj + 1
         !   write (fmiout, *) MergedSepBlockTrajID(j,:)
         !end do
         !write (fmiout, *) "And these are our StoSel (updated) modes: ", StoSelMode

         ! ----------------------------- !!! --------------------------------------------------

         ! Now we have everything for figuring out in which of the separated bundles selection
         ! needs to happen, and translating for back from the separated bundles to the total B1

         ! Determine how many block bundles we need to create, based on merged blocktrajid:
         ! and based on whether we need to carry out SS in the merged block or not
         !BundlesMultDim = 0
         !do i = 1, NumBlockOuter
         !   if (StoSelMode(i) /= 0) then
         !      BundlesMultDim = BundlesMultDim + 1
         !   else
         !      write (fmiout, *) "No separated Bundle will be created for block ", i, &
         !                        "because StoSelMode is zero"
         !   end if
         !end do

         !nblockMerged = BundlesMultDim
         !write (fmiout, *) "NumBlockOuter = ", NumBlockOuter, "StoSelMode = ", StoSelMode
         !write (fmiout, *) "BundlesMultDim based on merged blocktrajid: ", BundlesMultDim
         !write (fmiout, *) "nblock based on merged blocktrajid: ", nblockMerged

         ! Determine how many CBFs we have per S-only and T-only block
         ! sepB creation and subsequent selection only needs to happen if
         ! we have more than one CBF in the block

         NumTrajSepBlock = 0
         NumCBFsBlock = 0
         BlockCBFID = 0
         call FMS_getNumCBFsBlock(B1, MergedBlockTrajID, NumTrajSepBlock, NumCBFsBlock, BlockCBFID)

         !write (fmiout, *) "BlockCBFID after merging: "
         !do i = 1, B1%NCBFs + 1
         !   write (fmiout, *) BlockCBFID(i, :)
         !end do
         !write (fmiout, *) "NumCBFsBlock: ", NumCBFsBlock

         IsSelection = .false.
         do iBlocks = 1, 2
            if (NumCBFsBlock(iBlocks) > 1 ) then
               IsSelection(iBlocks) = .true. ! Set StoSelMode of block to zero if we have
                                             ! no more than one CBF in the block
            end if
         end do
         !write (fmiout, *) "Singlet selection set to ", IsSelection(1), "Triplet selection set to ", IsSelection(2)
         !write (fmiout, *) "Number Singlet CBFs", NumCBFsBlock(1), "Number Triplet CBFs", NumCBFsBlock(2)

         allocate (sepB(2))
         allocate (SaveNorms(2))

         !write (fmiout, *) "-----------------------------------------------------------"
         !write (fmiout, *) "The B1 norm is: "
         !write (fmiout, *) "Norm: ", FMS_Norm(B1)
         !write (fmiout, *) "and the corresponding amplitudes are"
         !do j = 1, B1%NumTraj
         !   write (fmiout, *) "Trajectory ", j, "TrajID ", B1%Trajectory(j)%TrajID, &
         !                     "AMplitude: ", B1%Trajectory(j)%Amplitude
         !end do

         do iBlocks = 1, 2
            if (IsSelection(iBlocks)) then
               !write (fmiout, *) "Creating separated Bundle for block ", iBlocks, &
               !                  "with StoSelMode ", StoSelMode(iBlocks)

               ! (2) CREATE SEPARATED BUNDLES ------------------------------------------------------
               ! Create and fill temp separated bundles
               ! with corresponding trajs from B1

               call FMS_getMultBundles(B1, sepB, iBlocks, MergedBlockTrajID(:, iBlocks), BlockCBFID(:,iBlocks), &
                                   NumTrajSepBlock(iBlocks), NumCBFsBlock(iBlocks))

               !write (fmiout, *) "After creating the sep bundle norms are: "
               !write (fmiout, *) "Block ", iBlocks, "Norm: ", FMS_Norm(sepB(iBlocks))
               !write (fmiout, *) "and the corresponding amplitudes are"
               !do j = 1, NumTrajSepBlock(iBlocks)
               !   write (fmiout, *) "Trajectory ", j, "TrajID ", sepB(iBlocks)%Trajectory(j)%TrajID, & 
               !                     "AMplitude: ", sepB(iBlocks)%Trajectory(j)%Amplitude
               !end do
               !write (fmiout, *) "----------------------------"

               ! Set the norm of each temporary bundle to one (for SS)
               ! ---> by just normalizing
               ! and remember actual norm within total bundle (for continuing propagation after SS)
               ! ---> SaveNorms
               !write (fmiout, *) "What's the norm of the newly created bundle?"
               !write (fmiout, *) "Block ", iBlocks, "Norm: ", FMS_Norm(sepB(iBlocks))
               SaveNorms(iBlocks) = FMS_Norm(sepB(iBlocks))

               !write (fmiout, *) "Did we copy over to SaveNorms correctly?", SaveNorms(iBlocks)

               !write (fmiout, *) "----------------------------"
               !write (fmiout, *) "Normalizing..."

               !write (fmiout, *) "SaveNorms(iBlocks)", SaveNorms(iBlocks)
               !write (fmiout, *) "NumTrajSepBlock(iBlocks)", NumTrajSepBlock(iBlocks)
               !write (fmiout, *) "MergedBlockTrajID(:,iBlocks)", MergedBlockTrajID(:,iBlocks)
               do j = 1, NumTrajSepBlock(iBlocks)
                  sepB(iBlocks)%Trajectory(j)%Amplitude = sepB(iBlocks)%Trajectory(j)%Amplitude &
                                                               * dcmplx( 1.d0 / sqrt(SaveNorms(iBlocks)) )
               end do

               !write (fmiout, *) "After Normalizing prior to SS, the sep bundle norms are: "
               !write (fmiout, *) "Block ", iBlocks, "Norm: ", FMS_Norm(sepB(iBlocks))
               !write (fmiout, *) "and the corresponding amplitudes are"
               !do j = 1, NumTrajSepBlock(iBlocks)
               !   write (fmiout, *) "Trajectory ", j, "TrajID ", sepB(iBlocks)%Trajectory(j)%TrajID, & 
               !                     "AMplitude: ", sepB(iBlocks)%Trajectory(j)%Amplitude
               !end do
               !write (fmiout, *) "----------------------------"

               ! (3) PERFORM SS ON SEPARATED BUNDLES ------------------------------------------------

               !write (fmiout, *) "Performing stochastic selection for block Bundle", iBlocks, &
               !                  "with", sepB(iBlocks)%NCBFs, "CBFs"

               call perform_stochastic_selection(sepB(iBlocks), selectionTime)

               ! Set the norm of each temporary bundle back to SaveNorms
               ! and renormalize adjusting traj amplitudes
               !write (fmiout, *) "----------------------------"
               !write (fmiout, *) "After Selection, the sep bundle norms are: "

               !write (fmiout, *) "Block ", iBlocks, "Norm: ", FMS_Norm(sepB(iBlocks))
               !write (fmiout, *) "and the corresponding amplitudes are"
               !do j = 1, NumTrajSepBlock(iBlocks)
               !   write (fmiout, *) "Trajectory ", j, "TrajID ", sepB(iBlocks)%Trajectory(j)%TrajID, & 
               !                     "AMplitude: ", sepB(iBlocks)%Trajectory(j)%Amplitude
               !end do
               !
               !write (fmiout, *) "----------------------------"
               !write (fmiout, *) "Renormalizing..."

               !write (fmiout, *) "SaveNorms(iBlocks)", SaveNorms(iBlocks)
               !write (fmiout, *) "NumTrajSepBlock(iBlocks)", NumTrajSepBlock(iBlocks)
               !write (fmiout, *) "MergedBlockTrajID(:,iBlocks)", MergedBlockTrajID(:,iBlocks)

               ! (4) WHO DIED AND CLEANING UP ---------------------------------------------------------

               ! Move indices of ForceKill trajs further in gliForceKill
               ! to avoid they are overwritten in triplet selection step.
               ! Important if singlet and triplet selection both happen in
               ! one selection step!
               if (iBlocks == 1) then
                  if (any(gliForceKill(:) /= 0)) then
                     do j = 1, NumTrajSepBlock(1)
                        gliForceKill( NumTrajSepBlock(2)+j ) = gliForceKill(j)
                        write (fmiout, *) "gliForceKill element ", NumTrajSepBlock(2)+j, &
                                          "will be overwritten to ", gliForceKill(j)
                     end do
                  end if
               end if

               RNorm = FMS_Norm(sepB(iBlocks))
               !write (fmiout, *) "RNorm", RNorm
               do j = 1, NumTrajSepBlock(iBlocks)
                  sepB(iBlocks)%Trajectory(j)%Amplitude = sepB(iBlocks)%Trajectory(j)%Amplitude &
                                                               * dcmplx( sqrt(SaveNorms(iBlocks) / RNorm ) )
               end do

               !write (fmiout, *) "After Renormalizing prior to copying back, the sep bundle norms are: "

               !write (fmiout, *) "Block ", iBlocks, "Norm: ", FMS_Norm(sepB(iBlocks))
               !write (fmiout, *) "and the corresponding amplitudes are"
               !do j = 1, NumTrajSepBlock(iBlocks)
               !   write (fmiout, *) "Trajectory ", j, "TrajID ", sepB(iBlocks)%Trajectory(j)%TrajID, & 
               !                     "AMplitude: ", sepB(iBlocks)%Trajectory(j)%Amplitude
               !end do
               !write (fmiout, *) "----------------------------"

               call FMS_copy_sepB_to_original_bundle(B1, sepB(iBlocks), iBlocks, NumTrajSepBlock(iBlocks), &
                                                        MergedBlockTrajID(:,iBlocks)) !, MergedSepBlockTrajID)
               !write (fmiout, *) "States were copied back without error"

               !write (fmiout, *) "(After copying back to original Bundle)"
               !write (fmiout, *) "Number of Dead Trajs is now: ", B1%NumDeadTraj


               ! Destroy the temporary block bundle
               call sepB(iBlocks)%destroy()

            !else if (iBlocks == 1) then
            !   write (fmiout, *) "No separated Bundle created for the singlet block &
            !                      because no more than one CBF."
            !else if (iBlocks == 2) then
            !   write (fmiout, *) "No separated Bundle created for the triplet block &
            !                      because no more than one CBF."
            end if
         end do

         deallocate(SaveNorms, sepB, StoSelMode)

         !write (fmiout, *) "After copy_MultBundles_to_original_Bundle, the B1 norm is: "
         !write (fmiout, *) "Norm: ", FMS_Norm(B1)
         !write (fmiout, *) "and the corresponding amplitudes are"
         !do j = 1, B1%NumTraj
         !   write (fmiout, *) "Trajectory ", j, "has ID", B1%Trajectory(j)%TrajID, "and AMplitude: ", &
         !                     B1%Trajectory(j)%Amplitude
         !end do

         !write (fmiout, *) "-----------------------------------------------------------"

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "Check that Centroids of B1 remain unaffected"
         !write (fmiout, *) "Part 2: B1 Centroids after"
         !do iTraj = 1, (((B1%NCBFs - 1) * B1%NCBFs) / 2)
         !   write (fmiout, *) "Centroid number", iTraj
         !   write (fmiout, *) "Is centroid to trajectories", B1%Centroids(iTraj)%CentID
         !   write (fmiout, *) "And has position: ", B1%Centroids(iTraj)%Particle(1)%get_pos()
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "Check that TrajIDs of B1 remain unaffected"
         !write (fmiout, *) "(At the end of StochasticCollapse)"
         !do iTraj = 1, B1%NumTraj
         !   write (fmiout, *) "Traj number", iTraj, "has TrajID", B1%Trajectory(iTraj)%TrajID
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "These trajs are in B1 (very end of FMS_StochasticCollapse)"
         !do i = 1, B1%NumTraj
         !   write (fmiout, *) "TrajID: ", B1%Trajectory(i)%TrajID, "CBF: ", B1%Trajectory(i)%CBF
         !   write (fmiout, *) "Ms: ", B1%Trajectory(i)%Ms
         !   write (fmiout, *) "Amp: ", B1%Trajectory(i)%Amplitude
         !end do
         !write (fmiout, *) " ----------------------------------------------------"


         ! Create BundlesMultDim many temp bundles and fill them up
         ! with corresponding trajs from B1
         !write (fmiout, *) "This is the blocktrajid to be used in getMultBundles: ", blocktrajid
         !call getMultBundles(B1, sepB, BundlesMultDim, MergedBlockTrajID)

         ! Set the norm of each temporary bundle to one (for SS)
         ! ---> by just normalizing
         ! and remember actual norm within total bundle (for continuing propagation after SS)
         ! ---> SaveNorms

         !allocate (SaveNorms(2))
         !!allocate (SaveTrajAmps(B1%NumTraj + 1, 2))
         !write (fmiout, *) "What's the norm of the newly created bundles?"
         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "Block ", i, "Norm: ", FMS_Norm(sepB(i))
         !      SaveNorms(i) = FMS_Norm(sepB(i))
         !      !call FMS_Branching(sepB(i), SaveTrajAmps(i,:))
         !   end if
         !end do
         !write (fmiout, *) "Did we copy over to SaveNorms correctly?", SaveNorms

         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      do iTraj = 1, NumTrajSepBlock(i)
         !         write(fmiout, *) "Amplitude of traj", sepB(i)%Trajectory(iTraj)%TrajID , &
         !                          sepB(i)%Trajectory(iTraj)%Amplitude
         !         !SaveTrajAmps(iTraj, i) = sepB(i)%Trajectory(iTraj)%Amplitude
         !      end do
         !   end if
         !end do

         !write (fmiout, *) "-----------------------------------------------------------"
         !write (fmiout, *) "The B1 norm is: "
         !write (fmiout, *) "Norm: ", FMS_Norm(B1)
         !write (fmiout, *) "and the corresponding amplitudes are"
         !do iTraj = 1, B1%NumTraj
         !   write (fmiout, *) "Trajectory ", iTraj, "TrajID ", B1%Trajectory(iTraj)%TrajID, &
         !                     "AMplitude: ", B1%Trajectory(iTraj)%Amplitude
         !end do

         !write (fmiout, *) "----------------------------"
         !write (fmiout, *) "Normalizing..."

         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "SaveNorms(i)", SaveNorms(i)
         !      write (fmiout, *) "NumTrajSepBlock(i)", NumTrajSepBlock(i)
         !      write (fmiout, *) "MergedBlockTrajID(:,i)", MergedBlockTrajID(:,i)
         !      !call FMS_RenormalizeBlockAmplitudes(sepB(i), SaveNorms(i), NumTrajSepBlock(i), MergedBlockTrajID(:,i))
         !      do j = 1, NumTrajSepBlock(i)
         !         !iTraj = MergedBlockTrajID(j,i)
         !         !sepB(i)%Trajectory(j)%Amplitude = sepB(i)%Trajectory(j)%Amplitude &
         !         !                                             * dcmplx( 1.d0 / sqrt(SaveNorms(i)) )
         !         sepB(i)%Trajectory(j)%Amplitude = sepB(i)%Trajectory(j)%Amplitude &
         !                                                      * dcmplx( 1.d0 / sqrt(SaveNorms(i)) )
         !      end do
         !   end if
         !end do
         
         !write (fmiout, *) "After Normalizing prior to SS, the sep bundle norms are: "
         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "Block ", i, "Norm: ", FMS_Norm(sepB(i))
         !      write (fmiout, *) "and the corresponding amplitudes are"
         !      do iTraj = 1, NumTrajSepBlock(i)
         !         write (fmiout, *) "Trajectory ", iTraj, "TrajID ", sepB(i)%Trajectory(iTraj)%TrajID, & 
         !                           "AMplitude: ", sepB(i)%Trajectory(iTraj)%Amplitude
         !      end do
         !   end if
         !end do
         !write (fmiout, *) "----------------------------"

         ! (3) Perform SS on the separated bundles
         !do i = 1, nblockMerged
         !do i = 1, 2

         !   if (IsSelection(i) .eqv. .true.) then

         !      !write (fmiout, *) "StoSelMode for block ", i, "is ", StoSelMode(i)
         !      ! Do the stochastic selection
         !      write (fmiout, *) "Performing stochastic selection for block Bundle", i, &
         !                        "with", sepB(i)%NCBFs, "CBFs"

         !      !call perform_stochastic_selection(sepB(i), selectionTime, IsSelection)
         !      call perform_stochastic_selection(sepB(i), selectionTime)

         !      ! Adjust gliForceKill only if selection actually happened
         !      if (IsSelection(i)) then

         !         write (fmiout, *) "Number of dead trajs: ", B1%NumDeadTraj
         !         do j = 1, B1%NumDeadTraj
         !            write (fmiout, *) "TrajID of dead traj", j, ": ", B1%DeadTraj(j)%TrajID
         !         end do

         !      end if

         !      ! Adjust gliForceKill only if selection actually happened
         !      !if (i>1 .and. IsSelection) then
         !      !   !write (fmiout, *) "Adjusting gliForceKill using MergedSepBlockTrajID (this is the full matrix): "
         !      !   !do j = 1, B1%NumTraj + 1
         !      !   !   write (fmiout, *) MergedSepBlockTrajID(j,:)
         !      !   !end do

         !      !   !write (fmiout, *) "and MergedBlockTrajID (this is the full matrix): "
         !      !   !do j = 1, B1%NumTraj + 1
         !      !   !   write (fmiout, *) MergedBlockTrajID(j,:)
         !      !   !end do

         !      !   write (fmiout, *) "gliForcKill for block Bundle", i
         !      !   write (fmiout, *) gliForceKill
         !      !   !write (fmiout, *) "adjusting to B1 indiced..."

         !      !   SaveForceKill = gliForceKill

         !      !   do iTraj = 1, B1%NumTraj + 1                       ! Why did I ever start this one from 2??
         !      !      do j = 1, B1%NumTraj + 1
         !      !         if (MergedSepBlockTrajID(iTraj,j) == 0) cycle
         !      !         do DeadID = 1, size(SaveForceKill)
         !      !            !write (fmiout, *) "many numbers (?)", SaveForceKill(DeadID)
         !      !            if (MergedSepBlockTrajID(iTraj,j) == SaveForceKill(DeadID)) then
         !      !               write (fmiout, *) "Match for gliForceKill: ", SaveForceKill(DeadID)
         !      !               write (fmiout, *) "with MergedSepblock: ", MergedSepBlockTrajID(iTraj,j)
         !      !               write (fmiout, *) "i and j are: ", iTraj, j
         !      !               write (fmiout, *) "Corresponding MergedBlock: ", MergedBlockTrajID(iTraj,j)
         !      !               gliForceKill(DeadID) = MergedBlockTrajID(iTraj,j)
         !      !            end if
         !      !         end do
         !      !      end do
         !      !   end do

         !      !   write (fmiout, *) "Adjusted gliForcKill for block Bundle", i, "and obtained"
         !      !   write (fmiout, *) gliForceKill
         !      !end if

         !      if (IsSelection(i)) then
         !         write (fmiout, *) "(Before copying back to original Bundle)"
         !         write (fmiout, *) "Number of Dead Trajs is now: ", B1%NumDeadTraj
         !      end if

         !      !if (B1%NumDeadTraj > 0) then

         !      !   do iTraj = 1, B1%NumDeadTraj
         !      !      write (fmiout, *) "IDs of dead traj ", iTraj, ": ", B1%DeadTraj(iTraj)%TrajID
         !      !   end do

         !      !   !if (B1%NumDeadTraj >= 2) then
         !      !   !   call FMS_DieError("That's enough for now...")
         !      !   !end if

         !      !end if

         !   !else
         !   !   if (any(blocktrajid(:,i) /= 0)) then
         !   !      write (fmiout, *) "No stochastic selection because block Bundle", i, "contains a mix of S and T"
         !   !   else
         !   !      write (fmiout, *) "No stochastic selection because block Bundle", i, "was merged into other block"
         !   !   end if

         !   end if

         !end do


         ! Set the norm of each temporary bundle back to SaveNorms
         ! and renormalize adjusting traj amplitudes
         ! ---> call FMS_Renormalize
         !write (fmiout, *) "----------------------------"
         !write (fmiout, *) "After Selection, the sep bundle norms are: "
         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "Block ", i, "Norm: ", FMS_Norm(sepB(i))
         !      write (fmiout, *) "and the corresponding amplitudes are"
         !      do iTraj = 1, NumTrajSepBlock(i)
         !         write (fmiout, *) "Trajectory ", iTraj, "TrajID ", sepB(i)%Trajectory(iTraj)%TrajID, & 
         !                           "AMplitude: ", sepB(i)%Trajectory(iTraj)%Amplitude
         !      end do
         !   end if
         !end do

         !write (fmiout, *) "Renormalizing..."
         !do i = 1,  BundlesMultDim
         !   call FMS_Renormalize(sepB(i), SaveTrajAmps(i,:), SaveNorms(i))
         !end do
         !write (fmiout, *) "----------------------------"
         !write (fmiout, *) "After Renormalizing, the sep bundle norms are: "
         !do i = 1, BundlesMultDim
         !   write (fmiout, *) "Block ", i, "Norm: ", FMS_Norm(sepB(i))
         !end do
         !write (fmiout, *) "and the corresponding amplitudes are"
         !do i = 1, BundlesMultDim
         !   write (fmiout, *) "Block: ", i
         !   do iTraj = 1, sepB(i)%NumTraj
         !      write (fmiout, *) "Trajectory ", iTraj, "TrajID ", sepB(i)%Trajectory(iTraj)%TrajID, & 
         !                        "AMplitude: ", sepB(i)%Trajectory(iTraj)%Amplitude
         !   end do
         !end do
         !write (fmiout, *) "----------------------------"

         ! (4) Set appropriate amplitudes to zero in the original bundle
         ! If selection happened, copy information which trajs died

         !if (any(gliForceKill(:) /= 0)) then ! TODO: use something else here
                                             ! otherwise copy_MultBundles_to_original_bundle always gets called as soon as 
                                             ! we have a number in gliForceKill
                                             ! How to determine if kill happened?
                                             ! ---
                                             ! For now, go with setting gliForceKill back to zero at the start of each
                                             ! StoSel (but not sure if that breaks things somewhere...?)


         !write (fmiout, *) "blocktrajid right before copy_MultBundles_to_B1: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) blocktrajid(i,:)
         !end do

         !write (fmiout, *) "MergedBlockTrajID: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) MergedBlockTrajID(i,:)
         !end do

         !write (fmiout, *) "MergedSepBlockTrajID: "
         !do i = 1, B1%NumTraj + 1
         !   write (fmiout, *) MergedBlockTrajID(i,:)
         !end do

         !write (fmiout, *) "----------------------------"
         !write (fmiout, *) "Renormalizing..."

         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "SaveNorms(i)", SaveNorms(i)
         !      write (fmiout, *) "NumTrajSepBlock(i)", NumTrajSepBlock(i)
         !      write (fmiout, *) "MergedBlockTrajID(:,i)", MergedBlockTrajID(:,i)
         !      !call FMS_RenormalizeBlockAmplitudes(sepB(i), SaveNorms(i), NumTrajSepBlock(i), MergedBlockTrajID(:,i))
         !      RNorm = FMS_Norm(sepB(i))
         !      write (fmiout, *) "RNorm", RNorm
         !      do j = 1, NumTrajSepBlock(i)
         !         !iTraj = MergedBlockTrajID(j,i)
         !         !sepB(i)%Trajectory(j)%Amplitude = sepB(i)%Trajectory(j)%Amplitude &
         !         !                                             * dcmplx( sqrt(SaveNorms(i)) )
         !         sepB(i)%Trajectory(j)%Amplitude = sepB(i)%Trajectory(j)%Amplitude &
         !                                                      * dcmplx( sqrt(SaveNorms(i) / RNorm ) )
         !      end do
         !   end if
         !end do
         
         !write (fmiout, *) "After Renormalizing prior to copying back, the sep bundle norms are: "
         !do i = 1, 2
         !   if (IsSelection(i)) then
         !      write (fmiout, *) "Block ", i, "Norm: ", FMS_Norm(sepB(i))
         !      write (fmiout, *) "and the corresponding amplitudes are"
         !      do iTraj = 1, NumTrajSepBlock(i)
         !         write (fmiout, *) "Trajectory ", iTraj, "TrajID ", sepB(i)%Trajectory(iTraj)%TrajID, & 
         !                           "AMplitude: ", sepB(i)%Trajectory(iTraj)%Amplitude
         !      end do
         !   end if
         !end do
         !write (fmiout, *) "----------------------------"

         !do iBlocks = 1, 2
         !   if (IsSelection(iBlocks)) then
         !      call copy_MultBundles_to_original_bundle(B1, sepB(iBlocks), iBlocks, NumTrajSepBlock(iBlocks), &
         !                                               MergedBlockTrajID(:,iBlocks)) !, MergedSepBlockTrajID)
         !      ! 18/05/26: Don't we need to normalise here ???!!?!?!
         !      write (fmiout, *) "States were copied back without error"
         !   end if
         !end do

         !do iBlocks = 1, 2
         !   if (IsSelection(iBlocks)) then
         !      write (fmiout, *) "(After copying back to original Bundle)"
         !      write (fmiout, *) "Number of Dead Trajs is now: ", B1%NumDeadTraj

         !      write (fmiout, *) "After copy_MultBundles_to_original_Bundle, the B1 norm is: "
         !      write (fmiout, *) "Norm: ", FMS_Norm(B1)
         !      write (fmiout, *) "and the corresponding amplitudes are"
         !      do iTraj = 1, B1%NumTraj
         !         write (fmiout, *) "Trajectory ", iTraj, "has ID", B1%Trajectory(iTraj)%TrajID, "and AMplitude: ", &
         !                           B1%Trajectory(iTraj)%Amplitude
         !      end do
         !   end if
         !end do

         !write (fmiout, *) "-----------------------------------------------------------"
         
         ! Destroy the temporary block bundles
         !do i = 1, nblockMerged
         !do iBlocks = 1, 2
         !   if (IsSelection(iBlocks)) then
         !      call sepB(iBlocks)%destroy()
         !   end if
         !end do
         !deallocate(SaveNorms, sepB, StoSelMode)

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "Check that Centroids of B1 remain unaffected"
         !write (fmiout, *) "Part 2: B1 Centroids after"
         !do iTraj = 1, (((B1%NCBFs - 1) * B1%NCBFs) / 2)
         !   write (fmiout, *) "Centroid number", iTraj
         !   write (fmiout, *) "Is centroid to trajectories", B1%Centroids(iTraj)%CentID
         !   write (fmiout, *) "And has position: ", B1%Centroids(iTraj)%Particle(1)%get_pos()
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "Check that TrajIDs of B1 remain unaffected"
         !write (fmiout, *) "(At the end of StochasticCollapse)"
         !do iTraj = 1, B1%NumTraj
         !   write (fmiout, *) "Traj number", iTraj, "has TrajID", B1%Trajectory(iTraj)%TrajID
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !write (fmiout, *) " ----------------------------------------------------"
         !write (fmiout, *) "These trajs are in B1 (very end of FMS_StochasticCollapse)"
         !do i = 1, B1%NumTraj
         !   write (fmiout, *) "TrajID: ", B1%Trajectory(i)%TrajID, "CBF: ", B1%Trajectory(i)%CBF
         !   write (fmiout, *) "Ms: ", B1%Trajectory(i)%Ms
         !   write (fmiout, *) "Amp: ", B1%Trajectory(i)%Amplitude
         !end do
         !write (fmiout, *) " ----------------------------------------------------"

         !call FMS_DieError("Now, please go back and reconsider your life choices!")

         !else

         !   call perform_stochastic_selection(B1, selectionTime)

         !end if

         !call FMS_DieError("check that this works...")

      else

         call perform_stochastic_selection(B1, selectionTime)

      end if

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
      integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj) :: Coupled
! Number of blocks
      integer(kind=DefInt) :: nblock
! Matrix containing IDs of TBFs belonging to a block.
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1) :: blocktrajid
! Array containing number of TBFs belonging to a block.
      integer(kind=DefInt), dimension(B1%NumTraj + 1) :: ntrajblock
! Population of all blocks combined (not necessarily 1)
      real(kind=DefReal) :: totBlockpop
! Array containing the populations of each block
      real(kind=DefReal), dimension(B1%NumTraj + 1) :: Blockpop
! Index of selected block
      integer(kind=DefInt) :: iblockslct
! AIMSWISS: Current selection time
      real(kind=DefReal), intent(in) :: selectionTime
!! sep SS: Did selection happen?
!      logical, optional, intent(inout) :: IsSelection

      !write (fmiout, *) "Number of traj used for dimension of Coupled: ", B1%NumTraj
      !write (fmiout, *) "meaning Coupled dimension is: ", "( ", B1%NumTraj, ", ", B1%NumTraj, " )"

      !write (fmiout, *) "CentIDs in perform_stochastic_selection (iTraj, iCBF:"
      !do iTraj = 1, B1%NumTraj
      !   write (fmiout, *) B1%Trajectory(iTraj)%TrajID, B1%Trajectory(iTraj)%CBF
      !end do
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! First, decompose the hamiltonian into a block diagonal representation
! work out the coupling matrix
      Coupled = 0
      call FMS_BuildCoupled(B1, Coupled, selectionTime)

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Next, group trajectories into coupled blocks
      blocktrajid = 0
      ntrajblock = 0
      nblock = 0
      call FMS_GroupIntoBlocks(B1, Coupled, blocktrajid, ntrajblock, nblock)
      !write (fmiout, *) "Right before performing SS:"
      !write (fmiout, *) "Number of Blocks: ", nblock
      !write (fmiout, *) "Coupled is: ", Coupled

! There should be at least one block of trajectories
      if (nblock == 0) then
         write (fmiOut, *) 'Number of trajectory basis blocks is zero'
         call FMS_DieError("ERROR in FMS_StochasticCollapse")
      end if
! If just a single block then no more work to be done
      !if (nblock == 1) write (fmiout, *) "Only one block. No selection!"
      if (nblock == 1) return
! At this point, more than one uncoupled block of trajectories,
! so we want to stochastically collapse to one of them
!      if (present(IsSelection)) then
!         IsSelection = .true.
!      end if

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Compute populations of each block
      totBlockpop = 0.0d0
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
                                          ntrajblock(iblockslct), &
                                          blocktrajid(:, iblockslct))

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! Remove all TBFs that do not belong selected block
      call FMS_KillOtherBlocks(B1, nblock, ntrajblock, iblockslct, blocktrajid)
      !write (fmiout, *) "Leaving perform_stochastic_selection"

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
      child_i%SWISS%ParentOverlap = abs(overlap(parent_s, child_s))**2

      write (fmiOut, '(a,i0,a,f0.2)') 'SWISS: Trajectory ', parent_s%TrajID, ' and '// &
         'its child will decohere at t = ', child_i%SWISS%SelectionTime
      write (fmiOut, '(a,f5.3)') 'Their current absolute overlap is ', child_i%SWISS%ParentOverlap

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
      type(T_Trajectory), intent(in) :: parent, child
      integer(kind=DefInt), intent(in) :: npart
      real(kind=DefReal) :: decoherenceTime
      integer(kind=DefInt) :: ipart, jdim, ishifted, istart, ndim, &
                              parentSt, childSt
      real(kind=DefReal) :: forceElementDiff2, decoherenceRate2, &
                            parentForceElement, childForceElement, &
                            summand2, width

! Calculation of the squared decoherence rate
      parentSt = parent%StateID
      childSt = child%StateID
! stochastic selection is currently only implemented for
! multi-state spawning (may change in the future)
      if (childSt == parentSt) then
         write (fmiOut, *) "SWISS: States don't differ, decoherence not possible"
         call FMS_DieError("ERROR in FMS_CalculateDecoherenceTime")
      end if

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
                                 childForceElement)**2
            summand2 = forceElementDiff2 / (4.d0 * width)
            decoherenceRate2 = decoherenceRate2 + summand2
         end do
      end do

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
         if (B1%Trajectory(i)%SWISS%SelectionTime <= B1%CurrentTime) then
            performSelection = .true.
            selectionTime = B1%Trajectory(i)%SWISS%SelectionTime
            write (fmiout, *) 'SWISS: Time to select! The current time is ', &
               B1%CurrentTime, ' and the selection time is ', &
               selectionTime
            exit
         end if
      end do

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

      !write (fmiout, *) "---------------------------------------------------------------------------"
      !write (fmiout, *) " Inside FMS_BuildCoupled "
      !write (fmiout, *) "The number of traj we are using is: ", B1%NumTraj

      ntraj = B1%NumTraj

      do i = 1, ntraj
         Coupled(i, i) = 1
         do j = i, ntraj
            !  X.Z. Keep different trajectories in the same CBF together.
            if (B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf) then
               Coupled(i, j) = 1
               Coupled(j, i) = 1
               !write (fmiout, *) "Coupled(i,j) set to 1 because trajs have the same CBF"
               !write (fmiout, *) "CBFi, CBFj", B1%Trajectory(i)%cbf, B1%Trajectory(j)%cbf
               cycle
            end if
         end do
      end do

      if (glzStoSwiss) then

         call FMS_BuildCoupled_SWISS(B1, Coupled, selectionTime)

      else if (glzStoOlap) then

         call FMS_BuildCoupled_OSS(B1, Coupled)

      else

         !do i = 1, ntraj
         !   write (fmiout, *) "Building coupled for a triplet?", B1%Trajectory(i)%triplet
         !end do

         call FMS_BuildCoupled_ESS(B1, Coupled)
         !write (fmiOut, *) "Coupled matrix is:", Coupled, ",   Size is:", size(Coupled)

      end if

      !write (fmiout, *) "---------------------------------------------------------------------------"

   end subroutine FMS_BuildCoupled

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_BuildCoupled_ESS(B1, Coupled)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
         intent(inout) :: Coupled
      integer(kind=DefInt) :: ntraj, i, j

      ntraj = B1%NumTraj

      write (fmiout, *) "FMS_BuildCoupled_ESS is getting called :S"
      !write (fmiout, *) "Calculating coupling for trajs in bundle"
      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Are traj indices still fine here?"
      !write (fmiout, *) "Number trajs:", B1%NumTraj
      !write (fmiout, *) "Number CBFs:", B1%NCBFs
      !write (fmiout, *) " ----------------------------------------------------"

      do i = 2, ntraj
         do j = 1, i
            !write (fmiout, *) "ID, CBF, is triplet?"
            !write (fmiout, *) "i", B1%Trajectory(i)%TrajID, B1%Trajectory(i)%CBF, B1%Trajectory(i)%triplet
            !write (fmiout, *) "j", B1%Trajectory(j)%TrajID, B1%Trajectory(j)%CBF, B1%Trajectory(j)%triplet
            if (B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf) cycle
            if (abs(FMS_bH(B1, i, j)) > gldStochaThresh) then
               Coupled(i, j) = 1
               Coupled(j, i) = 1
            end if
         end do
      end do

! Uncover all indirect connections between TBFs
!!!!!!!!!!!!!!! Look at thattt:
      call FMS_ConvergeCoupled(ntraj, Coupled)

   end subroutine FMS_BuildCoupled_ESS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_BuildCoupled_OSS(B1, Coupled)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
         intent(inout) :: Coupled

! Overlap between TBF i and j
      complex(kind=DefComp) :: S_ij
      integer(kind=DefInt) :: ntraj, i, j

      !write (fmiout, *) "FMS_BuildCoupled_OSS is getting called"

      ntraj = B1%NumTraj

      do i = 2, ntraj
         do j = 1, i
            if (B1%Trajectory(i)%cbf == B1%Trajectory(j)%cbf) cycle
            S_ij = overlap(B1%Trajectory(i), &
                           B1%Trajectory(j))
            if (abs(S_ij) > gldStochaThresh) then
               Coupled(i, j) = 1
               Coupled(j, i) = 1
            end if
         end do
      end do

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
      logical :: success

      ntraj = B1%NumTraj
      success = .true.
      do i = 2, ntraj
         do j = 1, i
            if (B1%Trajectory(i)%ParentID == B1%Trajectory(j)%TrajID) then
               ! Add connection for all parent-child TBF pairs that haven't
               ! decohered!
               if (.not. (abs(B1%Trajectory(i)%SWISS%SelectionTime &
                              - selectionTime) <= FPZero)) then
                  Coupled(i, j) = 1
                  Coupled(j, i) = 1
               else
                  S_ik = abs(overlap(B1%Trajectory(i), B1%Trajectory(j)))**2
                  S_ik_max = swissThresh * B1%Trajectory(i)%SWISS%ParentOverlap
                  write (fmiout, '(a,i0,a,i0,a)') 'Check if overlap of the pair ', &
                     B1%Trajectory(i)%TrajID, ' and ', B1%Trajectory(j)%TrajID, &
                     ' is not too large'
                  if (S_ik > S_ik_max) then
                     success = .false.
                  end if
                  call FMS_WriteSelectionLog(B1%Trajectory(j), &
                                             B1%Trajectory(i), &
                                             S_ik, S_ik_max, &
                                             success)
                  success = .true.
               end if
            end if
         end do
      end do

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
      integer(kind=DefInt), dimension(ntraj, ntraj), &
         intent(inout) :: Coupled
! Coupling matrix for TBF basis of the previous iteration
      integer(kind=DefInt), dimension(ntraj, ntraj) :: Coupled_prev
      do
         Coupled_prev = Coupled

         Coupled = min(matmul(Coupled, Coupled), 1)

         if (all(Coupled == Coupled_prev)) exit
      end do

   end subroutine FMS_ConvergeCoupled

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_GroupIntoBlocks(B1, Coupled, blocktrajid, &
                                  ntrajblock, nblock)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), dimension(B1%NumTraj, B1%NumTraj), &
         intent(in) :: Coupled
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1), &
         intent(inout) :: blocktrajid
      integer(kind=DefInt), dimension(B1%NumTraj + 1), &
         intent(inout) :: ntrajblock
      integer(kind=DefInt), intent(inout) :: nblock

! Boolean array containing IDs of already considered TBFs
      logical, dimension(B1%NumTraj + 1) :: ztrajdone

      integer(kind=DefInt) :: ntraj, i, j, l, itraj

      !write (fmiout, *) "========================================="
      !write (fmiout, *) "What happens in GroupIntoBlocks?"

      ntraj = B1%NumTraj

      ztrajdone(:) = .false.

! Loop over ntraj possible blocks
      do i = 1, ntraj + 1

!   Add the first unsorted trajectory to this current block (for first block,
!   will always be first trajectory)
         do j = 1, ntraj
            if (ztrajdone(j)) cycle
            blocktrajid(1, i) = j
            ztrajdone(j) = .true.
            !write (fmiout, *) "The hitherto unsorted traj ", j, "was added to the block as the first traj"
            exit
         end do

         itraj = blocktrajid(1, i) !first trajectory in block

         !write (fmiout, *) "This is blocktrajid inside the i loop, i = ", i
         !do j = 1, B1%NumTraj + 1
         !   write (fmiout, *) blocktrajid(j,:)
         !end do

!   If no trajectories added, then exit loop
         if (itraj == 0) then
            nblock = i - 1
            !write (fmiout, *) "nblock is now ", nblock
            exit
         end if

!   Add other coupled trajectories to this current block
         l = 1
         do j = 1, ntraj
            !write (fmiout, *) "Now sorting traj j = ", j
            if (ztrajdone(j)) cycle !if this trajectory has already been sorted skip
            if (Coupled(j, itraj) == 1) then
               !write (fmiout, *) "traj ", j, "is in the same block as traj ", itraj
               l = l + 1
               blocktrajid(l, i) = j
               ztrajdone(j) = .true.
            end if
         end do
         ntrajblock(i) = l

      end do

      !write (fmiout, *) "This is the finished blocktrajid: "
      !do j = 1, B1%NumTraj + 1
      !   write (fmiout, *) blocktrajid(j,:)
      !end do

      !write (fmiout, *) "These trajectories are in the first block: ", blocktrajid(:,1)
      !write (fmiout, *) "and these are in the second block: ", blocktrajid(:,2)

      !write (fmiout, *) "========================================="

   end subroutine FMS_GroupIntoBlocks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_ComputeBlockPopulations(B1, blocktrajid, ntrajblock, &
                                          nblock, totBlockpop, Blockpop)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1), &
         intent(in) :: blocktrajid
      integer(kind=DefInt), dimension(B1%NumTraj + 1), &
         intent(in) :: ntrajblock
      integer(kind=DefInt), intent(in) :: nblock
      real(kind=DefReal), intent(inout) :: totBlockpop
      real(kind=DefReal), dimension(B1%NumTraj + 1), &
         intent(inout) :: Blockpop
! Temporaty FMS Bundle needed to store blocks
      type(T_TrajectoryBundle) :: BTemp

      integer(DefInt) :: i, j, jtraj, nstate, npart, ncbf

      nstate = B1%NumStates
      npart = B1%Trajectory(1)%NumParticles
      ncbf = B1%NCBFs

! Compute populations of each block:
!    Do this by building a temporary bundle of each block then use existing subroutine norm(bundle)
      do i = 1, nblock
         call BTemp%create(numtraj=ntrajblock(i), &
                           numdeadtraj=0, &
                           numstates=nstate, &
                           numparticles=npart, &
                           ncbfs=ncbf)

         do j = 1, ntrajblock(i)
            jtraj = blocktrajid(j, i)
            BTemp%Trajectory(j) = B1%Trajectory(jtraj)
         end do

         Blockpop(i) = FMS_Norm(BTemp)

         call BTemp%destroy()

         totBlockpop = totBlockpop + Blockpop(i)
      end do

   end subroutine FMS_ComputeBlockPopulations

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_SelectBlock(ntraj, nblock, ntrajblock, totBlockpop, &
                              Blockpop, iblockslct)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      integer(kind=DefInt), intent(in) :: ntraj
      integer(kind=DefInt), intent(in) :: nblock
      real(kind=DefReal), intent(in) :: totBlockpop
      real(kind=DefReal), dimension(ntraj + 1), &
         intent(in) :: Blockpop
      integer(kind=DefInt), dimension(ntraj + 1), &
         intent(in) :: ntrajblock
      integer(kind=DefInt), intent(inout) :: iblockslct
! Character denoting selected block
      character(len=1) :: cselect
! Random number used for Monte Carlo procedure
      real(kind=DefReal) :: drand
! Dummy variable for storing relative population of block
      real(kind=DefReal) :: Blockpopnorm
! Cumulative population of blocks used in Monte Carlo procedure
      real(kind=DefReal) :: Blockpopnormsum

      integer(kind=DefInt) :: i

!     Use random number to select which block to collapse to
      drand = fms_ranb(i4zero)
!       find which block was randomly selected ---------------- eq (14) SSAIMS paper
      Blockpopnormsum = 0.0d0
      do i = 1, nblock
         Blockpopnorm = Blockpop(i) / totBlockpop
         Blockpopnormsum = Blockpopnormsum + Blockpopnorm
         if (drand <= Blockpopnormsum) then
            iblockslct = i
            exit
         end if
      end do

      write (fmiOut, *) 'Stochastically collapsing'
      write (fmiOut, *) 'Block   NTraj   Pop     Selected'
      do i = 1, nblock
         cselect = ''
         if (i == iblockslct) cselect = '*'
         write (fmiOut, '(2(I3,5X),f8.4,X,A1)') i, ntrajblock(i), Blockpop(i), &
            cselect
      end do

   end subroutine FMS_SelectBlock

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_RenormalizeBlockAmplitudes(B1, slctBlockpop, slctNtrajblock, &
                                             slctBlocktrajid)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      type(T_TrajectoryBundle), intent(inout) :: B1

      real(kind=DefReal), intent(in) :: slctBlockpop
      integer(kind=DefInt), intent(in) :: slctNtrajblock
      integer(kind=DefInt), dimension(B1%NumTraj + 1), &
         intent(in) :: slctBlocktrajid
! Norm of molecular wavefunction before selection
      real(kind=DefReal) :: NormInit

! Renormalization factor
      complex(kind=DefComp) :: cnorm

      integer(kind=DefInt) :: j, jtraj

!     Renormalize selected block to total initial Norm (need not be 1)
      NormInit = FMS_Norm(B1)
      !write (fmiout, *) "Inside FMS_RenormalizeBlockAmplitude: "
      !write (fmiout, *) "This is NormInit", NormInit
      !write (fmiout, *) "This is slctBlockpop", slctBlockpop
      !write (fmiout, *) "This is slctNtrajblock", slctNtrajblock
      !write (fmiout, *) "This is slctBlockTrajID", slctBlocktrajid
      cnorm = dcmplx(sqrt(NormInit / slctBlockpop))
      !write (fmiout, *) "cnorm = ", cnorm
      !write (fmiout, *) "Renormalising the amplitudes (Inside perform_stochastic_selection): "
      do j = 1, slctNtrajblock
         jtraj = slctBlocktrajid(j)
         !write (fmiout, *) "Before (jtraj, Amp): ", jtraj, B1%Trajectory(jtraj)%Amplitude
         B1%Trajectory(jtraj)%Amplitude = B1%Trajectory(jtraj)%Amplitude * cnorm
         !write (fmiout, *) "After (jtraj, Amp): ", jtraj, B1%Trajectory(jtraj)%Amplitude
      end do

   end subroutine FMS_RenormalizeBlockAmplitudes

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine FMS_KillOtherBlocks(B1, nblock, ntrajblock, iblockslct, blocktrajid)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      type(T_TrajectoryBundle), intent(inout) :: B1
      integer(kind=DefInt), intent(in) :: nblock
      integer(kind=DefInt), dimension(B1%NumTraj + 1), &
         intent(in) :: ntrajblock
      integer(kind=DefInt), intent(in) :: iblockslct
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1), &
         intent(in) :: blocktrajid

      integer(kind=DefInt) :: i, j, l, jtraj
! Mark all other trajectories for death
      l = 0
      !write (fmiout, *) "-----------------------------------------------"
      !write (fmiout, *) "Marking trajectories for death"
      do i = 1, nblock
         !write (fmiout, *) "Current block: ", i
         !write (fmiout, *) "was it selected? ", i == iblockslct
         if (i == iblockslct) cycle
         do j = 1, ntrajblock(i)
            !write (fmiout, *) "We have a j! ", j
            l = l + 1
            jtraj = blocktrajid(j, i)
            gliForceKill(l) = B1%Trajectory(jtraj)%TrajID
            !write (fmiout, *) "Trajectory will be killed (TrajID, CBF): ", B1%Trajectory(jtraj)%TrajID, B1%Trajectory(jtraj)%CBF
            B1%Trajectory(jtraj)%Amplitude = dcmplx(0.0d0, 0.0d0)
         end do
      end do
      !write (fmiout, *) "Array of TrajIDs for dead/ non dead Trajs: ", gliForceKill, "size is: ", size(gliForceKill)
      !write (fmiout, *) "-----------------------------------------------"

! Remove any dead trajectory amplitudes
      do jtraj = 1, B1%NumDeadTraj
         B1%DeadTraj(jtraj)%Amplitude = dcmplx(0.0d0, 0.0d0)
         !write (fmiout, *) "Killing this trajectory: ", jtraj
         !write (fmiout, *) "Out of this many total dead trajs: ", B1%NumDeadTraj
         !write (fmiout, *) "Was it a triplet?", B1%DeadTraj(jtraj)%triplet
      end do

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
      integer(kind=DefInt) :: id_c, state_c, &
                              id_p, state_p
      integer(DefInt) :: iunit
      logical :: file_exists

      currentTime = parent%get_time()
      locSelectionTime = child%SWISS%SelectionTime

      id_p = parent%TrajID
      state_p = parent%StateID
      id_c = child%TrajID
      state_c = child%StateID

      if (l_select) then
         file_name = 'Select.log'
      else
         file_name = 'FailSelect.log'
      end if

      call FMS_OpenFile(file_name, iunit, file_exists)

      if (.not. file_exists) then
         write (iunit, '(A)') '#PredictedSTime ActualSTime'// &
            '  CID  CSt  PID  PSt  PredictedS  ActualS'
      end if

! write the selection log
101   format(2(1x, f9.2), 4(1x, i4), 2(1x, f9.5))
      write (iunit, 101) locSelectionTime, currentTime, &
         id_c, state_c, id_p, state_p, &
         predictedS, currentS
      close (unit=iunit)

   end subroutine FMS_WriteSelectionLog

   subroutine FMS_GroupIntoMultiplicityBlocks(B1, StoSelMode, NumBlockOuter, blocktrajid)

      ! - merge the total coupled matrix according to
      !   - Singlet only --> perform StoSel
      !   - Triplet only --> perform StoSel
      !   - Singlet and Triplet in the same block --> no StoSel on these mixed blocks

      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1), intent(in) :: blocktrajid
      integer(kind=DefInt), intent(in) :: NumBlockOuter
      integer(kind=DefInt), intent(inout), dimension(NumBlockOuter) :: StoSelMode(:)
      integer(kind=DefInt) :: iBlocks, iTraj
      logical :: isSing, isTrip, isSingTrip

      do iBlocks = 1, NumBlockOuter
         isSing = .false.
         isTrip = .false.
         isSingTrip = .false.
         do iTraj = 1, B1%NumTraj
            if (any(blocktrajid(:,iBlocks) == iTraj)) then
               if (B1%Trajectory(iTraj)%triplet) then
                  isTrip = .true.
               else
                  isSing = .true.
               end if

            end if
         end do
         if (isTrip .eqv. .true. .and. isSing .eqv. .true.) then
            isSingTrip = .true.
            StoSelMode(iBlocks) = 0 ! no selection because Sing and Trip are coupled
         else if (isSing .eqv. .true.) then
            StoSelMode(iBlocks) = 1 ! selection for Sing
         else if (isTrip .eqv. .true.) then
            StoSelMode(iBlocks) = 3 ! selection for Trip
         else
            call FMS_DieError("No valid Selection mode could be determined")
         end if

         !write (fmiout, *) "For block ", i, "StoSel mode is: "
         !write (fmiout, *) "isSing is ", isSing, ", isTrip is ", isTrip, ", and... isSingTrip is", isSingTrip
         !write (fmiout, *) "This is for the trajs ", blocktrajid(:,i)
         !write (fmiout, *) "This is the current StoSel mode: ", StoSelMode(i)
      end do

   end subroutine FMS_GroupIntoMultiplicityBlocks

   subroutine FMS_getMergedBlockTrajID(B1, NumBlockOuter, StoSelMode, ntrajblock, blocktrajid, &
                                   MergedBlockTrajID)
      
      type(T_TrajectoryBundle), intent(in) :: B1
      integer(kind=DefInt), intent(in) :: NumBlockOuter
      integer(kind=DefInt), dimension(NumBlockOuter), intent(in) :: StoSelMode(:)
      integer(kind=DefInt), dimension(B1%NumTraj + 1), intent(in) :: ntrajblock
      integer(kind=DefInt), dimension(B1%NumTraj + 1, B1%NumTraj + 1), intent(in) :: blocktrajid
      integer(kind=DefInt), dimension(B1%NumTraj + 1, 2), intent(inout) :: MergedBlockTrajID
      integer(kind=DefInt) :: i, j, currStoSelMode, skip_positions

      do i = 1, NumBlockOuter

         !write (fmiout, *) "Entering loop for i = ", i
         !write (fmiout, *) "Part of blocktrajid we are using: ", blocktrajid(i,:)
         !write (fmiout, *) "General question: what is NumTrajBlock???", NumTrajBlock

         if (StoSelMode(i) == 0) cycle

         currStoSelMode = StoSelMode(i)

         !write (fmiout, *) "In outer block ", i, "with StoSelMode ", currStoSelMode

         if (currStoSelMode == 1) then

            skip_positions = 1
            do while ( MergedBlockTrajID(skip_positions, 1) /= 0 )
               !write (fmiout, *) "skip_positions: ", skip_positions, "ID", MergedBlockTrajID(skip_positions, 1)
               skip_positions = skip_positions + 1
            end do
            skip_positions = skip_positions - 1

            do j = 1, ntrajblock(i)
               MergedBlockTrajID(j+skip_positions, 1) = blocktrajid(j,i)
            end do

         else if (currStoSelMode == 3) then

            skip_positions = 1
            do while ( MergedBlockTrajID(skip_positions, 2) /= 0 )
               skip_positions = skip_positions + 1
            end do
            skip_positions = skip_positions - 1

            do j = 1, ntrajblock(i)
               MergedBlockTrajID(j+skip_positions, 2) = blocktrajid(j,i)
            end do

         end if

         !do j = i+1, NumBlockOuter
         !   write (fmiout, *) "Entering loop for j = ", j

         !   if (StoSelMode(j) == currStoSelMode) then

         !      ! Copy over block j TrajIDs to i column of blocktrajid
         !      ! Set the NumTrajBlock(i:j) positions in the i column to the indices in the jth column
         !      !blocktrajid(NumTrajBlock(i)+1:NumTrajBlock(i)+NumTrajBlock(j), i) = blocktrajid(1:NumTrajBlock(j),j)
         !      !write (fmiout, *) "Copied over ", blocktrajid(1:NumTrajBlock(j),j), " to ",&
         !      !                  blocktrajid( sum(NumTrajBlock(1:j)) : NumTrajBlock(i) + NumTrajBlock(j), i)

         !      ! Count how many positions need to be skipped
         !      skip_positions = 0
         !      do k = 1, i
         !         if (StoSelMode(k) == currStoSelMode) then
         !            skip_positions = skip_positions + ntrajblock(k)
         !         end if
         !      end do
         !      !write (fmiout, *) "This many positions ", skip_positions, "need 2 be skipped in col", i

         !      !prev_add_traj = 0
         !      !write (fmiout, *) "++++++++++++++++++++"
         !      do k = 1, ntrajblock(j)
         !         !write (fmiout, *) "Number of trajs added to column", i, "in a previous step: ", prev_add_traj
         !         !write (fmiout, *) "Number of trajs in previous blocks: ", NumTrajBlock(1:j-1), "in sum: ", &
         !         !                  sum(NumTrajBlock(1:j-1))
         !         !write (fmiout, *) "Now replace element ", skip_positions+k, "in ", &
         !         !                   i, "th column by", blocktrajid(k,j)
         !         blocktrajid(skip_positions+k, i) = blocktrajid(k,j)
         !         !prev_add_traj = prev_add_traj + 1
         !      end do
         !      !write (fmiout, *) "++++++++++++++++++++"

         !      ! Set j column to zero
         !      blocktrajid(:,j) = 0

         !      !write (fmiout, *) "Part of blocktrajid was modified to: ", blocktrajid(i,:)
         !      !write (fmiout, *) "The ith column is now: ", blocktrajid(:,i)

         !   end if

         !end do

      end do

   end subroutine FMS_getMergedBlockTrajID

   subroutine FMS_getNumCBFsBlock(B1, MergedBlockTrajID, NumTrajSepBlock, NumCBFsBlock, BlockCBFID)

      ! the original bundle
      type(T_TrajectoryBundle), intent(in) :: B1
      ! 2D array holding the TrajIDs sorted according to StoSel blocks they belong to
      integer(kind=DefInt), intent(in), dimension(B1%NumTraj + 1, 2) :: MergedBlockTrajID

      ! Getting info on Traj, CBF indices for each StoSel block:
      !  - Traj information
      !    1D array holding the number of trajs per StoSel block
      integer(kind=DefInt), intent(inout), dimension(2) :: NumTrajSepBlock
      !  - CBF information
      !    2D array holding the CBFIDs sorted according to StoSel blocks they belong to
      integer(kind=DefInt), intent(inout), dimension(B1%NCBFs + 1, 2) :: BlockCBFID
      !    1D array holding the number of trajs per StoSel block
      integer(kind=DefInt), intent(inout), dimension(2) :: NumCBFsBlock
      integer(kind=DefInt) :: CBFcount, CBFcurr
      ! Iterating over StoSel blocks, trajs:
      integer(kind=DefInt) :: iblocks, iTraj

      ! Determine number of trajs and CBFs in each StoSel block

      !NumTrajSepBlock = 0
      !NumCBFsBlock = 0
      BlockCBFID = 0
      CBFcurr = 0
      !CBFcount = 0

      ! Count how many trajs and CBFs we have in the current block
      !
      ! - Construct NumTrajSepBlock: will be needed 4 allocating temp Bundle and 
      !   knowing how many trajs 2 add to temp Bundles
      !
      ! - Construct BlockCBFID: equivalent of MergedBlockTrajID for CBFs, 2 keep track
      !   of CBFs corresponding 2 trajs
      ! - and NumCBFsBlock (will be needed 4 allocating temp Bundle and 
      !   making sure CBF IDs in temp Bundles start from 1)
      !
      ! - CBFcurr, CBFcount: used to make sure that multiplets are copied
      !   over correctly, and not separated

      do iblocks = 1, 2

         CBFcount = 0
         do iTraj =1, B1%NumTraj

            !write (fmiout, *) "Counting NCBFs in block ", iblocks
            !write (fmiout, *) "The current CBF is: ", CBFcurr
            !write (fmiout, *) "This is MergedBlockTrajID", MergedBlockTrajID

            ! count traj if its ID is in the current block
            ! (done in this way because we can not assume the TrajIDs to be
            !  in ascending order)
            if (any(MergedBlockTrajID(:, iblocks) == iTraj)) then
               NumTrajSepBlock(iblocks) = NumTrajSepBlock(iblocks) + 1
               if (B1%Trajectory(iTraj)%CBF /= CBFcurr) then
                  CBFcurr = B1%Trajectory(iTraj)%CBF
                  CBFcount = CBFcount + 1
                  BlockCBFID(CBFcount, iblocks) = CBFcurr
                  NumCBFsBlock(iblocks) = NumCBFsBlock(iblocks) + 1
               end if
            end if

         end do

         ! useful printing --> keep!
         !write (fmiout, *) "Number of Trajs and CBFs in block ", iblocks
         !write (fmiout, *) NumTrajSepBlock(iblocks), NumCBFsBlock(iblocks)

      end do

      !write (fmiout, *) "This is BlockCBFID :"
      !do iblocks = 1, BundlesMultDim
      !   write (fmiout, *) BlockCBFID(:, iblocks)
      !end do

   end subroutine FMS_getNumCBFsBlock

   subroutine FMS_getMultBundles(B1, sepB, iBlocks, MergedBlockTrajID, BlockCBFID, iNumTraj, iNumCBF)

      ! Number of separated bundles we need to create
      ! is determined based on merged blocktrajid
      ! so this number is also how many Selections we will perform
      !integer(kind=DefInt), intent(in) :: BundlesMultDim
      ! the original bundle
      type(T_TrajectoryBundle), intent(in) :: B1
      ! 1D array holding the TrajIDs sorted according to StoSel blocks they belong to
      integer(kind=DefInt), intent(in), dimension(B1%NumTraj + 1) :: MergedBlockTrajID
      ! 1D array holding the CBFIDs sorted according to StoSel blocks they belong to
      integer(kind=DefInt), intent(in), dimension(B1%NCBFs + 1) :: BlockCBFID
      ! array holding the separated bundles
      type(T_TrajectoryBundle), intent(inout) :: sepB(2)

      integer(kind=DefInt), intent(in) :: iBlocks, iNumTraj, iNumCBF
      integer(kind=DefInt) :: iTraj, currTrajID, currCBFID
      ! Number of particles
      integer(kind=DefInt) :: npart
      
      !! Renormalization of the separated Bundles:
      !real(kind=DefReal), dimension(B1%NumStates) :: OldNorms
      !real(kind=DefReal) :: DNorm0

      npart = B1%Trajectory(1)%NumParticles

      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Total number of trajs and CBFs in unseparated Bundle: "
      !write (fmiout, *) B1%NumTraj, B1%NCBFs
      !write (fmiout, *) " ----------------------------------------------------"

      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Number of trajs and CBFs that will be in unseparated Bundle: "
      !write (fmiout, *) "Number S, Number T trajs:", NumSingTraj, NumTripTraj
      !write (fmiout, *) "Number CBFs S, Number CBFs T trajs:", nSingCBF, nTripCBF
      !write (fmiout, *) " ----------------------------------------------------"

      ! (2) Create temporary bundles, separated according 2 selection mode:
      !
      ! - S only, or T only (as for S/T combined, no SS should happen
      !                      no temp Bundles are created 4 this case)
      !do iblocks = 1, BundlesMultDim

      call sepB(iBlocks)%create(numtraj=iNumTraj, &
                                       numdeadtraj=0, &                 !!!!! TODO: Is it fine setting this to 0???
                                       numstates=B1%NumStates, &        ! yes, I think it is because we are not taking
                                       numparticles=npart, &            ! into account dead trajectories in any of
                                       ncbfs=iNumCBF)                   ! the previous steps

      !end do

      ! useful printing --> keep!
      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Block number ", iblocks
      !write (fmiout, *) "Number trajs, Number CBFs:", sepB(iblocks)%NumTraj, sepB(iblocks)%NCBFs
      !!do iblocks = 1, BundlesMultDim
      !!   write (fmiout, *) "Block number ", iblocks
      !!   write (fmiout, *) "Number trajs, Number CBFs:", sepB(iblocks)%NumTraj, sepB(iblocks)%NCBFs
      !!end do
      !write (fmiout, *) " ----------------------------------------------------"

      ! (3) Copy over traj info from B1 to the temp bundles:
      !
      ! - Create trajs and copy over matching info from B1 according 2 MergedBlockTrajID column
      ! - Copy over matching Centroids (i.e. Centroid for the correct two trajs in the sepB)

      !write (fmiout, *) "The MergedBlockTrajID we are using for copying trajs to sep Bundles: "
      !do iTraj = 1, iNumTraj
      !   write (fmiout, *) "Traj number ", iTraj, "with ID ", MergedBlockTrajID(iTraj), &
      !                     "should have CBFID ", B1%Trajectory(MergedBlockTrajID(iTraj))%CBF
      !end do
      !write (fmiout, *) "The CBF indices: "
      !do iCBF = 1, iNumCBF
      !   write (fmiout, *) "CBF number ", iCBF, "with ID ", BlockCBFID(iCBF)
      !end do

      do iTraj = 1, iNumTraj

         ! create the traj in the temp separated bundle and set it 2 corresponding traj in B1
         currTrajID = MergedBlockTrajID(iTraj)
         currCBFID = B1%Trajectory(currTrajID)%CBF

         call sepB(iBlocks)%Trajectory(iTraj)%create(npart, B1%NumStates)
         sepB(iBlocks)%Trajectory(iTraj) = B1%Trajectory(currTrajID)
         sepB(iBlocks)%Trajectory(iTraj)%TrajID = B1%Trajectory(currTrajID)%TrajID ! should give the same as...
         !sepB(iblocks)%Trajectory(currTrajID)%TrajID = currTrajID                      ! this one here
                                                                                       
      end do

      !write (fmiout, *) "What we filled into sepB of block ", iBlocks
      !do iTraj = 1, iNumTraj

      !   write (fmiout, *) "Traj Number ", iTraj, "has ID", sepB(iBlocks)%Trajectory(iTraj)%TrajID
      !   write (fmiout, *) "has CBF ID ", sepB(iBlocks)%Trajectory(iTraj)%CBF

      !end do

      !write (fmiout, *) " ----------------------------------------------------"
      !write (fmiout, *) "Original centroids in B1"
      !write (fmiout, *) "Does B1 have smth like positions for Centroids?", B1%Centroids(1)%Particle(1)%get_pos()
      !do iTraj = 1, (((B1%NCBFs - 1) * B1%NCBFs) / 2)
      !   write (fmiout, *) "Centroid number ", iTraj
      !   write (fmiout, *) "is centroid to trajectories ", B1%Centroids(iTraj)%CentID
      !   write (fmiout, *) "And has position: ", B1%Centroids(iTraj)%Particle(1)%get_pos()
      !end do

      !do iblocks = 1, BundlesMultDim

      !   CenCount = 1
      !   CBFcount = 1

      !   do iTraj = 1, NumTrajSepBlock(iblocks)

      !      addTraj = MergedBlockTrajID(iTraj,iblocks)
      !      ! iTraj --> ID traj will have in temp bundle
      !      ! addTraj --> ID of corresponding traj in B1

      !      write (fmiout, *) "iTraj ", iTraj
      !      write (fmiout, *) "represents traj", addTraj
      !      write (fmiout, *) "The BlockCBFID for current block: ", BlockCBFID(:, iblocks)
      !      write (fmiout, *) "has actual CBF", B1%Trajectory(addTraj)%CBF, "the BlockCBFID CBF is:", BlockCBFID(CBFcount,iblocks)

      !      if (addTraj == 0) cycle ! dont try to add traj if there is no corresponding traj in BlockTrajID
      !                              ! TODO: do we actually need this?

      !      ! increase CBFcount if CBFcurr has changed from previous iteration
      !      ! (BlockCBFID(CBFcount,iblocks) is the CBF of previous the iteration)
      !      if (B1%Trajectory(addTraj)%CBF /= BlockCBFID(CBFcount,iblocks)) then
      !         CBFcount = CBFcount + 1
      !      end if

      !      ! create the traj in the temp separated bundle and set it 2 corresponding traj in B1
      !      call sepB(iblocks)%Trajectory(iTraj)%create(npart, B1%NumStates)
      !      sepB(iblocks)%Trajectory(iTraj) = B1%Trajectory(addTraj)
      !      sepB(iblocks)%Trajectory(iTraj)%TrajID = B1%Trajectory(addTraj)%TrajID

      !      ! check if CBF IDs need to be reset to 1 (the case if we are not in the first separated block)
      !      ! to make sure CBF IDs also start at 1
      !      ! TODO: change such that we always make sure CBFIDs start from one, not only when we are
      !      !       in blocks after the first block
      !      if (iblocks>1) then
      !         BlockCBFCount = 0
      !         do jblocks = 1, iblocks-1
      !            BlockCBFCount = BlockCBFCount + NumCBFBlock(jblocks)
      !         end do
      !         !write (fmiout, *) "BlockCBFcount (all CBFs added to previous blocks: ", BlockCBFCount
      !         sepB(iblocks)%Trajectory(iTraj)%CBF = B1%Trajectory(addTraj)%CBF - BlockCBFCount
      !         !write (fmiout, *) "CBF receives the number (must start at 1)", B1%Trajectory(addTraj)%CBF - BlockCBFCount
      !      else
      !         sepB(iblocks)%Trajectory(iTraj)%CBF = B1%Trajectory(addTraj)%CBF
      !      end if

      !      ! Next, copy over matching centroids from B1 to temp bundle
      !      ! ---------------------------------------------------------------------------------------
      !      ! Careful here, figuring this out has been a mess but we (we=Vera) pray that it works now
      !      ! - do we really cover all possible CBF ID pairs, i.e. catch all centroids?
      !      ! ---------------------------------------------------------------------------------------
      !      !write (fmiout, *) "___________________________"

      !      ! - start from the 1st CBF ID (i.e. 1st row, current column of BlockCBFID)
      !      !   and go over all CBF IDs as CBFcount increases
      !      CBF_i = BlockCBFID(CBFcount, iblocks)
      !      !write (fmiout, *) "CBF_i is, CBF_j will be", CBF_i, CBF_i + 1
      !      !write (fmiout, *) "loop limit will be:", B1%NCBFs

      !      ! - loop over possible CBF ID pairs and check 4 matching CentID pair in B1
      !      is_addCen = .false.
      !      do CBF_j = CBF_i + 1, B1%NCBFs ! what if CBF IDs are not in ascending order, does this still work??? :S

      !         if (CenCount > (((sepB(iblocks)%NCBFs - 1) * sepB(iblocks)%NCBFs) / 2) ) cycle ! dont add more Cens
      !                                                                                                      ! if max of cens in
      !                                                                                                      ! temp bundle is reached

      !         !write (fmiout, *) "CBF_i and CBF_j", CBF_i, CBF_j
      !         if (any(CBF_i == BlockCBFID(:,iblocks)) .and. any(CBF_j == BlockCBFID(:,iblocks)) .and. CBF_i /= CBF_j) then

      !             !write (fmiout, *) "dream come true", CBF_j, CBF_i
      !             !write (fmiout, *) "CenID is", CenID

      !             ! For this CentID pair, find out corresponding ID of Centroid in B1 (CenID)
      !             do Cen_TrajID = 1, (((B1%NCBFs - 1) * B1%NCBFs) / 2)

      !                !write(fmiout, *) "Is Cen ", Cen_TrajID, "the correct one?"
      !                !write (fmiout, *) B1%Centroids(Cen_TrajID)%CentID
      !                Cen_i = B1%Centroids(Cen_TrajID)%CentID(1); Cen_j = B1%Centroids(Cen_TrajID)%CentID(2)

      !                if (CBF_i == Cen_i .and. CBF_j == Cen_j .or. CBF_j == Cen_i .and. CBF_i == Cen_j) then
      !                   CenID = Cen_TrajID
      !                   !write (fmiout, *) "The ID of the centroid we need is ", CenID
      !                   is_addCen = .true.
      !                   !cycle
      !                end if

      !             end do

      !             if (is_addCen .eqv. .true.) then
      !                ! Use the corresponding ID of Centroid in B1 (CenID) to copy over matching info
      !                sepB(iblocks)%Centroids(CenCount) = B1%Centroids(CenID)
      !                sepB(iblocks)%Centroids(CenCount)%CentID = [CBF_j, CBF_i]
      !                sepB(iblocks)%Centroids(CenCount)%TrajID = CenCount
      !                !write (fmiout, *) "We added Cen", CenID, "of B1 to block Bundle ", iblocks, "as Cen number", CenCount

      !                if (iblocks>1) then
      !                   !BlockCBFCount = 0
      !                   !do jblocks = 1, iblocks-1
      !                   !   BlockCBFCount = BlockCBFCount + NumCBFBlock(jblocks)
      !                   !end do

      !                   ! still in the same block, so reuse BlockCBFCount from above

      !                   sepB(iblocks)%Centroids(CenCount)%CentID = [CBF_j - BlockCBFCount, &
      !                                                                      CBF_i - BlockCBFCount]
      !                   !write (fmiout, *) "Cen receives the number (must start at 1)"
      !                end if
      !                CenCount = CenCount + 1
      !                cycle
      !             end if

      !         !else
      !         !   write (fmiout, *) "Index pair is not the one needed for the CentIDs in this block"

      !         end if

      !      end do

      !      !write (fmiout, *) "___________________________"

      !   end do
      !
      !end do

      !write (fmiout, *) "Do we have Centroids after copying them?", glzCentroids
      !do iblocks = 1, BundlesMultDim
      !   if (sepB(iblocks)%NCBFs > 1) then
      !      do iTraj = 1, (((sepB(iblocks)%NCBFs - 1) * sepB(iblocks)%NCBFs) / 2)
      !         write (fmiout, *) "Centroid number ", iTraj
      !         write (fmiout, *) "is centroid to trajectories ", sepB(iblocks)%Centroids(iTraj)%CentID
      !         write (fmiout, *) "And has position: ", sepB(iblocks)%Centroids(iTraj)%Particle(1)%get_pos()
      !      end do
      !   else
      !      write (fmiout, *) "This block only has one CBF: ", iblocks, ". Therefore no centroids."
      !   end if
      !end do

      !write (fmiout, *) " ----------------------------------------------------"

      !write (fmiout, *) "What we filled into the separated bundle: "
      !do iblocks = 1, BundlesMultDim
      !   do iTraj = 1, sepB(iblocks)%NumTraj

      !         write (fmiout, *) "Info for trajectory with ID ", sepB(iblocks)%Trajectory(iTraj)%TrajID
      !         write (fmiout, *) "Triplet?", sepB(iblocks)%Trajectory(iTraj)%triplet
      !         write (fmiout, *) "Ms", sepB(iblocks)%Trajectory(iTraj)%Ms
      !         write (fmiout, *) "CBF", sepB(iblocks)%Trajectory(iTraj)%CBF

      !   end do
      !end do
      !write (fmiout, *) " ----------------------------------------------------"

   end subroutine FMS_getMultBundles

! Copy over trajectories after the selection from sepB
! to the original bundle. Note that we must copy over
! also the trajectories marked for death.
   subroutine FMS_copy_sepB_to_original_bundle(B1, sepB, iBlocks, NumTrajiBlocks, &
                                                           MergedBlockTrajID) !, MergedSepBlockTrajID)
      integer(kind=DefInt), intent(in) :: iBlocks
      integer(kind=DefInt), intent(in) :: NumTrajiBlocks
      type(T_TrajectoryBundle), intent(inout) :: B1
      type(T_TrajectoryBundle), intent(in) :: sepB
      integer(kind=DefInt), dimension(B1%NumTraj + 1, 2), intent(in) :: &
                                                   MergedBlockTrajID !, MergedSepBlockTrajID
      integer(kind=DefInt) :: iTraj, jTraj
      integer(kind=DefInt) :: add_DeadTraj, iDeadTraj
                             !iSingTraj, iTripTraj, iCBF
      !integer(kind=DefInt) :: TrajIDSing, CBFIDSing, TrajIDTrip, CBFIDTrip

      ! I dont actually think we need to copy back the temp separated bundles
      ! All we need to do is to adjust gliForceKill to contain the TrajIDs for B1
      ! not those referring to the temp separated bundles

      !do i = 1, B1%NumTraj
      !   write (fmiout, *) "These trajs are in the original B1 (before copying back)"
      !   write (fmiout, *) "TrajID: ", B1%Trajectory(i)%TrajID, "DeadTime: ", B1%Trajectory(i)%DeadTime
      !   write (fmiout, *) "IsDead? (checking if DeadTime < CurrentTime)", B1%Trajectory(i)%DeadTime < B1%CurrentTime
      !end do

      ! How to get the information on which Trajs are dead?
      !write (fmiout, *) "These trajs are dead in the original B1 (before copying back)"
      !do i = 1, B1%NumDeadTraj
      !   write (fmiout, *) "TrajID: ", B1%DeadTraj(i)%TrajID, "DeadTime: ", B1%DeadTraj(i)%DeadTime
      !end do

      ! Try a different approach to what is (or will be) commented out further down

      !write (fmiout, *) "i (the how manyth block)= ", iBlocks
      do iTraj = 1, NumTrajiBlocks

         !write (fmiout, *) "iTraj (the how manyth traj)= ", iTraj
         !write(fmiout, *) "TrajID in sepB: ", sepB%Trajectory(iTraj)%TrajID

         add_DeadTraj = 0

         do iDeadTraj = 1, B1%NumDeadTraj
            if (B1%DeadTraj(iDeadTraj)%TrajID <= sepB%Trajectory(iTraj)%TrajID) then
               add_DeadTraj = add_DeadTraj + 1
            end if
         end do

         !write (fmiout, *) "add_DeadTraj is now ", add_DeadTraj

         do jTraj = 1, B1%NumTraj + B1%NumDeadTraj
            if (sepB%Trajectory(iTraj)%TrajID == jTraj) then

               !write (fmiout, *) "It's a match! for ", sepB%Trajectory(iTraj)%TrajID, &
               !                  "(sepB TrajID)", "and", jTraj, "(B1 TrajID)"

               !write (fmiout, *) "The amplitude ", B1%Trajectory(jTraj - add_DeadTraj)%Amplitude, &
               !                  "(the ", jTraj - add_DeadTraj, "th traj in B1)", &
               !                  "will be overwritten to", sepB%Trajectory(iTraj)%Amplitude
               B1%Trajectory(jTraj - add_DeadTraj)%Amplitude = sepB%Trajectory(iTraj)%Amplitude

            end if
         end do

         !do jTraj = 1, B1%NumTraj + B1%NumDeadTraj

            !add_DeadTraj = 0

            !if (sepB%Trajectory(iTraj)%TrajID == jTraj) then
            !   write (fmiout, *) "It's a match! for ", sepB(i)%Trajectory(iTraj)%TrajID, &
            !                     "(sepB TrajID)", "and", jTraj, "(B1 TrajID)"

            !   do iDeadTraj = 1, B1%NumDeadTraj
            !      if (B1%DeadTraj(iDeadTraj)%TrajID < jTraj) then
            !         add_DeadTraj = add_DeadTraj + 1
            !      end if
            !   end do

            !   write (fmiout, *) "The amplitude ", B1%Trajectory(jTraj - add_DeadTraj)%Amplitude, &
            !                     "will be overwritten to", sepB(i)%Trajectory(iTraj)%Amplitude
            !   B1%Trajectory(jTraj - add_DeadTraj)%Amplitude = sepB(i)%Trajectory(iTraj)%Amplitude
            !end if

         !end do

      end do

      !DeadCount = 0
      !DeadTrajIDs = 0
      !do iTraj = 1, B1%NumTraj ! Does this loop always account for all dead and alive trajs???
      !                         ! Even after several kills in the Bundle?
      !                         ! Do we need to add B1%NumDeadTraj somehow asp?
      !   !write (fmiout, *) " ------------------------"
      !   !write (fmiout, *) "Current iTraj: ", iTraj

      !   !if (any(B1%Trajectory(:)%TrajID == iTraj)) then
      !   !   write(fmiout, *) "This traj is alive"
      !   !end if

      !   if (any(B1%DeadTraj(:)%TrajID == iTraj)) then
      !      !write(fmiout, *) "This traj is dead"
      !      ! This TrajID needs to be remembered and 'left blank' when temp separated Bundles are copied back
      !      DeadCount = DeadCount + 1
      !      DeadTrajIDs(DeadCount) = iTraj
      !   end if

      !   if (any(B1%DeadTraj(:)%TrajID == iTraj) .and. any(B1%Trajectory(:)%TrajID == iTraj)) then
      !      !write(fmiout, *) "This traj is dead and alive"
      !      call FMS_DieError("Traj is dead and alive. This shouldn't happen")
      !   end if

      !   !write (fmiout, *) " ------------------------"
      !end do
      !!write (fmiout, *) "These are our dead trajectories: ", DeadTrajIDs

      !hitDead = .false.

      !do i = 1, BundlesMultDim
      !   if (StoSelMode(i) == 1 .or. StoSelMode(i) == 3) then
      !      !write (fmiout, *) "Killing happened in block Bundle ", i
      !      !write (fmiout, *) "Trajs need to be copied back to B1 to keep track of who died"

      !      do iTraj = 1, sepB(i)%NumTraj
      !      !   write (fmiout, *) " ------------------------"
      !      !   write (fmiout, *) "Copying back Traj ", iTraj
      !      !   write (fmiout, *) "Ms", sepB(i)%Trajectory(iTraj)%Ms
      !      !   write (fmiout, *) "CBF", sepB(i)%Trajectory(iTraj)%CBF
      !      !   write (fmiout, *) "What is the amplitude?", sepB(i)%Trajectory(iTraj)%Amplitude
      !         write (fmiout, *) " ------------------------"
      !         write (fmiout, *) "Is this traj still the same in B1?", B1%Trajectory(iTraj)%TrajID
      !         write (fmiout, *) "Ms", B1%Trajectory(iTraj)%Ms
      !         write (fmiout, *) "CBF", B1%Trajectory(iTraj)%CBF
      !         write (fmiout, *) "What is the amplitude?", B1%Trajectory(iTraj)%Amplitude
      !         write (fmiout, *) " ------------------------"
      !      end do
      !      !write (fmiout, *) " ----------------------------------------------------"
      !      !write (fmiout, *) "Now, actually copy and check again: "
      !      !write (fmiout, *) " ----------------------------------------------------"

      !      ! Copy back trajectories, adjust temp bundle TrajID, CBF back to refer to B1 
      !      !write (fmiout, *) "Copying trajs for block bundle ", i
      !      
      !      ! 20/05/26: Try something new -----------------------------------------------------

      !      write (fmiout, *) "i (the how manyth block)= ", i

      !      do j = 1, sepB(i)%NumTraj

      !         write (fmiout, *) "j (the how manyth traj)= ", j
      !         add_DeadTraj = 0

      !         do iDeadTraj = 1, DeadCount

      !            if (MergedBlockTrajID(j,i) == DeadTrajIDs(iDeadTraj)) then
      !               add_DeadTraj = iDeadTraj
      !               write (fmiout, *) "Hit the dead traj ", DeadTrajIDs(iDeadTraj)
      !            end if

      !         end do

      !         write (fmiout, *) "add_DeadTraj is ", add_DeadTraj

      !         !B1TrajID = MergedBlockTrajID(j,i)+add_DeadTraj
      !         !sepBTrajID = MergedSepBlockTrajID(j,i)
      !         !write (fmiout, *) "Are these really the same as TrajIDs? -----------"
      !         !write (fmiout, *) "B1TrajID", B1TrajID, "actual TrajID", &
      !         !                  B1%Trajectory(MergedBlockTrajID(j,i)+add_DeadTraj)%TrajID
      !         !write (fmiout, *) "sepBTrajID", sepBTrajID, "actual TrajID", &
      !         !                  sepB(i)%Trajectory(MergedSepBlockTrajID(j,i))%TrajID
      !         !write (fmiout, *) "--------------"

      !         B1TrajID = B1%Trajectory(MergedBlockTrajID(j,i)+add_DeadTraj)%TrajID
      !         sepBTrajID = sepB(i)%Trajectory(MergedSepBlockTrajID(j,i))%TrajID
      !         write (fmiout, *) "Are these really the same as TrajIDs? -----------"
      !         write (fmiout, *) "B1TrajID", B1TrajID, "actual TrajID", &
      !                           B1%Trajectory(MergedBlockTrajID(j,i)+add_DeadTraj)%TrajID
      !         write (fmiout, *) "sepBTrajID", sepBTrajID, "actual TrajID", &
      !                           sepB(i)%Trajectory(MergedSepBlockTrajID(j,i))%TrajID
      !         write (fmiout, *) "--------------"

      !         write (fmiout, *) "Attempting to copy amplitude of sep Bundle traj ", sepBTrajID, &
      !                           "over to traj ", B1TrajID, "in B1"
      !         !write (fmiout, *) "The corresponding TrajIDs are: ", B1%Trajectory(MergedBlockTrajID(j,i)+add_DeadTraj)%TrajID, &
      !         !                  "in B1 and ", sepB(i)%Trajectory(MergedSepBlockTrajID(j,i))%TrajID, "in sep Bundle"
      !         write (fmiout, *) "The amplitude ", B1%Trajectory(B1TrajID)%Amplitude, &
      !                           "will be overwritten to", sepB(i)%Trajectory(sepBTrajID)%Amplitude
      !         B1%Trajectory(B1TrajID)%Amplitude = &
      !                 sepB(i)%Trajectory(sepBTrajID)%Amplitude

      !      end do

      !      ! 20/05/26 end --------------------------------------------------------------------

      !      !do iTraj = 1, sepB(i)%NumTraj
      !      !   if (i>1) then
      !      !      B1%Trajectory(iTraj+TrajCount) = sepB(i)%Trajectory(iTraj) !TODO: figure out what original TrajID, CBF ID was

      !      !      if (any(sepB(i)%Trajectory(iTraj)%TrajID + sepB(i-1)%NumTraj == DeadTrajIDs(:))) then
      !      !         write (fmiout, *) "We hit the dead traj (iTraj)", &
      !      !                           sepB(i)%Trajectory(iTraj)%TrajID + sepB(i-1)%NumTraj
      !      !         write (fmiout, *) "Adding the DeadCount ", DeadCount, "from now on"
      !      !         hitDead = .true.
      !      !      end if

      !      !      if (hitDead) then
      !      !         B1%Trajectory(iTraj+TrajCount)%TrajID = sepB(i)%Trajectory(iTraj)%TrajID + sepB(i-1)%NumTraj + &
      !      !         DeadCount
      !      !         !B1%Trajectory(iTraj+TrajCount)%CBF = sepB(i)%Trajectory(iTraj)%CBF
      !      !         ! TODO: 30/04/26: Do we have to adjust CBF IDs just like for TrajID??
      !      !         B1%Trajectory(iTraj+TrajCount)%Amplitude = sepB(i)%Trajectory(iTraj)%Amplitude  ! used to break stuff
      !      !      else
      !      !         B1%Trajectory(iTraj+TrajCount)%TrajID = sepB(i)%Trajectory(iTraj)%TrajID + sepB(i-1)%NumTraj
      !      !         B1%Trajectory(iTraj+TrajCount)%Amplitude = sepB(i)%Trajectory(iTraj)%Amplitude  ! used to break stuff
      !      !         !B1%Trajectory(iTraj+TrajCount)%CBF = sepB(i)%Trajectory(iTraj)%CBF
      !      !      end if
      !      !      B1%Trajectory(iTraj+TrajCount)%CBF = sepB(i)%Trajectory(iTraj)%CBF + sepB(i-1)%NCBFs
      !      !      write (fmiout, *) "Traj of index ", iTraj, "was copied over and received ID, ", &
      !      !                        B1%Trajectory(iTraj+TrajCount)%TrajID

      !      !   ! copy over amplitudes for first block as well:
      !      !   ! as it's the first block we don't need to worry about:
      !      !   ! - how many trajs we had in previous blocks (for sure)
      !      !   ! - how many dead trajs (not too sure, but we pray for it each night)
      !      !   else if (i==1) then
      !      !      B1%Trajectory(iTraj)%Amplitude = sepB(i)%Trajectory(iTraj)%Amplitude
      !      !   end if

      !      !end do

      !      ! 23/04/26: Definitely go back here and check that copying back is done correctly
      !      !           also apart from amplitudes

      !      do iTraj = 1, sepB(i)%NumTraj
      !         !write (fmiout, *) " ------------------------"
      !         !write (fmiout, *) "Copying back Traj ", iTraj
      !         !write (fmiout, *) "Ms", sepB(i)%Trajectory(iTraj)%Ms
      !         !write (fmiout, *) "CBF", sepB(i)%Trajectory(iTraj)%CBF
      !         !write (fmiout, *) "What is the amplitude?", sepB(i)%Trajectory(iTraj)%Amplitude
      !         write (fmiout, *) " ------------------------"
      !         write (fmiout, *) "Is this traj still the same in B1?", B1%Trajectory(iTraj)%TrajID
      !         write (fmiout, *) "Ms", B1%Trajectory(iTraj)%Ms
      !         write (fmiout, *) "CBF", B1%Trajectory(iTraj)%CBF
      !         write (fmiout, *) "What is the amplitude?", B1%Trajectory(iTraj)%Amplitude
      !         write (fmiout, *) " ------------------------"
      !      end do
      !   end if
      !end do

   end subroutine FMS_copy_sepB_to_original_bundle

   subroutine fill_state_bundles(B1, BundleSS)
      type(T_TrajectoryBundle), intent(in) :: B1
      type(T_TrajectoryBundle), intent(inout) :: BundleSS(B1%NumStates)
      integer(kind=DefInt) :: iState, iTraj, nStateTraj, iSSTraj, npart

      npart = B1%Trajectory(1)%NumParticles

      do iState = 1, B1%NumStates

         ! Count the number of trajectories for a given state
         nStateTraj = 0
         do iTraj = 1, B1%NumTraj
            if (B1%Trajectory(iTraj)%StateID == iState) then
               nStateTraj = nStateTraj + 1
            end if
         end do

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
         end do
      end do

   end subroutine fill_state_bundles

! Copy over trajectories after the selection from BundleSS
! to the original bundle. Note that we must copy over
! also the trajectories marked for death.
   subroutine copy_state_bundles_to_original_bundle(B1, BundleSS)
      type(T_TrajectoryBundle), intent(inout), target :: B1
      type(T_TrajectoryBundle), intent(in), target :: BundleSS(B1%NumStates)
      type(T_TrajectoryBundle), pointer :: BSS_i
      type(T_Trajectory), pointer :: T_i
      integer(kind=DefInt) :: iState, iTrj, iSSTraj

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
            end do

         end do

      end do

   end subroutine copy_state_bundles_to_original_bundle

end module SelectionModule
