!     Copyright Todd J. Martinez and Raphael D. Levine, 1994

!>
!! Mark any trajectories to be removed as 'dead', then remove them
!! and rescale the bundle
!!
!! Trajectories are marked for deletion if they are on IgnoreState,
!! or have population less than PopToSpawn.
!! Dead trajectories are removed from the Bundle%Trajectory array
!! and appended to the Bundle%DeadTraj array.
!! The bundle matrices are rescaled accordingly.
!<
subroutine FMS_RemoveDead(B1)
   use GlobalModule, only: defInt, defReal, fmiOut, FPZero, &
                           glIgnoreState, glzCentroids, gliForceKill
   use BundleModule
   use BundleCalcsModule, only: FMS_bH, FMS_Mulliken
   use SpawnModule, only: spawn_couple, spawn_params
   implicit none
   type(T_TrajectoryBundle), intent(inout) :: B1
   type(T_TrajectoryBundle) :: BTemp

   integer(kind=DefInt) :: iTraj, nDead, iLive, iDead, jTraj
   integer(kind=DefInt) :: iCent, jCent, jDead, jLive
   integer(kind=DefInt) :: iCBF, jCBF, nCBF, n2CBF

   integer(DefInt), dimension(B1%NumTraj, B1%NumTraj) :: Coupled, Coupled_prev

   real(DefReal) :: Population(B1%NumTraj)
   real(kind=DefReal) :: pop

   integer(DefInt) :: ntraj, i, j, TrajID, StateID, nstate, n

   logical :: NotCoupled, OnIgnoreState, PopBelowThresh, &
              MarkForDeath, OverSpawnThresh, ForceDead

   ntraj = B1%NumTraj
   nstate = B1%NumStates

   ! First, decompose the hamiltonian into a block diagonal representation
   ! work out the coupling matrix
   Coupled = 0
   do i = 2, ntraj
      do j = 1, i
         if (abs(FMS_bH(B1, i, j)) > FPZero) then
            Coupled(i, j) = 1
            Coupled(j, i) = 1
         end if
      end do
   end do

   ! This matrix has the properties that it is non-zero is
   !    Coup^1   : directy connected
   !    Coup^2   : connected by 1 or less common trajectory
   !    Coup^3   : connected by 2 or less common trajectories
   ! We will iterate the matrix multiplication to convergence
   do
      Coupled_prev = Coupled

      Coupled = min(matmul(Coupled, Coupled), 1)

      if (all(Coupled == Coupled_prev)) exit
   end do

   ! Replace the diagonal with zero
   do i = 1, ntraj
      Coupled(i, i) = 0
   end do

   Population = abs(FMS_Mulliken(B1))

   ndead = 0
   do i = 1, ntraj

      TrajID = B1%Trajectory(i)%TrajID
      StateID = B1%Trajectory(i)%StateID

      OverSpawnThresh = .false.
      do n = 1, nstate
         OverSpawnThresh = OverSpawnThresh .or. spawn_couple(B1%Trajectory(i), n) > spawn_params%CFThresh
      end do

      ! Update DeadTime for those who will not be killed
      NotCoupled = all(Coupled(:, i) == 0)
      OnIgnoreState = (StateID == glIgnoreState)
      PopBelowThresh = (Population(i) < spawn_params%PopToSpawn)

      ! GAIMS changed
      pop = 0.d0
      do n = 1, ntraj
         if (B1%Trajectory(n)%CBF == B1%Trajectory(i)%CBF) then
            pop = pop + Population(n)
         end if
      end do
      PopBelowThresh = (pop < spawn_params%PopToSpawn)

      MarkForDeath = ((OnIgnoreState .and. NotCoupled .and. .not. OverSpawnThresh) &
                      .or. (PopBelowThresh .and. NotCoupled))
      ForceDead = any(gliForceKill(:) == TrajID)
      if (ForceDead) then
         B1%Trajectory(i)%DeadTime = -9999.0d0
      end if

      MarkForDeath = (MarkForDeath .or. ForceDead)

      if (.not. MarkForDeath) then
         B1%Trajectory(i)%DeadTime = B1%CurrentTime
      else
         write (fmiOut, '(a,i0,a)') 'Traj ', TrajID, ' was marked for death'
      end if

      ! Workout who is getting killed
      if (B1%Trajectory(i)%is_dead()) then
         ndead = ndead + 1
         if (ForceDead) then
            write (fmiOut, '(a,i0,a,i0)') '** Force Killing trajectory ', TrajID, ' on state ', StateID
         else if (OnIgnoreState) then
            write (fmiOut, '(a,i0,a,i0)') '** Killing trajectory ', TrajID, ' on state ', StateID
         else
            write (fmiOut, '(a,i0,a,f6.5)') '** Killing trajectory ', TrajID, ' pop < ', spawn_params%PopToSpawn
         end if
      end if
   end do

   ! Remove dead trajectories
   if (ndead > 0) then
      call BTemp%create(numtraj=B1%NumTraj - nDead, &
                        numdeadtraj=B1%NumDeadTraj + nDead, &
                        numstates=B1%NumStates, &
                        numparticles=B1%NumParticles, &
                        ncbfs=B1%NCBFs)
      BTemp%CurrentTime = B1%CurrentTime

      ! Copy old dead trajectories to new bundle
      do iTraj = 1, B1%NumDeadTraj
         call BTemp%DeadTraj(iTraj)%copy_from(B1%DeadTraj(iTraj))
         do jTraj = 1, B1%NumDeadTraj
            BTemp%DeadH(iTraj, jTraj) = B1%DeadH(iTraj, jTraj)
         end do
      end do

      iLive = 0
      iDead = 0
      nCBF = 0
      write (fmiOut, '(a,/,a)') 'Reconstructing trajectory bundle. Living trajectories:'

      do iTraj = 1, B1%NumTraj
         if (.not. B1%Trajectory(iTraj)%is_dead()) then
            ! Alive: copy trajectory to new bundle
            iLive = iLive + 1
            call BTemp%Trajectory(iLive)%copy_from(B1%Trajectory(iTraj))

            ! Get the CBF identifier. Here we use the fact that
            ! the 3 triplet trajecories are always stored in
            ! successive order starting with Ms=2
            if (BTemp%Trajectory(iLive)%Ms == 2) then
               nCBF = nCBF + 1
            end if
            if (B1%Trajectory(iTraj)%triplet) then
               write (fmiOut, '(4X,I4," S=1 Ms=",I0," CBF ",I0)') iTraj, BTemp%Trajectory(iLive)%Ms - 2, nCBF
            else
               write (fmiOut, '(4X,I4," S=0 Ms=",I0," CBF ",I0)') iTraj, 0, nCBF
            end if
            BTemp%Trajectory(iLive)%CBF = nCBF

            ! Copy centroids over too
            if (glzCentroids) then
               jLive = 0
               if (BTemp%Trajectory(iLive)%Ms == 2) then
                  n2CBF = 0
                  do jTraj = 1, iTraj - 1
                     if (.not. B1%Trajectory(jTraj)%is_dead() .and. (B1%Trajectory(jTraj)%Ms == 2)) then
                        n2CBF = n2CBF + 1
                        iCBF = B1%Trajectory(iTraj)%CBF
                        jCBF = B1%Trajectory(jTraj)%CBF
                        iCent = ((iCBF - 2) * (iCBF - 1)) / 2 + jCBF
                        jCent = ((nCBF - 2) * (nCBF - 1)) / 2 + n2CBF
                        call BTemp%Centroids(jCent)%copy_from(B1%Centroids(iCent))
                        BTemp%Centroids(jCent)%CentID(1) = nCBF
                        BTemp%Centroids(jCent)%CentID(2) = n2CBF
                     end if
                  end do
               end if
            end if

         else
            ! Dead: add the trajectory to the graveyard
            iDead = iDead + 1
            call BTemp%DeadTraj(B1%NumDeadTraj + iDead)%copy_from(B1%Trajectory(iTraj))

            ! Set dead time of recently killed trajectory to be current time
            ! to make it easier for user to restart killed trajectories!
            BTemp%DeadTraj(B1%NumDeadTraj + iDead)%DeadTime = BTemp%CurrentTime
            ! Copy new hamiltonian elements into dead hamiltonian
            jDead = 0
            do jTraj = 1, iTraj
               if (B1%Trajectory(jTraj)%is_dead()) then
                  jDead = jDead + 1
                  BTemp%DeadH(B1%NumDeadTraj + iDead, B1%NumDeadTraj + jDead) = FMS_bH(B1, iTraj, jTraj)
                  BTemp%DeadH(B1%NumDeadTraj + jDead, B1%NumDeadTraj + iDead) = FMS_bH(B1, jTraj, iTraj)
                  BTemp%DeadH(B1%NumDeadTraj + jDead, 1:B1%NumDeadTraj) = (0.d0, 0.d0)
                  BTemp%DeadH(1:B1%NumDeadTraj, B1%NumDeadTraj + jDead) = (0.d0, 0.d0)
               end if
            end do
         end if
      end do

      BTemp%NCBFs = nCBF
      call B1%copy_from(BTemp)
      call BTemp%destroy()

   end if

end subroutine FMS_RemoveDead
