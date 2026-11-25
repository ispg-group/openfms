!     Copyright Todd J. Martinez and Raphael D. Levine, 1994
!>
!!    Free memory, close interfaces, and quit succesfully.
!!    This is for succesful exits - for errors, use FMS_DieError.
!!    \param B1 (Optional) trajectory bundle to shut down cleanly
!!    \param T1 (Optional) trajectory to shut down cleanly
!!    @ingroup ESP
!<
subroutine FMS_Shutdown(B1, T1, terminate)
   use GlobalModule
   use BundleModule
   use TrajectoryModule
   implicit none

   type(T_TrajectoryBundle), optional, intent(inout) :: B1
   type(T_Trajectory), optional, intent(inout) :: T1
   logical, optional, intent(in) :: terminate
   integer(kind=DefInt) :: NumTraj

!     Free memory
   if (present(B1)) NumTraj = B1%NumTraj
   if (present(B1)) call B1%destroy()
   if (present(T1)) call T1%destroy()

!     End MPI Comm
   if (gliModel == glModelType%TC .and. associated(tc_finalize)) then
      call tc_finalize()
   end if

   ! Write timing info
   call cpu_time(glTimers%final)
   write (fmiOut, *)
   write (fmiOut, *) '== FMS DONE =='
   write (fmiOut, *)
   write (fmiOut, *) '-------------- Timing Info: -------------------'
   write (fmiOut, 1) '# electronic structure calls:  ', glNESCalls
   write (fmiOut, 1) '# rejected timesteps:          ', glnStepsRejected
   write (fmiOut, 2) 'Spawning - Forward propagation: ', glTimers%forwardprop
   write (fmiOut, 2) 'Spawning - Backward propagation:', glTimers%backprop
   if (GlzIterInv) then
      write (fmiOut, 2) 'Conjugate gradient time:       ', glTimers%cg
   else
      write (fmiOut, 2) 'S inversion time:              ', glTimers%inv
   end if
   write (fmiOut, 2) 'BuildHS time:                  ', glTimers%bhs
   write (fmiOut, 2) 'Overlap time:                  ', glTimers%olap
   write (fmiOut, 2) 'V matrix element time:         ', glTimers%pot
   write (fmiOut, 2) 'Electronic structure time:     ', glTimers%elecstruct
   write (fmiOut, 2) 'Total runtime:                 ', glTimers%final - glTimers%initial
   if (present(B1)) then
      write (fmiOut, 1) 'Final # of trajectories:       ', NumTraj
   end if
   flush (fmiOut)
   close (fmiOut)

1  format(A31, 1x, I8)
2  format(A31, 1x, F11.2)

!     Quit succesfully
!     This may be a little confusing, but want an option to return
   if (present(terminate)) then
      if (terminate) stop
   else
      stop
   end if

end subroutine FMS_Shutdown
