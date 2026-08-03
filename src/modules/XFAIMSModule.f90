module XFAIMSModule
   use GlobalModule
   use TrajectoryModule
   implicit none

   !> XFAIMS parameters
   !> Default values are set in read_namelist.f90 as for the other parameters.
   type :: t_xfaims_params
      real(DefReal) :: f0
      real(DefReal) :: freq
      real(DefReal) :: sigma
      real(DefReal) :: CEP
      real(DefReal) :: t0
      real(DefReal) :: polx
      real(DefReal) :: poly
      real(DefReal) :: polz
      real(DefReal) :: CoupFieldTimeStep
      real(DefReal) :: sp_spwn_i
      real(DefReal) :: sp_spwn_f
      logical :: onespawnonly
      logical :: IgnoreStateAferField
!   contains
!      procedure, public :: initialize => initialize_xfaims_params
   end type t_xfaims_params

   type(t_xfaims_params) :: xfaims_params

   private
   public :: xfaims_params, print_xfaims_params, activate_xfaims

contains

   subroutine print_xfaims_params()
      implicit none

      character(len=*), parameter :: divider = &
                                     ' -------------------------------------------------------'

      write (fmiOut, *)
      write (fmiOut, '(a)') divider
      write (fmiOut, '(a)') ' XFAIMS parameters'
      write (fmiOut, '(a)') divider

      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Field amplitude f0_xfr:', xfaims_params%f0, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Field frequency freq_xfr:', xfaims_params%freq, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Pulse width sigma_xfr:', xfaims_params%sigma, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8)') 'Carrier-envelope phase:', xfaims_params%CEP
      write (fmiOut, '(1x,a,t34,f16.4,a)') 'Pulse center t0_xfr:', xfaims_params%t0, ' a.u.'

      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization x:', xfaims_params%polx
      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization y:', xfaims_params%poly
      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization z:', xfaims_params%polz

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Field coupling timestep:', xfaims_params%CoupFieldTimeStep, ' a.u.'

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Spawning interval start:', xfaims_params%sp_spwn_i, ' a.u.'

      ! XFAIMS is only activated at multiples of the global timestep in openfms.F90.
      if (modulo(xfaims_params%sp_spwn_i, gldTimeStep) /= 0.0) then
         write (fmiOut, '(1x,a)') &
            'WARNING: XFAIMS will be initiated at the closest following multiple of TimeStep.'
      end if

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Spawning interval end:', xfaims_params%sp_spwn_f, ' a.u.'

      ! XFAIMS is only deactivated at multiples of the global timestep in openfms.F90.
      if (modulo(xfaims_params%sp_spwn_f, gldTimeStep) /= 0.0) then
         write (fmiOut, '(1x,a)') &
            'WARNING: XFAIMS will be terminated at the closest following multiple of TimeStep.'
      end if

      write (fmiOut, '(1x,a,t34,l16)') &
         'Only one field spawning:', xfaims_params%onespawnonly

      write (fmiOut, '(1x,a,t34,l16)') &
         'Ignore state after field:', xfaims_params%IgnoreStateAferField

      write (fmiOut, '(a)') divider
      write (fmiOut, *)

      flush (fmiOut)

   end subroutine print_xfaims_params

   subroutine activate_xfaims(time)
      real(DefReal), intent(in) :: time

      ! This is the point where we really decide if we are in the field region or not and set the corresponding timestep.
      ! Hence, we set the global variable determinig the field coupling here.
      if (glzxfaims) then
         if (time >= xfaims_params%sp_spwn_i .and. time < xfaims_params%sp_spwn_f) then
            ! Activate XFAIMS
            if (.not. glzXFActive) then
               glzXFActive = .true.
               write (fmiOut, *)
               write (fmiOut, *) '----------------------------------------'
               write (fmiOut, *) 'Entering field coupling region => XFAIMS'
               write (fmiOut, *) '----------------------------------------'
            end if
         end if
         if (time >= xfaims_params%sp_spwn_f) then
            ! Deactivate XFAIMS
            if (glzXFActive) then
               glzXFActive = .false.
               write (fmiOut, *)
               write (fmiOut, *) '---------------------------------------'
               write (fmiOut, *) 'Leaving field coupling region => XFAIMS'
               write (fmiOut, *) '---------------------------------------'
            end if
         end if
      end if

   end subroutine activate_xfaims

end module XFAIMSModule
