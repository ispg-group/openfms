module XFAIMSModule
   use GlobalModule
   use TrajectoryModule
   implicit none

   private
   public :: print_xfaims_params, activate_xfaims

contains

   subroutine print_xfaims_params()
      implicit none

      character(len=*), parameter :: divider = &
                                     ' -------------------------------------------------------'

      write (fmiOut, *)
      write (fmiOut, '(a)') divider
      write (fmiOut, '(a)') ' XFAIMS parameters'
      write (fmiOut, '(a)') divider

      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Field amplitude f0_xfr:', f0_xf, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Field frequency freq_xfr:', freq_xf, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8,a)') 'Pulse width sigma_xfr:', sigma_xf, ' a.u.'
      write (fmiOut, '(1x,a,t34,es16.8)') 'Carrier-envelope phase:', CEP_xf
      write (fmiOut, '(1x,a,t34,f16.4,a)') 'Pulse center t0_xfr:', t0_xf, ' a.u.'

      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization x:', polx_xf
      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization y:', poly_xf
      write (fmiOut, '(1x,a,t34,f16.8)') 'Polarization z:', polz_xf

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Field coupling timestep:', gldCoupFieldTimeStep, ' a.u.'

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Spawning interval start:', sp_spwn_i, ' a.u.'

      ! XFAIMS is only activated at multiples of the global timestep in openfms.F90.
      if (modulo(sp_spwn_i, gldTimeStep) /= 0.0) then
         write (fmiOut, '(1x,a)') &
            'WARNING: XFAIMS will be initiated at the closest following multiple of TimeStep.'
      end if

      write (fmiOut, '(1x,a,t34,f16.4,a)') &
         'Spawning interval end:', sp_spwn_f, ' a.u.'

      ! XFAIMS is only deactivated at multiples of the global timestep in openfms.F90.
      if (modulo(sp_spwn_f, gldTimeStep) /= 0.0) then
         write (fmiOut, '(1x,a)') &
            'WARNING: XFAIMS will be terminated at the closest following multiple of TimeStep.'
      end if

      write (fmiOut, '(1x,a,t34,l16)') &
         'Only one field spawning:', onespawnonly_xf

      write (fmiOut, '(1x,a,t34,l16)') &
         'Ignore state after field:', gldIgnoreStateAferField

      write (fmiOut, '(a)') divider
      write (fmiOut, *)

      flush (fmiOut)

   end subroutine print_xfaims_params

   subroutine activate_xfaims(time)
      real(DefReal), intent(in) :: time

      ! This is the point where we really decide if we are in the field region or not and set the corresponding timestep.
      ! Hence, we set the global variable determinig the field coupling here.
      if (glzxfaims) then
         if (time >= sp_spwn_i .and. time < sp_spwn_f) then
            ! Activate XFAIMS
            if (.not. glzXFActive) then
               glzXFActive = .true.
               write (fmiOut, *)
               write (fmiOut, *) '----------------------------------------'
               write (fmiOut, *) 'Entering field coupling region => XFAIMS'
               write (fmiOut, *) '----------------------------------------'
            end if
         end if
         if (time >= sp_spwn_f) then
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
