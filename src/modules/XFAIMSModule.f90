module XFAIMSModule
   use GlobalModule
   use TrajectoryModule
   implicit none

   private
   public :: print_xfaims_params

contains

   subroutine print_xfaims_params()
      implicit none

      character(len=*), parameter :: divider = &
                                     ' -------------------------------------------------------'

      write (fmiOut, *)
      write (fmiOut, '(a)') divider
      write (fmiOut, '(a)') ' XFAIMS parameters'
      write (fmiOut, '(a)') divider

      write (fmiOut, '(1x,a,t34,es16.8)') 'Field amplitude f0_xfr:', f0_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Field frequency freq_xfr:', freq_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Pulse center t0_xfr:', t0_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Pulse width sigma_xfr:', sigma_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Carrier-envelope phase:', CEP_xf

      write (fmiOut, '(1x,a,t34,es16.8)') 'Polarization x:', polx_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Polarization y:', poly_xf
      write (fmiOut, '(1x,a,t34,es16.8)') 'Polarization z:', polz_xf

      write (fmiOut, '(1x,a,t34,es16.8)') &
         'Field coupling timestep:', gldCoupFieldTimeStep

      write (fmiOut, '(1x,a,t34,es16.8)') &
         'Spawning interval start:', sp_spwn_i

      write (fmiOut, '(1x,a,t34,es16.8)') &
         'Spawning interval end:', sp_spwn_f

      write (fmiOut, '(1x,a,t34,l1)') &
         'Only one field spawning:', onespawnonly_xf

      write (fmiOut, '(1x,a,t34,l1)') &
         'Ignore state after field:', gldIgnoreStateAferField

      write (fmiOut, '(a)') divider
      write (fmiOut, *)

      flush (fmiOut)

   end subroutine print_xfaims_params

end module XFAIMSModule
