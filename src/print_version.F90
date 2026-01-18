! Copyright Todd J. Martinez and Raphael D. Levine, 1994
subroutine print_version(unit)
   implicit none
   integer, intent(in) :: unit
   character(len=*), parameter :: FMS_VERSION = "v0.1-alpha"
   ! GIT_VER is defined in Makefile
   character(len=*), parameter :: GIT_COMMIT=GIT_VER

   write(unit, *) 'version: '//FMS_VERSION
   write(unit, *) 'git commit: '//trim(GIT_COMMIT)
end subroutine print_version
