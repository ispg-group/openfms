!  Copyright Todd J. Martinez and Raphael D. Levine, 1994
!> @brief Parameters and methods for the template-driven electronic
!!  structure interface
!!
!! The template interface allows FMS to interface with arbitrary
!! electronic structure packages through command-line interaction. The
!! user must supply a template instructing FMS how to write an input
!! file for the ESP, and additional templates instructing FMS how to
!! read the necessary quantities from the resulting output.
!<
module TemplateModule
      use GlobalModule
      implicit none
      save
      private
      public :: FMS_TemplateInit, FMS_TemplateClean
      public :: FMS_RunCommand, FMS_WriteInput, FMS_ReadOutput
!      integer (kind=DefInt) :: tpiState,tpjState,tpkState
      integer (kind=DefInt), public :: tpigradsign

      contains

!------------------------------------------------------------------------
      subroutine FMS_TemplateInit()
      write(6,*)'Template Interface disabled'
      stop 1

      end subroutine FMS_TemplateInit

      subroutine FMS_RunCommand(ixState)
      integer (kind=DefInt), intent(in) :: ixState
      write(6,*)'[FMS_RunCommand] : this subroutine is not active'
      stop 1
      end subroutine FMS_RunCommand
!!    \see FMS_ReadOutput
!<
      subroutine FMS_WriteInput(par,ixstate,cfTemplate,Orb,CIVec)
      use ElecStrucModule

      real (kind=DefReal), optional :: Orb(esnBasis,esnBasis), CIVec(:,:)

      real (kind=DefReal) :: par(:)
      character(len=*), intent(in) :: cfTemplate

      integer (kind=DefInt) :: ixstate
      write(6,*)'[FMS_WriteInput] : this subroutine is not active'
      stop 1
      end subroutine FMS_WriteInput
!<
      subroutine FMS_ReadOutput(par,ixState,cfTemplate)
      real (kind=DefReal) :: par(:)
      character(len=*), intent(in) :: cfTemplate
      integer (kind=DefInt) :: ixState

      write(6,*)'[FMS_ReadOutput] : this subroutine is not active'
      stop 1
      end subroutine FMS_ReadOutput

!     Clear old input and output files so they won't be accidentally
!     re-used
      subroutine FMS_TemplateClean()
      write(6,*)'[FMS_TemplateClean] : this subroutine is not active'
      stop 1
      end subroutine FMS_TemplateClean

end module TemplateModule
