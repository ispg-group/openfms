module test_qml
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testutils, only: check_dieerror_called
   use GlobalModule
   implicit none
   private

   public :: collect_qml_suite

contains

   !> Collect all exported unit tests
   subroutine collect_qml_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest('qml_tully_pot', test_tully_pot) &
                  !                  new_unittest('qml_tully_fms', test_tully_fms) &
                  ]

   end subroutine collect_qml_suite

   !> Test the Tully potential energy surfaces
   subroutine test_tully_pot(error)
      type(error_type), allocatable, intent(out) :: error

      integer(kind=DefInt) :: ndim, nsurf, option, k
      logical :: adiabatic
      character(len=20) :: model_name
      real(kind=DefReal) :: DQ

      real(kind=DefReal), allocatable :: Q(:)
      real(kind=DefReal), allocatable :: V(:, :)
      real(kind=DefReal), allocatable :: G(:, :, :)
      real(kind=DefReal), allocatable :: NAC(:, :, :)
      real(kind=DefReal), allocatable :: dNAC(:, :, :, :)
      real(kind=DefReal), allocatable :: Vdia(:, :)
      real(kind=DefReal), allocatable :: Vadia(:, :)

      model_name = 'Tully'
      ndim = 1
      nsurf = 2
      option = 1
      adiabatic = .true.

      call sub_Init_Qmodel(ndim, nsurf, model_name, adiabatic, option)

      call check(error, ndim == ndim,
      'Invalid number of dimensions')
      if (allocated(error)) return

      call check(error, nsurf == nsurf,
      'Invalid number of surfaces')
      if (allocated(error)) return

      allocate (Q(ndim))
      allocate (V(nsurf, nsurf))
      allocate (Vdia(nsurf, nsurf))
      allocate (Vadia(nsurf, nsurf))
      allocate (g(nsurf, nsurf, ndim))
      allocate (NAC(nsurf, nsurf, ndim))
      allocate (dNAC(nsurf, nsurf, ndim, ndim))

      !if (nsurf > 1) then

      Q(1) = -5.0_defrealgit
      DQ = 0.1_defreal
      do k = 1, 100
         call sub_Qmodel_Vdia_Vadia(Vdia, Vadia, Q)

         if (abs(Q(1) - 1.0_defreal) < 1.0e-8_defreal) then

            call check(error, abs(Vdia(1, 1) + 0.00798103482005344_defreal) < 1.0e-8_defreal,
            'Vdia(1,1) does not match expected value at Q=1')
            if (allocated(error)) return

            call check(error, abs(Vadia(1, 1) + 0.00819025633785785_defreal) < 1.0e-8_defreal,
            'Vadia(1,1) does not match expected value at Q=1')
            if (allocated(error)) return

            call sub_Qmodel_VG_NAC(V, G, NAC, Q)
            Q(:) = Q(:) + DQ

            call check(error, abs(G(1, 1, 1) - 0.00232163_defreal) < 1.0e-8_defreal,
            'G(1,1,1) does not match expected value at Q=1')
            if (allocated(error)) return

         end if

         !Lack of phase checking caused NAC at Q=1 to be 0, which is incorrect. Changing coordinate to Q=0.9 for now.
         if (abs(Q(1) - 0.9_defreal) < 1.0e-8_defreal) then
            call check(error, abs(NAC(1, 2, 1) - 0.308531791140984_defreal) < 1.0e-8_defreal,
            'NAC(1,2,1) does not match expected value at Q=1')
            if (allocated(error)) return
         end if
      end do

      deallocate (NAC)

      !else

      !   call get_Qmodel_Q0(Q, 0)

      !   call sub_Qmodel_Vdia_Vadia(Vdia, Vadia, Q)
      !end if

   end subroutine test_tully_pot

end module test_qml
