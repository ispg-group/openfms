! Copyright Todd J. Martinez and Raphael D. Levine, 1994
module RandomModule
   use GlobalModule, only: DefInt, DefReal, fmiOut
   implicit none
   private
   public :: FMS_ranb, initialize_fortran_prng

contains
!>
!!    @brief Random number generator
!! \param IXT if <> 0, generator is reseeded.  To get further numbers in
!!     a pseudo-random sequence, ixt should be set to zero.
!!    @ingroup numerics
!<
   function fms_ranb(ixt)

      real*8 :: FMS_RANB
      integer*4, intent(in) :: ixt
      integer*4, save :: FMS_IX
      integer*4 :: I1, I2, I3, I4, I5, I6, I7, I8, I9
      data I2/16807/, I4/32768/, I5/65536/, I3/2147483647/

      if (ixt /= 0) then
         FMS_ix = ixt
         write (fmiOut, '(a, i0)') 'RanB seeded with ', ixt
      end if

      I6 = FMS_IX / I5
      I7 = (FMS_IX - I6 * I5) * I2
      I8 = I7 / I5
      I9 = I6 * I2 + I8
      I1 = I9 / I4
      FMS_IX = (((I7 - I8 * I5) - I3) + (I9 - I1 * I4) * I5) + I1
      if (FMS_IX < 0.0d0) FMS_IX = FMS_IX + I3
      FMS_RANB = (4.656613d-10) * dble(FMS_IX)

   end function fms_ranb

   ! This code was taken from ABIN
   ! Move this into a module?
   ! Initializing PRNG subroutine random_number() defined by Fortran standard
   ! https://stackoverflow.com/questions/51893720/correctly-setting-random-seeds-for-repeatability
   subroutine initialize_fortran_prng(user_seed)
      use, intrinsic :: iso_fortran_env, only: int64
      integer, intent(in) :: user_seed
      integer, allocatable :: seeds(:)
      integer :: i, seed_size
      integer(int64) :: s
      real(defReal) :: drans(100)

      call random_seed(size=seed_size)
      allocate (seeds(seed_size))

      ! We're use a simple PRNG defined below for the initial seed state
      s = int(user_seed, kind(s))
      do i = 1, seed_size
         s = lcg(s)
         seeds(i) = int(s)
      end do

      call random_seed(put=seeds)
      ! Prime the prng by discarding first 100 values
      call random_number(drans)
   end subroutine initialize_fortran_prng

   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
   ! Taken from:
   ! https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
   integer function lcg(s)
      use, intrinsic :: iso_fortran_env, only: int64
      ! TODO: I think this should be intent(in)?
      integer(int64), intent(inout) :: s

      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
   end function lcg

end module RandomModule
