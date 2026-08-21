module PropVVRingMod
   implicit none
contains

      subroutine PropVVRing(T1,itraj,TimeStep,skip0)
      use RingModule, only: ringprop, ring_centvec, FMS_rings
      use GlobalModule, only: DefReal, DefInt
      use TrajectoryModule, only: T_Trajectory
      use ParticleModule, only: T_Particle
      use TrajectoryCalcsModule, only: FMS_UpdatePES, FMS_Forces, FMS_PotentialT

      Type (T_Trajectory) T1
      integer, intent(in) :: itraj
      real(Kind=DefReal) TimeStep
      logical,optional :: skip0
      logical skip0here
      real(8), allocatable :: r_cent(:,:)
      integer :: iparticle

      if (present(skip0)) then
        skip0here=skip0
      else
        skip0here=.false.
      endif

      !note: since ringprop uses a temporary trajectory structure,
      !the evecs of T1 are not updated
      call ringprop(FMS_rings(itraj),TimeStep,force_func,skip0here,.true.)
      ! Copy ring centroid positions back to T1%Particle
      allocate(r_cent(FMS_rings(itraj)%ndim, FMS_rings(itraj)%natom))
      r_cent = ring_centvec(FMS_rings(itraj)%r, FMS_rings(itraj)%ndim, FMS_rings(itraj)%natom, FMS_rings(itraj)%n)
      do iparticle = 1, FMS_rings(itraj)%natom
        call T1%Particle(iparticle)%set_pos(r_cent(:, iparticle))
      end do
      deallocate(r_cent)
      ! Update electronic structure at centroid geometry
      call FMS_UpdatePES(T1, zForce=.true.)
      return

      contains
      function force_func(r_vec,ndim,natom,n)
      implicit none
      type (T_Trajectory) Ttemp
      integer,intent(in) ::  ndim,natom,n
      real(8), dimension(ndim,natom,n),intent(in) :: r_vec
      real(8), dimension(ndim,natom,n) :: force_func
      real(8), dimension(ndim*natom) :: forcevec
      real(8) :: V

      integer ibead,iparticle,idm,iindex
      !Create temporary trajectory structure, into which we will stuff
      !ring bead positions to get the force

      Ttemp=T1

      do ibead=1,n
        do iparticle=1,natom
          do idm=1,ndim
            call Ttemp%Particle(IParticle)%set_pos(idm, r_vec(idm,iparticle,ibead))
          enddo
        enddo
        if (FMS_rings(itraj)%method.ne.'GAUSSAVE') then
          forcevec=FMS_Forces(Ttemp)
          iindex=1
          do iparticle=1,natom
            do idm=1,ndim
              force_func(idm,iparticle,ibead)=forcevec(iindex)
              iindex=iindex+1
            enddo
          enddo
        else
         V=FMS_PotentialT(Ttemp)
          do iparticle=1,natom
            do idm=1,ndim
              force_func(idm,iparticle,ibead) = -2.0d0 * &
                 TTemp%Particle(IParticle)%Width * &
                 (TTemp%Particle(IParticle)%get_pos(idm) - &
                  T1%Particle(IParticle)%get_pos(idm)) * V
            enddo
          enddo
        endif
      enddo
      call Ttemp%destroy()

      end function force_func


      end subroutine PropVVRing

end module PropVVRingMod
