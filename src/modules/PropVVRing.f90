module PropVVRingMod
   implicit none
contains

      subroutine PropVVRing(T1,itraj,TimeStep,skip0)
      use RingModule, only: ringprop, ring_centvec, FMS_rings
      use GlobalModule, only: DefReal
      use TrajectoryModule, only: T_Trajectory
      use TrajectoryCalcsModule, only: FMS_UpdatePES, FMS_PotentialT, &
                                 FMS_PhaseDot, FMS_GetForce

      Type (T_Trajectory) T1
      integer, intent(in) :: itraj
      real(Kind=DefReal) TimeStep
      logical,optional :: skip0
      logical skip0here
      real(8), allocatable :: r_cent(:,:), p_cent(:,:)
      integer :: iparticle
      real(kind=DefReal) :: g1_0, g1_1, g2_0, g2_1

      if (present(skip0)) then
        skip0here=skip0
      else
        skip0here=.false.
      endif

      g1_0 = FMS_PhaseDot(T1)
      g2_0 = 2.d0 * dot_product(FMS_GetForce(T1), T1%get_vel())

      !note: since ringprop uses a temporary trajectory structure,
      !the evecs of T1 are not updated
      call ringprop(FMS_rings(itraj),TimeStep,force_func,skip0here,.true.)
      ! Copy ring centroid positions back to T1%Particle
      allocate(r_cent(FMS_rings(itraj)%ndim, FMS_rings(itraj)%natom))
      r_cent = ring_centvec(FMS_rings(itraj)%r, FMS_rings(itraj)%ndim, FMS_rings(itraj)%natom, FMS_rings(itraj)%n)
      p_cent = ring_centvec(FMS_rings(itraj)%p, FMS_rings(itraj)%ndim, FMS_rings(itraj)%natom, FMS_rings(itraj)%n)
      do iparticle = 1, FMS_rings(itraj)%natom
        call T1%set_pos(iparticle, r_cent(:, iparticle))
        call T1%set_mom(iparticle, p_cent(:, iparticle))
      end do
      deallocate(r_cent, p_cent)
      ! Update electronic structure at centroid geometry
      call FMS_UpdatePES(T1, zForce=.true.)

      g1_1 = FMS_PhaseDot(T1)
      g2_1 = 2.d0 * dot_product(FMS_GetForce(T1), T1%get_vel())

      T1%Phase = T1%Phase &
                + 0.5d0 * (g1_0 + g1_1) * TimeStep &
                - 0.125d0 * (g2_0 - g2_1) * TimeStep**2

      call T1%rescale_phases()
      return

      contains
      function force_func(r_vec,ndim,natom,n)
      implicit none
      type (T_Trajectory) Ttemp
      integer,intent(in) ::  ndim,natom,n
      real(8), dimension(ndim,natom,n),intent(in) :: r_vec
      real(8), dimension(ndim,natom,n) :: force_func
      real(8) :: V

      integer ibead,iparticle,idm
      !Create temporary trajectory structure, into which we will stuff
      !ring bead positions to get the force

      Ttemp=T1

      do ibead=1,n
        if (FMS_rings(itraj)%method.ne.'GAUSSAVE') then
          call FMS_EvaluateRingBead(Ttemp,r_vec(:,:,ibead),T1%StateID, &
             force_func(:,:,ibead),FMS_rings(itraj)%V_ext(ibead))
        else
         do iparticle=1,natom
           do idm=1,ndim
             call Ttemp%Particle(iparticle)%set_pos(idm,r_vec(idm,iparticle,ibead))
           enddo
         enddo
         V=FMS_PotentialT(Ttemp)
          FMS_rings(itraj)%V_ext(ibead)=V
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

      subroutine FMS_EvaluateRingBead(Ttemp,r_bead,iState,force_bead,potential)
      use TrajectoryModule, only: T_Trajectory
      use TrajectoryCalcsModule, only: FMS_Forces, FMS_PotentialT
      implicit none
      type(T_Trajectory), intent(inout) :: Ttemp
      real(8), intent(in) :: r_bead(:,:)
      integer, intent(in) :: iState
      real(8), intent(out) :: force_bead(:,:)
      real(8), intent(out) :: potential
      integer :: iparticle

      do iparticle = 1, size(r_bead, 2)
        call Ttemp%set_pos(iparticle, r_bead(:, iparticle))
      end do

      force_bead=reshape(FMS_Forces(Ttemp,iState),shape(force_bead))
      potential=FMS_PotentialT(Ttemp,iState)

      end subroutine FMS_EvaluateRingBead

end module PropVVRingMod
