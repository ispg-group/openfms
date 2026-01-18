! Copyright Todd J. Martinez and Raphael D. Levine, 1994
!--------------------------------------------------------------------------------------------------
! PROGRAM: Quantics Interface
!> @author Cris Sanz Sanz, Graham Worth
!> @author Marin Sapunar, Ruđer Boškovrć Institute
!> @author Daniel Hollas, University of Bristol, adapted for OpenFMS
!> @date June, 2017
!
! Interface to Quantics operator library
!--------------------------------------------------------------------------------------------------
module QuanticsModule
      use GlobalModule, only: FMS_DieError
#ifdef Quantics
      use rddvrmod, only: nspfdof, spfdof, dop, ndof, long, maxdim, nddstate
      use operdef, only: alloc_operdef, hopsdim, ldddb
      use rdopermod, only: operinfo, rdoper
      use iorst, only: rstinfo
      use dirdyn, only: ndoftsh,dercpdim,ndofddpes,ndofdd,&
                  dbnrec,nactdim,natmtsh,ldbsave,&
                  lupdhes,lnactdb,lddrddb,ddtrajnum,num_gp
      use dirdyn, only: alloc_dirdyn, alloc_dddb, atnam
      use directdyn, only: getddpes
      use potevalmod, only: calcdiab, calcdiabder
      use psidef, only: qcentdim,gwpdim,zcent,vdimgp,dimgp,ndimgp,zgp,nsgp,totgp,&
                        sbaspar,rsbaspar
      use openmpmod, only: lompqc

      use dd_db, only: dddb_gp, getdbnrec, preparedb
      use dbcootrans, only: ltshtrans, lddtrans ! , tshtransb
      !use op2lib, only: subvxxdo1
      !use xvlib, only: mvtxdd1, mvxxdd1


      implicit none
      private
      public :: run_quantics

      real(dop), allocatable :: qcoo(:),qcoo1(:),xgp(:)
      real(dop), allocatable :: pesdia(:,:)
      real(dop), allocatable :: derdia(:,:,:)
      real(dop), allocatable :: derad(:,:,:)
      complex(dop), allocatable :: rotmatz(:,:)
      integer(long), allocatable :: point(:)
      real(dop), allocatable :: hops(:)
      integer :: gdof

contains

      subroutine run_quantics(step, xyz0, cstate, en, gra, nadvec)
      use rddvrmod, only: ilog, ldd, basis, ldbsmall, rpbaspar
      integer, intent(in) :: step
      real(dop), intent(in) :: xyz0(:, :)
      integer, intent(in) :: cstate
      real(dop), intent(out) :: gra(:, :)
      real(dop), intent(out) :: en(:)
      real(dop), intent(out) :: nadvec(:, :, :)

      integer :: i, n, m, f
      real(dop), allocatable :: xyz(:, :)

      real(dop) :: time
      logical(kind=4), save :: initialized=.false.

      allocate(xyz, source=xyz0)

      open(ilog,file='quantics.log',status='unknown',position='append')

      if (.not. initialized) then
         call initialize_quantics()
         initialized = .true.
      end if

      gra = 0.0_dop
      nadvec = 0.0_dop
      ! DH change
      en = 0.0_dop

      ! reform xyz -> qcoo (Quantics dynamical coordinates)
      if (ltshtrans) then
         print * , 'ltshtrans option not implemented!'
         stop 1
         !call FMS_DieError('ltshtrans option not implemented!')
         ! call subvxxdo1(xyz,tshxcoo0,ndoftsh)
         ! call mvxxdd1(tshtransb,xyz,qcoo,maxdim,ndoftsh,gdof)
      else
         f=0
         do n=1,natmtsh
            do i=1,3
               f=f+1
               qcoo(f)=xyz(i,n)
            enddo
         enddo
      endif

! need to add frozen coordinates to qcoo
! (it assumes coordinates are the centre of a GWP)
      qcoo1 = 0.0
      m=1
      do n=1,nspfdof(m)
         f=spfdof(n,m)
         qcoo1(f) = qcoo(n)
      enddo
! Add in any frozen coordinates
      do f=1,ndof
         if (basis(f) == 19) qcoo1(f) = rpbaspar(1,f)
      enddo

! Initialise local DBs. Need to be in Cartesians.
      if (ldd .and. ldbsmall) then
         if (lddtrans) then
            print *, 'lddtrans option not implemented!'
            ! call FMS_DieError('lddtrans option not implemented!')
            ! call ddq2x(qcoo1,xgp)
         else
            xgp=qcoo1
         endif

         num_gp = 1
         call dddb_gp(dbnrec,xgp,num_gp)
      endif

! Perform calculation of QC.
      time=0.0d0
      if (ldd) call getddpes(time,qcoo,1,1)

! PES matrix in adiabatic (en) and diabatic (pesdia) representations
! rotmatz is the ADT matrix (as a complex)
      call calcdiab(hops,en,pesdia,rotmatz,point,qcoo1,1)

! Matrix of gradients in adiabatic (derad) and diabatic (derdia).
      call calcdiabder(hops,derad,derdia,rotmatz,qcoo1,1)

! Convert forces to Cartesian and extract forces / nact
      call extrgra(cstate,en,gra,nadvec,derad)

      close(ilog)

      end subroutine run_quantics

!#######################################################################

      subroutine initialize_quantics()
         use rddvrmod, only: c5, allocmemory, dbmemdim, ltraj, ldd, rlaenge, olaenge, &
            & idvr, irst, ioper, ddpath, ddname, dlaenge, ldbsmall, dname, feb, macheps, &
            & oname, rname, gdim, operfile, gwpm, &
            & rdmemdim, alloc_dvrdat, alloc_grddat, dvrinfo

         real(dop), external :: dlamch
         logical(kind=4)     :: check, lerr, linwf
         character(len=c5)   :: filename, string
         integer :: ilbl, jlbl
         integer :: chkdvr, chkgrd, chkpsi, chkprp

         macheps = dlamch('P')

         string='../..'
         ilbl=5
         call abspath(string,ilbl)
         dname = string
         dlaenge = index(dname,' ')-1
         oname = string
         olaenge = index(oname,' ')-1
         rname = string
         rlaenge = index(rname,' ')-1

! turn off parallelisation of QC calcs (omp threads do separate trajs).
         lompqc=.false.

!-----------------------------------------------------------------------
! get array dimensions
!-----------------------------------------------------------------------
         inquire(irst,opened=check)
         if (check) close(irst)
         filename=rname(1:rlaenge)//'/restart'
         ilbl=index(filename,' ')-1
         open(irst,file=filename(1:ilbl),form='unformatted',status='old')
         call rdmemdim(irst)
         close(irst)

!-----------------------------------------------------------------------
! For DD calculations, DB needs to be read into memory
!-----------------------------------------------------------------------
         if (ldd) then
            ldbsave = .true.
            if (lddrddb) then
               lupdhes = .true.
            else
               lupdhes = .false.
            endif

! if not using a DB, do not allocate large memory
            if (.not. (lddrddb)) then
               dbmemdim = 1
               ldbsmall = .false.
            endif
         endif

!-----------------------------------------------------------------------
! Allocate memory
!-----------------------------------------------------------------------
         allocmemory=0
         call alloc_dvrdat()
         call alloc_grddat()
         call alloc_operdef()
         if (ldd .or. ltraj) then
            allocate(gwpdim(1,1))
            allocate(zcent(1,1))
            allocate(vdimgp(1,1))
            allocate(dimgp(1,1))
            allocate(ndimgp(1,1))
            allocate(zgp(1))
            allocate(nsgp(1))
            allocate(rsbaspar(sbaspar,maxdim,1))
            call alloc_dirdyn()
         endif

!-----------------------------------------------------------------------
! Read system / DVR information
!-----------------------------------------------------------------------
         filename=dname(1:dlaenge)//'/dvr'
         ilbl=index(filename,' ')-1
         open(idvr,file=filename,form='unformatted',status='old')
         chkdvr=1
         call dvrinfo(lerr,chkdvr)
         close(idvr)

!-----------------------------------------------------------------------
! Read data from oper file
!-----------------------------------------------------------------------
         ddpath = ' '
         filename=oname(1:olaenge)//'/oper'
         ilbl=index(filename,' ')-1
         open(ioper,file=filename,form='unformatted',status='old')
         chkdvr=1
         chkdvr=2
         chkgrd=1
         call operinfo(lerr,chkdvr,chkgrd)

!-----------------------------------------------------------------------
! read in coordinate transformation information
!-----------------------------------------------------------------------
         if (lddtrans .or. ltshtrans) then
            call FMS_DieError('lddtrans/ltshtrans option not implemented!')
            !call alloc_dbcootrans()
            !call rdddtrans(ioper)
         endif

         close(ioper)

!-----------------------------------------------------------------------
! Open DB and find out what it contains before allocating memory
!-----------------------------------------------------------------------
         if (ldddb) then
            call preparedb(1)
            call getdbnrec(dbnrec)
            call alloc_dddb()
         endif

!-----------------------------------------------------------------------
! Read data needed by the operator
!-----------------------------------------------------------------------
         operfile=oname(1:olaenge)//'/oper'
         allocate(hops(hopsdim))
         chkdvr=2
         chkgrd=1
         call rdoper(hops,chkdvr,chkgrd)

!-----------------------------------------------------------------------
! DD needs to read restart file for info on how Shepard Interpolation is
! being done
!-----------------------------------------------------------------------
         if (ldd) then
           filename=rname(1:rlaenge)//'/restart'
           ilbl=index(filename,' ')-1
           open(irst,file=filename,form='unformatted',status='old')
           chkdvr=0
           chkgrd=0
           chkpsi=0
           chkprp=1
           call rstinfo(linwf,lerr,chkdvr,chkgrd,chkpsi,chkprp)
           close(irst)
         endif

!-----------------------------------------------------------------------
! get no. of states and no. of dynamical coordinates
!-----------------------------------------------------------------------
         nddstate = gdim(feb)
         gdof = nspfdof(1)

!-----------------------------------------------------------------------
! qcentdim is needed in getddpes as dimension of Ndof (effectively 1GWP)
!-----------------------------------------------------------------------
         if (ldd) then
            gwpdim(1,1) = 1
            zcent(1,1) = 1
            gwpm(1)=.true.
            qcentdim = nspfdof(1)

!needed for getddpes (No. of configurations is 1)
            vdimgp(1,1) = 1
            dimgp(1,1) = 1
            ndimgp(1,1) = 1
            nsgp(1)=1
            totgp = 1
            zgp(1) = 1

! no. of NACTS
            if (nddstate > 1) then
               lnactdb = .true.
               nactdim = nddstate*(nddstate-1)/2
            endif

! get trajectory number from directory name
            string=ddname
            call dirpath(string,ilbl)
            jlbl=ilbl-1
            ilbl=jlbl
            do
               if (string(ilbl-1:ilbl-1) == '.') exit
               ilbl=ilbl-1
            enddo
            read(string(ilbl:jlbl),*) ddtrajnum

         endif

         ! Allocate memory
         allocate(qcoo(ndoftsh))
         allocate(qcoo1(maxdim))
         allocate(xgp(maxdim))
         allocate(pesdia(nddstate,nddstate))
         allocate(rotmatz(nddstate,nddstate))
         allocate(point(maxdim))
         allocate(derdia(nddstate,nddstate,maxdim))
         allocate(derad(nddstate,nddstate,maxdim))

      end subroutine initialize_quantics

      subroutine extrgra(sta,en,gra,nadvec,derad)
      integer(long), intent(in)  :: sta
      integer(long)              :: s,s1,f,n,m
      real(dop), dimension(ndoftsh,nddstate,nddstate), intent(out) :: nadvec
      real(dop), dimension(ndof)                                   :: qnadvec
      real(dop), dimension(ndoftsh),intent(out)                    :: gra
      real(dop), dimension(ndof)                                   :: qgra
      real(dop), dimension(nddstate,nddstate,maxdim), intent(in)   :: derad
      real(dop), dimension(nddstate), intent(in)                   :: en
      real(dop) :: ediff


! extract gradient for present state and transform to Cartesian
! removing frozen coordinates
      qgra(:) = 0.0
      m = 1  ! only 1 mode in TSH
      do n=1,nspfdof(m)
         f=spfdof(n,m)
         qgra(n) = derad(sta,sta,f)
      enddo

      gra(:) = 0.0
      if (ltshtrans) then
         call FMS_DieError('ltshtrans option not implemented!')
         !call mvtxdd1(tshtransb,qgra,gra,maxdim,nspfdof(1),ndoftsh)
      else
         do f=1,ndoftsh
            gra(f)=qgra(f)
         enddo
      endif

! extract nacts and transform to Cartesian
! set  for frozen coordinates to 0
      qnadvec(:) = 0.0_dop
      nadvec(:,:,:) = 0.0_dop
      m = 1  ! only 1 mode in TSH
      do s=1,nddstate
         do s1=s+1,nddstate
            do n=1,nspfdof(m)
               f=spfdof(n,m)
               qnadvec(n) = derad(s1,s,f)
            enddo

            if (ltshtrans) then
               call FMS_DieError('ltshtrans option not implemented!')
               !call mvtxdd1(tshtransb,qnadvec,nadvec(1,s1,s),maxdim,&
               !     nspfdof(1),ndoftsh)
            else
               do f=1,ndoftsh
                  nadvec(f,s1,s)=qnadvec(f)
               enddo
            endif

            ediff = en(s) - en(s1)
            if (abs(ediff) < 1.0d-6) then
               if (ediff < 0.0) then
                  ediff=-1.0d-6
               else
                  ediff=1.0d-6
               endif
            endif
            nadvec(:,s1,s) = nadvec(:,s1,s) / ediff

            nadvec(:,s,s1) = -nadvec(:,s1,s)
         enddo
      enddo

      end subroutine extrgra

!#######################################################################
#else

      implicit none
      private
      public :: run_quantics
      contains

      subroutine run_quantics(step, xyz0, cstate, en, gra, nadvec)
      integer, intent(in) :: step
      real(8), intent(in) :: xyz0(:, :)
      integer, intent(in) :: cstate
      real(8), intent(out) :: gra(:, :)
      real(8), intent(out) :: en(:)
      real(8), intent(out) :: nadvec(:, :, :)

      gra = 0.0D0; en = 0.0D0; nadvec = 0.0D0

      call FMS_DieError('Quantics interface not available, ' // &
                      & 'recompile with ESP=quantics in CONFIGFMS')
      end subroutine run_quantics

#endif

end module QuanticsModule
