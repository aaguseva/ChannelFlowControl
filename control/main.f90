!c/********************************************************************/
!c/*                                                                  */
!c/*    This code solves incompressible Navier-Stokes equations       */
!c/*    for a three-dimensional channel flow for vorticity and        */
!c/*    laplacian of wall-normal velocity, with flow control at       */
!c/*    the wall through adjusting boundary conditions                */
!c/*                                                                  */
!c/*    Contributors: A. Guseva, O. Flores and J. Jiménez             */ 
!c/*                                                                  */
!c/*     REFERENCES                                                   */
!c/*                                                                  */   
!c/*       - Guseva, A. and Jiménez, J. (2022) Linear instability     */  
!c/*         and resonance effects in large-scale opposition          */    
!c/*         flow control. J. Fluid Mech., 935, A35.                  */    
!c/*                                                                  */
!c/*       - Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statis    */
!c/*         tics in fully developped channel flow at low Reynolds    */
!c/*         numbers, J. Fluid Mech., 177, 133-166                    */
!c/*                                                                  */
!c/********************************************************************/
      program chanddt

      use statistics
      use ctes
      use mpi_var

      implicit none
      include "mpif.h"
      include "ctes3D"

      integer nbuffsize,nwkasize,j
      real*4, allocatable:: vor(:),phi(:)
      real*4 u00(my),w00(my)
      real*4, allocatable::  wk(:)
      integer ihv,ihg,iphiwk,ivorwk,irf0u,irf0w,iu00wk,iw00wk,idvordy,ichwk,nstr
      integer nspsize
      real*4, allocatable:: sp(:)
      real(4), dimension(:), allocatable ::chwk
      real(4),dimension(:),allocatable:: u,v,w,dudy,dvdy,dwdy
      real(4),dimension(:),allocatable:: p_rhs, press
      real(4),dimension(:),allocatable:: tauxy,tauyz,tauxz


!c                              /*   initializes everything    */
      call mpi_initialize()
!c--------------- initializes commons and things
      call initcr(myid)
!c--------------  allocates spectra
      jspbb = jspb
      if (jspe.lt.jspb) jspbb=jspe
      nspsize = (mx1+1)*(nz1+1)*8*(jspe-jspbb+1)
      allocate(sp(nspsize))
!c--------------- allocates buffers
      nbuffsize = mx*max(mmy*mz,mmz*my)
      nwkasize  = 6*nbuffsize + 4*my
!c ----------  storage for the pressure

      allocate(vor(nbuffsize))
      allocate(phi(nbuffsize))
      allocate(wk(nwkasize))


      allocate(u(nbuffsize),v(nbuffsize),w(nbuffsize))
      allocate(dudy(nbuffsize),dvdy(nbuffsize),dwdy(nbuffsize))
      allocate(tauxy(nbuffsize),tauyz(nbuffsize),tauxz(nbuffsize))
      allocate(p_rhs(nbuffsize),press(nbuffsize))


      allocate(chwk(nbuffsize))
      nacumsp = 0
      sp(1:nspsize)=0.
!c               /*   read input data           */
      if (myid==0) write(*,*) 'calling getfil'
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call getfil(vor,phi,u00,w00,wk,myid)
!c      write(*,*) myid,' MEM: ',(nspsize+2*nbuffsize+nwkasize)*4/1024**2,
!c     .                ' MB + ',(3*mx*mz+6*(mgalx+2)*mgalz)*4/1024**2,
!c     .                ' MB '
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c       Cuando recorto el tama�o de la caja por un factor gamma,
!c       tambien aplico este factor a vor y phi:
!c       Lnew = gamma*Lold -> vor_new = vor_old*gamma
!c                         -> phi_new = phi_old*gamma^2
!c      do j=1,nbuffsize
!c         vor(j) = vor(j)*2.
!c         phi(j) = phi(j)*4.
!c      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                    /* start time advancement  */
      irf0u   = 1
      irf0w   = irf0u   + my
      iu00wk  = irf0w   + my
      iw00wk  = iu00wk  + my
      ihv     = iw00wk  + my
      ihg     = ihv     + nbuffsize
      iphiwk  = ihg     + nbuffsize
      ivorwk  = iphiwk  + nbuffsize
      idvordy = ivorwk  + nbuffsize
      ichwk   = idvordy + nbuffsize


      call cross1(vor,phi,u00,w00,wk(irf0u),wk(irf0w),wk(iu00wk),wk(iw00wk), &
                  wk(ihv),wk(ihg),wk(iphiwk),wk(idvordy),wk(ivorwk), &
                  wk(idvordy),wk(ichwk),sp,myid, &
                  u,v,w,dudy,dvdy,dwdy,tauxy,tauyz,tauxz,press,p_rhs)

!c                    /* finalize procedure      */
       call mpi_fin()

end program chanddt

!c ====================================================================
!c
!c                   get data field from file
!c                                             jjs  27/1/01
!c ====================================================================

      subroutine getfil(vor,phi,u00,w00,wk1,myid)

      use running
      use ctes

      implicit none
      include "ctes3D"
      include "mpif.h"

      integer istat(MPI_STATUS_SIZE),ierr
      integer myid,master,i,pe,k
      real*4  a0e
      real*4 vor(mx,mz,*),phi(mx,mz,*),u00(*),w00(*),wk1(*)
      real*4  Ree,alpe,bete,ce
      integer j,mxe,mye,mze,ntotr,iproc,mym,mmy2
      real*4 dumvar(1:20) ! dummy variable to match code input - output

      master = 0

!c      ---------------  zero everything, fitf
         vor(1:mx*mz*mmy,1,1) = 0.
         phi(1:mx*mz*mmy,1,1) = 0.

         u00(1:my) = 0.
         w00(1:my) = 0.

!c      ---------------  read from file and distribute

      if (myid.eq.master) then     ! ----- this is done by the master
         open(iinp,file=filinp,status='unknown',form='unformatted',access = 'stream')
         read(iinp) time,Ree,alpe,bete,a0e,mxe,mye,mze,pe,ce
         write(*,*)
         write(*,*) 'reading input file ...'
         write(*,*) 'time=',time,mxe,mye,mze,'mesh:',ce,pe
         write(*,*)  Ree, alpe, bete, a0e

         read(iinp) (dumvar(j),j=1,20)  ! reading boundaries into a dummy variable
         ntotr=2*mxe*mze
!c     ------------- 00 modes ------------------
         write(*,*) 'reading 0 modes'
         read(iinp) (wk1(j),j=1,2*mye)
         mym = min(my,mye)
         do j=1,mym
            u00(j) = wk1(2*j-1)
            w00(j) = wk1(2*j)
         enddo


         do iproc=1,numerop-1
            call MPI_SEND(mxe,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mye,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mze,1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(u00,my,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(w00,my,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,ierr)
         enddo
!c         ------------  other modes,  master node -----
         mmy2 = min(je,mye)-jb+1
         write(*,*) 'master reads its data'
         do j=1,mmy2
            read(iinp) (wk1(i),i=1,ntotr)
            call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
         enddo
!c         ------------  other modes,  distribute to slaves --
         do iproc=1,numerop-1
            mmy2 = min(jend(iproc),mye)-jbeg(iproc)+1
            do j=1,mmy2
               read(iinp) (wk1(i),i=1,ntotr)
               call MPI_SEND(wk1,ntotr,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,ierr)
            enddo
         enddo
         close(iinp)
      else          ! -------- this is done by all slaves
!c         ------------  receive 00 mode ----------
         call MPI_RECV(mxe,1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mye,1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mze,1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(u00,my,MPI_REAL,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(w00,my,MPI_REAL,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         ntotr=2*mxe*mze
!c         ------------  receive other modes ----------
         mmy2 = min(je,mye)-jb+1
         do j=1,mmy2
            call MPI_RECV(wk1,ntotr,MPI_REAL,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
         enddo
      endif
!
!c ---- temporal patch: zeroing (0,beta) modes of phi
      if (ifor.eq.1) then
         phi(1:2,1:mz,1:min(je,mye)-jb+1) = 0.
      endif

!c ----- before going any further exchange cuts in y with cuts in z!
      call chikj2jik(vor,vor,wk1,wk1,myid)
      call chikj2jik(phi,phi,wk1,wk1,myid)
    end subroutine getfil

!c/********************************************************************/
      subroutine initcr(myid)

      use statistics
      use running
      use ctes
      use boundary


      implicit none
      include "mpif.h"
      include "ctes3D"

      integer myid,nn1
      parameter (nn1=40000)        !!! buffer for chebishev

      integer istat(MPI_STATUS_SIZE),ierr
      real*4 pi
      integer mybu,mzbu,i,j,iproc,jj,k
      real*8 d11(my,5),d12(my,5),d21(my,5),d22(my,5)

      NAMELIST /parameters_setup/ Re, alp, bet, a0, imesh, gamma
      NAMELIST /parameters_timestep/ Delt, var_dt, CFL, nstep, nimag, nhist, &
                                     ntimes, istart
      NAMELIST /parameters_control/ alp_u, alp_v, alp_w, nplane_c, t_lag, x_lag, f_wU, &
                                    f_lse, savefields, dtc, mxmin, mzmin, mxmax, mzmax
      NAMELIST /parameters_file/ filout, filinp, filstt, id22
      NAMELIST /parameters_var/ ifor, nstart


       Pi=(4.*atan(1.))

       if (myid.eq.0) then
         read(*,nml=parameters_setup)
         read(*,nml=parameters_timestep)
         read(*,nml=parameters_control)
         read(*,nml=parameters_file)
         read(*,nml=parameters_var)
       endif


      call MPI_BCAST(Re,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(alp,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(bet,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(a0,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(Delt,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(var_dt,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(CFL,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gamma,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)


      call MPI_BCAST(alp_u,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(alp_v,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(alp_w,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nplane_c,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(t_lag,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(x_lag,1, MPI_REAL4,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_wU,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_lse,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(savefields,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dtc,1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(mxmin,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(mzmin,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(mxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(mzmax,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

      call MPI_BCAST(imesh,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nimag,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nhist,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ifor,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(id22,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nstart,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ntimes,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)


      iinp=32
      ispf=35
      istati=0
      isn=33
      iout=31

!c ---------------  compute y coordinates, pointers and
!c ---------------  modes values for FOURIER
      call malla(my,imesh,gamma,y,fmap,y2)
      call genexp
      call pointers(jbeg,jend,kbeg,kend)
      jb=jbeg(myid)
      je=jend(myid)
      kb=kbeg(myid)
      ke=kend(myid)

      do j = 0,numerop-1
         if ((nplane_c.ge.jbeg(j)).and.(nplane_c.le.jend(j))) then
           id_up = j
         else if (((my-(nplane_c - 1)).ge.jbeg(j)).and.((my-(nplane_c - 1)).le.jend(j))) then
           id_down = j
         endif
      enddo


      mmz = ke-kb+1
      mmy = je-jb+1

!c    ------------  initializes fast fourier transforms and CFDiff ----
      call cfti(mgalz)
      call rfti(mgalx)
      call prederiv1m(my,d11(1,1),d11(1,2),d11(1,3),d11(1,4),d11(1,5),&
                      d12(1,1),d12(1,2),d12(1,3),d12(1,4),d12(1,5))
      call prederiv2 (my,d21(1,1),d21(1,2),d21(1,3),d21(1,4),d21(1,5),&
                     d22(1,1),d22(1,2),d22(1,3),d22(1,4),d22(1,5),y2(1))

      do j=1,my
         do i=1,5
           prem1(i,j)= d11(j,i)
           dt12(i,j) = d12(j,i)
           dt21(i,j) = d21(j,i)
           dt22(i,j) = d22(j,i)
         enddo
      enddo

      call bandec5(prem1,my)

      trp(0) = 0.25*(y(2)-y(1))
      trp2(0) = 25d-2*(y2(2)-y2(1))
      do j=2,my-1
         trp(j-1) = 0.25*(y(j+1)-y(j-1))
         trp2(j-1) = 25d-2*(y2(j+1)-y2(j-1))
      enddo
      trp(my-1) = 0.25*(y(my)-y(my-1))
      trp2(my-1)= 25d-2*(y2(my)-y2(my-1))

      do j=2,my-1
         hy(j) = (y(j+1)-y(j-1))/(2.*2.5)
      enddo
      hy(1)  = (y(2)-y(1))/2.5
      hy(my) = (y(my)-y(my-1))/2.5

!c --------------  prepare spectra -----
      do j=1,nspec
         jsptot(j)           =jspecy(j)
         jsptot(2*nspec+2-j) = my-jspecy(j)+1
      enddo
      jsptot(nspec+1) = (my+1)/2

      do i=0,numerop-1
        jspbeg(i)=2*nspec+2
        !DEC$ NOVECTOR
         do j=2*nspec+1,1,-1
            if (jsptot(j).ge.jbeg(i)) jspbeg(i) = j
         enddo

         jspend(i)=0
         do j=1,2*nspec+1
            if (jsptot(j).le.jend(i)) jspend(i) = j
         enddo
      enddo

      do i=0,numerop-1
         do j=1,2*nspec+1
            do jj=jspbeg(i),jspend(i)
               if (j.eq.jj) jspiproc(j) = i
            enddo
         enddo
      enddo

      do j=1,my
         jsp(j) = 0
         do jj = 1,2*nspec+1
            if (jsptot(jj).eq.j) jsp(j) = 1
         enddo
      enddo

      jspb = jspbeg(myid)
      jspe = jspend(myid)

!c ------------------ Re/dt --------------------------------
      dtr = Re/Delt
!c ------------------ MPI Datatypes ------------------------
      do iproc=0,numerop-1
         if (iproc.ne.myid) then
             mzbu = kend(iproc)-kbeg(iproc)+1
             call MPI_TYPE_VECTOR(mmy,mx*mzbu,mx*mz,MPI_REAL,myslice(iproc),ierr)
             call MPI_TYPE_COMMIT(myslice(iproc),ierr)
         endif
      enddo
!      write(*,*) 'MPI_datatypes'
!c --------------  initialize stats -------------
    call stats_init()
!c --------------  write header for output -------------
    call io_write_header()

    end subroutine initcr
!c ==============================================================
      subroutine assign(work,vor,phi,mx,mz,mxe,mze)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          single jjs 4/01/01, rewritten jjs 28/01/01
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer   mxm,mzm,klen,kini1,kini2,i,k,k1,k2,mx,mz,mxe,mze
      real*4    vor(mx,mz),phi(mx,mz)
      real*4    work(2,mxe,mze)

      mxm=min(mx,mxe)
      mzm=min(mz,mze)

      klen=(mzm+1)/2
      kini1=mze - klen + 1
      kini2=mz - klen + 1

      do k=1,klen
         do i=1,mxm
            vor(i,k)=work(1,i,k)
            phi(i,k)=work(2,i,k)
         enddo
      enddo

      do k=1,klen-1
        k1 = k + kini1
        k2 = k + kini2
        do i=1,mxm
          vor(i,k2)=work(1,i,k1)
          phi(i,k2)=work(2,i,k1)
        enddo
      enddo

    end subroutine assign



  subroutine makepoiselle(u00,w00)
  use ctes
  implicit none
  include "ctes3D"
  real*4 u00(my),w00(my)
  integer j
  do j=1,my
     u00(j) = 1.-y(j)**2
     w00(j) = 0.
  enddo

  end subroutine makepoiselle

 subroutine mpi_initialize()
   use mpi_var
   implicit none
   include 'mpif.h'
   include "ctes3D"

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  if (numprocs.ne.numerop) then
     write(*,*) 'wrong number of processors',numprocs
     write(*,*) 'compiled with', numerop
     stop
  endif

  if (myid==0) write(*,*) 'MPI enviroment initialized ..'

  end subroutine mpi_initialize

  subroutine mpi_fin()
  use mpi_var
  implicit none
  include 'mpif.h'

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  write (*,*) 'process',myid,'over'
  call MPI_FINALIZE(ierr)

  end subroutine mpi_fin
