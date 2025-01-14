!! hvhg - calculate nonlinear term in physical space

subroutine cross1(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,&
                 hv,hg,phiwk,spwk,vorwk,dvordy,chwk,sp,myid, &
                 u,v,w,dudy,dvdy,dwdy,tauxy,tauyz,tauxz, &
                 press,p_rhs)

use statistics
use running
use ctes
use planes
use boundary

implicit none
include "mpif.h"
include "ctes3D"


integer myid,iproc,leng,leng1,leng2,istep
integer irun,rkstep,i,ii,k,kk,j,jj,i1,k1,k2,ipo1,ipo2,ipo3,ipo4,ipo
real*4  r1, dtr1, du0,duL,bcbr,bctr,tmp
real*8 rk,rk2
complex*8 bcb,bct,bcbdv,bctdv,bcbo,bcto
real*4  reynota,massu0,massw0,H00u,H00w,massu,massw
real*4  massu1,massu2,massw1,massw2,cteu,ctew

integer istat(MPI_STATUS_SIZE),ierr

character*200 fname

real*4, dimension(0:2*my-1,0:mx1,kb:ke) :: hg,vorwk,hv,phiwk,&
                                                 phi,vor,dvordy
real*4 chwk(*),work(20*my)
real*4 u00(0:*),w00(0:*),rf0u(0:*),u00wk(0:*),rf0w(0:*),w00wk(0:*)

real*4 dkx2,dkz2
complex*8 dkx,dkz,dkk

real*4 sp(0:mx1,0:nz1,8,jspbb:jspe),spwk(0:mx1,0:nz1,8,1:*)
character*4 ext,ext1

real*8 iter_time,write_time,laps_time,comm_time


complex(4), dimension(0:my1,0:mx1,kb:ke) :: u,v,w,dudy,dvdy,dwdy
complex(4), dimension(0:my1,0:mx1,kb:ke) :: tauxy,tauyz,tauxz
complex(4), dimension(0:my1,0:mx1,kb:ke) :: press, p_rhs
complex(4), dimension(2,0:mx1,kb:ke) :: bcpress

complex(4),dimension(:),allocatable:: g0, g1, dg0dy, dg1dy
complex(4) fbot, ftop
integer nbuffsize
character*200 flnm

real*8, dimension(0:5,0:6) :: minphi ! first kz, than kx

comm_time = 0D0
commtimer=0.0D0
transtimer=0.0D0
totaltimer=0.0D0

ihist    = 0
istati   = 0

irun   = 0   ! first time step is special in tim3rkp
icfl   = 1   ! first time step always needs a step size

if(myid.eq.0) then
   write(ext1,'(i4.4)') id22
   fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
	    write (*,*) fname
   open(39,file=fname,status='unknown')
endif
!c========================================================================
!c                 THIS IS THE TIME LOOP
!c========================================================================

do istep=1,nstep

  call io_calc_timers_io(istep)
!c    =========== write image to a file  =================  */
    if (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then
        call io_calc_wtime(write_time)
!c           /* procedure to receive all the slices and stats  */
        call stats_send(myid)
        call io_write_state(vor,phi,u00,w00,chwk,dvordy)
        call io_write_stat()
        call io_calc_spec(sp,spwk)
        call io_write_cf(write_time)

    endif

!c     ==================================  finished writing image */

!c                                   /*     time stepper             */

!c/********************************************************************/
!c/*      da un paso en el tiempo. Runge - Kutta  para terminos       */
!c/*      convectivos euler inverso para la parte viscosa.            */
!c/*                                                                  */
!c/*       Resuelve:    Gt  + Hg = 1/Re G"                            */
!c/*                    V"t + Hv = 1/Re V""                           */
!c/*                                                                  */
!c/*   input:                                                         */
!c/*     vor: vorticidad segun la direccion y (n) 		     */
!c/*     phi: laplaciana velocidad segun la direccion y  (n)          */
!c/*     vorwk: copia de vor para el calc. de los term. no lineales   */
!c/*     phiwk: copia de phi para el calc. de los term. no lineales   */
!c/*     hg: Hg                                                       */
!c/*     hv: Hv                                                       */
!c/*     dvordy: area de trabajo de dimension mx*nxymax               */
!c/*     chwk: area de trabajo para los chz2y                         */
!c/*                                                                  */
!c/*  output:                                                         */
!c/*     vor: vorticidad segun la direccion y (n+1)                   */
!c/*     phi: laplaciana velocidad segun la direccion y  (n+1)        */
!c/*      hg: contiene  v (velocidad segun y)                         */
!c/*      hv: contiene  dv/dy                                         */
!c!/*..................................................................*/
!c/*  MODULE FOR MPI SP2                                              */
!c/*..................................................................*/
!c!/*                                                                  */
!c/*   updated jjs 07/01/01                                           */
!c/*   in jik form jca						     */
!c/*   to CFdiff by of						     */
!c*/                                                                  */
!c*/                                                                  */
!c/********************************************************************/

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   this is done only for the first step
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(irun.eq.0) then

         irun = 1
!c	/*   Reinicia el tiempo del codigo
         if (istart.eq.1) then
             time = 0.
         endif
!c	/*   Todos necesitan el tiempo!!!!
         call  MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

         nbuffsize = mx*max(mmy*mz,mmz*my)
         allocate(g0(nbuffsize),g1(nbuffsize),dg0dy(nbuffsize),dg1dy(nbuffsize))
         call var_precompute(g0, g1, dg0dy, dg1dy)

!c!        /*   calcula la v a partir de phi */
         do k=kb,ke
            k1 = icx(k-1)
!c            k1 = k-1
            do i=0,mx1
               rk = bet2(k1)+alp2(i)
               bcb = 0.0!vwall(i,k-1)
               bct = -bcb
               call Lapvdv(phi(0,i,k),hg(0,i,k),hv(0,i,k),rk,bcb,bct)
            enddo
         enddo
!c	/*     prepara phiwk,vorwk,u00wk,w00wk */

        u00wk(0:my1)=u00(0:my1)
        w00wk(0:my1)=w00(0:my1)

        call deryr(u00wk,rf0u,1,my)
        call deryr(w00wk,rf0w,1,my)

        vorwk(0:2*my-1,0:mx1,kb:ke)=vor(0:2*my-1,0:mx1,kb:ke)
        phiwk(0:2*my-1,0:mx1,kb:ke)=phi(0:2*my-1,0:mx1,kb:ke)

!c				computes mass
        massu0=0.
        massw0=0.
        do j=0,my1
           massu0 = massu0 + trp(j)*u00(j)
        enddo
        if(myid.eq.0) write(*,*) 'readed massu0=',massu0
        massu0 = .8987636566162
        if(myid.eq.0) then
          write(*,*) 'la masa esta puesta a capon!!!!!!!!!!!!!!!!!!!'
          write (*,*) 'mass:',massu0,massw0
        endif

      endif
!cccccccccc    end special first step ccccccccccccccccccccc


!c                                    /*  Runge-Kutta third order  */
!c-----!

      do rkstep=1,3

        call var_bc_press(phiwk, bcpress,rkstep)
!c-----------------------  computes d (ome_2) / dy
!c-----------------------  dvordy      : d (vorwk) / dy -- F-F
         call deryr2(vorwk,dvordy,(mx1+1)*mmz,my,chwk)
!c c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
!c all arrays are  sent from y-x --> x-z
!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
         call chjik2ikj(phiwk,phiwk,chwk,chwk,myid)
         call chjik2ikj(hv,hv,chwk,chwk,myid)
         call chjik2ikj(hg,hg,chwk,chwk,myid)
         call chjik2ikj(dvordy,dvordy,chwk,chwk,myid)
         call chjik2ikj(vorwk,vorwk,chwk,chwk,myid)
!c----------------------- storage for 0 modes
        ipo1 = 1    + my
        ipo2 = ipo1 + my
        ipo3 = ipo2 + my
        ipo4 = ipo3 + my
!c-----------------------  all nodes compute 0  mode of vorticity
!c  work(1):du00;  work(ipo1):dw00;  work(ipo2):u00;  work(ipo3):w00;
        do j=0,my1
           work(1+j)   =rf0u(j)
           work(ipo1+j)=rf0w(j)
           work(ipo2+j)=u00wk(j)
           work(ipo3+j)=w00wk(j)
        enddo

        call hvhg(phiwk,vorwk,hv,hg,rf0u,rf0w,dvordy,work,sp,myid,rkstep,&
                  u1r,u2r,u3r,o1r,o2r,o3r,u1r,u2r,u3r,o1r,o2r,o3r,&
                  u,v,w,dudy,dvdy,dwdy,tauxy,tauyz,tauxz,press,p_rhs,bcpress,g0,g1,dg0dy,dg1dy,istep)

        H00u=0.
        H00w=0.
        do j=0,my1
           H00u = trp(j)*rf0u(j)
           H00w = trp(j)*rf0w(j)
        enddo

!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
!c  at this point: dvordy, rhv and rhg are the outs
!c  they must be trasformed from xy-xz before completion
!c  dvordy: dH1/dx + dH3/dz, hg: -dH3/dx + dH1/dz, hv: d^2H2/dx^2 + d^2H2/dz^2
!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c

         call chikj2jik(dvordy,dvordy,chwk,chwk,myid)
         call chikj2jik(hv,hv,chwk,chwk,myid)
         call chikj2jik(hg,hg,chwk,chwk,myid)!!

!c------------------------ computes dvordy = d (dH1/dx + dH3/dz) / dy
         call deryr2(dvordy,dvordy,(mx1+1)*mmz,my,chwk)
!c--- computes  rhv =d^2H2/dx^2 + d^2H2/dz^2 - d(dH1/dx + dH3/dz)/dy
         hv(0:2*my-1,0:mx1,kb:ke) = hv(0:2*my-1,0:mx1,kb:ke) - dvordy(0:2*my-1,0:mx1,kb:ke)
!c-----------
!c----------- advances in time : nonlinear terms explicitly

         r1=c(rkstep)*Deltat

         do j=0,my1
            u00wk(j)=u00(j)+r1*rf0u(j)
            w00wk(j)=w00(j)+r1*rf0w(j)
         enddo


         vorwk(0:2*my-1,0:mx1,kb:ke)=vor(0:2*my-1,0:mx1,kb:ke)+r1*hg(0:2*my-1,0:mx1,kb:ke)
         phiwk(0:2*my-1,0:mx1,kb:ke)=phi(0:2*my-1,0:mx1,kb:ke)+r1*hv(0:2*my-1,0:mx1,kb:ke)

         dtr1=dtr/c(rkstep)

        rf0u(0:my1)=-dtr1*u00wk(0:my1)
        rf0w(0:my1)=-dtr1*w00wk(0:my1)

        hg(0:2*my-1,0:mx1,kb:ke)=-dtr1*vorwk(0:2*my-1,0:mx1,kb:ke)
        hv(0:2*my-1,0:mx1,kb:ke)=-dtr1*phiwk(0:2*my-1,0:mx1,kb:ke)


!!c------------  updating boundary condition for the viscous time step
         if (f_wU==1) t_lag = r1 ! in phase with U(nplane_c) if the flag is set to one

         do k=kb,ke
            k1 = icx(k-1)
            do i=0,mx1
               rk = bet2(k1)+alp2(i)
               rk2 = dtr1+rk

               if ((mod(istep,dtc)==0).and.(i.ge.mxmin).and.(i.le.mxmax).and.(k1.ge.mzmin).and.(k1.le.mzmax).and.((i+k1).ne.0)) then
                   fbot = exp(-xalp(i)*((u00(nplane_c-1)-a0)*t_lag + x_lag)) 
                   ftop = exp(-xalp(i)*((u00(my1-(nplane_c-1))-a0)*t_lag + x_lag)) 
                   bcbo = fbot*(- xalp(i)*alp_w*w_bc(1,i,k-1) + xbet(k-1)*alp_u*u_bc(1,i,k-1))
                   bcto = ftop*(- xalp(i)*alp_w*w_bc(2,i,k-1) + xbet(k-1)*alp_u*u_bc(2,i,k-1))
                   bcb = fbot*(alp_v * v_bc(1,i,k-1))
                   bct = ftop*(alp_v * v_bc(2,i,k-1))
                   bcbdv = fbot*(- xalp(i)*alp_u*u_bc(1,i,k-1) - xbet(k-1)*alp_w*w_bc(1,i,k-1))
                   bctdv = ftop*(- xalp(i)*alp_u*u_bc(2,i,k-1) - xbet(k-1)*alp_w*w_bc(2,i,k-1))
               else
                 bcbo = 0.
                 bcto = 0.
                 bcb = 0.
                 bct = 0.
                 bcbdv = 0.
                 bctdv = 0.
               endif
               call lapsov(phiwk(0,i,k),hg(0,i,k),hv(0,i,k),hv(0,i,k),&
                          vorwk(0,i,k),hg(0,i,k),&
                          rk2,rk,bcb,bct,bcbdv,bctdv,bcbo,bcto)


            enddo
         enddo


!!c                        /* Kx = Kz = 0 modes         */
!!c			Se calculan en dos fases:


!!c	Velocidad sin correccion de presion-masa
        bcbr = 0.   !! Estas no tienen por que ser nulas!!
        bctr = 0.
        rk = dtr1
 
        call Lapv1(rf0u,u00wk,rk,bcbr,bctr)
        call Lapv1(rf0w,w00wk,rk,bcbr,bctr)

!!c	Velocidad correccion presion-massa, cc homogeneas:
        bcbr = 0.
        bctr = 0.
        do j=0,my1
           rf0u(j) = -dtr1
           rf0w(j) = -dtr1
        enddo

        rk = dtr1
        call Lapv1(rf0u,rf0u,rk,bcbr,bctr)
        call Lapv1(rf0w,rf0w,rk,bcbr,bctr)
 


!!c!!	Calculo masa de cada flujo:

        massu1=0.
        massu2=0.
        massw1=0.
        massw2=0.
        do j=0,my1
           massu1 = massu1 + trp(j)*u00wk(j)
           massw1 = massw1 + trp(j)*w00wk(j)
           massu2 = massu2 + trp(j)*rf0u(j)
           massw2 = massw2 + trp(j)*rf0w(j)
        enddo

        cteu =(massu0 - massu1)/massu2
        ctew =(massw0 - massw1)/massw2

        do j=0,my1
           u00wk(j) = u00wk(j) + cteu*rf0u(j)
           w00wk(j) = w00wk(j) + ctew*rf0w(j)
        enddo

        call deryr(u00wk,rf0u,1,my)
        call deryr(w00wk,rf0w,1,my)


end do             !!! finalizado subpaso del RK3

      do j=0,my1
         u00(j)=u00wk(j)
         w00(j)=w00wk(j)
      enddo

      do k=kb,ke
         do i=0,mx1
            do j= 0,2*my-1
               vor(j,i,k)=vorwk(j,i,k)
               phi(j,i,k)=phiwk(j,i,k)
            enddo
         enddo
      enddo

!!c-------------------- send WxL, WzL
      if (ihist.eq.1) then
         if (myid.eq.numerop-1) then
            chwk (1) = WxL
            chwk (2) = WzL
            chwk (3) = uvL
            call MPI_SEND(chwk,3,MPI_REAL,0,0,MPI_COMM_WORLD,ierr)
         endif
         if (myid.eq.0) then
             call MPI_RECV(chwk,3,MPI_REAL,numerop-1,0,MPI_COMM_WORLD,istat,ierr)
             WxL =  chwk(1)
             WzL =  chwk(2)
             uvL =  chwk(3)
         endif
      endif


!!c                                   /*     write history record     */

      if(myid.eq.0) then
          reynota=0.5*re*re*(abs(wz0/re+uv0)+abs(wzl/re+uvL))
          massu = massu1 +cteu*massu2
          massw = massw1 +ctew*massw2
          if (ihist.eq.1) then
             tmp=my1/2

325   format(i5,9(d14.6))
             write(*,325) istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,u00(floor(tmp))*Re/sqrt(reynota),&
                          massw,uv0,uvL

324   format(17(d22.14))
             write(39,324) time,-1.*Wz0,WzL,sqrt(reynota),ener,u00(floor(tmp))*Re/sqrt(reynota),massw

         endif
      endif

!!c				/* time:
      time=time+Deltat

      if(istati.eq.1) istati = 0
      if(icfl.eq.1)   icfl   = 0 ! if have not applied control this time step
      if(ihist.eq.1)  ihist  = 0
      if(mod(istep,dtc)==0) icfl = 1 !if applied control recalculate timestep

      if (myid.eq.0) then
        totaltimer=totaltimer+MPI_WTIME()
        print *,istep,MPI_WTIME()+iter_time-commtimer+comm_time,commtimer-comm_time,MPI_WTIME()+iter_time
        comm_time = commtimer
      end if

    end do ! Time-stepping loop

    call io_write_footer()

    end subroutine cross1

!c!/********************************************************************
!c/*                                                                  */
!c/*         computes the forcing terms (convective terms)            */
!c/*                                                                  */
!c/*    input:                                                        */
!c/*      phi: Delt(v)   (F-F-Tch)                                    */
!c/*      vor: vorticity (F-F-Tch)                                    */
!c/*      rhg: velocity along y axis (F-F-phys)                       */
!c/*      rhv: dv/dy                 (F-F-phys)                       */
!c/*      work: area de trabajo de dimension al menos  max(mgalx,4*my) */
!c!/*                                                                  */
!c/*   output:                                                        i/
!c/*     rhv: nonlinear term for the phi equation    (F-F-Tch)       */
!c/*     rhg: non linear term for vorticity equation (F-F-Tch)       */
!c/*     rf0u: non linear term for the evolution of Kx=Kz=0  u        */
!c!/*     rf0w: non linear term for the evolution of Kx=Kz=0  w        */
!c/*                                                                  */
!c/*..................................................................*/
!c/*  MODULE FOR MPI SP2                                              */
!c/*..................................................................*
!c/*                                                                  */
!c/*    updated jjs 22/12/00     (INCOMPLETE)                         */
!c/*    single  jjs  4/01/01                                          */
!c/*    low storage : 24/01/01 (incomplete)
!c/********************************************************************/
      subroutine hvhg(phic,ome2c,rhvc,rhgc,rf0u,rf0w,ome1c,work2,sp,myid,rkstep,&
                      u1r,u2r,u3r,o1r,o2r,o3r,u1c,u2c,u3c,o1c,o2c,o3c, &
                      u,v,w,dudy,dvdy,dwdy,tauxy,tauyz,tauxz,press,p_rhs,bcpress,g0,g1,dg0dy,dg1dy,istep)

      use running
      use ctes
      use statistics
      use boundary


      implicit none
      include "ctes3D"
      include "mpif.h"

      real*4 sp(0:mx1,0:nz1,8,jspbb:jspe)
      real*4 uner(9)
      complex*8, dimension(0:mx1,0:mz1,jb:*) :: phic,ome1c,ome2c,rhgc,rhvc

      real*4 rf0u(*),rf0w(*),work2(*)

!!c--------------- 6 * (mgalx+2)  * mgalz planes

      real*4, dimension(mgalx+2,mgalz) :: u1r, u2r, u3r, o1r, o2r, o3r
      complex*8, dimension(0:mx1,0:mz1) :: u1c, u2c, u3c, o1c, o2c, o3c

      complex*8 dk
      real*4 dk2
      integer myid,rkstep,istep
      integer ipo1,ipo2,ipo3
      integer i,j,k,jj,iy,kk,jndex
      integer istat(MPI_STATUS_SIZE),ierr
      integer mmyr
      integer iproc,pproc
      integer pnodes

      real*4 cflx,cfly,cflz,hxalp,hzalp,hyy,cfl0,reigmx1
      integer uggg
      real*8 aa
      complex*16 cc

      real*4 temp
      character*4 ext1, ext2
      character*200 flnm
      integer ilocalu,ilocalw,icount

      complex*8, dimension(0:mx1,0:mz1,jb:je) :: u,v,w,dudy,dvdy,dwdy
      complex*8, dimension(0:mx1,0:mz1,jb:je) :: tauxy,tauyz,tauxz
      complex*8, dimension(0:mx1,0:mz1,jb:je) :: press, p_rhs
      complex(4), dimension(2,0:mx1,kb:ke) :: bcpress
      integer ntot
      real(4),dimension(:),allocatable:: wk1,wk2,wk3

      complex*8, dimension(0:mx1,0:mz1,my,3) :: LSE
      complex*8, dimension(0:mx1,0:mz1,3) :: mes_state
      complex*8, dimension(0:mx1,0:mz1,jb:je) :: v_predict
      complex*8, dimension(0:mx1,0:mz1,my) :: v_pred_2

      complex(4), dimension(0:my1,0:mx1,kb:ke):: g0, g1, dg0dy, dg1dy

      ipo1 = 1    + my
      ipo2 = ipo1 + my
      ipo3 = ipo2 + my

!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
!c  at this point:
!c  rhv is dv / dy -- F-F-P
!c  rhg is v -- F-F-P
!c  phi is nabla^2(v) -- F-F-P
!c  ome1 is d(omega_2)/dy -- F-F-P
!c  ome2 is omega_2 --F-F-P
!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c

!c---------------------- start operating by planes
    !  w_bc = (0., 0.)
!c---------------- Initialize variables out of the y loop

      do kk=1,9
          uner(kk) = 0.
      enddo

      iy = jspb

      cflx = 0.
      cfly = 0.
      cflz = 0.
      cfl0 = 0.

      hxalp=alp*mx*0.5
      hzalp=bet*mz*0.5


      DO J = JB,JE

!c---------------------- computes non 0 modes of ome1, ome3

          do k=1,mz1
             dk  = 1.0/xbet(k)
             o3c(0,k) = -ome1c(0,k,j) * dk
             o1c(0,k) = - phic(0,k,j) * dk
         enddo

          do k=0,mz1
             dk  = xbet(k)
             dk2 = bet2(k)
             do i=1,mx1
                temp = 1.0/(alp2(i) + dk2)
                o3c(i,k) = ( ome1c(i,k,j)*dk-phic(i,k,j)*xalp(i))*temp
                o1c(i,k) = ( ome1c(i,k,j)*xalp(i) + phic(i,k,j)*dk)*temp
            enddo
         enddo

!c---------------------- computes non 0 modes of u,w
          do k=1,mz1
             dk  = 1.0/xbet(k)
             u3c(0,k) = -rhvc (0,k,j) * dk
             u1c(0,k) = ome2c (0,k,j) * dk
         enddo

          do k=0,mz1
             dk  = xbet(k)
             dk2 = bet2(k)
             do i=1,mx1
                temp = 1.0/(alp2(i) + dk2)
                u3c(i,k) = ( rhvc(i,k,j)*dk + ome2c(i,k,j)*xalp(i))*temp
                u1c(i,k) = ( rhvc(i,k,j)*xalp(i) - ome2c(i,k,j)*dk)*temp
            enddo
         enddo

!c---------------------- all nodes, 0 modes, u,w,ome1,ome3

         jj = j-1
         o3c(0,0) = -work2(      j)
         o1c(0,0) =  work2(ipo1+jj)
         u1c(0,0) =  work2(ipo2+jj)
         u3c(0,0) =  work2(ipo3+jj)

!c -------------- copy v and omega_2 into their planes
         do k=0,mz1
            do i=0,mx1
                u2c(i,k) = rhgc(i,k,j)
                o2c(i,k) = ome2c(i,k,j)   !!! warning, ome2c(0,0,j) must be zero
         enddo
      enddo

 ! Copy u and w at the specified nplane_c to use later as boundary condition

  if (j==nplane_c) then
    u_bc(1,:,:) = u1c
    v_bc(1,:,:) = u2c
    w_bc(1,:,:) = u3c
  endif
  if (j==(my-(nplane_c-1))) then
     u_bc(2,:,:) = u1c
     v_bc(2,:,:) = u2c
     w_bc(2,:,:) = u3c
  endif

!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
!c  at this point:
!c  3-D arrays in F-F-P  ---------------------------------
!c  rhv is dv / dy  -----  rhg is v -----phi is nabla^2(v)
!c  ome1 is d(omega_2)/dy ---------------ome2 is omega_2
!c  3-D arrays  in F-F-P  --------------------------------
!c
!c  2-D arrays    --------------
!c  u1 is u ------ u2 is v ------u3 is w
!c  o1 is omega_1 ------ o2 is omega_2 ------ o3 is omega_3
!c  all variables in Fourierx -- Fourier z -- Physical y
!c  2-D arrays    --------------
!c  everything in ONE x-z plane
!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c

u(:,:,j) = u1c ! Fill in big matrices for derivatives and pressure
v(:,:,j) = u2c
w(:,:,j) = u3c

!c======
!c                      /*     statistics of  v, omega_1, omega_3      */
!c                      /*     statistics of  u,w, omega_2             */

if (istati.eq.1 .and. rkstep.eq.1) then
    call stats_addfou(sp,uner,u1c,u2c,u3c,o1c,o2c,o3c,rhvc,iy,j)
endif
!
        if (rkstep.eq.1) then

!c------------- compute vorticity & Re stress  at walls
            if (j.eq.1) then
                Wx0 = o1c(0,0)
                Wz0 = o3c(0,0)
            endif
            if (j.eq.my) then
                WxL = o1c(0,0)
                WzL = o3c(0,0)
            endif

        endif

!c                             /*      Move everything to PPP          */

!c ------- substract umax/2 to u00 to increase dt !!!!
         u1c(0,0) = u1c(0,0) - a0

         call fourxz(u1c,u1r,1,1)         ! u
         call fourxz(u2c,u2r,1,1)         ! v
         call fourxz(u3c,u3r,1,1)         ! w
         call fourxz(o1c,o1r,1,1)         !omega_1
         call fourxz(o2c,o2r,1,1)         !omega_2
         call fourxz(o3c,o3r,1,1)         !omega_3


!c ----------- triple products

         if (rkstep.eq.1) then


            if (istati.eq.1) then
            do k = 1,mgalz
               do i=1,mgalx
                  aa = u2r(i,k)
                  uuv(j) = uuv(j) +aa*u1r(i,k)**2
                  wwv(j) = wwv(j) +aa*u3r(i,k)**2
                  vvv(j) = vvv(j) +aa**3
               enddo
            enddo
            endif

            if (j.eq.1) then
                uv0=0.
                do k=1,mgalz
                   do i=1,mgalx
                     uv0 = uv0 + u1r(i,k)*u2r(i,k)
                   enddo
                enddo
                uv0= uv0/(mgalz*mgalx)
            endif
            if (j.eq.my) then
                uvL=0.
                do k=1,mgalz
                   do i=1,mgalx
                     uvL = uvL + u1r(i,k)*u2r(i,k)
                   enddo
                enddo
                uvL= uvL/(mgalz*mgalx)
            endif

         endif

         uggg=0

         if (icfl.eq.1.and.rkstep.eq.1) then

!c                            /*      estimates maximum time step     */
!c/********************************************************************/
!c/*                                                                  */
!c/*    estimates spatial eigenvalues and computes maximum            */
!c/*    time step.                                                    */
!c/*                                                                  */
!c/*  input :                                                         */
!c/*    rhv,rhg,phi : velocities in (phys-phys-phys) plane            */
!c/*                                                                  */
!c/********************************************************************/
            hyy = hy(j)
            do k=1,mgalz
               do i=1,mgalx
                   cflx = max(cflx,abs(u1r(i,k)) )
                   cfl0 = max(cfl0,abs(u2r(i,k)) )
                   cflz = max(cflz,abs(u3r(i,k)) )
               enddo
            enddo

            cfly = max(cfly,cfl0/hyy)
         endif
!c        /* rhg= H1 = v.omega_3 - w.omega_2 (F-F-P)  */
!c        /* phi= H2 = w.omega_1 - u.omega_3 (F-F-P)  */
!c        /* rhv= H3 = u.omega_2 - v.omega_1 (F-F-P)  */
!c/********************************************************************/
!c/*         computes u X  omega                                      */
!c/********************************************************************/
         do k=1,mgalz
            do i=1,mgalx
               aa = u2r(i,k)
               u2r(i,k) = u2r(i,k)*o3r(i,k)-u3r(i,k)*o2r(i,k)
               u3r(i,k) = u3r(i,k)*o1r(i,k)-u1r(i,k)*o3r(i,k)
               u1r(i,k) = u1r(i,k)*o2r(i,k)-aa*o1r(i,k)
            enddo
         enddo
!c--------------------- at this point
!c--------------------- u1 : H3,  u3 : H2, u2 : H1
!c--------------------- back to F-F-T
         call fourxz(u1c,u1r,-1,1)
         call fourxz(u2c,u2r,-1,1)
         call fourxz(u3c,u3r,-1,1)
!c                          /* saves coefficients for Kx = Kz = 0    */
          rf0u(j)=real(u2c(0,0))
          rf0w(j)=real(u1c(0,0))
!c                         /*   u2 = - dH3/dx + dH1/dz      */
!c                         /*   o1 = dH1/dx + dH3/dz      */

         do k=0,mz1
            do i=0,mx1
               o1c(i,k) = xalp(i)*u2c(i,k)+xbet(k)*u1c(i,k)
               u2c(i,k) = -xalp(i)*u1c(i,k)+xbet(k)*u2c(i,k)
            enddo
         enddo

!c                         /*   rhv = d^2 H2/dx^2 + d^H2/dz^2      */
         do k=0,mz1
            do i=0,mx1
               u1c(i,k) = - u3c(i,k)*(alp2(i)+bet2(k))
            enddo
         enddo

!c --------------------- copy planes into outputs
         do k=0,mz1
            do i=0,mx1
                rhvc(i,k,j)=u1c(i,k)
                rhgc(i,k,j)=u2c(i,k)
                ome1c(i,k,j)=o1c(i,k)
            enddo
         enddo


!c -------------------------- finishes the y loop
      ENDDO

if ((savefields==1).and.(rkstep==1).and.(mod(istep-1,nimag).eq.0).and.(istep.ne.1)) then
    ext2 = 'u'
    call io_phys_field(u,u1c,u1r,ext2)
    ext2 = 'v'
    call io_phys_field(v,u1c,u1r,ext2)
    ext2 = 'w'
    call io_phys_field(w,u1c,u1r,ext2)

    call var_duvw_dy(u,v,w,dudy,dvdy,dwdy)
    call var_stress(u,v,w,dudy,dwdy,tauxy,tauyz,tauxz)

    ext2 = 'txy'
    call io_phys_field(tauxy,u1c,u1r,ext2)
    ext2 = 'tyz'
    call io_phys_field(tauyz,u1c,u1r,ext2)

    ntot=(mgalx+2)*mgalz
    allocate(wk1(ntot),wk2(ntot),wk3(ntot))
    call  var_presrhs(u,v,w,dudy,dvdy,dwdy,p_rhs,wk1,wk2,wk3,wk1,wk2,wk3)
    deallocate(wk1,wk2,wk3)
    call var_pressure(press,p_rhs,bcpress,g0,g1,dg0dy,dg1dy)

    ext2 = 'pneu'
    call io_phys_field(press,u1c,u1r,ext2)
endif


call MPI_Barrier(MPI_COMM_WORLD,ierr);
call MPI_BCAST(u_bc(1,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_up, MPI_COMM_WORLD, ierr)
call MPI_BCAST(v_bc(1,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_up, MPI_COMM_WORLD, ierr)
call MPI_BCAST(w_bc(1,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_up, MPI_COMM_WORLD, ierr)

call MPI_BCAST(u_bc(2,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_down, MPI_COMM_WORLD, ierr)
call MPI_BCAST(v_bc(2,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_down, MPI_COMM_WORLD, ierr)
call MPI_BCAST(w_bc(2,:,:),(mx1+1)*(mz1+1), MPI_COMPLEX, id_down, MPI_COMM_WORLD, ierr)

!c ----------------------- some things have to be done after the y loop

!c ---------------------- adds up the total energy

!c/********************************************************************/
!c/*        computes the energy
!c/********************************************************************/
       if (istati.eq.1.and.rkstep.eq.1.and.ihist.eq.1) then
          call MPI_ALLREDUCE(uner,ener,9,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
          do i=1,9
             ener(i)=sqrt(abs(ener(i)))
          enddo
       endif
!c-------------------- computes Deltat
      if (icfl.eq.1.and.rkstep.eq.1) then
         cflx = cflx*hxalp
         cflz = cflz*hzalp
         cfl0 = max(cflx,max(cfly,cflz))
         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
         if (reigmx1.lt.1e-1)  uggg=1
         if (var_dt==1) then
            Deltat=CFL/reigmx1
         else
            Deltat=Delt
         endif
         dtr=Re/Deltat
         if (uggg.ne.0) then
            write(*,*) 'UGGG', uggg,myid,ihist,jb,je,kb,ke
         endif
      endif

      do iproc=0,numerop-1
         if (iproc.ne.myid) then
           mmyr=jend(iproc)-jbeg(iproc)+1
           call MPI_SENDRECV(rf0u(jb),mmy,MPI_REAL,iproc,0,rf0u(jbeg(iproc)),&
                           mmyr,MPI_REAL,iproc,0,MPI_COMM_WORLD,istat,ierr)

           call MPI_SENDRECV(rf0w(jb),mmy,MPI_REAL,iproc,0,rf0w(jbeg(iproc)),&
                           mmyr,MPI_REAL,iproc,0,MPI_COMM_WORLD,istat,ierr)
         endif
      enddo

!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
!c  at this point: ome1, rhv and rhg are the outs
!c  they must be trasformed from xy-xz before completion
!c  o1: dH1/dx + dH3/dz
!c  u2: -dH3/dx + dH1/dz
!c  u1: d^2H2/dx^2 + d^2H2/dz^2
!c c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c

end subroutine hvhg
