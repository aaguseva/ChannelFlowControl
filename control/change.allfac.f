       subroutine chikj2jik(xy,xz,wk1,wk2,myid)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV               c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
       use ctes

       implicit none
       include 'mpif.h'
       include "ctes3D"
       integer myid,icontrol,iresult,ilocal,icount

       real*4 xy(mx,mz,*),xz(2*my,mx1+1,kb:*),wk1(*),
     &        wk2(mx,kb:ke,*)

       integer istat(MPI_STATUS_SIZE),ierr
       integer iproc,nsetotr,ipoxz,mzbu,ntoto
       integer i,j,k,mmz2

c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c NOTE:
cc  wk2 dimensioned at least (mx*mzp*myp)
c/cccccccccccccccccccccccccccccccccccccccccccccccc/

cc-------------------------- checks NANS 1/2 -----------------------
c       icontrol=0
c       icount = 0
c       ilocal = 0
c
c       do j=1,je-jb+1
c          do k=1,mz
c             do i=1,mx
c                if (xy(i,k,j).ne.xy(i,k,j)) then
c                   icontrol=1
c                   icount = icount + 1
c                   ilocal = 1
c                endif
c             enddo
c          enddo
c       enddo
c
c       if (ilocal.eq.1) then
c          write(*,*) 'NaN in chikj2jik before MPI_SENDRECV,myid=,',
c     .                myid,'   count=',icount
c          call system("hostname")
c       endif
cc--------------------------end checks NANS 1/2 -----------------------


       mmz2=mmz*mx

       if (myid.eq.0) then
         commtimer = commtimer-MPI_WTIME()
       endif

       do iproc=0,numerop-1

          if(iproc.ne.myid)then

             nsetotr=mmz2*(jend(iproc)-jbeg(iproc)+1)
             ipoxz=1+mmz2*(jbeg(iproc)-1)

             call MPI_SENDRECV(xy(1,kbeg(iproc),1),1,myslice(iproc),
     .            iproc,0,wk1(ipoxz),nsetotr,
     .            MPI_REAL,iproc,0,
     .            MPI_COMM_WORLD,istat,ierr)

          else

            mzbu=kend(iproc)-kbeg(iproc)+1
            ipoxz=1+mmz2*(jb-1)
            call pacy2z(wk1(ipoxz),xy,mzbu,kbeg(iproc),iproc)

          endif

       enddo

       if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
       endif


       if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
       endif

c       do j=1,my
c          do k=kb,ke
c             do i=1,mx1+1
c                xz(2*j-1,i,k) = wk2(2*i-1,k,j)
c                xz(2*j  ,i,k) = wk2(2*i  ,k,j)
c             enddo
c          enddo
c       enddo
       do k=kb,ke
             do j=1,my
          do i=1,mx1+1
                xz(2*j-1,i,k) = wk2(2*i-1,k,j)
                xz(2*j  ,i,k) = wk2(2*i  ,k,j)
             enddo
          enddo
       enddo


       if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
       endif

cc-----------------------------checks NANS 2/2-----------------------
c       icount=0
c       ilocal=0
c       do k=kb,ke
c          do i=1,mx1+1
c             do j=1,2*my
c                if (xz(j,i,k).ne.xz(j,i,k)) then
c                   icontrol=1
c                   ilocal=1
c                   icount = icount + 1
c                endif
c             enddo
c          enddo
c       enddo
c
c       if (ilocal.eq.1) then
c         write(*,*) 'NaN in chikj2jik after MPI_SENDRECV,myid=,',
c     .               myid,'   count:',icount
c         call system("hostname")
c       endif
c
c       call MPI_ALLREDUCE(icontrol,iresult,1,MPI_INTEGER,
c     .                         MPI_MAX,MPI_COMM_WORLD,ierr)
c
c       if (iresult.eq.1) then
c          write(*,*) 'Finishing process ... ',myid
c          stop
c       endif
cc--------------------------end checks NANS 2/2 -----------------------

       end


       subroutine chjik2ikj(xz,xy,wk1,wk2,myid)
       use ctes
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
c MODULE FOR SP2 MPI uses SENDRECV               c/
c sends a block (mx,mz,jb:je) to a (mx,kb:ke,my) c/
c/cccccccccccccccccccccccccccccccccccccccccccccccc/
       implicit none
       include "mpif.h"
       include "ctes3D"

       real*4 xz(2*my,mx1+1,kb:*),xy(mx,mz,*),wk1(*),wk2(mx,kb:ke,*)
       integer myid,icontrol,ilocal,icount,iresult

       integer istat(MPI_STATUS_SIZE),ierr
       integer iproc,nsetots,ipoxz,mzbu,ntoto
       integer i,j,k,mmy2,mmz2

       mmz2=mmz*mx
       mmy2=mmy*mx

cc-------------------------- checks NANS 1/2 -----------------------
c       icontrol=0
c       icount = 0
c       ilocal = 0
c
c       do k=kb,ke
c          do i=1,mx1+1
c             do j=1,2*my
c                if (xz(j,i,k).ne.xz(j,i,k)) then
c                   icontrol=1
c                   icount = icount + 1
c                   ilocal = 1
c                endif
c             enddo
c          enddo
c       enddo
c
c       if (ilocal.eq.1) then
c          write(*,*) 'NaN in chjik2ikj before MPI_SENDRECV,myid=,',
c     .                myid,'   count=',icount
c          call system("hostname")
c       endif
cc--------------------------end checks NANS 1/2 -----------------------

       if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
       endif


c       do k=kb,ke
c          do j=1,my
c             do i=1,mx1+1
c                wk2(2*i-1,k,j) = xz(2*j-1,i,k)
c                wk2(2*i  ,k,j) = xz(2*j  ,i,k)
c             enddo
c          enddo
c       enddo
      do k=kb,ke
         do j=1,my
            do i=1,mx1+1
               wk2(2*i-1,k,j) = xz(2*j-1,i,k)
               wk2(2*i  ,k,j) = xz(2*j  ,i,k)
            enddo
         enddo
      enddo

      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
         commtimer = commtimer-MPI_WTIME()
      endif



       do iproc=0,numerop-1

c           write(*,*) 'myid,iproc',myid,iproc
          if(iproc.ne.myid)then

             nsetots=(jend(iproc)-jbeg(iproc)+1)*mmz2
             ipoxz=1+mmz2*(jbeg(iproc)-1)

             call MPI_SENDRECV(wk1(ipoxz),nsetots,
     .            MPI_REAL,iproc,0,xy(1,kbeg(iproc),1),
     .            1,myslice(iproc),
     .            iproc,0,MPI_COMM_WORLD,istat,ierr)

          else

           mzbu=kend(iproc)-kbeg(iproc)+1
           ipoxz=1+mmz2*(jb-1)
           call unpacz2y(wk1(ipoxz),xy,mzbu,kbeg(iproc),iproc)

          endif

      enddo

      if (myid.eq.0) then
        commtimer = commtimer+MPI_WTIME()
      endif

cc-------------------------- checks NANS 2/2 -----------------------
c       icount = 0
c       ilocal = 0
c
c       do j=1,je-jb+1
c          do k=1,mz
c             do i=1,mx
c                if (xy(i,k,j).ne.xy(i,k,j)) then
c                   icontrol=1
c                   icount = icount + 1
c                   ilocal = 1
c                endif
c             enddo
c          enddo
c       enddo
c
c       if (ilocal.eq.1) then
c          write(*,*) 'NaN in chjik2ikj after MPI_SENDRECV,myid=,',
c     .                myid,'   count=',icount
c          call system("hostname")
c       endif
c
c       call MPI_ALLREDUCE(icontrol,iresult,1,MPI_INTEGER,
c     .                         MPI_MAX,MPI_COMM_WORLD,ierr)
c
c       if (iresult.eq.1) then
c          write(*,*) 'Finishing process ... ',myid
c          stop
c       endif
cc--------------------------end checks NANS 2/2 -----------------------

      end


      subroutine unpacz2y(xyi,xyo,mzbu,kb1,iproc)
      use ctes
      implicit none
      include "ctes3D"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
!      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
!     .               kbeg(0:numerop-1),kend(0:numerop-1),
!     .               jb,je,kb,ke,mmy,mmz
!      save /point/
      integer mzbu,kb1,iproc
      real*4 xyi(mx,kb1:kb1+mzbu-1,*),xyo(mx,mz,*)
      integer  i,j,k

c     kb1 = kbeg(iproc)
      do j=1,mmy
         do k=kbeg(iproc),kend(iproc)
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end


      subroutine pacy2z(xyo,xyi,mzbu,kb1,iproc)
      use ctes
      implicit real*4(a-h,o-z)
      include "ctes3D"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c IN xzo OUT xzi                                               c
c  unpack from kb(iproc) till ke(iproc)                        c
c  and    from jb(myid)  till je(myid)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
!     .               kbeg(0:numerop-1),kend(0:numerop-1),
!     .               jb,je,kb,ke,mmy,mmz
!
!      save /point/
      integer mzbu,kb1,iproc
      dimension xyi(mx,mz,*),xyo(mx,kb1:kb1+mzbu-1,*)
      integer i,j,k

c     kb1 = kbeg(iproc)
      do j=1,mmy
         do k=kbeg(iproc),kend(iproc)
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end

      subroutine pointers(jb,je,kb,ke)
c ===============================================
c      jjs,  aug/2001  (bug in old version)
c ===============================================
      implicit none

      include "ctes3D"
      integer jb(numerop),je(numerop),kb(numerop),ke(numerop),n,n1,n2

      n1=my/numerop
      n2=my-numerop*n1

      jb(1)=1
      do n=1,n2
         je(n)  = jb(n)+n1
         jb(n+1)= je(n)+1
      enddo
      do n=n2+1,numerop-1
         je(n)=jb(n)+n1-1
         jb(n+1)= je(n)+1
      enddo
      je(numerop)=jb(numerop)+n1-1


      n1=mz/numerop
      n2=mz-numerop*n1

      kb(1)=1
      do n=1,n2
         ke(n)  = kb(n)+n1
         kb(n+1)= ke(n)+1
      enddo
      do n=n2+1,numerop-1
         ke(n)=kb(n)+n1-1
         kb(n+1)= ke(n)+1
      enddo
      ke(numerop)=kb(numerop)+n1-1

      end


c       ------- OLD VERSION! IT HAS A BUG!
c       -----------------------------------------
c       subroutine pointers(jb,je,kb,ke)
c       include "ctes3D"
c       integer jb(numerop),je(numerop),kb(numerop),ke(numerop)
c
c
c       n1=my/numerop
c       n2=my-numerop*n1
c
c       jb(1)=1
c       je(1)=jb(1)+n1
c       do n=2,n2
c          jb(n)=je(n-1)+1
c          je(n)=jb(n)+n1
c       enddo
c       do n=n2+1,numerop
c          jb(n)=je(n-1)+1
c          je(n)=jb(n)+n1-1
c       enddo
c
c
c       n1=mz/numerop
c       n2=mz-numerop*n1
c
c       kb(1)=1
c       ke(1)=kb(1)+n1
c       do n=2,n2
c          kb(n)=ke(n-1)+1
c          ke(n)=kb(n)+n1
c       enddo
c       do n=n2+1,numerop
c          kb(n)=ke(n-1)+1
c          ke(n)=kb(n)+n1-1
c       enddo
c
c       end
