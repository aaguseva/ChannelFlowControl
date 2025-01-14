subroutine stats_init()
  use statistics
  implicit none
  include "ctes3D"

  um(1:my)  = 0.
  vm(1:my)  = 0.
  wm(1:my)  = 0.
  up(1:my)  = 0.
  vp(1:my)  = 0.
  wp(1:my)  = 0.
  uvr(1:my) = 0.
  uwr(1:my) = 0.
  vwr(1:my) = 0.
  w1m(1:my) = 0.
  w2m(1:my) = 0.
  w3m(1:my) = 0.
  w1p(1:my) = 0.
  w2p(1:my) = 0.
  w3p(1:my) = 0.
  ep(1:my) = 0.
  uuv(1:my) = 0.
  wwv(1:my) = 0.
  vvv(1:my) = 0.

  nacum = 0
end subroutine stats_init

subroutine stats_send(myid)
use statistics
use ctes
!use point
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          single jjs 4/01/01
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
include "mpif.h"
include "ctes3D"

integer istat(MPI_STATUS_SIZE)

 if(myid.ne.0) then
    ipo=jb
    leng=mmy
    call MPI_SEND(um(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(vm(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(wm(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(up(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(vp(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(wp(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w1m(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w2m(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w3m(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w1p(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w2p(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(w3p(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(uvr(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(uwr(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(vwr(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(ep(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(uuv(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(wwv(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(vvv(ipo),leng,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
 else
    do iproc=1,numerop-1
       ipo=jbeg(iproc)
       leng=jend(iproc)-jbeg(iproc)+1
       call MPI_RECV(um(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(vm(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(wm(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(up(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(vp(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(wp(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w1m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w2m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w3m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w1p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w2p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(w3p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(uvr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(uwr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(vwr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(ep(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(uuv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(wwv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
       call MPI_RECV(vvv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
   enddo
endif

end subroutine stats_send

subroutine stats_addfou(sp,uner,u1c,u2c,u3c,o1c,o2c,o3c,rhvc,iy,j)
  use statistics
  use running
  use ctes
  use mpi_var

  implicit none
  include "mpif.h"
  include "ctes3D"
  complex(4), dimension(0:mx1,0:mz1):: u1c,u2c,u3c,o1c,o2c,o3c
  complex(4),dimension(0:mx1,0:mz1,jb:je):: rhvc
  real*4 sp(0:mx1,0:nz1,8,jspbb:jspe)
  real*4 uner(9), hyy
  real*8 aa
  complex*16 cc
  integer i,j,k,kk,iy


  !c -----------------  spectra
  !c ---------------- only if my plane contains spectra information

  if (jsp(j).eq.1) then
      do kk = 0,mz1
        k = icx(kk)
  !c                  write(*,*) j,iy,jspb,jspe
        sp(0,k,1,iy) = sp(0,k,1,iy)+u1c(0,kk)*conjg(u1c(0,kk))
        sp(0,k,2,iy) = sp(0,k,2,iy)+u2c(0,kk)*conjg(u2c(0,kk))
        sp(0,k,3,iy) = sp(0,k,3,iy)+u3c(0,kk)*conjg(u3c(0,kk))
        sp(0,k,4,iy) = sp(0,k,4,iy)+real(u1c(0,kk)*conjg(u2c(0,kk)))
        sp(0,k,5,iy) = sp(0,k,5,iy)+aimag(u1c(0,kk)*conjg(u2c(0,kk)))
        sp(0,k,6,iy) = sp(0,k,6,iy)+o1c(0,kk)*conjg(o1c(0,kk))
        sp(0,k,7,iy) = sp(0,k,7,iy)+o2c(0,kk)*conjg(o2c(0,kk))
        sp(0,k,8,iy) = sp(0,k,8,iy)+o3c(0,kk)*conjg(o3c(0,kk))
        do i = 1,mx1
           sp(i,k,1,iy) = sp(i,k,1,iy)+2.*u1c(i,kk)*conjg(u1c(i,kk))
           sp(i,k,2,iy) = sp(i,k,2,iy)+2.*u2c(i,kk)*conjg(u2c(i,kk))
           sp(i,k,3,iy) = sp(i,k,3,iy)+2.*u3c(i,kk)*conjg(u3c(i,kk))
           sp(i,k,4,iy) =sp(i,k,4,iy)+2.*real(u1c(i,kk)*conjg(u2c(i,kk)))
           sp(i,k,5,iy) =sp(i,k,5,iy)+2.*aimag(u1c(i,kk)*conjg(u2c(i,kk)))
           sp(i,k,6,iy) = sp(i,k,6,iy)+2.*o1c(i,kk)*conjg(o1c(i,kk))
           sp(i,k,7,iy) = sp(i,k,7,iy)+2.*o2c(i,kk)*conjg(o2c(i,kk))
           sp(i,k,8,iy) = sp(i,k,8,iy)+2.*o3c(i,kk)*conjg(o3c(i,kk))
        enddo
     enddo
  !c ----------------  update spectra y index
  iy = iy + 1
  endif
  !c---------- just add 1 to nacum once !!!!
  if (j.eq.je)  nacumsp = nacumsp +1

  !c                    ----- intensities ------------
  !c                    ----- & dissipation ------------
  do kk=1,9
      ener(kk) = 0.
  enddo
  do k=0,mz1
  !c              intensities ----------------
      aa = u1c(0,k)*conjg(u1c(0,k))
      up(j) = up(j) + aa
      ener(1)=ener(1) + aa

      aa= u2c(0,k)*conjg(u2c(0,k))
      vp(j) = vp(j) + aa
      ener(2)=ener(2) + aa

      aa= u3c(0,k)*conjg(u3c(0,k))
      wp(j) = wp(j) + aa
      ener(3)=ener(3) + aa

      aa= real(u1c(0,k)*conjg(u2c(0,k)))
      uvr(j)= uvr(j) + aa
      ener(4)=ener(4) + aa

      aa= real(u1c(0,k)*conjg(u3c(0,k)))
      uwr(j)= uwr(j) + aa
      ener(5)=ener(5) + aa

      aa= real(u3c(0,k)*conjg(u2c(0,k)))
      vwr(j)= vwr(j) + aa
      ener(6)=ener(6) + aa

      aa= o1c(0,k)*conjg(o1c(0,k))
      w1p(j)= w1p(j) + aa
      ener(7)=ener(7) + aa

      aa = o2c(0,k)*conjg(o2c(0,k))
      w2p(j)= w2p(j) + aa
      ener(8)=ener(8) + aa

      aa = o3c(0,k)*conjg(o3c(0,k))
      w3p(j)= w3p(j) + aa
      ener(9)=ener(9) + aa
  !c              dissipation  ----------------
      aa =  bet2(k) *( u1c(0,k)*conjg(u1c(0,k)) +u2c(0,k)*conjg(u2c(0,k)) +&
            u3c(0,k)*conjg(u3c(0,k)) ) +rhvc(0,k,j)*conjg(rhvc(0,k,j))

      cc = o1c(0,k) + xbet(k)*u2c(0,k)
      aa = aa + cc*conjg(cc)
      cc = o3c(0,k)
      aa = aa + cc*conjg(cc)

      ep(j) = ep(j) + aa

      do i=1,mx1
          aa = 2.*u1c(i,k)*conjg(u1c(i,k))
          up(j) = up(j) + aa
          ener(1)=ener(1) + aa

          aa= 2.*u2c(i,k)*conjg(u2c(i,k))
          vp(j) = vp(j) + aa
          ener(2)=ener(2) + aa

          aa= 2.*u3c(i,k)*conjg(u3c(i,k))
          wp(j) = wp(j) + aa
          ener(3)=ener(3) + aa

          aa= 2.*real(u1c(i,k)*conjg(u2c(i,k)))
          uvr(j)= uvr(j) + aa
          ener(4)=ener(4) + abs(aa)

          aa= 2.*real(u1c(i,k)*conjg(u3c(i,k)))
          uwr(j)= uwr(j) + aa
          ener(5)=ener(5) + abs(aa)

          aa= 2.*real(u3c(i,k)*conjg(u2c(i,k)))
          vwr(j)= vwr(j) + aa
          ener(6)=ener(6) + abs(aa)

          aa= 2.*o1c(i,k)*conjg(o1c(i,k))
          w1p(j)= w1p(j) + aa
          ener(7)=ener(7) + aa

          aa = 2.*o2c(i,k)*conjg(o2c(i,k))
          w2p(j)= w2p(j) + aa
          ener(8)=ener(8) + aa

          aa = 2.*o3c(i,k)*conjg(o3c(i,k))
          w3p(j)= w3p(j) + aa
          ener(9)=ener(9) + aa

  !c              dissipation  ----------------

          aa = (alp2(i) + bet2(k))*( u1c(i,k)*conjg(u1c(i,k)) +&
                u2c(i,k)*conjg(u2c(i,k)) +u3c(i,k)*conjg(u3c(i,k)) ) +&
                rhvc(i,k,j)*conjg(rhvc(i,k,j) )

          cc = o1c(i,k) + xbet(k)*u2c(i,k)
          aa = aa + cc*conjg(cc)
          cc = o3c(i,k) - xalp(i)*u2c(i,k)
          aa = aa + cc*conjg(cc)

          ep(j) = ep(j) + 2.*aa

        enddo
    enddo

  !cc --------------- add this plane energy

    hyy = hy(j)
    do kk = 1,9
       uner(kk) = uner(kk) + ener(kk)*hyy
    enddo

  !c ------------ means

    um(j) = um(j)+u1c(0,0)
    vm(j) = vm(j)+u2c(0,0)
    wm(j) = wm(j)+u3c(0,0)
    w1m(j)= w1m(j)+o1c(0,0)
    w2m(j)= w2m(j)+o2c(0,0)
    w3m(j)= w3m(j)+o3c(0,0)

  !c-------------- update nacum just once !!!

    if (j.eq.je)   nacum = nacum+1

    if (myid.eq.0) then
       Wx0a=Wx0a+o1c(0,0)
       Wz0a=Wz0a+o3c(0,0)
    endif


end subroutine stats_addfou
