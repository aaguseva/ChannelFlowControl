
subroutine var_duvw_dy(u,v,w,dudy,dvdy,dwdy)
  use ctes
  use mpi_var

  implicit none
  include "ctes3D"

  complex(4),dimension(0:my1,0:mx1,kb:ke):: u, v, w, dudy, dvdy, dwdy
  real(4), dimension(:), allocatable ::chwk1
  integer*8 nbuffsize

  nbuffsize = mx*max(mmy*mz,mmz*my)
  allocate(chwk1(nbuffsize))
  call chikj2jik(u,u,chwk1,chwk1,myid)
  call chikj2jik(v,v,chwk1,chwk1,myid)
  call chikj2jik(w,w,chwk1,chwk1,myid)

  call deryr2(u,dudy,(mx1+1)*mmz,my,chwk1)
  call deryr2(v,dvdy,(mx1+1)*mmz,my,chwk1)
  call deryr2(w,dwdy,(mx1+1)*mmz,my,chwk1)

  call chjik2ikj(u,u,chwk1,chwk1,myid)
  call chjik2ikj(v,v,chwk1,chwk1,myid)
  call chjik2ikj(w,w,chwk1,chwk1,myid)

  call chjik2ikj(dudy,dudy,chwk1,chwk1,myid)
  call chjik2ikj(dvdy,dvdy,chwk1,chwk1,myid)
  call chjik2ikj(dwdy,dwdy,chwk1,chwk1,myid)

  deallocate(chwk1)

end subroutine var_duvw_dy


subroutine var_stress(u,v,w,dudy,dwdy,tauxy,tauyz,tauxz)
  use ctes
  use mpi_var
  implicit none
  include "mpif.h"
  include "ctes3D"

  complex(4),dimension(0:mx1,0:mz1,jb:je)::u,v,w,dudy,dwdy,tauxy,tauyz,tauxz
  integer i,j,k

       !----------------Shear stresses------------------!
  do j=jb,je
    do k=0,mz1
        do i=0,mx1
               tauyz(i,k,j)   = dwdy(i,k,j)+v(i,k,j)*xbet(k)      ! tau_yz: dwdy+dvdz
               tauxz(i,k,j)   = u(i,k,j)*xbet(k)+w(i,k,j)*xalp(i) ! tau_xz: dudz+dwdx
               tauxy(i,k,j)   = v(i,k,j)*xalp(i)+dudy(i,k,j)      ! tau_xy: dvdx+dudy
          enddo
     enddo
  enddo

endsubroutine var_stress


subroutine var_pressure(presjik,pres,bcpress,g0,g1,dg0dy,dg1dy)
  use ctes
  use mpi_var
  use planes

  implicit none
  include "mpif.h"
  include "ctes3D"

  complex(4), dimension(0:my1,0:mx1,kb:ke):: presjik, presj, pres, dpdy
  complex(4), dimension(0:my1,0:mx1,kb:ke):: g0, g1, dg0dy, dg1dy
  complex(4), dimension(2,0:mx1,kb:ke) ::  bcpress
  real(4), dimension(:), allocatable ::chwk1
  complex(4) :: pcoef0, pcoef1
  complex*8 bcb,bct
  real*8 rK
  integer nbuffsize, i, k, k1
  character*4 ext1, ext2
  complex(4), dimension(0:my1,0:mx1,kb:ke) ::  coef0, coef1

  nbuffsize = mx*max(mmy*mz,mmz*my)
  allocate(chwk1(nbuffsize))
  call chikj2jik(pres,pres,chwk1,chwk1,myid)

  do k=kb,ke
     do i=0,mx1
        if ((k.eq.1).and.(i.eq.0)) then
            presjik(0:my1,i,k)=(0.0,0.0)
        else
          k1 = icx(k-1)
          rK = bet2(k1)+alp2(i)

          bcb = (0.0, 0.0)
          bct = (0.0, 0.0)
          call Lapvdv(pres(0,i,k),presj(0,i,k),dpdy(0,i,k),rk,bcb,bct)

          pcoef1 = (bcpress(1,i,k) - dpdy(0,i,k))*dg0dy(my1,i,k)/dg0dy(0,i,k)
          pcoef1 = (bcpress(2,i,k) - dpdy(my1,i,k)) - pcoef1
          pcoef1 = pcoef1/(dg1dy(my1,i,k) - (dg0dy(my1,i,k)/dg0dy(0,i,k))*dg1dy(0,i,k))
          pcoef0 = (bcpress(1,i,k) - dpdy(0,i,k) - dg1dy(0,i,k)*pcoef1)/dg0dy(0,i,k)

          presjik(0:my1,i,k) = presj(0:my1,i,k) + pcoef0*g0(0:my1,i,k) + pcoef1*g1(0:my1,i,k)

        endif

      enddo
   enddo
   call chjik2ikj(presjik,presjik,chwk1,chwk1,myid)
   deallocate(chwk1)

end subroutine var_pressure

subroutine var_presrhs(u,v,w,dudy,dvdy,dwdy,pres,wk1,wk2,wk3,wk1r,wk2r,wk3r)
  use ctes

  implicit none
  include "mpif.h"
  include "ctes3D"

  integer i,k,kk,j,jj,iplan,nvecy
  complex(4), dimension(0:mx1,0:mz1,jb:je):: u,v,w,dudy,dvdy,dwdy,pres
  complex(4), dimension(0:mx1,0:mz1):: wk1,wk2,wk3
  real(4), dimension(mgalx+2 ,mgalz):: wk1r,wk2r,wk3r
  ! Input: data should be transposed to xz planes
  do j = jb,je
    ! dudy and dvdx
    do i=0,mx1
       wk1(i,:) = xalp(i)*v(i,:,j) !dvdx(nplany,0:mz1,mx1+1)
    enddo
    wk3 = dudy(:,:,j)  ! 4
    call fourxz(wk1,wk1,1,1)    !
    call fourxz(wk3,wk3,1,1)
    wk1r = 2d0*wk1r*wk3r   ! 2 and 4

   !dvdz and dwdy
   do k=0,mz1
      wk2(:,k) = xbet(k)*v(:,k,j) !uyzx(nplany,0:mz1,mx1+1)
   enddo
   wk3 = dwdy(:,:,j)
   call fourxz(wk2,wk2,1,1)    !
   call fourxz(wk3,wk3,1,1)
   wk1r = wk1r + 2d0*wk2r*wk3r

   !dudz and dwdx
   do i=0,mx1
      wk2(i,:) = xalp(i)*w(i,:,j)
   enddo
   do k=0,mz1
      wk3(:,k) = xbet(k)*u(:,k,j)
   enddo
   call fourxz(wk2,wk2,1,1)
   call fourxz(wk3,wk3,1,1)
   wk1r = wk1r + 2d0*wk2r*wk3r

   !(dvdy)^2 and (dwdz)^2
   wk2 = dvdy(:,:,j) !uyzx(nplany,0:mz1,mx1+1)
   do k=0,mz1
      wk3(:,k) = xbet(k)*w(:,k,j)
   enddo
   call fourxz(wk2,wk2,1,1)
   call fourxz(wk3,wk3,1,1)
   wk1r = wk1r + wk2r*wk2r + wk3r*wk3r

    !(dudx)^2
    do i=0,mx1
       wk2(i,:) = xalp(i)*u(i,:,j)
    enddo
    call fourxz(wk2,wk2,1,1)
    wk1r =  -(wk1r + wk2r*wk2r) ! RHS for each j is here

    call fourxz(wk1,wk1,-1,1)     !   back to FFP
    pres(:,:,j) = wk1             ! copy RHS to pressure
  end do

end subroutine var_presrhs

subroutine var_precompute(g0, g1, dg0dy, dg1dy)
  use ctes
  implicit none
  include "mpif.h"
  include "ctes3D"

  complex(4), dimension(0:my1,0:mx1,kb:ke):: g0, g1, dg0dy, dg1dy, zero_rhs
  complex*8 bcb,bct
  real*8 rK
  integer  i, k, k1

  zero_rhs = (0.0,0.0)

  do k=kb,ke
     do i=0,mx1
        if ((k.eq.1).and.(i.eq.0)) then
            g0(0:my1,i,k)=(0.0,0.0)
            g1(0:my1,i,k)=(0.0,0.0)
            dg0dy(0:my1,i,k)=(0.0,0.0)
            dg1dy(0:my1,i,k)=(0.0,0.0)
        else
          k1 = icx(k-1)
          rK = bet2(k1)+alp2(i)
          bcb = (1.0, 0.0)
          bct = (0.0, 0.0)
          call Lapvdv(zero_rhs(0,i,k),g0(0,i,k),dg0dy(0,i,k),rk,bcb,bct)
          bcb = (0.0, 0.0)
          bct = (1.0, 0.0)
          call Lapvdv(zero_rhs(0,i,k),g1(0,i,k),dg1dy(0,i,k),rk,bcb,bct)
        endif

      enddo
   enddo

end subroutine var_precompute

subroutine var_bc_press(phi, bcpress,rkstep)
    use ctes
    use planes
    use mpi_var
    implicit none

    include "ctes3D"

    complex(4) bcpress(2,0:mx1,kb:ke), phi(0:my1,0:mx1,kb:ke)
    real(4), dimension(:), allocatable ::chwk1
    integer nbuffsize, rkstep
    character*4 ext2

    bcpress(1,:,:) =  phi(0,:,:)/re
    bcpress(2,:,:) =  phi(my1,:,:)/re

end subroutine var_bc_press
