      subroutine genexp
      !use wave
      use ctes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      computes wavenumbers and fourier indices
c
c      updated j.j.s.     22/12/00
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include "ctes3D"

      real zero
      integer i,j,k

      zero = 0e0

      do 10 k=0,nz1
         xbet(k) = cmplx(zero,bet*k)
         icx(k) = k
 10   continue

      do 20 k=nz1+1,mz1
         xbet(k) = cmplx(zero ,-bet*(mz1+1-k))
 20   continue

      do 30 k=1,nz1
         icx(mz-k) = k
 30   continue

      do 40 i=0,mx1
         iax(2*i+1) = i
         iax(2*i+2) = i
 40   continue

      do i=0,mx1
         xalp(i) = cmplx(zero ,alp*i)
      enddo

      do i=0,mx1
         alp2(i) = -xalp(i)**2
      enddo

      do j=0,mz1
         bet2(j) = -xbet(j)**2
      enddo

      end
