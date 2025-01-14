
subroutine io_write_state(vor,phi,u00,w00,buff1,buff2)
  use running
  use ctes
  use planes
  use boundary
  use mpi_var

  implicit none
  include "mpif.h"
  include "ctes3D"

  character*200 fil
  character*4 ext1
  integer i,j,k,ntotr,leng

  real(4), dimension(mx,mz,jb:je):: vor,phi,buff1,buff2
  real(4) u00(*),w00(*)

  !      ---------------  copy everything
  call chjik2ikj(vor,vor,buff1,buff1,myid)
  call chjik2ikj(phi,phi,buff1,buff1,myid)

  buff1(:,:,jb:je) = vor(:,:,jb:je)
  buff2(:,:,jb:je) = phi(:,:,jb:je)

  if(id22.gt.9999.or.id22.lt.0) then
     write(*,*) 'number of images out of range'
     stop
  endif
  !------- collect data in the master, and write it to file
  if(myid==0) then

     write(ext1,'(i4.4)') id22
     fil = trim(filout)//'.'//ext1

     open (iout,file=fil,status='unknown',form='unformatted',access='stream')
     rewind(iout)
     write(iout) time,Re,alp,bet,a0,mx,my,mz,imesh,gamma, &
                 mxmin,mzmin,mxmax,mzmax,alp_u,alp_v,alp_w,t_lag,x_lag,(j,j=1,11)
     write(iout) (u00(j),w00(j),j=1,my)
     do j=jb,je
        write(iout) ((buff1(i,k,j),buff2(i,k,j),i=1,mx),k=1,mz)
     enddo

     do iproc=1,numerop-1
        leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
        call MPI_RECV(buff1,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        call MPI_RECV(buff2,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        do j=jb,jb+jend(iproc)-jbeg(iproc)
           write(iout) ((buff1(i,k,j),buff2(i,k,j),i=1,mx),k=1,mz)
        enddo

     enddo
     close(iout)

  else
     call MPI_SEND(buff1,mx*mz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
     call MPI_SEND(buff2,mx*mz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
  endif

  ! ----- before going any further exchange cuts in y with cuts in z!
  call chikj2jik(vor,vor,buff1,buff1,myid)
  call chikj2jik(phi,phi,buff1,buff1,myid)

end subroutine io_write_state

subroutine io_calc_spec(sp,spwk)
    use statistics
    use mpi_var
    use ctes

    implicit none
    include "mpif.h"
    include "ctes3D"

    real*4 sp(0:mx1,0:nz1,8,jspbb:jspe),spwk(0:mx1,0:nz1,8,1:*)
    integer i, j, jj
    integer leng, leng1, leng2

  !c------- upper half sends spectra
              !DEC$ NOVECTOR
              do j=2*nspec+1,nspec+1,-1
                 jj = 2*nspec+2 - j
                 if (myid.eq.jspiproc(j).and.myid.ne.jspiproc(jj)) then
                    leng = (mx1+1)*(nz1+1)*8
                    call MPI_SEND(sp(0,0,1,j),leng,MPI_REAL,jspiproc(jj),0,MPI_COMM_WORLD,ierr)
                 endif
  !c------- lower half receives upper half spectra and computes average
                 if (myid.eq.jspiproc(jj)) then
                    leng = (mx1+1)*(nz1+1)*8
                    leng1 = (mx1+1)*(nz1+1)*3
                    leng2 = (mx1+1)*(nz1+1)*5
                    if (jspiproc(j).ne.myid) then
                       call MPI_RECV(spwk(0,0,1,1),leng,MPI_REAL,jspiproc(j),0,MPI_COMM_WORLD,istat,ierr)
                    else
                      do i=0,leng-1
                         spwk(i,0,1,1) = sp(i,0,1,j)
                      enddo
                    endif
  !c------------- velocity spectra are symmetric
                    do i = 0,leng1-1
                       spwk(i,0,1,1) =.5*(spwk(i,0,1,1) + sp(i,0,1,jj))
                    enddo
  !c------------- velocity cospectra are skew-symmetric
                    do i = leng1,leng2-1
                       spwk(i,0,1,1) =.5*(-spwk(i,0,1,1) + sp(i,0,1,jj))
                   enddo
  !c------------- vorticity spectra are symmetric
                    do i = leng2,leng-1
                      spwk(i,0,1,1) =.5*(spwk(i,0,1,1) + sp(i,0,1,jj))
                    enddo
  !c------- everybody sends data to the master
                    if(myid.ne.0) then
  !c-------  only lower half sends averaged spectra to the master
                       leng = (mx1+1)*(nz1+1)*8
                       call MPI_SEND(spwk,leng,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
                    else
  !c------  the master first writes its stuff
                       call io_write_spec(spwk,jj)
                    endif
                 endif
                 if (myid.eq.0.and.jspiproc(jj).ne.myid) then
  !c------ then receives everything from everybody
                    leng = (mx1+1)*(nz1+1)*8
                    call MPI_RECV(spwk,leng,MPI_REAL,jspiproc(jj),jspiproc(jj),MPI_COMM_WORLD,istat,ierr)
  !c------- and writes it
                    call io_write_spec(spwk, jj)
                 endif
              enddo
end subroutine io_calc_spec

subroutine io_write_spec(sp,jjsp)
  use statistics
  use running
  use ctes


  implicit none
  include "mpif.h"
  include "ctes3D"

  real*4 sp(0:mx1,0:nz1,8)
  integer i,j,k, jjsp
  character*4 ext1
  character*200 fnamespe


  if (nstart.ne.0.and.nacumsp.ne.0) then
!c                                 /*       write spectra        */
     if (jjsp.eq.1) then

        write(ext1,'(i4.4)') id22
        fnamespe=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.spe'
        open (ispf,file=fnamespe,status='unknown',form='unformatted', access='stream')
        rewind(ispf)
        write(ispf) time,Re, alp, bet, mx,my,mz,-(nspec+1),nacumsp
        write(ispf) (jsptot(j), j=1,nspec+1)

        do k=1,8
           write(ispf) (sp(i,0,k),i=0,(mx1+1)*(nz1+1) - 1)
        enddo

     else
     
        do k=1,8
           write(ispf) (sp(i,0,k),i=0,(mx1+1)*(nz1+1) - 1)
        enddo

     endif


   endif
end subroutine io_write_spec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine inoutphys takes a variable transposed to x-z planes!
! transforms it to physical space                               !
! and saves in a file filout.ext1.ext2                          !
! careful, transpose data if they are in kb:ke!                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_phys_field(varspec,wk1c,wk1r,ext2)

  use running
  use ctes
  use mpi_var

  implicit none
  include "mpif.h"
  include "ctes3D"

  ! ----------------------- workspaces  -------------------------------!

  integer i,k,j
  complex(4), dimension(0:mx1,0:mz1,jb:je):: varspec ! var in spectral space
  complex(4), dimension(0:mx1,0:mz1):: wk1c ! work arrays
  real(4), dimension(mgalx+2 ,mgalz):: wk1r
  real(4), dimension(mgalx+2,mgalz,jb:je):: varphys ! var in phys space to print to file
  integer leng ! length of buffer for MPI sends receives leng = (jend(iproc)-jbeg(iproc)+1)*(mgalx+2)*mgalz
  character*200 flnm
  character*4 ext1, ext2

  !------------------------------------------------------------------!
  !---Fourier transforms of the planes------!
  !---IMPORTANT!!! fourxz expects the same buffer as input and output parameter
  !---When calling subroutine, pass the same buffer to both wk1 AND wk1r
  !---They have the same memory address and
  !---If one changed, the second is changed too


 do j = jb,je    ! nvecy planes
     wk1c = varspec(:,:,j)
     call fourxz(wk1c,wk1c,1,1)    !
     varphys(:,:,j)=wk1r
 end do

!  ------- collect data in the master, and write it to file
  if(myid==0) then

     write(ext1,'(i4.4)') id22

     flnm = trim(filout)//'.'//trim(ext1)//'.'//trim(ext2)
     write(*,*) 'filename', flnm, 'ext1', ext1
     open (iout,file=flnm,status='unknown',access="stream")
     rewind(iout)
     write(iout) time,Re,alp,bet,a0,mgalx,my,mgalz
     do j=jb,je
        write(iout) ((varphys(i,k,j),i=1,mgalx),k=1,mgalz)
     enddo
     do iproc=1,numerop-1
        leng=(jend(iproc)-jbeg(iproc)+1)*(mgalx+2)*mgalz
        call MPI_RECV(varphys,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        do j=jb,jb+jend(iproc)-jbeg(iproc)
           write(iout) ((varphys(i,k,j),i=1,mgalx),k=1,mgalz)
        enddo
     enddo
     close(iout)
  else
     call MPI_SEND(varphys,(mgalx+2)*mgalz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
  endif

endsubroutine io_phys_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine inoutspec takes a variable transposed to x-z planes!
! saves in in spectral space as given by the code               !
! and saves in a file filout.ext1.ext2                          !
! careful, transpose data if they are in kb:ke!                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_spec_field(varspec,ext2)

  use running
  use ctes
  use mpi_var

  implicit none
  include "mpif.h"
  include "ctes3D"

  ! ----------------------- workspaces  -------------------------------!

  integer i,k,j
  complex(4), dimension(0:mx1,0:mz1,jb:je):: varspec, varaux ! var in spectral space
  integer leng ! length of buffer for MPI sends receives leng = (jend(iproc)-jbeg(iproc)+1)*(mgalx+2)*mgalz
  character*200 flnm
  character*4 ext1, ext2

 do j = jb,je    ! nvecy planes
     varaux(:,:,j)=varspec(:,:,j)
 end do

!  ------- collect data in the master, and write it to file
  if(myid==0) then

     write(ext1,'(i4.4)') id22

     flnm = trim(filout)//'.'//trim(ext1)//'.'//trim(ext2)
     write(*,*) 'filename', flnm, 'ext1', ext1
     open (iout,file=flnm,status='unknown',access="stream")
     rewind(iout)
     write(iout) time,Re,alp,bet,a0,mx1+1,my,mz1+1
     do j=jb,je
        write(iout) ((varaux(i,k,j),i=0,mx1),k=0,mz1)
     enddo
     do iproc=1,numerop-1
        leng=(jend(iproc)-jbeg(iproc)+1)*(mx1+1)*(mz1+1)
        call MPI_RECV(varaux,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
        do j=jb,jb+jend(iproc)-jbeg(iproc)
           write(iout) ((varaux(i,k,j),i=0,mx1),k=0,mz1)
        enddo
     enddo
     close(iout)
  else
     call MPI_SEND(varaux,(mx1+1)*(mz1+1)*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
  endif

endsubroutine io_spec_field

subroutine io_write_header()
  use mpi_var
  use boundary
  use running
  use ctes

  implicit none
  include 'mpif.h'
  include 'ctes3D'

  if(myid.eq.0) then
     write(*,*) '**********************************************************'
     write(*,*) 'Reynolds number and box size'
     write(*,'(a10,f8.2)') 'Re =',Re
     write(*,'(a10,f8.3,a10,f6.3)') 'alp =',alp,'  bet =',bet
     write(*,*) 'Numerical resolution'
     write(*,'(a10,i5,a10,i5,a10,i5)') 'mgalx =',mgalx,'mgalz =',mgalz,'my =',my
     write(*,'(a10,i5,a10,i5,a10,i5)') 'mx =',mx,'mz =',mz
     write(*,*) 'Timestep and file writing'
     write(*,'(a10,i6,a10,i6,a10,i5)') 'nstep =',nstep,'nimag =',nimag,'nhist =',nhist
     write(*,'(a10,e12.4,a10,f5.2)') 'Delt =',Delt,'  CFL =',CFL
     write(*,*) 'Disturbances at the boundary'
     write(*,'(a10,3f8.3,a10,2f8.3)') 'opposition',alp_u,alp_v,alp_w,'phase',t_lag, x_lag
     write(*,'(a10,2i4,a10,2i4)') 'mxwall',mxmin, mxmax,'mzwall',mzmin, mzmax
     write(*,*) 'Input/output'
     write(*,'(a,a)') '    reading from:  ',filinp
     write(*,'(a,a)') '        write in:  ',filout
     write(*,*) '**********************************************************'
  endif
end subroutine io_write_header


subroutine io_write_footer()
  use ctes
  use running
  use mpi_var

  implicit none
  include "mpif.h"

  if (myid.eq.0) then
    print *,"Total time: ",totaltimer
    print *,"Trans. time: ",transtimer
    print *,"Comm. time: ",commtimer
    print *,"Comm/Total: ",commtimer/totaltimer
    print *,"Trans/Total: ",transtimer/totaltimer
    print *,"Aver time per step: ",totaltimer/nstep
  end if

end subroutine io_write_footer


subroutine io_write_stat()
use statistics
use mpi_var
use running
use ctes

implicit none
include "ctes3D"
character*4 ext1
character*200 fnamesta

real*4 timed,Reed,alped,beted,a0ed
real*8 fac
integer j
  !c                                 /*       write statistics       */
if (nstart.ne.0.and.nacum.ne.0) then

    if (myid==0) then

        !write(*,*) 'stat esc', nstart,nacum
        timed = time
        Reed  = Re
        alped = alp
        beted = bet
        a0ed  = a0

        fac = 1./dble(mgalx*mgalz)

        write(ext1,'(i4.4)') id22
        fnamesta=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.sta'
        open (isn,file=fnamesta,status='unknown',form='unformatted', access='stream')
        rewind(isn)
        write(isn) nacum,Wx0a/nacum,Wz0a/nacum
        write(isn) my,timed,Reed,alped,beted,a0ed
        write(isn) (um(j), vm(j), wm(j), up(j), vp(j), wp(j), &
                 w1m(j),w2m(j),w3m(j),w1p(j),w2p(j),w3p(j),uvr(j),uwr(j), &
                 vwr(j),ep(j),uuv(j)*fac,wwv(j)*fac,vvv(j)*fac,j=1,my)
        write(isn) (y(j),j=1,my)
        write(isn) (fmap(j),j=1,my)
    endif
endif

end subroutine io_write_stat

subroutine io_calc_wtime(write_time)
  use running
  use mpi_var
  implicit none
  include "mpif.h"
  real*8 write_time

  if (myid.eq.0) then
      write_time = -MPI_WTIME()
      write(*,*) "time is ", time
  endif

end subroutine io_calc_wtime

subroutine io_calc_timers_io(istep)
use mpi_var
use running
use ctes
use statistics

implicit none
include "mpif.h"
real*8 iter_time
integer istep

if (myid.eq.0) then
  	totaltimer = totaltimer-MPI_WTIME()
  	iter_time=-MPI_WTIME()
endif

if (mod(istep-1,nhist) .eq. 0) then
    ihist=1
    icfl= 1
endif

if (mod(istep-1,ntimes) .eq.0 .and.nstart.ne.0) then
    istati=1
endif

end subroutine io_calc_timers_io

subroutine io_write_cf(write_time)
  use running
  use mpi_var

  implicit none
  include "mpif.h"
  real*8 write_time
  character*4 ext1
  character*200 fname

  if (myid.eq.0) then
      write(*,*) 'time write:',MPI_WTIME()+write_time
      id22 = id22+1
      close(39)
      write(ext1,'(i4.4)') id22
      fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
      write (*,*) fname
      open(39,file=fname,status='unknown')
  endif

end subroutine io_write_cf



