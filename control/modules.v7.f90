module running
  save
  real*4  Delt,Deltat,dtr
  real*4  time
  integer nimag,nstep,nhist,ihist,istart,icfl, var_dt
  real*4  ener(9),Wx0,Wz0,WxL,WzL, uv0, uvL
  real*8  Wx0a,Wz0a

  !  ------- file things ----
  integer iinp,iout,istt,ihre,id22,isn,ispf
  character*200 filinp,filout,filstt

end module running

module ctes
  private
  include 'ctes3D'
  save
  real*8, public :: commtimer,transtimer,totaltimer
  integer, public:: myslice(0:numerop-1)
  integer, dimension(0:numerop-1), public:: jbeg,jend,kbeg,kend
  integer, public :: jb,je,kb,ke,mmy,mmz
  real*8, dimension(my), public ::  fmap,y2
  real*4, public ::  Re,alp,bet,a0,CFL
  real*4, dimension(my), public::y,hy
  complex*8, public:: xalp(0:mx1), xbet(0:mz1)
  real*4, public ::  alp2(0:mx1),bet2(0:mz1)
  integer, public :: icx(0:mz1), iax(mx)
  integer, dimension(0:numerop-1), public :: jspbeg, jspend
  integer, public  :: jspb,jspe
  integer, public :: ifor
  real*4, public:: gamma
  integer, public:: imesh
  ! precomputed diffenetiation matrices
  real(8),public:: prem1(5,my)
  real(8),public:: dt12(5,my),dt21(5,my),dt22(5,my)
  ! Runge-Kutta
  real*4, public :: c(3)
  parameter(c = (/ 1./3., 1./2., 1. /))
end module ctes


module boundary
  private
  include 'ctes3D'
  save
  ! wavenumbers for control
  integer, public:: mxmin,mzmin, mxmax, mzmax
  ! velocity at the boundaries
  complex*8, dimension(2,0:mx1,0:mz1), public :: u_bc, v_bc, w_bc
  real*4, public :: alp_u, alp_v, alp_w
  real*4, public :: t_lag, x_lag
  integer, public :: nplane_c, id_up, id_down, f_wU,f_lse,dtc,savefields
end module boundary


module statistics
  private
  include 'ctes3D'
  save
  ! rms and mean quantities
  real*8, dimension(my), public :: um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p
  real*8, dimension(my), public:: uvr,uwr,vwr,ep,uuv,wwv,vvv
  integer, public :: istati,ntimes,nacum,nstart
  ! parameters for spectra
  integer, dimension(2*nspec+1), public :: jsptot, jspiproc
  integer, public  :: nacumsp,jsp(my),jspbb
  ! mass flux calculations
  real*4, public:: trp(0:my1)
  real*8, public :: trp2(0:my1)
end module statistics

module planes
  private
  include 'ctes3D'
  save
  real*4, dimension(mgalx+2,mgalz), public :: u1r,u2r,u3r,o1r,o2r,o3r
end module planes


module mpi_var
  private
  include 'mpif.h'
  save
  integer, public :: istat(MPI_STATUS_SIZE),ierr
  integer, public :: master,myid,iproc,numprocs
end module mpi_var
