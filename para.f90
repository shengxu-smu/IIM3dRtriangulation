module para
  implicit none
  !-----------------
  !constants                       
  double precision, parameter:: pi=3.14159265358979323d0
  double precision, parameter:: e =2.71828182845904524d0
  double precision, parameter:: ZERO=0.0d0
  double precision, parameter:: HALF=0.5d0
  double precision, parameter:: ONE=1.0d0
  double precision, parameter:: TWO=2.0d0
  double precision, parameter:: eps1=1.0d-8
  double precision, parameter:: DPeps=1.0d-15
  double precision, parameter:: Re=10.0d0,Re1=1.0d0/Re
  !-----------------
  integer :: i,j,k
  integer :: nstart,nend
  integer :: nx,ny,nz
  integer :: nstep,icfl,iread,iwrite,itout,iplot,ianimation
  double precision :: cflc,cflv,dtcfl,dtfix,tout,t,t0,dt,dt1,dt2

  double precision :: x0,y0,z0,xl,yl,zl,dx,dy,dz,hdx,hdy,hdz
  double precision :: dx1,dx2,dy1,dy2,dz1,dz2,betay,betaz
  double precision :: epsilonx,epsilony,epsilonz,eps
  !-------------------------------
  integer, parameter :: isingular=1,imove=0 
  integer, parameter :: lw_ubc=1,le_ubc=2,ls_ubc=1,ln_ubc=1,lb_ubc=1,lt_ubc=1
  integer, parameter :: lw_vbc=1,le_vbc=2,ls_vbc=1,ln_vbc=1,lb_vbc=1,lt_vbc=1
  integer, parameter :: lw_wbc=1,le_wbc=2,ls_wbc=1,ln_wbc=1,lb_wbc=1,lt_wbc=1
  integer, parameter :: lw_pbc=2,le_pbc=2,ls_pbc=2,ln_pbc=2,lb_pbc=2,lt_pbc=2
  
end module para

