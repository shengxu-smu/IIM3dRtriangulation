!----------------------------------------------------------------------
subroutine run_input
  !---------------------
  use para
  !--------------------------
  namelist /domainpara/ x0,y0,z0,xl,yl,zl,nx,ny,nz,tstart,tend,&
       epsilonx,epsilony,epsilonz,eps
  !-----------------------------------------------
  open(unit=8,file='DAT/input.run',status='old')
  rewind 8
  read(8,*)
  read(8,*)nstep,icfl,iread,iwrite,itout,iplot,ianimation
  read(8,*)
  read(8,*)
  read(8,*)cflc,cflv,dtcfl,dtfix,tout
  close(8)
  !-----------------------------------------------------
  open(88,file='./DAT/domainpara.txt')
  read(88,domainpara)
  close(88)
  !---------------------------------------------------
  
  open(unit=61,file='DAT/fbyfs.dat',status='unknown')
  open(unit=62,file='DAT/fbypp.dat',status='unknown')
  open(unit=63,file='DAT/xst.dat',status='unknown')
  open(unit=64,file='DAT/yst.dat',status='unknown')
  open(unit=65,file='DAT/zst.dat',status='unknown')
  open(unit=66,file='DAT/qt.dat',status='unknown')
  !------------------------------------------------------
  return

end subroutine run_input

!-----------------------------------------------------------------------
