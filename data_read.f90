subroutine data_read
  !----------------------------------------------------------------
  use para
  use Lagrange
  use field
  !----------------------------------------------------------------
  open(unit=11,file='DAT/old.out',form='unformatted',status='old')
  rewind 11
  read(11) nstart,t,t0
  read(11) xsc,ysc,zsc,phi,theta,psi
  read(11) xs,ys,zs,us,vs,ws
  read(11) x,y,z,xe,ye,ze,xc,yc,zc
  read(11) u,v,w,p
  close(11)
  return
  !---------------------------------------------------------------
end subroutine data_read
