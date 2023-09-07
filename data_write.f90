subroutine data_write

  use para
  use Lagrange
  use field

  !-------------------
  open(unit=12,file='DAT/new.out',form='unformatted',status='unknown')
  
  rewind 12
  write(12) nend,t,t0
  write(12) xsc,ysc,zsc,phi,theta,psi
  write(12) xs,ys,zs,us,vs,ws
  write(12) x,y,z,xe,ye,ze,xc,yc,zc
  write(12) u,v,w,p
  close(12)
  return
  
  

end subroutine data_write
