!-----------------------------------------------------------------------

subroutine animation
  use para
  use field
  use Lagrange
  
  integer :: m,ms,m2
  ms =1

  call vorticity

  DO m=1,ms
     write(63,200) (xs(m2,m),m2=nv4obj(m-1)+1,nv4obj(m)) 
     write(64,200) (ys(m2,m),m2=nv4obj(m-1)+1,nv4obj(m))
     write(65,200) (zs(m2,m),m2=nv4obj(m-1)+1,nv4obj(m))
     
  ENDDO

  do k=1,nz
     do j=1,ny
        write(66,200)(o(i,j,k),i=1,nx)
     enddo
  enddo

200 format(1x,1000e16.6e4)

  return
end subroutine animation


!c-----------------------------------------------------------------------
