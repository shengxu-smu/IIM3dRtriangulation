! this subroutine is to discretize the domain using the MAC scheme
! Modified: NOV. 3rd,2021
!-----------------------------------------------------------------------
subroutine mesh
  
  use para
  use field

  dx=(xl-x0)/dble(nx)
  hdx=0.5d0*dx
  dx1=1.0d0/dx
  dx2=dx1*dx1

  dy=(yl-y0)/dble(ny)
  hdy=0.5d0*dy
  dy1=1.0d0/dy
  dy2=dy1*dy1

  dz=(zl-z0)/dble(nz)
  hdz=0.5d0*dz
  dz1=1.0d0/dz
  dz2=dz1*dz1
  
  betay=dx*dx*dy2
  betaz=dx*dx*dz2
 
  do i=-2,2*nx
     x(i)=x0+hdx*dble(i)
  enddo
  do j=-2,2*ny
     y(j)=y0+hdy*dble(j)
  enddo
  do k=-2,2*nz
     z(k)=z0+hdz*dble(k)
  enddo
  
  xc(0)=x(-2)
  do i=0,nx
     xc(i+1)=x(2*i)
     xe(i)=xc(i)+hdx
  enddo
  yc(0)=y(-2)
  do j=0,ny
     yc(j+1)=y(2*j)
     ye(j)=yc(j)+hdy
  enddo
  zc(0)=z(-2)
  do k=0,nz
     zc(k+1)=z(2*k)
     ze(k)=zc(k)+hdz
  enddo
  
  return
end subroutine mesh
!-----------------------------------------------------------------------
