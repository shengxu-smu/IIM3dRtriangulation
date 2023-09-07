subroutine allocateEuler


  use para
  use Euler

  ! u,v,w
  allocate(u(0:nx,0:ny+1,0:nz+1))
  allocate(v(0:nx+1,0:ny,0:nz+1))
  allocate(w(0:nx+1,0:ny+1,0:nz))
  ! p
  allocate(p(0:nx+1,0:ny+1,0:nz+1))
  
  ! io indicator
  allocate(iop(0:nx+1,0:ny+1,0:nz+1))
  allocate(iou(0:nx+1,0:ny+1,0:nz+1))
  allocate(iov(0:nx+1,0:ny+1,0:nz+1))
  allocate(iow(0:nx+1,0:ny+1,0:nz+1))

end subroutine allocateEuler
