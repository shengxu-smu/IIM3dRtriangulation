subroutine allocate_vfft

  use para
  use vfft
  
  allocate(wsavei(3*nx+15),wsavej(3*ny+15),wsavek(3*nz+15))
  allocate(ri(ny,nx),wi(ny,nx))
  allocate(rj(nx,ny),wj(nx,ny))
  allocate(rk(nx,nz),wk(nx,nz))

end subroutine allocate_vfft
