subroutine allocate_field

  use para
  use field

  ! grid
  !allocate(x(-2:2*nx+1),y(-2:2*ny+1),z(-2:2*nz+1))
  allocate(x(-2:2*nx),y(-2:2*ny),z(-2:2*nz))
  allocate(xe(0:nx),ye(0:ny),ze(0:nz))
  allocate(xc(0:nx+1),yc(0:ny+1),zc(0:nz+1))
  ! velocity,pressure,d, and o.... (d:divergence of the velocity,o:vorticity)
  allocate(u(0:nx,0:ny+1,0:nz+1))
  allocate(v(0:nx+1,0:ny,0:nz+1))
  allocate(w(0:nx+1,0:ny+1,0:nz))
  allocate(p(0:nx+1,0:ny+1,0:nz+1))
  allocate(d(0:nx+1,0:ny+1,0:nz+1))
  allocate(o(0:nx+1,0:ny+1,0:nz+1))
  ! 
  allocate(un(0:nx,0:ny+1,0:nz+1))
  allocate(vn(0:nx+1,0:ny,0:nz+1))
  allocate(wn(0:nx+1,0:ny+1,0:nz))
  allocate(dn(0:nx+1,0:ny+1,0:nz+1))

  allocate(um(0:nx,0:ny+1,0:nz+1))
  allocate(vm(0:nx+1,0:ny,0:nz+1))
  allocate(wm(0:nx+1,0:ny+1,0:nz))

  !------
  allocate(fcee(0:nx+1,0:ny+1,0:nz+1))
  allocate(fcec(0:nx+1,0:ny+1,0:nz+1))
  allocate(fece(0:nx+1,0:ny+1,0:nz+1))
  allocate(fecc(0:nx+1,0:ny+1,0:nz+1))
  allocate(fcce(0:nx+1,0:ny+1,0:nz+1))
  allocate(fccc(0:nx+1,0:ny+1,0:nz+1))
  allocate(feec(0:nx+1,0:ny+1,0:nz+1))
  !-------
  allocate(data1(0:nx+1,0:ny+1,0:nz+1))
  allocate(data2(0:nx+1,0:ny+1,0:nz+1))
  allocate(data3(0:nx+1,0:ny+1,0:nz+1))
  allocate(data4(0:nx+1,0:ny+1,0:nz+1))
  allocate(data5(0:nx+1,0:ny+1,0:nz+1))
  allocate(data6(0:nx+1,0:ny+1,0:nz+1))
  allocate(data7(0:nx+1,0:ny+1,0:nz+1))
  allocate(data8(0:nx+1,0:ny+1,0:nz+1))
  allocate(data9(0:nx+1,0:ny+1,0:nz+1))
  !-------
  ! io indicator                                                              
  allocate(iop(0:nx+1,0:ny+1,0:nz+1))
  allocate(iou(0:nx+1,0:ny+1,0:nz+1))
  allocate(iov(0:nx+1,0:ny+1,0:nz+1))
  allocate(iow(0:nx+1,0:ny+1,0:nz+1))
  !--------------------------------
  allocate(pw(ny,nz),pe(ny,nz))
  allocate(ps(nx,nz),pn(nx,nz))
  allocate(pb(nx,ny),pt(nx,ny))
  !----------  
end subroutine allocate_field
