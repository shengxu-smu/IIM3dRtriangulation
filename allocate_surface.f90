subroutine allocate_surface

  use Lagrange

  
  allocate(xsc(nobj),ysc(nobj),zsc(nobj)) ! the center of mass at time t
  allocate(xsc0(nobj),ysc0(nobj),zsc0(nobj)) ! center of mass at refrence time t0  
  allocate(xsct(nobj),ysct(nobj),zsct(nobj)) ! velocity,
  allocate(xsctt(nobj),ysctt(nobj),zsctt(nobj))  !acceleration
  allocate(phi(nobj),theta(nobj),psi(nobj))  !Euler angles
  allocate(phi0(nobj),theta0(nobj),psi0(nobj)) !Euler angles - at t0
  allocate(phit(nobj),thetat(nobj),psit(nobj)) 
  allocate(phitt(nobj),thetatt(nobj),psitt(nobj))
  allocate(omegax(nobj),omegay(nobj),omegaz(nobj))   !angular velocity
  allocate(omegaxt(nobj),omegayt(nobj),omegazt(nobj))  !angular acceleration

  ! the coordinates of the sufrace points (the tirangular panel's vertices 
  allocate(xs(nvertices,nobj))
  allocate(ys(nvertices,nobj))
  allocate(zs(nvertices,nobj))

  allocate(xss(nvertices,nobj))
  allocate(yss(nvertices,nobj))
  allocate(zss(nvertices,nobj))
  ! the velocity at vertices on the surface 
  allocate(us(nvertices,nobj))
  allocate(vs(nvertices,nobj))
  allocate(ws(nvertices,nobj))

  allocate(pJC0(nvertices))

  allocate(dudnp_T(3,npanels))
  allocate(dvdnp_T(3,npanels))
  allocate(dwdnp_T(3,npanels))
  allocate(dpdnPJC_T(3,npanels))

end subroutine allocate_surface
