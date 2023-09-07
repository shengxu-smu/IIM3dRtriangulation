!-----------------------------------------------------------------------
subroutine surface_move
  
  use para
  use Lagrange
  
  double precision :: si,ci,sii,cii,siii,ciii,s1,c1,s2,c2,s3,c3
  double precision :: xt,yt,zt,xb,yb,zb,xn,yn,zn,gacobi
  double precision :: xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3
  double precision :: xx(5),yy(5),zz(5)
  double precision :: x1(5),y1(5),z1(5),x2(5),y2(5),z2(5),x3(5),y3(5),z3(5)
  !-----------------
  integer :: n,m,ms
  !ms =1
  !------------------
  t0=t0+0.5d0*dt

  DO m=1,nobj
     
     si=dsin(phi(m))
     ci=dcos(phi(m))

     sii=dsin(theta(m))
     cii=dcos(theta(m))

     siii=dsin(psi(m))
     ciii=dcos(psi(m))

     xsc(m)=xsc0(m)
     ysc(m)=ysc0(m)
     zsc(m)=zsc0(m)

     xsct(m)=0.0d0
     ysct(m)=0.0d0
     zsct(m)=0.0d0

     xsctt(m)=0.0d0
     ysctt(m)=0.0d0
     zsctt(m)=0.0d0

     phi(m)=phi0(m)
     theta(m)=theta0(m)
     psi(m)=psi0(m)


     phit(m)=0.0d0
     thetat(m)=0.0d0
     psit(m)=0.0d0

     phitt(m)=0.0d0
     thetatt(m)=0.0d0
     psitt(m)=0.0d0

     s1=dsin(phi(m))
     c1=dcos(phi(m))

     s2=dsin(theta(m))
     c2=dcos(theta(m))

     s3=dsin(psi(m))
     c3=dcos(psi(m))

     omegax(m)=s3*s2*phit(m)+c3*thetat(m)
     omegay(m)=s3*thetat(m)-c3*s2*phit(m)
     omegaz(m)=c2*phit(m)+psit(m)

     omegaxt(m)=s3*s2*phitt(m)+s3*c2*thetat(m)*phit(m)+  &
                c3*psit(m)*s2*phit(m)+    &
                c3*thetatt(m)-s3*psit(m)*thetat(m)
     omegayt(m)=s3*thetatt(m)+c3*psit(m)*thetat(m)  &
                -c3*s2*phitt(m)-c3*c2*thetat(m)*phit(m)  &
                +s3*psit(m)*phit(m)
     omegazt(m)=c2*phitt(m)-s2*thetat(m)*phit(m)+psitt(m)
     !-------------------------------------- 
     do iv = nv4obj(m-1)+1,nv4obj(m)

        xx1=c1*vertex(1,iv)-s1*vertex(2,iv)
        yy1=s1*vertex(1,iv)+c1*vertex(2,iv)
        zz1=vertex(3,iv)

        xx2=xx1
        yy2=c2*yy1-s2*zz1
        zz2=s2*yy1+c2*zz1

        !XX3 = R3*XX2                                                           
        xx3=c3*xx2-s3*yy2
        yy3=s3*xx2+c3*yy2
	zz3=zz2
        ! Xs = Xsc + XX3                                                        
        xs(iv,m)=xsc(m)+xx3
        ys(iv,m)=ysc(m)+yy3
        zs(iv,m)=zsc(m)+zz3

        !Us = d(Xsc)/dt + Omega x (XX3)                                         
        us(iv,m)=xsct(m)+omegay(m)*zz3-omegaz(m)*yy3
        vs(iv,m)=ysct(m)+omegaz(m)*xx3-omegax(m)*zz3
        ws(iv,m)=zsct(m)+omegax(m)*yy3-omegay(m)*xx3

        vertex(1,iv) = xs(iv,m)
        vertex(2,iv) = ys(iv,m)
        vertex(3,iv) = zs(iv,m)
     enddo
  enddo
       
  return
end subroutine surface_move
!-----------------------------------------------------------------------
