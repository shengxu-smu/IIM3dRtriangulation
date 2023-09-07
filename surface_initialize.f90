!----------------------------------------------------------------------------
!     the refrence:Xu,S., Singular forces in the immersed interface method
!      for rigid objects in 3D, Applied Mathematics Letters 22 (2009) 827833
! ---------------------------------------------------------------------------
subroutine surface_initialize
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  use para
  use Lagrange
  use field
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=
  integer :: m,iv
  double precision:: s1,c1,s2,c2,s3,c3
  double precision:: xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO m=1,nobj
     !-----------------------------------------------------------
     !the motion of the refrence point(fixed with respect to S)
     !Xsc = (xsc,ysc,zsc) - the new position
     !Xsc0 = (xsc0,ysc0,zsc0)-the position at some refrence time t 
     xsc(m)=xsc0(m)
     ysc(m)=ysc0(m) 
     zsc(m)=zsc0(m)   
     ! velocity at Xsc:  d(Xsc)/dt
     xsct(m)=0.0d0
     ysct(m)=0.0d0  
     zsct(m)=0.0d0
     !accelaration  d2(Xsc)/dtdt
     xsctt(m)=0.0d0
     ysctt(m)=0.0d0  
     zsctt(m)=0.0d0   
     !the Euler angles
     phi(m)=phi0(m)
     theta(m)=theta0(m)  
     psi(m)=psi0(m)
     !
     phit(m)=0.0d0
     thetat(m)=0.0d0  
     psit(m)=0.0d0
     !         
     phitt(m)=0.0d0
     thetatt(m)=0.0d0 
     psitt(m)=0.0d0
     !----------------   
     s1=dsin(phi(m))
     c1=dcos(phi(m))
     
     s2=dsin(theta(m))
     c2=dcos(theta(m))
     
     s3=dsin(psi(m))
     c3=dcos(psi(m))
     !----------------------------------------
     ! the angular velocity Omega=(omegax,omegay,omegaz),which is given in equation (14)
     omegax(m)=s3*s2*phit(m)+c3*thetat(m)
     omegay(m)=s3*thetat(m)-c3*s2*phit(m)
     omegaz(m)=c2*phit(m)+psit(m)
     ! The angular acceleration
     omegaxt(m)= s3*s2*phitt(m)+s3*c2*thetat(m)*phit(m)+&
                 c3*psit(m)*s2*phit(m)+&
                 c3*thetatt(m)-s3*psit(m)*thetat(m)
     omegayt(m)= s3*thetatt(m)+c3*psit(m)*thetat(m)&
                 -c3*s2*phitt(m)-c3*c2*thetat(m)*phit(m)&
                 +s3*psit(m)*phit(m)
     omegazt(m)=c2*phitt(m)-s2*thetat(m)*phit(m)+psitt(m)
     !----------------------------------------------------
     !------------------------------------------------------------------
     do iv = nv4obj(m-1)+1,nv4obj(m)           
        !compute  rhs of equation 11: multiply the rotation matrix R by
        !the difference between x0-xc0 
                      
        !XX1 = R1*(x0-xc0)        ** !xc0=0
        xx1=c1*vertex(1,iv)-s1*vertex(2,iv)
        yy1=s1*vertex(1,iv)+c1*vertex(2,iv)
        zz1=vertex(3,iv)
        ! XX2 = R2*XX1
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
     enddo
     !----------------------------------------------------------
  ENDDO
      
  return
end subroutine surface_initialize
!-----------------------------------------------------------------------
