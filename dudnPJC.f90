subroutine dudnPJC(s0,u_s0,n,ds,dudn)
  !---------------------------------------------------------------------------------
  use para
  use field
  use Lagrange
  !---------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------
  !interface
  !----------------------------------------------------------------------------------
  ! declaration:
  double precision, dimension(3), intent(in) :: s0,u_s0,n
  double precision, intent(in) :: ds
  double precision, dimension(3), intent(out) :: dudn
  
  integer, parameter:: nn=2,nijk=3
  double precision :: uu(nn),vv(nn),ww(nn),hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
  integer :: id,jd,kd,inn,ic,ie,jc,je,kc,ke,iu,ju,ku,iv,jv,kv,iw,jw,kw
  double precision, dimension(3) :: s,foo
  !----------------------------------------------------------------------------------
  dudn = 0.0d0
  id = int(sign(1.0d0,n(1)))    ! x-direction                
  jd = int(sign(1.0d0,n(2)))    ! y-direction            
  kd = int(sign(1.0d0,n(3)))    ! z-direction
  !---------------------------------------------------------------------------------
  do inn=1,nn
     !s = s0+ds*n*dble(inn)
     s(1) = s0(1)+ds*n(1)*dble(inn)
     s(2) = s0(2)+ds*n(2)*dble(inn)
     s(3) = s0(3)+ds*n(3)*dble(inn)

     i = int((s(1)-x0)/hdx)  ! index                                                      
     j = int((s(2)-y0)/hdy)
     k = int((s(3)-z0)/hdz)

     ic = int((s(1)-xc(0))*dx1)
     ie = int((s(1)-xe(0))*dx1)
     jc = int((s(2)-yc(0))*dy1)
     je = int((s(2)-ye(0))*dy1)
     kc = int((s(3)-zc(0))*dz1)
     ke = int((s(3)-ze(0))*dz1)

     iu = ie
     ju = jc
     ku = kc
     
     iv = ic
     jv = je
     kv = kc
     
     iw = ic
     jw = jc
     kw = ke

     if(id.lt.0.0d0) then
        iu = iu+1
        iv = iv+1
        iw = iw+1
     endif
     
     if(jd.lt.0.0d0) then
        ju = ju+1
        jv = jv+1
        jw = jw+1
     endif

     if(kd.lt.0.0d0) then
        ku = ku+1
        kv = kv+1
        kw = kw+1
     endif
     !------------------------------------------------------------------------
     !-----------------------u------------------------------------------------
     do k=1,nijk
        do j=1,nijk
           do i=1,nijk 
              hi(i)=xe(iu+id*(i-1))  ! x coordinate
              fi(i)=u(iu+id*(i-1),ju+jd*(j-1),ku+kd*(k-1))
           enddo
           call interpolate(hi,fi,nijk,s(1),fj(j),foo)     
           hj(j)=yc(ju+jd*(j-1))
        enddo
        call interpolate(hj,fj,nijk,s(2),fk(k),foo)
        hk(k)=zc(ku+kd*(k-1))
     enddo
     call interpolate(hk,fk,nijk,s(3),uu(inn),foo)
     !----------------------v-----------------------------------------------       
     do k=1,nijk
        do j=1,nijk
           do i=1,nijk
              hi(i)=xc(iv+id*(i-1))
              fi(i)=v(iv+id*(i-1),jv+jd*(j-1),kv+kd*(k-1))
           enddo
           call interpolate(hi,fi,nijk,s(1),fj(j),foo)
           hj(j)=ye(jv+jd*(j-1))
        enddo
        call interpolate(hj,fj,nijk,s(2),fk(k),foo)
        hk(k)=zc(kv+kd*(k-1))
     enddo
     call interpolate(hk,fk,nijk,s(3),vv(inn),foo)
     !----------------------w------------------------------------------------        
     do k=1,nijk
        do j=1,nijk
           do i=1,nijk
              hi(i)=xc(iw+id*(i-1))
              fi(i)=w(iw+id*(i-1),jw+jd*(j-1),kw+kd*(k-1))
           enddo
           call interpolate(hi,fi,nijk,s(1),fj(j),foo)
           hj(j)=yc(jw+jd*(j-1))
        enddo
        call interpolate(hj,fj,nijk,s(2),fk(k),foo)
        hk(k)=ze(kw+kd*(k-1))
     enddo
     call interpolate(hk,fk,nijk,s(3),ww(inn),foo)
  enddo
  !---------------------------------------------------------------------------------
  if(nn.eq.1) then                                   
     dudn(1) = (uu(1)-u_s0(1))/ds
     dudn(2) = (vv(1)-u_s0(2))/ds
     dudn(3) = (ww(1)-u_s0(3))/ds
  elseif(nn.eq.2) then
     dudn(1) = (-3.0d0*u_s0(1) + 4.0d0*uu(1) - uu(2))/(2.0d0*ds)
     dudn(2) = (-3.0d0*u_s0(2) + 4.0d0*vv(1) - vv(2))/(2.0d0*ds)
     dudn(3) = (-3.0d0*u_s0(3) + 4.0d0*ww(1) - ww(2))/(2.0d0*ds)
  elseif(nn.eq.3) then
     !dudn(1) = (18.0d0*uu(1) - 9.0d0*uu(2) + 2.0d0*uu(3) - 11.0d0*u_s0(1))/(6.0d0*ds)
     !dudn(2) = (18.0d0*vv(1) - 9.0d0*vv(2) + 2.0d0*vv(3) - 11.0d0*u_s0(2))/(6.0d0*ds)
     !dudn(3) = (18.0d0*ww(1) - 9.0d0*ww(2) + 2.0d0*ww(3) - 11.0d0*u_s0(3))/(6.0d0*ds)
  endif

     
end subroutine dudnPJC
!----------------------------------------------------------
