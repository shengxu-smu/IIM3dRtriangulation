!----------------------------------------------------------------------------------------------------
!This subroutine computes the force of the fluid acting on the panel's vertices.
! For a vertex S0 at the surface, we have:
!           F = p|+ + (1/Re)*du/dn|+  -------> 1
!  where p|+ = (1/2)*(p(S1)+p(S(-1) +[p]+[dp/dn]*ds) ----->2
!    and du/dn|+ is known
!    The point S1 and S2 are ds away from S0 in the normal direction outside and inside the object
!    p @ S1 and S2 is computed using the trilinear interpolation
!    [p] and [dp/dn] have been computed in the subroutine JC4uvwp
!-----------------------------------------------------------------
subroutine Force4Panel(A,B,C,n,dudn_A,dudn_B,dudn_C,dpdn_A,dpdn_B,&
     dpdn_C,pJ_A,pJ_B,pJ_C,ds,T_fx,T_fy,T_fz)

  !-------------
  use para
  use field
  use Lagrange

  implicit none
  !------------------
  double precision,dimension(3), intent(in):: A,B,C,n,dudn_A,dudn_B,dudn_C
  double precision, intent(in):: dpdn_A,dpdn_B,dpdn_C,ds,pJ_A,pJ_B,pJ_C
  double precision, dimension(3), intent(out)::T_fx,T_fy,T_fz
  double precision::Fx_A,Fy_A,Fz_A
  double precision::Fx_B,Fy_B,Fz_B
  double precision::Fx_C,Fy_C,Fz_C  
  double precision, dimension(3):: Sp,Sm,foo,nm ! S @ surface point + ds, and S @ surface point minus ds
  integer :: id,jd,kd,ic,jc,kc
  integer :: idm,jdm,kdm
  double precision :: P_Sp,P_Sm
  integer, parameter:: nijk=3
  double precision :: hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
  !-----------------------------------------------------------------------------
  !intialize the force to be zero
  Fx_A = 0.0d0
  Fy_A = 0.0d0
  Fz_A = 0.0d0

  Fx_B = 0.0d0
  Fy_B = 0.0d0
  Fz_B = 0.0d0
  
  Fx_C = 0.0d0
  Fy_C = 0.0d0
  Fz_C = 0.0d0
  !----------------------------------------------------------------------------
  ! direction
  id = int(sign(1.0d0,n(1))) 
  jd = int(sign(1.0d0,n(2)))
  kd = int(sign(1.0d0,n(3)))
  nm=-n
  idm = int(sign(1.0d0,nm(1)))
  jdm = int(sign(1.0d0,nm(2)))
  kdm = int(sign(1.0d0,nm(3)))
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------------- 
  ! For a point A:
  ! 1. compute p(S1) 
  ! 1.1 compute S1 
  Sp(1) = A(1)+ds*n(1)
  Sp(2) = A(2)+ds*n(2)
  Sp(3) = A(3)+ds*n(3)
  ! 1.2 indices
  ic = int((Sp(1)-xc(0))*dx1)
  jc = int((Sp(2)-yc(0))*dy1)
  kc = int((Sp(3)-zc(0))*dz1)
  if(id.lt.0.0d0) then
     ic = ic+1
  endif
  if(jd.lt.0.0d0) then
     jc = jc+1
  endif
  if(kd.lt.0.0d0) then
     kc = kc+1
  endif
  ! 1.3 interpolate
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+id*(i-1))  ! x coordinate
           fi(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sp(1),fj(j),foo)
        hj(j)=yc(jc+jd*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sp(2),fk(k),foo)
     hk(k)=zc(kc+kd*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sp(3),P_Sp,foo)
  !------------------------------------
  ! 2. compute p(S-1)
  ! 2.1 compute S-1     
  Sm(1) = A(1)+ds*nm(1)
  Sm(2) = A(2)+ds*nm(2)
  Sm(3) = A(3)+ds*nm(3)
  ! 2.2 indices 
  ic = int((Sm(1)-xc(0))*dx1)
  jc = int((Sm(2)-yc(0))*dy1)
  kc = int((Sm(3)-zc(0))*dz1)
  if(idm.lt.0.0d0) then
     ic = ic+1
  endif
  if(jdm.lt.0.0d0) then
     jc = jc+1
  endif
  if(kdm.lt.0.0d0) then
     kc = kc+1
  endif
  ! 2.3 interpolate                                                                                       
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+idm*(i-1))  ! x coordinate                                                                                                      
           fi(i)=p(ic+idm*(i-1),jc+jdm*(j-1),kc+kdm*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sm(1),fj(j),foo)
        hj(j)=yc(jc+jdm*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sm(2),fk(k),foo)
     hk(k)=zc(kc+kdm*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sm(3),P_Sm,foo)
  !---------------------------------------------------------
  Fx_A = -0.50d0*(P_Sp+P_Sm+pJ_A-dpdn_A*ds)*n(1)+Re1*dudn_A(1)
  Fy_A = -0.50d0*(P_Sp+P_Sm+pJ_A-dpdn_A*ds)*n(2)+Re1*dudn_A(2)
  Fz_A = -0.50d0*(P_Sp+P_Sm+pJ_A-dpdn_A*ds)*n(3)+Re1*dudn_A(3)
  !---------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------------- 
  ! @ B
  Sp(1) = B(1)+ds*n(1)
  Sp(2) = B(2)+ds*n(2)
  Sp(3) = B(3)+ds*n(3)
  ! 1.2 indices                                                                                                                                       
  ic = int((Sp(1)-xc(0))*dx1)
  jc = int((Sp(2)-yc(0))*dy1)
  kc = int((Sp(3)-zc(0))*dz1)
  if(id.lt.0.0d0) then
     ic = ic+1
  endif
  if(jd.lt.0.0d0) then
     jc = jc+1
  endif
  if(kd.lt.0.0d0) then
     kc = kc+1
  endif
  ! 1.3 interpolate                                                                                                                                   
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+id*(i-1))  ! x coordinate                                                                                                      
           fi(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sp(1),fj(j),foo)
        hj(j)=yc(jc+jd*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sp(2),fk(k),foo)
     hk(k)=zc(kc+kd*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sp(3),P_Sp,foo)
  ! 2.1 compute S-1                                                                                                                                   
  Sm(1) = B(1)+ds*nm(1)
  Sm(2) = B(2)+ds*nm(2)
  Sm(3) = B(3)+ds*nm(3)
  ! 2.2 indices                                                                                                                                       
  ic = int((Sm(1)-xc(0))*dx1)
  jc = int((Sm(2)-yc(0))*dy1)
  kc = int((Sm(3)-zc(0))*dz1)
  if(idm.lt.0.0d0) then
     ic = ic+1
  endif
  if(jdm.lt.0.0d0) then
     jc = jc+1
  endif
  if(kdm.lt.0.0d0) then
     kc = kc+1
  endif
  ! 2.3 interpolate                                                                                                                                   
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+idm*(i-1))  ! x coordinate                                                                                                      
           fi(i)=p(ic+idm*(i-1),jc+jdm*(j-1),kc+kdm*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sm(1),fj(j),foo)
        hj(j)=yc(jc+jdm*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sm(2),fk(k),foo)
     hk(k)=zc(kc+kdm*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sm(3),P_Sm,foo)
  !---------------------------------------------------------                                                                                          
  Fx_B = -0.50d0*(P_Sp+P_Sm+pJ_B-dpdn_B*ds)*n(1)+Re1*dudn_B(1)
  Fy_B = -0.50d0*(P_Sp+P_Sm+pJ_B-dpdn_B*ds)*n(2)+Re1*dudn_B(2)
  Fz_B = -0.50d0*(P_Sp+P_Sm+pJ_B-dpdn_B*ds)*n(3)+Re1*dudn_B(3)
  !-----------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------------- 
  ! @ C
  Sp(1) = C(1)+ds*n(1)
  Sp(2) = C(2)+ds*n(2)
  Sp(3) = C(3)+ds*n(3)
  ! 1.2 indices                                                                                                                                       
  ic = int((Sp(1)-xc(0))*dx1)
  jc = int((Sp(2)-yc(0))*dy1)
  kc = int((Sp(3)-zc(0))*dz1)
  if(id.lt.0.0d0) then
     ic = ic+1
  endif
  if(jd.lt.0.0d0) then
     jc = jc+1
  endif
  if(kd.lt.0.0d0) then
     kc = kc+1
  endif
  ! 1.3 interpolate                                                                                                                                   
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+id*(i-1))  ! x coordinate                                                                                                      
           fi(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sp(1),fj(j),foo)
        hj(j)=yc(jc+jd*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sp(2),fk(k),foo)
     hk(k)=zc(kc+kd*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sp(3),P_Sp,foo)
  ! 2.1 compute S-1
  Sm(1) = C(1)+ds*nm(1)
  Sm(2) = C(2)+ds*nm(2)
  Sm(3) = C(3)+ds*nm(3)
  ! 2.2 indices
  ic = int((Sm(1)-xc(0))*dx1)
  jc = int((Sm(2)-yc(0))*dy1)
  kc = int((Sm(3)-zc(0))*dz1)
  if(idm.lt.0.0d0) then
     ic = ic+1
  endif
  if(jdm.lt.0.0d0) then
     jc = jc+1
  endif
  if(kdm.lt.0.0d0) then
     kc = kc+1
  endif
  ! 2.3 interpolate                                                                                                                                   
  do k=1,nijk
     do j=1,nijk
        do i=1,nijk
           hi(i)=xc(ic+idm*(i-1))  ! x coordinate                                                                               
           fi(i)=p(ic+idm*(i-1),jc+jdm*(j-1),kc+kdm*(k-1))
        enddo
        call interpolate(hi,fi,nijk,Sm(1),fj(j),foo)
        hj(j)=yc(jc+jdm*(j-1))
     enddo
     call interpolate(hj,fj,nijk,Sm(2),fk(k),foo)
     hk(k)=zc(kc+kdm*(k-1))
  enddo
  call interpolate(hk,fk,nijk,Sm(3),P_Sm,foo)
  !---------------------------------------------------------                                                                                          
  Fx_C = -0.50d0*(P_Sp+P_Sm+pJ_C-dpdn_C*ds)*n(1)+Re1*dudn_C(1)
  Fy_C = -0.50d0*(P_Sp+P_Sm+pJ_C-dpdn_C*ds)*n(2)+Re1*dudn_C(2)
  Fz_C = -0.50d0*(P_Sp+P_Sm+pJ_C-dpdn_C*ds)*n(3)+Re1*dudn_C(3)  
  !-----------------------------------------------------------
  T_fx = (/Fx_A,Fx_B,Fx_C/)
  T_fy = (/Fy_A,Fy_B,Fy_C/)
  T_fz = (/Fz_A,Fz_B,Fz_C/)

  return
  
end subroutine Force4Panel
