!-----------------------------------------------------------------------------------
! This subroutine is to compute the fluid force and torque acting on the boundary 
!----------------------------------------------------------------------------------
subroutine time_output
  use para
  use field
  use Lagrange

  !------------------------------------------------------------
  implicit none
  !----------------------------------------------------------
  integer, parameter:: ms=1
  integer ::m,ip,np,np1,np2,nv,nv1,nv2,nA,nB,nC
  double precision,dimension(ms) :: S_Fx,S_Fy,S_Fz,TQx,TQy,TQz
  double precision :: M_Fx,M_Fy,M_Fz,M_TQx,M_TQy,M_TQz
  double precision :: panel_area,ds,pJ_A,pJ_B,pJ_C 
  double precision, dimension(1:3) :: A,B,C,n,ABxAC  
  double precision, dimension(1:3) :: V_fx,V_fy,V_fz
  double precision, dimension(1:3) :: d_A,d_B,d_C,F_A,F_B,F_C,G_A,G_B,G_C
  
  interface
     function crossproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision, dimension(3) :: crossproduct
     end function crossproduct
     function dotproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision :: dotproduct
     end function dotproduct
     subroutine Force4Panel(A,B,C,n,dudn_A,dudn_B,dudn_C,dpdn_A,dpdn_B,dpdn_C,pJ_A,pJ_B,pJ_C,ds,T_fx,T_fy,T_fz)
       double precision,dimension(3), intent(in):: A,B,C,n,dudn_A,dudn_B,dudn_C
       double precision, intent(in):: dpdn_A,dpdn_B,dpdn_C,ds,pJ_A,pJ_B,pJ_C
       double precision, dimension(3), intent(out)::T_fx,T_fy,T_fz
       
       double precision::Fx_A,Fy_A,Fz_A
       double precision::Fx_B,Fy_B,Fz_B
       double precision::Fx_C,Fy_C,Fz_C  
       double precision, dimension(3):: Sp,Sm,foo ! S @ surface point + ds, and S @ surface point minus ds   
       integer :: id,jd,kd,ic,jc,kc
       double precision :: P_Sp,P_Sm
       integer, parameter:: nijk=3
       double precision :: hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
     end subroutine Force4Panel

  end interface
  !-------------------------------------------------------------------------------------
  ds = 1.01d0*sqrt(dx*dx+dy*dy+dz*dz)
  S_Fx = 0.0d0
  S_Fy = 0.0d0
  S_Fz = 0.0d0
  TQx = 0.0d0
  TQy = 0.0d0
  TQz = 0.0d0
  !--------------------------------------------------------------------------------------
  DO m=1,ms
     nv1 = nv4obj(m-1)
     nv2 = nv4obj(m)
     nv =nv2-nv1
     np1 = np4obj(m-1)
     np2 = np4obj(m)
     np =np2-np1
     do ip = np1+1,np2
        nA = panel(1,ip)+nv1
        nB = panel(2,ip)+nv1
        nC = panel(3,ip)+nv1
        A = vertex(:,nA)
        B = vertex(:,nB)
        C = vertex(:,nC)
        ABxAC = crossproduct(B-A,C-A)
        !-----------------------------------------------
        ! normal vector tot eh panel
        n = ABxAC
        n = n/sqrt(dotproduct(n,n))
        !------------------------------------------------
        ! area of the pane;
        panel_area = 0.50d0*dsqrt(dotproduct(ABxAC,ABxAC))
        !!------------------------------------------------------
        ! JCs [p] at the vertices
        pJ_A=pJC0(nA)
        pJ_B=pJC0(nB)
        pJ_C=pJC0(nC)
        !-----------------------------------------------------------------------
        !-----Force------------------------------------------------------------
        ! 1. compute Fx,Fy,&Fz @ the vertices A,B,&C
        ! need:
        !      - Cooedinates A,B, and C and the nurmal vector:
        !
        call Force4Panel(A,B,C,n,(/dudnp_T(1,ip),dvdnp_T(1,ip),dwdnp_T(1,ip)/),&
             (/dudnp_T(2,ip),dvdnp_T(2,ip),dwdnp_T(2,ip)/),&
             (/dudnp_T(3,ip),dvdnp_T(3,ip),dwdnp_T(3,ip)/),&
             dpdnPJC_T(1,ip),dpdnPJC_T(2,ip),dpdnPJC_T(3,ip),&
             pJ_A,pJ_B,pJ_C,ds,V_fx,V_fy,V_fz)
        ! 2. interpolate force from the vertices to the centroid ( geometric center) 
        M_Fx = (1.0d0/3.0d0)*(V_fx(1)+V_fx(2)+V_fx(3))
        M_Fy = (1.0d0/3.0d0)*(V_fy(1)+V_fy(2)+V_fy(3))
        M_Fz = (1.0d0/3.0d0)*(V_fz(1)+V_fz(2)+V_fz(3))
        ! 3. multiply the force @ the centroid by the panel area, and update the result
        S_Fx(m)  = S_Fx(m)  + M_Fx*panel_area
        S_Fy(m)  = S_Fy(m)  + M_Fy*panel_area
        S_Fz(m)  = S_Fz(m)  + M_Fz*panel_area
        !---------------------------------------------------------
        !---------------------------------------------------------
        !---------- Torque---------------------------------------
        ! A: (A-xc) X (Fx_A,Fy_A,Fz_A)
        d_A = A-(/xsc(m),ysc(m),zsc(m)/) ! distance from the center of the object
        F_A = (/V_fx(1),V_fy(1),V_fz(1)/)! the force at the vertex
        G_A=crossproduct( d_A,F_A ) ! torque @ A
        ! B: (B-xc) X (Fx_B,Fy_B,Fz_B)
	d_B = B-(/xsc(m),ysc(m),zsc(m)/)
        F_B = (/V_fx(2),V_fy(2),V_fz(2)/)
        G_B=crossproduct( d_B,F_B )
        ! C: (C-xc) X (Fx_C,Fy_C,Fz_C)
        d_C = C-(/xsc(m),ysc(m),zsc(m)/)
        F_C = (/V_fx(3),V_fy(3),V_fz(3)/)
        G_C=crossproduct( d_C,F_C )
        ! interpolate @ the geometric center of the panel
        M_TQx =  (1.0d0/3.0d0)*(G_A(1)+G_B(1)+G_C(1)) 
        M_TQy =  (1.0d0/3.0d0)*(G_A(2)+G_B(2)+G_C(2))
        M_TQz =  (1.0d0/3.0d0)*(G_A(3)+G_B(3)+G_C(3))
        ! multiply by the panel area and update the result
        TQx(m) = TQx(m) + M_TQx*panel_area
        TQy(m) = TQy(m) + M_TQy*panel_area
        TQz(m) = TQz(m) + M_TQz*panel_area
     enddo
  ENDDO
  
  write(61,100) t,(S_Fx(m),S_Fy(m),S_Fz(m),TQx(m),TQy(m),TQz(m),m=1,ms)
  !print*,t,(S_Fx(m),TQx(m),m=1,ms)

  
100 format(1x,200e16.6e4)

  
  return
end subroutine time_output
  

!-----------------------------------------------------------------------
