subroutine JC4uvwp
  !-----------------------------------------------------------------------
  ! Description:
  ! This subroutine computes the principal and Cartesian jump conditions at
  ! the vertices of the each triangle panel and then interpolate the jump
  ! condition at the intersection points between the grid-lines and the panel.
  !----------------------------------------------------------------------
  !use
  use para
  use field
  use Lagrange
  !-----------------------------------------------------------------------
  !Declaration:
  implicit none
  !---------------------------------------------------------------
  integer :: iobj,ip,nv,np,nA,nB,nC
  double precision, dimension(1:3) :: Xcc,Ucc,dotUcc,Omegac,dotOmegac
  double precision, dimension(1:3) :: A,B,C,M1,M2,M3,M,n
  double precision, dimension(1:3) :: uPJC_A,uPJC_B,uPJC_C
  double precision, dimension(1:3) :: u_A,u_B,u_C,u_M1,u_M2,u_M3
  double precision, dimension(1:3) :: dudnn,dudnp,dudnp_A,dudnp_B,dudnp_C
  double precision, dimension(1:3) :: dudnp_M1,dudnp_M2,dudnp_M3
  double precision, dimension(1:3) :: dudnPJC_A,dudnPJC_B,dudnPJC_C
  double precision, dimension(1:3) :: dudnPJC_M1,dudnPJC_M2,dudnPJC_M3
  double precision, dimension(1:3) :: d2udn2p,dudnnPJC_A,dudnnPJC_B,dudnnPJC_C
  double precision, dimension(1:3) :: dudnnPJC_M1,dudnnPJC_M2,dudnnPJC_M3
  double precision :: ds,KK
  double precision :: l1,l2,den,numx,numy,numz
  double precision, dimension(1:3) :: tau,beta
  double precision, dimension(1:3) :: dudxJC_A,dudyJC_A,dudzJC_A
  double precision, dimension(1:3) :: dudxJC_B,dudyJC_B,dudzJC_B
  double precision, dimension(1:3) :: dudxJC_C,dudyJC_C,dudzJC_C

  double precision, dimension(1:7) :: rhs
  double precision, dimension(1:6) :: jumpd2u,jumpd2p
  
  double precision, dimension(1:3) :: dudxxJC_A,dudyyJC_A,dudzzJC_A     ![d2u/dxidxj] @ A,B,&C 
  double precision, dimension(1:3) :: dudxyJC_A,dudxzJC_A,dudyzJC_A
  double precision, dimension(1:3) :: dudxxJC_B,dudyyJC_B,dudzzJC_B
  double precision, dimension(1:3) :: dudxyJC_B,dudxzJC_B,dudyzJC_B
  double precision, dimension(1:3) :: dudxxJC_C,dudyyJC_C,dudzzJC_C
  double precision, dimension(1:3) :: dudxyJC_C,dudxzJC_C,dudyzJC_C
  
  double precision :: dpdxJC_A,dpdyJC_A,dpdzJC_A    ![dp/dxj] @ A,B,&C 
  double precision :: dpdxJC_B,dpdyJC_B,dpdzJC_B
  double precision :: dpdxJC_C,dpdyJC_C,dpdzJC_C
  double precision :: dpdxJC_M1,dpdyJC_M1,dpdzJC_M1
  double precision :: dpdxJC_M2,dpdyJC_M2,dpdzJC_M2
  double precision :: dpdxJC_M3,dpdyJC_M3,dpdzJC_M3

  double precision :: dpdnPJC_A,dpdnPJC_B,dpdnPJC_C

  double precision :: dpdxxJC_A,dpdxyJC_A,dpdxzJC_A   ![d2p/dxidxj] @ A,B,&C!  
  double precision :: dpdyyJC_A,dpdyzJC_A,dpdzzJC_A
  double precision :: dpdxxJC_B,dpdxyJC_B,dpdxzJC_B
  double precision :: dpdyyJC_B,dpdyzJC_B,dpdzzJC_B
  double precision :: dpdxxJC_C,dpdxyJC_C,dpdxzJC_C
  double precision :: dpdyyJC_C,dpdyzJC_C,dpdzzJC_C
  double precision :: dpdnnJC_A,dpdnnJC_B,dpdnnJC_C    ![d2p/dn2] @ A,B,&C       
  double precision :: dudxp,dudyp,dudzp      !dui/dxj|+                                   
  double precision :: dvdxp,dvdyp,dvdzp
  double precision :: dwdxp,dwdyp,dwdzp
  double precision :: dudxn,dudyn,dudzn      !dui/dxj|-              
  double precision :: dvdxn,dvdyn,dvdzn
  double precision :: dwdxn,dwdyn,dwdzn
  double precision :: Sp  
  double precision, dimension(1:3) :: R1,R2,R3,numx2,numy2,numz2
  double precision :: dpdt_A,dpdt_B,dpdt_C          ! dp/dtau              
  double precision :: dpdt_M1,dpdt_M2,dpdt_M3       
  double precision :: dpdb_A,dpdb_B,dpdb_C          ! dp/dbeta             
  double precision :: dpdb_M1,dpdb_M2,dpdb_M3    
  double precision :: rhs1,rhs2
  double precision, dimension(nvertices) :: rhs4pJC,pJC
  
  double precision,dimension(36) :: A_JC,B_JC,C_JC
  double precision :: ave4rhs,tol
  integer ::cnt
  double precision:: rhs_A,rhs_B,rhs_C
  !-----------------------------------------------------------------------
  double precision, dimension(:), allocatable :: xx,mrhs,r,Ap,pp
  double precision :: ceps,err,dotrr,alpha,betaa
  integer:: iter,nv1,nv2,np1,np2,iv
  !-----------------------------------------------------------------------
  double precision, dimension(9) :: Adum,Adup,Bdum,Bdup,Cdum,Cdup
  !-----------------------------------------------------------------------
  interface
     function crossproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision, dimension(3) :: crossproduct
     end function crossproduct
     !----------------------------------------------------
     function dotproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision :: dotproduct
     end function dotproduct
     !---------------------------------------------------
     subroutine cjc4d2u(tau,beta,rhs,jumpd2u)
       double precision, dimension(7), intent(in) :: rhs
       double precision, dimension(3), intent(in) :: tau,beta
       double precision, dimension(6), intent(out) :: jumpd2u
       double precision, dimension(7) :: dd
       integer :: iv
       integer, dimension(7) :: id
       integer, dimension(6) :: ix
       double precision, dimension(6):: xx
       double precision :: fo1,fo2,fo3,bo1,bo2,bo3
     end subroutine cjc4d2u
     !-----------------------------------------------------------------
     subroutine dudnPJC(s0,u_s0,n,ds,dudn)
       double precision, dimension(3), intent(in) :: s0,u_s0,n
       double precision, intent(in) :: ds
       double precision, dimension(3), intent(out) :: dudn
       integer, parameter:: nn=2,nijk=3
       double precision :: uu(nn),vv(nn),ww(nn),hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
       integer :: id,jd,kd,inn,ic,ie,jc,je,kc,ke,iu,ju,ku,iv,jv,kv,iw,jw,kw
       double precision, dimension(3) :: s,foo
     end subroutine dudnPJC
     !------------------------------------------------------------------
     subroutine dudnnPJC(s0,u_s0,n,ds,dudnn)
       double precision, dimension(3), intent(in) :: s0,u_s0,n
       double precision, intent(in) :: ds
       double precision, dimension(3), intent(out) :: dudnn
       
       integer, parameter:: nn=3,nijk=4
       double precision :: uu(nn),vv(nn),ww(nn),hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
       integer :: id,jd,kd,inn,ic,ie,jc,je,kc,ke,iu,ju,ku,iv,jv,kv,iw,jw,kw
       double precision, dimension(3) :: s,foo
     end subroutine dudnnPJC
     !--------------------------------------------------------------------------------
     subroutine IrrP_JC2(A,B,C,tau,beta,n,A_JC,B_JC,C_JC,ip,iobj,Xcc,Ucc,Omegac,&
          Adup,Adum,Bdup,Bdum,Cdup,Cdum)
       integer :: s
       integer :: ic0,ic1,ie0,ie1,jc0,jc1,je0,je1,kc0,kc1,ke0,ke1
       double precision, dimension(1:3) :: A2,B2,C2,P2,Q
       double precision :: coord3
       integer :: kc,ke,je,jc,ie,ic  
       double precision, dimension(3), intent(in) :: A,B,C,tau,beta,n,Xcc,Ucc,Omegac
       double precision,dimension(36), intent(in) ::A_JC,B_JC,C_JC
       double precision,dimension(9), intent(in) :: Adup,Adum,Bdup,Bdum,Cdup,Cdum
       double precision, dimension(36) :: Q_JC
       integer, intent(in) :: ip,iobj
     end subroutine IrrP_JC2
     !---------------------------------------------------------------------
     subroutine IrrP_pressureJC2(nrow,pJC)
       integer, intent(in) :: nrow
       double precision,dimension(nrow),intent(in):: pJC
       integer :: in,ip,iobj,nA,nB,nC
       integer :: ic,iff,jc,jf,kc,kf
       double precision, dimension(1:3) :: Q,A,B,C
       double precision :: fA,fB,fC,fQ
     end subroutine IrrP_pressureJC2
     subroutine CartJump1st4u(A,B,C,n,dudnPJC,dudxJC,dudyJC,dudzJC)
       double precision,dimension(3), intent(in) :: A,B,C,n,dudnPJC
       double precision,dimension(3), intent(out) :: dudxJC,dudyJC,dudzJC
       double precision,dimension(3):: tau,beta
       double precision :: l1,l2,den,numx,numy,numz,d1,d2,d3
     end subroutine CartJump1st4u
     subroutine CartJump2nd4u(A,B,C,n,dudnnPJC_A,dudxJC_A,dudyJC_A,dudzJC_A,&
          dudxJC_B,dudyJC_B,dudzJC_B,dudxJC_C,dudyJC_C,dudzJC_C,&
          dudxxJC_A,dudxyJC_A,dudxzJC_A,dudyyJC_A,dudyzJC_A,dudzzJC_A)


       double precision,dimension(3), intent(in) :: A,B,C,n,dudnnPJC_A
       double precision,dimension(3), intent(in) :: dudxJC_A,dudyJC_A,dudzJC_A
       double precision,dimension(3), intent(in) :: dudxJC_B,dudyJC_B,dudzJC_B
       double precision,dimension(3), intent(in) :: dudxJC_C,dudyJC_C,dudzJC_C
       
       double precision,dimension(3), intent(out) ::dudxxJC_A,dudxyJC_A,dudxzJC_A
       double precision,dimension(3), intent(out) ::dudyyJC_A,dudyzJC_A,dudzzJC_A
       
       double precision,dimension(3) :: tau,beta
       double precision :: l1,l2
       double precision,dimension(7) :: rhs
       double precision,dimension(6) :: jumpd2u
     end subroutine CartJump2nd4u
     subroutine CartJump2nd4p(A,B,C,n,Omegac,dudnP_A,&
          dpdx_A,dpdy_A,dpdz_A,dpdx_B,dpdy_B,dpdz_B,dpdx_C,dpdy_C,dpdz_C,&
          dpdnnPJC_A,dpdxx_A,dpdxy_A,dpdxz_A,dpdyy_A,dpdyz_A,dpdzz_A,&
          Adum,Adup)

       double precision,dimension(3), intent(in) :: A,B,C,n
       double precision,dimension(3), intent(in) :: Omegac,dudnP_A
       double precision, intent(in) :: dpdx_A,dpdy_A,dpdz_A
       double precision, intent(in) :: dpdx_B,dpdy_B,dpdz_B
       double precision, intent(in) :: dpdx_C,dpdy_C,dpdz_C

       double precision ,intent(out):: dpdnnPJC_A
       double precision ,intent(out):: dpdxx_A,dpdxy_A,dpdxz_A
       double precision ,intent(out):: dpdyy_A,dpdyz_A,dpdzz_A
       double precision, dimension(9), intent(out) :: Adum,Adup
       double precision,dimension(3) :: tau,beta,R1,R2,R3,numx,numy,numz
       double precision :: l1,l2,den,Sp
       double precision :: dudxn,dudyn,dudzn,dudxp,dudyp,dudzp
       double precision :: dvdxn,dvdyn,dvdzn,dvdxp,dvdyp,dvdzp
       double precision :: dwdxn,dwdyn,dwdzn,dwdxp,dwdyp,dwdzp
       double precision,dimension(7) :: rhs
       double precision,dimension(6) :: jumpd2p

       double precision ::R1i,R1j,R1k
       double precision ::R2i,R2j,R2k
       double precision ::R3i,R3j,R3k

       double precision ::unumx,unumy,unumz
       double precision ::vnumx,vnumy,vnumz
       double precision ::wnumx,wnumy,wnumz

     end subroutine CartJump2nd4p
     subroutine rhs4p(A,B,C,M1,M2,M3,&
          dpdxJC_A,dpdyJC_A,dpdzJC_A,dpdxJC_B,dpdyJC_B,dpdzJC_B,&
          dpdxJC_C,dpdyJC_C,dpdzJC_C,dpdxJC_M1,dpdyJC_M1,dpdzJC_M1,&
          dpdxJC_M2,dpdyJC_M2,dpdzJC_M2,dpdxJC_M3,dpdyJC_M3,dpdzJC_M3,&
          rhs_A,rhs_B,rhs_C)
       
       double precision, dimension(3), intent(in) :: A,B,C,M1,M2,M3
       double precision, intent(in) :: dpdxJC_A,dpdyJC_A,dpdzJC_A
       double precision, intent(in) :: dpdxJC_B,dpdyJC_B,dpdzJC_B
       double precision, intent(in) :: dpdxJC_C,dpdyJC_C,dpdzJC_C
       double precision, intent(in) :: dpdxJC_M1,dpdyJC_M1,dpdzJC_M1
       double precision, intent(in) :: dpdxJC_M2,dpdyJC_M2,dpdzJC_M2
       double precision, intent(in) :: dpdxJC_M3,dpdyJC_M3,dpdzJC_M3
       double precision, intent(out) :: rhs_A,rhs_B,rhs_C 
       double precision, dimension(3):: tau,beta
       double precision :: l1,l2,rhs1,rhs2,rhs3
       double precision :: dpdt_A,dpdt_B,dpdt_C,dpdt_M1,dpdt_M2,dpdt_M3
       double precision :: dpdb_A,dpdb_B,dpdb_C,dpdb_M1,dpdb_M2,dpdb_M3
     end subroutine rhs4p
     
     subroutine CG4p(nv1,nv2,np1,np2,rhs,xx)
       integer, intent(in):: nv1,nv2,np1,np2
       double precision, dimension(nv1+1:nv2),intent(in) :: rhs
       double precision, dimension(nv1+1:nv2),intent(out) :: xx
       integer:: i,iv,ip,nv,nA,nB,nC,iter
       double precision, dimension(nv1+1:nv2)::r,pp,ADiag,mrhs,zz,Ap
       double precision :: eps, err,dotrr,dotrz,alpha,beta
     end subroutine CG4p

  end interface
  !------------------------------------------------------------------------
  !internal
  !-------------------------------------------------------------------------
  rhs4pJC = 0.0d0
  !-------------------------------------------------------------------------
  ds = 1.01d0*sqrt(dx*dx+dy*dy+dz*dz) 
  !-----------------------------------------------------------------------
  do iobj =1,nobj
     nv1 = nv4obj(iobj-1)
     nv2 = nv4obj(iobj)
     np1 = np4obj(iobj-1)
     np2 = np4obj(iobj)
     
     nv = nv2-nv1
     np = np2-np1
     !xsc and xsc0
     Xcc = (/xsc(iobj),ysc(iobj),zsc(iobj)/)
     Ucc = (/xsct(iobj),ysct(iobj),zsct(iobj)/)
     dotUcc = (/xsctt(iobj),ysctt(iobj),zsctt(iobj)/)
     Omegac = (/omegax(iobj),omegay(iobj),omegaz(iobj)/)
     dotOmegac = (/omegaxt(iobj),omegayt(iobj),omegazt(iobj)/)

     do ip = np1+1,np2
        !-------------------------------------------------------------------
        ! index of the panel's vertices
        nA = panel(1,ip)+nv1
        nB = panel(2,ip)+nv1
        nC = panel(3,ip)+nv1
        ! vertices of the panel       
        A = vertex(:,nA)
        B = vertex(:,nB)
        C = vertex(:,nC)
        !-------------------------------------------------------------------
        M1 = 0.5*(A+B)  !middle point of AB   
        M2 = 0.5*(B+C)  !middle point of BC                
        M3 = 0.5*(A+C)  !middle point of AC      
        M = 0.5*(M1+M3)  !middle point of DF and AE
        !------------------------------------------------------------------
        ! normal vector to the panel 
        n = crossproduct(B-A,C-A)
        n = n/sqrt(dotproduct(n,n))
        !------------------------------------------------------------------
        !==================================================================
        !======== Jump Conditions for velocity=============================
        !==================================================================
        !--------PJC for velocity  [u] ------------------------------------
        ! [u] = 0 
        uPJC_A = 0.0d0
        uPJC_B = 0.0d0
        uPJC_C = 0.0d0
        ! ------PJC for the 1st order derivative [du/dn]-------------------
        ! 1. velocity at the vertices A,B &C (surface velocity)
        ! 2. du/dn |-, du/dn|+ at A,B,and C
        ! 3. [du/dn] = du/dn|+  -  du/dn|- @ A,B &C
        u_A  = Ucc  + crossproduct(Omegac,(A-Xcc))
        u_B  = Ucc  + crossproduct(Omegac,(B-Xcc))
        u_C  = Ucc  + crossproduct(Omegac,(C-Xcc))
        u_M1 = Ucc  + crossproduct(Omegac,(M1-Xcc))
        u_M2 = Ucc  + crossproduct(Omegac,(M2-Xcc))
        u_M3 = Ucc  + crossproduct(Omegac,(M3-Xcc))
        !------------------------------------------------------------------
        ! PJC for 1st normal derivative of u
        dudnn = crossproduct(Omegac,n)
        call dudnPJC(A,u_A,n,ds,dudnp_A)
        dudnPJC_A = dudnp_A - dudnn
        call dudnPJC(B,u_B,n,ds,dudnp_B)
        dudnPJC_B = dudnp_B - dudnn
        call dudnPJC(C,u_C,n,ds,dudnp_C)
        dudnPJC_C = dudnp_C - dudnn
        call dudnPJC(M1,u_M1,n,ds,dudnp_M1)
        dudnPJC_M1 = dudnp_M1 - dudnn
        call dudnPJC(M2,u_M2,n,ds,dudnp_M2)
        dudnPJC_M2 = dudnp_M2 - dudnn
        call dudnPJC(M3,u_M3,n,ds,dudnp_M3)
        dudnPJC_M3 = dudnp_M3 - dudnn
        !------PJCs for the 2nd normal derivative for u [d2u/dn2]----------------
        call dudnnPJC(A,u_A,n,ds,d2udn2p)
        KK = 4.0d0
        dudnnPJC_A = d2udn2p + KK*dudnPJC_A 
        call dudnnPJC(B,u_B,n,ds,d2udn2p)
        dudnnPJC_B = d2udn2p + KK*dudnPJC_B
        call dudnnPJC(C,u_C,n,ds,d2udn2p)
        dudnnPJC_C = d2udn2p + KK*dudnPJC_C
        call dudnnPJC(M1,u_M1,n,ds,d2udn2p)
        dudnnPJC_M1 = d2udn2p + KK*dudnPJC_M1
        call dudnnPJC(M2,u_M2,n,ds,d2udn2p)
        dudnnPJC_M2 = d2udn2p + KK*dudnPJC_M2
        call dudnnPJC(M3,u_M3,n,ds,d2udn2p)
        dudnnPJC_M3 = d2udn2p + KK*dudnPJC_M3
        !-------------------------------------------------------------------
        !------------------------------------------------------------------
        ! CJCs for 1st derivative of velocity [du/dxi]
        ! A:
        call CartJump1st4u(A,B,C,n,dudnPJC_A,dudxJC_A,dudyJC_A,dudzJC_A)
        ! B:
        call CartJump1st4u(B,C,A,n,dudnPJC_B,dudxJC_B,dudyJC_B,dudzJC_B)
        ! C:
        call CartJump1st4u(C,A,B,n,dudnPJC_C,dudxJC_C,dudyJC_C,dudzJC_C)
        ! CJCs for the 2nd derivatives for velocity
        !A:
        call CartJump2nd4u(A,B,C,n,dudnnPJC_A,dudxJC_A,dudyJC_A,dudzJC_A,&
             dudxJC_B,dudyJC_B,dudzJC_B,dudxJC_C,dudyJC_C,dudzJC_C,&
             dudxxJC_A,dudxyJC_A,dudxzJC_A,dudyyJC_A,dudyzJC_A,dudzzJC_A)
        !B:
        call CartJump2nd4u(B,C,A,n,dudnnPJC_B,dudxJC_B,dudyJC_B,dudzJC_B,&
             dudxJC_C,dudyJC_C,dudzJC_C,dudxJC_A,dudyJC_A,dudzJC_A,&
             dudxxJC_B,dudxyJC_B,dudxzJC_B,dudyyJC_B,dudyzJC_B,dudzzJC_B)
        !C:
        call CartJump2nd4u(C,A,B,n,dudnnPJC_C,dudxJC_C,dudyJC_C,dudzJC_C,&
             dudxJC_A,dudyJC_A,dudzJC_A,dudxJC_B,dudyJC_B,dudzJC_B,&
             dudxxJC_C,dudxyJC_C,dudxzJC_C,dudyyJC_C,dudyzJC_C,dudzzJC_C)
        !=================================================================
        !====================================================================
        !---------- Jump conditions for pressure----------------------------
        !===================================================================
        ! CJCs for the 1st derivatives of pressure [dp/dxi]
        !-----@ A----------------------------------------------------------
        dpdxJC_A = Re1*dudnnPJC_A(1) - dotOmegac(2)*(A(3)-Xcc(3))+dotOmegac(3)*(A(2)-Xcc(2))
        dpdyJC_A = Re1*dudnnPJC_A(2) - dotOmegac(3)*(A(1)-Xcc(1))+dotOmegac(1)*(A(3)-Xcc(3))
        dpdzJC_A = Re1*dudnnPJC_A(3) - dotOmegac(1)*(A(2)-Xcc(2))+dotOmegac(2)*(A(1)-Xcc(1))
        !----@ B------------------------------------------------------------
        dpdxJC_B = Re1*dudnnPJC_B(1) - dotOmegac(2)*(B(3)-Xcc(3))+dotOmegac(3)*(B(2)-Xcc(2))
        dpdyJC_B = Re1*dudnnPJC_B(2) - dotOmegac(3)*(B(1)-Xcc(1))+dotOmegac(1)*(B(3)-Xcc(3))
        dpdzJC_B = Re1*dudnnPJC_B(3) - dotOmegac(1)*(B(2)-Xcc(2))+dotOmegac(2)*(B(1)-Xcc(1))
        !---@ C----------------------------------------------------------
        dpdxJC_C = Re1*dudnnPJC_C(1) - dotOmegac(2)*(C(3)-Xcc(3))+dotOmegac(3)*(C(2)-Xcc(2))
        dpdyJC_C = Re1*dudnnPJC_C(2) - dotOmegac(3)*(C(1)-Xcc(1))+dotOmegac(1)*(C(3)-Xcc(3))
        dpdzJC_C = Re1*dudnnPJC_C(3) - dotOmegac(1)*(C(2)-Xcc(2))+dotOmegac(2)*(C(1)-Xcc(1))
        !---@ M1----------------------------------------------------------                            
        dpdxJC_M1 = Re1*dudnnPJC_M1(1) - dotOmegac(2)*(M1(3)-Xcc(3))+dotOmegac(3)*(M1(2)-Xcc(2))
        dpdyJC_M1 = Re1*dudnnPJC_M1(2) - dotOmegac(3)*(M1(1)-Xcc(1))+dotOmegac(1)*(M1(3)-Xcc(3))
        dpdzJC_M1 = Re1*dudnnPJC_M1(3) - dotOmegac(1)*(M1(2)-Xcc(2))+dotOmegac(2)*(M1(1)-Xcc(1))
        !---@ M2----------------------------------------------------------                                 
        dpdxJC_M2 = Re1*dudnnPJC_M2(1) - dotOmegac(2)*(M2(3)-Xcc(3))+dotOmegac(3)*(M2(2)-Xcc(2))
        dpdyJC_M2 = Re1*dudnnPJC_M2(2) - dotOmegac(3)*(M2(1)-Xcc(1))+dotOmegac(1)*(M2(3)-Xcc(3))
        dpdzJC_M2 = Re1*dudnnPJC_M2(3) - dotOmegac(1)*(M2(2)-Xcc(2))+dotOmegac(2)*(M2(1)-Xcc(1))
        !---@ M3----------------------------------------------------------                                
        dpdxJC_M3 = Re1*dudnnPJC_M3(1) - dotOmegac(2)*(M3(3)-Xcc(3))+dotOmegac(3)*(M3(2)-Xcc(2))
        dpdyJC_M3 = Re1*dudnnPJC_M3(2) - dotOmegac(3)*(M3(1)-Xcc(1))+dotOmegac(1)*(M3(3)-Xcc(3))
        dpdzJC_M3 = Re1*dudnnPJC_M3(3) - dotOmegac(1)*(M3(2)-Xcc(2))+dotOmegac(2)*(M3(1)-Xcc(1))
        !====================================================================
        !====================================================================
        ! PJC for 1st derivative of p at the verticies
        dpdnPJC_A = dotproduct( (/dpdxJC_A,dpdyJC_A,dpdzJC_A/) ,n)
        dpdnPJC_B = dotproduct( (/dpdxJC_B,dpdyJC_B,dpdzJC_B/) ,n)
        dpdnPJC_C = dotproduct( (/dpdxJC_C,dpdyJC_C,dpdzJC_C/) ,n)
        !===========================================================================================
        ! Jump condition for the 2nd order derivative of pressure
        !--------------------------------------------------------------------------------------------
        call CartJump2nd4p(A,B,C,n,Omegac,dudnp_A,&
             dpdxJC_A,dpdyJC_A,dpdzJC_A,dpdxJC_B,dpdyJC_B,dpdzJC_B,dpdxJC_C,dpdyJC_C,dpdzJC_C,&
             dpdnnJC_A,dpdxxJC_A,dpdxyJC_A,dpdxzJC_A,dpdyyJC_A,dpdyzJC_A,dpdzzJC_A,&
             Adum,Adup)
        call CartJump2nd4p(B,C,A,n,Omegac,dudnp_B,&
             dpdxJC_B,dpdyJC_B,dpdzJC_B,dpdxJC_C,dpdyJC_C,dpdzJC_C,dpdxJC_A,dpdyJC_A,dpdzJC_A,&
             dpdnnJC_B,dpdxxJC_B,dpdxyJC_B,dpdxzJC_B,dpdyyJC_B,dpdyzJC_B,dpdzzJC_B,&
             Bdum,Bdup)
        call CartJump2nd4p(C,A,B,n,Omegac,dudnp_C,&
             dpdxJC_C,dpdyJC_C,dpdzJC_C,dpdxJC_A,dpdyJC_A,dpdzJC_A,dpdxJC_B,dpdyJC_B,dpdzJC_B,&
             dpdnnJC_C,dpdxxJC_C,dpdxyJC_C,dpdxzJC_C,dpdyyJC_C,dpdyzJC_C,dpdzzJC_C,&
             Cdum,Cdup)
        !============================================================================================
        !---------------------------------------------------------------------------------------------
        ! RHS for the  system for [p]
        call rhs4p(A,B,C,M1,M2,M3,&
             dpdxJC_A,dpdyJC_A,dpdzJC_A,dpdxJC_B,dpdyJC_B,dpdzJC_B,&
             dpdxJC_C,dpdyJC_C,dpdzJC_C,dpdxJC_M1,dpdyJC_M1,dpdzJC_M1,&
             dpdxJC_M2,dpdyJC_M2,dpdzJC_M2,dpdxJC_M3,dpdyJC_M3,dpdzJC_M3,&
             rhs_A,rhs_B,rhs_C)
        rhs4pJC(nA) = rhs4pJC(nA) +rhs_A
        rhs4pJC(nB) = rhs4pJC(nB) +rhs_B
        rhs4pJC(nC) = rhs4pJC(nC) +rhs_C
        !============================================================================================
        !============================================================================================
        !============================================================================================
        ! store JC
        ! interpolate
        A_JC(1:3)   = (/ dudxJC_A(1),dudyJC_A(1),dudzJC_A(1)  /)
        A_JC(4:6)   = (/ dudxJC_A(2),dudyJC_A(2),dudzJC_A(2)  /)
        A_JC(7:9)   = (/ dudxJC_A(3),dudyJC_A(3),dudzJC_A(3)  /)
        A_JC(10:12) = (/ dpdxJC_A   ,dpdyJC_A   ,dpdzJC_A     /)
        A_JC(13:18) = (/ dudxxJC_A(1),dudxyJC_A(1),dudxzJC_A(1),  &
                         dudyyJC_A(1),dudyzJC_A(1),dudzzJC_A(1) /)
        A_JC(19:24) = (/ dudxxJC_A(2),dudxyJC_A(2),dudxzJC_A(2),  &
                         dudyyJC_A(2),dudyzJC_A(2),dudzzJC_A(2) /)
        A_JC(25:30) = (/ dudxxJC_A(3),dudxyJC_A(3),dudxzJC_A(3),  &
                         dudyyJC_A(3),dudyzJC_A(3),dudzzJC_A(3) /)
        A_JC(31:36) = (/ dpdxxJC_A,dpdxyJC_A,dpdxzJC_A,  &
                         dpdyyJC_A,dpdyzJC_A,dpdzzJC_A /)

        B_JC(1:3)   = (/ dudxJC_B(1),dudyJC_B(1),dudzJC_B(1)  /)
        B_JC(4:6)   = (/ dudxJC_B(2),dudyJC_B(2),dudzJC_B(2)  /)
        B_JC(7:9)   = (/ dudxJC_B(3),dudyJC_B(3),dudzJC_B(3)  /)
        B_JC(10:12) = (/ dpdxJC_B   ,dpdyJC_B   ,dpdzJC_B     /)
        B_JC(13:18) = (/ dudxxJC_B(1),dudxyJC_B(1),dudxzJC_B(1),  &
                         dudyyJC_B(1),dudyzJC_B(1),dudzzJC_B(1) /)
        B_JC(19:24) = (/ dudxxJC_B(2),dudxyJC_B(2),dudxzJC_B(2),  &
                         dudyyJC_B(2),dudyzJC_B(2),dudzzJC_B(2) /)
        B_JC(25:30) = (/ dudxxJC_B(3),dudxyJC_B(3),dudxzJC_B(3),  &
                         dudyyJC_B(3),dudyzJC_B(3),dudzzJC_B(3) /)
        B_JC(31:36) = (/ dpdxxJC_B,dpdxyJC_B,dpdxzJC_B,  &
                         dpdyyJC_B,dpdyzJC_B,dpdzzJC_B /)

        C_JC(1:3)   = (/ dudxJC_C(1),dudyJC_C(1),dudzJC_C(1)  /)
        C_JC(4:6)   = (/ dudxJC_C(2),dudyJC_C(2),dudzJC_C(2)  /)
        C_JC(7:9)   = (/ dudxJC_C(3),dudyJC_C(3),dudzJC_C(3)  /)
        C_JC(10:12) = (/ dpdxJC_C  ,dpdyJC_C   ,dpdzJC_C      /)
        C_JC(13:18) = (/ dudxxJC_C(1),dudxyJC_C(1),dudxzJC_C(1),  &
                         dudyyJC_C(1),dudyzJC_C(1),dudzzJC_C(1) /)
        C_JC(19:24) = (/ dudxxJC_C(2),dudxyJC_C(2),dudxzJC_C(2),  &
                         dudyyJC_C(2),dudyzJC_C(2),dudzzJC_C(2) /)
        C_JC(25:30) = (/ dudxxJC_C(3),dudxyJC_C(3),dudxzJC_C(3),  &
                         dudyyJC_C(3),dudyzJC_C(3),dudzzJC_C(3) /)
        C_JC(31:36) = (/ dpdxxJC_C,dpdxyJC_C,dpdxzJC_C,  &
                         dpdyyJC_C,dpdyzJC_C,dpdzzJC_C /)       
        !================================================================
        tau = B-A
        l1 = sqrt(dotproduct(tau,tau))
        tau = tau/l1
        
        beta = C-A
        l2 = sqrt(dotproduct(beta,beta))
        beta = beta/l2

        call IrrP_JC2(A,B,C,tau,beta,n,A_JC,B_JC,C_JC,ip,iobj,Xcc,Ucc,Omegac,&
             Adup,Adum,Bdup,Bdum,Cdup,Cdum)

        !---------------------------------
        dudnp_T(:,ip)   = (/dudxJC_A(1)*n(1)+dudyJC_A(1)*n(2)+dudzJC_A(1)*n(3)&
             ,dudxJC_B(1)*n(1)+dudyJC_B(1)*n(2)+dudzJC_B(1)*n(3)&
             ,dudxJC_C(1)*n(1)+dudyJC_C(1)*n(2)+dudzJC_C(1)*n(3)/)
        dvdnp_T(:,ip)   = (/dudxJC_A(2)*n(1)+dudyJC_A(2)*n(2)+dudzJC_A(2)*n(3)&
             ,dudxJC_B(2)*n(1)+dudyJC_B(2)*n(2)+dudzJC_B(2)*n(3)&
             ,dudxJC_C(2)*n(1)+dudyJC_C(2)*n(2)+dudzJC_C(2)*n(3)/)

        dwdnp_T(:,ip)   = (/dudxJC_A(3)*n(1)+dudyJC_A(3)*n(2)+dudzJC_A(3)*n(3)&
             ,dudxJC_B(3)*n(1)+dudyJC_B(3)*n(2)+dudzJC_B(3)*n(3)&
             ,dudxJC_C(3)*n(1)+dudyJC_C(3)*n(2)+dudzJC_C(3)*n(3)/)
        dpdnPJC_T(:,ip) = (/dpdnPJC_A,dpdnPJC_B,dpdnPJC_C/)
     enddo
  enddo
  !==============================
  do iobj =1,nobj
     nv1 = nv4obj(iobj-1)
     nv2 = nv4obj(iobj)
     np1 = np4obj(iobj-1)
     np2 = np4obj(iobj)

     ave4rhs = sum(rhs4pJC(nv1+1:nv2))/(nv2-nv1)
     rhs4pJC(nv1+1:nv2) = rhs4pJC(nv1+1:nv2)-ave4rhs
     ave4rhs = sum(rhs4pJC(nv1+1:nv2))/(nv2-nv1)
     call CG4p(nv1,nv2,np1,np2,rhs4pJC(nv1+1:nv2),pJC(nv1+1:nv2))
     pJC0(nv1+1:nv2)= pJC(nv1+1:nv2)
  enddo
  !===============================================================
  !===============================================================
  !==============================================
  ! interpolate [p] @ intersection points
  call IrrP_pressureJC2(nvertices,pJC)
  
  return
     
end subroutine JC4uvwp
