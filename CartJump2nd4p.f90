subroutine CartJump2nd4p(A,B,C,n,Omegac,dudnP_A,&
     dpdx_A,dpdy_A,dpdz_A,dpdx_B,dpdy_B,dpdz_B,dpdx_C,dpdy_C,dpdz_C,&
     dpdnnPJC_A,dpdxx_A,dpdxy_A,dpdxz_A,dpdyy_A,dpdyz_A,dpdzz_A,&
     Adum,Adup)


  implicit none

  double precision,dimension(3), intent(in) :: A,B,C,n
  double precision,dimension(3), intent(in) :: Omegac,dudnP_A
  double precision, intent(in) :: dpdx_A,dpdy_A,dpdz_A
  double precision, intent(in) :: dpdx_B,dpdy_B,dpdz_B
  double precision, intent(in) :: dpdx_C,dpdy_C,dpdz_C

  double precision, intent(out):: dpdnnPJC_A
  double precision, intent(out):: dpdxx_A,dpdxy_A,dpdxz_A
  double precision, intent(out):: dpdyy_A,dpdyz_A,dpdzz_A
  double precision, dimension(9), intent(out) :: Adum,Adup

  double precision:: signx,signy,signz
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

  !-----------------
  interface
     function crossproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision, dimension(3) :: crossproduct
     end function crossproduct
     !------------------------                                          
     function dotproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision :: dotproduct
     end function dotproduct
     !------
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
  end interface
  !============
  !--------------------------------------------
  tau = B-A
  l1 = sqrt(dotproduct(tau,tau))
  tau = tau/l1
  !--------------------------------------------
  beta = C-A
  l2 = sqrt(dotproduct(beta,beta))
  beta = beta/l2
  !---------------------------------------------
  !Print*,Omegac
  ! Step.1 :  compute the principal jump condition [d2p/dnn] = [Sp]
  !   we need to compute:
  !                   dui/dxj|-
  !                   dui/dxj|+
  ! dui/dxj|- :
  
  dudxn =  0.0d0
  dudyn = -Omegac(3)
  dudzn =  Omegac(2)
  
  dvdxn =  Omegac(3)
  dvdyn =  0.0d0
  dvdzn = -Omegac(1)

  dwdxn = -Omegac(2)
  dwdyn =  Omegac(1)
  dwdzn =  0.0d0      

  Adum = (/dudxn,dudyn,dudzn,dvdxn,dvdyn,dvdzn,dwdxn,dwdyn,dwdzn/)
  
  ! dui/dxj|+   solve C1 * dui/dxj|+ = RHS
  !         |taux    tauy    tayz |            |(Omega x tau)|i   |
  !   C1 =  |betax   betay   betaz|      RHS = |(Omega x beta)|i  |            
  !         |nx      ny      nz   |            | d2u/dnn_A|i      |
  !dpdnnPJC_A = Sp
  R1 = crossproduct(Omegac,tau)
  R2 = crossproduct(Omegac,beta)
  R3 = dudnP_A
  
  R1i = R1(1)
  R1j = R1(2)
  R1k = R1(3)

  R2i = R2(1)
  R2j = R2(2)
  R2k = R2(3)

  R3i = R3(1)
  R3j = R3(2)
  R3k = R3(3)

  unumx = R1i*( beta(2)*n(3)    - beta(3)*n(2))&
        - R2i*( tau(2) *n(3)    - tau(3) *n(2))&
        + R3i*( tau(2) *beta(3) - tau(3) *beta(2))
  unumy =  R2i*( tau(1) *n(3)    - tau(3) *n(1)    )&
         - R1i*( beta(1)*n(3)    - beta(3)*n(1)    )&
         - R3i*( tau(1) *beta(3) - tau(3) *beta(1) )
  unumz =  R1i*( beta(1)*n(2)    - beta(2)*n(1)    )&
         - R2i*( tau(1) *n(2)    - tau(2) *n(1)    )&
         + R3i*( tau(1) *beta(2) - tau(2) *beta(1) )

  vnumx = R1j*( beta(2)*n(3)    - beta(3)*n(2))&
        - R2j*( tau(2) *n(3)    - tau(3) *n(2))&
        + R3j*( tau(2) *beta(3) - tau(3) *beta(2))
  vnumy =  R2j*( tau(1) *n(3)    - tau(3)*n(1)   ) &
         - R1j*( beta(1)*n(3)    - beta(3)*n(1)    )&
         - R3j*( tau(1) *beta(3) - tau(3) *beta(1) )
  vnumz =  R1j*( beta(1)*n(2)    - beta(2)*n(1)    )&
         - R2j*( tau(1) *n(2)    - tau(2) *n(1)    )&
         + R3j*( tau(1) *beta(2) - tau(2) *beta(1) )

  wnumx = R1k*( beta(2)*n(3)    - beta(3)*n(2))&
        - R2k*( tau(2) *n(3)    - tau(3) *n(2))&
        + R3k*( tau(2) *beta(3) - tau(3) *beta(2))
  wnumy =  R2k*( tau(1) *n(3)    - tau(3) *n(1)    )&
         - R1k*( beta(1)*n(3)    - beta(3)*n(1)    )&
         - R3k*( tau(1) *beta(3) - tau(3) *beta(1) )
  wnumz =  R1k*( beta(1)*n(2)    - beta(2)*n(1)    )&
         - R2k*( tau(1) *n(2)    - tau(2) *n(1)    )&
         + R3k*( tau(1) *beta(2) - tau(2) *beta(1) )!

  den = tau(1)*(beta(2)*n(3)-beta(3)*n(2)) &                         
       -tau(2)*(beta(1)*n(3)-beta(3)*n(1)) &
       +tau(3)*(beta(1)*n(2)-beta(2)*n(1))

  dudxp = unumx/den
  dudyp = unumy/den
  dudzp = unumz/den

  dvdxp = vnumx/den
  dvdyp = vnumx/den
  dvdzp = vnumx/den
  
  dwdxp = wnumx/den
  dwdyp = wnumx/den
  dwdzp = wnumx/den
  Adup = (/dudxp,dudyp,dudzp,dvdxp,dvdyp,dvdzn,dwdxp,dwdyp,dwdzp/)
  
  Sp = 2.0d0*( dudxp*dvdyp - dudyp*dvdxp +&
               dudxp*dwdzp - dudzp*dwdxp +&
               dvdyp*dwdzp - dvdzp*dwdyp - &
               dudxn*dvdyn + dudyn*dvdxn - &
               dudxn*dwdzn +  dudzn*dwdxn -&
               dvdyn*dwdzn + dvdzn*dwdyn)  

  
  dpdnnPJC_A = Sp
  !------------    
  rhs(1) = (dpdx_B-dpdx_A)/l1
  rhs(2) = (dpdy_B-dpdy_A)/l1                
  rhs(3) = (dpdz_B-dpdz_A)/l1                
  rhs(4) = (dpdx_C-dpdx_A)/l2                 
  rhs(5) = (dpdy_C-dpdy_A)/l2               
  rhs(6) = (dpdz_C-dpdz_A)/l2
  rhs(7) = dpdnnPJC_A
  
  call  cjc4d2u(tau,beta,rhs,jumpd2p)

  dpdxx_A=jumpd2p(1)
  dpdxy_A=jumpd2p(2)
  dpdxz_A=jumpd2p(3)
  dpdyy_A=jumpd2p(4)
  dpdyz_A=jumpd2p(5)
  dpdzz_A=jumpd2p(6)

  return
  
end subroutine CartJump2nd4p
