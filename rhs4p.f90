subroutine rhs4p(A,B,C,M1,M2,M3,&
     dpdxJC_A,dpdyJC_A,dpdzJC_A,dpdxJC_B,dpdyJC_B,dpdzJC_B,&
     dpdxJC_C,dpdyJC_C,dpdzJC_C,dpdxJC_M1,dpdyJC_M1,dpdzJC_M1,&
     dpdxJC_M2,dpdyJC_M2,dpdzJC_M2,dpdxJC_M3,dpdyJC_M3,dpdzJC_M3,&
     rhs_A,rhs_B,rhs_C)


  implicit none
  
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
  end interface
  !============
  !--A-------------
  tau = B-A
  l1 = sqrt(dotproduct(tau,tau))
  tau = tau/l1
  
  beta = C-A
  l2 = sqrt(dotproduct(beta,beta))
  beta = beta/l2
  
  dpdt_A  = dotproduct( (/dpdxJC_A,dpdyJC_A,dpdzJC_A/) ,tau)
  dpdt_B  = dotproduct( (/dpdxJC_B,dpdyJC_B,dpdzJC_B/) ,tau)
  dpdt_M1 = dotproduct( (/dpdxJC_M1,dpdyJC_M1,dpdzJC_M1/),tau)

  dpdb_A  = dotproduct( (/dpdxJC_A,dpdyJC_A,dpdzJC_A/) ,beta)
  dpdb_C  = dotproduct( (/dpdxJC_C,dpdyJC_C,dpdzJC_C/) ,beta)
  dpdb_M3 = dotproduct( (/dpdxJC_M3,dpdyJC_M3,dpdzJC_M3/),beta)

  rhs1 = (l1/6.0d0)*(dpdt_A + 4.0d0*dpdt_M1 + dpdt_B )
  rhs2 = (l2/6.0d0)*(dpdb_A + 4.0d0*dpdb_M3 + dpdb_C )
  rhs_A = rhs1 + rhs2             
  !print*,rhs1,rhs2,rhs_A
  !------------
  !----B-----------------
  tau  = C-B
  l1   = sqrt(dotproduct(tau,tau))
  tau  = tau/l1

  beta = A-B
  l2 = sqrt(dotproduct(beta,beta))
  beta = beta/l2


  dpdt_C  = dotproduct( (/dpdxJC_C,dpdyJC_C,dpdzJC_C/) ,tau)
  dpdt_B  = dotproduct( (/dpdxJC_B,dpdyJC_B,dpdzJC_B/) ,tau)
  dpdt_M2 = dotproduct( (/dpdxJC_M2,dpdyJC_M2,dpdzJC_M2/),tau)

  dpdb_B  = dotproduct( (/dpdxJC_B,dpdyJC_B,dpdzJC_B/) ,beta)
  dpdb_A  = dotproduct( (/dpdxJC_A,dpdyJC_A,dpdzJC_A/) ,beta)
  dpdb_M1 = dotproduct( (/dpdxJC_M1,dpdyJC_M1,dpdzJC_M1/),beta)

  rhs1 = (l1/6.0d0)*(dpdt_C + 4.0d0*dpdt_M2 + dpdt_B )
  rhs2 = (l2/6.0d0)*(dpdb_B + 4.0d0*dpdb_M1 + dpdb_A )
  rhs_B = rhs1 + rhs2

  !----
  !--C----------------
  tau = A-C
  l1 = sqrt(dotproduct(tau,tau))
  tau = tau/l1

  beta = B-C
  l2 = sqrt(dotproduct(beta,beta))
  beta = beta/l2

  dpdt_A  = dotproduct( (/dpdxJC_A,dpdyJC_A,dpdzJC_A/) ,tau)
  dpdt_C  = dotproduct( (/dpdxJC_C,dpdyJC_C,dpdzJC_C/) ,tau)
  dpdt_M3 = dotproduct( (/dpdxJC_M3,dpdyJC_M3,dpdzJC_M3/),tau)

  dpdb_B  = dotproduct( (/dpdxJC_B,dpdyJC_B,dpdzJC_B/) ,beta)
  dpdb_C  = dotproduct( (/dpdxJC_C,dpdyJC_C,dpdzJC_C/) ,beta)
  dpdb_M2 = dotproduct( (/dpdxJC_M2,dpdyJC_M2,dpdzJC_M2/),beta)

  rhs1 = (l1/6.0d0)*(dpdt_A + 4.0d0*dpdt_M3 + dpdt_C )
  rhs2 = (l2/6.0d0)*(dpdb_B + 4.0d0*dpdb_M2 + dpdb_C )
  rhs_C = rhs1 + rhs2

  return
  
  
end subroutine rhs4p
