subroutine CartJump2nd4u(A,B,C,n,dudnnPJC_A,dudxJC_A,dudyJC_A,dudzJC_A,&
     dudxJC_B,dudyJC_B,dudzJC_B,dudxJC_C,dudyJC_C,dudzJC_C,&
     dudxxJC_A,dudxyJC_A,dudxzJC_A,dudyyJC_A,dudyzJC_A,dudzzJC_A)

  ! This subroutine computes [d2u/dxidxj] @ a vertex A in the panel ABC
  ! using [d2u/dnn] @ A, and [du/dxj] @ A,B,C
  !=======================================
  implicit none
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
  !
  tau  = B-A
  l1   = sqrt(dotproduct(tau,tau))
  tau  = tau/l1
  
  beta = C-A
  l2   = sqrt(dotproduct(beta,beta))
  beta = beta/l2

  !--u-----------------------------------------
  rhs(1) = (dudxJC_B(1)-dudxJC_A(1))/l1
  rhs(2) = (dudyJC_B(1)-dudyJC_A(1))/l1                
  rhs(3) = (dudzJC_B(1)-dudzJC_A(1))/l1                
  rhs(4) = (dudxJC_C(1)-dudxJC_A(1))/l2                 
  rhs(5) = (dudyJC_C(1)-dudyJC_A(1))/l2               
  rhs(6) = (dudzJC_C(1)-dudzJC_A(1))/l2                                      
  rhs(7) = dudnnPJC_A(1)
  call  cjc4d2u(tau,beta,rhs,jumpd2u)
  !du2ndJC_A = jumpd2u
  dudxxJC_A(1)=jumpd2u(1)
  dudxyJC_A(1)=jumpd2u(2)
  dudxzJC_A(1)=jumpd2u(3)
  dudyyJC_A(1)=jumpd2u(4)
  dudyzJC_A(1)=jumpd2u(5)
  dudzzJC_A(1)=jumpd2u(6)
  !----v--------------------
  rhs(1) = (dudxJC_B(2)-dudxJC_A(2))/l1                                       
  rhs(2) = (dudyJC_B(2)-dudyJC_A(2))/l1
  rhs(3) = (dudzJC_B(2)-dudzJC_A(2))/l1
  rhs(4) = (dudxJC_C(2)-dudxJC_A(2))/l2
  rhs(5) = (dudyJC_C(2)-dudyJC_A(2))/l2
  rhs(6) = (dudzJC_C(2)-dudzJC_A(2))/l2
  rhs(7) = dudnnPJC_A(2)
  call  cjc4d2u(tau,beta,rhs,jumpd2u)
  dudxxJC_A(2)=jumpd2u(1)
  dudxyJC_A(2)=jumpd2u(2)
  dudxzJC_A(2)=jumpd2u(3)
  dudyyJC_A(2)=jumpd2u(4)
  dudyzJC_A(2)=jumpd2u(5)
  dudzzJC_A(2)=jumpd2u(6)
  !--w----------------------------------
  rhs(1) = (dudxJC_B(3)-dudxJC_A(3))/l1
  rhs(2) = (dudyJC_B(3)-dudyJC_A(3))/l1
  rhs(3) = (dudzJC_B(3)-dudzJC_A(3))/l1
  rhs(4) = (dudxJC_C(3)-dudxJC_A(3))/l2
  rhs(5) = (dudyJC_C(3)-dudyJC_A(3))/l2
  rhs(6) = (dudzJC_C(3)-dudzJC_A(3))/l2
  rhs(7) = dudnnPJC_A(3)
  call  cjc4d2u(tau,beta,rhs,jumpd2u)
  dudxxJC_A(3)=jumpd2u(1)
  dudxyJC_A(3)=jumpd2u(2)
  dudxzJC_A(3)=jumpd2u(3)
  dudyyJC_A(3)=jumpd2u(4)
  dudyzJC_A(3)=jumpd2u(5)
  dudzzJC_A(3)=jumpd2u(6)
  !-------------------------------
  return
  
end subroutine CartJump2nd4u
!===================================================
