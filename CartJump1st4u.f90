subroutine CartJump1st4u(A,B,C,n,dudnPJC,dudxJC,dudyJC,dudzJC)
  !============================================================
  ! This subroutine computes the jump condition of du/dxj @ a vertex A in the panel ABC
  !===============================================================================
  implicit none
  double precision,dimension(3), intent(in) :: A,B,C,n,dudnPJC
  double precision,dimension(3), intent(out) :: dudxJC,dudyJC,dudzJC
  double precision,dimension(3):: tau,beta
  double precision :: l1,l2,den,numx,numy,numz,d1,d2,d3
  !----------------------------------------------------------------------------
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
  tau = B-A
  l1 = sqrt(dotproduct(tau,tau))
  tau = tau/l1

  beta = C-A
  l2 = sqrt(dotproduct(beta,beta))
  beta = beta/l2

  dudxJC(:) = n(1)*dudnPJC(:)
  dudyJC(:) = n(2)*dudnPJC(:)                                                
  dudzJC(:) = n(3)*dudnPJC(:)  

  return
  !==============================================================
end subroutine CartJump1st4u
!================================================================
