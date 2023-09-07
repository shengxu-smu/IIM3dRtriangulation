! 
subroutine CG4p(nv1,nv2,np1,np2,rhs,xx)
  !------------------------------------------------------------------
  ! The conjugate gradient to solve the spd linear system  A*([p]) =(rhs)
  ! The original system is semi-definite, so we adjust the 1st row and column of the matrix 
  ! Input:
  !        nv1,nv2,np1,np2
  !               -nv1,nv2: the numbers of vertices
  !               -np1,np2: the numbers of panels
  !        rhs    - the rhs of the linear system for the pressure for
  !                 current object.
  !                 rhs(nv1+1:nv2)
  ! Output:
  !        xx     - the resulting jump condition [p] @ vertices
  !-----------------------------------------------------------
  !use
  use Lagrange
  !----------------------------------------------------------------------
  implicit none
  integer, intent(in):: nv1,nv2,np1,np2
  double precision, dimension(nv1+1:nv2),intent(in) :: rhs
  double precision, dimension(nv1+1:nv2),intent(out) :: xx
  integer:: i,iv,ip,nv,nA,nB,nC,iter
  double precision, dimension(nv1+1:nv2)::r,pp,ADiag,mrhs,zz,Ap,Ax
  double precision :: eps, err,dotrr,dotrz,alpha,beta
  interface
     function dotproduct2(A,B,nv1,nv2) 
       integer, intent(in)::nv1,nv2
       double precision, dimension(nv1+1:nv2), intent(in) :: A,B
       double precision :: dotproduct2
       integer:: i
       double precision::AdotB
     end function dotproduct2
  end interface
  !------------------------------------------------------
  !-----------------------------------------------------
  xx = pJC0(nv1+1:nv2) ! initial guess is [p] from the previous step
  !mrhs = rhs
  ! for Ax=b ---> r = b-Ax
  Ax = 0.0d0
  do ip = np1+1,np2
     nA = panel(1,ip)+nv1
     nB = panel(2,ip)+nv1
     nC = panel(3,ip)+nv1
     if (nA .eq.nv1+1) then
        Ax(nA) = Ax(nA) - 2.0d0*xx(nA)
        Ax(nB) = Ax(nB) - 2.0d0*xx(nB) + xx(nC)
        Ax(nC) = Ax(nC) - 2.0d0*xx(nC) + xx(nB)
     elseif(nB.eq.nv1+1) then
        Ax(nA) = Ax(nA) - 2.0d0*xx(nA) + xx(nC)
        Ax(nB) = Ax(nB) - 2.0d0*xx(nB)
        Ax(nC) = Ax(nC) - 2.0d0*xx(nC) + xx(nA)
     elseif (nC.eq.nv1+1) then
        Ax(nA) = Ax(nA) - 2.0d0*xx(nA) + xx(nB)
        Ax(nB) = Ax(nB) - 2.0d0*xx(nB) + xx(nA)
        Ax(nC) = Ax(nC) - 2.0d0*xx(nC)
     else
        Ax(nA) = Ax(nA) - 2.0d0*xx(nA) + xx(nB) + xx(nC)
        Ax(nB) = Ax(nB) - 2.0d0*xx(nB) + xx(nA) + xx(nC)
        Ax(nC) = Ax(nC) - 2.0d0*xx(nC) + xx(nA) + xx(nB)
     endif
  enddo
  mrhs = rhs-Ax
  mrhs(nv1+1) =0.0d0
  r = mrhs
  pp = r
  eps = 10e-8
  err = 100.d0
  iter = 0
  dotrr = dotproduct2(r,r,nv1,nv2)
  if (abs(dotrr).le.eps) then
     print*,'CG: ||r||_2 =  ',dotrr
     return
  endif
  do while (err > eps)
     iter = iter + 1
     dotrr = dotproduct2(r,r,nv1,nv2)    
     !------------------------------
     !Ap
     Ap = 0.0d0     
     do ip = np1+1,np2
        ! indices of the current panel's vertices
        nA = panel(1,ip)+nv1
        nB = panel(2,ip)+nv1
        nC = panel(3,ip)+nv1
        if (nA .eq.nv1+1) then
           Ap(nA) = Ap(nA) - 2.0d0*pp(nA)
           Ap(nB) = Ap(nB) - 2.0d0*pp(nB) + pp(nC)
           Ap(nC) = Ap(nC) - 2.0d0*pp(nC) + pp(nB)
        elseif(nB.eq.nv1+1) then
           Ap(nA) = Ap(nA) - 2.0d0*pp(nA) + pp(nC)
           Ap(nB) = Ap(nB) - 2.0d0*pp(nB) 
           Ap(nC) = Ap(nC) - 2.0d0*pp(nC) + pp(nA) 
        elseif (nC.eq.nv1+1) then
           Ap(nA) = Ap(nA) - 2.0d0*pp(nA) + pp(nB) 
           Ap(nB) = Ap(nB) - 2.0d0*pp(nB) + pp(nA) 
           Ap(nC) = Ap(nC) - 2.0d0*pp(nC) 
        else
           Ap(nA) = Ap(nA) - 2.0d0*pp(nA) + pp(nB) + pp(nC)
           Ap(nB) = Ap(nB) - 2.0d0*pp(nB) + pp(nA) + pp(nC)
           Ap(nC) = Ap(nC) - 2.0d0*pp(nC) + pp(nA) + pp(nB)
        endif
     enddo
     !-------------------------------------------------------
     alpha = dotrr/dotproduct2(pp,Ap,nv1,nv2)    
     xx = xx + alpha*pp
     r = r-alpha*Ap
     beta = dotproduct2(r,r,nv1,nv2)/dotrr
     pp = r + beta*pp
     err = dsqrt(dotproduct2(r,r,nv1,nv2))
     !---
  enddo
  !--------------------------------------------------------------------
  !print*,'niter',iter,'residual',err
  return
end subroutine CG4p
