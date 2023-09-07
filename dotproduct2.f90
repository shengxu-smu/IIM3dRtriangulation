!=============================================================================
!     March 2022
!=============================================================================
function dotproduct2(A,B,nv1,nv2)
  !------------------------------------------------------------------------------
  ! This function is to compute the inner product between 2 vectors in 3d
  ! the input:
  !           a =(a1,a2,a3,....) 
  !           b = (b1,b2,b3,.....)
  ! Output:
  !           dotproduct=a.b = a1*b1 +a2*b3 +a3*b3+......
  !-----------------------------------------------------------------------------  
  !  implicit none
  integer, intent(in)::nv1,nv2
  double precision, dimension(nv1+1:nv2), intent(in) :: A,B
  double precision :: dotproduct2
  integer:: i
  double precision::AdotB
  !----
  AdotB =0.0d0
  do i=nv1+1,nv2
     AdotB =AdotB +A(i)*B(i)
  enddo
  dotproduct2 = AdotB
  return
end function dotproduct2
 
