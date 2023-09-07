!=============================================================================
!     September 2020
!=============================================================================
function dotproduct(a,b)
!------------------------------------------------------------------------------
! This function is to compute the inner product between 2 vectors in 3d
! the input:
!           a =(a1,a2,a3) 
!           b = (b1,b2,b3)
! Output:
!           dotproduct=a.b = a1*b1 +a2*b3 +a3*b3
!-----------------------------------------------------------------------------
 
!  implicit none
  double precision, dimension(3), intent(in) :: a,b
  double precision :: dotproduct

  dotproduct = a(1)*b(1) + a(2)*b(2) +a(3)*b(3)
  
end function dotproduct
 
