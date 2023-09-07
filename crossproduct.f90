!==========================================================================
! Summer 2020
!===========================================================================
function crossproduct(a,b) 
  !-------------------------------------------------------------------------
  !   Description:     
  !     This function is to compute the cross product of 2 vecters a and b 
  !     Input:
  !           a and b:  3d arrays a(a1,a2,a3) and b(b1,b2,b3)
  !     Output:
  !            crossproduct  =( a2*b3-a3*b2 , a3*b1-a1*b3 , a1*b2-a2*b1 ) 
  !--------------------------------------------------------------------------

  !--------------Declaration--------------------------------
  implicit none
  double precision, dimension(3), intent(in) :: a,b
  double precision, dimension(3) :: crossproduct
  !---------------------------------------------------------
  ! find the cross product 
  crossproduct(1) = a(2)*b(3) - a(3)*b(2)
  crossproduct(2) = a(3)*b(1) - a(1)*b(3)
  crossproduct(3) = a(1)*b(2) - a(2)*b(1) 

  return
  
end function  crossproduct
