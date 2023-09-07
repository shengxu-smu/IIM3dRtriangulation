!===========================================================
! November 2020
!=============================================================
function coord4intersection(A,B,C,gridline,i) result(coord3)
  !------------------------------------------------------------
  ! Description:
  !             This function computes the 3rd coordinate of the
  !             intersection point between the gridline and the
  !             panel ABC
  !
  !             For a panel ABC,and a point p such that                                      
  !             A =(x1,y1,z1), B = (x2,y2,z2), C = (x3,y3,z3))                               
  !             p = (xint,yint,zint)                                                         
  !             then use the dot product:                                                    
  !             vec{AP}.vec{N}= n1(xint-x1)+n2(yint-y1)+n3(zint-z1)=0                        
  !             where N=(n1,n2,n3) the normal vector to the panel ABC 
  ! Input:
  !       A,B,C      - the coordinates of the vertices ABC  
  !                    Each one consists of 3 entries corresponding to
  !                    x,y,and z coordinates
  !
  !       gridline   - the projected gridline coordinates into xy,
  !                    xz, or yz plane. It has 2 entries depending
  !                    on which plane it was projected into
  !
  !       i          - this is a flag input to determine the direction
  !                    of the gridline
  !                    direction i=1:
  !                             gridline parallel to x-axis-> compute xint                                              
  !                    direction i=2:
  !                             gridline parallel to y-axis-> compute yint                                               
  !                    direction i=3:
  !                             gridline parallel to z-axis-. compute zint 
  !
  ! Output:
  !       coord3     - the coordinate of xint,yint,or zint depending on i
  !-----------------------------------------------------------------------
  implicit none
  ! declaration
  double precision:: coord3
  double precision, dimension(1:3), intent(in):: A,B,C
  double precision, dimension(1:2), intent(in):: gridline
  integer, intent(in):: i
  integer:: m,n
  double precision, dimension(1:3):: normal
  double precision:: eps
  !------------------------------------------------------------------------
  interface
     function crossproduct(a,b)
       double precision, dimension(1:3) :: crossproduct
       double precision, dimension(1:3),intent(in) :: a, b
     end function crossproduct
  end interface
  !-----------------------------------------------------------------------
  eps=1.0d-12
  ! compute the normal vector to the panel
  normal=crossproduct(B-A,C-A)
  normal=normal/dsqrt(dot_product(normal,normal))
  ! select which direction
  select case(i)
  case(1)
     m=2
     n=3
  case(2)
     m=1
     n=3
  case(3)
     m=1
     n=2
  end select
  !---------------------------------------------------
  if(dabs(normal(i))>eps) then 
     coord3=A(i)-(normal(m)*(gridline(1)-A(m))+normal(n)*(gridline(2)-A(n)))/(normal(i))
  else
     coord3=(A(i)+B(i)+C(i))/3.0d0
  end if
  !---------------------------------------------------------------------------------------
end function coord4intersection
!------------------------------------------------------------------------------------------
