!=============================================================================================
! August 2020
!============================================================================================
subroutine pointinside(a,b,c,p,s)
  !-------------------------------------------------------------------------------------------
  ! Description: This subroutine is to test whether a grid point P lies inside a triangle ABC or not
  !              using sameside test
  ! Input:
  !       a =(a1,a2,a3): the coordinates of the point A
  !       b = (b1,b2,b3): the coordinates of the point B
  !       c = (c1,c2,c3): the coordinates of the point C
  !       p = (p1,p2,p3): the coordinates of the point P
  ! Output:
  !       s: the result of the test:
  !          if s = 1 then P is located inside the traingle
  !          if s = 0 then P is located on an edge or a vertex of the triangle
  !          if s = -1 then P is located outside the triangle
  !-----------------------------------------------------------------------------------------
  !-----------------Declarations-------------------------
  implicit none
  integer :: proceed
  double precision, dimension(3),intent(in) :: a,b,c,p
  double precision, dimension(3) :: vec1,vec2,vec3,vec4,vec5,vec6
  integer, intent(out) :: s
  !---------------------------------------------------
  !interface 
  interface
     function crossproduct(a,b)
       double precision, dimension(3), intent(in)::a,b
       double precision, dimension(3) :: crossproduct
     end function crossproduct
     function dotproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision :: dotproduct
     end function dotproduct
  end interface

  !---------------------------------------------------
  ! initilize s to be -1
  s = -1
  !---------------------------------------------------------------
  ! First: test the side vec{BA}
  ! call the crossproduct function to compute:
  ! vec1 = vec(BP} x vec{BA}
  ! vec2 = vec{BC} x vec{BA}
  ! find the sign of the dot product of vec1 and vec2
  vec1 = crossproduct(p-b,a-b)
  vec2 = crossproduct(c-b,a-b)
  
  if ( dotproduct(vec1,vec2) >  0.0d0 ) then
     proceed = 1
  elseif ( dotproduct(vec1,vec2) .eq. 0.0d0 ) then
     s = 0
     proceed = 1
  else
     s = -1
     return
  endif  
  !---------------------------------------------------------------                                      
  ! Second: if the first test is successful, then test  side vec{AC}
  ! call the crossproduct function to compute:                                       
  ! vec3 = vec(AP} x vec{AC}
  ! vec4 = vec{AB} x vec{AC}
  ! Then find the sign of the dotproduct of vec3 and vec4
  
  if ( proceed .eq. 1 ) then 
     vec3 = crossproduct(p-a,c-a)
     vec4 = crossproduct(b-a,c-a)
     
     if ( dotproduct(vec3,vec4) >  0.0d0 ) then
        proceed = 1           
     elseif ( dotproduct(vec3,vec4) .eq. 0.0d0 ) then
        s = 0
        proceed = 1
     else
        s = -1
        return
     endif  
  endif
  !---------------------------------------------------------------
  ! Third: if the second test is successful, then test  side vec{CB}
  ! call the crossproduct function to compute:
  ! vec5 = vec(CP} x vec{CB}
  ! vec6 = vec{CA} x vec{CB}
  ! Then find the sign of the dotproduct of vec5 and vec6  
  if ( proceed .eq. 1) then 

     vec5 = crossproduct(p-c,b-c)
     vec6 = crossproduct(a-c,b-c)

     if ( ( dotproduct(vec5,vec6) >  0.0d0) .and. (s.ne.0) ) then
        s = 1
     elseif ( (dotproduct(vec5,vec6) > 0.0d0  ) .and. (s .eq.0)  ) then
        s = 0
     elseif ( ( dotproduct(vec5,vec6) .eq. 0.0d0)  ) then
        s = 0
     else
        s = -1
     endif
  endif

  return
 
end subroutine pointinside
