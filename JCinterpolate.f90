!-----------------------------------------------------------------
! May 3rd
!----------------------------------------------------------------

subroutine JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
  !'''''''''''''''''''''''''''''''''''''''''''''''''''
  ! This subroutine is to interpolate the jump condition from the vertices of the triangle ABC
  ! to the intersection point Q between the cartesian grid line and the triangle panel
  !
  ! Input:
  !       A - the coordinate (x,y,z) of the vertex A
  !       B - the coordinate (x,y,z) of the vertex B
  !       C - the coordinate (x,y,z) of the vertex C
  !       Q - the coordinate (x,y,z) of the intersection point Q
  !       fA,fB,fC - the jump condition of the 1st or 2nd derivatives at the vertices in the direction x,y,or z 
  ! Output:
  !       fQ - the jump condition at Q
  !''''''''''''''''''''''''''''''''''''''''''''''''''''''
  !======= Declarations =========
  implicit none
  interface
     function crossproduct(a,b)
       double precision, dimension(3), intent(in) :: a,b
       double precision, dimension(3) :: crossproduct
     end function crossproduct
     function dotproduct(a,b)
        double precision, dimension(3), intent(in) :: a,b
        double precision :: dotproduct
      end function dotproduct
  end interface
  
  integer :: i

  double precision, dimension(1:3), intent(in) :: A,B,C,Q
  double precision, intent(in) :: fA,fB,fC
  double precision, intent(out) :: fQ

  double precision, dimension(1:3) :: vec1,vec2,vec3,vec4
  double precision :: normvec1, normvec2, normvec3, normvec4

  !======= Internals =========
  ! find vec1 = ||vec{QB} x vec{QC}|| by calling the crossprod function
  vec1 =  crossproduct(B-Q,C-Q)
  normvec1 = sqrt(dotproduct(vec1,vec1))
  ! find vec2 = ||vec{AQ} x vec{AC}|| by calling the crossprod function
  vec2 =  crossproduct(Q-A,C-A)
  normvec2 = sqrt(dotproduct(vec2,vec2))
  ! find vec3 = ||vec{AB} x vec{AQ}|| by calling the crossprod function
  vec3 =  crossproduct(B-A,Q-A)
  normvec3 = sqrt(dotproduct(vec3,vec3))
  ! find vec4 = ||vec{AB} x vec{AC}|| by calling the crossprod function.                                        
  vec4 =  crossproduct(B-A,C-A)
  normvec4 = sqrt(dotproduct(vec4,vec4))

  !find the [f] at Q  through interpolation
  fQ = (1.d0/normvec4)*(normvec1*fA+normvec2*fB+normvec3*fC)


  return
end subroutine JCinterpolate
