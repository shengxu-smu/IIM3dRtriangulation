! -*- Mode: Fortran90; -*-
!--------------------------------------------------------------------------
! Sheng Xu (original author)
! Department of Mathematics
! SMU
! (modified and commented by Dale Pearson on August 24, 2012)
!==========================================================================

subroutine c2_solver(a,b,c,x,y,z,r,icheck,v)

!--------------------------------------------------------------------------
! Subroutine Description: This subroutine will solve a overdetermined 7x6
! linear system C_2*q = r analytically via Gaussian Elimination. This 
! technique is discussed in [1]. This program will start at equation (99)
! on page 2085.
! 
! Inputs:
! a,b,c,x,y,z - the components of tau and beta. the values of a,b and c take
! on the tau values. a,b and c could be either tau_1, tau_2 or tau_3 and
! similarly x,y or z could be either b_1, b_2 or b_3 (see pg. 2085 in [1]). 
! These are double precision (real*8) scalars.  
! 
! r - the right hand side of the equation. This will be a 1D array of size 7.
!
! icheck - a flag that ensures the user consistency among the variables. If 
!          icheck = 1, then the subroutine will output the extra equation and 
!          the appropriate variable to the screen. Note: If you do not want 
!          to use this feature, then input any integer not equal to 1.
!        
! Output:
!
! v - the solution array. In this case this will be the six second order jump
!     conditions. This will be a 1D array of size 6.
!
! References
!
! [1] Sheng Xu, Z. Jane Wang, A 3D immersed interface method for fluid-solid
!     interaction, Compt. Methods Appl Mech Engrg. 197 (2008) 2068-2086. 
!---------------------------------------------------------------------------

!======= Inclusions ===========

!======= Declarations =========
  implicit none
  integer, intent(in) :: icheck
  double precision :: a,b,c,x,y,z,r(7),v(6)
  double precision :: s2,s3,m1,m2,m3,e,f,g,h

  ! input protect the program from erroneous input by stopping the program 
  ! (in this case a)
  if(a.eq.0.0d0) then
     write(*,*)' !! wrong arrangement in matrix !!'
     stop
  endif
 
  ! declare the variables that we will need to perform Gaussian elimination.
  s2 = a*a+c*c
  s3 = a*a+b*b
  m1 = b*z-c*y
  m2 = c*x-a*z
  m3 = a*y-b*x
  e = 2.d0*b*c*(a*x+b*y)-s3*(c*y+b*z)
  f = s2*(a*x+b*y)-s3*(a*x+c*z)
  g = -s3*m2-2.0d0*b*c*m3
  h = -s2*m3

  ! test if the tau and beta vectors are parallel, if they are then stop the program
  if((m1*m1+m2*m2+m3*m3).eq.0.0d0) then
     write(*,*)' !! parallel vectors !!'
     stop
  endif
  
  ! perform Gaussian elimination on (99):   
  
  ! 1st column
  r(2)=r(2)-a*r(1)
  r(3)=r(3)-x*r(1)

  ! 2nd column
  r(2)=a*r(2)-b*r(4)
  r(3)=a*r(3)-y*r(4)
  r(5)=a*r(5)-x*r(4)
  
  ! 3rd column
  r(2)=r(2)-c*r(6)
  r(3)=r(3)-z*r(6)
  r(7)=a*r(7)-x*r(6)

  ! 4th column
  r(3)=s3*r(3)-(a*x+b*y)*r(2)
  r(5)=s3*r(5)+m3*r(2)

  ! find the pivot by determining if m3 or m2 is greater. 
  ! If m3 >= m2, then we first find v(6) and subsequently 
  ! find v(5). This part also contains the code for icheck
  ! (explained above)
  if(abs(m3).ge.abs(m2)) then
     r(3)=m3*r(3)-e*r(7)
     r(5)=m3*r(5)-g*r(7)
     v(6)=r(5)/(m3*h+m2*g)
     v(5)=(r(7)+m2*v(6))/m3
     if(icheck.eq.1) then
        if((m3*f+m2*e).ne.0.0d0) then
           write(*,*)v(6),r(5)/(m3*h+m2*g)
        else
           write(*,*)m3*f+m2*e,r(3)
        endif
     endif
     ! If m3 < m2, then we first find v(5) and subsequently 
     ! find v(6).
  else
     r(3)=m2*r(3)+f*r(7)
     r(5)=m2*r(5)+h*r(7)
     v(5)=r(5)/(m3*h+m2*g)
     v(6)=-(r(7)-m3*v(5))/m2
     if(icheck.eq.1) then
        if((m3*f+m2*e).ne.0.0d0) then
           write(*,*)v(5),r(5)/(m3*h+m2*g)
        else
           write(*,*)m3*f+m2*e,r(3)
        endif
     endif
  endif
  
  ! back solve the rest of the variables (v4 to v1) 
  v(4)=-(r(2)+2.0d0*b*c*v(5)+s2*v(6))/s3
  v(3)=(r(6)-b*v(5)-c*v(6))/a
  v(2)=(r(4)-b*v(4)-c*v(5))/a
  v(1)=r(1)-v(4)-v(6)

  return ! exits out the subroutine
end subroutine c2_solver
