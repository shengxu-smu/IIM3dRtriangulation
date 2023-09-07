  !================================================================
  ! March 2021
  !=================================================================
subroutine onesidedFD(su,uu,ds,nn,flag,du)
  !----------------------------------------------------------------
  ! Description
  !         one sided finite difference
  ! Input:
  !      su - su(3) velocity at vertex (u,v,w)  
  !      uu - u(nn,3) velocity at pints s1,s2,..
  !                  columns:  velocity components u,v,w
  !                  rows: velocity for points
  !      nn - number of points in the n direction
  !      flag -order of derivative
  !            1 fist derivative
  !            2 second derivative
  !      ds - the distance ds =sqrt(dx^2+dy^2+dz^2)
  ! Output:
  !       du = du/dn|+   OR  d2u/dn2|+
  !---------------------------------------------------------------
  ! DECLARATION:
  implicit none
  integer, intent(in) :: nn,flag
  double precision, dimension(nn,3), intent(in) :: uu
  double precision, dimension(3),intent(in) :: su
  double precision, intent(in):: ds
  double precision, dimension(3), intent(out) :: du
  double precision, dimension(3) :: u1,u2,u3,u4
  integer :: i
  !--------------------------------------------------------------
  if (flag.eq.1) then
     if(nn.eq.1)then  ! du/dn = (u(s1)-u(s0))/dsn
        u1 = uu(1,:)
        du = (u1-su)/ds   
     endif
     if(nn.eq.2)then  !  du/dn=(4u(s1)-u(s2)-3u(s0))/2dsn
        u1 = uu(1,:)
        u2 = uu(2,:)
        du = (-3.0d0*su + 4.0d0*u1 - u2)/(2.0d0*ds)
     endif
     if(nn.eq.3)then !If 3 pts,du/dn=(18u(s1)-9u(s2)+2u(s3)-11u(s0))/6dsn    
        u1 = uu(1,:)
        u2 = uu(2,:)
        u3 = uu(3,:)
        du = (18.0d0*u1 - 9.0d0*u2 + 2.0d0*u3 - 11.0d0*su)/(6.0d0*ds)
     endif
  elseif(flag.eq.2)then
     if(nn.eq.1)then       
        ! error
     endif
     if(nn.eq.2)then  ! d2u/dnn  O(ds)                     
        u1 = uu(1,:)
        u2 = uu(2,:)
        du = ( su  -2.0d0*u1 + u2)/(ds*ds)
     endif
     if(nn.eq.3)then !  d2u/dnn   O(ds^2)
        u1 = uu(1,:)
        u2 = uu(2,:)
        u3 = uu(3,:)
        du = ( 2.0d0*su -5.0d0*u1 + 4.0d0*u2 -u3)/(ds*ds)
     endif
     if(nn.eq.4)then !  d2u/dnn   O(ds^3)                    
        u1 = uu(1,:)
        u2 = uu(2,:)
        u3 = uu(3,:)
        u4 = uu(4,:)
        du = ( 35.0d0*su -104.0d0*u1 + 114.0d0*u2 -56.0d0*u3 +11.0d0*u4)/(12.0d0*ds*ds)
     endif
  endif
  
  
end subroutine onesidedFD
