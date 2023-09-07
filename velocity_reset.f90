!-----------------------------------------------------------------------
! update the velocity inside the object
subroutine velocity_reset
  use para
  use field
  use Lagrange
  !---------------------------------------------------------------------
  integer :: myobj
  double precision :: s1,s2,s3,c1,c2,c3
  !---------------------------------------------------------------------
  DO m=1,nobj
     myobj = objects(m)
     !------------------------------------
     !u
     do i = 0,nx+1
        do j = 0,ny+1
           do k = 0,nz+1
              if (iou(i,j,k).eq.myobj) then
                 u(i,j,k) = xsct(m)+omegay(m)*(zc(k)-zsc(m))&
                                   -omegaz(m)*(yc(j)-ysc(m))
              endif
           enddo
        enddo
     enddo
     !----------------------------------------------------------
     !v
     do i = 0,nx+1
        do j = 0,ny+1
           do k = 0,nz+1       
              if (iov(i,j,k).eq.myobj) then
                 v(i,j,k) = ysct(m)-omegax(m)*(zc(k)-zsc(m)) &
                                   +omegaz(m)*(xc(i)-xsc(m))
              endif
           enddo
        enddo
     enddo
     !-------------------------------------------------------------
     !w
     do i = 0,nx+1
        do j = 0,ny+1
           do k = 0,nz+1
              if (iow(i,j,k).eq.myobj) then
                 w(i,j,k) = zsct(m)+omegax(m)*(yc(j)-ysc(m))&
                                   -omegay(m)*(xc(i)-xsc(m))
              endif
           enddo
        enddo
     enddo
     !-----------------------------------------------------------------
  ENDDO

  return
  !---------------------------------------------------------------------  
end subroutine velocity_reset
!-----------------------------------------------------------------------
