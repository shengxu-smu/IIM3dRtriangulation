subroutine field_initial
  !------------------------------------------
  use para
  use field
  !------------------------------------------
  do k=0,nz+1
     do j=0,ny+1
        do i=0,nx
           u(i,j,k)=1.0d0
        enddo
     enddo
  enddo
  !
  do k=0,nz+1
     do j=0,ny
        do i=0,nx+1
           v(i,j,k)=0.0d0
        enddo
     enddo
  enddo
  !
  do k=0,nz
     do j=0,ny+1
        do i=0,nx+1
           w(i,j,k)=0.0d0
        enddo
     enddo
  enddo
  !
  do k=0,nz+1
     do j=0,ny+1
        do i=0,nx+1
           p(i,j,k)=0.0d0
        enddo
     enddo
  enddo
  !
  do k=0,nz+1
     do j=0,ny+1
        do i=0,nx+1
           d(i,j,k)=0.0d0
        enddo
     enddo
  enddo
  !
  return
  !
end subroutine field_initial
!-------------------------------------------------
