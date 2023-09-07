!-----------------------------------------------------------------------

subroutine u_strain
  use para
  use field
  use Lagrange
  
  integer :: ie,je,ke
  integer :: n,m,ms
  ms =1
  !---------
  do k=1,nz
     do j=1,ny
        do i=1,nx
           data1(i,j,k)=dx1*(fecc(i,j,k)-fecc(i-1,j,k))
           data2(i,j,k)=dy1*(fcec(i,j,k)-fcec(i,j-1,k))
           data3(i,j,k)=dz1*(fcce(i,j,k)-fcce(i,j,k-1))
        enddo
     enddo
  enddo
  
  IF(isingular.eq.1) THEN
     DO m=1,ms

        do n=1,nicjc
           i=icjc(1,n,m)
           j=icjc(2,n,m)
           ke=icjc(3,n,m)
           k=ke+1
           data3(i,j,k)=data3(i,j,k)+udz(n,m)
        enddo

        do n=1,nickc
           i=ickc(1,n,m)
           k=ickc(2,n,m)
           je=ickc(3,n,m)
           j=je+1
           data2(i,j,k)=data2(i,j,k)+udy(n,m)
        enddo
        do n=1,njckc
           j=jckc(1,n,m)
           k=jckc(2,n,m)
           ie=jckc(3,n,m)
           i=ie+1
           data1(i,j,k)=data1(i,j,k)+udx(n,m)
        enddo

        
     ENDDO
  ENDIF
  
  !---------
  
  return
end subroutine u_strain


!-----------------------------------------------------------------------



