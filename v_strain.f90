!-----------------------------------------------------------------------
!
subroutine v_strain
  use para
  use Lagrange
  use field
  
  integer:: ie,je,ke

  integer :: n,m,ms
  ms =1
  !---------------

  do k=1,nz
     do j=1,ny
        do i=1,nx
           data4(i,j,k)=dx1*(fecc(i,j,k)-fecc(i-1,j,k))
           data5(i,j,k)=dy1*(fcec(i,j,k)-fcec(i,j-1,k))
           data6(i,j,k)=dz1*(fcce(i,j,k)-fcce(i,j,k-1))
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
           data6(i,j,k)=data6(i,j,k)+vdz(n,m)
        enddo

        do n=1,nickc
           i=ickc(1,n,m)
           k=ickc(2,n,m)
           je=ickc(3,n,m)
           j=je+1
           data5(i,j,k)=data5(i,j,k)+vdy(n,m)
        enddo
        do n=1,njckc
           j=jckc(1,n,m)
           k=jckc(2,n,m)
           ie=jckc(3,n,m)
           i=ie+1
           data4(i,j,k)=data4(i,j,k)+vdx(n,m)
        enddo

     ENDDO
  ENDIF

  !---------------
  
  return
end subroutine v_strain


!-----------------------------------------------------------------------



