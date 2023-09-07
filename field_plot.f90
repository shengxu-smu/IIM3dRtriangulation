!-----------------------------------------------------------------
subroutine field_plot
  !---------------------------------------------------------------
  use para
  use field
  !---------------------------------------------------------------
  integer:: i0,j0,k0  ! section
  ! gridlines c
  open(unit=81,file='DAT/xc.dat',status='unknown')
  open(unit=82,file='DAT/yc.dat',status='unknown')
  open(unit=83,file='DAT/zc.dat',status='unknown')
  do i=1,nx
     write(81,100)xc(i)
  enddo
  do j=1,ny
     write(82,100)yc(j)
  enddo
  do k=1,nz
     write(83,100)zc(k)
  enddo
  close(81)
  close(82)
  close(83)
  !--------------------------------------------------------------
  call vorticity
  !--------------------------------------------------------------
  i0=int(nx/2)!2-x0/dx
  j0=int(ny/2)!2-y0/dy
  k0=int(nz/2)!2-z0/dz
  !--------------------------------------------------------------
  open(unit=84,file='DAT/xu.dat',status='unknown')
  open(unit=85,file='DAT/yu.dat',status='unknown')
  open(unit=86,file='DAT/zu.dat',status='unknown')
  
  open(unit=87,file='DAT/xv.dat',status='unknown')
  open(unit=88,file='DAT/yv.dat',status='unknown')
  open(unit=89,file='DAT/zv.dat',status='unknown')
  
  open(unit=90,file='DAT/xw.dat',status='unknown')
  open(unit=91,file='DAT/yw.dat',status='unknown')
  open(unit=92,file='DAT/zw.dat',status='unknown')
  
  open(unit=93,file='DAT/xp.dat',status='unknown')
  open(unit=94,file='DAT/yp.dat',status='unknown')
  open(unit=95,file='DAT/zp.dat',status='unknown')

  
  do i=1,nx
     write(84,100) data1(i,j0,k0)
     write(87,100) data4(i,j0,k0)
     write(90,100) data7(i,j0,k0)
     write(93,100) p(i,j0,k0)
  enddo
  do j=1,ny
     write(85,100) data1(i0,j,k0)
     write(88,100) data4(i0,j,k0)
     write(91,100) data7(i0,j,k0)
     write(94,100)p(i0,j,k0)
  enddo
  do k=1,nz
     write(86,100) data1(i0,j0,k)
     write(89,100) data4(i0,j0,k)
     write(92,100) data7(i0,j0,k)
     write(95,100) p(i0,j0,k)
  enddo
  do i=84,95
     close(i)
  enddo
  
100 format(1x,e16.6e4)
!-------------------------------------------------
  open(unit=101,file='DAT/u.dat',status='unknown')
  open(unit=102,file='DAT/v.dat',status='unknown')
  open(unit=103,file='DAT/w.dat',status='unknown')
  open(unit=104,file='DAT/o1.dat',status='unknown')
  open(unit=105,file='DAT/o2.dat',status='unknown')
  open(unit=106,file='DAT/o3.dat',status='unknown')
  open(unit=107,file='DAT/p.dat',status='unknown')
  open(unit=108,file='DAT/d.dat',status='unknown')
  open(unit=109,file='DAT/q.dat',status='unknown')
  do k=1,nz
     do j=1,ny

        write(101,200)(data1(i,j,k),i=1,nx)
        write(102,200)(data4(i,j,k),i=1,nx)
        write(103,200)(data7(i,j,k),i=1,nx)

        write(104,200)(fecc(i,j,k),i=1,nx)
        write(105,200)(fcec(i,j,k),i=1,nx)
        write(106,200)(fcce(i,j,k),i=1,nx)

        write(107,200)(p(i,j,k),i=1,nx)
        write(108,200)(d(i,j,k),i=1,nx)
        write(109,200)(o(i,j,k),i=1,nx)
     enddo
  enddo
  do i=101,109
     close(i)
  enddo
  
200 format(1x,1000e16.6e4)
!------------------------------------------------------------------
  
  return
end subroutine field_plot


!-----------------------------------------------------------------------
