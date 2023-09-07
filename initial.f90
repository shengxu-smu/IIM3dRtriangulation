!-----------------------------------------------------------------------
subroutine initial
  use para
  use Lagrange
  use field
  use vfft
  
  integer :: n,m,ms,iv,nA,nB,nC
  !----------------------------------------------
  call vcosti(nx,wsavei)
  call vcosti(ny,wsavej)
  call vcosti(nz,wsavek)
  !-----------------------------------------------
  t=0.0d0
  t0=0.0d0
  !-----------------------------------------------
  call surface_initialize 
  
  open(unit=24,file='DAT/xs0.dat',status='unknown')
  open(unit=25,file='DAT/ys0.dat',status='unknown')
  open(unit=26,file='DAT/zs0.dat',status='unknown')
  DO m=1,nobj                                                              
     do iv=nv4obj(m-1)+1,nv4obj(m)
        write(24,*) xs(iv,m)                                       
        write(25,*) ys(iv,m)                                                  
        write(26,*) zs(iv,m)                                         
     enddo
  ENDDO
  close(24)
  close(25)
  close(26)
  !-----------------------------
  call field_initial
  !---------------------------------
  call velocity_reset
  
100 format(1x,1000e16.6e4)

  return
end subroutine initial
!-----------------------------------------------------------------------
