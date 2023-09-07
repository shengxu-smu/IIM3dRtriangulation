!---------------------------------------------------------
subroutine surface_plot
  use para
  use Lagrange
  integer :: m,ms,iv,ip,jp
  ms =1

  !--------
  write(*,*)
  write(*,*)'even numbers:'
  write(*,*)'niejc = ',niejc
  write(*,*)'nicje = ',nicje
  write(*,*)'nicjc = ',nicjc
  write(*,*)'niekc = ',niekc
  write(*,*)'nicke = ',nicke
  write(*,*)'nickc = ',nickc
  write(*,*)'njekc = ',njekc
  write(*,*)'njcke = ',njcke
  write(*,*)'njckc = ',njckc

  open(unit=71,file='DAT/xs.dat',status='unknown')
  open(unit=72,file='DAT/ys.dat',status='unknown')
  open(unit=73,file='DAT/zs.dat',status='unknown')
  open(unit=74,file='DAT/us.dat',status='unknown')
  open(unit=75,file='DAT/vs.dat',status='unknown')
  open(unit=76,file='DAT/ws.dat',status='unknown')
  open(unit=77,file='DAT/nv.dat',status='unknown')
  open(unit=78,file='DAT/np.dat',status='unknown')
  open(unit=79,file='DAT/panel.dat',status='unknown')
  write(77,*) nv4obj(0)
  write(78,*) np4obj(0)
  DO m=1,ms
     write(77,*) nv4obj(m)
     do iv=1+nv4obj(m-1),nv4obj(m)

        write(71,*) xs(iv,m)
        write(72,*) ys(iv,m)
        write(73,*) zs(iv,m)
        write(74,*) us(iv,m)
        write(75,*) vs(iv,m)
        write(76,*) ws(iv,m)
     enddo
  ENDDO
  DO m=1,ms
     write(78,*) np4obj(m)
     do ip=1+np4obj(m-1),np4obj(m)
        write(79,*) (panel(jp,ip),jp=1,3)
     enddo
  ENDDO
  

  Do i=71,79
     close(i)
  ENDDo
  
200 format(1x,1000e16.6e4)
  
  return
end subroutine surface_plot


!c-----------------------------------------------------------------------
