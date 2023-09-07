!-----------------------------------------------------------------------
subroutine jc_pressure
 
  use para
  use Lagrange

  double precision::signx,signy,signz
  integer :: n,m,ms
  ms =1
  !--------
  DO m=1,ms
     
     do n=1,nicjc
        ! the sign of the normal vector
        signx=sign(1.0d0,ficjc(10,n,m))
        signy=sign(1.0d0,ficjc(11,n,m))
        signz=sign(1.0d0,ficjc(12,n,m))      

        ficjc(22,n,m)=signx*ficjc(22,n,m)
        ficjc(23,n,m)=signy*ficjc(23,n,m)
        ficjc(24,n,m)=signz*ficjc(24,n,m)

        ficjc(43,n,m)=signx*ficjc(43,n,m)
        ficjc(46,n,m)=signy*ficjc(46,n,m)
        ficjc(48,n,m)=signz*ficjc(48,n,m)

     enddo

     do n=1,nickc        
        signx=sign(1.0d0,fickc(10,n,m))
        signy=sign(1.0d0,fickc(11,n,m))
        signz=sign(1.0d0,fickc(12,n,m))

        fickc(22,n,m)=signx*fickc(22,n,m)
        fickc(23,n,m)=signy*fickc(23,n,m)
        fickc(24,n,m)=signz*fickc(24,n,m)

        ficKc(43,n,m)=signx*ficKc(43,n,m)
        ficKc(46,n,m)=signy*ficKc(46,n,m)
        ficKc(48,n,m)=signz*ficKc(48,n,m)
     enddo

     do n=1,njckc
        
        signx=sign(1.0d0,fjckc(10,n,m))
        signy=sign(1.0d0,fjckc(11,n,m))
        signz=sign(1.0d0,fjckc(12,n,m))
        
        fjckc(22,n,m)=signx*fjckc(22,n,m)
        fjckc(23,n,m)=signy*fjckc(23,n,m)
        fjckc(24,n,m)=signz*fjckc(24,n,m)
        
        fjckc(43,n,m)=signx*fjckc(43,n,m)
        fjckc(46,n,m)=signy*fjckc(46,n,m)
        fjckc(48,n,m)=signz*fjckc(48,n,m)

        !if(fjckc(22,n,m).ne.0.0d0)  print*,fjckc(43:48,:,:)
     enddo
     
  ENDDO

  !print*,fjckc(43:48,:,:)
  !pause
  return
end subroutine jc_pressure
!-----------------------------------------------------------------------
