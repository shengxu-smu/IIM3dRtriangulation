!-----------------------------------------------------------------------
! Still working on this
!-----------------------------------------------------------------------
subroutine IrrP_JC2(A,B,C,tau,beta,n,A_JC,B_JC,C_JC,ip,iobj,Xcc,Ucc,Omegac,&
     Adup,Adum,Bdup,Bdum,Cdup,Cdum)
  !====================================================================
  ! For each panel:
  !              1. find the intersection points
  !                 9 type of intersection points
  !              2. interpolate the jump condition at this intersection points
  !              3. compute the jumps contribution to the fd stencil
  !              4. Compute the jump contribution to the field-interpolation
  !===========================================================================
  !use
  use para
  use field
  use Lagrange
  !==========================================================================
  implicit none
  !============================================================================
  !interface
  interface
     subroutine JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
       integer :: i
       double precision, dimension(1:3), intent(in) :: A,B,C,Q
       double precision, intent(in) :: fA,fB,fC
       double precision, intent(out) :: fQ
       double precision, dimension(1:3) :: vec1,vec2,vec3,vec4
       double precision :: normvec1, normvec2, normvec3, normvec4
     end subroutine JCinterpolate
     subroutine JCinterpolate2(A,B,C,Q,fA,fB,fC,fQ)
       double precision, dimension(1:3), intent(in) :: A,B,C,Q
       double precision, dimension(1:36), intent(in) :: fA,fB,fC
       double precision, dimension(1:36), intent(out) :: fQ
       double precision, dimension(1:3) :: vec1,vec2,vec3,vec4
       double precision :: normvec1, normvec2, normvec3, normvec4
     end subroutine JCinterpolate2
     function coord4intersection(A,B,C,gridline,i) result(coord3)
         double precision:: coord3
         double precision, dimension(1:3), intent(in):: A,B,C
         double precision, dimension(1:2), intent(in):: gridline
         integer, intent(in):: i
         integer:: m,n
         double precision, dimension(1:3):: normal
         double precision:: eps
       end function coord4intersection
       subroutine pointinside(a,b,c,p,s)
         integer :: proceed
         double precision, dimension(3),intent(in) :: a,b,c,p
         double precision, dimension(3) :: vec1,vec2,vec3,vec4,vec5,vec6
         integer, intent(out) :: s
       end subroutine pointinside
       function crossproduct(a,b)
         double precision, dimension(3), intent(in) :: a,b
         double precision, dimension(3) :: crossproduct
       end function crossproduct
    end interface
  !==========================================================================
  integer :: s,ii
  integer :: ic0,ic1,ie0,ie1,jc0,jc1,je0,je1,kc0,kc1,ke0,ke1
  double precision, dimension(1:3) :: A2,B2,C2,P2,Q
  double precision :: coord3
  integer :: kc,ke,je,jc,ie,ic
  double precision, dimension(3), intent(in) :: A,B,C,tau,beta,n,Xcc,Ucc,Omegac
  double precision,dimension(36), intent(in) ::A_JC,B_JC,C_JC
  double precision,dimension(9), intent(in) :: Adup,Adum,Bdup,Bdum,Cdup,Cdum
  double precision, dimension(36) :: Q_JC
  integer, intent(in) :: ip,iobj
  double precision, dimension(3):: u_Q
  !=========================================================================
  !==========================================================================
  ic0 = floor((min(A(1),B(1),C(1))-xc(0))/dx)-1
  ic1 = floor((max(A(1),B(1),C(1))-xc(0))/dx)+1
  ie0 = floor((min(A(1),B(1),C(1))-xe(0))/dx)-1
  ie1 = floor((max(A(1),B(1),C(1))-xe(0))/dx)+1
  
  jc0 = floor((min(A(2),B(2),C(2))-yc(0))/dy)-1
  jc1 = floor((max(A(2),B(2),C(2))-yc(0))/dy)+1
  je0 = floor((min(A(2),B(2),C(2))-ye(0))/dy)-1
  je1 = floor((max(A(2),B(2),C(2))-ye(0))/dy)+1

  kc0 = floor((min(A(3),B(3),C(3))-zc(0))/dz)-1
  kc1 = floor((max(A(3),B(3),C(3))-zc(0))/dz)+1
  ke0 = floor((min(A(3),B(3),C(3))-ze(0))/dz)-1
  ke1 = floor((max(A(3),B(3),C(3))-ze(0))/dz)+1
  !========================================================== 
  ! project the panel into xy plane
  ! find intersection point
  ! interpolate 1st and 2nd derivatives
  !-----------------------------------------------------------
  ! project to xy-plane                                                   
  A2 = (/A(1),A(2),0.0d0/)
  B2 = (/B(1),B(2),0.0d0/)
  C2 = (/C(1),C(2),0.0d0/)
  !write(*,*) nobj
  !-----------------------------------------------------------
  ! icjc
  do j=jc0,jc1
     do i=ic0,ic1
        P2 = (/xc(i),yc(j),0.0d0/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge                   
           P2(1) = P2(1)+eps ! shift the grid line x-direction                 
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(2) = P2(2)+eps ! shift y-direction                            
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(2)/),3)
           Q = (/xc(i),yc(j),coord3/)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           !----------------------------------
           kc = floor((coord3-zc(0))/dz) ! interect at kc:kc+1
           ke = floor((coord3-ze(0))/dz)  ! intersection @ kf:kf+1!!
           
           nicjc = nicjc+1
           
           icjc(1,nicjc,1) = i
           icjc(2,nicjc,1) = j
           icjc(3,nicjc,1) = ke
           icjc(4,nicjc,1) = kc
           icjc(5,nicjc,1) = ip
           icjc(6,nicjc,1) = iobj

           ficjc(1:3,nicjc,1) = Q
           ficjc(4:6,nicjc,1) = tau
           ficjc(7:9,nicjc,1) = beta
           ficjc(10:12,nicjc,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           ficjc(13:48,nicjc,1) = Q_JC
           ficjc(50:55,nicjc,1) =  0.0d0!ficjc(43:48,nicjc,1)

           ficjc(56,nicjc,1) = u_Q(1)
           ficjc(57,nicjc,1) = u_Q(2)
           ficjc(58,nicjc,1) = u_Q(3)

           !u
           call JCinterpolate(A,B,C,Q,Adup(1),Bdup(1),Cdup(1),ficjc(59,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(1),Bdum(1),Cdum(1),ficjc(60,nicjc,1))
           
           call JCinterpolate(A,B,C,Q,Adup(2),Bdup(2),Cdup(2),ficjc(61,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(2),Bdum(2),Cdum(2),ficjc(62,nicjc,1))

           call JCinterpolate(A,B,C,Q,Adup(3),Bdup(3),Cdup(3),ficjc(63,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(3),Bdum(3),Cdum(3),ficjc(64,nicjc,1))
           !v
           call JCinterpolate(A,B,C,Q,Adup(4),Bdup(4),Cdup(4),ficjc(65,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(4),Bdum(4),Cdum(4),ficjc(66,nicjc,1))

           call JCinterpolate(A,B,C,Q,Adup(5),Bdup(5),Cdup(5),ficjc(67,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(5),Bdum(5),Cdum(5),ficjc(68,nicjc,1))

           call JCinterpolate(A,B,C,Q,Adup(6),Bdup(6),Cdup(6),ficjc(69,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(6),Bdum(6),Cdum(6),ficjc(70,nicjc,1))
           ! w
           call JCinterpolate(A,B,C,Q,Adup(7),Bdup(7),Cdup(7),ficjc(71,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(7),Bdum(7),Cdum(7),ficjc(72,nicjc,1))

           call JCinterpolate(A,B,C,Q,Adup(8),Bdup(8),Cdup(8),ficjc(73,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(8),Bdum(8),Cdum(8),ficjc(74,nicjc,1))

           call JCinterpolate(A,B,C,Q,Adup(9),Bdup(9),Cdup(9),ficjc(75,nicjc,1))
           call JCinterpolate(A,B,C,Q,Adum(9),Bdum(9),Cdum(9),ficjc(76,nicjc,1))
           
        endif
     enddo
  enddo
  !----------------------------------------------------------
  ! icje
  do j=je0,je1
     do i =ic0,ic1
        P2 = (/xc(i),ye(j),0.0d0/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge          
           P2(1) = P2(1)+eps ! shift the grid line x-direction                     
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(2) = P2(2)+eps ! shift y-direction        
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(2)/),3)
           Q = (/xc(i),ye(j),coord3/)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           kc = floor((coord3-zc(0))/dz)  ! intersect @ kc:kc+1
           ke = floor((coord3-ze(0))/dz)  ! intersection @ kf:kf+1
           !---------------
           nicje = nicje +1

           icje(1,nicje,1) = i
           icje(2,nicje,1) = j
           icje(3,nicje,1) = ke
           icje(4,nicje,1) = kc
           icje(5,nicje,1) = ip
           icje(6,nicje,1) = iobj


           ficje(1:3,nicje,1) = Q
           ficje(4:6,nicje,1) = tau
           ficje(7:9,nicje,1) = beta
           ficje(10:12,nicje,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           ficje(13:48,nicje,1) = Q_JC
           ficje(50:55,nicje,1) =  0.0d0!ficje(43:48,nicje,1)

           
           ficje(56,nicje,1) = u_Q(1)
           ficje(57,nicje,1) = u_Q(2)
           ficje(58,nicje,1) = u_Q(3)

           !v                                                                                                                                              
           call JCinterpolate(A,B,C,Q,Adup(6),Bdup(6),Cdup(6),ficje(59,nicje,1))
           call JCinterpolate(A,B,C,Q,Adum(6),Bdum(6),Cdum(6),ficje(60,nicje,1))
           ! w                                                                              
           call JCinterpolate(A,B,C,Q,Adup(9),Bdup(9),Cdup(9),ficje(61,nicje,1))
           call JCinterpolate(A,B,C,Q,Adum(9),Bdum(9),Cdum(9),ficje(62,nicje,1))

           

        endif
     enddo !end i                 
  enddo !end j    
  !----------------------------------------------------------
  ! iejc                                                   
  do j=jc0,jc1
     do i =ie0,ie1
        P2 = (/xe(i),yc(j),0.0d0/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge               
           P2(1) = P2(1)+eps ! shift the grid line x-direction    
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(2) = P2(2)+eps ! shift y-direction
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(2)/),3)
           Q = (/xe(i),yc(j),coord3/)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           
           kc = floor((coord3-zc(0))/dz) ! interect at kc:kc+1         
           ke = floor((coord3-ze(0))/dz)  ! intersection @ kf:kf+1

           niejc = niejc+1

           iejc(1,niejc,1) = i
           iejc(2,niejc,1) = j
           iejc(3,niejc,1) = ke
           iejc(4,niejc,1) = kc
           iejc(5,niejc,1) = ip
           iejc(6,niejc,1) = iobj
           
           fiejc(1:3,niejc,1) = Q
           fiejc(4:6,niejc,1) = tau
           fiejc(7:9,niejc,1) = beta
           fiejc(10:12,niejc,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           fiejc(13:48,niejc,1) = Q_JC
           fiejc(50:55,niejc,1) = 0.0d0! fiejc(43:48,niejc,1)


           fiejc(56,niejc,1) = u_Q(1)
           fiejc(57,niejc,1) = u_Q(2)
           fiejc(58,niejc,1) = u_Q(3)

           !u                                                                                                                                              
           call JCinterpolate(A,B,C,Q,Adup(3),Bdup(3),Cdup(3),fiejc(59,niejc,1))
           call JCinterpolate(A,B,C,Q,Adum(3),Bdum(3),Cdum(3),fiejc(60,niejc,1))
           !v                                                                                                                                             
           call JCinterpolate(A,B,C,Q,Adup(9),Bdup(9),Cdup(9),fiejc(61,niejc,1))
           call JCinterpolate(A,B,C,Q,Adum(9),Bdum(9),Cdum(9),fiejc(62,niejc,1))
           
        endif
     enddo !end i
  enddo !end j      
  !==========================================================
  ! project the panel into xz plane                                     
  ! find intersection point                                               
  ! interpolate                  
  A2 = (/A(1),0.0d0,A(3)/)
  B2 = (/B(1),0.0d0,B(3)/)
  C2 = (/C(1),0.0d0,C(3)/)
  !---------------------------------------------------------
  ! ickc
  do k=kc0,kc1
     do i=ic0,ic1
        P2 = (/xc(i),0.0d0,zc(k)/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge 
           P2(1) = P2(1)+eps ! shift the grid line x-direction                   
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction   
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
           Q = (/ xc(i) ,coord3 ,zc(k)  /)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           
           je = floor((coord3-ye(0))/dy) ! intersect @ jf:jf+1
           jc = floor((coord3-yc(0))/dy)!

           !-----
           nickc =nickc+1

           ickc(1,nickc,1) = i
           ickc(2,nickc,1) = k
           ickc(3,nickc,1) = je
           ickc(4,nickc,1) = jc
           ickc(5,nickc,1) = ip
           ickc(6,nickc,1) =iobj
           
           fickc(1,nickc,1) = Q(1)!Qx
           fickc(2,nickc,1) = Q(3)!Qz
           fickc(3,nickc,1) = Q(2)!Qy
           fickc(4:6,nickc,1) = tau
           fickc(7:9,nickc,1) = beta
           fickc(10:12,nickc,1) = n
           
           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)

           fickc(13:48,nickc,1) = Q_JC
           fickc(50:55,nickc,1) =  0.0d0!fickc(43:48,nickc,1)



           fickc(56,nickc,1) = u_Q(1)
           fickc(57,nickc,1) = u_Q(2)
           fickc(58,nickc,1) = u_Q(3)

           !u
           call JCinterpolate(A,B,C,Q,Adup(1),Bdup(1),Cdup(1),fickc(59,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(1),Bdum(1),Cdum(1),fickc(60,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(2),Bdup(2),Cdup(2),fickc(61,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(2),Bdum(2),Cdum(2),fickc(62,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(3),Bdup(3),Cdup(3),fickc(63,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(3),Bdum(3),Cdum(3),fickc(64,nickc,1))
           !v                                                                                                                                              
           call JCinterpolate(A,B,C,Q,Adup(4),Bdup(4),Cdup(4),fickc(65,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(4),Bdum(4),Cdum(4),fickc(66,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(5),Bdup(5),Cdup(5),fickc(67,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(5),Bdum(5),Cdum(5),fickc(68,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(6),Bdup(6),Cdup(6),fickc(69,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(6),Bdum(6),Cdum(6),fickc(70,nickc,1))
           ! w                                                                                                                                             
           call JCinterpolate(A,B,C,Q,Adup(7),Bdup(7),Cdup(7),fickc(71,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(7),Bdum(7),Cdum(7),fickc(72,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(8),Bdup(8),Cdup(8),fickc(73,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(8),Bdum(8),Cdum(8),fickc(74,nickc,1))
           call JCinterpolate(A,B,C,Q,Adup(9),Bdup(9),Cdup(9),fickc(75,nickc,1))
           call JCinterpolate(A,B,C,Q,Adum(9),Bdum(9),Cdum(9),fickc(76,nickc,1))
        endif
     enddo
  enddo
  !-------------------------------------------------------------------
  ! icke
  do k = ke0,ke1
     do i=ic0,ic1
        P2 = (/xc(i),0.0d0,ze(k)/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge                     
           P2(1) = P2(1)+eps ! shift the grid line x-direction           
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction                            
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
           Q = (/ xc(i) ,coord3 ,ze(k)  /)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           je = floor((coord3-ye(0))/dy) ! intersect @ jf:jf+1  !
           jc = floor((coord3-yc(0))/dy)  ! intersect at jc:jc+1

           nicke = nicke+1

           icke(1,nicke,1) = i
           icke(2,nicke,1) = k
           icke(3,nicke,1) = je
           icke(4,nicke,1) = jc
           icke(5,nicke,1) = ip
           icke(6,nicke,1) = iobj!!
           
           ficke(1,nicke,1) = Q(1)!Qx
           ficke(2,nicke,1) = Q(3)!Qz
           ficke(3,nicke,1) = Q(2)!Qy   

           ficke(4:6,nicke,1) = tau
           ficke(7:9,nicke,1) = beta
           ficke(10:12,nicke,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           ficke(13:48,nicke,1) = Q_JC
           ficke(50:55,nicke,1) =  0.0d0! ficke(43:48,nicke,1)

           ficke(56,nicke,1) = u_Q(1)
           ficke(57,nicke,1) = u_Q(2)
           ficke(58,nicke,1) = u_Q(3)
           !v
           call JCinterpolate(A,B,C,Q,Adup(5),Bdup(5),Cdup(5),ficke(59,nicke,1))
           call JCinterpolate(A,B,C,Q,Adum(5),Bdum(5),Cdum(5),ficke(60,nicke,1))
           ! w
           call JCinterpolate(A,B,C,Q,Adup(8),Bdup(8),Cdup(8),ficke(61,nicke,1))
           call JCinterpolate(A,B,C,Q,Adum(8),Bdum(8),Cdum(8),ficke(62,nicke,1))
        endif
     enddo
  enddo
  !-----------------------------------
  ! iekc
  do k=kc0,kc1
     do i=ie0,ie1
        P2 = (/xe(i),0.0d0,zc(k)/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge
           P2(1) = P2(1)+eps ! shift the grid line x-direction
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
           Q = (/ xe(i) ,coord3 ,zc(k)  /)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           
           je = floor((coord3-ye(0))/dy) ! intersect @ je:je+1 
           jc = floor((coord3-yc(0))/dy)
          
           niekc =niekc+1
           
           iekc(1,niekc,1) = i
           iekc(2,niekc,1) = k
           iekc(3,niekc,1) = je
           iekc(4,niekc,1) = jc
           iekc(5,niekc,1) = ip
           iekc(6,niekc,1) = iobj

           !fiekc(1:3,niekc,1) = Q
           fiekc(1,niekc,1) = Q(1)!Qx
           fiekc(2,niekc,1) = Q(3)!Qz
           fiekc(3,niekc,1) = Q(2)!Qy   

           fiekc(4:6,niekc,1) = tau
           fiekc(7:9,niekc,1) = beta
           fiekc(10:12,niekc,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           fiekc(13:48,niekc,1) = Q_JC
           fiekc(50:55,niekc,1) =  0.0d0!fiekc(43:48,niekc,1)

           fiekc(56,niekc,1) = u_Q(1)
           fiekc(57,niekc,1) = u_Q(2)
           fiekc(58,niekc,1) = u_Q(3)

           
           !u
           
           call JCinterpolate(A,B,C,Q,Adup(2),Bdup(2),Cdup(2),fiekc(59,niekc,1))
           call JCinterpolate(A,B,C,Q,Adum(2),Bdum(2),Cdum(2),fiekc(60,niekc,1))
           !v
           
           call JCinterpolate(A,B,C,Q,Adup(5),Bdup(5),Cdup(5),fiekc(61,niekc,1))
           call JCinterpolate(A,B,C,Q,Adum(5),Bdum(5),Cdum(5),fiekc(62,niekc,1))


        endif
     enddo
  enddo
  !========================================================== 
  ! project the panel into yz  plane    
  ! find intersection point         
  ! interpolate 
  !========================================================== 
  A2 = (/0.0d0,A(2),A(3)/)
  B2 = (/0.0d0,B(2),B(3)/)
  C2 = (/0.0d0,C(2),C(3)/)
  !-------------------------
  ! jckc  
  do k=kc0,kc1
     do j=jc0,jc1
        P2 = (/0.0d0,yc(j),zc(k)/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge 
           P2(2) = P2(2)+eps ! shift the grid line y-direction  
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction   
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
           Q  = (/ coord3,yc(j),zc(k) /)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           
           ie = floor((coord3-xe(0))/dx)  ! intersect @ iff:iff+1
           ic = floor((coord3-xc(0))/dx) 
           !--------------------------------
           !------------------------------
           njckc =njckc+1
           
           jckc(1,njckc,1) = j
           jckc(2,njckc,1) = k
           jckc(3,njckc,1) = ie
           jckc(4,njckc,1) = ic
           jckc(5,njckc,1) = ip
           jckc(6,njckc,1) = iobj!

           !fjckc(1:3,njckc,1) = Q
           fjckc(1,njckc,1) = Q(2)!Qy
           fjckc(2,njckc,1) = Q(3)!Qz
           fjckc(3,njckc,1) = Q(1)!Qx   

           fjckc(4:6,njckc,1) = tau
           fjckc(7:9,njckc,1) = beta
           fjckc(10:12,njckc,1) = n!!

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           fjckc(13:48,njckc,1) = Q_JC
           fjckc(50:55,njckc,1) =  0.0d0!fjckc(43:48,njckc,1)

           fjckc(56,njckc,1) = u_Q(1)
           fjckc(57,njckc,1) = u_Q(2)
           fjckc(58,njckc,1) = u_Q(3)

           !u
           call JCinterpolate(A,B,C,Q,Adup(1),Bdup(1),Cdup(1),fjckc(59,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(1),Bdum(1),Cdum(1),fjckc(60,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(2),Bdup(2),Cdup(2),fjckc(61,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(2),Bdum(2),Cdum(2),fjckc(62,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(3),Bdup(3),Cdup(3),fjckc(63,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(3),Bdum(3),Cdum(3),fjckc(64,njckc,1))
           !v
           call JCinterpolate(A,B,C,Q,Adup(4),Bdup(4),Cdup(4),fjckc(65,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(4),Bdum(4),Cdum(4),fjckc(66,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(5),Bdup(5),Cdup(5),fjckc(67,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(5),Bdum(5),Cdum(5),fjckc(68,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(6),Bdup(6),Cdup(6),fjckc(69,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(6),Bdum(6),Cdum(6),fjckc(70,njckc,1))
           ! w
           call JCinterpolate(A,B,C,Q,Adup(7),Bdup(7),Cdup(7),fjckc(71,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(7),Bdum(7),Cdum(7),fjckc(72,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(8),Bdup(8),Cdup(8),fjckc(73,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(8),Bdum(8),Cdum(8),fjckc(74,njckc,1))
           call JCinterpolate(A,B,C,Q,Adup(9),Bdup(9),Cdup(9),fjckc(75,njckc,1))
           call JCinterpolate(A,B,C,Q,Adum(9),Bdum(9),Cdum(9),fjckc(76,njckc,1))
        endif
     enddo
  enddo
  !-------------------------------------------------
  ! jcke
  do k=ke0,ke1
     do j=jc0,jc1
        P2 = (/0.0d0,yc(j),ze(k)/)
        call pointinside(A2,B2,C2,P2,s)
        if( s == 0) then ! the grid line through panel edge 
           P2(2) = P2(2)+eps ! shift the grid line y-direction  
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction 
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
           Q = (/coord3,yc(j),ze(k)/)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)

           
           ic = floor((coord3-xc(0))/dx)  ! intersect @ ic:ic+1
           ie= floor((coord3-xe(0))/dx)  ! intersect @ iff:iff+1 
           
           njcke =njcke+1
           !-----------
           jcke(1,njcke,1) = j
           jcke(2,njcke,1) = k
           jcke(3,njcke,1) = ie
           jcke(4,njcke,1) = ic
           jcke(5,njcke,1) = ip
           jcke(6,njcke,1) = iobj!

           !fjcke(1:3,njcke,1) = Q
           fjcke(1,njcke,1) = Q(2)
           fjcke(2,njcke,1) = Q(3)
           fjcke(3,njcke,1) = Q(1)
           fjcke(4:6,njcke,1) = tau
           fjcke(7:9,njcke,1) = beta
           fjcke(10:12,njcke,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           fjcke(13:48,njcke,1) = Q_JC
           fjcke(50:55,njcke,1) =  0.0d0!fjcke(43:48,njcke,1)
           !------------
           
           fjcke(56,njcke,1) = u_Q(1)
           fjcke(57,njcke,1) = u_Q(2)
           fjcke(58,njcke,1) = u_Q(3)
           !u
           call JCinterpolate(A,B,C,Q,Adup(1),Bdup(1),Cdup(1),fjcke(59,njcke,1))
           call JCinterpolate(A,B,C,Q,Adum(1),Bdum(1),Cdum(1),fjcke(60,njcke,1))
           ! w                                                                                                                                             
           call JCinterpolate(A,B,C,Q,Adup(7),Bdup(7),Cdup(7),fjcke(61,njcke,1))
           call JCinterpolate(A,B,C,Q,Adum(7),Bdum(7),Cdum(7),fjcke(62,njcke,1))           
        endif
     enddo
  enddo
  !-------------------------------------------------
  ! jekc     
  do k=kc0,kc1
     do j=je0,je1
        P2 = (/0.0d0,ye(j),zc(k)/)
        call pointinside(A2,B2,C2,P2,s)!
        if( s == 0) then ! the grid line through panel edge   
           P2(2) = P2(2)+eps ! shift the grid line y-direction  
           call pointinside(A2,B2,C2,P2,s)
           if(s==0) then
              P2(3) = P2(3)+eps ! shift z-direction   
              call pointinside(A2,B2,C2,P2,s)
           endif
        endif
        if(s==1) then
           coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
           Q      = (/ coord3,ye(j),zc(k)  /)
           u_Q = Ucc + crossproduct(Omegac,Q-Xcc)
           
           ic     = floor((coord3-xc(0))/dx)    ! intersect @ ic:ic+1
           ie     = floor((coord3-xe(0))/dx)
           !------------------------------------------------------------
           !---------------------------------------------------------
           njekc =njekc+1
             
           jekc(1,njekc,1) = j
           jekc(2,njekc,1) = k
           jekc(3,njekc,1) = ie
           jekc(4,njekc,1) = ic
           jekc(5,njekc,1) = ip
           jekc(6,njekc,1) = iobj


           !fjekc(1:3,njekc,1) = Q
           fjekc(1,njekc,1) = Q(2)!Qy
           fjekc(2,njekc,1) = Q(3)!Qz
           fjekc(3,njekc,1) = Q(1)!Qx

           fjekc(4:6,njekc,1) = tau
           fjekc(7:9,njekc,1) = beta
           fjekc(10:12,njekc,1) = n

           call JCinterpolate2(A,B,C,Q,A_JC,B_JC,C_JC,Q_JC)
           fjekc(13:48,njekc,1) = Q_JC
           fjekc(50:55,njekc,1) = 0.0d0!fjekc(43:48,njekc,1)
           !---------------
           fjekc(56,njekc,1) = u_Q(1)
           fjekc(57,njekc,1) = u_Q(2)
           fjekc(58,njekc,1) = u_Q(3)

           !u                                                                                                
           call JCinterpolate(A,B,C,Q,Adup(1),Bdup(1),Cdup(1),fjekc(59,njekc,1))
           call JCinterpolate(A,B,C,Q,Adum(1),Bdum(1),Cdum(1),fjekc(60,njekc,1))
           !v                                                                                                  
           call JCinterpolate(A,B,C,Q,Adup(4),Bdup(4),Cdup(4),fjekc(61,njekc,1))
           call JCinterpolate(A,B,C,Q,Adum(4),Bdum(4),Cdum(4),fjekc(62,njekc,1))
           !
        endif
     enddo
  enddo

  return 
end subroutine IrrP_JC2
!----------------------------------------------------------------------------------------
