!===================================================================
! September 2020
!===================================================================
subroutine raycrossing
  !-----------------------------------------------------------------
  !  This subroutine is to :
  !          1. find the intersection points between the grid lines
  !             and the interface
  !          2. determin which points insid an object and which one
  !             is outside
  !--------------------------------------------------------------------
  use para
  use field
  use Lagrange
  !--------------------
  implicit none

  integer :: s
  integer :: myobj, mypanel
  integer :: nA,nB,nC
  
  integer :: ic0,ic1,ie0,ie1,jc0,jc1,je0,je1,kc0,kc1,ke0,ke1
  integer :: kc,ke

  double precision, dimension(3) :: A,B,C
  double precision, dimension(3) :: A2,B2,C2,P2
  double precision :: coord3
  !-------------------
  interface
     subroutine pointinside(a,b,c,p,s)
       integer :: proceed
       double precision, dimension(3),intent(in) :: a,b,c,p
       double precision, dimension(3) :: vec1,vec2,vec3,vec4,vec5,vec6
       integer, intent(out) :: s
     endsubroutine pointinside
     function coord4intersection(A,B,C,gridline,i) result(coord3)
       double precision:: coord3
       double precision, dimension(1:3), intent(in):: A,B,C
       double precision, dimension(1:2), intent(in):: gridline
       integer, intent(in):: i
     end function coord4intersection
  end interface
  !--------------------------------------------
  ! initialize the io pointer to be zero
  iou = 0
  iov = 0
  iow = 0
  iop = 0
  ! Initialize the numbers of intersection points to be zero 
  nicjc = 0
  nicje = 0
  niejc = 0

  nickc = 0
  nicke = 0
  niekc = 0

  njckc = 0
  njcke = 0
  njekc = 0
  !--------------------------------------------
  !--------------------------------------------
  do obj =1,nobj  ! loop through objects
     !---------------------------------
     myobj = objects(obj) 
     !-------------------------------------------------------
     do mypanel = np4obj(obj-1)+1,np4obj(obj) !loop through panels
        !------------------------------------------------------
        ! the vertices of the current panel
        nA = panel(1,mypanel)+nv4obj(obj-1)
        nB = panel(2,mypanel)+nv4obj(obj-1)
        nC = panel(3,mypanel)+nv4obj(obj-1)

        A = vertex(:,nA)
        B = vertex(:,nB)
        C = vertex(:,nC)
        !-------------------------------------------------------
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
        !====================================================================
        ! project to xy-plane
        A2 = (/A(1),A(2),0.0d0/)
        B2 = (/B(1),B(2),0.0d0/)
        C2 = (/C(1),C(2),0.0d0/)
        !-------------------------------------------------------------------
        ! icjc
        do j = jc0,jc1
           do i = ic0,ic1
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
                 !find the z coordinate of the intersection point
                 coord3 = coord4intersection(A,B,C,(/P2(1),P2(2)/),3) 
                 kc = floor((coord3-zc(0))/dz) ! interect at kc:kc+1
                 !update the 3d array IO
                 do k = 0,kc
                    iop(i,j,k) = myobj-iop(i,j,k) ! io array for p @ xc,yc,zc
                 enddo
                 ke = floor((coord3-ze(0))/dz)  ! intersection @ ke:ke+1
                 do k = 0,ke
                    iow(i,j,k) = myobj-iow(i,j,k) ! io array for w @ xc,yc,ze
                 enddo
                 nicjc = nicjc+1
              endif
           enddo !endi
        enddo !endj
        !------------------------------------------------------------------------------------------
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
                 kc = floor((coord3-zc(0))/dz)  ! intersect @ kc:kc+1
                 do k = 0,kc
                    iov(i,j,k) = myobj-iov(i,j,k)  ! the io array for the velocity component v
                 enddo
                 nicje = nicje +1
              endif
           enddo !end i
        enddo !end j
        !----------------------------------------------------------------
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
                 kc = floor((coord3-zc(0))/dz)    ! intersect @ kc:kc+1
                 do k = 0,kc
                    iou(i,j,k) = myobj-iou(i,j,k)   ! the io array for the velocity component u
                 enddo
                 niejc = niejc+1
              endif
           enddo !end i
        enddo !end j
        !====================================================================================
        ! project into xz-plane
        A2 = (/A(1),0.0d0,A(3)/)
        B2 = (/B(1),0.0d0,B(3)/)
        C2 = (/C(1),0.0d0,C(3)/)
        !-------------------------------------------------------------
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
                 !coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
                 nickc =nickc+1
              endif
           enddo
        enddo
        !-----------------------------------------------------------------
        ! icke     
         do k=ke0,ke1
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
                 !coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
                 nicke = nicke+1
              endif
           enddo
        enddo
        !------------------------------------------------------------------
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
                 !coord3 = coord4intersection(A,B,C,(/P2(1),P2(3)/),2)
                 niekc = niekc+1
              endif
           enddo
        enddo
        !===================================================================
        !project into yz-plane
        A2 = (/0.0d0,A(2),A(3)/)
        B2 = (/0.0d0,B(2),B(3)/)
        C2 = (/0.0d0,C(2),C(3)/)
        !------------------------------------------------------------------
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
                 !coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
                 njckc =njckc+1
              endif
           enddo
        enddo
        !-------------------------------------------------------------------
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
                 !coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
                 njcke =njcke+1
              endif
           enddo
        enddo
        !------------------------------------------------------------------
        ! jekc
        do k=kc0,kc1
           do j=je0,je1
              P2 = (/0.0d0,ye(j),zc(k)/)
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
                 !coord3 = coord4intersection(A,B,C,(/P2(2),P2(3)/),1)
                 njekc =njekc+1
              endif
           enddo
        enddo      
     enddo  ! end panel
  enddo ! end objects

  return
  
end subroutine raycrossing
!--------------------------------------------------------------------------------
