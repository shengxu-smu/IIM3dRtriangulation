subroutine IrrP_pressureJC2(nrow,pJC)
  !--------------------------------------------------
  !Description
  !--------------------------------------------------
  !Use
  use para
  use field
  use Lagrange
  !---------------------------------------------------
  !Declaration:
  implicit none
  integer, intent(in) :: nrow
  double precision,dimension(nrow),intent(in):: pJC
  integer :: in,ip,iobj,nA,nB,nC
  integer :: ic,iff,jc,jf,kc,kf
  double precision, dimension(1:3) :: Q,A,B,C
  double precision :: fA,fB,fC,fQ
  !---------------------------------------------------
  interface
     subroutine JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
       integer :: i
       double precision, dimension(1:3), intent(in) :: A,B,C,Q
       double precision, intent(in) :: fA,fB,fC
       double precision, intent(out) :: fQ
       double precision, dimension(1:3) :: vec1,vec2,vec3,vec4
       double precision :: normvec1, normvec2, normvec3, normvec4
     end subroutine JCinterpolate
  end interface
  !----------------------------------------------------
  !====================================================
  !1.icjc
  do in =1,nicjc !loop through intersection
     ! the panel and object 
     ip   = icjc(5,in,1)
     iobj = icjc(6,in,1)
     ! the coord of the intersection
     Q = ficjc(1:3,in,1)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     ficjc(49,in,1) = fQ
  enddo
  !----------------------------------------------------
  !2.icje
  do in =1,nicje !loop through intersection 
     ! location of the intersection 
     ip   = icje(5,in,1)
     iobj = icje(6,in,1)
     ! the coord of the intersection 
     Q = ficje(1:3,in,1)
     ! the indices of the panel's vertices     
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices  
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices  
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     ficje(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !3.iejc
  do in =1,niejc !loop through intersection
     ! location of the intersection 
     ip   = iejc(5,in,1)
     iobj = iejc(6,in,1)
     ! the coord of the intersection
     Q = fiejc(1:3,in,1) 
     ! the indices of the panel's vertices 
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices 
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices  
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fiejc(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !4.ickc
  do in =1,nickc !loop through intersection
     ! location of the intsection pt
     ip   = ickc(5,in,1)
     iobj = ickc(6,in,1)
     ! the coord of the intersection 
     Q = (/ fickc(1,in,1), fickc(3,in,1), fickc(2,in,1) /)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices 
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices 
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices                                                                 
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fickc(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !5.icke
  do in =1,nicke !loop through intersection
     ! location of the intsection pt
     ip   = icke(5,in,1)
     iobj = icke(6,in,1)
     ! the coord of the intersection
     Q = (/ ficke(1,in,1), ficke(3,in,1), ficke(2,in,1) /)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     ficke(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !6.iekc
  do in =1,niekc !loop through intersection 
     ! location of the intsection pt 
     ip   = iekc(5,in,1)
     iobj = iekc(6,in,1)
     ! the coord of the intersection 
     Q = (/ fiekc(1,in,1), fiekc(3,in,1), fiekc(2,in,1) /)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices                                     
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fiekc(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !7.jckc
  do in =1,njckc !loop through intersection
     ! location of the intersection 
     ip   = jckc(5,in,1)
     iobj = jckc(6,in,1)
    ! the coord of the intersection
     Q = (/ fjckc(3,in,1), fjckc(1,in,1), fjckc(2,in,1) /)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices 
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fjckc(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !8.jcke
  do in =1,njcke !loop through intersection
     ! location of the intersection
     ip   = jcke(5,in,1)
     iobj = jcke(6,in,1)
    ! the coord of the intersection
     Q =  (/ fjcke(3,in,1), fjcke(1,in,1), fjcke(2,in,1) /)
     ! the indices of the panel's vertices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fjcke(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  !9.jekc
  do in =1,njekc !loop through intersection
     ! location of the intersection 
     ip   = jekc(5,in,1)
     iobj = jekc(6,in,1)
     ! the coord of the intersection
     Q = (/ fjekc(3,in,1), fjekc(1,in,1), fjekc(2,in,1) /)
     ! the indices of the panel's vetrices
     nA = panel(1,ip)+nv4obj(iobj-1)
     nB = panel(2,ip)+nv4obj(iobj-1)
     nC = panel(3,ip)+nv4obj(iobj-1)
     ! the panel vertices
     A = vertex(:,nA)
     B = vertex(:,nB)
     C = vertex(:,nC)
     ! the Jump conditions at vertices 
     fA = pJC(nA)
     fB = pJC(nB)
     fC = pJC(nC)
     ! interpolate JCs at Q from the vertices 
     call JCinterpolate(A,B,C,Q,fA,fB,fC,fQ)
     fjekc(49,in,1) = fQ
  enddo
  !---------------------------------------------------- 
  
end subroutine IrrP_pressureJC2
