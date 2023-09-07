module Lagrange
  implicit none
  !-------------------------------------------------
  !objects' info.
  integer :: nobj,obj,nvertices,npanels
  !integer :: iv,m,ms
  integer, dimension(:), allocatable :: objects,nv4obj,np4obj
  integer, dimension(:,:), allocatable :: panel
  double precision, dimension(:,:), allocatable :: vertex,curvature
  double precision, dimension(:), allocatable:: xsc0,ysc0,zsc0,phi0,theta0,psi0
  double precision, dimension(:), allocatable:: xsc,ysc,zsc,phi,theta,psi
  double precision, dimension(:), allocatable:: xsct,ysct,zsct,phit,thetat,psit
  double precision, dimension(:), allocatable:: xsctt,ysctt,zsctt,phitt,thetatt,psitt
  double precision, dimension(:), allocatable:: omegax,omegay,omegaz,omegaxt,omegayt,omegazt
  !-------------------------------------------------------------
  ! panel-gridline intersections
  !
  ! the number of intersection points
  integer :: nicjc,nicje,niejc
  integer :: nickc,nicke,niekc
  integer :: njckc,njcke,njekc
  ! the location of intersections
  double precision, dimension(:,:,:), allocatable :: icjc,icje,iejc
  double precision, dimension(:,:,:), allocatable :: ickc,icke,iekc
  double precision, dimension(:,:,:), allocatable :: jckc,jcke,jekc
  ! the coordinates and JCs @ intersection
  double precision, dimension(:,:,:), allocatable :: ficjc,ficje,fiejc
  double precision, dimension(:,:,:), allocatable :: fickc,ficke,fiekc
  double precision, dimension(:,:,:), allocatable :: fjckc,fjcke,fjekc
  !---------------------------------------------------------------------
  double precision, dimension(:,:),allocatable :: xs,ys,zs,xss,yss,zss ! 1:number of vertices of object, 2: the object 
  double precision, dimension(:,:),allocatable :: fst,snd,us,vs,ws
  !----------------------------------------------------------------------------------
  integer, dimension(:,:),allocatable:: kuiejc,kvicje,kwicjc 
  integer, dimension(:,:),allocatable:: kwiejc,kwicje,kuicje
  integer, dimension(:,:),allocatable:: juiekc,jwicke,jvickc  
  integer, dimension(:,:),allocatable:: jviekc,jvicke,jwiekc
  integer, dimension(:,:),allocatable:: ivjekc,iwjcke,iujckc   
  integer, dimension(:,:),allocatable:: iujekc,iujcke,ivjcke
  double precision, dimension(:,:),allocatable:: ukiejc,vkicje,wkicjc
  double precision, dimension(:,:),allocatable:: wkiejc,wkicje,ukicje
  double precision, dimension(:,:),allocatable:: ujiekc,wjicke,vjickc
  double precision, dimension(:,:),allocatable:: vjiekc,vjicke,wjiekc
  double precision, dimension(:,:),allocatable:: vijekc,wijcke,uijckc
  double precision, dimension(:,:),allocatable:: uijekc,uijcke,vijcke
  !-----------------------------------------------------------------
  double precision, dimension(:,:),allocatable:: udx,vdx,wdx
  double precision, dimension(:,:),allocatable:: udy,vdy,wdy
  double precision, dimension(:,:),allocatable:: udz,vdz,wdz
  !------------------------------------------------------------------
  double precision, dimension(:,:),allocatable:: uudx,uvdy,uwdz,pdx
  double precision, dimension(:,:,:),allocatable:: udxx,udyy,udzz
  double precision, dimension(:,:),allocatable:: vudx,vvdy,vwdz,pdy
  double precision, dimension(:,:,:),allocatable:: vdxx,vdyy,vdzz
  double precision, dimension(:,:),allocatable:: wudx,wvdy,wwdz,pdz
  double precision, dimension(:,:,:),allocatable:: wdxx,wdyy,wdzz
  double precision, dimension(:,:,:),allocatable:: pdxx,pdyy,pdzz
  !----------------------------------------------------------------------------------
  double precision, dimension(:), allocatable::pJC0 ! initial guess for the JC of P
  !----------------------------------------------------------------------------------
  double precision, dimension(:,:), allocatable::dudnp_T,dvdnp_T,dwdnp_T,dpdnPJC_T
  !----------------------------------------------------------------------------------
end module Lagrange
