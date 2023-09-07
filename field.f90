module field
  
  implicit none

  double precision,dimension(:),allocatable  :: x,y,z,xe,ye,ze,xc,yc,zc

  double precision,dimension(:,:),allocatable :: pw,pe,ps,pn,pb,pt
  double precision,dimension(:,:,:),allocatable :: u,v,w,p,d,o
  double precision,dimension(:,:,:),allocatable :: un,vn,wn,dn,um,vm,wm
  double precision,dimension(:,:,:),allocatable :: fcee,fcec,fece,fecc,fcce,fccc,feec
  
  double precision,dimension(:,:,:),allocatable :: data1,data2,data3,data4,data5,data6,data7,data8,data9

  integer*2, dimension(:,:,:), allocatable:: iou,iov,iow,iop
  
end module field
