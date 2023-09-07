!-------------------------------------------------------------------------------
!  December 2020
!-------------------------------------------------------------------------------
subroutine loadobject
  !-----------------------------------------------------------------------------
  ! Description:
  !  This subroutine is to load the objects information for the Input txt files
  !  DIRECTION: /INPUT/OBJECT
  !  Files:
  !       1- objlist:
  !           nobj  - the number of object
  !           the objects list
  !       2. vobj1,vobj2,....
  !           nvertex
  !           ndimension
  !           the coordinates of vertices 
  !       3. pobj1,pobj2,....
  !           npanel
  !           npoint
  !           the vertices of each panel
  !       4. cobj1,cobj2, .....
  !           curvature of each vertex 
  !---------------------------------------------------------------------
  use Lagrange
  use para
  !---------------------------------------------------------------------
  ! Declaration:
  !  integer :: i
  integer :: object,nvertex,ndimension,npanel,npoint
  integer :: myobj,mypanel,myvertex,nk
  double precision :: h
  
  character(5) filesuffix

  namelist /objectinfo/nobj
  namelist /vertexinfo/nvertex,ndimension
  namelist /panelinfo/npanel,npoint,h
  namelist /curvatureinfo/nk
  !---------------------------------------------------------------------
  ! load the number of object and objects list from the tx file objlist
  open(100,file='./DAT/OBJECT/objlist')
  read(100,objectinfo)
  allocate(objects(nobj))
  read(100,*) (objects(i),i=1,nobj)
  close(100)
  !---------------------------------------------------------------------
  allocate(nv4obj(0:nobj))
  allocate(np4obj(0:nobj))
  !---------------------------------------------------------------------
  nvertices = 0
  npanels = 0
  nv4obj(0) = 0
  np4obj(0) = 0
  !---------------------------------------------------------------------
  do obj = 1,nobj

     myobj = objects(obj)
     write(filesuffix,'(i0)') myobj

     open(101,file='./DAT/OBJECT/vobj'//filesuffix)
     read(101,vertexinfo)
     nvertices = nvertices + nvertex
     nv4obj(obj) = nv4obj(obj-1) + nvertex
     close(101)

     open(102,file='./DAT/OBJECT/pobj'//filesuffix)
     read(102,panelinfo)
     np4obj(obj) = np4obj(obj-1) + npanel
     npanels = npanels + npanel
     close(102)
     
  enddo
  !--------------------------------------------------------
  allocate(vertex(3,nvertices))
  allocate(panel(3,npanels))
  allocate(curvature(3,npanels))

  do obj =1,nobj
     
     myobj = objects(obj)
     write(filesuffix,'(i0)') myobj

     open(101,file='./DAT/OBJECT/vobj'//filesuffix)
     read(101,vertexinfo)
     read(101,*)((vertex(i,myvertex),i=1,3),myvertex = nv4obj(obj-1)+1,nv4obj(obj))
     close(101)

     open(102,file='./DAT/OBJECT/pobj'//filesuffix)
     read(102,panelinfo)
     read(102,*)((panel(i,mypanel),i=1,3),mypanel=np4obj(obj-1)+1,np4obj(obj))
     close(102)
     
     open(103,file='./DAT/OBJECT/cobj'//filesuffix)
     read(103,curvatureinfo)
     read(103,*)((curvature(i,mypanel),i=1,3),mypanel=np4obj(obj-1)+1,np4obj(obj))
     close(103) 

  enddo
  !-------------------------------------------------------------------------------  
end subroutine loadobject
