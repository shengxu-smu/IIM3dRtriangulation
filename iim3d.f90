program main
  !------------------------------------------------------------------
  ! use:
  use para
  use field
  use Lagrange
  use vfft
  !------------------------------------------------------------------------
  implicit none
  ! declaration:
  integer :: n
  character*8 :: day1,day2
  integer*4 :: time1(3),time2(3)
  double precision :: tmp,tstart,tend
  real*4 etime,total,cost(2)
  !-----------------------------------------------------------------------
  !interface:  
  !=======================================================================
  ! initialization:
  call run_input             !1. input
  call allocate_field        !2. domaian?
  call allocate_vfft
  call mesh                  !3. discretize the domaian(MAC scheme)
  call loadobject            !4. load the objects information
  call allocate_surface
  !--------------------------------------------------
  if(nobj.gt.1) then
     write(*,*)'  !! note: multiple objects !!'
  endif
  !------------------------------------------
  xsc0 = 0.0d0
  ysc0 = 0.0d0
  zsc0 = 0.0d0  
  phi0 = 0.0d0
  theta0 = 0.0d0
  psi0 = 0.0d0
  vertex = 0.50d0*vertex
  !-------
  call raycrossing          !6. find the intersection -----------
  call allocate_Lagrange
  !----------------------------
  call initial
  nstart = 0
  !=====================
  if(iread.eq.1) then
         call data_read
  endif
  !====================================================================
  nstart=nstart+1
  nend=nstart+nstep-1
  tstart=t
  tmp=t
  call date_and_time(day1)
  call itime(time1)
  !--------------------------------------------------------------------
  do n=nstart,nend
     call cfl               ! compute the time step using CFL condition
     write(*,*)
     write(*,*)'! n = ',n,' t = ',t
     print*, dt
     call rk4                ! time integration using RK 
     if(itout.eq.1.and.mod(n,1).eq.0) then
        call time_output
     endif
     if(ianimation.eq.1) then
        if(tmp.ge.t-dt.and.tmp.lt.t)then
           call animation
           tmp=tmp+tout
        endif
     endif
     write(*,*) '========================================'
  enddo
  tend=t
  call date_and_time(day2)
  call itime(time2)
  !-------------------------------------------------------------------
  !------------------------------------------------------------------
  if(iwrite.eq.1) then
     call data_write
  endif
  if(iplot.eq.1) then
     call surface_plot
     call field_plot
  endif

  !------------------------------------------------------------------------
  ! output and error analysis
  total=etime(cost)
  call run_output(tstart,tend,total,day1,time1,day2,time2)
  write(*,*)
  write(*,*)'start @'
  write(*,*)'date: ', day1
  write(*,2000) time1
  write(*,*)'end @'
  write(*,*)'date: ', day2
  write(*,2000) time2
  write(*,*)'user time   = ',cost(1)/3600.0d0,' hours'
  write(*,*)'system time = ',cost(2)/3600.0d0,' hours'
  write(*,*)'total time  = ',total/3600.0d0,' hours'
2000 format (' time: ',i2.2, ':', i2.2, ':', i2.2 )

  stop

end program main
