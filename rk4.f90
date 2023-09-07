subroutine rk4
  !--------------------------------------------------------
  use para
  use field
  use Lagrange
  !--------------------------------------------------------
  integer :: krk,modify
  double precision ::  fac
  integer*4 :: time1(3),time2(3),times1(3),times2(3)
  !-------------------------------------------------------
  modify = 0
  !======================================================
  !----Step 1 of RK4-------------------------------------
  krk=1
  fac=0.5d0 
  if(isingular.eq.1) then
     call resetLagrange
     call singular_call
  else
     call u_interpolate
     call u_strain
     call v_interpolate
     call v_strain
     call w_interpolate
     call w_strain
  endif
  
  call old_save
  call pressure(fac)
  call rhs
  
  !--------
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
           v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
           w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
           !
           um(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
           vm(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
           wm(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
        enddo
     enddo
  enddo
  !----------------------------------------------------
  if(imove.eq.1) then
     call surface_move
     call deallocate_Lagrange
     call raycrossing
     call allocate_Lagrange
  endif 
  !----------------------------------------------------
  call velocity_reset
  call ubc
  call vbc
  call wbc
  !=====================================================
  !=====================================================
  !----Step 2 of RK4------------------------------------- 
  krk=2
  fac=0.5d0
  if(isingular.eq.1) then
     call resetLagrange
     call singular_call
  else
     call u_interpolate
     call u_strain
     call v_interpolate
     call v_strain
     call w_interpolate
     call w_strain
  endif
  
  call pressure(fac)
  call rhs

  
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
           v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
           w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
           um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
           vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
           wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
        enddo
     enddo
  enddo

  !-------------------
  call velocity_reset
  call ubc
  call vbc
  call wbc
  !======================================================================
   !----Step 3 of RK4------------------------------------- 
  krk=3
  fac=1.0d0
  if(isingular.eq.1) then
     call resetLagrange
     call singular_call
  else
     call u_interpolate
     call u_strain
     call v_interpolate
     call v_strain
     call w_interpolate
     call w_strain
  endif

  call pressure(fac)
  call rhs
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
           v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
           w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
           um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
           vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
           wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
        enddo
     enddo
  enddo

  if(imove.eq.1) then
     call surface_move
     call deallocate_Lagrange
     call raycrossing
     call allocate_Lagrange
  endif


  !----------------------------------------------------
  call velocity_reset
  call ubc
  call vbc
  call wbc
  !=====================================================================
  !----Step 4 of RK4------------------------------------- 
  krk=4
  fac=1.0d0  
  if(isingular.eq.1) then
     call resetLagrange
     call singular_call
  else
     call u_interpolate
     call u_strain
     call v_interpolate
     call v_strain
     call w_interpolate
     call w_strain
  endif
  
  call pressure(fac)
  call rhs

  do k=1,nz
     do j=1,ny
        do i=1,nx
           u(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
           v(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
           w(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
        enddo
     enddo
  enddo

  call velocity_reset
  call ubc
  call vbc
  call wbc
  !=======================  
2000 format (A,' time: ',i2.2, ':', i2.2, ':', i2.2 )

  return
  
end subroutine rk4
