!----------------------------------------------------------------------
subroutine singular_call
  !-----------------------------------
  use para
  use field
  use Lagrange
  !-----------------------------------
  integer*4 :: time1(3),time2(3)
  !-----------------------------------
  call JC4uvwp
  !---------------------
  call jc_velocity
  !-------------------
  call correction_interpolate
  call correction_strain
  
  call u_interpolate
  call u_strain
  
  call v_interpolate
  call v_strain
  
  call w_interpolate
  call w_strain
  !------------------------------------
  !-----------------------------------
  
  call jc_pressure 
  call correction_difference

2000 format (A,' time: ',i2.2, ':', i2.2, ':', i2.2 )
  return
end subroutine singular_call


!-----------------------------------------------------------------------
