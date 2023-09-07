subroutine cjc4d2u(tau,beta,rhs,jumpd2u)

  !
  implicit none
  !----
  interface
     subroutine c2_solver(a,b,c,x,y,z,r,icheck,v)
       integer, intent(in) :: icheck
       double precision :: a,b,c,x,y,z,r(7),v(6)
       double precision :: s2,s3,m1,m2,m3,e,f,g,h
     end subroutine c2_solver
  end interface
  !======================
  double precision, dimension(7), intent(in) :: rhs
  double precision, dimension(3), intent(in) :: tau,beta
  double precision, dimension(6), intent(out) :: jumpd2u
  double precision, dimension(7) :: dd
  integer :: iv
  integer, dimension(7) :: id
  integer, dimension(6) :: ix
  double precision, dimension(6):: xx
  double precision :: fo1,fo2,fo3,bo1,bo2,bo3
  !============================

  ! if tau(1).ne.0, then we have equation (96) on 2085      
   if(abs(tau(1)).ge.abs(tau(2)).and.abs(tau(1)).ge.abs(tau(3))) then

      id(1) = 7
      id(2) = 1
      id(3) = 4
      id(4) = 2
      id(5) = 5
      id(6) = 3
      id(7) = 6

      ix(1) = 1
      ix(2) = 2
      ix(3) = 3
      ix(4) = 4
      ix(5) = 5
      ix(6) = 6

      fo1 = tau(1)
      fo2 = tau(2)
      fo3 = tau(3)

      bo1 = beta(1)
      bo2 = beta(2)
      bo3 = beta(3)

   endif

   ! if tau(2).ne.0, then we have equation (97) on 2085
   if(abs(tau(2)).ge.abs(tau(1)).and.abs(tau(2)).ge.abs(tau(3))) then

      id(1) = 7
      id(2) = 2
      id(3) = 5
      id(4) = 3
      id(5) = 6
      id(6) = 1
      id(7) = 4

      ix(1) = 6
      ix(2) = 3
      ix(3) = 5
      ix(4) = 1
      ix(5) = 2
      ix(6) = 4

      fo1 = tau(2)
      fo2 = tau(3)
      fo3 = tau(1)

      bo1 = beta(2)
      bo2 = beta(3)
      bo3 = beta(1)
   endif
   
   ! if tau(3).ne.0, then we have equation (98) on 2085
   if(abs(tau(3)).ge.abs(tau(1)).and.abs(tau(3)).ge.abs(tau(2))) then

      id(1) = 7
      id(2) = 3
      id(3) = 6
      id(4) = 2
      id(5) = 5
      id(6) = 1
      id(7) = 4

      ix(1) = 6
      ix(2) = 5
      ix(3) = 3
      ix(4) = 4
      ix(5) = 2
      ix(6) = 1

      fo1 = tau(3)
      fo2 = tau(2)
      fo3 = tau(1)

      bo1 = beta(3)
      bo2 = beta(2)
      bo3 = beta(1)
   endif
   
   ! arrange the rhs in the proper order according to id indices            
   do iv = 1,7
      dd(iv) = rhs(id(iv))
   enddo
   
   ! call the c2_solver to find the jump condtions at the current vertex    
   call c2_solver(fo1,fo2,fo3,bo1,bo2,bo3,dd,0,xx)
   ! place the jump conditions in the proper order according to the ix indices
   do iv = 1,6
      jumpd2u(iv) = xx(ix(iv))
   enddo


  
end subroutine cjc4d2u
