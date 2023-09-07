subroutine deallocate_Lagrange

  use Lagrange

  !---------------------------------------
  !locations of intersections                                                                     
  deallocate(icjc)   
  deallocate(icje)    
  deallocate(iejc)

  deallocate(ickc)    
  deallocate(icke) 
  deallocate(iekc) 

  deallocate(jckc) 
  deallocate(jcke) 
  deallocate(jekc)                         
  !--------------------------------------------------------------------------
  deallocate(ficjc) 
  deallocate(ficje)     
  deallocate(fiejc)

  deallocate(fickc)    
  deallocate(ficke)       
  deallocate(fiekc)                                                                           

  deallocate(fjckc) 
  deallocate(fjcke)   
  deallocate(fjekc) 
  !---------------------------------------------------------------------
  deallocate(kuiejc,kvicje,kwicjc)
  deallocate(kwiejc,kwicje,kuicje)

  deallocate(juiekc,jwicke,jvickc)
  deallocate(jviekc,jvicke,jwiekc)

  deallocate(ivjekc,iwjcke,iujckc)
  deallocate(iujekc,iujcke,ivjcke)


  deallocate(ukiejc,vkicje,wkicjc)
  deallocate(wkiejc,wkicje,ukicje)
  deallocate(ujiekc,wjicke,vjickc)
  deallocate(vjiekc,vjicke,wjiekc)
  deallocate(vijekc,wijcke,uijckc)
  deallocate(uijekc,uijcke,vijcke)

  deallocate(udz,vdz,wdz)
  deallocate(udy,vdy,wdy)
  deallocate(udx,vdx,wdx)


  deallocate(uwdz,vwdz,wwdz,pdz)
  deallocate(udzz,vdzz,wdzz,pdzz)

  deallocate(uvdy,wvdy,vvdy,pdy)
  deallocate(udyy,wdyy,vdyy,pdyy)

  deallocate(vudx,wudx,uudx,pdx)
  deallocate(vdxx,wdxx,udxx,pdxx)
  

  return
end subroutine deallocate_Lagrange
