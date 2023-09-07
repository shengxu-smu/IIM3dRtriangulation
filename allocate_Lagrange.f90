subroutine allocate_Lagrange

  use Lagrange

  !---------------------------------------
  !locations of intersections                                                                     
  allocate(icjc(6,nicjc,1)) ! 1: ic, 2: jc, 3: ke, 4: kc, 5: panel, 6: object                        
  allocate(icje(6,nicje,1)) ! 1: ic, 2: je, 3: ke, 4: kc, 5: panel, 6: object                        
  allocate(iejc(6,niejc,1)) ! 1: ie, 2: jc, 3: ke, 4: kc, 5: panel, 6: object                        

  allocate(ickc(6,nickc,1)) ! 1: ic, 2: kc, 3: je, 4: jc, 5: panel, 6: object                        
  allocate(icke(6,nicke,1)) ! 1: ic, 2: ke, 3: je, 4: jc, 5: panel, 6: object                        
  allocate(iekc(6,niekc,1)) ! 1: ie, 2: kc, 3: je, 4: jc, 5: panel, 6: object                        

  allocate(jckc(6,njckc,1)) ! 1: jc, 2: kc, 3: ie, 4: ic, 5: panel, 6: object                        
  allocate(jcke(6,njcke,1)) ! 1: jc, 2: ke, 3: ie, 4: ic, 5: panel, 6: object
  allocate(jekc(6,njekc,1)) ! 1: je, 2: kc, 3: ie, 4: ic, 5: panel, 6: object                        
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  !---------------------------------------------------------------                                 
  !   1: i, 2: j, 3: z-coordinate, 4-5: tau, 7-9: beta, 10-12:n
  !   13-15: [du/dxi], 16-18:[dv/dxi], 19-21:[dw/dxi], 22-24:[dp/dxi]                              
  !   25-30: [du/dxidxj], 31-36: [dv/dxidxj], 37-42: [dw/dxidxj], 43-48: [dp/dxidxj]               
  !   49 : [p]                                                                                     
  !   50-55: [dp/dxidxj]                                                                           
  allocate(ficjc(76,nicjc,1))  !  1. ic, 2: jc, 3: z-coordinate,                  
  allocate(ficje(62,nicje,1))  !  1. ic, 2: je, 3: z-coordinate,                                  
  allocate(fiejc(62,niejc,1))  !  1. ie, 2: jc, 3: z-coordinate,                                     
  !   4-5: tau, 7-9: beta, 10-12:n                                                                 
  !   13-15: [du/dxi], 16-18:[dv/dxi], 19-21:[dw/dxi], 22-24:[dp/dxi]                              
  !   25-30: [du/dxidxj], 31-36: [dv/dxidxj], 37-42: [dw/dxidxj], 43-48: [dp/dxidxj]               
  !   49 : [p]                                                                                     
  !   50-55: [dp/dxidxj]                                                            
  allocate(fickc(76,nickc,1))      !   1: ic, 2: kc, 3: y-coordinate,                              
  allocate(ficke(62,nicke,1))      !   1: ic, 2: ke, 3: y-coordinate,                              
  allocate(fiekc(62,niekc,1))      !   1: ie, 2: kc, 3: y-coordinate,
  
  !   4-5: tau, 7-9: beta, 10-12:n                                                                 
  !   13-15: [du/dxi], 16-18:[dv/dxi], 19-21:[dw/dxi], 22-24:[dp/dxi]                              
  !   25-30: [du/dxidxj], 31-36: [dv/dxidxj], 37-42: [dw/dxidxj], 43-48: [dp/dxidxj]               
  !   49 : [p]                                                                                     
  !   50-55: [dp/dxidxj]                                                                           
  allocate(fjckc(76,njckc,1))  ! 1: jc, 2: kc, 3: x-coordinate                                     
  allocate(fjcke(62,njcke,1))  ! 1: jc, 2: ke, 3: x-coordinate                                     
  allocate(fjekc(62,njekc,1))  ! 1: je, 2: kc, 3: x-coordinate
  !---------------------------------------------------------------------
  allocate(kuiejc(niejc,1),kvicje(nicje,1),kwicjc(nicjc,1))
  allocate(kwiejc(niejc,1),kwicje(nicje,1),kuicje(nicje,1))
  allocate(juiekc(niekc,1),jwicke(nicke,1),jvickc(nickc,1))
  allocate(jviekc(niekc,1),jvicke(nicke,1),jwiekc(niekc,1))
  allocate(ivjekc(njekc,1),iwjcke(njcke,1),iujckc(njckc,1))
  allocate(iujekc(njekc,1),iujcke(njcke,1),ivjcke(njcke,1))


  allocate(ukiejc(niejc,1),vkicje(nicje,1),wkicjc(nicjc,1))
  allocate(wkiejc(niejc,1),wkicje(nicje,1),ukicje(nicje,1))
  allocate(ujiekc(niekc,1),wjicke(nicke,1),vjickc(nickc,1))
  allocate(vjiekc(niekc,1),vjicke(nicke,1),wjiekc(niekc,1))
  allocate(vijekc(njekc,1),wijcke(njcke,1),uijckc(njckc,1))
  allocate(uijekc(njekc,1),uijcke(njcke,1),vijcke(njcke,1))

  allocate(udz(nicjc,1),vdz(nicjc,1),wdz(nicjc,1))
  allocate(udy(nickc,1),vdy(nickc,1),wdy(nickc,1))
  allocate(udx(njckc,1),vdx(njckc,1),wdx(njckc,1))


  allocate(uwdz(niejc,1),vwdz(nicje,1),wwdz(nicjc,1),pdz(nicjc,1))
  allocate(udzz(2,niejc,1),vdzz(2,nicje,1),wdzz(2,nicjc,1),pdzz(2,nicjc,1))

  allocate(uvdy(niekc,1),wvdy(nicke,1),vvdy(nickc,1),pdy(nickc,1))
  allocate(udyy(2,niekc,1),wdyy(2,nicke,1),vdyy(2,nickc,1),pdyy(2,nickc,1))

  allocate(vudx(njekc,1),wudx(njcke,1),uudx(njckc,1),pdx(njckc,1))
  allocate(vdxx(2,njekc,1),wdxx(2,njcke,1),udxx(2,njckc,1),pdxx(2,njckc,1))  
  !-------------------------------------------------------------------------
end subroutine allocate_Lagrange
!----------------------------------------------------------------------------
