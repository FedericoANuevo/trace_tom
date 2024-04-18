pro test_index

  n1 = 3 & i0 = 1
  n2 = 5 & j0 = 3
  n3 = 7 & k0 = 6

  array           = fltarr(n1,n2,n3)
  array(i0,j0,k0) = -666.

  ind_wh = where(array eq -666.)
  ixdex  = 0L
  index  = i0 + n1*j0 + n1*n2*k0

  print,'index from where'  ,ind_wh 
  print,'index from formula',index

; Recover the 3D-index   
  i = index    mod n1
  j = index/n1 mod n2
  k = index/(n1*n2)
  
  help,i ,j ,k  
  help,i0,j0,k0  
  
  STOP
  return
end
