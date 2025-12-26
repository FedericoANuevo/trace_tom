pro function_power_law, x, A, f, pder

  Ninst = n_elements(A)/2

  if Ninst eq 1 then begin
          f = A[0] * x^(-A[1])
     df_dA0 = x^(-A[1])
     df_dA1 = -A[0] * x^(-A[1]) * alog(x)
     pder   = [ [df_dA0] , [df_dA1] ]
  endif
  
  if Ninst eq 2 then begin
          f = A[0] * x^(-A[1]) + A[2] * x^(-A[3])
     df_dA0 = x^(-A[1])
     df_dA1 = -A[0] * x^(-A[1]) * alog(x)
     df_dA2 = x^(-A[3])
     df_dA3 = -A[2] * x^(-A[3]) * alog(x)
     pder   = [ [df_dA0] , [df_dA1] , [df_dA2] , [df_dA3] ]
  endif
  
  if Ninst eq 3 then begin
     f = A[0] * x^(-A[1]) + A[2] * x^(-A[3]) + A[4] * x^(-A[5]) 
     df_dA0 = x^(-A[1])
     df_dA1 = -A[0] * x^(-A[1]) * alog(x)
     df_dA2 = x^(-A[3])
     df_dA3 = -A[2] * x^(-A[3]) * alog(x)
     df_dA4 = x^(-A[5])
     df_dA5 = -A[4] * x^(-A[5]) * alog(x)
     pder   = [ [df_dA0] , [df_dA1] , [df_dA2] , [df_dA3] , [df_dA4] , [df_dA5] ]
  endif

end
