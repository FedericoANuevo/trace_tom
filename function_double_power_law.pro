pro function_double_power_law, x, A, f, pder

  f = A[0] * x^(-A[1]) + A[2] * x^(-A[3])

  df_dA0 = x^(-A[1])
  df_dA1 = -A[0] * x^(-A[1]) * alog(x)
  df_dA2 = x^(-A[3])
  df_dA3 = -A[2] * x^(-A[3]) * alog(x)
  
  pder  = [ [df_dA0] , [df_dA1] , [df_dA2] , [df_dA3] ]

end
