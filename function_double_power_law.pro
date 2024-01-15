pro function_double_power_law, r, A, f, dfdA

  common radcrits, radcritA, radcritB

  f = A[0] * (r/radcritA)^(-A[1]) + A[2] * (r/radcritB)^(-A[3])

  df_dA0 = (r/radcritA )^(-A[1])
  df_dA1 = -A[0] * (r/radcritA )^(-A[1]) * alog(r/radcritA )
  df_dA2 = (r/radcritB)^(-A[3])
  df_dA3 = -A[2] * (r/radcritB)^(-A[3]) * alog(r/radcritB)
  
  df_dA  = [ [df_dA0] , [df_dA1] , [df_dA2] , [df_dA3] ]

  return
end
