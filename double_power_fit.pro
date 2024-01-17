pro double_power_fit, rad1, rad2, radsamp, Nesamp, A, chisq
  
  A = fltarr(4)
  
  radcrit = rad1 + (rad2-rad1)/2. ; a dynamic value based on inspection of specific cases
  
        iA = where(radsamp lt radcrit)
        iB = where(radsamp ge radcrit)
  radsampA = radsamp(iA)
  radsampB = radsamp(iB)
   NesampA =  Nesamp(iA)
   NesampB =  Nesamp(iB)
  
  linear_fit, alog(radsampA), alog(NesampA), AN, r2N, /linfit_idl
  A[0] = exp(AN[0]) ; cm-3
  A[1] =    -AN[1]  ; dimensionless exponent of power law

  linear_fit, alog(radsampB), alog(NesampB), AN, r2N, /linfit_idl
  A[2] = exp(AN[0]) ; cm-3
  A[3] =    -AN[1]  ; dimensionless exponent of power law

  weights = 0.*radsamp + 1. ; No weights
 ;weights = 1./Nesamp
  
  yfit_func = CURVEFIT(radsamp, Nesamp, weights, A, SIGMA, function_name='function_double_power_law', status=status, iter=iter, chisq=chisq, itmax=100, /double)

  test = finite(A,/NAN)
  if test(0) eq 1 then stop ; sucede si el determinante es nulo
  
  return
end
