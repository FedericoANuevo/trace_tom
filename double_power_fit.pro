pro double_power_fit, radmin, radmax, radsamp, Nesamp, A, chisq

  common radcrits, radcritA, radcritB
  
  A = fltarr(4)
  
  radcritA = radmin 
 ;radcritB = radmin + (radmax-radmin)/2.
  radcritB = 1.2
  
        iA = where(radsamp lt radcritB)
        iB = where(radsamp ge radcritB)
  radsampA = radsamp(iA)
  radsampB = radsamp(iB)
   NesampA =  Nesamp(iA)
   NesampB =  Nesamp(iB)
  
  linear_fit, alog(radsampA/radcritA), alog(NesampA), AN, r2N, /linfit_idl
  A[0] = exp(AN[0]) ; cm-3
  A[1] =    -AN[1]  ; dimensionless exponent of power law

  linear_fit, alog(radsampB/radcritB), alog(NesampB), AN, r2N, /linfit_idl
  A[2] = exp(AN[0]) ; cm-3
  A[3] =    -AN[1]  ; dimensionless exponent of power law

  yfit_func = CURVEFIT(radsamp, Nesamp, weights, A, SIGMA, FUNCTION_NAME='function_double_power_law',STATUS=status,iter=iter,chisq=chisq,itmax=2000)

  test = finite(A,/NAN)
  if test(0) eq 1 then stop ; sucede si el determinante es nulo

  return
end
