pro double_power_fit, radsamp, Nesamp, A, chisq, noweight=noweight, weighted=weighted
  
     rad1 = min(radsamp)
     rad2 = max(radsamp)
  radcrit = rad1 + (rad2-rad1)/3. ; The value "2" is based on inspection of specific cases.
                                  ; Tests indicate it provides better results than other
                                  ; reasonable values (f = 3, 4).

        iA = where(radsamp lt radcrit)
        iB = where(radsamp ge radcrit)
  radsampA = radsamp(iA)
  radsampB = radsamp(iB)
   NesampA =  Nesamp(iA)
   NesampB =  Nesamp(iB)

; Initial estimate for A:
  A = fltarr(4)
  linear_fit, alog(radsampA), alog(NesampA), AN, r2N, /linfit_idl
  A[0] = exp(AN[0]) ; cm-3
  A[1] =    -AN[1]  ; dimensionless exponent of power law
  linear_fit, alog(radsampB), alog(NesampB), AN, r2N, /linfit_idl
  A[2] = exp(AN[0]) ; cm-3
  A[3] =    -AN[1]  ; dimensionless exponent of power law

; Set weights:
  if keyword_set(noweight) then weights = 0.*Nesamp + 1.
  if keyword_set(weighted) then weights =    Nesamp

; Fit:
  yfit_func = CURVEFIT(radsamp, Nesamp, weights, A, SIGMA, function_name='function_double_power_law', status=status, iter=iter, chisq=chisq, itmax=100, /double)

  test = finite(A,/NAN)
  if test(0) eq 1 then stop ; sucede si el determinante es nulo
  
  return
end
