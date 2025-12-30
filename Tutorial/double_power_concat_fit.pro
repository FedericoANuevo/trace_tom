pro double_power_concat_fit, rad_concat, Ne_concat, inst_concat, A, chisq, weighted=weighted
  
 ;rad1 = min(rad_concat)
 ;rad2 = max(rad_concat)
 ;radcrit = rad1 + (rad2-rad1)/2. ; The value "2" is based on inspection of specific cases.
                                  ; Tests indicate it provides better results than other
                                  ; reasonable values (f = 3, 4).
  radcrit = 2.0
  
  iA = where(rad_concat lt radcrit)
  iB = where(rad_concat ge radcrit)
  rad_concatA = rad_concat(iA)
  rad_concatB = rad_concat(iB)
  Ne_concatA  =  Ne_concat(iA)
  Ne_concatB  =  Ne_concat(iB)

; Initial estimate for A:
  A = fltarr(4)
  linear_fit, alog(rad_concatA), alog(Ne_concatA), AN, r2N, /linfit_idl
  A[0] = exp(AN[0])             ; cm-3
  A[1] =    -AN[1]              ; dimensionless exponent of power law
  linear_fit, alog(rad_concatB), alog(Ne_concatB), AN, r2N, /linfit_idl
  A[2] = exp(AN[0])             ; cm-3
  A[3] =    -AN[1]              ; dimensionless exponent of power law

; Set weights:
  if NOT keyword_set(weighted) then weights = 0.*Ne_concat + 1.
  if     keyword_set(weighted) then weights = 1./(Ne_concat / mean(Ne_concat))^2
; Do not consider data points above Rmax_fit
  Rmax_fit = 5.5
  index = where(rad_concat ge Rmax_fit) & if index(0) ne -1 then weights(index)=0.

; Fit:
  yfit_func = CURVEFIT(rad_concat, Ne_concat, weights, A, SIGMA,$
                       function_name='function_double_power_law',$
                       status=status, iter=iter, chisq=chisq, itmax=100, /double)

 ;test = finite(A,/NAN)
 ;if test(0) eq 1 then stop ; sucede si el determinante es nulo
  
  return
end
