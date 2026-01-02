pro double_power_concat_fit, x_concat, y_concat, inst_concat, A, chisq, weighted=weighted
  
 ;rad1 = min(x_concat)
 ;rad2 = max(x_concat)
 ;radcrit = rad1 + (rad2-rad1)/2. ; The value "2" is based on inspection of specific cases.
                                  ; Tests indicate it provides better results than other
                                  ; reasonable values (f = 3, 4).
  xcrit = 2.0
  
  iA = where(x_concat lt xcrit)
  iB = where(x_concat ge xcrit)
  x_concatA = x_concat(iA)
  x_concatB = x_concat(iB)
  y_concatA = y_concat(iA)
  y_concatB = y_concat(iB)

; Initial estimate for A:
  A = fltarr(4)
  linear_fit, alog(x_concatA), alog(y_concatA), AN, r2N, /linfit_idl
  A[0] = exp(AN[0])             ; cm-3
  A[1] =    -AN[1]              ; dimensionless exponent of power law
  linear_fit, alog(x_concatB), alog(y_concatB), AN, r2N, /linfit_idl
  A[2] = exp(AN[0])             ; cm-3
  A[3] =    -AN[1]              ; dimensionless exponent of power law

; Set weights:
  if NOT keyword_set(weighted) then weights = 0.*y_concat + 1.
  if     keyword_set(weighted) then weights = 1./(y_concat / mean(y_concat))^2
; Do not consider data points above xmax_fit
  xmax_fit = 5.5
  index = where(x_concat ge xmax_fit) & if index(0) ne -1 then weights(index)=0.

; Fit:
  yfit_func = CURVEFIT(x_concat, y_concat, weights, A, SIGMA,$
                       function_name='function_double_power_law',$
                       status=status, iter=iter, chisq=chisq, itmax=100, /double)

 ;test = finite(A,/NAN)
 ;if test(0) eq 1 then stop ; sucede si el determinante es nulo
  
  return
end
