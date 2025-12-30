pro power_concat_fit, Inst_list, rad_concat, Ne_concat, inst_concat, A, chisq, weighted=weighted

; Compute an initial set of values for the fit parameters.
  Ninst = n_elements(Inst_list)
  A     = fltarr(2*Ninst)
  for inst = 0,Ninst-1 do begin
     index = where(inst_concat eq inst_list(inst))
     rad_fit = rad_concat(index)
      Ne_fit =  Ne_concat(index)
     linear_fit, alog(rad_fit), alog(Ne_fit), AN, r2N, /linfit_idl
    A[2*inst  ] = exp(AN[0]) ; cm-3
    A[2*inst+1] =    -AN[1]  ; dimensionless exponent of power law
  endfor

; Set weights:
  if NOT keyword_set(weighted) then weights = 0.*Ne_concat + 1.
  if     keyword_set(weighted) then weights =  1./  (Ne_concat / mean(Ne_concat))^2
; Do not consider data points above Rmax_fit
  Rmax_fit = 5.5
  index = where(rad_concat ge Rmax_fit) & if index(0) ne -1 then weights(index)=0.
  
; Fit:
  yfit_func = CURVEFIT(rad_concat, Ne_concat, weights, A, SIGMA,$
                       function_name='function_power_law',$
                       status=status, iter=iter, chisq=chisq, itmax=100, /double)

 ;test = finite(A,/NAN)
 ;if test(0) eq 1 then stop ; sucede si el determinante es nulo
  
  return
end
