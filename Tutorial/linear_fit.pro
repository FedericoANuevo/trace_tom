; Performs linear fit to data in different selectable ways.

pro linear_fit, xsamp, ysamp, A, r2, chisqr = chisqr, $
                linfit_idl = linfit_idl, theil_sen = theil_sen

  if keyword_set(linfit_idl) then A =   linfit(xsamp, ysamp, prob=prob, /double)
  if keyword_set(theil_sen)  then A = theilsen(xsamp, ysamp)


  chisqr = total((ysamp - (A[0]+A[1]*xsamp))^2 )/(n_elements(ysamp)-n_elements(A))
  
  r2 = r2_function(ysamp,A[0]+A[1]*xsamp)

  return
end
