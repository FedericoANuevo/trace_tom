pro linear_fit, xsamp, ysamp, A, r2, $
                linfit_idl = linfit_idl, theil_sen = theil_sen

  if keyword_set(linfit_idl) then A =   linfit(xsamp, ysamp, prob=prob, /double)
  if keyword_set(theil_sen)  then A = theilsen(xsamp, ysamp)

  r2 = r2_function(ysamp,A[0]+A[1]*xsamp)

  return
end
