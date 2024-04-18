function r2_function,y,yfit
  SStot = total( (y-mean(y))^2 )
  SSres = total( (y-yfit   )^2 )
  r2    = 1.-SSres/SStot
  return,r2
end
