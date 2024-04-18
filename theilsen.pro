; Performs theil-sen estimation, avoiding pairs with same abscissa.
function theilsen, x, y
  x = float(x)
  y = float(y)
  num = n_elements(x)
  n=float(num)*(num-1)/2
  theil=fltarr(n)
  t=0.
  for i=0,num-2 do begin
     for j=i+1,num-1 do begin
        if x[j] ne x[i] then theil[t] = (y[j]-y[i]) / (x[j]-x[i])
        if x[j] eq x[i] then theil[t] = -666.
        t=t+1
     endfor
  endfor
  theil     = theil( where(theil ne -666.) )
  slope     = median(theil    , /even)
  intercept = median(y-slope*x, /even)
  results   = [intercept, slope]
  return, results
end
