;
; PURPOSE: This code interpolates yv(xv) of a fieldline into yi(xi),
; where xi is a common fine grid.
;
; INPUTS:
; xv and yv: FLTARRAYS with the yv(xv) data points.
; instrumental flag.
;
; OUTPUTS:
; xi and yi, the fine common grid and interpolated values.
;
; HISTORY: v1.0: AMV, November 2023, IAFE.
;
pro interpol_fl, xv=xv, yv=yv, xi=xi, yi=yi, $
                 aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                 mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

; Set parameters for fine common grid, dependening on instrument:
  if keyword_set(aia) then begin
     xi_min = 1.02
     xi_max = 1.25
     Nxi    = 50
  endif
  if keyword_set(mk4) then begin
     xi_min = 1.15
     xi_max = 1.50
     Nxi    = 100
  endif
  if keyword_set(lascoc2) then begin
     xi_min = 2.50
     xi_max = 6.00
     Nxi    = 450
  endif
  
; Set-up fine common grid:
  dxi = (xi_max-xi_min)/Nxi  
  xi  = xi_min + dxi/2. + dxi * findgen(Nxi)

; Interpol into fine common grid:
  yi  = INTERPOL( yv, xv, xi, /spline )

  return
end
