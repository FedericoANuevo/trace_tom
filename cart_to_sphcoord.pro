;
; PURPOSE: This code takes in a 3D cartesian positio vector
; V (1D 3-element array) 
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;


; Esta rutina da las coordenas esféricas de un vector V en cartesianas
; V=(x,y,z) > SPHCOORD=(r,th,ph)
; V es un vector (de entrada)  de 3 componentes x, y, z [Rsun]  
; sphcoord es el vector salida de 3 componentes r[Rsun], th [RAD], ph [RAD]


pro cart_to_sphcoord,V,sphcoord
  x=V[0]*1d
  y=V[1]*1d
  z=V[2]*1d
  r  = sqrt(x^2+y^2+z^2)
  th = acos(z/r)
; devuelve siempre un valor de ph entre 0 y 2*!pi
; x = 3. & r = x/cos(30.*!dtor) & y = r*sin(30.*!dtor)
; x = +abs(x) & y = +abs(y) 
  ph =          atan(y/x)
  if x gt 0 and y lt 0. then ph = 2*!dpi + atan(y/x)
  if x lt 0             then ph =   !dpi - atan(y/abs(x))
  sphcoord=[r,th,ph]
  return
end
