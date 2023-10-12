pro synth_radial_fieldline, lat0 = lat0 , lon0 = lon0, dr = dr, dir = dir, outfile = outfile

  if not keyword_set(dir)     then  dir     = '/data1/DATA/fieldlines_judit/'
  if not keyword_set(outfile) then  outfile = 'synth_radial_fieldline.txt'  
  if not keyword_set(lat0)    then  lat0    = 70.
  if not keyword_set(lon0)    then  lon0    =  0.
  if not keyword_set(dr)      then  dr      =  0.005

  th0  = (90. - lat0)*!dtor
  ph0  = lon0        *!dtor
  rmin = 1.
  rmax = 10.
  N    = fix((rmax-rmin)/dr)
  
  r_l  = rmin + dr*findgen(N) + dr/2
  th_l = fltarr(N) + th0
  ph_l = fltarr(N) + ph0

  x_l = r_l * sin(th_l) * cos(ph_l)
  y_l = r_l * sin(th_l) * sin(ph_l)
  z_l = r_l * cos(th_l) 

; Write output file
  openw,2,dir+outfile
  printf,2,'  X [Rs]            Y [Rs]            Z [Rs]            '
    for i = 0,N-1 do begin
     printf,2,x_l(i),y_l(i),z_l(i), FORMAT = '(3(E18.10))' 
  endfor
  close,2


     

  
  
  return
end
