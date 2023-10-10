pro read_fieldline
  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

; Path to use de x-tools  
 !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH  

; Read the DEMT 
  root_dir = '/data1/' &  dir  =root_dir+'DATA/ldem_files/'
  file ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' 
  read_ldem,file,/ldem,/gauss1,dir=dir


; directory and file-name of field line  
  dir  = '~/Downloads/'
  file = 'test2.dat';'test.dat'
  outfile = file+'.demt.dat'
; inicializo las variables
  xx = ''
  x  = 0.
  y  = 0.
  y  = 0.

; filas de datos a leer
  N1  = 10L
  N2  = 703L;10063L

; arrays 1D con la geometria de la linea 
  N = N2-N1+1
  x_l  = fltarr(N)
  y_l  = fltarr(N)
  z_l  = fltarr(N)
  r_l  = fltarr(N)
  th_l = fltarr(N)
  ph_l = fltarr(N)
; Para guardar Ne y Te a lo largo de la linea  
  Ne_l = fltarr(N)  
  Te_l = fltarr(N)
  
  
  
  openr,1,dir+file              ; abro el ASCII donde esta la linea de campo
  openw,2,dir+outfile
  for i = 1L, N1-1 do begin
     readf,1,xx
  endfor
  for i = N1,N2 do begin
     readf,1,x,y,z
     V = [x,y,z]
     cart_to_sphcoord,V,sphcoord
     r0  = sphcoord[0]
     th0 = sphcoord[1]
     ph0 = sphcoord[1]
     determindex,r0,th0,ph0,irad,ilat,ilon
     if irad ne -1 and ilat ne -1 and ilon ne -1 then begin
        Ne_l(i-N1) =  N_e(irad,ilat,ilon)
        Te_l(i-N1) =  Tm (irad,ilat,ilon)
        index      = irad + nr*ilat + nr*nth*ilon
     endif else begin
        Ne_l(i-N1) = -1.
        Te_l(i-N1) = -1.
        index      = -1
     endelse
     printf,2,x,y,z, Ne_l(i-N1),Te_l(i-N1),index
     x_l(i-N1) = x
     y_l(i-N1) = y
     z_l(i-N1) = z
     
     r_l (i-N1) = r0
     th_l(i-N1) = th0
     ph_l(i-N1) = ph0
  endfor
  close,1
  close,2
  

  STOP
  return
end
