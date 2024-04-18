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
  outfile = file+'.demt.out'

  readcol,dir+file,x_l,y_l,z_l
  N    = n_elements(x_l)
; Geometria de la linea en coord. esfericas
  r_l  = fltarr(N)
  th_l = fltarr(N)
  ph_l = fltarr(N)
; Para guardar Ne y Te a lo largo de la linea  
  Ne_l  = fltarr(N)  
  Te_l  = fltarr(N)
  ind_l = intarr(N)*0L
  
  for i = 0,N-1 do begin
     x = x_l(i)
     y = y_l(i)
     z = z_l(i)
     V = [x,y,z]
     cart_to_sphcoord,V,sphcoord
     r0  = sphcoord[0]
     th0 = sphcoord[1]
     ph0 = sphcoord[1]
     determindex,r0,th0,ph0,irad,ilat,ilon
     if irad ne -1 and ilat ne -1 and ilon ne -1 then begin
        Ne_l (i) =  N_e(irad,ilat,ilon)
        Te_l (i) =  Tm (irad,ilat,ilon)
        ind_l(i) =  irad + nr*ilat + nr*nth*ilon
     endif else begin
        Ne_l (i) = -1.
        Te_l (i) = -1.
        ind_l(i) = -1
     endelse
  endfor


  openw,1,dir+outfile
  for i = 0,N-1 do begin
     printf,1,x_l(i),y_l(i),z_l(i), Ne_l(i),Te_l(i),ind_l(i)
  endfor
  close,1
  
  STOP
  return
end
