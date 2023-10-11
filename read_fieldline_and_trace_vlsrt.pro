;
; PURPOSE: This code reads a single AWSoM fieldline tracing file,
;          traces th VL-SRT 3D N_e along it,
;          and produces an expanded fieldline tracing file including
;          the tomographic electron density as a new column.
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;

pro read_fieldline_and_trace_vlsrt,instr,dir,file,rad,lat,lon,nr,nt,np,N_e
; Define name of output file
  outfile = file+'_'+instr+'.out'  
; Read the ASCII fieldline file with readcol (SolarSoft)  
  readcol,dir+file,x_l,y_l,z_l,FORMAT='D'  
; elements of 1D vector with field-line coord. 
  N    = n_elements(x_l)
; Define the fieldline spherical coordinate arrays  
  r_l  = dblarr(N)
  th_l = dblarr(N)
  ph_l = dblarr(N)
; Define the fieldline tomographic variables arrays
  Ne_l  = fltarr(N)  
  Te_l  = fltarr(N)
  ind_l = lonarr(N)*0L
; Trace tomographic products along the fieldline  
  for i = 0,N-1 do begin
     x = x_l(i)
     y = y_l(i)
     z = z_l(i)
     V = [x,y,z]
     cart_to_sphcoord,V,sphcoord
     r0  = sphcoord[0]
     th0 = sphcoord[1]
     ph0 = sphcoord[1]
     determindex,r0,th0,ph0,irad,ilat,ilon,rad,lat,lon
     if irad ne -1 and ilat ne -1 and ilon ne -1 then begin
        Ne_l (i) =  N_e (irad,ilat,ilon)
        ind_l(i) =  irad + nr*ilat + nr*nt*ilon
     endif else begin
        Ne_l (i) = -1.
        ind_l(i) = -1
     endelse
  endfor
; Write output file
  openw,2,dir+outfile
  printf,2,'  X [Rs]            Y [Rs]            Z [Rs]            Ne [cm^ -3]           Cell-3D-index'
  for i = 0,N-1 do begin
     printf,2,x_l(i),y_l(i),z_l(i), Ne_l(i),ind_l(i),$
            FORMAT = '(4(E18.10)," ",I9)' 
  endfor
  close,2
  return
end
