;
; PURPOSE: This code reads a single AWSoM fieldline tracing file
;          (which is esentially an ASCII table),
;          traces DEMT 3D products (N_e, T_e) along it,
;          and produces an expanded fieldline tracing file including
;          the tomographic producs as new columns in the table.
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;

pro read_fieldline_and_trace_demt,instr,dir,file,rad,lat,lon,nr,nt,np,N_e,T_e,W_T,ldem_flag,csv=csv,trace_Bs=trace_Bs
; Define name of output file
  outfile = file+'_'+instr+'.out'
  if not keyword_set(csv) then begin
;    Read the ASCII fieldline file with readcol (SolarSoft) 
     if NOT keyword_set(trace_Bs) then $
     readcol,dir+file,x_l,y_l,z_l,FORMAT='D,D,D'
     if     keyword_set(trace_Bs) then $
     readcol,dir+file,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,FORMAT='D,D,D,D,D,D,D,D'
  endif else begin
     data = read_csv(dir+file,header = header)
     x_l  = reform(data.field2)
     y_l  = reform(data.field3)
     z_l  = reform(data.field4)
  endelse
  
; elements of 1D vector with field-line coord. 
  N    = n_elements(x_l)
; Define the fieldline spherical coordinate arrays  
  r_l  = dblarr(N)
  th_l = dblarr(N)
  ph_l = dblarr(N)
; Define the fieldline tomographic variables arrays
  Ne_l  = fltarr(N)  
  Te_l  = fltarr(N)
  WT_l  = fltarr(N)
  ldem_flag_l = fltarr(N)
  ind_l = lonarr(N)
; Trace tomographic products along the fieldline  
  for i = 0,N-1 do begin
     x = x_l(i)
     y = y_l(i)
     z = z_l(i)
     V = [x,y,z]
     cart_to_sphcoord,V,sphcoord
     r0  = sphcoord[0] & r_l (i) = r0
     th0 = sphcoord[1] & th_l(i) = th0
     ph0 = sphcoord[2] & ph_l(i) = ph0
     determindex,r0,th0,ph0,irad,ilat,ilon,rad,lat,lon
     if irad ne -1 and ilat ne -1 and ilon ne -1 then begin
        Ne_l (i) =  N_e (irad,ilat,ilon)
        Te_l (i) =  T_e (irad,ilat,ilon)
        WT_l (i) =  W_T (irad,ilat,ilon)
        ldem_flag_l(i) = ldem_flag(irad,ilat,ilon)
        ind_l(i) =  irad + nr*ilat + nr*nt*ilon
     endif else begin
        Ne_l (i) = -1.
        Te_l (i) = -1.
        WT_l (i) = -1.
        ldem_flag_l(i) = -1.
        ind_l(i) = -1
     endelse
  endfor
; Write output file
  openw,2,dir+outfile
  if NOT keyword_set(trace_Bs) then begin
     printf,2,'  X [Rs]            Y [Rs]            Z [Rs]                  r [Rs]        lat [deg]          lon [deg]        Ne [cm^ -3]       Te [K]               WT [K]          LDEM-flag         Cell-3D-index'
     for i = 0,N-1 do begin
        printf,2,x_l(i),y_l(i),z_l(i),r_l(i),90. -th_l(i)/!dtor, ph_l(i)/!dtor, Ne_l(i),Te_l(i),WT_l(i),ldem_flag_l(i),ind_l(i),$
               FORMAT = '(10(E18.10)," ",I9)'
     endfor
  endif
  if     keyword_set(trace_Bs) then begin
     printf,2,'  X [Rs]            Y [Rs]            Z [Rs]            s [Rs]            B_r [G]           B_th [G]          B_ph [G]          B [G]             r [Rs]            lat [deg]         lon [deg]         Ne [cm^ -3]       Te [K]            WT [K]            LDEM-flag         Cell-3D-index'
     for i = 0,N-1 do begin
        printf,2,x_l(i),y_l(i),z_l(i),s_l(i),Br_l(i),Bth_l(i),Bph_l(i),B_l(i),r_l(i),90. -th_l(i)/!dtor, ph_l(i)/!dtor, Ne_l(i),Te_l(i),WT_l(i),ldem_flag_l(i),ind_l(i),$
               FORMAT = '(15(E18.10)," ",I9)'
     endfor
  endif
  close,2
  return
end
