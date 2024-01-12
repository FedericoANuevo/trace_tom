;
; PURPOSE: This code merges the tracing of tomographic products based
; on data from different instruments along given AWSoM fieldlines, and
; store all the information in a pointer structure.
;
; It also performs, if requested (/fits) analytical fits to the
; various tomographic products along individual field lines, including
; the fits as part of the structure.
;
; INPUTS:
; fir_fl and fl_list: STRINGS. directory where field lines are
; located, and filename of list of names of field lines filenames.
;
; FLAGS: one per allowed instrument, use those for which you want to
; merge results, provided in any order.
; Also /fits.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
;          V1.1 AMV, IAFE, November 2023. Added Sampling, expanded to
;          all instruments and performed overall polishing.
;          V2.0 AMV, IAFE, January 2024. Added fit to tomographic resuls.
;

pro merge_trace_struct, dir_fl = dir_fl, fl_list = fl_list, $
                        aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                        mk4 = mk4, kcor = kcor, lascoc2 = lascoc2, $
                        struture_filename = structure_filename, $
                        fits = fits
  
  common all_data, trace_data, $
     N_fl, Npt_max, Npt_v,$
     x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A
  
  if not keyword_set(dir_fl) or not keyword_set(fl_list) then STOP

; Set up filename for output structure:
  if not keyword_set(structure_filename) then structure_filename = fl_list
  structure_filename = structure_filename + '-tracing-structure-merge'

; Read the list with the field-lines  
  N_fl     = 0
  filename = ''
  openr,1,dir_fl+fl_list
  readf,1,N_fl

; Maximum number of point along the fieldline  
  Npt_max = 10100 ; Ask JUDIT for an optimal value.
; Default value in all arrays.
  default = -678.
  
; Define needed arrays:
  x_A   = dblarr(N_fl,Npt_max) + default
  y_A   = dblarr(N_fl,Npt_max) + default
  z_A   = dblarr(N_fl,Npt_max) + default
  rad_A = dblarr(N_fl,Npt_max) + default
  lat_A = dblarr(N_fl,Npt_max) + default
  lon_A = dblarr(N_fl,Npt_max) + default
  Npt_v = intarr(N_fl)      
  if keyword_set(aia) then begin
     Ne_aia_A    = fltarr(N_fl,Npt_max) + default
     Tm_aia_A    = fltarr(N_fl,Npt_max) + default
     index_aia_A = intarr(N_fl,Npt_max) + default
     index_sampling_aia_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(euvia) then begin
     Ne_euvia_A    = fltarr(N_fl,Npt_max) + default
     Tm_euvia_A    = fltarr(N_fl,Npt_max) + default
     index_euvia_A = intarr(N_fl,Npt_max) + default
     index_sampling_euvia_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(euvib) then begin
     Ne_euvib_A    = fltarr(N_fl,Npt_max) + default
     Tm_euvib_A    = fltarr(N_fl,Npt_max) + default
     index_euvib_A = intarr(N_fl,Npt_max) + default
     index_sampling_euvib_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(eit) then begin
     Ne_eit_A    = fltarr(N_fl,Npt_max) + default
     Tm_eit_A    = fltarr(N_fl,Npt_max) + default
     index_eit_A = intarr(N_fl,Npt_max) + default
     index_sampling_eit_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(mk4) then begin
     Ne_mk4_A    = fltarr(N_fl,Npt_max) + default
     index_mk4_A = intarr(N_fl,Npt_max) + default
     index_sampling_mk4_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(kcor) then begin
     Ne_kcor_A    = fltarr(N_fl,Npt_max) + default
     index_kcor_A = intarr(N_fl,Npt_max) + default
     index_sampling_kcor_A = intarr(N_fl,Npt_max) + default
  endif
  if keyword_set(lascoc2) then begin
     Ne_c2_A    = fltarr(N_fl,Npt_max) + default
     index_c2_A = intarr(N_fl,Npt_max) + default
     index_sampling_c2_A = intarr(N_fl,Npt_max) + default
  endif
           
  for i_fl = 0,N_fl-1 do begin ; Start loop in fieldlines  
     initialized = 'no'
     readf,1,filename
     outfile       = filename + '_merge'
     
   ; AIA
     if keyword_set(aia)   then begin
        file_aia   = filename+'_aia.out'
        outfile    = outfile +'_aia'
        if i_fl eq 0 then structure_filename = structure_filename + '_aia'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_aia_l  ,Tm_aia_l  ,index_aia_l  ,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_aia_l, index_l=index_aia_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_aia_l],[Tm_aia_l],[index_aia_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [AIA, cm^-3]   Tm [AIA, K]        AIA-3DInd/SampInd'
           output_format = '(8(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
        endif
        Ne_aia_A   (i_fl,0:N_l-1)          = Ne_aia_l
        Tm_aia_A   (i_fl,0:N_l-1)          = Tm_aia_l
        index_aia_A(i_fl,0:N_l-1)          = index_aia_l
        index_sampling_aia_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; EUVI-A
     if keyword_set(euvia)  then begin
        file_euvia = filename+'_euvia.out'
        outfile   =  outfile +'_euvia'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvia'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,index_euvia_l,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvia_l, index_l=index_euvia_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_euvia_l],[Tm_euvia_l],[index_euvia_l],[index_sampling_l]]
           header_str = header_str +'         Ne [EUVIA, cm^-3]  Tm [EUVIA, K]     EUVIA-3DInd      EUVIA-Samp'
           output_format = output_format + ',"  ",2(E18.10),"  ",2(I9)' 
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_euvia_l],[Tm_euvia_l],[index_euvia_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EUVIA, cm^-3] Tm [EUVIA, K]    EUVIA-3DInd      EUVIA-Samp'
           output_format = '(8(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
        endif
        Ne_euvia_A   (i_fl,0:N_l-1) = Ne_euvia_l
        Tm_euvia_A   (i_fl,0:N_l-1) = Tm_euvia_l
        index_euvia_A(i_fl,0:N_l-1) = index_euvia_l    
        index_sampling_euvia_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; EUVI-B
     if keyword_set(euvib)  then begin
        file_euvib = filename+'_euvib.out'
        outfile    = outfile +'_euvib'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvib'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,index_euvib_l,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvib_l, index_l=index_euvib_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_euvib_l],[Tm_euvib_l],[index_euvib_l],[index_sampling_l]]
           header_str = header_str +'         Ne [EUVIB, cm^-3]  Tm [EUVIB, K]     EUVIB-3DInd      EUVIB-Samp'
           output_format = output_format + ',"  ",2(E18.10),"  ",2(I9)'
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_euvib_l],[Tm_euvib_l],[index_euvib_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EUVIB, cm^-3] Tm [EUVIB, K]    EUVIB-3DInd      EUVIB-Samp'
           output_format = '(8(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
        endif
        Ne_euvib_A   (i_fl,0:N_l-1) = Ne_euvib_l
        Tm_euviB_A   (i_fl,0:N_l-1) = Tm_euvib_l
        index_euvib_A(i_fl,0:N_l-1) = index_euvib_l
        index_sampling_euvib_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; EIT
     if keyword_set(eit)    then begin
        file_euvib = filename+'_eit.out'
        outfile    = outfile +'_eit'
        if i_fl eq 0 then structure_filename = structure_filename + '_eit'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,index_eit_l,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_eit_l, index_l=index_eit_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_eit_l],[Tm_eit_l],[index_eit_l],[index_sampling_l]]
           header_str = header_str +'         Ne [EIT, cm^-3]  Tm [EIT, K]     EIT-3DInd      EIT-Samp'
           output_format = output_format + ',"  ",2(E18.10),"  ",2(I9)'
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_eit_l],[Tm_eit_l],[index_eit_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EIT, cm^-3]   Tm [EIT, K]        EIT-3DInd      EIT-Samp'
           output_format = '(8(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
        endif
        Ne_eit_A   (i_fl,0:N_l-1) = Ne_eit_l
        Tm_eit_A   (i_fl,0:N_l-1) = Tm_eit_l
        index_eit_A(i_fl,0:N_l-1) = index_eit_l
        index_sampling_eit_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; Mk4
     if keyword_set(mk4)    then begin
        file_mk4   = filename+'_mk4.out'
        outfile    = outfile +'_mk4'
        if i_fl eq 0 then structure_filename = structure_filename + '_mk4'
        readcol,dir_fl+file_mk4,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_mk4_l             ,index_mk4_l  ,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_mk4_l, index_l=index_mk4_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_mk4_l],[index_mk4_l],[index_sampling_l]] 
           header_str = header_str +'   Ne [Mk4, cm^-3]     Mk4-3DInd/SampInd'
           output_format = output_format + ',"  ",E18.10,"  ",2(I9)'
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_mk4_l],[index_mk4_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [Mk4, cm^-3]    Mk4-3DInd/SampInd'
           output_format = '(7(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
        endif
        Ne_mk4_A   (i_fl,0:N_l-1) = Ne_mk4_l
        index_mk4_A(i_fl,0:N_l-1) = index_mk4_l
        index_sampling_mk4_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; KCOR
     if keyword_set(kcor)    then begin
        file_kcor   = filename+'_kcor.out'
        outfile     = outfile +'_kcor'
        if i_fl eq 0 then structure_filename = structure_filename + '_kcor'
        readcol,dir_fl+file_kcor,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_kcor_l           ,index_kcor_l ,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_kcor_l, index_l=index_kcor_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_kcor_l],[index_kcor_l],[index_sampling_l]] 
           header_str = header_str +'   Ne [KCOR, cm^-3]    KCOR-3DInd/SampInd'
           output_format = output_format + ',"  ",E18.10,"  ",2(I9)'
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_kcor_l],[index_kcor_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [KCOR, cm^-3]   KCOR-3DInd/SampInd'
           output_format = '(7(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
        endif
        Ne_kcor_A   (i_fl,0:N_l-1) = Ne_kcor_l
        index_kcor_A(i_fl,0:N_l-1) = index_kcor_l        
        index_sampling_kcor_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; LASCO-C2
     if keyword_set(lascoc2) then begin
        file_c2    = filename+'_lascoc2.out'
        outfile    = outfile +'_lascoc2'
        if i_fl eq 0 then structure_filename = structure_filename + '_lascoc2'
        readcol,dir_fl+file_c2 ,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_c2_l              ,index_c2_l   ,FORMAT='D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_c2_l, index_l=index_c2_l, index_sampling_l=index_sampling_l
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_c2_l],[index_c2_l],[index_sampling_l]] 
           header_str = header_str +'   Ne [C2, cm^-3]      C2-3DInd/SampInd'
           output_format = output_format + ',"  ",E18.10,"  ",2(I9)'
        endif
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_c2_l],[index_c2_l],[index_sampling_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [C2, cm^-3]     C2-3DInd/SampInd'
           output_format = '(7(E18.10)," ",2(I9)'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
        endif
        Ne_c2_A   (i_fl,0:N_l-1) = Ne_c2_l
        index_c2_A(i_fl,0:N_l-1) = index_c2_l    
        index_sampling_c2_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; Close output filename
     outfile = outfile + '.out'
   ; Close output format string
     output_format = output_format + ')'
   ; Transpose output_columns array
     output_columns = transpose(output_columns)
   ; Write output file
      openw,2,dir_fl+outfile
     printf,2,header_str
     for i = 0,n_elements(x_l)-1 do printf,2,output_columns(*,i),FORMAT=output_format
     close,2
  endfor  ; End loop in fieldlines    
  close,1

; Create a pointer structure to store traced information  
  trace_data = { N_fl:        ptr_new(N_fl)                          ,$
                 Npt_max:     ptr_new(Npt_max)                       ,$
                 Npt_v:       ptr_new(Npt_v)                         ,$
                 x:           ptr_new(x_A)                           ,$
                 y:           ptr_new(y_A)                           ,$
                 z:           ptr_new(z_A)                           ,$
                 rad:         ptr_new(rad_A)                         ,$
                 lat:         ptr_new(lat_A)                         ,$
                 lon:         ptr_new(lon_A)                         ,$
                 Ne_aia:      ptr_new()                              ,$
                 Tm_aia:      ptr_new()                              ,$
                 index_aia:   ptr_new()                              ,$
                 Ne_euvia:    ptr_new()                              ,$
                 Tm_euvia:    ptr_new()                              ,$
                 index_euvia: ptr_new()                              ,$
                 Ne_euvib:    ptr_new()                              ,$
                 Tm_euvib:    ptr_new()                              ,$
                 index_euvib: ptr_new()                              ,$
                 Ne_eit:      ptr_new()                              ,$
                 Tm_eit:      ptr_new()                              ,$
                 index_eit:   ptr_new()                              ,$
                 Ne_mk4:      ptr_new()                              ,$
                 index_mk4:   ptr_new()                              ,$
                 Ne_kcor:     ptr_new()                              ,$
                 index_kcor:  ptr_new()                              ,$
                 Ne_c2:       ptr_new()                              ,$
                 index_c2:    ptr_new()                              ,$
                 index_sampling_aia:   ptr_new()                     ,$
                 index_sampling_euvia: ptr_new()                     ,$
                 index_sampling_euvib: ptr_new()                     ,$
                 index_sampling_eit:   ptr_new()                     ,$
                 index_sampling_mk4:   ptr_new()                     ,$                        
                 index_sampling_kcor:  ptr_new()                     ,$
                 index_sampling_c2:    ptr_new()                      }

; Store traced data
  if keyword_set(aia) then begin
     trace_data.Ne_aia     = ptr_new(    Ne_aia_A)
     trace_data.Tm_aia     = ptr_new(    Tm_aia_A)
     trace_data.index_aia  = ptr_new( index_aia_A)
     trace_data.index_sampling_aia = ptr_new(index_sampling_aia_A)
  endif
  if keyword_set(euvia) then begin
     trace_data.Ne_euvia     = ptr_new(    Ne_euvia_A)
     trace_data.Tm_euvia     = ptr_new(    Tm_euvia_A)
     trace_data.index_euvia  = ptr_new( index_euvia_A)
     trace_data.index_sampling_euvia = ptr_new(index_sampling_aia_A)
  endif
  if keyword_set(euvib) then begin
     trace_data.Ne_euvib     = ptr_new(    Ne_euvib_A)
     trace_data.Tm_euvib     = ptr_new(    Tm_euvib_A)
     trace_data.index_euvib  = ptr_new( index_euvib_A)
     trace_data.index_sampling_euvib = ptr_new(index_sampling_euvib_A)
  endif
  if keyword_set(eit) then begin
     trace_data.Ne_eit     = ptr_new(    Ne_eit_A)
     trace_data.Tm_eit     = ptr_new(    Tm_eit_A)
     trace_data.index_eit  = ptr_new( index_eit_A)
     trace_data.index_sampling_eit = ptr_new(index_sampling_eit_A)
  endif
  if keyword_set(mk4) then begin
     trace_data.Ne_mk4     = ptr_new(    Ne_mk4_A)
     trace_data.index_mk4  = ptr_new( index_mk4_A)
     trace_data.index_sampling_mk4 = ptr_new(index_sampling_mk4_A)
  endif
  if keyword_set(kcor) then begin
     trace_data.Ne_kcor     = ptr_new(    Ne_kcor_A)
     trace_data.index_kcor  = ptr_new( index_kcor_A)
     trace_data.index_sampling_kcor = ptr_new(index_sampling_kcor_A)
  endif
  if keyword_set(lascoc2) then begin
     trace_data.Ne_c2     = ptr_new(    Ne_c2_A)
     trace_data.index_c2  = ptr_new( index_c2_A)
     trace_data.index_sampling_c2 = ptr_new(index_sampling_c2_A)
  endif

; Perform fits if requested:
  if keyword_set(fits) then $
  fit_trace_data, aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                  mk4 = mk4, kcor = kcor, lascoc2 = lascoc2
  
 ; Save structure in dir_fl:
  save, trace_data, filename = dir_fl + structure_filename + '.sav'

  return
end
