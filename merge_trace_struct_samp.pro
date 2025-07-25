;
; PURPOSE: This code merges the tracing of tomographic products based
; on data from different instruments along given AWSoM/FDIPS-PFSS fieldlines, and
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
;
; HISTORY: V1.0 FAN, IAFE, July 2025.

pro merge_trace_struct_samp, fl_dir=fl_dir, fl_list=fl_list, $
                             aia=aia, euvia=euvia, euvib=euvib, eit=eit, $
                             mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2, $
                             struture_filename=structure_filename,$
                             trace_Bs=trace_Bs

  if not keyword_set(fl_dir) or not keyword_set(fl_list) then STOP

; Set up filename for output structure:
  if not keyword_set(structure_filename) then structure_filename = fl_list
  structure_filename = structure_filename + '-tracing-structure-merge-samp'

; Read the list with the field-lines  
  N_fl     = 0L
  filename = ''
  openr,1,fl_dir+fl_list
  readf,1,N_fl

; Maximum number of point along the fieldline
; FEDE: Tal vez podemos optimizar el valor de esta
; variable. Crucial para optimizar el uso de RAM.  
  Npt_max_samp = 150
  
; Default value in all arrays.
  default = -678.
  
; Define needed arrays:
; 1st Line Geometry and Bfield:

  Footpoint_Rad_A = fltarr(N_fl)
  Footpoint_Lon_A = fltarr(N_fl)
  Footpoint_Lat_A = fltarr(N_fl)
  Apex_Rad_A      = fltarr(N_fl)
  Apex_Lon_A      = fltarr(N_fl)
  Apex_Lat_A      = fltarr(N_fl)

; 2nd Tomography Products: 
  if keyword_set(aia) then begin
     Npt_v_aia            = intarr(N_fl)
     Ne_aia_A             = fltarr(N_fl,Npt_max_samp) + default
     Tm_aia_A             = fltarr(N_fl,Npt_max_samp) + default
     WT_aia_A             = fltarr(N_fl,Npt_max_samp) + default
     ldem_flag_aia_A      = fltarr(N_fl,Npt_max_samp) + default
     rad_aia_A            = fltarr(N_fl,Npt_max_samp) + default
     lat_aia_A            = fltarr(N_fl,Npt_max_samp) + default
     lon_aia_A            = fltarr(N_fl,Npt_max_samp) + default
     if keyword_set(trace_Bs) then begin
        s_aia_A            = fltarr(N_fl,Npt_max_samp) + default
        Br_aia_A           = fltarr(N_fl,Npt_max_samp) + default
        Bth_aia_A          = fltarr(N_fl,Npt_max_samp) + default
        Bph_aia_A          = fltarr(N_fl,Npt_max_samp) + default
        B_aia_A            = fltarr(N_fl,Npt_max_samp) + default
     endif
  endif
  
  if keyword_set(euvia) then begin
     Npt_v_euvia            = intarr(N_fl)
     Ne_euvia_A             = fltarr(N_fl,Npt_max_samp) + default
     Tm_euvia_A             = fltarr(N_fl,Npt_max_samp) + default
     WT_euvia_A             = fltarr(N_fl,Npt_max_samp) + default
     ldem_flag_euvia_A      = fltarr(N_fl,Npt_max_samp) + default
      rad_euvia_A           = fltarr(N_fl,Npt_max_samp) + default
      lat_euvia_A           = fltarr(N_fl,Npt_max_samp) + default
      lon_euvia_A           = fltarr(N_fl,Npt_max_samp) + default
      if keyword_set(trace_Bs) then begin
         s_euvia_A           = fltarr(N_fl,Npt_max_samp) + default
         Br_euvia_A          = fltarr(N_fl,Npt_max_samp) + default
         Bth_euvia_A         = fltarr(N_fl,Npt_max_samp) + default
         Bph_euvia_A         = fltarr(N_fl,Npt_max_samp) + default
         B_euvia_A           = fltarr(N_fl,Npt_max_samp) + default
      endif
   endif
  
  if keyword_set(euvib) then begin
     Npt_v_euvib            = intarr(N_fl)
     Ne_euvib_A             = fltarr(N_fl,Npt_max_samp) + default
     Tm_euvib_A             = fltarr(N_fl,Npt_max_samp) + default
     WT_euvib_A             = fltarr(N_fl,Npt_max_samp) + default
     ldem_flag_euvib_A      = fltarr(N_fl,Npt_max_samp) + default
      rad_euvib_A           = fltarr(N_fl,Npt_max_samp) + default
      lat_euvib_A           = fltarr(N_fl,Npt_max_samp) + default
      lon_euvib_A           = fltarr(N_fl,Npt_max_samp) + default
      if keyword_set(trace_Bs) then begin
         s_euvib_A           = fltarr(N_fl,Npt_max_samp) + default
         Br_euvib_A          = fltarr(N_fl,Npt_max_samp) + default
         Bth_euvib_A         = fltarr(N_fl,Npt_max_samp) + default
         Bph_euvib_A         = fltarr(N_fl,Npt_max_samp) + default
         B_euvib_A           = fltarr(N_fl,Npt_max_samp) + default
      endif
   endif

  if keyword_set(eit) then begin
     Npt_v_eit            = intarr(N_fl)
     Ne_eit_A             = fltarr(N_fl,Npt_max_samp) + default
     Tm_eit_A             = fltarr(N_fl,Npt_max_samp) + default
     WT_eit_A             = fltarr(N_fl,Npt_max_samp) + default
     ldem_flag_eit_A      = fltarr(N_fl,Npt_max_samp) + default
      rad_eit_A           = fltarr(N_fl,Npt_max_samp) + default
      lat_eit_A           = fltarr(N_fl,Npt_max_samp) + default
      lon_eit_A           = fltarr(N_fl,Npt_max_samp) + default
      if keyword_set(trace_Bs) then begin
        s_eit_A           = fltarr(N_fl,Npt_max_samp) + default
        Br_eit_A          = fltarr(N_fl,Npt_max_samp) + default
        Bth_eit_A         = fltarr(N_fl,Npt_max_samp) + default
        Bph_eit_A         = fltarr(N_fl,Npt_max_samp) + default
        B_eit_A           = fltarr(N_fl,Npt_max_samp) + default
     endif
   endif
  
  if keyword_set(mk4) then begin
     Npt_v_mk4            = intarr(N_fl)
     Ne_mk4_A             = fltarr(N_fl,Npt_max_samp) + default
     rad_mk4_A            = fltarr(N_fl,Npt_max_samp) + default
     lat_mk4_A            = fltarr(N_fl,Npt_max_samp) + default
     lon_mk4_A            = fltarr(N_fl,Npt_max_samp) + default
     if keyword_set(trace_Bs) then begin
        s_mk4_A           = fltarr(N_fl,Npt_max_samp) + default
        Br_mk4_A          = fltarr(N_fl,Npt_max_samp) + default
        Bth_mk4_A         = fltarr(N_fl,Npt_max_samp) + default
        Bph_mk4_A         = fltarr(N_fl,Npt_max_samp) + default
        B_mk4_A           = fltarr(N_fl,Npt_max_samp) + default
     endif
  endif
  
  if keyword_set(kcor) then begin
     Npt_v_kcor            = intarr(N_fl)
     Ne_kcor_A             = fltarr(N_fl,Npt_max_samp) + default
     rad_kcor_A            = fltarr(N_fl,Npt_max_samp) + default
     lat_kcor_A            = fltarr(N_fl,Npt_max_samp) + default
     lon_kcor_A            = fltarr(N_fl,Npt_max_samp) + default
     if keyword_set(trace_Bs) then begin
        s_kcor_A           = fltarr(N_fl,Npt_max_samp) + default
        Br_kcor_A          = fltarr(N_fl,Npt_max_samp) + default
        Bth_kcor_A         = fltarr(N_fl,Npt_max_samp) + default
        Bph_kcor_A         = fltarr(N_fl,Npt_max_samp) + default
        B_kcor_A           = fltarr(N_fl,Npt_max_samp) + default
     endif
  endif
  if keyword_set(ucomp) then begin
     Npt_v_ucomp            = intarr(N_fl)
     Ne_ucomp_A             = fltarr(N_fl,Npt_max_samp) + default
     rad_ucomp_A            = fltarr(N_fl,Npt_max_samp) + default
     lat_ucomp_A            = fltarr(N_fl,Npt_max_samp) + default
     lon_ucomp_A            = fltarr(N_fl,Npt_max_samp) + default
     if keyword_set(trace_Bs) then begin
        s_ucomp_A           = fltarr(N_fl,Npt_max_samp) + default
        Br_ucomp_A          = fltarr(N_fl,Npt_max_samp) + default
        Bth_ucomp_A         = fltarr(N_fl,Npt_max_samp) + default
        Bph_ucomp_A         = fltarr(N_fl,Npt_max_samp) + default
        B_ucomp_A           = fltarr(N_fl,Npt_max_samp) + default
     endif
  endif
  
  if keyword_set(lascoc2) then begin
     Npt_v_c2            = intarr(N_fl)
     Ne_c2_A             = fltarr(N_fl,Npt_max_samp) + default
     rad_c2_A            = fltarr(N_fl,Npt_max_samp) + default
     lat_c2_A            = fltarr(N_fl,Npt_max_samp) + default
     lon_c2_A            = fltarr(N_fl,Npt_max_samp) + default
     if keyword_set(trace_Bs) then begin
        s_c2_A           = fltarr(N_fl,Npt_max_samp) + default
        Br_c2_A          = fltarr(N_fl,Npt_max_samp) + default
        Bth_c2_A         = fltarr(N_fl,Npt_max_samp) + default
        Bph_c2_A         = fltarr(N_fl,Npt_max_samp) + default
        B_c2_A           = fltarr(N_fl,Npt_max_samp) + default
     endif
  
  endif
           
  for i_fl = 0L,N_fl-1 do begin ; Start loop in fieldlines  
     initialized = 'no'
     readf,1,filename
     
   ; AIA
     if keyword_set(aia)   then begin
        file_aia   = filename+'_aia.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_aia'
        if NOT keyword_set(trace_Bs) then $
        readcol,fl_dir+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_aia_l,Tm_aia_l,$
                WT_aia_l, ldem_flag_aia_l, index_aia_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
        readcol,fl_dir+file_aia,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_aia_l,Tm_aia_l,$
                WT_aia_l, ldem_flag_aia_l, index_aia_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_aia_l, index_l=index_aia_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           index_foot = where( rad_l eq min(rad_l))
           index_apex = where( rad_l eq max(rad_l))
           footpoint_rad_A(i_fl) = rad_l(index_foot[0])
           footpoint_lat_A(i_fl)=  lat_l(index_foot[0])
           footpoint_lon_A(i_fl) = lon_l(index_foot[0])
           apex_rad_A(i_fl) = rad_l(index_apex[0])
           apex_lat_A(i_fl)=  lat_l(index_apex[0])
           apex_lon_A(i_fl) = lon_l(index_apex[0])
        endif
        indsamp = where(index_sampling_l eq 1)
        if indsamp[0] ne -1 then begin
           Nl_samp          = n_elements(indsamp)
           Npt_v_aia(i_fl)  = Nl_samp 
           Ne_aia_A            (i_fl,0:Nl_samp-1) = Ne_aia_l(indsamp)
           Tm_aia_A            (i_fl,0:Nl_samp-1) = Tm_aia_l(indsamp)
           WT_aia_A            (i_fl,0:Nl_samp-1) = WT_aia_l(indsamp)
           ldem_flag_aia_A     (i_fl,0:Nl_samp-1) = ldem_flag_aia_l(indsamp)
           rad_aia_A           (i_fl,0:Nl_samp-1) = rad_l(indsamp)
           lat_aia_A           (i_fl,0:Nl_samp-1) = lat_l(indsamp)
           lon_aia_A           (i_fl,0:Nl_samp-1) = lon_l(indsamp)
           if keyword_set(trace_Bs) then begin
               s_aia_A          (i_fl,0:Nl_samp-1) =   S_l (indsamp)
              Br_aia_A          (i_fl,0:Nl_samp-1) =  Br_l (indsamp)
             Bth_aia_A          (i_fl,0:Nl_samp-1) = Bth_l (indsamp)
             Bph_aia_A          (i_fl,0:Nl_samp-1) = Bph_l (indsamp)
               B_aia_A          (i_fl,0:Nl_samp-1) =   B_l (indsamp)
           endif
        endif
     endif
; Nota del coder: Hay que repetir lo que se codeó con AIA para los
; demás instrumentos ...
     
   ; EUVI-A
     if keyword_set(euvia)  then begin
        file_euvia = filename+'_euvia.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvia'
        if NOT keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,$
                WT_euvia_l, ldem_flag_euvia_l, index_euvia_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvia,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,$
                WT_euvia_l, ldem_flag_euvia_l, index_euvia_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvia_l, index_l=index_euvia_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_euvia_A            (i_fl,0:N_l-1) = Ne_euvia_l
        Tm_euvia_A            (i_fl,0:N_l-1) = Tm_euvia_l
        WT_euvia_A            (i_fl,0:N_l-1) = WT_euvia_l
        ldem_flag_euvia_A     (i_fl,0:N_l-1) = ldem_flag_euvia_l
        index_euvia_A         (i_fl,0:N_l-1) = index_euvia_l    
        index_sampling_euvia_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; EUVI-B
     if keyword_set(euvib)  then begin
        file_euvib = filename+'_euvib.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_euvib'
        if NOT keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvib,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
        readcol,fl_dir+file_euvib,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,$
                WT_euvib_l, ldem_flag_euvib_l, index_euvib_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_euvib_l, index_l=index_euvib_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_euvib_A            (i_fl,0:N_l-1) = Ne_euvib_l
        Tm_euvib_A            (i_fl,0:N_l-1) = Tm_euvib_l
        WT_euvib_A            (i_fl,0:N_l-1) = WT_euvib_l
        ldem_flag_euvib_A     (i_fl,0:N_l-1) = ldem_flag_euvib_l
        index_euvib_A         (i_fl,0:N_l-1) = index_euvib_l
        index_sampling_euvib_A(i_fl,0:N_l-1) = index_sampling_l
     endif
   ; EIT
     if keyword_set(eit)    then begin
        file_eit = filename+'_eit.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_eit'
        if NOT keyword_set(trace_Bs) then $
        readcol,fl_dir+file_eit,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,$
                WT_eit_l, ldem_flag_eit_l, index_eit_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
        readcol,fl_dir+file_eit,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,$
                WT_eit_l, ldem_flag_eit_l, index_eit_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_eit_l, index_l=index_eit_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_eit_A            (i_fl,0:N_l-1) = Ne_eit_l
        Tm_eit_A            (i_fl,0:N_l-1) = Tm_eit_l
        WT_eit_A            (i_fl,0:N_l-1) = WT_eit_l
        ldem_flag_eit_A     (i_fl,0:N_l-1) = ldem_flag_eit_l
        index_eit_A         (i_fl,0:N_l-1) = index_eit_l
        index_sampling_eit_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; Mk4
     if keyword_set(mk4)    then begin
        file_mk4   = filename+'_mk4.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_mk4'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_mk4,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_mk4_l,$
                   index_mk4_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_mk4,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_mk4_l,$
                   index_mk4_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_mk4_l, index_l=index_mk4_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_mk4_A   (i_fl,0:N_l-1) = Ne_mk4_l
        index_mk4_A(i_fl,0:N_l-1) = index_mk4_l
        index_sampling_mk4_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; KCOR
     if keyword_set(kcor)    then begin
        file_kcor   = filename+'_kcor.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_kcor'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_kcor,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_kcor,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_kcor_l,$
                   index_kcor_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_kcor_l, index_l=index_kcor_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_kcor_A   (i_fl,0:N_l-1) = Ne_kcor_l
        index_kcor_A(i_fl,0:N_l-1) = index_kcor_l        
        index_sampling_kcor_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; UCOMP
       if keyword_set(ucomp)    then begin
        file_ucomp   = filename+'_ucomp.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_ucomp'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_ucomp,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_ucomp_l,$
                   index_ucomp_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_ucomp,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_ucomp_l,$
                   index_ucomp_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_ucomp_l, index_l=index_ucomp_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_ucomp_A   (i_fl,0:N_l-1) = Ne_ucomp_l
        index_ucomp_A(i_fl,0:N_l-1) = index_ucomp_l        
        index_sampling_ucomp_A(i_fl,0:N_l-1) = index_sampling_l
     endif

   ; LASCO-C2
     if keyword_set(lascoc2) then begin
        file_c2    = filename+'_lascoc2.out'
        if i_fl eq 0 then structure_filename = structure_filename + '_lascoc2'
        if NOT keyword_set(trace_Bs) then $
           readcol,fl_dir+file_c2,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_c2_l,$
                   index_c2_l  ,FORMAT='D,D,D,D,D,D'
        if     keyword_set(trace_Bs) then $
           readcol,fl_dir+file_c2,x_l,y_l,z_l,s_l,Br_l,Bth_l,Bph_l,B_l,rad_l,lat_l,lon_l,Ne_c2_l,$
                   index_c2_l  ,FORMAT='D,D,D,D,D,D,D,D,D,D,D'
        sample_fl, Ne_l=Ne_c2_l, index_l=index_c2_l, index_sampling_l=index_sampling_l
        if initialized eq 'no' then begin
           N_l         = n_elements(x_l)
           Npt_v(i_fl) = N_l  
           initialized = 'yes'
           x_A   (i_fl,0:N_l-1) = x_l
           y_A   (i_fl,0:N_l-1) = y_l
           z_A   (i_fl,0:N_l-1) = z_l
           rad_A (i_fl,0:N_l-1) = rad_l
           lat_A (i_fl,0:N_l-1) = lat_l
           lon_A (i_fl,0:N_l-1) = lon_l  
           if keyword_set(trace_Bs) then begin
             s_A (i_fl,0:N_l-1) =   s_l
            Br_A (i_fl,0:N_l-1) =  Br_l
           Bth_A (i_fl,0:N_l-1) = Bth_l
           Bph_A (i_fl,0:N_l-1) = Bph_l
             B_A (i_fl,0:N_l-1) =   B_l
           endif
        endif
        Ne_c2_A   (i_fl,0:N_l-1) = Ne_c2_l
        index_c2_A(i_fl,0:N_l-1) = index_c2_l    
        index_sampling_c2_A(i_fl,0:N_l-1) = index_sampling_l
     endif
                              
  endfor                        ; End loop in fieldlines    
  close,1
  STOP
; POINTER-STRUCTURE  
; Create a pointer structure to store field line extraction information  
  trace_data = { N_fl:                ptr_new(N_fl)                ,$
                 Npt_max_samp:        ptr_new(Npt_max_samp)        ,$
                 Npt_v_samp:          ptr_new(Npt_v)               ,$
                 footpoint_rad:       ptr_new(footpoint_rad_A)     ,$
                 footpoint_lat:       ptr_new(footpoint_lat_A)     ,$
                 footpoint_lon:       ptr_new(footpoint_lon_A)     ,$
                 apex_rad:            ptr_new(apex_rad_A)          ,$
                 apex_lat:            ptr_new(apex_lat_A)          ,$
                 apex_lon:            ptr_new(apex_lon_A)           }
  undefine,Npt_max_samp
  undefine,Npt_v_samp
  undefine,footpoint_rad_A
  undefine,footpoint_lat_A
  undefine,footpoint_lon_A
  undefine,apex_rad_A
  undefine,apex_lat_A
  undefine,apex_lon_A


if keyword_set(trace_Bs) then begin
; Read-in the leg-label of all fieldlines.
  leg_label_A = lonarr(N_fl)
  leglab=0L
  openr,1,fl_dir+'legs-label.dat'
  readf,1,N_fl
  for i_fl=0L,N_fl-1 do begin
     readf,1,leglab
     leg_label_A(i_fl)=leglab
  endfor
  close,1
; Read-in the leg-length of all fieldlines.
  leg_length_A = fltarr(N_fl)
  leglength=0.
  openr,1,fl_dir+'legs-length.dat'
  readf,1,N_fl
  for i_fl=0L,N_fl-1 do begin
     readf,1,leglength
     leg_length_A(i_fl)=leglength
  endfor
  close,1
; Read-in the leg footpoint Bfield of all fieldlines.
  leg_footbfield_A = fltarr(N_fl,3)
  xx='' & Br=0. & Bth=0. & Bph=0.
  openr,1,fl_dir+'legs-footpoint-Bfield.dat'
  readf,1,N_fl
  readf,1,xx
  for i_fl=0L,N_fl-1 do begin
     readf,1,Br,Bth,Bph
     leg_footbfield_A(i_fl,*)=[Br,Bth,Bph]
  endfor
  close,1
; Add arrays to structure
  trace_data = create_struct( trace_data ,$
               'leg_label'      , ptr_new(leg_label_A)     ,$
               'leg_length'     , ptr_new(leg_length_A)    ,$
               'leg_footbfield' , ptr_new(leg_footbfield_A) )
  undefine,leg_label_A
  undefine,leg_length_A
  undefine,leg_footbfield_A
endif

; Store into the structure traced tomographic resuls
  if keyword_set(aia) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_aia' , ptr_new(            Ne_aia_A) ,$
                   'Tm_aia' , ptr_new(            Tm_aia_A) ,$
                   'WT_aia' , ptr_new(            WT_aia_A) ,$
            'ldem_flag_aia' , ptr_new(     ldem_flag_aia_A) ,$
                  'rad_aia' , ptr_new(           rad_aia_A) ,$
                  'lat_aia' , ptr_new(           lat_aia_A) ,$
                  'lon_aia' , ptr_new(           lon_aia_A) )
     undefine,Ne_aia_A
     undefine,Tm_aia_A
     undefine,WT_aia_A
     undefine,ldem_flag_aia_A
     undefine,rad_aia_A
     undefine,lat_aia_A
     undefine,lon_aia_A
     if keyword_set(trace_Bs) then begin
        trace_data = create_struct( trace_data ,$
                                    's_aia', ptr_new(s_aia_A), $
                                    'B_aia', ptr_new(B_aia_A), $
                                    'Br_aia', ptr_new(Br_aia_A), $
                                    'Bth_aia', ptr_new(Bth_aia_A), $
                                    'Bph_aia', ptr_new(Bph_aia_A) )
        undefine,s_aia_A
        undefine,B_aia_A
        undefine,Br_aia_A
        undefine,Bth_aia_A
        undefine,Bph_aia_A
     endif
  endif
; Nota del coder: Hay que repetir lo que se codeó con AIA para los
; demás instrumentos ...
  
  if keyword_set(euvia) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_euvia' , ptr_new(            Ne_euvia_A) ,$
                 'Tm_euvia' , ptr_new(            Tm_euvia_A) ,$
                 'WT_euvia' , ptr_new(            WT_euvia_A) ,$
          'ldem_flag_euvia' , ptr_new(     ldem_flag_euvia_A) ,$
              'index_euvia' , ptr_new(         index_euvia_A) ,$
     'index_sampling_euvia' , ptr_new(index_sampling_euvia_A)  )
     undefine,Ne_euvia_A
     undefine,Tm_euvia_A
     undefine,WT_euvia_A
     undefine,ldem_flag_euvia_A
     undefine,index_euvia_A
     undefine,index_sampling_euvia_A
  endif
  if keyword_set(euvib) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_euvib' , ptr_new(            Ne_euvib_A) ,$
                 'Tm_euvib' , ptr_new(            Tm_euvib_A) ,$
                 'WT_euvib' , ptr_new(            WT_euvib_A) ,$
          'ldem_flag_euvib' , ptr_new(     ldem_flag_euvib_A) ,$
              'index_euvib' , ptr_new(         index_euvib_A) ,$
     'index_sampling_euvib' , ptr_new(index_sampling_euvib_A)  )
     undefine,Ne_euvib_A
     undefine,Tm_euvib_A
     undefine,WT_euvib_A
     undefine,ldem_flag_euvib_A
     undefine,index_euvib_A
     undefine,index_sampling_euvib_A
  endif
  if keyword_set(eit) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_eit' , ptr_new(            Ne_eit_A) ,$
                   'Tm_eit' , ptr_new(            Tm_eit_A) ,$
                   'WT_eit' , ptr_new(            WT_eit_A) ,$
            'ldem_flag_eit' , ptr_new(     ldem_flag_eit_A) ,$
                'index_eit' , ptr_new(         index_eit_A) ,$
       'index_sampling_eit' , ptr_new(index_sampling_eit_A)  )
     undefine,Ne_eit_A
     undefine,Tm_eit_A
     undefine,WT_eit_A
     undefine,ldem_flag_eit_A
     undefine,index_eit_A
     undefine,index_sampling_eit_A
  endif
  if keyword_set(mk4) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_mk4' , ptr_new(            Ne_mk4_A) ,$
                'index_mk4' , ptr_new(         index_mk4_A) ,$
       'index_sampling_mk4' , ptr_new(index_sampling_mk4_A)  )
     undefine,Ne_mk4_A
     undefine,index_mk4_A
     undefine,index_sampling_mk4_A
  endif
  if keyword_set(kcor) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_kcor' , ptr_new(            Ne_kcor_A) ,$
              'index_kcor' , ptr_new(         index_kcor_A) ,$
     'index_sampling_kcor' , ptr_new(index_sampling_kcor_A)  )
     undefine,Ne_kcor_A
     undefine,index_kcor_A
     undefine,index_sampling_kcor_A
  endif
  if keyword_set(ucomp) then begin
     trace_data = create_struct( trace_data ,$
                 'Ne_ucomp' , ptr_new(            Ne_ucomp_A) ,$
              'index_ucomp' , ptr_new(         index_ucomp_A) ,$
     'index_sampling_ucomp' , ptr_new(index_sampling_ucomp_A)  )
     undefine,Ne_ucomp_A
     undefine,index_ucomp_A
     undefine,index_sampling_ucomp_A
  endif
  if keyword_set(lascoc2) then begin
     trace_data = create_struct( trace_data ,$
                   'Ne_c2' , ptr_new(            Ne_c2_A) ,$
                'index_c2' , ptr_new(         index_c2_A) ,$
       'index_sampling_c2' , ptr_new(index_sampling_c2_A)  )
     undefine,Ne_c2_A
     undefine,index_c2_A
     undefine,index_sampling_c2_A
  endif

; Perform fits:
  fit_trace_data, aia=aia, euvia=euvia, euvib=euvib, eit=eit,$
                  mk4=mk4, kcor=kcor, ucomp=ucomp, lascoc2=lascoc2,$
                  fl_dir=fl_dir, trace_data = trace_data 
  
 ; Save structure in fl_dir:
  save, trace_data, filename = fl_dir + structure_filename + '.sav'

  print,'Output in:'
  print,fl_dir + structure_filename + '.sav'
  
  return
end
