;
; PURPOSE: This code merges the tracing of tomographic products based
; on data from different instruments along the same fieldline.
;
; INPUTS:
; fir_fl and fl_list: STRINGS. directory where field lines are
; located, and filename of list of names of field lines filenames.
;
; FLAGS: one per allowed instrument, use those for which you want to
; merge results, in any order.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
;

pro merge_trace, dir_fl = dir_fl, fl_list = fl_list, $
                 aia = aia, euvia = euvia, euvib = euvib, eit = eit, $
                 mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

  if not keyword_set(fl_list) or not keyword_set(fl_list) then STOP
  
  ; Read the list with the field-lines  
  N_fl     = 0
  filename = ''
  openr,1,dir_fl+fl_list
  readf,1,N_fl
  for i_fl = 0,N_fl-1 do begin
     initialized = 'no'
     readf,1,filename
     outfile       = filename + '_merge'

   ; AIA
     if keyword_set(aia)   then begin
        file_aia   = filename+'_aia.out'
        outfile    = outfile +'_aia'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_aia_l  ,Tm_aia_l  ,index_aia_l  ,FORMAT='D,D,D,D,D,D'
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_aia_l],[Tm_aia_l],[index_aia_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [AIA, cm^-3]   Tm [AIA, K]        AIA-3Dind'
           output_format = '(8(E18.10)," ",I9'
        endif
     endif

   ; EUVI-A
     if keyword_set(euvia)  then begin
        file_euvia = filename+'_euvia.out'
        outfile   =  outfile +'_euvia'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvia_l,Tm_euvia_l,index_euvia_l,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_euvia_l],[Tm_euvia_l],[index_euvia_l]]
           header_str = header_str +'         Ne [EUVIA, cm^-3]  Tm [EUVIA, K]     EUVIA-3Dind'
           output_format = output_format + ',"  ",2(E18.10),"  ",I9' 
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_euvia_l],[Tm_euvia_l],[index_euvia_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EUVIA, cm^-3] Tm [EUVIA, K]    EUVIA-3Dind'
           output_format = '(8(E18.10)," ",I9'
        endif
     endif

   ; EUVI-B
     if keyword_set(euvib)  then begin
        file_euvib = filename+'_euvib.out'
        outfile    = outfile +'_euvib'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_euvib_l,Tm_euvib_l,index_euvib_l,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_euvib_l],[Tm_euvib_l],[index_euvib_l]]
           header_str = header_str +'         Ne [EUVIB, cm^-3]  Tm [EUVIB, K]     EUVIB-3Dind'
           output_format = output_format + ',"  ",2(E18.10),"  ",I9'
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_euvib_l],[Tm_euvib_l],[index_euvib_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EUVIB, cm^-3] Tm [EUVIB, K]    EUVIB-3Dind'
           output_format = '(8(E18.10)," ",I9'
        endif        
     endif

   ; EIT
     if keyword_set(eit)    then begin
        file_euvib = filename+'_eit.out'
        outfile    = outfile +'_eit'
        readcol,dir_fl+file_aia,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_eit_l,Tm_eit_l,index_eit_l,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_eit_l],[Tm_eit_l],[index_eit_l]]
           header_str = header_str +'         Ne [EIT, cm^-3]  Tm [EIT, K]     EIT-3Dind'
           output_format = output_format + ',"  ",2(E18.10),"  ",I9'
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_eit_l],[Tm_eit_l],[index_eit_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [EIT, cm^-3]   Tm [EIT, K]        EIT-3Dind'
           output_format = '(8(E18.10)," ",I9'
        endif        
     endif

   ; Mk4
     if keyword_set(mk4)    then begin
        file_mk4   = filename+'_mk4.out'
        outfile    = outfile +'_mk4'
        readcol,dir_fl+file_mk4,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_mk4_l             ,index_mk4_l  ,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_mk4_l],[index_mk4_l]] 
           header_str = header_str +'  Ne [Mk4, cm^-3]     Mk4-3Dind'
           output_format = output_format + ',"  ",E18.10,"  ",I9'
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_mk4_l],[index_mk4_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [Mk4, cm^-3]     Mk4-3Dind'
           output_format = '(7(E18.10)," ",I9'
        endif        
     endif

   ; KCOR
     if keyword_set(kcor)    then begin
        file_kcor   = filename+'_kcor.out'
        outfile     = outfile +'_kcor'
        readcol,dir_fl+file_kcor,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_kcor_l           ,index_kcor_l ,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_kcor_l],[index_kcor_l]] 
           header_str = header_str +'  Ne [KCOR, cm^-3]     KCOR-3Dind'
           output_format = output_format + ',"  ",E18.10,"  ",I9'
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_kcor_l],[index_kcor_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [KCOR, cm^-3]     KCOR-3Dind'
           output_format = '(7(E18.10)," ",I9'
        endif
     endif

   ; LASCO-C2
     if keyword_set(lascoc2) then begin
        file_c2    = filename+'_lascoc2.out'
        outfile    = outfile +'_lascoc2'
        readcol,dir_fl+file_c2 ,x_l,y_l,z_l,rad_l,lat_l,lon_l,Ne_c2_l              ,index_c2_l   ,FORMAT='D,D,D,D,D,D'
        if initialized eq 'yes' then begin
           output_columns = [[[output_columns]],[Ne_c2_l],[index_c2_l]] 
           header_str = header_str +'  Ne [C2, cm^-3]     C2-3Dind'
           output_format = output_format + ',"  ",E18.10,"  ",I9'
        endif
        if initialized eq 'no' then begin
           initialized = 'yes'
           output_columns = [[x_l],[y_l],[z_l],[rad_l],[lat_l],[lon_l],[Ne_c2_l],[index_c2_l]]
           header_str = '  X [Rs]            Y [Rs]            Z [Rs]            RAD [Rs]          LAT [deg]         LON [deg]         Ne [C2, cm^-3]     C2-3Dind'
           output_format = '(7(E18.10)," ",I9'
        endif
     endif

   ; Close output filename
     outfile       = outfile + '.out'
   ; Close output format string
     output_format = output_format + ')'
   ; Transpose output_columns array
     output_columns = transpose(output_columns)
   ; Write output file
     openw,2,dir_fl+outfile
     printf,2,header_str
     N = n_elements(x_l)
     for i = 0,N-1 do printf,2,output_columns(*,i),FORMAT = output_format
     close,2
  endfor
     close,1
  return
end
