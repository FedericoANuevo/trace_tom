pro expand_csv_for_sam_and_yeimy
  
    common data, N_fl, Npt_max, Npt_v, x_A, y_A, z_A, rad_A, lat_A, lon_A,$
       Ne_aia_A, Tm_aia_A, WT_aia_A, ldem_flag_aia_A, index_aia_A, index_sampling_aia_A,$
       Ne_euvia_A, Tm_euvia_A,  WT_euvia_A, ldem_flag_euvia_A, index_euvia_A, index_sampling_euvia_A,$
       Ne_euvib_A, Tm_euvib_A, WT_euvib_A, ldem_flag_euvib_A, index_euvib_A, index_sampling_euvib_A,$
       Ne_eit_A, Tm_eit_A, WT_eit_A, ldem_flag_eit_A, index_eit_A, index_sampling_eit_A,$
       Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
       Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
       Ne_c2_A, index_c2_A, index_sampling_c2_A,$
       rad_fit_aia_A, Ne_fit_aia_A, Tm_fit_aia_A, fitflag_aia_A,scN_fit_aia_A,scT_fit_aia_A,$
       rad_fit_euvia_A, Ne_fit_euvia_A, Tm_fit_euvia_A, fitflag_euvia_A,scN_fit_euvia_A,scT_fit_euvia_A,$
       rad_fit_euvib_A, Ne_fit_euvib_A, Tm_fit_euvib_A, fitflag_euvib_A,scN_fit_euvib_A,scT_fit_euvib_A,$
       rad_fit_eit_A, Ne_fit_eit_A, Tm_fit_eit_A, fitflag_eit_A,scN_fit_eit_A,scT_fit_eit_A,$
       rad_fit_c2_A, Ne_fit_c2_A, fitflag_c2_A,scN_fit_c2_A,$
       rad_fit_mk4_A, Ne_fit_mk4_A, fitflag_mk4_A,scN_fit_mk4_A,$
       N0_fit_aia_A,lN_fit_aia_A,T0_fit_aia_A,dTdr_fit_aia_A,$
       N0_fit_euvia_A,lN_fit_euvia_A,T0_fit_euvia_A,dTdr_fit_euvia_A,$
       N0_fit_euvib_A,lN_fit_euvib_A,T0_fit_euvib_A,dTdr_fit_euvib_A,$
       N0_fit_eit_A,lN_fit_eit_A,T0_fit_eit_A,dTdr_fit_eit_A,$
       N1_fit_aia_A,N2_fit_aia_A,p1_fit_aia_A,p2_fit_aia_A,$
       N1_fit_euvia_A,N2_fit_euvia_A,p1_fit_euvia_A,p2_fit_euvia_A,$
       N1_fit_euvib_A,N2_fit_euvib_A,p1_fit_euvib_A,p2_fit_euvib_A,$
       N1_fit_eit_A,N2_fit_eit_A,p1_fit_eit_A,p2_fit_eit_A,$
       N0_fit_mk4_A,lN_fit_mk4_A,$
       N1_fit_mk4_A,N2_fit_mk4_A,p1_fit_mk4_A,p2_fit_mk4_A,$
       N1_fit_c2_A,N2_fit_c2_A,p1_fit_c2_A,p2_fit_c2_A,$
       lN_fit_c2_A,$
       fit_F_Ne_aia,fit_F_Ne_mk4,fit_F_Ne_c2,$
       fit_F_Ne_euvia,fit_F_Ne_euvib,fit_F_eit_c2

; Directory and File-name 
  dir  = '/media/Data1/data1/DATA/flines_Sam-Yeimy/'
  file = 'April2024_PFSS-tracing-structure-merge_aia_lascoc2.sav'  
; Load the traced-data-structure (stored in common data )
  load_traced_data_structure,dir=dir,structure_filename=file,/aia,/lascoc2

  for ifl = 0,N_fl-1 do begin

     id_l = indgen(Npt_v(ifl))
     x_l = reform(x_A(ifl,0:Npt_v(ifl)-1))
     y_l = reform(y_A(ifl,0:Npt_v(ifl)-1))
     z_l = reform(z_A(ifl,0:Npt_v(ifl)-1))
     tmp  = reform(index_sampling_aia_A(ifl,0:Npt_v(ifl)-1))
     Ne_aia_l = reform(Ne_aia_A(ifl,0:Npt_v(ifl)-1))
     Tm_aia_l = reform(Tm_aia_A(ifl,0:Npt_v(ifl)-1))
     index    = where(tmp le 0)
     Ne_aia_l(index) = -1.
     Tm_aia_l(index) = -1.
     tmp     = reform(index_sampling_c2_A(ifl,0:Npt_v(ifl)-1))
     Ne_c2_l = reform(Ne_c2_A(ifl,0:Npt_v(ifl)-1))
     index    = where(tmp le 0)
     Ne_c2_l(index) = -1
     if ifl lt   10                 then index_str = '000' + strmid(string(ifl),7,1)
     if ifl ge   10 and ifl lt  100 then index_str = '00'  + strmid(string(ifl),6,2)
     if ifl ge  100 and ifl lt 1000 then index_str =  '0'  + strmid(string(ifl),5,3)
     if ifl ge 1000                 then index_str =         strmid(string(ifl),4,4)
     
     output_file = dir+'fline_'+index_str+'_tom.csv'
     table_hdr = ['x [Rsun]','y [Rsun]','z [Rsun]','Ne AIA [cm-3]','Te AIA [K]','Ne C2 [cm-3]']
;    La Keyword table_header solo se puede usar con IDL 8.0     
     write_csv,output_file,id_l,x_l,y_l,z_l,Ne_aia_l,Tm_aia_l,Ne_c2_l;,table_header = table_hdr
;    STOP
  endfor

  
  return
end
  



  
