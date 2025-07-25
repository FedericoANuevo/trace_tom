;
; PURPOSE: This code is the main routine, where the AWSoM fieldline
; files are provided, as well as all the needed parameters and files
; to read the tomographic producs to be traced.
;
; CALLING SEQUENCE EXAMPLES:
; tom_trace_call, /demt , /trace_Bs
; tom_trace_call, /ucomp, /trace_Bs
; tom_trace_call, /kcor , /trace_Bs
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;          V1.1 AMV, CLaSP, May-2025, centralized and clearer dir structure.
;                    Also added /nfs flag.
;          V1.2 AMV&FAN, CLaSP, May-2025, added trace_Bs flag (also in
;                    called sub-routines´
;

pro tom_trace_call,demt=demt,lasco=lasco,kcor_mk4=kcor_mk4,ucomp=ucomp,$
                   nfs1=nfs1,nfs2=nfs2,trace_Bs=trace_Bs

;===============================================================================================
; Define PROJECT_NAME, a string suffix to construct the full PATHS to the required files.
  PROJECT_NAME = 'CR2254'
; PROJECT_NAME = 'CR2261'

; Define field_line_geometry_suffix_dir
  field_line_geometry_suffix_dir='_aunifgrid_multirad_2x2deg_HMI-PolFil/'
; field_line_geometry_suffix_dir='_aunifgrid_2.50Rs_2x2deg_HMI-PolFil/'
  
; Provide FL_LIST, the file which informs the number of field lines and the
; filenames of the ASCII files containing the geometry of each line.
;
; fl_list = 'fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt'
; fl_list = 'fdips_field_150X180X360_mrmqs220831t1302c2261_000.ubdat_fline-filenames_list.txt'
; fl_list = 'fdips_field_150X180X360_hmi.Synoptic_Mr.2254.ubdat_fline-filenames_list.txt'
  fl_list = 'fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2254.ubdat_fline-filenames_list.txt'
; fl_list = 'fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2261.ubdat_fline-filenames_list.txt'
;===============================================================================================

; --------------------This block should not require edits.---------------------------
; Set  FL_DIR, where the field-lines geometry files should be located,
; and TOM_DIR, where the 3D tomography products to trace should be located.
  base_dir = '/data1/'
  if keyword_set(nfs1) then base_dir = '/data/Data1/data1/'
  if keyword_set(nfs2) then base_dir = '/data/Data2/data1/'
  if not keyword_set(field_line_geometry_suffix_dir) then field_line_geometry_suffix_dir='/'
   fl_dir = base_dir+'DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry'+field_line_geometry_suffix_dir
  tom_dir = base_dir+'DATA/trace_tom_files/'+PROJECT_NAME+'/tomography_3Dproducts/' 
;------------------------------------------------------------------------------------
    
; Set TOM_FILE, the filename (or filename suffix in case of DEMT)
; containing the tomographic 3D products, the instrument label, 
; and all parameters of the tomography computational ball.
; Then call tom_trace.
   
  if keyword_set(demt) then begin
    ;tom_file = 'CR2082_HOLLOW_compound2.dat' & instr = 'euvia'
    ;tom_file = 'LDEM.CR2099_aia_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound1.dat'    & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound2.dat'    & instr = 'aia'
    ;tom_file = 'LDEM.April-2024_aia_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & instr = 'aia'
     if PROJECT_NAME eq 'CR2254' then $
        tom_file = 'LDEM.feb-mar_2022_segment2_aia_Hollow_3Bands_ucomp_comparison_gauss1_lin_Norm-median_singlStart'
     if PROJECT_NAME eq 'CR2261' then $
        tom_file = 'LDEM.chip07_aia_Hollow_3Bands_ucomp_comparison_exp2_gauss1_lin_Norm-median_singlStart'
     instr    = 'aia'
     nr       = 21   ;30
     nt       = 60   ;90
     np       = 2*nt
     rmin     = 1.09 ;1.0
     rmax     = 1.30 ;1.3
     Irmin    = 1.09 ;1.02
     Irmax    = 1.25 ;1.25
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,trace_Bs=trace_Bs
  endif

  if keyword_set(lasco) then begin
    ;tom_file = 'x_LascoC2pB_CR2099_shifted_std-grid_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L6.e-6'
    ;tom_file = 'x_LascoC2pB_CR2082_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L8.2e-6'
    ;tom_file = 'x_LascoC2pB_April-2024_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L1.1e-5'
     instr    = 'lascoc2'
     nr       = 60
     nt       = 60
     np       = 2*nt
     rmin     = 2.5
     rmax     = 8.5
     Irmin    = 2.5
     Irmax    = 6.0
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,/csv,trace_Bs=trace_Bs
  endif

  if keyword_set(kcor_mk4) then begin
    ;tom_file = 'x_Mk4_CR2099_shifted_Rmin1.15_Rmax1.85_IRmin1.15_IRmax1.50_70x90x180_BF2_L5.e-6' & instr = 'mk4'
     if PROJECT_NAME eq 'CR2254' then $
        tom_file = 'x_kcor_2022_Feb-Mar_segment2_Rmin1.09_Rmax1.75_IRmin1.09_IRmax1.50_66x60x120_BF2_L2.6e-5_CORRECT-FACT-0.79'
     if PROJECT_NAME eq 'CR2261' then $
        tom_file = 'x_kcor_20220823-20220905_Rmin1.09_Rmax1.75_IRmin1.09_IRmax1.50_66x60x120_BF2_L2.5e-5_CORRECT-FACT-0.79'
     instr    = 'kcor'
     nr       = 66   ;70
     nt       = 60   ;90
     np       = 2*nt
     rmin     = 1.09 ;1.15
     rmax     = 1.75 ;1.85
     Irmin    = 1.09 ;1.15
     Irmax    = 1.50 ;1.50
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,trace_Bs=trace_Bs
  endif

  if keyword_set(ucomp) then begin
     if PROJECT_NAME eq 'CR2254' then $
        tom_file = 'Ne_ratio_ucomp_1074-1079_second_target_segment2.dat'
     if PROJECT_NAME eq 'CR2261' then $
        tom_file = 'Ne_ratio_ucomp_1074-1079.dat'
     instr    = 'ucomp'
     nr       = 16
     nt       = 60
     np       = 2*nt
     rmin     = 1.09
     rmax     = 1.25
     Irmin    = 1.09
     Irmax    = 1.20
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,trace_Bs=trace_Bs
  endif
  
  return
end
