;
; PURPOSE: This code is the main routine, where the AWSoM fieldline
; files are provided, as well as all the needed parameters and files
; to read the tomographic producs to be traced.
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;

; CALLING SEQUENCE EXAMPLES:
; tom_trace_call,/demt
; tom_trace_call,/lasco
; tom_trace_call,/mlso

pro tom_trace_call,demt=demt,lasco=lasco,mlso=mlso

; base_dir = '/data1/DATA/fieldlines_judit/'
  base_dir = '/data1/DATA/'
; base_dir = '/media/Data1/data1/DATA/'
; dir      = 'radial_synth_fieldlines/' & fl_list='list_synth.txt'
; dir      = 'CR2099/map1/' & fl_list  = 'list.map1.txt' 
; dir      = 'CR2099/map7/' & fl_list  = 'list.map7.txt'
; dir      = 'CR2099/map12/'& fl_list  = 'list.map12.txt'
; dir      = 'CR2082/map1/' & fl_list  = 'list.map1.txt'
; dir      = 'flines_Sam-Yeimy/' & fl_list  = 'list.txt'
; dir      = 'CR2082/map1_new/' & fl_list  = 'list.map1.new.txt'
; dir      = 'CR2082/map7_new/' & fl_list  = 'list.map7.new.txt'
; dir      = 'CR2099/map7_new/' & fl_list  = 'list.map7.new.txt'
; dir      = 'CR2099/map1_new/' & fl_list  = 'list.map1.new.txt' 
  dir      = 'fl_fdips/CR2254/' & fl_list  = 'fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt'
  fl_dir   = base_dir + dir

  if keyword_set(demt) then begin
     tom_dir  = '/data1/DATA/ldem_files/'
    ;tom_file = 'CR2082_HOLLOW_compound2.dat' & instr = 'euvia'
    ;tom_file = 'LDEM.CR2099_aia_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound1.dat'    & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound2.dat'    & instr = 'aia'
    ;tom_file = 'LDEM.April-2024_aia_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & instr = 'aia'
     tom_file = 'LDEM.feb-mar_2022_segment2_aia_Hollow_3Bands_ucomp_comparison_gauss1_lin_Norm-median_singlStart'
     nr       = 21;30
     nt       = 60;90
     np       = 2*nt
     rmin     = 1.09;1.0
     rmax     = 1.3;1.3
     Irmin    = 1.09;1.02
     Irmax    = 1.25;1.25
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax;,/csv
  endif

  if keyword_set(lasco) then begin
     instr    = 'lascoc2'
     tom_dir  = '/data1/tomography/bindata/'
 ;   tom_file = 'x_LascoC2pB_CR2099_shifted_std-grid_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L6.e-6'
 ;   tom_file = 'x_LascoC2pB_CR2082_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L8.2e-6'
     tom_file = 'x_LascoC2pB_April-2024_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L1.1e-5'
     nr       = 60
     nt       = 60
     np       = 2*nt
     rmin     = 2.5
     rmax     = 8.5
     Irmin    = 2.5
     Irmax    = 6.0
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,/csv
  endif

  if keyword_set(mlso) then begin
;    instr    = 'mk4'
;    tom_dir  = '/data1/tomography/bindata/'
;    tom_file = 'x_Mk4_CR2099_shifted_Rmin1.15_Rmax1.85_IRmin1.15_IRmax1.50_70x90x180_BF2_L5.e-6'
     instr    = 'kcor'
     tom_dir  = '/data1/tomography/bindata/'
     tom_file = 'x_kcor_2022_Feb-Mar_segment2_Rmin1.09_Rmax1.75_IRmin1.09_IRmax1.50_66x60x120_BF2_L2.6e-5'
     nr       = 66;70
     nt       = 60;90
     np       = 2*nt
     rmin     = 1.09;1.15
     rmax     = 1.75;1.85
     Irmin    = 1.09;1.15
     Irmax    = 1.50;1.50
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax
  endif
  


 if keyword_set(ucomp) then begin
     instr    = 'ucomp'
     tom_dir  = '/data1/DATA/ldem_files/'
     tom_file = 'Ne_ratio_ucomp_1074-1079_second_target_segment2.dat'
     nr       = 16
     nt       = 60
     np       = 2*nt
     rmin     = 1.09
     rmax     = 1.20
     Irmin    = 1.09
     Irmax    = 1.25
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax
  endif
  
  return
end
