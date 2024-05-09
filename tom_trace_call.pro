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

  base_dir = '/data1/DATA/fieldlines_judit/'
; dir      = 'radial_synth_fieldlines/' & fl_list='list_synth.txt'
; dir      = 'CR2099/map1/' & fl_list  = 'list.map1.txt' 
; dir      = 'CR2099/map7/' & fl_list  = 'list.map7.txt'
; dir      = 'CR2099/map12/'& fl_list  = 'list.map12.txt'
  dir      = 'CR2082/map1/' & fl_list  = 'list.map1.txt'
  fl_dir   = base_dir + dir



  
  if keyword_set(demt) then begin
     tom_dir  = '/data1/DATA/ldem_files/'
     tom_file = 'CR2082_HOLLOW_compound2.dat' & instr = 'euvia'
    ;tom_file = 'LDEM.CR2099_aia_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound1.dat'    & instr = 'aia'
    ;tom_file = 'CR2099_AIA_compound2.dat'    & instr = 'aia'
     nr       = 30
     nt       = 90
     np       = 2*nt
     rmin     = 1.0
     rmax     = 1.3
     Irmin    = 1.02
     Irmax    = 1.25
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax
  endif

  if keyword_set(lasco) then begin
     instr    = 'lascoc2'
     tom_dir  = '/data1/tomography/bindata/'
     tom_file = 'x_LascoC2pB_CR2099_shifted_std-grid_Rmin2.5_Rmax8.5_IRmin2.5_IRmax6.0_60x60x120_BF4_L6.e-6'
     nr       = 60
     nt       = 60
     np       = 2*nt
     rmin     = 2.5
     rmax     = 8.5
     Irmin    = 2.5
     Irmax    = 6.0
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax
  endif

  if keyword_set(mlso) then begin
     instr    = 'mk4'
     tom_dir  = '/data1/tomography/bindata/'
     tom_file = 'x_Mk4_CR2099_shifted_Rmin1.15_Rmax1.85_IRmin1.15_IRmax1.50_70x90x180_BF2_L5.e-6'
     nr       = 70
     nt       = 90
     np       = 2*nt
     rmin     = 1.15
     rmax     = 1.85
     Irmin    = 1.15
     Irmax    = 1.50
     tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
               nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax
  endif
  
  return
end
