;
; PURPOSE: This code reads into memory tomographic producs for either
;          VL tomography (N_e) or EUV tomography (N_e, T_e).
;          It then loops across AWSoM fieldline files and, for each one,
;          it calls the routine that actually reads each fieldline
;          file and it traces the tomographic producs along them.
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
; HISTORY: V1.0.1 FAN & AMV, CLaSP, May 2025.
; Write the tom. grid parameters in an ASCII file
;

pro tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
              nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax,csv=csv,$
              trace_Bs=trace_Bs

; Write in an ascii file the tom. grid parameters  
   openw,1,fl_dir+'tom.grid.'+instr+'.dat'
  printf,1,'      nr      nt     np       rmin         rmax         Irmin        Irmax'
  printf,1,nr,nt,np,rmin,rmax,Irmin,Irmax
   close,1

; Read Tom results to trace
  if instr eq 'aia' or instr eq 'euvia' or instr eq 'euvib' or instr eq 'eit' then $
     read_demt ,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e,T_e,W_T,ldem_flag
  if instr eq 'lascoc2' or instr eq 'mk4' or instr eq 'kcor' or instr eq 'ucomp' then $
     read_vlsrt,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e

; Read the list with the field-lines  
  N=0L
  filename=''
  openr,1,fl_dir+fl_list
  readf,1,N
  for i = 0L,N-1 do begin
     readf,1,filename
;    Read the fieldline and trace Tom. results along the line
     if instr eq 'aia' or instr eq 'euvia' or instr eq 'euvib' or instr eq 'eit' then $
        read_fieldline_and_trace_demt,instr,fl_dir,filename,rad,lat,lon,nr,nt,np,N_e,T_e,W_T,ldem_flag,csv=csv,trace_Bs=trace_Bs
     if instr eq 'lascoc2' or instr eq 'mk4' or instr eq 'kcor' or instr eq 'ucomp' then $
        read_fieldline_and_trace_vlsrt,instr,fl_dir,filename,rad,lat,lon,nr,nt,np,N_e,csv=csv,trace_Bs=trace_Bs
  endfor
  close,1
  
  return
end
