;
; PURPOSE: This code reads into memory tomographic producs for either
;          VL tomography (N_e) or EUV tomography (N_e, T_e).
;          It then loops across AWSoM fieldline files and, for each one,
;          it calls the routine that actually reads each fieldline
;          file and it traces the tomographic producs along them.
;
; HISTORY: V1.0 FAN & AMV, CLaSP, October 2023.
;

pro tom_trace,instr=instr,tom_dir=tom_dir,tom_file=tom_file,fl_dir=fl_dir,fl_list=fl_list,$
              nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,Irmin=Irmin,Irmax=Irmax

; Read Tom results to trace
  if instr eq 'aia' or instr eq 'euvia' or instr eq 'euvib' or instr eq 'eit' then $
     read_demt ,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e,T_e,W_T,ldem_flag
  if instr eq 'lascoc2' or instr eq 'mk4' or instr eq 'kcor' then $
     read_vlsrt,tom_dir,tom_file,nr,nt,np,rmin,rmax,Irmin,Irmax,rad,lat,lon,N_e
; Read the list with the field-lines  
  N=0
  filename=''
  openr,1,fl_dir+fl_list
  readf,1,N
  for i = 0,N-1 do begin
     readf,1,filename
;    Read the fieldline and trace Tom. results along the line
     if instr eq 'aia' or instr eq 'euvia' or instr eq 'euvib' or instr eq 'eit' then $
        read_fieldline_and_trace_demt,instr,fl_dir,filename,rad,lat,lon,nr,nt,np,N_e,T_e,W_T,ldem_flag
     if instr eq 'lascoc2' or instr eq 'mk4' or instr eq 'kcor' then $
        read_fieldline_and_trace_vlsrt,instr,fl_dir,filename,rad,lat,lon,nr,nt,np,N_e
  endfor
  close,1
  
  return
end
