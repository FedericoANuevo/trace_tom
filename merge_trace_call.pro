;
; PURPOSE: This code is the main routine to merge tracings, where the AWSoM fieldline
; files are provided, as well as flags are set for every instrument
; whose already traced tomography is to me merged.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
; 

pro merge_trace_call
  
  dir_fl  = '/data1/DATA/fieldlines_judit/'
  fl_list = 'list_synth.txt'

  merge_trace, dir_fl = dir_fl, fl_list = fl_list, /aia, /mk4, /lascoc2
  
  return
end
