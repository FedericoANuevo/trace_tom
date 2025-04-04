;
; PURPOSE: This code is the main routine to merge tracings, where the AWSoM fieldline
; files are provided, as well as flags are set for every instrument
; whose already traced tomography is to me merged.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
; 

pro merge_trace_struct_call

; dir_fl  = '/data1/DATA/fieldlines_judit/radial_synth_fieldlines/' & fl_list = 'list_synth.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2099/map1/'             & fl_list = 'list.map1.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2099/map7/'             & fl_list = 'list.map7.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2099/map12/'            & fl_list = 'list.map12.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2099/map7_new/'         & fl_list = 'list.map7.new.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2099/map1_new/'         & fl_list = 'list.map1.new.txt' 
; merge_trace_struct, dir_fl = dir_fl, fl_list = fl_list,/aia,/mk4,/lascoc2
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2082/map1/'             & fl_list = 'list.map1.txt'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2082/map1_new/'         & fl_list = 'list.map1.new.txt' &  structure_filename = 'CR2082_AWSoM-map1'
; dir_fl  = '/data1/DATA/fieldlines_judit/CR2082/map7_new/'         & fl_list = 'list.map7.new.txt' &  structure_filename = 'CR2082_AWSoM-map7'
; merge_trace_struct, dir_fl = dir_fl, fl_list = fl_list,/euvia,/lascoc2,struture_filename = structure_filename
  dir_fl  = '/data1/DATA/flines_Sam-Yeimy/'         & fl_list = 'list.txt' &  structure_filename = 'April204_PFSS'
  merge_trace_struct, dir_fl = dir_fl, fl_list = fl_list,/aia,/lascoc2,struture_filename = structure_filename
  
  return
end
