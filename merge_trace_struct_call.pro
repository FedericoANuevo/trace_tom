;
; PURPOSE: This code is the main routine to merge tracings, where the AWSoM fieldline
; files are provided, as well as flags are set for every instrument
; whose already traced tomography is to me merged.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
; 

pro merge_trace_struct_call, nfs1=nfs1, nfs2=nfs2, opcl=opcl

; Define PROJECT_NAME, a string suffix to construct the full PATHS to the required files.
  PROJECT_NAME = 'CR2254'

; Provide FL_LIST, the file which informs the number of field lines and the
; filenames of the ASCII files containing the geometry of each line.
  fl_list = 'fdips_field_150x180x360_mrmqs220221t2004c2254_000.ubdat_fline-filenames_list.txt'

; --------------------This block should not require edits.---------------------------
; Set  FL_DIR, where the field-lines geometry files should be located,
; and TOM_DIR, where the 3D tomography products to trace should be located.
  base_dir = '/data1/'
  if keyword_set(nfs1) then base_dir = '/data/Data1/data1/'
  if keyword_set(nfs2) then base_dir = '/data/Data2/data1/'
   fl_dir = base_dir+'DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry/'
;------------------------------------------------------------------------------------

 ; merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, opcl=opcl, /aia, /kcor, /ucomp, struture_filename=structure_filename
   merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, opcl=opcl, /aia               , struture_filename=structure_filename
  return
end

; fl_dir  = '/data1/DATA/fieldlines_judit/radial_synth_fieldlines/' & fl_list = 'list_synth.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2099/map1/'             & fl_list = 'list.map1.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2099/map7/'             & fl_list = 'list.map7.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2099/map12/'            & fl_list = 'list.map12.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2099/map7_new/'         & fl_list = 'list.map7.new.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2099/map1_new/'         & fl_list = 'list.map1.new.txt' 
; merge_trace_struct, fl_dir = fl_dir, fl_list = fl_list,/aia,/mk4,/lascoc2
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2082/map1/'             & fl_list = 'list.map1.txt'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2082/map1_new/'         & fl_list = 'list.map1.new.txt' &  structure_filename = 'CR2082_AWSoM-map1'
; fl_dir  = '/data1/DATA/fieldlines_judit/CR2082/map7_new/'         & fl_list = 'list.map7.new.txt' &  structure_filename = 'CR2082_AWSoM-map7'
; merge_trace_struct, fl_dir = fl_dir, fl_list = fl_list,/euvia,/lascoc2,struture_filename = structure_filename
; fl_dir  = '/data1/DATA/flines_Sam-Yeimy/'         & fl_list = 'list.txt' &  structure_filename = 'April204_PFSS'
