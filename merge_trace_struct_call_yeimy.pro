;
; PURPOSE: This code is the main routine to merge tracings, where the AWSoM fieldline
; files are provided, as well as flags are set for every instrument
; whose already traced tomography is to me merged.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
; 

pro merge_trace_struct_call_yeimy, nfs1=nfs1, nfs2=nfs2, trace_Bs=trace_Bs

; Define PROJECT_NAME, a string suffix to construct the full PATHS to the required files.
  PROJECT_NAME = 'April24'

; Provide FL_LIST, the file which informs the number of field lines and the
; filenames of the ASCII files containing the geometry of each line.
; fl_list = 'fdips_field_150X180X360_mrbqs240414t1304c2283_230_shift.ubdat_fline-filenames_list.txt'
  fl_list = 'list_yeimy-fl.txt'
; Define field_line_geometry_suffix_dir
  field_line_geometry_suffix_dir='_yeimy/'
  if not keyword_set(field_line_geometry_suffix_dir) then $
  field_line_geometry_suffix_dir='/'

; --------------------This block should not require edits.---------------------------
; Set  FL_DIR, where the field-lines geometry files should be located,
; and TOM_DIR, where the 3D tomography products to trace should be located.
  base_dir = '/data1/'
  if keyword_set(nfs1) then base_dir = '/data/Data1/data1/'
  if keyword_set(nfs2) then base_dir = '/data/Data2/data1/'
   fl_dir = base_dir+'DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry'+field_line_geometry_suffix_dir
;------------------------------------------------------------------------------------

   merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, trace_Bs = trace_Bs, /aia ,/lascoc2
   
  return
end

