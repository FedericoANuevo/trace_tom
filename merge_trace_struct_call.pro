;
; PURPOSE:
; This code is the main routine to merge tracings. Here the AWSoM fieldline
; files are provided, as well as flags are set for every instrument
; for which already traced tomographic prodcuts are to me merged.
;
; HISTORY: V1.0 AMV & FAN, CLaSP, October 2023.
; 
; merge_trace_struct_call, /trace_Bs

pro merge_trace_struct_call, nfs1=nfs1, nfs2=nfs2, trace_Bs=trace_Bs

;===============================================================================================
; Define PROJECT_NAME, a string suffix to construct the full PATHS to the required files.
; PROJECT_NAME = 'CR2254'
; PROJECT_NAME = 'CR2261'
  PROJECT_NAME = 'April24'
  
; Define field_line_geometry_suffix_dir
; field_line_geometry_suffix_dir = '_aunifgrid_multirad_3x3deg_HMI-PolFil/'                                    
; field_line_geometry_suffix_dir = '_aunifgrid_2.50Rs_1x1deg_HMI-PolFil/'
  field_line_geometry_suffix_dir = '_equatorial-ring/'                                    
  
; Provide FL_LIST, the file which informs the number of field lines and the
; filenames of the ASCII files containing the geometry of each line.
; fl_list = 'fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2254_prep.ubdat_fline-filenames_list.txt'
; fl_list = 'fdips_field_150X180X360_hmi.Synoptic_Mr_polfil.2261_prep.ubdat_fline-filenames_list.txt'
  fl_list = 'Bfield_AWSoM_April24.ubdat_fline-filenames_list.txt'  
;===============================================================================================


; --------------------This block should not require edits.---------------------------
; Set  FL_DIR, where the field-lines geometry files should be located,
; and TOM_DIR, where the 3D tomography products to trace should be located.
  base_dir = '/data1/'
  if keyword_set(nfs1) then base_dir = '/data/Data1/data1/'
  if keyword_set(nfs2) then base_dir = '/data/Data2/data1/'
  if not keyword_set(field_line_geometry_suffix_dir) then field_line_geometry_suffix_dir='/'
  fl_dir = base_dir+'DATA/trace_tom_files/'+PROJECT_NAME+'/field_lines_geometry'+field_line_geometry_suffix_dir
;------------------------------------------------------------------------------------

;  merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, trace_Bs=trace_Bs, /aia, /kcor, /ucomp, structure_filename=structure_filename
   merge_trace_struct, fl_dir=fl_dir, fl_list=fl_list, trace_Bs=trace_Bs, /aia, /lascoc2     , structure_filename=structure_filename
  return
end
