;
; PURPOSE: This routine loads into memory a specified traced-data structure.
;
; INPUTS:
; dir and structure_filename: STRINGS; dir where the SAV file is
; located and its filename.
;
; OUTPUT:
; trace_data: STRUCTURE; containing the tracing of tomographic
; products along field lines, as well as their geometry.
;
; HISTORY: V1.0, AMV, November 2023, IAFE.
;
pro load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data
  restore, filename = dir + structure_filename
  return
end
