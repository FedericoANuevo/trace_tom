;
; PURPOSE: This routine loads into memory a specified traced-data structure.
;
; INPUTS:
; dir and structure_filename: STRINGS; dir where the SAV file is
; located and its filename.
;
; OUTPUTS:
; trace_data: STRUCTURE; containing the tracing of tomographic
; products along field lines, as well as their geometry.
;
; HISTORY: V1.0, AMV, November 2023, IAFE.
;
pro load_traced_data_structure, dir=dir, structure_filename=structure_filename, trace_data=trace_data
  common data, N_fl, Npt_max
  restore, filename = dir + structure_filename
  N_fl    = *trace_data.N_fl
  Npt_max = *trace_data.Npt_max
  return
end
