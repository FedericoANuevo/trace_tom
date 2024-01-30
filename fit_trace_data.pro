;
; PURPOSE: Performs analytical fits to tomographic products along each
; traced field line and stores results in the trace_data structure.
;
; HISTORY: V1.0, AMV, January 2023, IAFE.
;          V1.1, AMV & FAN, January 2023, IAFE. Added Theil-Sen and R² metric.
;          V1.2, AMV, January 2023, IAFE. Added Double-Power-Law for Mk4.
;                                         Added fits' parameters to output structure.
;                                         Added 2-PoweLaw fit for C2.
;          V1.3, AMV, January 2023, IAFE. Added DPL fit for EUV.
;                                         Added fitted function name
;                                         and automated parameter
;                                         storing in output structure.
;          V1.4, AMV, January 2023, IAFE. Added <lambda_N> for all
;                                         functional fits.
;
pro fit_trace_data, aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                    mk4 = mk4, kcor = kcor, lascoc2 = lascoc2

    common all_data, trace_data, $
     N_fl, Npt_max, Npt_v,$
     x_A, y_A, z_A, rad_A, lat_A, lon_A,$
     Ne_aia_A, Tm_aia_A, index_aia_A, index_sampling_aia_A,$
     Ne_euvia_A, Tm_euvia_A, index_euvia_A, index_sampling_euvia_A,$
     Ne_euvib_A, Tm_euvib_A, index_euvib_A, index_sampling_euvib_A,$
     Ne_eit_A, Tm_eit_A, index_eit_A, index_sampling_eit_A,$
     Ne_mk4_A, index_mk4_A, index_sampling_mk4_A,$
     Ne_kcor_A, index_kcor_A, index_sampling_kcor_A,$
     Ne_c2_A, index_c2_A, index_sampling_c2_A

    common radcrits, radcritA, radcritB

    default = -678.
    
    if keyword_set(aia) then begin
       radmin = 1.0 & radmax = 1.3
       drad_fit = 0.01
       Npt_fit = round((radmax-radmin)/drad_fit)
       fitflag_aia_A = fltarr(N_fl        ) + default
        N0_fit_aia_A = fltarr(N_fl        ) + default
        lN_fit_aia_A = fltarr(N_fl        ) + default
        N1_fit_aia_A = fltarr(N_fl        ) + default
        N2_fit_aia_A = fltarr(N_fl        ) + default
        p1_fit_aia_A = fltarr(N_fl        ) + default
        p2_fit_aia_A = fltarr(N_fl        ) + default
       scN_fit_aia_A = fltarr(N_fl        ) + default
        T0_fit_aia_A = fltarr(N_fl        ) + default
      dTdr_fit_aia_A = fltarr(N_fl        ) + default
       scT_fit_aia_A = fltarr(N_fl        ) + default
        Ne_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_aia_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_aia_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_aia_A(ifl,*))
          ind_samp_aia = where(tmp eq 1)
          radsamp = reform(rad_A(ifl,ind_samp_aia)) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /aia
          if covgflag eq 'yes' then begin
             fitflag_aia_A(ifl) = +1.
             Nesamp = reform(Ne_aia_A(ifl,ind_samp_aia))
             Tmsamp = reform(Tm_aia_A(ifl,ind_samp_aia))
              goto,skip_aia_isohthermal_hydrostatic
                 fit_F_Ne_aia  = 'IHS'
                 linear_fit, 1./radsamp   , alog(Nesamp), AN, r2N, /linfit_idl
                 scN_fit_aia_A(ifl)   = r2N             
                  N0_fit_aia_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
                  lN_fit_aia_A(ifl)   = 1./AN[1]         ; Rsun
                  Ne_fit_aia_A(ifl,*) = N0_fit_aia_A(ifl) * exp(-(1/lN_fit_aia_A(ifl))*(1.-1./rad_fit_aia_A   )) ; cm-3
                  print,lN_fit_aia_A(ifl)
              skip_aia_isohthermal_hydrostatic:
             ;goto,skip_aia_double_power_law
                  fit_F_Ne_aia  = 'DPL'
                  double_power_fit, radsamp, Nesamp, A, chisq, /weighted
                 scN_fit_aia_A(ifl)   = chisq
                  N1_fit_aia_A(ifl)   = A[0] ; cm-3
                  p1_fit_aia_A(ifl)   = A[1] ; dimensionless exponent of power law
                  N2_fit_aia_A(ifl)   = A[2] ; cm-3
                  p2_fit_aia_A(ifl)   = A[3] ; dimensionless exponent of power law
                  Ne_fit_aia_A(ifl,*) = A[0] * rad_fit_aia_A^(-A[1]) + A[2] * rad_fit_aia_A^(-A[3]) ; cm-3
                  dNe_dr              = - A[1]*A[0] * rad_fit_aia_A^(-A[1]-1) - A[3]*A[2] * rad_fit_aia_A^(-A[3]-1) ; cm-3 / Rsun
                  indsamp = where(rad_fit_aia_A ge min(radsamp) and rad_fit_aia_A le max(radsamp) AND dNe_dr lt 0.)
                 ;indsamp = where(rad_fit_aia_A ge min(radsamp) and rad_fit_aia_A le max(1.2    ) AND dNe_dr lt 0.)
                  lN_fit_aia_A(ifl)   =  int_tabulated(rad_fit_aia_A(indsamp),abs(dNe_dr(indsamp)/reform(Ne_fit_aia_A(ifl,indsamp)))^(-1)) / (max(rad_fit_aia_A(indsamp))-min(rad_fit_aia_A(indsamp)))
              skip_aia_double_power_law:
             ;Linear fit to Te(r)
                  linear_fit,    radsamp-1.,      Tmsamp , AT, r2T, /theil_sen
                  scT_fit_aia_A(ifl)   = r2T
                   T0_fit_aia_A(ifl)   = AT[0]            ; K
                 dTdr_fit_aia_A(ifl)   = AT[1]            ; K/Rsun
                   Tm_fit_aia_A(ifl,*) = T0_fit_aia_A(ifl) + dTdr_fit_aia_A(ifl)       *      (rad_fit_aia_A-1.)  ; K
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                           ,$
                                  'fitflag_aia',ptr_new( fitflag_aia_A) ,$
                                 'fit_F_Ne_aia',ptr_new(  fit_F_Ne_aia) ,$
                                  'scN_fit_aia',ptr_new( scN_fit_aia_A) ,$
                                  'scT_fit_aia',ptr_new( scT_fit_aia_A) ,$
                                  'rad_fit_aia',ptr_new( rad_fit_aia_A) ,$
                                   'Ne_fit_aia',ptr_new(  Ne_fit_aia_A) ,$
                                   'Tm_fit_aia',ptr_new(  Tm_fit_aia_A) ,$
                                   'T0_fit_aia',ptr_new(  T0_fit_aia_A) ,$
                                 'dTdr_fit_aia',ptr_new(dTdr_fit_aia_A) ,$
                                   'lN_fit_aia',ptr_new(  lN_fit_aia_A) )
       if fit_F_Ne_aia eq 'IHS' then $
          trace_data = create_struct( trace_data                        ,$
                                   'N0_fit_aia',ptr_new(  N0_fit_aia_A) )                                                                        
       if fit_F_Ne_aia eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_aia',ptr_new( N1_fit_aia_A) ,$
                                    'N2_fit_aia',ptr_new( N2_fit_aia_A) ,$
                                    'p1_fit_aia',ptr_new( p1_fit_aia_A) ,$
                                    'p2_fit_aia',ptr_new( p2_fit_aia_A) )

    endif                       ; AIA

    if keyword_set(mk4) then begin
       radmin = 1.1 & radmax = 1.5
       drad_fit = 0.01
       Npt_fit = round((radmax-radmin)/drad_fit)
       fitflag_mk4_A = fltarr(N_fl        ) + default
        N0_fit_mk4_A = fltarr(N_fl        ) + default
        lN_fit_mk4_A = fltarr(N_fl        ) + default
        N1_fit_mk4_A = fltarr(N_fl        ) + default
        N2_fit_mk4_A = fltarr(N_fl        ) + default
        p1_fit_mk4_A = fltarr(N_fl        ) + default
        p2_fit_mk4_A = fltarr(N_fl        ) + default
       scN_fit_mk4_A = fltarr(N_fl        ) + default
        Ne_fit_mk4_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_mk4_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_mk4_A(ifl,*))
          ind_samp_mk4 = where(tmp eq 1)
          radsamp = reform(rad_A(ifl,ind_samp_mk4)) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /mk4
          if covgflag eq 'yes' then begin
             fitflag_mk4_A(ifl) = +1.
             Nesamp = reform(Ne_mk4_A(ifl,ind_samp_mk4))
             goto,skip_mk4_isohthermal_hydrostatic
                 fit_F_Ne_mk4  = 'IHS'
                  linear_fit, 1./radsamp   ,alog(Nesamp), AN, r2N, /linfit_idl
                  scN_fit_mk4_A(ifl)   = r2N
                  N0_fit_mk4_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
                  lN_fit_mk4_A(ifl)   = 1./AN[1]         ; Rsun
                  Ne_fit_mk4_A(ifl,*) = N0_fit_mk4_A(ifl) * exp(-(1/lN_fit_mk4_A(ifl))*(1.-1./rad_fit_mk4_A)) ; cm-3
             skip_mk4_isohthermal_hydrostatic:
             goto,skip_mk4_single_power_law
                  fit_F_Ne_mk4  = 'SPL'
                  linear_fit, alog(radsamp), alog(Nesamp), AN, r2N, /linfit_idl
                 scN_fit_mk4_A(ifl)   = r2N
                  N1_fit_mk4_A(ifl)   = exp(AN[0]) ; cm-3
                  p1_fit_mk4_A(ifl)   =    -AN[1]  ; dimensionless exponent of power law
                  Ne_fit_mk4_A(ifl,*) = N1_fit_mk4_A(ifl) * rad_fit_mk4_A^(-p1_fit_mk4_A(ifl)) ; cm-3
                  indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp))
                  lN_fit_mk4_A(ifl)   = int_tabulated( rad_fit_mk4_A(indsamp), rad_fit_mk4_A(indsamp) / reform(p1_fit_mk4_A(ifl)) ) / (max(rad_fit_mk4_A(indsamp))-min(rad_fit_mk4_A(indsamp))) ; Rsun
             skip_mk4_single_power_law:
            ;goto,skip_mk4_double_power_law
                  fit_F_Ne_mk4  = 'DPL'
                  double_power_fit, radsamp, Nesamp, A, chisq, /weighted
                 scN_fit_mk4_A(ifl)   = chisq
                  N1_fit_mk4_A(ifl)   = A[0] ; cm-3
                  p1_fit_mk4_A(ifl)   = A[1] ; dimensionless exponent of power law
                  N2_fit_mk4_A(ifl)   = A[2] ; cm-3
                  p2_fit_mk4_A(ifl)   = A[3] ; dimensionless exponent of power law
                  Ne_fit_mk4_A(ifl,*) = A[0] * rad_fit_mk4_A^(-A[1]) + A[2] * rad_fit_mk4_A^(-A[3]) ; cm-3
                  dNe_dr              = - A[1]*A[0] * rad_fit_mk4_A^(-A[1]-1) - A[3]*A[2] * rad_fit_mk4_A^(-A[3]-1) ; cm-3 / Rsun
                  indsamp = where(rad_fit_mk4_A ge min(radsamp) and rad_fit_mk4_A le max(radsamp) AND dNe_dr lt 0.)
                  lN_fit_mk4_A(ifl)   = int_tabulated(rad_fit_mk4_A(indsamp),abs(dNe_dr(indsamp)/reform(Ne_fit_mk4_A(ifl,indsamp)))^(-1)) / (max(rad_fit_mk4_A(indsamp))-min(rad_fit_mk4_A(indsamp))) ; Rsun
            skip_mk4_double_power_law:
         endif             ; covgflag = 'yes'          
       endfor              ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_mk4',ptr_new(fitflag_mk4_A) ,$
                                 'fit_F_Ne_mk4',ptr_new( fit_F_Ne_mk4) ,$
                                  'scN_fit_mk4',ptr_new(scN_fit_mk4_A) ,$
                                  'rad_fit_mk4',ptr_new(rad_fit_mk4_A) ,$
                                   'Ne_fit_mk4',ptr_new( Ne_fit_mk4_A) ,$
                                   'lN_fit_mk4',ptr_new(  lN_fit_mk4_A) )
       if fit_F_Ne_mk4 eq 'IHS' then $
          trace_data = create_struct( trace_data                        ,$
                                   'N0_fit_mk4',ptr_new(  N0_fit_mk4_A) )                                                                        
       if fit_F_Ne_mk4 eq 'SPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_mk4',ptr_new( N1_fit_mk4_A) ,$
                                    'p1_fit_mk4',ptr_new( p1_fit_mk4_A) )
       if fit_F_Ne_mk4 eq 'DPL' then $
          trace_data = create_struct( trace_data                        ,$
                                    'N1_fit_mk4',ptr_new( N1_fit_mk4_A) ,$
                                    'N2_fit_mk4',ptr_new( N2_fit_mk4_A) ,$
                                    'p1_fit_mk4',ptr_new( p1_fit_mk4_A) ,$
                                    'p2_fit_mk4',ptr_new( p2_fit_mk4_A) )
    endif ; MK4

    if keyword_set(lascoc2) then begin
       radmin = 2.5 & radmax = 6.0
       drad_fit = 0.1
       Npt_fit = round((radmax-radmin)/drad_fit)
       fitflag_c2_A = fltarr(N_fl        ) + default
        N1_fit_c2_A = fltarr(N_fl        ) + default
        N2_fit_c2_A = fltarr(N_fl        ) + default
        p1_fit_c2_A = fltarr(N_fl        ) + default
        p2_fit_c2_A = fltarr(N_fl        ) + default
        lN_fit_c2_A = fltarr(N_fl        ) + default
       scN_fit_c2_A = fltarr(N_fl        ) + default
        Ne_fit_c2_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_c2_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_c2_A(ifl,*))
          ind_samp_c2 = where(tmp eq 1)
          radsamp = reform(rad_A(ifl,ind_samp_c2)) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /lascoc2
          if covgflag eq 'yes' then begin
             fitflag_c2_A(ifl) = +1.
             Nesamp = reform(Ne_c2_A(ifl,ind_samp_c2))
             goto,skip_c2_single_power_law
                 fit_F_Ne_c2 = 'SPL'
                 linear_fit, alog(radsamp), alog(Nesamp), AN, r2N, /linfit_idl
                 scN_fit_c2_A(ifl)  = r2N
                 N1_fit_c2_A(ifl)   = exp(AN[0]) ; cm-3
                 p1_fit_c2_A(ifl)   =    -AN[1]  ; dimensionless exponent of power law
                 Ne_fit_c2_A(ifl,*) = N1_fit_c2_A(ifl) * rad_fit_c2_A^(-p1_fit_c2_A(ifl)) ; cm-3
                 indsamp = where(rad_fit_c2_A ge min(radsamp) and rad_fit_c2_A le max(radsamp))
                 lN_fit_c2_A(ifl)   = int_tabulated( rad_fit_c2_A(indsamp), rad_fit_c2_A(indsamp) / reform(p1_fit_c2_A(ifl)) ) / (max(rad_fit_c2_A(indsamp))-min(rad_fit_c2_A(indsamp))) ; Rsun
             skip_c2_single_power_law:
            ;goto,skip_c2_double_power_law
                 fit_F_Ne_c2  = 'DPL'
                 double_power_fit, radsamp, Nesamp, A, chisq, /weighted
                 scN_fit_c2_A(ifl)   = chisq
                  N1_fit_c2_A(ifl)   = A[0] ; cm-3
                  p1_fit_c2_A(ifl)   = A[1] ; dimensionless exponent of power law
                  N2_fit_c2_A(ifl)   = A[2] ; cm-3
                  p2_fit_c2_A(ifl)   = A[3] ; dimensionless exponent of power law
                  Ne_fit_c2_A(ifl,*) = A[0] * rad_fit_c2_A^(-A[1]) + A[2] * rad_fit_c2_A^(-A[3]) ; cm-3
                  dNe_dr              = - A[1]*A[0] * rad_fit_c2_A^(-A[1]-1) - A[3]*A[2] * rad_fit_c2_A^(-A[3]-1) ; cm-3 / Rsun
                  indsamp = where(rad_fit_c2_A ge min(radsamp) and rad_fit_c2_A le max(radsamp) AND dNe_dr lt 0.)
                  lN_fit_c2_A(ifl)   = int_tabulated(rad_fit_c2_A(indsamp),abs(dNe_dr(indsamp)/reform(Ne_fit_c2_A(ifl,indsamp)))^(-1)) / (max(rad_fit_c2_A(indsamp))-min(rad_fit_c2_A(indsamp))) ; Rsun
             skip_c2_double_power_law:
           endif                ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                        ,$
                                  'fitflag_c2',ptr_new(fitflag_c2_A) ,$
                                 'fit_F_Ne_c2',ptr_new( fit_F_Ne_c2) ,$
                                  'scN_fit_c2',ptr_new(scN_fit_c2_A) ,$                                   
                                  'rad_fit_c2',ptr_new(rad_fit_c2_A) ,$
                                   'Ne_fit_c2',ptr_new( Ne_fit_c2_A) ,$
                                   'lN_fit_c2',ptr_new( lN_fit_c2_A) )
       if fit_F_Ne_c2 eq 'SPL' then $
          trace_data = create_struct( trace_data                     ,$
                                   'N1_fit_c2',ptr_new( N1_fit_c2_A) ,$
                                   'p1_fit_c2',ptr_new( p1_fit_c2_A) )
       if fit_F_Ne_c2 eq 'DPL' then $
          trace_data = create_struct( trace_data                     ,$
                                   'N1_fit_c2',ptr_new( N1_fit_c2_A) ,$
                                   'N2_fit_c2',ptr_new( N2_fit_c2_A) ,$
                                   'p1_fit_c2',ptr_new( p1_fit_c2_A) ,$
                                   'p2_fit_c2',ptr_new( p2_fit_c2_A) )
    endif; LASCOC2

    return
 end
