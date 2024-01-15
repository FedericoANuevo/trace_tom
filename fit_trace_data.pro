;
; PURPOSE: Performs analytical fits to tomographic products along each
; traced field line and stores results in the trace_data structure.
;
; HISTORY: V1.0, AMV, January 2023, IAFE.
;          V1.1, AMV & FAN, January 2023, IAFE. Added Theil-Sen and R² metric.
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
       r2N_fit_aia_A = fltarr(N_fl        ) + default
        T0_fit_aia_A = fltarr(N_fl        ) + default
      dTdr_fit_aia_A = fltarr(N_fl        ) + default
       r2T_fit_aia_A = fltarr(N_fl        ) + default
        Ne_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_aia_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_aia_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_aia_A(ifl,*))
          ind_samp_aia = where(tmp eq 1)
          radsamp = rad_A(ifl,ind_samp_aia) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /aia
          if covgflag eq 'yes' then begin
             fitflag_aia_A(ifl) = +1.
             Nesamp = Ne_aia_A(ifl,ind_samp_aia)
             Tmsamp = Tm_aia_A(ifl,ind_samp_aia)
             linear_fit, 1./radsamp   , alog(Nesamp), AN, r2N, /linfit_idl
             linear_fit,    radsamp-1.,      Tmsamp , AT, r2T, /theil_sen
             r2N_fit_aia_A(ifl)   = r2N
             r2T_fit_aia_A(ifl)   = r2T
              N0_fit_aia_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
              lN_fit_aia_A(ifl)   = 1./AN[1]         ; Rsun
              T0_fit_aia_A(ifl)   = AT[0]            ; K
            dTdr_fit_aia_A(ifl)   = AT[1]            ; K/Rsun
              Ne_fit_aia_A(ifl,*) = N0_fit_aia_A(ifl) * exp(-(1/lN_fit_aia_A(ifl))*(1.-1./rad_fit_aia_A   )) ; cm-3
              Tm_fit_aia_A(ifl,*) = T0_fit_aia_A(ifl) + dTdr_fit_aia_A(ifl)       *      (rad_fit_aia_A-1.)  ; K
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_aia',ptr_new(fitflag_aia_A) ,$
                                  'r2N_fit_aia',ptr_new(r2N_fit_aia_A) ,$
                                  'r2T_fit_aia',ptr_new(r2T_fit_aia_A) ,$
                                  'rad_fit_aia',ptr_new(rad_fit_aia_A) ,$
                                   'Ne_fit_aia',ptr_new( Ne_fit_aia_A) ,$
                                   'Tm_fit_aia',ptr_new( Tm_fit_aia_A) )
    endif; AIA

    if keyword_set(mk4) then begin
       radmin = 1.15 & radmax = 1.5
       drad_fit = 0.01
       Npt_fit = round((radmax-radmin)/drad_fit)
       fitflag_mk4_A = fltarr(N_fl        ) + default
        N0_fit_mk4_A = fltarr(N_fl        ) + default
        lN_fit_mk4_A = fltarr(N_fl        ) + default
         p_fit_mk4_A = fltarr(N_fl        ) + default
       r2N_fit_mk4_A = fltarr(N_fl        ) + default
        Ne_fit_mk4_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_mk4_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_mk4_A(ifl,*))
          ind_samp_mk4 = where(tmp eq 1)
          radsamp = rad_A(ifl,ind_samp_mk4) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /mk4
          if covgflag eq 'yes' then begin
             fitflag_mk4_A(ifl) = +1.
             Nesamp = Ne_mk4_A(ifl,ind_samp_mk4)
             goto,skip_isohthermal_hydrostatic
                 linear_fit, 1./radsamp   ,alog(Nesamp), AN, r2N, /linfit_idl
                 r2N_fit_mk4_A(ifl)   = r2N
                  N0_fit_mk4_A(ifl)   = exp(AN[0]+AN[1]) ; cm-3
                  lN_fit_mk4_A(ifl)   = 1./AN[1]         ; Rsun
                  Ne_fit_mk4_A(ifl,*) = N0_fit_mk4_A(ifl) * exp(-(1/lN_fit_mk4_A(ifl))*(1.-1./rad_fit_mk4_A   )) ; cm-3
             skip_isohthermal_hydrostatic:
            ;goto,single_power_law
                 radcrit = radmin
                 linear_fit, alog(radsamp/radcrit), alog(Nesamp), AN, r2N, /linfit_idl
                 r2N_fit_mk4_A(ifl)   = r2N
                  N0_fit_mk4_A(ifl)   = exp(AN[0]) ; cm-3
                   p_fit_mk4_A(ifl)   =    -AN[1]  ; dimensionless exponent of power law
                  Ne_fit_mk4_A(ifl,*) = N0_fit_mk4_A(ifl) * (rad_fit_mk4_A / radcrit)^(-p_fit_mk4_A(ifl)) ; cm-3
             single_power_law:
             goto,double_power_law
             double_power_fit, radmin, radmax, radsamp, Nesamp, A, chisq
                 r2N_fit_mk4_A(ifl)   = chisq
                  Ne_fit_mk4_A(ifl,*) = A[0] * (rad_fit_mk4_A / radcritA)^(-A[1]) + A[2] * (rad_fit_mk4_A / radcritB)^(-A[3]) ; cm-3
             double_power_law:
         endif             ; covgflag = 'yes'          
       endfor              ; field lines loop.
       trace_data = create_struct( trace_data                          ,$
                                  'fitflag_mk4',ptr_new(fitflag_mk4_A) ,$
                                  'r2N_fit_mk4',ptr_new(r2N_fit_mk4_A) ,$
                                  'rad_fit_mk4',ptr_new(rad_fit_mk4_A) ,$
                                   'Ne_fit_mk4',ptr_new( Ne_fit_mk4_A) )
    endif; MK4

    if keyword_set(lascoc2) then begin
       radmin = 2.5 & radmax = 6.0
       drad_fit = 0.1
       Npt_fit = round((radmax-radmin)/drad_fit)
       radcrit = radmin
       fitflag_c2_A = fltarr(N_fl        ) + default
        N0_fit_c2_A = fltarr(N_fl        ) + default
         p_fit_c2_A = fltarr(N_fl        ) + default
       r2N_fit_c2_A = fltarr(N_fl        ) + default
        Ne_fit_c2_A = fltarr(N_fl,Npt_fit) + default
       rad_fit_c2_A = radmin + drad_fit/2. + drad_fit * findgen(Npt_fit)
       for ifl=0,N_fl-1 do begin      
          tmp = reform(index_sampling_c2_A(ifl,*))
          ind_samp_c2 = where(tmp eq 1)
          radsamp = rad_A(ifl,ind_samp_c2) ; Rsun
          test_coverage, radsamp=radsamp, covgflag=covgflag, /lascoc2
          if covgflag eq 'yes' then begin
             fitflag_c2_A(ifl) = +1.
             Nesamp = Ne_c2_A(ifl,ind_samp_c2)
             linear_fit, alog(radsamp/radcrit), alog(Nesamp), AN, r2N, /linfit_idl
             r2N_fit_c2_A(ifl)   = r2N
              N0_fit_c2_A(ifl)   = exp(AN[0]) ; cm-3
               p_fit_c2_A(ifl)   =    -AN[1]  ; dimensionless exponent of power law
              Ne_fit_c2_A(ifl,*) = N0_fit_c2_A(ifl) * (rad_fit_c2_A / radcrit)^(-p_fit_c2_A(ifl)) ; cm-3
          endif                 ; covgflag = 'yes'          
       endfor                   ; field lines loop.
       trace_data = create_struct( trace_data                        ,$
                                  'fitflag_c2',ptr_new(fitflag_c2_A) ,$
                                  'r2N_fit_c2',ptr_new(r2N_fit_c2_A) ,$                                   
                                  'rad_fit_c2',ptr_new(rad_fit_c2_A) ,$
                                   'Ne_fit_c2',ptr_new( Ne_fit_c2_A) )
    endif; LASCOC2

    return
 end

