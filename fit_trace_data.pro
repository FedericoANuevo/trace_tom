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

    default = -678.
    
    if keyword_set(aia) then begin
       Npt_fit = 30 & rmin = 1.0 & rmax = 1.3
       fitflag_aia_A = fltarr(N_fl        ) + default
        N0_fit_aia_A = fltarr(N_fl        ) + default
        lN_fit_aia_A = fltarr(N_fl        ) + default
       r2N_fit_aia_A = fltarr(N_fl        ) + default

        Ne_fit_aia_A = fltarr(N_fl,Npt_fit) + default
        Tm_fit_aia_A = fltarr(N_fl,Npt_fit) + default
       dr_fit       = (rmax-rmin)/float(Npt_fit)
        r_fit_aia   = rmin + dr_fit/2. + dr_fit * findgen(Npt_fit)
       for ifl=0,Nfl-1 do begin      
          tmp = reform(index_sampling_aia_A(ifl,*))
          ind_samp_aia = where(tmp eq 1)
          rfit =    rad_A(ifl,ind_samp_aia) ; Rsun
          test_coverage, rfit=rfit, covgflag=covgflag, /aia
          if covgflag eq 'yes' then begin
             fitflag_aia_A(ifl) = +1.
             Nefit = Ne_aia_A(ifl,ind_samp_aia)
             Tmfit = Tm_aia_A(ifl,ind_samp_aia)
             A = linfit(1/rfit,alog(Nefit),prob=prob,/double)
              N0_fit_aia_A(ifl) = exp(A[0]+A[1])
              lN_fit_aia_A(ifl) = 1./A[1]
             r2N_fit_aia_A(ifl) = r2

          endif                 ; covgflag = 'yes'
          
          stop
       endfor                   ; field lines loop.
    endif

    

    return
 end



pro test_coverage, rfit=rfit, covgflag=covgflag, $
                   aia = aia, euvia = euvia, euvib = euvib, eit = eit,$
                   mk4 = mk4, kcor = kcor, lascoc2 = lascoc2
  covgflag = 'no'
  if keyword_set(aia) or keyword_set(euvi) or keyword_set(eit) then begin
     R1=1.10
     R2=1.15
     R3=1.20
     R4=1.25
     if (where(rfit lt R1)               )(0) ne -1 AND
        (where(rfit gt R1 and rfit le R2))(0) ne -1 AND
        (where(rfit gt R2 and rfit le R3))(0) ne -1 AND 
        (where(rfit gt R3 and rfit le R4 )(0) ne -1 THEN covgflag = 'yes'
     endif  
  return
end
