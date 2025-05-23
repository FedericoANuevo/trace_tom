pro determine_fit_critical_parameters, radsamp=radsamp, rad_fit=rad_fit, $
                                       radfit_min=radfit_min, radfit_max=radfit_max, range_fit=range_fit
 ;Determine radsamp_max
  radsamp_max = max(radsamp)
 ;Determine the min and max rad over which we will actually evaluate the fit.
  radfit_min =                  min(rad_fit)
  radfit_max = min([radsamp_max,max(rad_fit)]) 
 ;Determine range of rad_fit_[instrument] over which we will actually evaluate the fit.
  range_fit = where(rad_fit ge radfit_min AND rad_fit le radfit_max)
  return
end
