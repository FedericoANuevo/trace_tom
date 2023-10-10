pro determindex,r0,th0,ph0,irad,ilat,ilon
  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common structure2,pfss_data
                                ; Purpose: Given the position vector
                                ; coordinates (r0,th0,ph0)                                                                                                                                                                                                                 
                                ;          find the 3 1D-indexes of
                                ;          the tomographic grid cell                                                                                                                                                                                                                 
                                ;          that contains that
                                ;          position.                                                                                                                                                                                                                                       
 drad = rad(1) - rad (0)
 dlat = lat(1) - lat (0)
 dlon = lon(1) - lon (0)
 rad0 = r0
 lat0 = 90 - th0/!dtor
 lon0 =      ph0/!dtor
 irad = (where(rad0 ge rad-drad/2 AND rad0 lt rad+drad/2))(0)
 ilat = (where(lat0 ge lat-dlat/2 AND lat0 lt lat+dlat/2))(0)
 ilon = (where(lon0 ge lon-dlon/2 AND lon0 lt lon+dlon/2))(0)
 return
end
