pro planet_orbit, start_date, planet, orbit_x, orbit_y

  for i = 0, n_elements( planet ) - 1 do begin
 
    if ( planet eq 'mercury' ) then begin 
      planet_no = 1
      period = 0.2408 * 365.25
    endif
    if ( planet eq 'venus' ) then begin 
      planet_no = 2 
      period = 0.6152 * 365.25  
    endif
    if ( planet eq 'earth' ) then begin 
      planet_no = 3 
      period = 1.0 * 365.25 
    endif 
    if ( planet eq 'mars' ) then begin 
      planet_no = 4 
      period = 1.8809 * 365.25 
    endif  
    if ( planet eq 'jupiter' ) then begin 
      planet_no = 5 
      period = 11.862 * 365.25
    endif  
    if ( planet eq 'saturn' ) then begin 
      planet_no = 6 
      period = 29.458 * 365.25 
    endif 
    if ( planet eq 'uranus' ) then begin 
      planet_no = 7 
      period = 84.01 * 365.25
    endif  
    if ( planet eq 'neptune' ) then begin 
      planet_no = 8 
      period = 164.79 * 365.25  
    endif
    if ( planet eq 'pluto' ) then begin 
      planet_no = 9 
      period = 248.54 * 365.25 
    endif 
  endfor
  
  jd_struct = anytim2jd( start_date )
  jd = jd_struct.int + jd_struct.frac
    
  n_steps = sqrt( period )
    
  step_size = period / n_steps
  
  orbit_x = fltarr( n_steps + 1 ) 
  orbit_y = orbit_x
    
  for i = 0, n_steps do begin
    
    helio, jd, planet_no, rad, lon, lat
    lon = lon + 79. + 180.
    polrec, rad, lon, dx, dy, /degrees
    orbit_x[ i ] = dx
    orbit_y[ i ] = dy
    jd = jd + step_size
        
  endfor

end
