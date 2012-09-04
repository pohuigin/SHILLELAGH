;input heliocentric radius and longitude (corrected from vernal equinox and 180..)
;scale goes integers 0->50
;hel_lon is observed heliographic longitude of eruption
;tstart is cme lift off time
pro parkernplanets_blob, hel_rad, hel_lon, tstart, cme_scale=scale, cme_width=width, cme_color=color, verbose=verbose, $
	arr_cme_lon=hel_lon_arr, arr_cme_rad=hel_rad_arr;, earthlon=earthhlon

if n_elements(hel_lon) lt 1 then hel_lon=0.
if hel_rad le 0. then hel_rad=.05

; Convert CME lift off date to Julian date
st_jd_struct = anytim2jd( anytim(tstart,/vms) )
st_jd = st_jd_struct.int + st_jd_struct.frac
; Calculate earth's heliocentric latitude, longitude and distance from Sun
helio, st_jd, 3, earthhrad, earthhlon, earthhlat
vernal_equinox = 77.;79. ; the longitude of Capella in degrees 
earthhlon = earthhlon + vernal_equinox + 180.

;!!!!
;hel_lon=10.;-30.;-40.

polrec, hel_rad, hel_lon+earthhlon, pos_x, pos_y, /degrees

if keyword_set(verbose) then print,'posx',pos_x,'posy',pos_y

;Find rotation angle of CME plot symbol 
if pos_x lt 0 and pos_y ge 0 then theta=90.+tan((pos_x/pos_y))/!dtor
if pos_x lt 0 and pos_y lt 0 then theta=180.-tan((pos_x/pos_y))/!dtor
if pos_x ge 0 and pos_y lt 0 then theta=-90+tan((pos_x/pos_y))/!dtor
if pos_x ge 0 and pos_y ge 0 then theta=-90.+tan((pos_y/pos_x))/!dtor

;stop

if keyword_set(verbose) then print,'THETA=',theta
theta=theta

;if hel_rad le .15 then circle_sym, thick = 2, /fill else begin
	;Create CME plot symbol
	if n_elements(scale) lt 1 then scale=60
	n=220
	xyrcoord,[0,n,n],xx,yy,rr
	circimg=fltarr(n,n)
	;circimg[where(round(rr) eq round(scale))]=1.
	;circimg[where(round(rr) eq round(scale)-.1*scale)]=1.
	;circimg[n/2.-round(scale):n/2.-round(scale)+.1*scale,n/2.+1]=1.
	;circimg[n/2.+round(scale)-.1*scale:n/2.+round(scale),n/2.+1]=1.
	circimg[where(round(rr) ge (round(scale)-.2*scale) and round(rr) le round(scale))]=1.
	circimg[*,0:n/2]=0. & circimg=shift(circimg,0,-n/4.)
	circimg=rot(circimg,theta) ;180-
	w1=where(circimg eq 1) 
	xy1 = Find_Boundary(w1, XSize=n, YSize=n)
	wgood=round(findgen(49)/48.*n_elements(xy1[0,*]))
	x1=reform(xy1[0,*]) & y1=reform(xy1[1,*])
	x1=x1[wgood] & y1=y1[wgood]
	;x1=w1 mod n & y1=w1/n
	x1=(x1-n/2)/10. & y1=(y1-n/2)/10.
	
	;Plot CME
	if keyword_set(verbose) then print,'CME X=',pos_x,'Y=',pos_y
	usersym,x1,y1,/fill;thick=2
;endelse

if n_elements(color) lt 1 then color=6
;oplot, [pos_x,pos_x], [pos_y,pos_y], ps=8, color=color
plotsym,0,1,/fill

;Plot CME track
if n_elements(hel_lon_arr) gt 1 and n_elements(hel_rad_arr) gt 1 then begin
	polrec, hel_rad_arr, hel_lon_arr+earthhlon, pos_x_arr, pos_y_arr, /degrees
	oplot,pos_x_arr, pos_y_arr,ps=8,color=color
endif

oplot, [pos_x,pos_x], [pos_y,pos_y], ps=8, color=0

;Plot CME track
if n_elements(width) gt 0 then begin
	polrec, [hel_rad,hel_rad], [hel_lon+earthhlon+width/2.,hel_lon+earthhlon-width/2.], width_x, width_y, /degrees
	oplot,[0,width_x[0]], [0,width_y[0]],lines=2,color=color
	oplot,[0,width_x[1]], [0,width_y[1]],lines=2,color=color
	oplot,[width_x[0],pos_x,width_x[1]],[width_y[0],pos_y,width_y[1]],lines=2,color=color
endif

;Restore plotting parameters
set_line_color
circle_sym, thick = 2, /fill

end

;-------------------------------------------------------------------------->
; Programme to plot the Parker spiral and planets on a given date.

; Peter Gallagher (GSFC) - May 2004
; Paul Higgins (TCD) - Nov 2010 - added CME blob plotting 

pro parkernplanets_cme, date, foot_points = foot_points, planets = planets, $
                    all = all, inner = inner, outer = outer, $
		    vel_wind = invel_wind, fov = fov, perspective = perspective, $
		    cme_rad=cme_rad, cme_lon=cme_lon, cme_width=cme_width, cme_scale=cme_scale, cme_color=cme_color, cme_tstart=cme_tstart, $
		    arr_cme_rad=arr_cme_rad, arr_cme_lon=arr_cme_lon, _extra=_extra, verbose=verbose, $
		    earthlon=earthlon

  if ( n_elements( date ) eq 0 ) then get_utc, date, /vms
  
  if ( n_elements( foot_points ) eq 0 ) then foot_points = [ -90, -60, -30, 0, 30, 60, 90 ]

  if ( keyword_set( inner ) eq 0 ) then inner = 'inner'
  
  if ( keyword_set( vel_wind ) eq 0 ) then vel_wind = 500. else vel_wind = invel_wind ; km/s
  
  !p.background = 255
  
; Define planets to plot
  if keyword_set( all )   then planets = [ 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'neptune', 'pluto' ]
  if keyword_set( inner ) then planets = [ 'mercury', 'venus', 'earth', 'mars' ]
  if keyword_set( outer ) then planets = [ 'jupiter', 'saturn', 'neptune', 'pluto' ]

; Make sure planets string is lower case
  planets = strlowcase( planets )
     
; Convert time to VMS format   
  date = anytim( date, /vms )
  
; Assign a number to each planet
  planet_nos = intarr( n_elements( planets ) )
  
  for i = 0, n_elements( planets ) - 1 do begin
  
    if ( planets[ i ] eq 'mercury' ) then planet_nos[ i ] = 1 
    if ( planets[ i ] eq 'venus' )   then planet_nos[ i ] = 2 
    if ( planets[ i ] eq 'earth' )   then planet_nos[ i ] = 3 
    if ( planets[ i ] eq 'mars' )    then planet_nos[ i ] = 4 
    if ( planets[ i ] eq 'jupiter' ) then planet_nos[ i ] = 5 
    if ( planets[ i ] eq 'saturn' )  then planet_nos[ i ] = 6 
    if ( planets[ i ] eq 'uranus' )  then planet_nos[ i ] = 7 
    if ( planets[ i ] eq 'neptune' ) then planet_nos[ i ] = 8 
    if ( planets[ i ] eq 'pluto' )   then planet_nos[ i ] = 9 

  endfor

; Convert date to Julian date
  jd_struct = anytim2jd( date )
  jd = jd_struct.int + jd_struct.frac
  
; Calculate planetary heliocentric latitude, longitude and distance from Sun
  helio, jd, planet_nos, hel_rad, hel_lon, hel_lat
 
; Rotate coordinates relative to the Vernal Equinox
  vernal_equinox = 77.;79. ; the longitude of Capella in degrees 
  hel_lon = hel_lon - vernal_equinox; + 180.    
  
; Convert coordinates from polar to rectangular
  polrec, hel_rad, hel_lon, pos_x, pos_y, /degrees
  
; Output longitude of earth
  earthlon=hel_lon[2]
  
; Calculate the lat and long of the planet with perspective
  
  if keyword_set( perspective ) then begin
  
    if ( perspective eq 'mercury' ) then perspective = 1 
    if ( perspective eq 'venus' )   then perspective = 2 
    if ( perspective eq 'earth' )   then perspective = 3 
    if ( perspective eq 'mars' )    then perspective = 4 
    if ( perspective eq 'jupiter' ) then perspective = 5 
    if ( perspective eq 'saturn' )  then perspective = 6 
    if ( perspective eq 'uranus' )  then perspective = 7 
    if ( perspective eq 'neptune' ) then perspective = 8 
    if ( perspective eq 'pluto' )   then perspective = 9 

  endif else begin
  
    perspective = 3 ; Set default perspective to Earth
  
  endelse
  
  helio, jd, perspective, per_rad, per_lon, per_lat
  per_lon = per_lon + vernal_equinox + 180.
  
; Now plot planet positions 
  set_line_color
  circle_sym, thick = 2, /fill
  
; Calculate plot range
  max_orbit = ceil( max( hel_rad ) + 0.2 * max( hel_rad ) )
  
  if ( keyword_set( fov ) eq 0 ) then xrange = [  max_orbit, -max_orbit  ] else xrange = [ fov / 2., -fov / 2. ]
   
  plot, pos_x, pos_y, psym = 3, xrange = xrange, yrange = xrange, $
        xtitle = 'Heliocentric Distance (AU)', ytitle = 'Heliocentric Distance (AU)', $
	title = strmid( date, 0, 17 ) + ' (' + arr2str( fix(round(vel_wind)), /trim ) + ' km/sec)', /ys, /xs, color = 0, _extra=_extra

;Plot CME - phiggins
  if n_elements(cme_rad) gt 0 and n_elements(cme_lon) gt 0 and n_elements(cme_tstart) gt 0 then $
    parkernplanets_blob, cme_rad, cme_lon, cme_tstart, cme_width=cme_width, arr_cme_rad=arr_cme_rad, arr_cme_lon=arr_cme_lon, $
      cme_scale=cme_scale, cme_color=cme_color, verbose=verbose;, earthlon=earthlon
	  
; Plot planetary orbits and planet names
  
  for i = 0, n_elements( hel_rad ) - 1 do begin
      
    if ( planets[ i ] eq 'mercury' ) then color = 9
    if ( planets[ i ] eq 'venus' )   then color = 7
    if ( planets[ i ] eq 'earth' )   then color = 5
    if ( planets[ i ] eq 'mars' )    then color = 3
    if ( planets[ i ] eq 'jupiter' ) then color = 6
    if ( planets[ i ] eq 'saturn' )  then color = 2
    if ( planets[ i ] eq 'uranus' )  then color = 4
    if ( planets[ i ] eq 'neptune' ) then color = 10
    if ( planets[ i ] eq 'pluto' )   then color = 9
  
    planet_orbit, jd, planets[ i ], orbit_x, orbit_y
    oplot, orbit_x, orbit_y, psym = -3, color = 200
    
    plots, pos_x[ i ], pos_y[ i ], psym = 8, color = color, symsize = 2
    
    ;if ( hel_rad[ i ] gt 3 ) then begin
    xyouts, pos_x[ i ], pos_y[ i ] - 0.1, $
            strupcase( strmid( planets[ i ], 0 , 1 ) ) + strmid( planets[ i ], 1, 9 ), $
	    align = 0.5, /data, color = 0
    ;endif
  endfor
    
; Now plot the parker spiral (assume it is Archemedean)
  vel_wind = ( vel_wind / 150e6 ) * 24. * 60. * 60. ; AU per day

  rot_sun =  14.4 ; rotation rate of the Sun in degrees per day
  theta_spiral = -findgen( 360 * 6. )

  r_spiral = -( vel_wind / rot_sun ) * theta_spiral
  
  for j = 0, n_elements( foot_points ) - 1 do begin

    polrec, r_spiral, ( theta_spiral + per_lon ) + foot_points[ j ], spir_x, spir_y, /degrees
    oplot, spir_x, spir_y, color = 0, thick = 1

    if ( j eq 0 ) then index = where( r_spiral ge max_orbit * 0.65 )
    index = index[ 0 ]
   
    xyouts, spir_x[ index ], spir_y[ index ], string( foot_points[ j ] )+'!Uo', $
            /data, align = 0.5, color = 0, charsize = 1
  
  endfor

; Now plot a blue field line through Earth

    foot_point = 50
    polrec, r_spiral, ( theta_spiral + per_lon ) + foot_point, spir_x, spir_y, /degrees
    ;oplot, spir_x, spir_y, color = 5, thick = 2

    index = where( r_spiral ge max_orbit * 0.65 )
    index = index[ 0 ]
   
    ;xyouts, spir_x[ index ], spir_y[ index ], string( foot_point )+'!Uo', $
    ;       /data, align = 0.5, color = 0, charsize = 1
	    
; Now plot a blue field line through Mercury

  helio, jd, 1, per_rad, per_lon, per_lat
  per_lon = per_lon + vernal_equinox + 180.
 

    foot_point = 22
    polrec, r_spiral, ( theta_spiral + per_lon ) + foot_point, spir_x, spir_y, /degrees
    ;oplot, spir_x, spir_y, color = 9, thick = 2

    index = where( r_spiral ge max_orbit * 0.65 )
    index = index[ 0 ]
   
    ;xyouts, spir_x[ index ], spir_y[ index ], string( foot_point )+'!Uo', $
    ;        /data, align = 0.5, color = 0, charsize = 1
	    
  if keyword_set( satellites ) then begin
; Plot Cassini

  ;readcol, 'cassini_locations', cas_year, cas_doy, cas_rad, cas_lon, cas_lon, skip = 1
  
  restore,'cassini_locations.sav'
  
  cas_date = anytim( doy2utc( cas_doy, cas_year ) )
    
  dif_date = abs( cas_date - anytim( date ) )
  
  ind_date = where( dif_date eq min( dif_date ) )
     
  cas_rad = cas_rad[ ind_date ]
  cas_lon = cas_lon[ ind_date ]
  cas_rad = cas_rad[ 0 ]
  cas_lon = cas_lon[ 0 ] + 78. + 180.
   
  ;draw_circle, 0., 0., cas_rad, /data, color = 3, thick = 1, line = 2

  polrec, cas_rad, cas_lon, cas_x, cas_y, /deg
  plots, cas_x, cas_y, psym = 5, color = 0, symsize = 1
  xyouts, cas_x - 0.2, cas_y - 0.2, 'Cassini', $
          color = 0, charsize = 1

; Plot Galileo

  ;readcol, 'galileo_locations', gal_year, gal_doy, gal_rad, gal_lon, gal_lon, skip = 1
  
  ;save, f = 'galileo_locations.sav', gal_year, gal_doy, gal_rad, gal_lon, gal_lon
  
  restore,'galileo_locations.sav'
  
  gal_date = anytim( doy2utc( gal_doy, gal_year ) )
    
  dif_date = abs( gal_date - anytim( date ) )
  
  ind_date = where( dif_date eq min( dif_date ) )
     
  gal_rad = gal_rad[ ind_date ]
  gal_lon = gal_lon[ ind_date ]
  gal_rad = gal_rad[ 0 ]
  gal_lon = gal_lon[ 0 ] + 78. + 180.
   
  ;draw_circle, 0., 0., gal_rad, /data, color = 3, thick = 1, line = 2

  polrec, gal_rad, gal_lon, gal_x, gal_y, /deg
  plots, gal_x, gal_y, psym = 5, color = 0, symsize = 1
  xyouts, gal_x - 0.3, gal_y + 0.1, 'Galileo', $
          color = 0, charsize = 1

; Plot Ulysses

  ;readcol, 'ulysses_locations', uly_year, uly_doy, uly_rad, uly_lon, uly_lon, skip = 1
  
  ;save, f = 'ulysses_locations.sav', uly_year, uly_doy, uly_rad, uly_lon, uly_lon
  
  restore,'ulysses_locations.sav'
  
  uly_date = anytim( doy2utc( uly_doy, uly_year ) )
    
  dif_date = abs( uly_date - anytim( date ) )
  
  ind_date = where( dif_date eq min( dif_date ) )
     
  uly_rad = uly_rad[ ind_date ]
  uly_lon = uly_lon[ ind_date ]
  uly_rad = uly_rad[ 0 ]
  uly_lon = uly_lon[ 0 ] + 78. + 180.
   
  ;draw_circle, 0., 0., uly_rad, /data, color = 3, thick = 1, line = 2

  polrec, uly_rad, uly_lon, uly_x, uly_y, /deg
  plots, uly_x, uly_y, psym = 5, color = 0, symsize = 1
  xyouts, uly_x + 1., uly_y + 0.6, 'Ulysses', $
          color = 0, charsize = 1
  
;Plot STEREO
  thistim=anytim(date)
  readcol,orbitp+'STA_COHO1HR_MERGED_MAG_PLASMA_6262.txt',sta_dd,sta_tt,sta_rad,sta_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
  sta_tim=anytim(strmid(sta_dd,6,4)+'-'+strmid(sta_dd,3,2)+'-'+strmid(sta_dd,0,2)+'T'+sta_tt)
  readcol,orbitp+'STB_COHO1HR_MERGED_MAG_PLASMA_6262.txt',stb_dd,stb_tt,stb_rad,stb_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
  stb_tim=anytim(strmid(stb_dd,6,4)+'-'+strmid(stb_dd,3,2)+'-'+strmid(stb_dd,0,2)+'T'+stb_tt)
  wbest=(where(abs(sta_tim-thistim) eq min(abs(sta_tim-thistim))))[0]
  polrec, sta_rad[wbest], sta_hglon[wbest]+hel_lon[2], sta_pos_x, sta_pos_y, /degrees
  oplot,[sta_pos_x,sta_pos_x],[sta_pos_y,sta_pos_y],ps=8
  
  
  endif
  
; Finally plot   the Sun!
  plots, [ 0 , 0 ], [ 0, 0 ], psym = 8, color = 2, symsize = 3
  xyouts, 0, 0-0.1, 'Sun', align = 0.5, /data, color = 0


end
