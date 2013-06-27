;print;+
; NAME:
;       sw_mann
; PURPOSE:
;       Solve a special version of Parkers solarwind solution (Mann et al 2005)
;		using numerical root finding with addiont of dipolar magnetic field (Mann et al 2003)
;		which yield model Alphane speed 
;
; EXPLANATION:
;		Uses Paker solar wind soliution at given tempeateure (which sets the sound speed and sonic radius)
; 		to give velocity as function of r, in conjuntion with a derived mass flux constaint at 1 au to 
;		give density as a function of r. Model magnetic field (diploe fixt at 1 au) to give alphane speed.
;		
;        
; CALLING SEQUENCE:
;       sw_mann, r, vs=v_sol, [/tog, /plot, n=n, va=va, vms=sms]
;
; INPUTS:
;		r an arry of distances to solve for, if not give will sole in range 1 - 250 Rsun
;		/tog set print plot
;		/plot generate a standard plot with range 1-250 Rsun
;
; OUTPUTS:
;		vs= paramter to contain return sw velocity (km/s)
;		n=  paramter to contain return sw density  (cm^-3)
;		va= paramter to contain alphane velcity    (km/s)
;       vms= paramter to contain (max) fast magnetosonic velocityr (km/s)
;
;
; AUTHOR:
;		Shane Maloney, March 15 2011
;-

pro sw_mann, toggle=toggle, plot=plot, r, vs=vs, n=n, va=va, vms=vms
	
	; Aussme temp of source of SW
	T=1.4d * (1e6)
	
	; Set constants
	; SI
	G = 6.671d * (1e-11) 		;Gravitantional
	kb = 1.3807d * (1e-23)		;Boltzmann
	mp = 1.6726d * (1e-27)		;proton mass
	ms = 1.99d * (1e30)			;solar mass
	rs = 6.958d * (1e8)			;solar radius
	
	mu = 0.6d					;mean molecular wight of solar wind
	
	; cgs
	C = 6.3d*10.0d^34.0d		; mass flux constant at 1 au
	
	; Geometry
	; Angle of prop 0 from equator 
	delta = 0.0d
	delta = delta*(!dpi/180.0d)
	
	
	; Derived qunatites
	vc = sqrt((kb*T)/(mu*mp))      ;sound speed  m/s
	rc = (G*ms) / (2.0d*(vc^2.0d)) ;sonic radius m
	
	; Loop
	if n_params() ne 1 then begin
		r = linspc(1.0, 250.0d, 10000.0)
	endif
	r=r*rs
	v = dblarr(n_elements(r))
	n = v
	va = v
	vms = v
	
	; Set up common block to pass extra parms to newton
	common parms, rhs

	for i=0., n_elements(r)-1 do begin
	
		br = 2.20d * ((rs/r[i])^3.0d) * sqrt(1.0 + 3.0d*sin(delta)^2.0d) ;cgs
		
		rhs = 4.0d*alog(r[i]/rc)+4*(rc/r[i]) - 3.0d
		v[i] = newton(r[i]/rc, 'sw') *vc   ; m/s
		print,v[i]
		n[i] = c/ (v[i] * r[i]^2 * 100.0^3.0) ; cm^-3
		va[i] = br / sqrt(4*!dpi*n[i]*(mp*1000.0))/ 100.0d ; m/s have to convert kg to g ->mp*1000.0d
		vms[i] = sqrt(vc^2+va[i]^2) ;m/s
	endfor
	
	if keyword_set(plot) then begin
	
		set_line_color
		;!p.multi = [0, 1, 2]
		!p.charsize=1.5
		!p.thick=2.0
		
		if keyword_set(toggle) then toggle, f='mann_model.eps', /eps, xs=8, ys=6
		
		;DEVICE, SET_FONT='Times', /TT_FONT
		
		
		; Plot the data, omit right and top axes:  
		PLOT, r/rs, v/1000.0, $  
			;TITLE = 'Coronal Model', $  
			;XTITLE = 'Distance (R'+SunSymbol()+')', $  
			YTITLE ='Velocity (km s!u-1!n)', $  
			;XSTYLE=9, YSTYLE=9, XMARGIN=[12, 8], YMARGIN=[4, 4], /ylog, /xlog, $
			xr=[1, 200], YRANGE =[1e-2, 1e3], pos=[0.175, 0.50, .90, 0.90], $
			xtickn=[' ',' ',' ' ,' '], xthick=3, ythick=3, charth=3
		; Draw the top x-axis, supplying labels, etc.  
		; Make the characters smaller so they will fit:  
		AXIS, XAXIS=1, XTITLE='Distance (AU)', xr=[(1.0d*rs)/(1.49e11), (250.0d*rs)/(1.49e11)], /xs, xthick=3, ythick=3, charth=3
		; Draw the right y-axis. Scale the current y-axis minimum  
		; values from Fahrenheit to Celsius and make them  
		; the new min and max values. Set YSTYLE=1 to make axis exact.  
		AXIS, YAXIS=1, yr=[1, 1e10], YSTYLE = 1, YTITLE = 'Density (cm!u-3!n)', /ylog, /save, xthick=3, ythick=3, charth=3
	
		oplot, r/rs, n, lines=2
		
		plot, r/rs, v/1000, XTITLE = 'Distance (R'+SunSymbol()+')', $
			XSTYLE=9, YSTYLE=1, XMARGIN=[12, 8], YMARGIN=[4, 4], /xlog, xthick=3, ythick=3, charth=3, $
			xr=[1, 100], yr=[1e-2, 800], YTITLE ='Velocity (km s!u-1!n)', pos=[0.175, 0.1, .90, 0.50]
		;AXIS, XAXIS=1, XTITLE='Distance (AU)', xr=[min(1*rs)/(1.49e11), max(100*rs)/(1.49e11)],  /xs
		;verline, rc/rs
		;horline, vc/1000.0
		oplot, r/rs, va/1000.0, lines=2
		oplot, r/rs, vms/1000.0, lines=3
		
		if keyword_set(toggle) then toggle

		endif
		r=r/rs
		vs=v/1000.0d
		va=va/1000.0d
		vms=vms/1000.0d
end


; Function necessary to use newtown. Is simplifed version of Parker sol with all terms involving 
; v taken to one side
function sw, x
	common parms
	y =  (x^2)-(2.0d*alog(x))-rhs
	return, y
end

	