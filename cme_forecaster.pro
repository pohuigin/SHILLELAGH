;+
;NAME: CME_FORECASTER
;
;USAGE:	
;	IDL> cme_forecaster,'12-dec-2008 08:52',-5.,250.,50.,cmeheight=6.,cmemass=0.,trange=[-1./48.,6.],tbin=1./48.,/verb
;	IDL> cme_forecaster,/byrne_08,/test,tbin=.5,trange=[-2,10]
;	IDL> cme_forecaster,'1-mar-2003',20.,500.,60.
;	IDL> cme_forecaster,/byrne_08,outtim=tim_arr,outrad=rad_arr
;INPUT:
;	cmetim= 	time when CME erupts
;	cmelon=		heliographic longitude of erupting CME
;	cmevel= 	(mean) velocity of CME propagation
;	cmewidth=	radial expansion width of CME
;INPUT KEYWORDS:
;	cmemass=	a number between 1 and 3. 1 = 1E14g, 2 = 1E15g, 3 = 1E16g. Default is 1E15g.
;	cmeheight=	starting height for CME. default is 10 Rsun
;	coeffdrag=	a unitless constant dependent on the viscosity? CME surface texture? Default=1.
;	ballistic=	set to use ballistic CME model (no drag)
;	velwind=	solar wind velocity (does not affect CME currently). Only used in ballistic model.
;	tbin= 		interval between points
;	trange=		days from time of CME eruption to model propagation
;	verbose=	print run messages for debugging
;	test=		require hit-return to progress propagation plot
;	plot_write=	write png plots of propagation
;	plot_path=	path to write plots to
;	orbit_path=	path to read in orbit files of spacecraft etc
;OUTPUT KEYWORDS:
;	outtim=		CME time points (in days since eruption)  
;	outrad=		CME distance from sun in AU (149.6x10^6 kilometers)
;
;HISTORY:
;	Written - 08-NOV-2010 P.A. Higgins, P.T. Gallagher
;	Added drag model - 17-NOV-2010 P.A. Higgins
;
;	Contact: pohuigin@gmail.com, peter.gallagher@tcd.ie
;
;-

pro cme_forecaster, cmetim, cmelon, cmevel, cmewidth, $ ;CME parameters
	cmemass=cmemass, cmeheight=cmeheight, coeffdrag=coeffdrag, velwind=velwind, $
	res1tore=res1tore, verbose=verbose, test=test, $
	outtim=outtim, outlon=outlon, outrad=outrad, $ ;CME output parameters
	timearr=cmetimarr, distarr=cmedistarr, widtharr=cmewidtharr, $ ;measured/modelled CME parameter arrays
	ballistic=ballistic, no_forecast=no_forecast, plot_write=plot_write, plot_path=plot_path, orbit_path=orbit_path, $
	trange=intrange, tbin=intbin, $ ;parker movie parameters
	byrne_08=byrne_08, stereo=stereo, themis=themis, omni=omni

	if not keyword_set(plot_path) then plotp='~/science/procedures/cme_propagation/plots/' else plotp=plot_path
	if not keyword_set(orbit_path) then orbitp='~/science/procedures/cme_propagation/orbits/' else orbitp=orbit_path
	modelp='~/science/procedures/cme_propagation/sw_prop_save/'; else orbitp=orbit_path
	insitup='~/science/procedures/cme_propagation/sw_prop_insitu/'
	rootp='~/science/procedures/cme_propagation/'
	
	if n_elements(intrange) lt 1 then trange=[-2./12.,7.] else trange=intrange
	if n_elements(intbin) lt 1 then tbin=1./12. else tbin=intbin
	if n_elements(cmetim) gt 0 then cmetim=anytim(cmetim)
	
	rsun=6.955d5 ;in kilometers
	au_m=149.6d6 ;1 AU in kilometers
	l1_m=au_m-au_m*.01
	vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees 
	mass_h_kg=1.672622d-27 ;mass of hydrogen in kilograms
	mass_cme_kg=1d11 ;mass of CME in kilograms
	cmcub_p_kmcub=(1d5/1d)^3. ;cm cubed per km cubed
	if n_elements(cmemass) lt 1 then cme_mass=10. else cme_mass=10.^(cmemass) ;Mass of CME = CMEMASS x MASS_CME_KG
	print,'CME Mass = '+string(cme_mass*mass_cme_kg,form='(E12.2)')+' Kg'
	
	if n_elements(cmeheight) lt 1 then cme_height=10.*rsun else cme_height=cmeheight*rsun
	if n_elements(coeffdrag) lt 1 then dragcoeff=1. else dragcoeff=coeffdrag
	
	if n_elements(cmetim) lt 1 or n_elements(cmelon) lt 1 or n_elements(cmevel) lt 1 or n_elements(cmewidth) lt 1 then byrne_08=1.
	
	if keyword_set(byrne_08) then begin
		;if n_elements(cmetim) lt 1 then cmetim=anytim('12-dec-2008 09:00')
		;if n_elements(cmelon) lt 1 then cmelon=0.
		;if n_elements(cmevel) lt 1 then cmevel=400.
		;if n_elements(cmewidth) lt 1 then cmewidth=55.
		
		readcol,orbitp+'byrne_kins_meanspread.txt',d1,cme_height_arr,cme_top_a,d3,cme_bottom_a,d5,d6,d7,d8,d9,byrne_tim,form='F,F,F,F,F,F,F,F,F,F,A'
		cme_tim_arr=(anytim(byrne_tim))
		cme_r_arr=cme_height_arr*rsun;/au_m
		cme_vel_arr=deriv(cme_tim_arr,cme_r_arr)
		if (where(cme_top_a gt 180))[0] ne -1 then cme_top_a[where(cme_top_a gt 180)]=cme_top_a[where(cme_top_a gt 180)]-360.
		if (where(cme_bottom_a gt 180))[0] ne -1 then cme_bottom_a[where(cme_bottom_a gt 180)]=cme_bottom_a[where(cme_bottom_a gt 180)]-360.
		cme_width_arr=cme_top_a-cme_bottom_a
		
	endif
;	
;	if keyword_set(themis) then begin
;		readcol,orbitp+'THD_OR_SSC_6262.txt',thd_dd,thd_tt,thd_rr,dum1,thd_gsm_lon, dum2, dum3, dum4, thd_bow, thd_mag_pause,form='A,A,F,F,F,F,F,F,F',delim=' '
;		thd_tim=anytim(strmid(thd_dd,6,4)+'-'+strmid(thd_dd,3,2)+'-'+strmid(thd_dd,0,2)+'T'+thd_tt)
;		readcol,orbitp+'THE_OR_SSC_6262.txt',the_dd,the_tt,the_rr,dum1,the_gsm_lon, dum2, dum3, dum4, the_bow, the_mag_pause,form='A,A,F,F,F,F,F,F,F',delim=' '
;		the_tim=anytim(strmid(the_dd,6,4)+'-'+strmid(the_dd,3,2)+'-'+strmid(the_dd,0,2)+'T'+the_tt)
;	endif
;
;	if keyword_set(omni) then begin
;		readcol,orbitp+'OMNI_COHO1HR_MERGED_MAG_PLASMA_6262.txt',omn_dd,omn_tt,omn_hglon,omn_br,omn_bt,omn_bn,omn_vel,omn_rho,omn_temp,form='A,A,F,F,F,F,F,F,F',delim=' '
;		omn_tim=anytim(strmid(omn_dd,6,4)+'-'+strmid(omn_dd,3,2)+'-'+strmid(omn_dd,0,2)+'T'+omn_tt)
;		if (where(omn_hglon ge 180.))[0] ne -1 then omn_hglon[where(omn_hglon gt 180.)]=omn_hglon[where(omn_hglon gt 180.)]-360.
;	endif
	
	
	if n_elements(velwind) lt 1 then velwind=400.

	window, xsize=700, ysize=700		
	!p.multi=0
	!p.position=[.07,.05,.97,.95]
	!p.background = 255
	!p.color = 0   
	!p.charthick=1
	!p.charsize = 2
	
	;arr_cme_rad=0
	;arr_cme_lon=cmelon
	
	if n_elements(cmetimarr) gt 1 and n_elements(cmewidtharr) gt 1 and n_elements(cmedistarr) gt 1 then begin
		
	endif else begin
		timarr=anytim(cmetim)+(findgen((trange[1]-trange[0])/tbin)*tbin+trange[0])*24.*3600.
		timvmsarr=anytim(timarr,/vms)
	endelse 

;	if keyword_set(stereo) then begin
;		sta_pos=GET_STEREO_LONLAT( timvmsarr, 'A', system = 'HEE', /degrees )
;		stb_pos=GET_STEREO_LONLAT( timvmsarr, 'B', system = 'HEE', /degrees )
		;readcol,orbitp+'STA_COHO1HR_MERGED_MAG_PLASMA_6262.txt',sta_dd,sta_tt,sta_rad,sta_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
		;sta_tim=anytim(strmid(sta_dd,6,4)+'-'+strmid(sta_dd,3,2)+'-'+strmid(sta_dd,0,2)+'T'+stb_tt)
;		;readcol,orbitp+'STB_COHO1HR_MERGED_MAG_PLASMA_6262.txt',stb_dd,stb_tt,stb_rad,stb_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
		;stb_tim=anytim(strmid(stb_dd,6,4)+'-'+strmid(stb_dd,3,2)+'-'+strmid(stb_dd,0,2)+'T'+stb_tt)
;	endif

if keyword_set(ballistic) then begin
	etaearth=anytim(cmetim)+au_m/cmevel
	etal1=anytim(cmetim)+l1_m/cmevel
	;etathe=anytim(cmetim)+(au_m+thd_rr)
	;etathd
	print,'CME ETA @ EARTH: '+anytim(etaearth,/vms)
	print,'CME ETA @ L1: '+anytim(etal1,/vms)
	print,' '
endif
print,'CME properties are displayed in lower-left-hand corner of plot.'
if keyword_set(test) then begin & print,' ' & print,'Press ENTER/RETURN repeatedly.' & endif
;if keyword_set(omni) then begin
;	!p.multi=[0,1,2]
;	plot,(omn_tim-min(omn_tim))/3600./24.,omn_bn,ytitle='OMNI Magnetic Field (B!dR!n,B!dT!n,B!dN!n)',/xsty,xran=([etal1-24.*3600.*5.,etal1+24.*3600.*5.]-min(omn_tim))/24./3600.
;	oplot,(omn_tim-min(omn_tim))/3600./24.,omn_bt,color=5
;	oplot,(omn_tim-min(omn_tim))/3600./24.,omn_br,color=4
;	xyouts,.2,.65,'B!dN!n',color=0,/norm,chart=2 & xyouts,.25,.65,'B!dT!n',color=5,/norm,chart=2 & xyouts,.3,.65,'B!dR!n',color=4,/norm,chart=2
;	vline,([etal1-24.*3600.,etal1,etal1+24.*3600.]-min(omn_tim))/24./3600.,color=3,lines=2
;	plot,(omn_tim-min(omn_tim))/3600./24.,omn_hglon,/xsty,ytit='Heliographic Longitude', xtit='Days since 1-Dec-2008'
;	hline,[cmelon-cmewidth/2.-10.,cmelon-cmewidth/2.,cmelon+cmewidth/2.,cmelon+cmewidth/2.+10.],lines=2,color=3
;	vline,([etal1-24.*3600.,etal1,etal1+24.*3600.]-min(omn_tim))/24./3600.,color=3,lines=2
;	wbound=where(omn_hglon ge cmelon-cmewidth/2. and omn_hglon le cmelon+cmewidth/2.)
;	oplot,((omn_tim-min(omn_tim))/3600./24.)[wbound],omn_hglon[wbound],ps=4,color=3
;stop
;endif

;if keyword_set(themis) then begin
;	!p.multi=[0,1,2]
;	plot,(thd_tim-min(thd_tim))/3600./24.,thd_bow,ytitle='Magnetosphere', xtit='Days',/xsty
;	oplot,(thd_tim-min(thd_tim))/3600./24.,thd_mag_pause,color=5
;	oplot,(the_tim-min(the_tim))/3600./24.,the_bow,color=0
;	oplot,(the_tim-min(the_tim))/3600./24.,the_mag_pause,color=5
	
;stop
;endif
;!p.multi=0

;Calculate earth positions
earth_jd_struct = anytim2jd( anytim([cmetim,timarr],/vms) )
earth_jd = earth_jd_struct.int + earth_jd_struct.frac
helio, earth_jd, 3, earth_rad, earth_lon, earth_lat
earth_lon = earth_lon + vernal_equinox
earth_lon=sw_theta_shift(earth_lon)
cme_earthlon=earth_lon[0] & earth_lon=earth_lon[1:*]

;Model the CME propagation
	cme_rad=cme_height
	
	arr_cme_tim=fltarr(n_elements(timarr))
	arr_cme_rad=fltarr(n_elements(timarr))+cme_rad
	arr_cme_lon=fltarr(n_elements(timarr))
	arr_cme_width=fltarr(n_elements(timarr))
	arr_cme_area=fltarr(n_elements(timarr))
	arr_cme_vel=fltarr(n_elements(timarr))+cmevel
	arr_cme_acc=fltarr(n_elements(timarr))
	arr_sw_vel=fltarr(n_elements(timarr))
	arr_sw_rho=fltarr(n_elements(timarr))
	nextacc=0.
	for i=0, n_elements(timarr)-1 do begin
		if keyword_set(verbose) then print,'Time = '+anytim(timarr[i],/vms)
		thistim=timarr[i]
		arr_cme_tim[i]=thistim
		deltat=double(tbin*3600.*24.);(thistim-anytim(cmetim))

		;Calculate this earth position
		thisearth_lon=earth_lon[i]
		
		if cmelon gt 0 then modelcraft='stb'
		if cmelon eq 0 then modelcraft=['sta','stb']
		if cmelon lt 0 then modelcraft='sta'
		
		print,'thistim='+anytim(thistim,/vms)+' cmetim='+anytim(cmetim,/vms)
		
		if thistim ge cmetim then begin
			;initialize the acceleration to 0.
			if not keyword_set(ballistic) then begin
				thisacc=nextacc
				
				;Calculate HGI Longitude of CME from initial earth longitude
				cmetheta=cmelon+cme_earthlon
	
				;Read SW model data
				if not keyword_set(res1tore) then begin
					if i eq 0 then sw_propagate,anytim(thistim,/vms),spacecraft=modelcraft,trange=0.,tbin=0. $
						else sw_propagate,anytim(thistim,/vms),spacecraft=modelcraft,trange=0.,tbin=0.,/res1
				endif
				swff=file_search(modelp+'*.sav')
				timff=anytim(file2time(swff))
				wbest=where(abs(timff-thistim) eq min(abs(timff-thistim)))
				fftim=timff[wbest]
				swff=swff[wbest]
				restore,swff
				spirals=sw_arrays;[*,reverse(findgen((size(sw_arrays))[2]))]
				
				;Calculate a bounding arc that the CME passes through
				swslicewidth=.0125
				thisarcr12=[(arr_cme_rad[i-1]-swslicewidth*au_m) > 0, arr_cme_rad[i-1]+swslicewidth*au_m]
				thisarctheta12=[cmetheta-cmewidth/2., cmetheta+cmewidth/2.]
				if (where(spirals[*,0] ge thisarcr12[0] and spirals[*,0] le thisarcr12[1]))[0] eq -1 or (where(spirals[*,0] ge thisarctheta12[1] and spirals[*,1] le thisarctheta12[1]))[0] eq -1 then goto,goto_save_model_run
				thisarcprop=spirals[where(spirals[*,0] ge thisarcr12[0] and spirals[*,0] le thisarcr12[1]),*]
				thisarcprop=thisarcprop[where(thisarcprop[*,0] ge thisarctheta12[1] and thisarcprop[*,1] le thisarctheta12[1]),*]
				
				;Make arc for over plotting
				rarcarr=[fltarr(10)+thisarcr12[0],fltarr(10)+thisarcr12[1]] & rarcarr=[rarcarr,rarcarr[0]]
				thetaarcarr=findgen(10)/9.*(thisarctheta12[1]-thisarctheta12[0])+thisarctheta12[0] & thetaarcarr=[thetaarcarr,reverse(thetaarcarr),thetaarcarr[0]]
				
				;Average the properties within the bounding arc
				thisrhosw=mean(thisarcprop[*,3]*mass_h_kg*cmcub_p_kmcub) ;get from #/cm^3 to Kg/Km^3
				thisvelsw=mean(thisarcprop[*,2])  ;todo: need to translate to radial. only want velocity acting on cme...
				
				;Assume CME properties
				thismass=cme_mass*mass_cme_kg;1d15 ;in grams
				thisdrag=dragcoeff;1.
				
				thiscmevel=arr_cme_vel[i-1]+thisacc*double(deltat)
				thiscmerad=deltat*thiscmevel+arr_cme_rad[i-1]
				thiscmearea=2.*!pi*(thiscmerad)^2.*(1.-cos(cmewidth/2.*!dtor))
				
				if keyword_set(verbose) then begin
					print,'Acceleration = '+strtrim(thisacc,2)
					print,'CME Vel = '+strtrim(thiscmevel,2)
					print,'SW Vel = '+strtrim(thisvelsw,2)
					print,'Delta Vel = '+strtrim((thisacc*deltat),2)
				endif
				
	;			stop
				
				nextacc=(-.5d)*double(thisdrag)/double(thismass)*double(thisrhosw)*double(thiscmearea)*double(thiscmevel-thisvelsw)*double(abs(thiscmevel-thisvelsw))
				
				arr_cme_rad[i]=thiscmerad
				arr_cme_lon[i]=cmetheta
				arr_cme_width[i]=cmewidth
				arr_cme_area[i]=thiscmearea
				arr_cme_vel[i]=thiscmevel
				arr_cme_acc[i]=thisacc
				arr_sw_vel[i]=thisvelsw
				arr_sw_rho[i]=thisrhosw
			endif else begin
				cme_rad=deltat*cmevel;/au_m
				arr_cme_rad[i]=cme_rad
				arr_cme_lon[i]=cmelon
				arr_cme_width[i]=cmewidth
				arr_cme_vel[i]=thiscmevel
			endelse
		endif else begin
			if not keyword_set(res1tore) then begin
				if i eq 0 then sw_propagate,anytim(thistim,/vms),spacecraft=modelcraft,trange=0.,tbin=0. $
					else sw_propagate,anytim(thistim,/vms),spacecraft=modelcraft,trange=0.,tbin=0.,/res1
			endif
			swff=file_search(modelp+'*.sav')
			timff=anytim(file2time(swff))
			wbest=where(abs(timff-thistim) eq min(abs(timff-thistim)))
			fftim=timff[wbest]
			swff=swff[wbest]
			restore,swff
			spirals=sw_arrays
		endelse
	;	if keyword_set(plot_helio) then begin
			sw_plot_points, anytim(thistim,/vms), spirals, /magnetic, title_string='- '+strjoin(modelcraft)+' - B Field'
			setcolors,/sys,/sil,/quie
			if not keyword_set(ballistic) and thistim ge cmetim then begin
				oplot,rarcarr,thetaarcarr*!dtor,/polar,lines=0,color=0
				oplot,rarcarr,thetaarcarr*!dtor,/polar,lines=2,color=255
				plotsym,0,1,/fill
				oplot,[arr_cme_rad[i],arr_cme_rad[i]],[arr_cme_lon[i],arr_cme_lon[i]]*!dtor,ps=8,color=0
				oplot,[0,arr_cme_rad],[arr_cme_lon[i],arr_cme_lon[i]]*!dtor,/polar,lines=0,color=0
				oplot,[0,arr_cme_rad],[arr_cme_lon[i],arr_cme_lon[i]]*!dtor,/polar,lines=2,color=255
			endif
	;	endif
		
		;if thistim lt anytim(cmetim) then $
		;	parkernplanets, anytim(thistim,/vms),/all, vel_wind = velwind, verbose=verbose, earthlon=earthlon $;,/satellites $ ;,foot_points = [str2arr(round(parkerwest))], inner = inner
		;else begin
		;	parkernplanets, anytim(thistim,/vms), /all, vel_wind = velwind, verbose=verbose, $
		;		cme_rad=cme_rad, cme_lon=cmelon, cme_width=cmewidth, cme_scale=cme_scale, cme_color=cme_color, cme_tstart=cmetim, $
		;		arr_cme_rad=arr_cme_rad, arr_cme_lon=arr_cme_lon,earthlon=earthlon;,/satellites
		;endelse
		;im=tvrd(true=1)
		xyouts,.1,.07,'CME '+anytim(cmetim,/vms)+' rad='+strtrim(fix(round(arr_cme_rad[i])),2)+' vel='+strtrim(fix(round(arr_cme_vel[i])),2)+' hglon='+strtrim(fix(round(arr_cme_lon[i])),2)+' width='+strtrim(fix(round(arr_cme_width[i])),2),/norm,chars=1.4,charthick=3,color=0		
		xyouts,.1,.07,'CME '+anytim(cmetim,/vms)+' rad='+strtrim(fix(round(arr_cme_rad[i])),2)+' vel='+strtrim(fix(round(arr_cme_vel[i])),2)+' hglon='+strtrim(fix(round(arr_cme_lon[i])),2)+' width='+strtrim(fix(round(arr_cme_width[i])),2),/norm,chars=1.4,charthick=3,color=150
		
;		if keyword_set(stereo) then begin
;			polrec, sta_pos[0,i]/au_m, sta_pos[1,i]+earthlon, sta_x, sta_y, /degrees
;			oplot,[sta_x,sta_x], [sta_y,sta_y],ps=8,color=3
;			xyouts,sta_x,sta_y,'STA',/data
;			polrec, stb_pos[0,i]/au_m, stb_pos[1,i]+earthlon, stb_x, stb_y, /degrees
;			oplot,[stb_x,stb_x], [stb_y,stb_y],ps=8,color=4
;			xyouts,stb_x,stb_y,'STB',/data
;			if keyword_set(verbose) then print,'STA HEE lon',sta_pos[1,i],' STB HEE lon',stb_pos[1,i]
;		endif
		
		;gah=strtrim(arr2str(ga),2)
		
		if keyword_set(test) then begin & blank='' & read,blank & endif
	;	if keyword_set(plot_write) then begin
			if keyword_set(ballistic) then window_capture,file=plotp+'CME_parker_propagate_'+string(i,form='(I03)'),/png $
				else window_capture,file=plotp+'CME_shillelagh_propagate_'+string(i,form='(I03)'),/png
	;	endif
	endfor
	goto_save_model_run:
	if i ge n_elements(timarr)-1 then print,'MODEL RUN COMPLETED!' else print,'MODEL BOUNDS INSUFFICIENT FOR PROPAGATION!'
	save,arr_cme_rad,arr_cme_lon,arr_cme_width,arr_cme_area,arr_cme_vel,arr_cme_acc,arr_sw_vel,arr_sw_rho,file=rootp+'cme_model_drag_'+time2file(anytim(cmetim,/vms))+'hglon'+strtrim(round(cmetheta))
	
	plot,arr_cme_rad/rsun,arr_cme_vel,ps=4,xran=[0,70],/xsty,yran=[100,600],/ysty
	oplot,arr_cme_rad/rsun,arr_sw_vel,ps=4,color=!red
	if keyword_set(byrne08) then oplot,cme_height_arr,cme_vel_arr,ps=-4,color=!blue
	
stop

	outtim=(timarr-anytim(cmetim))/3600./24.
	outrad=arr_cme_rad

if keyword_set(verbose) then print,'CME ETA @ EARTH: '+anytim(etaearth,/vms)
if keyword_set(verbose) then print,'CME ETA @ L1: '+anytim(etal1,/vms)

return

end
