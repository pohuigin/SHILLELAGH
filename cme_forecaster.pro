;+
;NAME: CME_FORECASTER
;
;USAGE:	
;	IDL> cme_forecaster,/byrne_08,/test,tbin=.5,trange=[-2,10]
;	IDL> cme_forecaster,'1-mar-2003',20.,500.,60.
;	IDL> cme_forecaster,/byrne_08,outtim=tim_arr,outrad=rad_arr
;INPUT:
;	cmetim= 	time when CME erupts
;	cmelon=		heliographic longitude of erupting CME
;	cmevel= 	(mean) velocity of CME propagation
;	cmewidth=	radial expansion width of CME
;	velwind=	solar wind velocity (does not affect CME currently)
;OUTPUT KEYWORDS:
;	outtim=		CME time points (in days since eruption)  
;	outrad=		CME distance from sun in AU (149.6x10^6 kilometers)
;INPUT KEYWORDS:
;	tbin= 		interval between points
;	trange=		days from time of CME eruption to model propagation
;	verbose=	print run messages for debugging
;	test=		require hit-return to progress propagation plot
;	plot_write=	write png plots of propagation
;	plot_path=	path to write plots to
;	orbit_path=	path to read in orbit files of spacecraft etc
;
;	Written: 08-NOV-2010 P.A. Higgins, P.T. Gallagher
;	Contact: pohuigin@gmail.com, peter.gallagher@tcd.ie
;-

pro cme_forecaster, cmetim, cmelon, cmevel, cmewidth, velwind, $ ;CME parameters
	outtim=outtim, outlon=outlon, outrad=outrad, $ ;CME output parameters
	timearr=cmetimarr, distarr=cmedistarr, widtharr=cmewidtharr, $ ;measured/modelled CME parameter arrays
	no_forecast=no_forecast,verbose=verbose, test=test, plot_write=plot_write, plot_path=plot_path, orbit_path=orbit_path, $
	trange=intrange, tbin=intbin, $ ;parker movie parameters
	byrne_08=byrne_08, stereo=stereo, themis=themis, omni=omni

	if not keyword_set(plot_path) then plotp='~/science/procedures/cme_propagation/plots/' else plotp=plot_path
	if not keyword_set(orbit_path) then orbitp='~/science/procedures/cme_propagation/orbits/' else orbitp=orbit_path

	if n_elements(intrange) lt 1 then trange=[-2,10] else trange=intrange
	if n_elements(intbin) lt 1 then tbin=1 else tbin=intbin
	
	rsun=6.955d5 ;in kilometers
	au_m=149.6d6 ;1 AU in kilometers
	l1_m=au_m-au_m*.01
	
	if n_elements(cmetim) lt 1 or n_elements(cmelon) lt 1 or n_elements(cmevel) lt 1 or n_elements(cmewidth) lt 1 then byrne_08=1.
	
	if keyword_set(byrne_08) then begin
		if n_elements(cmetim) lt 1 then cmetim='12-dec-2008 09:00'
		if n_elements(cmelon) lt 1 then cmelon=0.
		if n_elements(cmevel) lt 1 then cmevel=400.
		if n_elements(cmewidth) lt 1 then cmewidth=55.
		
		readcol,orbitp+'byrne_kins_meanspread.txt',d1,cme_r_arr,cme_top_a,d3,cme_bottom_a,d5,d6,d7,d8,d9,byrne_tim,form='F,F,F,F,F,F,F,F,F,F,A'
		cme_tim_arr=(anytim(byrne_tim)-cmetim)/3600./24.
		cme_r_arr=cme_r_arr*6.955d8/au_m
		if (where(cme_top_a gt 180))[0] ne -1 then cme_top_a[where(cme_top_a gt 180)]=cme_top_a[where(cme_top_a gt 180)]-360.
		if (where(cme_bottom_a gt 180))[0] ne -1 then cme_bottom_a[where(cme_bottom_a gt 180)]=cme_bottom_a[where(cme_bottom_a gt 180)]-360.
		cme_width_arr=cme_top_a-cme_bottom_a
		
		stop
		
	endif
	
	if keyword_set(themis) then begin
		readcol,orbitp+'THD_OR_SSC_6262.txt',thd_dd,thd_tt,thd_rr,dum1,thd_gsm_lon, dum2, dum3, dum4, thd_bow, thd_mag_pause,form='A,A,F,F,F,F,F,F,F',delim=' '
		thd_tim=anytim(strmid(thd_dd,6,4)+'-'+strmid(thd_dd,3,2)+'-'+strmid(thd_dd,0,2)+'T'+thd_tt)
		readcol,orbitp+'THE_OR_SSC_6262.txt',the_dd,the_tt,the_rr,dum1,the_gsm_lon, dum2, dum3, dum4, the_bow, the_mag_pause,form='A,A,F,F,F,F,F,F,F',delim=' '
		the_tim=anytim(strmid(the_dd,6,4)+'-'+strmid(the_dd,3,2)+'-'+strmid(the_dd,0,2)+'T'+the_tt)
	endif

	if keyword_set(omni) then begin
		readcol,orbitp+'OMNI_COHO1HR_MERGED_MAG_PLASMA_6262.txt',omn_dd,omn_tt,omn_hglon,omn_br,omn_bt,omn_bn,omn_vel,omn_rho,omn_temp,form='A,A,F,F,F,F,F,F,F',delim=' '
		omn_tim=anytim(strmid(omn_dd,6,4)+'-'+strmid(omn_dd,3,2)+'-'+strmid(omn_dd,0,2)+'T'+omn_tt)
		if (where(omn_hglon ge 180.))[0] ne -1 then omn_hglon[where(omn_hglon gt 180.)]=omn_hglon[where(omn_hglon gt 180.)]-360.
	endif
	
	if n_elements(velwind) lt 1 then velwind=400.

	window,0, xsize=700, ysize=700		
	!p.multi=0
	!p.background = 255
	!p.color = 0   
	!p.charthick=1
	!p.charsize = 2
	
	arr_cme_rad=0
	arr_cme_lon=cmelon
	
	if n_elements(cmetimarr) gt 1 and n_elements(cmewidtharr) gt 1 and n_elements(cmedistarr) gt 1 then begin
		
	endif else begin
		timarr=anytim(cmetim)+(findgen((trange[1]-trange[0])/tbin)*tbin+trange[0])*24.*3600.
		timvmsarr=anytim(timarr,/vms)
	endelse 

	if keyword_set(stereo) then begin
		sta_pos=GET_STEREO_LONLAT( timvmsarr, 'A', system = 'HEE', /degrees )
		stb_pos=GET_STEREO_LONLAT( timvmsarr, 'B', system = 'HEE', /degrees )
		;readcol,orbitp+'STA_COHO1HR_MERGED_MAG_PLASMA_6262.txt',sta_dd,sta_tt,sta_rad,sta_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
		;sta_tim=anytim(strmid(sta_dd,6,4)+'-'+strmid(sta_dd,3,2)+'-'+strmid(sta_dd,0,2)+'T'+stb_tt)
		;readcol,orbitp+'STB_COHO1HR_MERGED_MAG_PLASMA_6262.txt',stb_dd,stb_tt,stb_rad,stb_hglon,form='A,A,F,F,F,F,F,F,F',delim=' '
		;stb_tim=anytim(strmid(stb_dd,6,4)+'-'+strmid(stb_dd,3,2)+'-'+strmid(stb_dd,0,2)+'T'+stb_tt)
	endif

etaearth=anytim(cmetim)+au_m/cmevel
etal1=anytim(cmetim)+l1_m/cmevel
;etathe=anytim(cmetim)+(au_m+thd_rr)
;etathd
print,'CME ETA @ EARTH: '+anytim(etaearth,/vms)
print,'CME ETA @ L1: '+anytim(etal1,/vms)
print,' '
print,'CME properties are displayed in lower-left-hand corner of plot.'
if keyword_set(test) then begin & print,' ' & print,'Press ENTER/RETURN repeatedly.' & endif
if keyword_set(omni) then begin
	!p.multi=[0,1,2]
	plot,(omn_tim-min(omn_tim))/3600./24.,omn_bn,ytitle='OMNI Magnetic Field (B!dR!n,B!dT!n,B!dN!n)',/xsty,xran=([etal1-24.*3600.*5.,etal1+24.*3600.*5.]-min(omn_tim))/24./3600.
	oplot,(omn_tim-min(omn_tim))/3600./24.,omn_bt,color=5
	oplot,(omn_tim-min(omn_tim))/3600./24.,omn_br,color=4
	xyouts,.2,.65,'B!dN!n',color=0,/norm,chart=2 & xyouts,.25,.65,'B!dT!n',color=5,/norm,chart=2 & xyouts,.3,.65,'B!dR!n',color=4,/norm,chart=2
	vline,([etal1-24.*3600.,etal1,etal1+24.*3600.]-min(omn_tim))/24./3600.,color=3,lines=2
	plot,(omn_tim-min(omn_tim))/3600./24.,omn_hglon,/xsty,ytit='Heliographic Longitude', xtit='Days since 1-Dec-2008'
	hline,[cmelon-cmewidth/2.-10.,cmelon-cmewidth/2.,cmelon+cmewidth/2.,cmelon+cmewidth/2.+10.],lines=2,color=3
	vline,([etal1-24.*3600.,etal1,etal1+24.*3600.]-min(omn_tim))/24./3600.,color=3,lines=2
	wbound=where(omn_hglon ge cmelon-cmewidth/2. and omn_hglon le cmelon+cmewidth/2.)
	oplot,((omn_tim-min(omn_tim))/3600./24.)[wbound],omn_hglon[wbound],ps=4,color=3
stop
endif

if keyword_set(themis) then begin
	!p.multi=[0,1,2]
	plot,(thd_tim-min(thd_tim))/3600./24.,thd_bow,ytitle='Magnetosphere', xtit='Days',/xsty
	oplot,(thd_tim-min(thd_tim))/3600./24.,thd_mag_pause,color=5
	oplot,(the_tim-min(the_tim))/3600./24.,the_bow,color=0
	oplot,(the_tim-min(the_tim))/3600./24.,the_mag_pause,color=5
	
stop
endif
!p.multi=0

	for i=0, n_elements(timarr)-1 do begin
		;days=trange[0]
		;parkerstart=anytim(cmetim)+(ga+days)*84600.
		;parkerwest=cmelon+(ga-days)*14.
		;help,parkerwest,inner,anytim(parkerstart,/vms),cmevel, velwind
		if keyword_set(verbose) then print,'Time = '+anytim(timarr[i],/vms)
		thistim=timarr[i]
		
		if thistim lt anytim(cmetim) then $
			parkernplanets, anytim(thistim,/vms),/all, vel_wind = velwind, verbose=verbose, earthlon=earthlon $;,/satellites $ ;,foot_points = [str2arr(round(parkerwest))], inner = inner
		else begin
			cme_rad=(thistim-anytim(cmetim))*cmevel/au_m
			arr_cme_rad=[arr_cme_rad,cme_rad]
			arr_cme_lon=[arr_cme_lon,cmelon]
			parkernplanets, anytim(thistim,/vms), /all, vel_wind = velwind, verbose=verbose, $
				cme_rad=cme_rad, cme_lon=cmelon, cme_width=cmewidth, cme_scale=cme_scale, cme_color=cme_color, cme_tstart=cmetim, $
				arr_cme_rad=arr_cme_rad, arr_cme_lon=arr_cme_lon,earthlon=earthlon;,/satellites
		endelse
		;im=tvrd(true=1)
		
		xyouts,1.8,1.8,'CME '+cmetim+' vel='+strtrim(fix(round(cmevel)),2)+' lon='+strtrim(fix(round(cmelon)),2)+' width='+strtrim(fix(round(cmewidth)),2),/data,chars=1.4
		
		if keyword_set(stereo) then begin
			polrec, sta_pos[0,i]/au_m, sta_pos[1,i]+earthlon, sta_x, sta_y, /degrees
			oplot,[sta_x,sta_x], [sta_y,sta_y],ps=8,color=3
			xyouts,sta_x,sta_y,'STA',/data
			polrec, stb_pos[0,i]/au_m, stb_pos[1,i]+earthlon, stb_x, stb_y, /degrees
			oplot,[stb_x,stb_x], [stb_y,stb_y],ps=8,color=4
			xyouts,stb_x,stb_y,'STB',/data
			if keyword_set(verbose) then print,'STA HEE lon',sta_pos[1,i],' STB HEE lon',stb_pos[1,i]
		endif
		
		;gah=strtrim(arr2str(ga),2)
		
		if keyword_set(test) then begin & blank='' & read,blank & endif
		if keyword_set(plot_write) then window_capture,file=plotp+'CME_parker_propagate_'+string(i,form='(I03)'),/png

	endfor

	outtim=(timarr-anytim(cmetim))/3600./24.
	outrad=arr_cme_rad

if keyword_set(verbose) then print,'CME ETA @ EARTH: '+anytim(etaearth,/vms)
if keyword_set(verbose) then print,'CME ETA @ L1: '+anytim(etal1,/vms)

return

end
