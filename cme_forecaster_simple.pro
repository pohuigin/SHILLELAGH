;----------------------------------------------------------------------------->
;ETA = Time of arrival for CME to fit to.
;TOLERANCE = Time in seconds for closeness of best fit. 

pro cme_forecaster_fit, ineta, inseed, intolerance

eta=anytim(ineta)
tol=intolerance

if n_elements(inseed) lt 1 then thiscoeff=inseed else thiscoeff=90.
cme_forecaster_simple,/no_plot, eta=thiseta0, coefficient=thiscoeff, /debug
cme_forecaster_simple,/no_plot, eta=thiseta, coefficient=thiscoeff*2., /debug
lastcoeff=thiscoeff
thiscoeff=thiscoeff*2.
lasteta=thiseta0
tdiff=thiseta-eta
coeffarr=[lastcoeff,thiscoeff]
tdiffarr=[thiseta0-eta,tdiff]

while abs(tdiff) gt tol do begin
	if tdiff gt 0 and lasteta-eta gt 0 then newcoeff=thiscoeff/2.
	if tdiff lt 0 and lasteta-eta lt 0 then newcoeff=thiscoeff*2.
	if tdiff lt 0 and lasteta-eta gt 0 then newcoeff=abs(lastcoeff-thiscoeff)/2.
	if tdiff gt 0 and lasteta-eta lt 0 then newcoeff=abs(lastcoeff-thiscoeff)/2.

	cme_forecaster_simple,/no_plot, eta=neweta, coefficient=newcoeff, /debug

	tdiff=neweta-eta
	lastcoeff=thiscoeff
	thiscoeff=newcoeff

	print,'COEFF:',thiscoeff
	print,'TDIFF:',tdiff
	
	coeffarr=[coeffarr,thiscoeff]
	tdiffarr=[tdiffarr,tdiff]

endwhile

	print,'COEFF:',thiscoeff
	print,'TDIFF:',tdiff
	print,'ETA:',anytim(eta,/vms)
	print,'FIT ETA:',anytim(neweta,/vms)
	
	save,neweta,newcoeff,eta,coeffarr,tdiffarr,tol,inseed,file='~/science/papers/arg_cme_flare_20110607/dragmodelruns'+time2file(systim(/utc))+'.sav' 
stop

end

;----------------------------------------------------------------------------->
;TLINEAR = Make time-steps linear. Otherwise they will become more sparse as R^-4
;COEFFICIENT = Input coefficient for drag model instead of each variable.

pro cme_forecaster_simple,no_plot=no_plot, verbose=verbose, tlinear=tlinear, $
	eta=outeta, coefficient=incoeff, debug=debug, shillelagh=shillelagh, $
	write_plot=write_plot

;CHOOSE PATHS
;----------------------------------------------------------------------------->
plotp='~/science/data/cme_propagation/cme_propagate_runs/'
propp='~/science/data/cme_propagation/sw_prop_save/'
fileruns='~/science/data/cme_propagation/cme_propagate_runs/ireland_transformer_failure.txt'
;----------------------------------------------------------------------------->

;CHOOSE INPUTS
;----------------------------------------------------------------------------->
;CME parameters (june 7th event)
cmetim=anytim('4-aug-2011 03:41:00.000') ;'7-jun-2011 06:50:00.000'
cmelon=36. ;54. ;in degrees
cmevel=1000. ;1100. ;in km/s
cmewidth=60. ;in degrees
cme_mass=1d12 & txtmass='11' ;in Kg
cme_height=15. ;in R Sun
dragcoeff=1.

;Solar wind parameters
velwind=500.
rhowind=5.
;----------------------------------------------------------------------------->

;PHYSICAL CONSTANTS
;----------------------------------------------------------------------------->
;SW Properties fall off as 1/r^alpha
alpha_b=2.5d ;between 2 and 3
alpha_rho=4.d ;!!!NOT USED: switched to using empirical model in Chen 1996 (Eqn. 17.)
alpha_t=1.d

;Physical Constants
au_km=149.6d6 ;1 AU in kilometers
vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees or solar ascending node?? 
nan=0/0.
r_sun=6.955d5 ;km
l1_km=au_km-au_km*.01
mass_h_kg=1.672622d-27 ;mass of hydrogen in kilograms
mass_cme_kg=1d11 ;mass of CME in kilograms
cmcub_p_kmcub=(1d5/1d)^3. ;cm cubed per km cubed
omegasun=360d/(25.2d*3600d*24d) ;in degrees/second from diff. rot. of latitudes (-10 -> +10)
;----------------------------------------------------------------------------->

if n_elements(incoeff) eq 1 then coeff=incoeff

if keyword_set(shillelagh) then begin
	ffprop=file_search(propp+'sw_prop*.sav')
	ttprop=anytim(file2time(ffprop))
	;ttprop=anytim(strmid((file2time((strmid(ffprop,79,13)),/cc)),0,17)+'00');anytim(file2time((strmid(ffprop,79,13))))
endif

txtvel=strtrim(fix(cmevel),2)
txtwidth=strtrim(fix(cmewidth),2)
txtheight=strtrim(fix(cme_height),2)
txtswvel=strtrim(fix(velwind),2)
txtswrho=strtrim(fix(rhowind),2)

;Model parameters
trange=[0.,4.]
npt=1024. ;x,y dimensions for model box

if keyword_set(tlinear) then begin
	tbin=1./24.
	timarr=anytim(cmetim)+(findgen((trange[1]-trange[0])/tbin)*tbin+trange[0])*24.*3600.
	nstep=n_elements(timarr)
endif else begin
	tbin=1d-4
	nstep=((trange[1]-trange[0])/tbin)^(1./2.)+1. ;100.
	;Generate CME propagation arrays
	timarr=(findgen(nstep))^(2.)
	timarr=(timarr/(timarr[1]-timarr[0])*tbin+trange[0])*24.*3600.+anytim(cmetim)
	;timarr=(timarr/(max(timarr)-1)*(trange[1]-trange[0])+trange[0])*24.*3600.+anytim(cmetim)
endelse

timvmsarr=anytim(timarr,/vms)

;Generate 2d SW space
xyrcoord,[2,npt,npt],xx,yy & xx=xx-npt/2. & yy=yy-npt/2.
recpol, (xx/npt*4.*au_km), (yy/npt*4.*au_km), swradarr, swthtarr, /degrees
spacebin=2.*au_km/npt
swthtarr=reform(swthtarr,npt^2.) ;+ vernal_equinox + 180.
swradarr=reform(swradarr,npt^2.)

swvelarr=fltarr(npt^2.)+velwind

swrhoarr0=fltarr(npt^2.)+rhowind
constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]
swrhoarr=calc_helio_props(swradarr, fltarr(npt^2.)+au_km, swrhoarr0, /density, constants=constants_arr)
;swrhoarr_mann=calc_helio_props(swradarr, fltarr(npt^2.)+au_km, swrhoarr0, /density, constants=constants_arr,/mann)

;[r, theta, vel, rho]
sw_arrays=fltarr(npt^2.,4) & sw_arrays[*,0]=swradarr & sw_arrays[*,1]=swthtarr
sw_arrays[*,2]=swvelarr & sw_arrays[*,3]=swrhoarr

;Put initial CME height into Km/s
cme_height=cme_height*r_sun

setcolors,/sys,/quiet,/silent
if not keyword_set(no_plot) then begin
	window, xsize=700, ysize=700		
	!p.multi=0
;	!p.position=[.07,.05,.97,.95]
	!p.background = 255
	!p.color = 0   
	!p.thick=2.
	!p.charsize = 2
endif

;Make some test plots
if not keyword_set(no_plot) then begin
	plot,swradarr/au_km,swrhoarr,/ylog,ps=0,chars=2,ytit='Rho',xtit='Rad [Au]'
	;oplot,swradarr/au_km,swrhoarr_mann,ps=0,col=!red
	xyouts,.3,.8,'Chen 1996, Eq. 17',/norm
	xyouts,.3,.9,'Mann 1999, Eq. 10',/norm
	if keyword_set(write_plot) then window_capture,file=plotp+'shillelagh/density_profile',/png
;stop
endif

;Calculate for ballistic
etaearth=anytim(cmetim)+(au_km-cme_height)/cmevel
etal1=anytim(cmetim)+(l1_km-cme_height)/cmevel
;etathe=anytim(cmetim)+(au_m+thd_rr)
;etathd
if keyword_set(verbose) then begin
	print,'BALLISTIC'
	print,'CME ETA @ EARTH: '+anytim(etaearth,/vms)
	print,'CME ETA @ L1: '+anytim(etal1,/vms)
	print,' '
endif

;Calculate earth positions
earth_jd_struct = anytim2jd( anytim([cmetim,timarr],/vms) )
earth_jd = earth_jd_struct.int + earth_jd_struct.frac
helio, earth_jd, 3, earth_rad, earth_lon, earth_lat
earth_lon = earth_lon + vernal_equinox;+180.
earth_lon=sw_theta_shift(earth_lon)
cme_earthlon=earth_lon[0] & earth_lon=earth_lon[1:*]

;Model the CME propagation
cme_rad=cme_height

arr_cme_tim=timarr;fltarr(n_elements(timarr))
arr_cme_rad=fltarr(n_elements(timarr))+cme_rad
arr_cme_lon=fltarr(n_elements(timarr))
arr_cme_width=fltarr(n_elements(timarr))
arr_cme_area=fltarr(n_elements(timarr))
arr_cme_vel=fltarr(n_elements(timarr))+cmevel
arr_cme_acc=fltarr(n_elements(timarr))
arr_sw_vel=fltarr(n_elements(timarr))
arr_sw_rho=fltarr(n_elements(timarr))
arr_cme_mass=fltarr(n_elements(timarr))
nextacc=0.

for i=1, n_elements(timarr)-1 do begin
	if keyword_set(verbose) then print,'Time = '+anytim(timarr[i],/vms)
	thistim=timarr[i]
	arr_cme_tim[i]=thistim
	;deltat=double(((tbin*(i+1.))^2.)*3600.*24.);(thistim-anytim(cmetim))
	if i lt n_elements(timarr)-1 then deltat=double((timarr[i+1]-timarr[i]))

	if keyword_set(shillelagh) then begin
		wffbest=(where(abs(ttprop-thistim) eq min(abs(ttprop-thistim))))[0]
		restore,ffprop[wffbest];,/verb
;print,ffprop[wffbest]
;		stop
	endif


	;Calculate this earth position
	thisearth_lon=earth_lon[i]

	if keyword_set(debug) then $
		if 1./tbin mod i eq 0 then print,'thistim='+anytim(thistim,/vms)+' cmetim='+anytim(cmetim,/vms)
	
	if thistim ge cmetim then begin
		;initialize the acceleration to 0.
		thisacc=nextacc
		
		;Calculate HGI Longitude of CME from initial earth longitude
		cmetheta=cmelon+cme_earthlon

		;Read SW model data
		spirals=sw_arrays;[*,reverse(findgen((size(sw_arrays))[2]))]
		
		;Calculate a bounding arc that the CME passes through
		;swslicewidth=.0125
		maxslicewidth=0.1
		swslicewidth=(arr_cme_rad[i]/au_km)^2.*(maxslicewidth) > (spacebin*2.)/au_km
;	help,swslicewidth
;stop

		thisarcr12=[(arr_cme_rad[i-1]-swslicewidth*au_km) > 0, arr_cme_rad[i-1]+swslicewidth*au_km]
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
		thismass=cme_mass;1d12 ;in Kg
		;thismass=(1.43d5*arr_cme_rad[i-1]-5d11) < 1d12
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

		if n_elements(coeff) ne 1 then coeff=double(thisdrag)/double(thismass)*double(thiscmearea)
		
		if keyword_set(verbose) then print,'COEFF:',coeff
		
		nextacc=(-.5d)*double(coeff)*double(thisrhosw)*double(thiscmevel-thisvelsw)*double(abs(thiscmevel-thisvelsw))
				
		arr_cme_rad[i]=thiscmerad
		arr_cme_lon[i]=cmetheta
		arr_cme_width[i]=cmewidth
		arr_cme_area[i]=thiscmearea
		arr_cme_vel[i]=thiscmevel
		arr_cme_acc[i]=thisacc
		arr_sw_vel[i]=thisvelsw
		arr_sw_rho[i]=thisrhosw
		arr_cme_mass[i]=thismass
		
;		print,'NEXT ACCELERATION',nextacc
;		print,'DELTA VEL',thisacc*double(deltat)
;		help,thisdrag,thismass,thisrhosw,thiscmearea,thiscmevel,thisvelsw
	endif else begin
	;	cme_rad=deltat*cmevel;/au_m
	;	arr_cme_rad[i]=cme_rad
		arr_cme_lon[i]=cmelon
		arr_cme_width[i]=cmewidth
	;	arr_cme_vel[i]=thiscmevel
	endelse

;stop

goto_save_model_run:

	;Do plotting
	if not keyword_set(no_plot) then begin
;stop
		parkernplanets_cme,timvmsarr[i],/inner,vel=velwind
		
		if keyword_set(shillelagh) then oplot,spirals[*,0]/au_km,spirals[*,1],ps=3,color=!black,/polar
		
		setcolors,/sys,/sil,/quie
		if thistim ge cmetim then begin
			oplot,rarcarr/au_km,thetaarcarr*!dtor,/polar,lines=0,color=0
			oplot,rarcarr/au_km,thetaarcarr*!dtor,/polar,lines=2,color=255
			plotsym,0,1,/fill
			oplot,[arr_cme_rad[i],arr_cme_rad[i]]/au_km,[arr_cme_lon[i],arr_cme_lon[i]]*!dtor,ps=8,color=0
			oplot,[0,arr_cme_rad]/au_km,[arr_cme_lon[i]/au_km,arr_cme_lon[i]]*!dtor,/polar,lines=0,color=0
			oplot,[0,arr_cme_rad]/au_km,[arr_cme_lon[i],arr_cme_lon[i]]*!dtor,/polar,lines=2,color=255
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
;		xyouts,.15,.03,'CME '+anytim(cmetim,/vms)+' rad='+strtrim(fix(round(arr_cme_rad[i])),2)+' vel='+strtrim(fix(round(arr_cme_vel[i])),2)+' hglon='+strtrim(fix(round(arr_cme_lon[i])),2)+' width='+strtrim(fix(round(arr_cme_width[i])),2),/norm,chars=1.4,charthick=6,color=150	
		polyfill,[0.,1.,1.,0.],[0.,0.,.07,.07],color=150,/fill,/norm
		xyouts,.15,.03,'CME '+anytim(cmetim,/vms)+' rad='+strtrim(long(round(arr_cme_rad[i]/r_sun)),2)+' vel='+strtrim(fix(round(arr_cme_vel[i])),2)+' hglon='+strtrim(fix(round(arr_cme_lon[i])),2)+' width='+strtrim(fix(round(arr_cme_width[i])),2),/norm,chars=1.4,charthick=2,color=0
	
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
			if keyword_set(write_plot) then window_capture,file=plotp+'drag_shillelagh/CME_shillelagh_simp_propagate_'+string(i,form='(I03)'),/png
	endif

;wait,.2
;stop

endfor

wearth=(where(abs(arr_cme_rad-au_km) eq min(abs(arr_cme_rad-au_km))))[0]
wl1=(where(abs(arr_cme_rad-l1_km) eq min(abs(arr_cme_rad-l1_km))))[0]

outeta=arr_cme_tim[wearth]-(arr_cme_rad[wearth]-au_km)/arr_cme_vel[wearth]
outl1=arr_cme_tim[wl1]-(arr_cme_rad[wl1]-l1_km)/arr_cme_vel[wl1]

print,'CME ETA @ EARTH: '+anytim(outeta,/vms)
print,'CME ETA @ L1: '+anytim(outl1,/vms)
print,' '

;if not keyword_set(no_plot) then begin

	tseriespts=findgen(nstep)/(nstep-1.)*(trange[1]-trange[0])*3600.*24.
	constacc=abs(((reverse(arr_cme_vel))[0]-arr_cme_vel[0])/(max(tseriespts[5:*]-min(tseriespts[5:*]))))
	
	window, xsize=1000, ysize=700		
	!p.multi=[0,2,2]
	!p.charsize=1.4
	utplot,arr_cme_tim-min(arr_cme_tim),arr_sw_rho/mass_h_kg/cmcub_p_kmcub,min(arr_cme_tim),ytit='SW DENS [#/cm^3]',color=!forest,/ylog,yran=[1d-2,1d6],/ysty
	vline,cmetim-min(arr_cme_tim),color=!red,/vlog
	vline,arr_cme_tim[wearth]-min(arr_cme_tim),color=!blue,/vlog
	utplot,arr_cme_tim-min(arr_cme_tim),arr_cme_rad/r_sun,min(arr_cme_tim),ytit='CME RAD [rsun]',ps=4,/ylog
	oplot,tseriespts,(tseriespts*arr_cme_vel[0]+arr_cme_rad[0])/r_sun,color=!gray,lines=2
	oplot,tseriespts,(tseriespts*arr_cme_vel[0]+constacc*tseriespts+arr_cme_rad[0])/r_sun,color=!gray
	vline,cmetim-min(arr_cme_tim),color=!red,/vlog
	vline,arr_cme_tim[wearth]-min(arr_cme_tim),color=!blue,/vlog
	utplot,arr_cme_tim-min(arr_cme_tim),arr_cme_vel,min(arr_cme_tim),ytit='CME VEL [km/s]',ps=4,/ylog,yran=[200,2000],/ysty
	oplot,tseriespts[5:*],tseriespts[5:*]*((reverse(arr_cme_vel))[0]-arr_cme_vel[0])/(max(tseriespts[5:*]-min(tseriespts[5:*])))+arr_cme_vel[0],color=!gray
	hline,arr_cme_vel[0],color=!gray,/vlog
	vline,cmetim-min(arr_cme_tim),color=!red,/vlog
	vline,arr_cme_tim[wearth]-min(arr_cme_tim),color=!blue,/vlog
	oplot,arr_cme_tim-min(arr_cme_tim),arr_sw_vel,color=!forest
	utplot,arr_cme_tim-min(arr_cme_tim),abs(arr_cme_acc),min(arr_cme_tim),ytit='CME ACC [km/s^2]',ps=4,/ylog
	vline,cmetim-min(arr_cme_tim),color=!red,/vlog
	vline,arr_cme_tim[wearth]-min(arr_cme_tim),color=!blue,/vlog
	hline,constacc,color=!gray
	xyouts,.02,.47,'ETA @ Earth:'+anytim(arr_cme_tim[wearth],/vms)+', L1:'+anytim(arr_cme_tim[wl1],/vms)+' Param: '+txtmass+'kg '+txtvel+'kmps '+txtwidth+'deg '+txtheight+'rsun '+txtswvel+'kmps '+txtswrho+'rho',/norm
	!p.multi=0
	
	if keyword_set(write_plot) then window_capture,file=plotp+'cme_evolution_'+txtmass+'kg_'+txtvel+'kmps_'+txtwidth+'deg_'+txtheight+'rsun_'+txtswvel+'kmps_'+txtswrho+'rho',/png

	spawn,'echo '''+systim(/ut)+''' >> '+fileruns
	spawn,'echo '''+txtmass+'kg '+txtvel+'kmps '+txtwidth+'deg '+txtheight+'rsun '+txtswvel+'kmps '+txtswrho+'rho'+''' >> '+fileruns
	spawn,'echo ''CME ETA @ EARTH: '+anytim(arr_cme_tim[wearth],/vms)+''' >> '+fileruns
	spawn,'echo '' '' >> '+fileruns 
stop

;endif

stop

end