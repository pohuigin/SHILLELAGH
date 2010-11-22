;EXAMPLE: IDL> sw_plot_points,'23-nov-2008',/vel,foot=[-20,30],chvel=550.,/save

pro sw_plot_points, time, array, velocityplot=velocityplot, densityplot=densityplot, $
	temperatureplot=temperatureplot, magneticfieldplot=magneticfieldplot, planets=planets, $
	footpoints=footpoints, chvelocity=inchvelocity, save_plot=save_plot, $
	title_string=title_string, cbarpos=cbarpos, pmulti=pmulti, no_stereo=no_stereo;,multiplot=multiplot

thistime=anytim(time)

if n_elements(title_string) lt 1 then title_string=''

;SW Properties fall off as 1/r^alpha
alpha_b=2.5d ;between 2 and 3
alpha_rho=1.d
alpha_t=1.d
;Physical Constants
au_km=149.6d6 ;1 AU in kilometers
vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees 
nan=0/0.
r_sun=6.955d5
omegasun=360d/(25.2d*3600d*24d) ;in degrees/second from diff. rot. of latitudes (-10 -> +10)
constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]

;list save files and find nearest to TIME
savp=sw_paths(/savep)
plotp=sw_paths(/plotinterpheliop)
if n_elements(array) gt 1 then spirals=array else begin
	swff=file_search(savp+'sw_prop*.sav')
	timff=anytim(file2time(swff))
	wbest=where(abs(timff-thistime) eq min(abs(timff-thistime)))
	thistime=timff[wbest]
	swff=swff[wbest]
	restore,swff
	spirals=sw_arrays;[*,reverse(findgen((size(sw_arrays))[2]))]

print,swff
endelse


;find the earth's position
earth_jd_struct = anytim2jd( anytim(thistime,/vms) )
earth_jd = earth_jd_struct.int + earth_jd_struct.frac
helio, earth_jd, 3, earth_rad, earth_lon, earth_lat
earth_lon = earth_lon + vernal_equinox
earth_lon=sw_theta_shift(earth_lon)
this_earthlon=earth_lon[0]

wxs=700 & wys=700
window,xs=wxs,ys=wys

;Plot the axes
setcolors,/sys,/sil,/qui
if n_elements(pmulti) gt 0 then !p.multi=pmulti
plot,/polar,spirals[*,0],spirals[*,1]*!dtor,ps=3,xr=[2.*au_km,-2.*au_km]*1.2,yr=[2.*au_km,-2.*au_km]*1.2,/iso,title=strmid(anytim(thistime,/vms),0,15)+' '+title_string,/nodata,chars=1.4

;Set up color dot plots
loadct,5,/silent
;loadct,0,/silent
plotsym,0,1,/fill
if keyword_set(magneticfieldplot) then begin & brange=[0.,3.] & propnum=5 & endif; then brange=[0.,15.]
if keyword_set(velocityplot) then begin & velrange=[100.,700.] & propnum=2 & endif;[100.,900.]
if keyword_set(densityplot) then begin & densrange=[0.,5.] & propnum=3 & endif;[0.,15.]
if keyword_set(temperatureplot) then begin & temprange=[3,8] & propnum=4 & endif;[1d6,20d6]
if n_elements(propnum) eq 0 then begin & brange=[0.,3.] & propnum=5 & endif

;Plot spiral extrapolation
;oplot,[spirals[*,0],spirals[*,0]], [spirals[*,1],spirals[*,1]]*!dtor,ps=3
nspiral=n_elements(spirals[*,0])

case propnum of 
	2: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((spirals[m,propnum])-velrange[0])/(velrange[1]-velrange[0])*255.
	3: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>densrange[0] < densrange[1])-densrange[0])/(densrange[1]-densrange[0])*255. < 255.
	4: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>temprange[0] < temprange[1])-temprange[0])/(temprange[1]-temprange[0])*255. < 255.
	5: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>brange[0] < brange[1])-brange[0])/(brange[1]-brange[0])*255. < 255.
endcase



;Plot color dots (non-extrapolated, only back propagated)
;if keyword_set(no_alpha) then for k=0l,nblob-1l do oplot,/polar,[sw[k,0,i],sw[k,0,i]],[sw[k,1,i],sw[k,1,i]]*!dtor,ps=8,color=((sw[k,5,i])-brange[0])/(brange[1]-brange[0])*255. $
;	else for k=0l,nblob-1l do oplot,/polar,[sw[k,0,i],sw[k,0,i]],[sw[k,1,i],sw[k,1,i]]*!dtor,ps=8,color=(alog10(sw[k,5,i])>0.)/3.*255. < 255.
	
;Plot color bar
;if j eq 2 then if keyword_set(no_alpha) then color_table, [0.,10.],/right,/top,title='|B| [nT]' else 
if n_elements(cbarpos) lt 1 then begin
	cbarpos0=!p.clip/float([wxs,wys,wxs,wys]) & cbarpos=cbarpos0
	cbarpos[0]=cbarpos0[2]-(cbarpos0[2]-cbarpos0[0])*.5
	cbarpos[2]=cbarpos0[2];-(cbarpos0[2]-cbarpos0[0])*.05;[.6,.85,.95,.9]
	cbarpos[3]=cbarpos0[3];-(cbarpos0[3]-cbarpos0[1])*.05
	cbarpos[1]=cbarpos0[3]-(cbarpos0[3]-cbarpos0[1])*.05
endif
stop
case propnum of 
	2: color_table, velrange,cbarpos[[0,2]],cbarpos[[1,3]],title='Velocity [Km/s]',/shadow
	3: color_table, densrange,cbarpos[[0,2]],cbarpos[[1,3]],title='LOG |p| [cm^-3]',/shadow
	4: color_table, temprange,cbarpos[[0,2]],cbarpos[[1,3]],title='LOG |T| [K]',/shadow
	5: color_table, brange,cbarpos[[0,2]],cbarpos[[1,3]],title='LOG |B| [nT]',/shadow
endcase

if n_elements(pmulti) gt 0 then !p.multi=pmulti
plot,/polar,spirals[*,0],spirals[*,1]*!dtor,/nodata,/noerase,xr=[2.*au_km,-2.*au_km]*1.2,yr=[2.*au_km,-2.*au_km]*1.2,/iso,chars=1.4

;Save the Stereo A position
;		if j eq 1 then wbest1=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))

;Plot space craft
;		if j eq 2 then begin
;			setcolors,/sys,/quie,/sile
;			plot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,/nodata,/noerase,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso
;Plot the Sun
plotsym,0,2,/fill
oplot,/polar,[0,0]*au_km,[0,0],ps=8,color=0;!black
plotsym,0,1,/fill
oplot,/polar,[0,0]*au_km,[0,0],ps=8,color=255

if keyword_set(no_stereo) then goto,skip_stereo
;Plot Stereo B
stb_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'B', system = 'HCI', /degrees )
;wbest2=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
plotsym,0,1.5,/fill
oplot,/polar,[stb_pos[0],stb_pos[0]],[stb_pos[1],stb_pos[1]]*!dtor,ps=8,color=0
plotsym,0,1,/fill
oplot,/polar,[stb_pos[0],stb_pos[0]],[stb_pos[1],stb_pos[1]]*!dtor,ps=8,color=!white
oplot,/polar,[0,stb_pos[0]],[0,stb_pos[1]]*!dtor,lines=1,color=!white
;oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!green
;plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!black
;Plot Stereo A
sta_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'A', system = 'HCI', /degrees )
;wbesta=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
plotsym,0,1.5,/fill
oplot,/polar,[sta_pos[0],sta_pos[0]],[sta_pos[1],sta_pos[1]]*!dtor,ps=8,color=0
plotsym,0,1,/fill
oplot,/polar,[sta_pos[0],sta_pos[0]],[sta_pos[1],sta_pos[1]]*!dtor,ps=8,color=!white
oplot,/polar,[0,sta_pos[0]],[0,sta_pos[1]]*!dtor,lines=1,color=!white
;oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!red
;plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!black
skip_stereo:

;Plot Earth
polrec, earth_rad[0]*au_km, this_earthlon, earthpos_x, earthpos_y, /degrees
plotsym,0,1.5,/fill
oplot,[earthpos_x,earthpos_x],[earthpos_y,earthpos_y],ps=8,color=0
plotsym,0,1,/fill
oplot,[earthpos_x,earthpos_x],[earthpos_y,earthpos_y],ps=8,color=!white
oplot,[0,earthpos_x],[0,earthpos_y],lines=1,color=!white
;oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!yellow
;plotsym,0,2,/fill & oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!black

;Plot Mars
if keyword_set(planets) then begin
	jd_struct = anytim2jd(anytim(thistime,/vms))
	jd = jd_struct.int + jd_struct.frac
	helio, jd, [1,2,4,5,6], planet_rad, planet_lon, planet_lat
	planet_lon=planet_lon + vernal_equinox
	polrec, planet_rad*au_km, planet_lon, planet_x,planet_y, /deg
	setcolors,/sys,/sile,/quie
	planetcolors=[!orange,!magenta,!red,!purple,!yellow]
	for i=0,n_elements(planet_rad)-1 do begin
		plotsym,0,1.5,/fill
		oplot,[planet_x[i],planet_x[i]],[planet_y[i],planet_y[i]],color=0,ps=8
		plotsym,0,1.,/fill
		oplot,[planet_x[i],planet_x[i]],[planet_y[i],planet_y[i]],color=planetcolors[i],ps=8
	endfor
endif

;oplot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,ps=3		

;Over plot selected spirals 
if n_elements(footpoints) gt 0 then begin
	if n_elements(inchvelocity) lt 1 then chvelocity=700. else chvelocity=inchvelocity
	for i=0,n_elements(footpoints)-1 do begin
		thetafootarr=findgen(360)-360.
		
		;delta R = R(blob) - delta theta (velocity(blob) / omega_sun)
		rfootarr=-thetafootarr/omegasun*chvelocity
		
		;theta = theta(blob) + delta theta
		polrec, rfootarr, thetafootarr+footpoints[i]+this_earthlon, footpos_x, footpos_y, /degrees
		oplot,footpos_x,footpos_y,color=0
		oplot,footpos_x,footpos_y,lines=2,color=255
	endfor

endif

if keyword_set(save_plot) then window_capture,file=plotp+'sw_properties_type'+strtrim(propnum,2)+'_'+time2file(anytim(thistime,/vms)),/png






end