pro stream_insitu,time,tau,plot_data=plot_data

trange=[anytim(time),anytim(time)+tau]
tdata=[anytim(time)-7.*3600.*24.,anytim(time)+tau+7.*3600.*24.]

sw_omni=sw_get_data(tdata, /omni);, sta=sta, stb=stb, ace=ace, wind=wind, $
	;constants=constants_arr, geoindex=geoindex
w0=where(abs(tdata[0]-sw_omni[6,*]) eq min(abs(tdata[0]-sw_omni[6,*])))
w1=where(abs(tdata[1]-sw_omni[6,*]) eq min(abs(tdata[1]-sw_omni[6,*])))
sw_omni=sw_omni[*,w0:w1]

sw_geo=sw_get_data(tdata, /geoindex);, sta=sta, stb=stb, ace=ace, wind=wind, $
	;constants=constants_arr, geoindex=geoindex
wg0=(where(abs(tdata[0]-sw_geo[0,*]) eq min(abs(tdata[0]-sw_geo[0,*]))))[0]
wg1=(where(abs(tdata[1]-sw_geo[0,*]) eq min(abs(tdata[1]-sw_geo[0,*]))))[0]
sw_geo=sw_geo[*,wg0:wg1]

if keyword_set(plot_data) then begin
	window,xs=900,ys=700
	!p.color=0
	!p.background=255
	setcolors,/sys,/sile,/quie
	!p.multi=[0,2,3]
;velocity
	utplot,sw_omni[6,*]-anytim(time),sw_omni[2,*],anytim(time,/vms),ytit='Velocity', chars=2
	vline,trange-anytim(time),lines=2,color=!red
;density
	utplot,sw_omni[6,*]-anytim(time),sw_omni[3,*],anytim(time,/vms),ytit='Density', chars=2
	vline,trange-anytim(time),lines=2,color=!red
;magnetic field
	utplot,sw_omni[6,*]-anytim(time),sw_omni[7,*],anytim(time,/vms),ytit='Magnetic Field R(bk)-T(rd)-N(bl)', chars=2
	oplot,sw_omni[6,*]-anytim(time),sw_omni[8,*],color=!red
	oplot,sw_omni[6,*]-anytim(time),sw_omni[9,*],color=!blue
	vline,trange-anytim(time),lines=2,color=!red
;inclination is it asin or acos????
	smth=11 & brad=smooth(sw_omni[7,*],smth) & btan=smooth(sw_omni[8,*],smth)
	utplot,sw_omni[6,*]-anytim(time),asin(btan/((brad)^2.+(btan)^2.))/!dtor,anytim(time,/vms),ytit='Inclination (data=bk, model=rd)', chars=2,ps=4
	model_spiral_inclination,inclinations,velocity=sw_omni[2,*]
	oplot,sw_omni[6,*]-anytim(time),inclinations,color=!red,ps=4
	vline,trange-anytim(time),lines=2,color=!red

	utplot,sw_geo[0,*]-anytim(time),sw_geo[1,*],anytim(time,/vms), ytit='Kp', chars=2
	vline,trange-anytim(time),lines=2,color=!red
	utplot,sw_geo[0,*]-anytim(time),sw_geo[2,*],anytim(time,/vms), ytit='DST', chars=2,yran=[-30,30],/ysty
;over plot radial field
	oplot,sw_geo[0,*]-anytim(time),sw_omni[7,*]*3.,color=!red
	axis,yaxis=1,color=!red, ytit='B radial',yran=[-10,10]
	hline,0,lines=2
	vline,trange-anytim(time),lines=2,color=!red
;	utplot,sw_geo[0,*]-anytim(time),sw_geo[4,*],anytim(time,/vms), ytit='Auroral Electrojet', chars=2
;	oplot,sw_geo[0,*]-anytim(time),sw_geo[3,*],color=!gray,lines=1
;	oplot,sw_geo[0,*]-anytim(time),sw_geo[5,*],color=!gray,lines=2
;	oplot,sw_geo[0,*]-anytim(time),sw_geo[6,*],color=!gray,lines=3
;	vline,trange-anytim(time),lines=2,color=!red
;	utplot,sw_geo[0,*]-anytim(time),sw_geo[7,*],anytim(time,/vms), ytit='Polar Cap (Thule St.)', chars=2
	!p.multi=0

;!!!want to plot dst or ae between each hemisphere!! does (-) CH give activity in one hemisphere, and (+) in the other???

stop

endif

end

;-------------------------------------------------------------------------->
pro spiral_plot, earthlon, _extra=_extra

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

window,xs=700,ys=700
!p.multi=0
!p.color=0
!p.background=255
setcolors,/sys,/sile,/quie
plot,fltarr(359)+au_km,findgen(359.)*!dtor, /polar,xrange=[au_km*1.2,-au_km*1.2],yrange=[au_km*1.2,-au_km*1.2],/xsty,/ysty, _extra=_extra

polrec,au_km,earthlon,earth_x,earth_y,/degrees
oplot,[0,earth_x],[0,earth_y],lines=2,color=!gray
plotsym,0,2,/fill
oplot,[earth_x,earth_x],[earth_y,earth_y],ps=8,color=!blue

end

;-------------------------------------------------------------------------->
pro spiral_forecast, westfootpoint, earthlon, time, velocity, plot_spiral=plot_spiral, $
	outtime=outtime, outeta=deltat

if n_elements(earthlon) lt 1 then earthlon=0

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

polrec,au_km,earthlon,x_earth0,y_earth0,/degrees

make_spiral, rfootarr, thetafootarr, footpos_x, footpos_y, $
	foot=westfootpoint, earthlon=earthlon, velocity=velocity
spiral_earth_crossing, rfootarr, thetafootarr, footpos_x, footpos_y, $
	thetaearth, x_earth=x_earth1, y_earth=y_earth1

deltatheta=vangle([x_earth1,y_earth1,0.],[x_earth0,y_earth0,0.])/!dtor
;deltat=(thetaearth-earthlon)/(omegasun-omegaearth)
deltat=(deltatheta)/(omegasun-omegaearth)
print,'Time to arrival: '+strtrim(deltat/3600./24.,2)+' days'
print,'Time of arrival: '+anytim(anytim(time)+deltat,/vms)
outtime=anytim(time)+deltat
if deltat lt 0 then print,'Stream has already passed Earth.'

if keyword_set(plot_spiral) then begin
	spiral_plot, earthlon, title='Forecast'
	oplot,footpos_x,footpos_y,color=!red
	make_spiral, rfootarr1, thetafootarr1, footpos_x1, footpos_y1, $
		foot=westfootpoint+deltat*omegasun, earthlon=earthlon, velocity=velocity
	spiral_earth_crossing, rfootarr1, thetafootarr1, footpos_x1, footpos_y1, $
		thetaearth
	oplot,footpos_x1,footpos_y1,color=!red,lines=2
;calculated earth position at predicted time of earth crossing
	polrec,au_km,earthlon+deltat*omegaearth,earth_calcx1,earth_calcy1,/degrees
;calculated earth crossing from modeled stream
	polrec,au_km,thetaearth,earth_x,earth_y,/degrees
	plotsym,0,2
	oplot,[earth_x,earth_x],[earth_y,earth_y],ps=8,color=!red
	plotsym,0,1,/fill
	oplot,[earth_calcx1,earth_calcx1],[earth_calcy1,earth_calcy1],ps=8,color=!blue
endif

end

;-------------------------------------------------------------------------->
pro spiral_earth_crossing, rfootarr, thetafootarr, footpos_x, footpos_y, $
	thetaearth, x_earth=x_earth, y_earth=y_earth, $
	wbest=wbest, x12=x12, y12=y12

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

	wbest=where(abs(rfootarr-au_km) eq min(abs(rfootarr-au_km)))
	if rfootarr[wbest] le au_km then begin
		x1=footpos_x[wbest] & x2=footpos_x[wbest-1]
		y1=footpos_y[wbest] & y2=footpos_y[wbest-1]
		thetaearth=(thetafootarr[wbest]+thetafootarr[wbest-1])/2.
	endif else begin
		x2=footpos_x[wbest] & x1=footpos_x[wbest+1]
		y2=footpos_y[wbest] & y1=footpos_y[wbest+1]
		thetaearth=(thetafootarr[wbest]+thetafootarr[wbest+1])/2.
	endelse
	
	polrec, au_km, thetaearth, x_earth, y_earth, /degrees

	x12=[x1,x2]
	y12=[y1,y2]

end

;-------------------------------------------------------------------------->
;footpoints = 2 element array in HG lon
;velocities = 2 element array in km/s
pro stream_duration, tauearth, dthetaearth, footpoints=footpoints, velocities=velocities, $
	constants=constants, earthlon=earthlon, plot_spiral=plot_spiral

if n_elements(earthlon) lt 1 then earthlon=0.
earthlon=0.

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

dthetafoot=abs(footpoints[1]-footpoints[0])
taufoot=dthetafoot/(omegasun-omegaearth)

make_spiral, rfootarr0, thetafootarr0, footpos_x0, footpos_y0, $
	foot=footpoints[0], earthlon=earthlon, velocity=velocities[0], $
	bin=bin, constants=constants_arr
spiral_earth_crossing, rfootarr0, thetafootarr0, footpos_x0, footpos_y0, $
	thetaearth0, x_earth=x_earth0, y_earth=y_earth0;, $
;	wbest=wbest, x12=x12, y12=y12

make_spiral, rfootarr1, thetafootarr1, footpos_x1, footpos_y1, $
	foot=footpoints[1], earthlon=earthlon, velocity=velocities[1], $
	bin=bin, constants=constants_arr
spiral_earth_crossing, rfootarr1, thetafootarr1, footpos_x1, footpos_y1, $
	thetaearth1, x_earth=x_earth1, y_earth=y_earth1;, $
;	wbest=wbest, x12=x12, y12=y12

dthetaearth=(abs(thetaearth1-thetaearth0))[0]
tauearth=dthetaearth/(omegasun-omegaearth)

help,dthetafoot,dthetaearth
print,'tau foot=',taufoot/3600./24.,'tau earth=',tauearth/3600./24.

if keyword_set(plot_spiral) then begin
	spiral_plot, earthlon,title='Stream Duration'
	oplot,footpos_x1,footpos_y1,color=!red
	oplot,footpos_x0,footpos_y0,color=!blue
	oplot,[0,x_earth0],[0,y_earth0],color=!gray,lines=2
	oplot,[0,x_earth1],[0,y_earth1],color=!gray,lines=2
endif

end

;-------------------------------------------------------------------------->
pro make_spiral, rfootarr, thetafootarr, footpos_x, footpos_y, $
	foot=infootpoint, earthlon=inearthlon, velocity=inchvelocity, $
	bin=inbin, constants=constants_arr

if n_elements(inbin) lt 1 then bin=.1 else bin=inbin
footpoint=infootpoint
earthlon=inearthlon
chvelocity=inchvelocity

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

	thetafootarr=footpoint+earthlon+(-findgen(360./float(bin))*float(bin))
	
	;delta R = R(blob) - delta theta (velocity(blob) / omega_sun)
	rfootarr=abs((thetafootarr-thetafootarr[0])/omegasun*chvelocity)
	thetafootarr=sw_theta_shift(thetafootarr)
	
	;theta = theta(blob) + delta theta
	polrec, rfootarr, thetafootarr, footpos_x, footpos_y, /degrees
end

;-------------------------------------------------------------------------->
;Analyze the properties of archimedian parker spirals
;inclination
;earth crossing HSSWSs
;expansion factors
pro model_spiral_inclination,inclinations,test=test,bin=bin,velocity=velocity,width=width

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

earthlon=0.

if n_elements(bin) lt 1 then bin=.1
if n_elements(velocity) gt 0 then begin
	chvelocity=velocity
	n=n_elements(velocity)
endif else begin
	n=100.
	chvelocity=findgen(n)*20.+300.
endelse


footpoints=fltarr(n)
inclinations=fltarr(n)

for i=0,n-1 do begin	
	make_spiral, rfootarr, thetafootarr, footpos_x, footpos_y, $
		foot=footpoints[i], earthlon=earthlon, velocity=chvelocity[i], $
		bin=bin, constants=constants_arr

;	thetafootarr=findgen(360./float(bin))*float(bin)-360.+footpoints[i]+earthlon
	
	;delta R = R(blob) - delta theta (velocity(blob) / omega_sun)
;	rfootarr=-thetafootarr/omegasun*chvelocity[i]
	
	;theta = theta(blob) + delta theta
;	polrec, rfootarr, thetafootarr, footpos_x, footpos_y, /degrees

	spiral_earth_crossing, rfootarr, thetafootarr, footpos_x, footpos_y, $
		thetaearth, x_earth=x_earth, y_earth=y_earth, $
		wbest=wbest, x12=x12, y12=y12

;	wbest=where(abs(rfootarr-au_km) eq min(abs(rfootarr-au_km)))
;	if rfootarr[wbest] le au_km then begin
;		x1=footpos_x[wbest] & x2=footpos_x[wbest-1]
;		y1=footpos_y[wbest] & y2=footpos_y[wbest-1]
;		thetaearth=(thetafootarr[wbest]+thetafootarr[wbest-1])/2.
;	endif else begin
;		x2=footpos_x[wbest] & x1=footpos_x[wbest+1]
;		y2=footpos_y[wbest] & y1=footpos_y[wbest+1]
;		thetaearth=(thetafootarr[wbest]+thetafootarr[wbest+1])/2.
;	endelse	
;	polrec, au_km, thetaearth, x_earth, y_earth, /degrees

	incearth=atan(y_earth/x_earth)/!dtor
	incspiral=atan((y12[1]-y12[0])/(x12[1]-x12[0]))/!dtor
	;if incspiral gt 0 and incearth lt 0. then incspiral=90.-incspiral
;	inclinations[i]=incearth-incspiral
	inclinations[i]=vangle([(x12[1]-x12[0]),(y12[1]-y12[0]),0],-[x_earth,y_earth,0])/!dtor
	
	if keyword_set(test) then begin
		thisvel=(chvelocity[i])[0]
		thisinc=(inclinations[i])[0]
		incspiral=(incspiral)[0]
		incearth=(incearth)[0]

		help,thisvel,incearth,incspiral,thisinc
		window,xs=900,ys=450
		setcolors,/sys,/quie,/sile
		!p.multi=[0,2,1]
		plot,fltarr(359)+au_km,findgen(359.)*!dtor, /polar,xrange=[au_km*1.2,-au_km*1.2],yrange=[au_km*1.2,-au_km*1.2],/xsty,/ysty
		oplot,[0,x_earth],[0,y_earth],ps=-1,color=!blue
		oplot,footpos_x,footpos_y,color=0,thick=2,lines=2
;		oplot,footpos_x,footpos_y,ps=4,color=255
		mm=(y2-y1)/(x2-x1)
		bb=y1-x1*mm
		oplot,[-2.*au_km,x1,x2,2.*au_km],[-2.*au_km*mm+bb,y1,y2,2.*au_km*mm+bb],ps=-1,color=!red
yrangebin=bin*10.
yrange=[y1-4.*(y2-y1)/float(yrangebin),y2+4.*(y2-y1)/float(yrangebin)]
;xrangebin=bin/200.
xrange=[x1-4.*(y2-y1)/float(yrangebin),x2+4.*(y2-y1)/float(yrangebin)]
;[x1-4.*(x2-x1)/float(xrangebin),x2+4.*(x2-x1)/float(xrangebin)]
		plot,fltarr(359)+au_km,findgen(359.)*!dtor, /polar,xrange=xrange,yrange=yrange,/xsty,/ysty,/iso
		oplot,[0,x_earth],[0,y_earth],ps=-1,color=!blue
		oplot,[0,x_earth],[0,y_earth],ps=4,color=!blue
		oplot,footpos_x,footpos_y,color=0,thick=2,lines=2
		;oplot,footpos_x,footpos_y,ps=4,color=255
		oplot,[-2.*au_km,x1,x2,2.*au_km],[-2.*au_km*mm+bb,y1,y2,2.*au_km*mm+bb],ps=-1,color=!red
		oplot,[-2.*au_km,x1,x2,2.*au_km],[-2.*au_km*mm+bb,y1,y2,2.*au_km*mm+bb],ps=1,color=!blue
		oplot,[-2.*au_km,x1,x2,2.*au_km],[-2.*au_km*mm+bb,y1,y2,2.*au_km*mm+bb],ps=4,color=!blue
xyouts,.16,.12,'velocity='+strtrim(fix(round(thisvel)),2)+' inclination='+strtrim(fix(round(thisinc)),2),/norm,chars=1.4
wshow
stop
	endif
endfor

if keyword_set(test) then begin
	plot,chvelocity,inclinations,yran=[0,50],xran=[500,1500],ps=4,ytit='Inclination',xtit='Flow Velocity',chars=1.4
	window_capture,file=sw_paths(/root)+'spiral_inclination_vs_flow_velocity',/png
	wshow
endif

;if (where(inclinations lt 0))[0] ne -1 then inclinations[where(inclinations lt 0)]=inclinations[where(inclinations lt 0)]+180.

end

;-------------------------------------------------------------------------->
;Case studies:
;model_parker_spirals, '09-mar-1998 01:18',[20,5], [550,400]
;model_parker_spirals, '11-nov-1997 01:21',[20,-5], [550,400]
;model_parker_spirals, '31-jun-2006 23:10',[30,-10], [600,400]
pro model_parker_spirals, time, footpoints, velocities

if n_elements(time) lt 1 then begin
	time='22-Aug-2010 22:05:43'
;[LEADING, FOLLOWING]
	footpoints=[10,-20] ;from on disk coronal hole feature extraction
	velocities=[700,300] ;velocities of coronal hole boundaries - estimated from shillelagh?
endif
earthlon=0. ;use helios to get earth longitude at TIME

spiral_forecast, footpoints[0], earthlon, time, velocities[0], /plot_spiral, outt=ftime,outeta=teta

stop

stream_duration, tau, dthetaearth, footpoints=footpoints, velocities=velocities, $
	earthlon=earthlon, /plot_spiral

stop

stream_insitu,ftime,tau,/plot_data
window_capture,sw_paths(/chforecast)+'coronal_hole_forecast'+time2file(time),/png

stop

;What is the inclination of a spiral to the earth-sun line. (infinite velocity => 0)
model_spiral_inclination,/test,velocity=velocity

stop

end