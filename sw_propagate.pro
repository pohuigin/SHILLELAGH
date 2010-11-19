;-------------------------------------------------------------------------->
;SHILLElagh (Solar-Heliospheric Image Longitude Latitude Environment ... lagh)
;A semi-empirical 2D model of the solar wind propagation through the heliosphere. 
;The inputs are OMNI, Stereo A and B insitu ascii data downloaded from CDAWeb. 
;The model covereage is 1-jan-2007 to the present for STEREO+OMNI, and 1-jan-1996 
;to the present for OMNI alone.
;Time, Radius, HGI Longitude, Density, Temperature, and Magnetic Field.
;
;EXAMPLE:
;	1.) Plot the solar wind outflow applying a 1/r^alpha fall-off in B-field. 
;		Hit return key to progress movie.
;	IDL>sw_propagate, '12-Dec-2008', trange=[-5,10], tbin=.5, /test,/plot_points
;	
;	2.) Plot the solar wind outflow applying a 1/r^alpha fall-off in B-field.
;		Allow model to run and output images to directory.
;	IDL>sw_propagate,'12-dec-2008',/no_alph,/plot_points
;NOTES:
;	This algorithm must be run in IDL 32 bit mode due to its utilization of SPICE.
;
;HISTORY:
;	12-Nov-2010 - P.A.Higgins - Written
;
;-------------------------------------------------------------------------->
;Set up the paths to the in situ data files.
function sw_paths, insitup=insitup, plotp=plotp, savep=savep, plotinterpheliop=plotinterpheliop

result=''

if keyword_set(insitup) then result='~/science/procedures/cme_propagation/sw_prop_insitu/'
if keyword_set(plotp) then result='~/science/procedures/cme_propagation/sw_prop_plots/'
if keyword_set(savep) then result='~/science/procedures/cme_propagation/sw_prop_save/'
if keyword_set(plotinterpheliop) then result='~/science/procedures/cme_propagation/sw_prop_interp_plots/'

return, result

end

;-------------------------------------------------------------------------->
;Interpolate the data points to a grid. Choose between Cylindrical and Heliocentric
;Returns a stack of maps
;RBIN= in AU
;RRANGE= in AU
;THETARANGE= in degrees
;THETABIN= in degrees
function sw_interp, sw_array, tim, constants=constants_arr, rbin=inrbin, $
	rrange=inrrange, thetabin=inthetabin, thetarange=inthetarange, $
	cylindrical=cylindrical, heliocentric=heliocentric, $
	plot_interp=plot_interp,save_maps=save_maps

;constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]

if not keyword_set(rbin) then rbin=.005*au_km else rbin=inrbin*au_km
if not keyword_set(rrange) then rrange=2.*au_km else rrange=inrrange*au_km
if not keyword_set(thetabin) then thetabin=.5 else thetabin=inthetabin
if not keyword_set(thetarange) then thetarange=[0,360] else thetarange=inthetarange ;minmax(sw_array[where(sw_array[*,0] le sqrt(2.*rrange^2.)),1]) else thetarange=inthetarange

if keyword_set(cylindrical) then begin
;Do theta in x direction, r in y direction
	nxy=[(thetarange[1]-thetarange[0])/thetabin,rrange/rbin]
	sw_array=sw_array[where(sw_array[*,0] le sqrt(2.*rrange^2.)),*]
	yrcoord=sw_array[*,0]/rbin
	xrcoord=(sw_array[*,1]+thetarange[0])/thetabin
	xrcoord=round(xrcoord)
	yrcoord=round(yrcoord)
	map_stack=fltarr(nxy[0],nxy[1],4)
	miss=nan;-1d31
	
	for i=2,5 do begin
		triangulate,xrcoord,yrcoord,triangles
		thisproj=trigrid(xrcoord,yrcoord,sw_array[*,i],triangles,nx=nxy[0],ny=nxy[1],/quint)
		thisproj=smooth(thisproj,[5,5])
		
		thisdat=fltarr(nxy[0],nxy[1])+miss
		thisdat[xrcoord,yrcoord]=sw_array[*,i]
		
		window,1,xs=700,ys=700
		!p.multi=[0,1,2]
		plot_image,thisdat > 267 < 595
		plot_image,thisproj > 267 < 595
		!p.multi=0
		
		stop
		
		;gaparr=fltarr(nxy[0],nxy[1])
		;gaparr[xrcoord,yrcoord]=sw_array[*,i]
		;gaparr=smooth(gaparr,[5,5])
		;wgap=where(gaparr eq 0)
		
		thisdat=fltarr(nxy[0],nxy[1])+miss
		thisdat[xrcoord,yrcoord]=sw_array[*,i]
		
		;FILL_MISSING, thisdat, miss,1,/extrap
		;FILL_MISSING, thisdat, miss,2,/extrap
		
		thisdat[wgap]=nan
		
		plot_image,alog10(thisdat)
		stop
		thisproj=thisdat
		;[nx,ny,[Vel, ]]
		map_stack[*,*,i-2]=thisproj
	endfor
endif

;goto,skipheliocentric

if keyword_set(heliocentric) then begin
	nxy=[2.*rrange/rbin,2.*rrange/rbin]
	sw_array=sw_array[where(sw_array[*,0] le sqrt(2.*rrange^2.)),*]
	polrec, sw_array[*,0], sw_array[*,1], pos_x, pos_y, /degrees
	xrcoord=(pos_x+rrange)/rbin
	yrcoord=(pos_y+rrange)/rbin
	xrcoord=round(xrcoord)
	yrcoord=round(yrcoord)
	map_stack=fltarr(nxy[0],nxy[1],4)
	miss=nan;-1d31
	if keyword_set(plot_interp) then window,1,xs=700,ys=700
	
	!p.multi=[0,2,2]
	for i=2,5 do begin
		triangulate,xrcoord,yrcoord,triangles
		thisproj=trigrid(xrcoord,yrcoord,sw_array[*,i],triangles,nx=nxy[0],ny=nxy[1],/quint)
;[v,rho,temp,bmag]
		
		if keyword_set(plot_interp) then begin
			;!p.multi=[5-i,2,2]
			;plot_image,thisdat > 267 < 595
			if i eq 2 then plot_image,thisproj > 267 < 595,/noerase
			if i eq 4 then plot_image,alog10(thisproj > 0) < 10 > 0,/noerase
			if i eq 3 then plot_image,alog10(thisproj > 0) < 3 > 0,/noerase
			if i eq 5 then plot_image,alog10(thisproj > 0) < 3 > 0,/noerase
			window_capture,file=sw_paths(/plotinterpheliop)+'frame_'+strtrim(tim,2),/png
		endif
		

		;thisdat=fltarr(nxy[0],nxy[1])+miss
		;thisdat[xrcoord,yrcoord]=sw_array[*,i]
		;thisproj=warp_tri(xrcoord,yrcoord,xrcoord,yrcoord,thisdat,output_size=nxy,/extrap)
		;thisproj = interpu_2d(thisdat, xrcoord, xrcoord, where(finite(thisdat) ne 1) mod nxy[0], where(finite(thisdat) ne 1)/nxy[0], inp_unc=inp_unc, ux_new=ux_new, uy_new=uy_new, out_unc=out_unc)
		;FILL_MISSING, thisdat, miss,nxy
		;thisproj=thisdat
		;F = InterpolateHeliosphere(R, F3D, xcgrid=xcgrid, xlgrid=xlgrid, rrgrid=rrgrid, 	$
 		;[/fillbad, /degrees])
		;R = KRIG2D(sw_array[*,i], xrcoord, yrcoord,GS=[1.,1.], BOUNDS=[0,0,nxy[0]-1.,nxy[1]-1.],expon=[.25,0.])

		map_stack[*,*,i-2]=thisproj
	endfor
	!p.multi=0

	
endif

;skipheliocentric:

save,map_stack,file=sw_paths(/sav)+'sw_properties_heliospheric_'+time2file(anytim(tim,/vms))+'.sav'

return,map_stack

end

;-------------------------------------------------------------------------->
;Propagate an array of plasma blobs along their spiral using given velocity
function sw_spiral, inblob, theta_range=inthetarange, theta_bin=inthetabin, r_range=inrrange, $
	constants=constants_arr

;constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]

if n_elements(inthetarange) lt 1 then thetarange=[-22.,22.] else thetarange=inthetarange
if n_elements(inthetabin) lt 1 then thetabin=1. else thetabin=inthetabin
if n_elements(inrrange) lt 1 then rrange=[-3.5*au_km,10.*au_km] else rrange=inrrange

;Get rid of small gap due to non integer theta range and large bin.
thetarange=thetarange+[-thetabin/2.,2.*thetabin]

blob_arr=inblob
blob_arr=blob_arr[where(blob_arr[*,0] le rrange[1] and blob_arr[*,0] gt rrange[0]),*]

blank=fltarr(1,6)
spiral_arr=blank

nblob=n_elements(blob_arr[*,0])
for i=0l,nblob-1l do begin
	thisblob=blob_arr[i,*]
	thisthetarange=thetarange
	
	;Account for blobs which have been propagated to -distance (imaginary distance) 
	;Start them at where they should have been emitted from the Sun.
	bloblt0=0.
	if thisblob[0] lt 0 then begin
		rbloblt0=abs(thisblob[0])
		deltatrot=abs(thisblob[0])/thisblob[2]
		thisblob[1]=thisblob[1]-omegasun*deltatrot
		thisblob[0]=0.
		bloblt0=1.
		thisthetarange[1]=thetarange[1];+omegasun*deltatrot
		thisthetarange[0]=thetarange[0]+omegasun*deltatrot
	endif
	
	;Find number of points along the spiral. Scale to the arc length to  be covered 
	npts=(thisthetarange[1]-thisthetarange[0])/thetabin;*(thisblob[0]/au_km)
	if npts lt 1 then continue
	
	;Properties of each point along spiral
	thisspiral=fltarr(npts,6)
	
	;Load SW properties for each point along spiral
	thisspiral[*,2]=thisblob[2]
	thisspiral[*,3]=thisblob[3]
	thisspiral[*,4]=thisblob[4]
	thisspiral[*,5]=thisblob[5]
	
	;Delta theta from blob starting point
	theta_spiral=findgen(npts)*thetabin+thisthetarange[0]
	
	;delta R = R(blob) - delta theta (velocity(blob) / omega_sun)
	if bloblt0 then r_spiral = thisblob[0]-( thisblob[2] / omegasun ) * theta_spiral $
		else r_spiral = thisblob[0]-( thisblob[2] / omegasun ) * theta_spiral
	thisspiral[*,0]=r_spiral
	
	;theta = theta(blob) + delta theta
	thisspiral[*,1]=thisblob[1]+theta_spiral

;Assume 1/R^alpha dependence of Solar Wind properties

;Filter out negative R values
	if bloblt0 then begin
		if rbloblt0 le 0 then continue
		if (where(thisspiral[*,0] gt 0))[0] eq -1 then continue else begin
			wgoodspiral=where(thisspiral[*,0] gt 0)
			if wgoodspiral[0] ne -1 then thisspiral=thisspiral[wgoodspiral,*]
		endelse
	
;B Field
		thisspiral[*,5]=thisspiral[*,5]*(rbloblt0/thisspiral[*,0])^alpha_b
;Temperature
		thisspiral[*,4]=thisspiral[*,4]*(rbloblt0/thisspiral[*,0])^alpha_t
;Density
		thisspiral[*,3]=thisspiral[*,3]*(rbloblt0/thisspiral[*,0])^alpha_rho
		
	endif else begin
;B Field
		thisspiral[*,5]=thisspiral[*,5]*(thisblob[0]/thisspiral[*,0])^alpha_b
;Temperature
		thisspiral[*,4]=thisspiral[*,4]*(thisblob[0]/thisspiral[*,0])^alpha_t
;Density
		thisspiral[*,3]=thisspiral[*,3]*(thisblob[0]/thisspiral[*,0])^alpha_rho
	endelse

	;Filter spiral points with R lt 0
	wgt0=where(thisspiral[*,0] gt 0)
	if wgt0[0] ne -1 then thisspiral=thisspiral[wgt0,*]
	
	spiral_arr=[spiral_arr,thisspiral]
	
;if i eq 50 then stop

endfor

outspiral=spiral_arr[1:*,*]
return,outspiral

end

;-------------------------------------------------------------------------->
;Make sure that longitudes run from 0-360
function sw_theta_shift, inarr

outarr=inarr
if (where(inarr ge 360))[0] ne -1 then outarr[where(inarr ge 360)]=outarr[where(inarr ge 360.)]-360.
if (where(inarr lt 0))[0] ne -1 then outarr[where(inarr lt 0)]=outarr[where(inarr lt 0.)]+360.

return, outarr

end

;-------------------------------------------------------------------------->

function sw_get_data, trange, omni=omni, sta=sta, stb=stb, ace=ace, wind=wind, $
	constants=constants_arr

orbitp=sw_paths(/insitu)
;constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]

if keyword_set(omni) then begin
	readcol,orbitp+'omni_'+strmid(time2file(trange[0],/date),0,4)+'.txt',omn_dd,omn_tt,omn_hglat,omn_hglon,omn_br,omn_bt,omn_bn,omn_bmag,omn_vel,omn_elev,omn_azim,omn_rho,omn_temp,form='A,A,F,F,F,F,F,F,F,F,F,F,F',delim=' '
	omn_tim=anytim(strmid(omn_dd,6,4)+'-'+strmid(omn_dd,3,2)+'-'+strmid(omn_dd,0,2)+'T'+omn_tt)
	omn_hglon=sw_theta_shift(omn_hglon)
	sc_vel_arr=omn_vel
	sc_tim_arr=omn_tim
	sc_r_arr=fltarr(n_elements(sc_tim_arr))+au_km
	sc_hgtheta_arr=omn_hglon ;need to convert to HAE or w/e.
	sc_bmag=omn_bmag;omn_br;(omn_br^2.+omn_bt^2.+omn_bn^2.)^.5
;todo: radial field
	sc_rho=omn_rho
	sc_temp=omn_temp
endif

if keyword_set(sta) then begin
	readcol,orbitp+'sta_'+strmid(time2file(trange[0],/date),0,4)+'.txt',sta_dd,sta_tt,sta_r,sta_hglat,sta_hglon,sta_br,sta_bt,sta_bn,sta_b,sta_v,sta_sw_lat,sta_sw_lon,sta_rho,sta_t,form='A,A,F,F,F,F,F,F,F,F,F,F,F,F',delim=' '
	sc_tim_arr=anytim(strmid(sta_dd,6,4)+'-'+strmid(sta_dd,3,2)+'-'+strmid(sta_dd,0,2)+'T'+sta_tt)
	sc_hgtheta_arr=sw_theta_shift(sta_hglon)
	sc_vel_arr=sta_v
	sc_r_arr=sta_r*au_km
	sc_bmag=sta_b;sta_br
	sc_rho=sta_rho
	sc_temp=sta_t
endif

if keyword_set(stb) then begin
	readcol,orbitp+'stb_'+strmid(time2file(trange[0],/date),0,4)+'.txt',stb_dd,stb_tt,stb_r,stb_hglat,stb_hglon,stb_br,stb_bt,stb_bn,stb_b,stb_v,stb_sw_lat,stb_sw_lon,stb_rho,stb_t,form='A,A,F,F,F,F,F,F,F,F,F,F,F,F',delim=' '
	sc_tim_arr=anytim(strmid(stb_dd,6,4)+'-'+strmid(stb_dd,3,2)+'-'+strmid(stb_dd,0,2)+'T'+stb_tt)
	sc_hgtheta_arr=sw_theta_shift(stb_hglon)
	sc_vel_arr=stb_v
	sc_r_arr=stb_r*au_km
	sc_bmag=stb_b;stb_br
	sc_rho=stb_rho
	sc_temp=stb_t
endif

;make an array of S=([r,theta,v,rho,temp,bmag,t],nblob)
sc_arr=fltarr(7,n_elements(sc_vel_arr))
sc_arr[0,*]=sc_r_arr
sc_arr[1,*]=sc_hgtheta_arr
sc_arr[2,*]=sc_vel_arr
sc_arr[3,*]=sc_rho
sc_arr[4,*]=sc_temp
sc_arr[5,*]=sc_bmag
sc_arr[6,*]=sc_tim_arr

;Check for bad data points
if (where(sc_arr[3,*] eq -1.*10^31.))[0] ne -1 then sc_arr[*,where(sc_arr[3,*] eq -1.*10^31.)]=nan
sc_arr=sc_arr[*,where(finite(sc_arr[3,*]) eq 1)]

return,sc_arr

end

;-------------------------------------------------------------------------->
;parker model solar wind speed vs. r from sun=> r-r_o = (v_sw/omega_sun)*(theta-theta_o)
;												v_sw=omega_sun*(r-r_o)/(theta-theta_o)
;intime = 	'DD-Mon-YYYY' to propagate to
;trange = 	[2 elements] days relative to input time
;tbin	=	days between extrapolation solutions

pro sw_propagate, intime, trange=intrange, tbin=intbin, no_alpha=no_alpha, full360=full360, $
	test=test, plot_points=plot_points, $
	interp_maps=interp_maps, save_maps=save_maps, plot_interp=plot_interp

time=anytim(intime)
if n_elements(intrange) lt 1 then trange=[time-5.*24.*3600.,time+10.*24.*3600.] else trange=time-intrange*24.*3600.
if n_elements(intbin) lt 1 then tbin=12.*3600. else tbin=intbin*24.*3600.

window,xs=750,ys=750

plotp=sw_paths(/plotp) ;'~/science/procedures/cme_propagation/sw_prop_movie/'

;SW Properties fall off as 1/r^alpha
alpha_b=2.5d ;between 2 and 3
alpha_rho=1.d
alpha_t=1.d
if keyword_set(no_alpha) then begin & alpha_b=0d & alpha_rho=0d & alpha_t=0d & endif

;Physical Constants
au_km=149.6d6 ;1 AU in kilometers
vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees 
nan=0/0.
r_sun=6.955d5
omegasun=360d/(25.2d*3600d*24d) ;in degrees/second from diff. rot. of latitudes (-10 -> +10)

constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]

;Read in situ data SC=([r,theta,v,rho,temp,bmag,t],nblob)
nspacecraft=3
sc_arr0=sw_get_data(trange,/omni,const=constants_arr)
sc_arr1=sw_get_data(trange,/sta,const=constants_arr)
sc_arr2=sw_get_data(trange,/stb,const=constants_arr)

ntime=round((trange[1]-trange[0])/tbin)
time_arr=findgen(ntime)*tbin+trange[0]

;Create SW propagation SW(r,theta,v,rho,temp,bmag,t_sc) array sw=[nblob,nparam=7,ntime]
nblob0=n_elements(sc_arr0[0,*])
sw0=fltarr(nblob0,7,ntime)
nblob1=n_elements(sc_arr1[0,*])
sw1=fltarr(nblob1,7,ntime)
nblob2=n_elements(sc_arr2[0,*])
sw2=fltarr(nblob2,7,ntime)

;Find range to solve spiral for each space craft: 1/2*Difference between earth and spacecraft + 90 degrees
wbestearth=(where(abs(time_arr[ntime/2.]-sc_arr0[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr0[6,*]))))[0]
;wbeststa=(where(abs(time_arr[ntime/2.]-sc_arr1[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr1[6,*]))))[0]
;wbeststb=(where(abs(time_arr[ntime/2.]-sc_arr2[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr2[6,*]))))[0]
stbhglon=(GET_STEREO_LONLAT( anytim(time_arr[ntime/2.],/vms), 'B', system = 'HEE', /degrees ))[1]
stahglon=(GET_STEREO_LONLAT( anytim(time_arr[ntime/2.],/vms), 'A', system = 'HEE', /degrees ))[1]

earthsta=stahglon;abs(sc_arr0[1,wbestearth]-sc_arr1[1,wbeststa])/2.
earthstb=stbhglon;abs(sc_arr0[1,wbestearth]-sc_arr2[1,wbeststb])/2.
if keyword_set(full360) then begin
	thetarange=ceil([[0.,360.], $
		[0,360.], $
		[0,360.]])
endif else begin
	thetarange=ceil([[earthstb/2., earthsta/2.], $
		[-earthsta/2.,180.-earthsta/2.], $
		[-1.*(180.+earthstb/2.),-1.*earthstb/2.]])
endelse

spiral_arrays=transpose([0.,0.,0.,0.,0.,0.])

;Run through each propagation time within TIMERANGE using TBIN
for i=0l,ntime-1l do begin
;Run through each spacecraft
	for j=0,nspacecraft-1 do begin
		ex=execute('sw=sw'+strtrim(j,2)+' & nblob=nblob'+strtrim(j,2)+' & sc_arr=sc_arr'+strtrim(j,2))
		thistime=time_arr[i]
;Fill time when blob measured is at spacecraft
		sw[*,6,i]=sc_arr[6,*]
		
;Fill velocity
		sw[*,2,i]=sc_arr[2,*]
		
;Fill B field
		sw[*,5,i]=sc_arr[5,*]
		
		
;Fill density
		sw[*,3,i]=sc_arr[3,*]
		
;Fill temperature
		sw[*,4,i]=sc_arr[4,*]
		
;Calculate radii
		sw[*,0,i]=sc_arr[0,*]-(sc_arr[6,*]-thistime)*sc_arr[2,*]

;Assume 1/R^alpha dependence of Solar Wind properties
sw=double(sw)
;B Field
		sw[*,5,i]=sw[*,5,i]*(au_km/abs(sw[*,0,i]))^alpha_b
;Temperature
		sw[*,4,i]=sw[*,4,i]*(au_km/abs(sw[*,0,i]))^alpha_t
;Density
		sw[*,3,i]=sw[*,3,i]*(au_km/abs(sw[*,0,i]))^alpha_rho

;Calculate earth longitude array
		if j eq 0 then begin
			earth_jd_struct = anytim2jd( anytim([thistime,reform(sw[*,6,i])],/vms) )
			earth_jd = earth_jd_struct.int + earth_jd_struct.frac
			helio, earth_jd, 3, earth_rad, earth_lon, earth_lat
			earth_lon = earth_lon + vernal_equinox
			earth_lon=sw_theta_shift(earth_lon)
			this_earthlon=earth_lon[0] & earth_lon=earth_lon[1:*]
		endif
		
;Calculate theta
		sw[*,1,i]=sc_arr[1,*] ;earth_lon
		sw[*,1,i]=sw_theta_shift(sw[*,1,i])
		
;Calculate spirals for each blob
		spirals=sw_spiral(sw[*,*,i],const=constants_arr, theta_range=thetarange[*,j])
		spiral_arrays=[spiral_arrays,spirals]
		
;Get rid of meaningless data  (R < R_sun)
		wbad=where(sw[*,0,i] lt r_sun)
		if wbad[0] ne -1 then sw[wbad,0,i]=nan
		
if not keyword_set(plot_points) then goto,skip_plot_points

;Plot the Sun
		if j eq 0 then begin
			setcolors,/sys,/sil,/qui
			plot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,ps=3,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso,title=anytim(thistime,/vms)
			plotsym,0,4,/fill
			oplot,/polar,[0,0]*au_km,[0,0],ps=8,color=!orange
		endif else oplot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,ps=3
		
		;Set up color dot plots
		loadct,5,/silent
		;loadct,0,/silent
		plotsym,0,1,/fill
		if j eq 0 then brange=[0.,15.];brange=[-10.,10.];minmax(sw[*,5,i])

		;Plot spiral extrapolation
		;oplot,[spirals[*,0],spirals[*,0]], [spirals[*,1],spirals[*,1]]*!dtor,ps=3
		nspiral=n_elements(spirals[*,0])
		if keyword_set(no_alpha) then for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((spirals[m,5])-brange[0])/(brange[1]-brange[0])*255. $
			else for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=(alog10(spirals[m,5])>0.)/3.*255. < 255.

		;Plot color dots (non-extrapolated, only back propagated)
		if keyword_set(no_alpha) then for k=0l,nblob-1l do oplot,/polar,[sw[k,0,i],sw[k,0,i]],[sw[k,1,i],sw[k,1,i]]*!dtor,ps=8,color=((sw[k,5,i])-brange[0])/(brange[1]-brange[0])*255. $
			else for k=0l,nblob-1l do oplot,/polar,[sw[k,0,i],sw[k,0,i]],[sw[k,1,i],sw[k,1,i]]*!dtor,ps=8,color=(alog10(sw[k,5,i])>0.)/3.*255. < 255.
			
;Plot color bar
		if j eq 2 then if keyword_set(no_alpha) then color_table, [0.,10.],/right,/top,title='|B| [nT]' else color_table, [0,3],/right,/top,title='LOG |B| [nT]'

		
;Save the Stereo A position
		if j eq 1 then wbest1=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))

;Plot space craft
		if j eq 2 then begin
			setcolors,/sys,/quie,/sile
;Plot Stereo B
			;stb_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'B', system = 'HCI', /degrees )
			wbest2=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
			plotsym,0,3,/fill
			;oplot,/polar,[stb_pos[0],stb_pos[0]],[stb_pos[1],stb_pos[1]],ps=8,color=!green
			oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!green
			plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!black
;Plot Stereo A
			;sta_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'A', system = 'HCI', /degrees )
			;wbesta=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
			plotsym,0,3,/fill
			;oplot,/polar,[sta_pos[0],sta_pos[0]],[sta_pos[1],sta_pos[1]],ps=8,color=!red
			oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!red
			plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!black
;Plot Earth
			plotsym,0,3,/fill
			oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!blue
			plotsym,0,2,/fill & oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!black
		endif
		
skip_plot_points:

	endfor	

	if keyword_set(plot_points) then window_capture,file=plotp+'frame_'+string(i,form='(I04)'),/png

	if keyword_set(interp_maps) then begin
		sw_arrays=spiral_arrays[1:*,*]
		sw_arrays[*,1]=sw_theta_shift(sw_arrays[*,1])
		mapstack=sw_interp(sw_arrays,thistime,/heliocentric,plot_interp=plot_interp,save_maps=save_maps,constants=constants_arr) ;,/cylindrical
	endif

	if keyword_set(test) then begin
		blah=''
		read,blah
	endif

endfor


end