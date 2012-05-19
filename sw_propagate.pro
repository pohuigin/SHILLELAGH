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
;
;	3.) Only use certain space craft to model the wind. Default assumes all.
;	IDL>sw_propagate,'12-dec-2008',spacecraft=['sta','omni']
;
;	4.) Only compute a single time.
;	IDL>sw_propagate,'12-dec-2008',trange=0
;
;	5.) Run a single STEREO instrument with a specific theta range
;	IDL>sw_propagate,'15-jun-2008',trange=0,space='sta', thetarange=[[0,0],[-60,-30],[0,0]]
;
;NOTES:
;	This algorithm must be run in IDL 32 bit mode due to its utilization of SPICE.
;
;HISTORY:
;	12-Nov-2010 - P.A.Higgins - Written
;
;TODO: 	1) 	correct all properties to be in solar equatorial plane instead of space craft plane
;		2) 	allow option for plotting points for different properties, not just B mag
;		3) 	read br, bn and bt into model- good for spotting coronal holes
;				bt should show very clear CME signal
;				bperp,spiral will tell the polarity of coronal holes
;		4) 	include properties like b rotation (clearest cme signal), beta sw, epsilon parameter (geo effective ness)
;		5) 	include SW inclination? B inclination (for current sheet topology)
;		6) 	AMDA property which goes negative when you cross the current sheet
;		7) 	allow adding more space craft in situ data slices.
;		8) 	allow option to plot 4 properties for CME detection or coronal hole detection
;				CMEs:[bmag, brotation or angle, density? temp?, epsilon?]
;				CHs:[btangent to spiral, density, velocity, temp?]
;		9) 	allow option to overplot individual spirals by inputting a list of heliographic longitudes
;		10) over plot colored arrows at the space craft to show the direction of the passing magnetic field
;		11) use the "los" bfield angle to correct the spiral
;				if b is along spiral, there should be no change
;				if b is perp spiral, it should pull the spiral to be in the direction of B
;		12) include mariner, voyager1, voyager 2, and map them to some radius range from saturn -> 2.5 or 3 AU 
;				and look at the interface between the inner and outer space craft data.
;		13) add the effect of the solar wind acceleration between 1 and 100 R sun (should go as r? r2? r3?) 
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

save,map_stack,file=sw_paths(/sav)+'sw_properties_heliospheric_'+time2file(anytim(tim,/vms))+'.sav',/compress

return,map_stack

end

;-------------------------------------------------------------------------->
;Propagate an array of plasma blobs along their spiral using given velocity
function sw_spiral, inblob, theta_range=inthetarange, theta_bin=inthetabin, r_range=inrrange, $
	constants=constants_arr, nparam=nspiralparam, earthlon=this_earthlon, no_alpha=no_alpha

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

;stop

;Get rid of small gap due to non integer theta range and large bin.
thetarange=thetarange+[-thetabin/2.,2.*thetabin]

blob_arr=inblob
if (where(blob_arr[*,0] le rrange[1] and blob_arr[*,0] gt rrange[0]))[0] eq -1 then return,nan
blob_arr=blob_arr[where(blob_arr[*,0] le rrange[1] and blob_arr[*,0] gt rrange[0]),*]

nblob=n_elements(blob_arr[*,0])
blank=fltarr(1,nspiralparam)
spiral_arr=blank
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
	thisspiral=fltarr(npts,nspiralparam)	
	
	;Load SW properties for each point along spiral
	;for j=2,n_elements(thisblob)-1 do thisspiral[*,j]=thisblob[j]
	thisspiral[*,2]=thisblob[2] ;velocity
	thisspiral[*,3]=thisblob[3] ;density
	thisspiral[*,4]=thisblob[4] ;temperature
	thisspiral[*,5]=thisblob[5] ;Bmagnitude
	;skip 6, the time column
	thisspiral[*,6]=thisblob[7] ;Bradial
	
	;Delta theta from blob starting point
	theta_spiral=findgen(npts)*thetabin+thisthetarange[0]+this_earthlon
	
	;delta R = R(blob) - delta theta (velocity(blob) / omega_sun)
	if bloblt0 then r_spiral = thisblob[0]-( thisblob[2] / omegasun ) * (theta_spiral-thisblob[1]) $
		else r_spiral = thisblob[0]-( thisblob[2] / omegasun ) * (theta_spiral-thisblob[1])
	thisspiral[*,0]=r_spiral
	
	;theta = theta(blob) + delta theta
;Center spiral at angle of earth NOT at angle of "blob"	
;	thisspiral[*,1]=thisblob[1]+theta_spiral
	thisspiral[*,1]=theta_spiral

;Assume 1/R^alpha dependence of Solar Wind properties

;Filter out negative R values
	if bloblt0 then begin
		if rbloblt0 le 0 then continue
		if (where(thisspiral[*,0] gt 0))[0] eq -1 then continue else begin
			wgoodspiral=where(thisspiral[*,0] gt 0)
			if wgoodspiral[0] ne -1 then thisspiral=thisspiral[wgoodspiral,*]
		endelse
	
		if not keyword_set(no_alpha) then begin
;B Field
			thisspiral[*,5]=thisspiral[*,5]*(rbloblt0/thisspiral[*,0])^alpha_b
;Temperature
			thisspiral[*,4]=thisspiral[*,4]*(rbloblt0/thisspiral[*,0])^alpha_t
;Density
			thisspiral[*,3]=calc_helio_props(thisspiral[*,0], rbloblt0, thisspiral[*,3], /density, constants=constants_arr);thisspiral[*,3]*(rbloblt0/thisspiral[*,0])^alpha_rho
		endif
	endif else begin
		if not keyword_set(no_alpha) then begin
;B Field
			thisspiral[*,5]=thisspiral[*,5]*(thisblob[0]/thisspiral[*,0])^alpha_b
;Temperature
			thisspiral[*,4]=thisspiral[*,4]*(thisblob[0]/thisspiral[*,0])^alpha_t
;Density
			thisspiral[*,3]=calc_helio_props(thisspiral[*,0], thisblob[0], thisspiral[*,3], /density, constants=constants_arr);thisspiral[*,3]*(thisblob[0]/thisspiral[*,0])^alpha_rho
		endif
	endelse

	;Filter spiral points with R lt 0
	wgt0=where(thisspiral[*,0] gt 0)
	if wgt0[0] ne -1 then thisspiral=thisspiral[wgt0,*]
	
	spiral_arr=[spiral_arr,thisspiral]
	
;if i eq 50 then stop

endfor

;stop

outspiral=spiral_arr[1:*,*]
return,outspiral

end

;-------------------------------------------------------------------------->
;parker model solar wind speed vs. r from sun=> r-r_o = (v_sw/omega_sun)*(theta-theta_o)
;												v_sw=omega_sun*(r-r_o)/(theta-theta_o)
;intime = 	'DD-Mon-YYYY' to propagate to
;trange = 	[2 elements] days relative to input time
;tbin	=	days between extrapolation solutions
;if you want to run the code for a list of specified times, 
;	set INTIME to the list of times and set TIMEREANGE=0.
;thetarange = if only running one SC, then you can set thetarange to a 2 element array of the +- angle (centered at 0, where 0=earth longitude) in degree to propagate
;outflist = output a list of the save files generated

pro sw_propagate, intime, trange=intrange, tbin=intbin, savpath=insavpath, no_alpha=no_alpha, full360=full360, $
	test=test, plot_points=plot_points, spacecraft=spacecraft, res1tore=res1tore, $
	interp_maps=interp_maps, save_maps=save_maps, plot_interp=plot_interp, filetag=infiletag, $
	custom=fcustom, nostereo=nostab, thetarange=inthetarange, outflist=fsavlist

if n_elements(insavpath) gt 0 then savpath=insavpath else savpath=sw_paths(/sav)
if n_elements(infiletag) gt 0 then filetag=infiletag else filetag=''
;initialise list of save files to output
fsavlist=''

time=anytim(intime)
if n_elements(intrange) lt 1 then trange=[time-5.*24.*3600.,time+10.*24.*3600.] else begin
	if n_elements(intrange) eq 1 then trange=intrange
	if n_elements(intrange) gt 1 then trange=time+intrange*24.*3600.
endelse
if n_elements(intbin) lt 1 then tbin=12.*3600. else tbin=intbin*24.*3600.

;if trange[0] ge anytim('1-jan-2007') then dostereo=1 else dostereo=0 
if time[0] ge anytim('1-jan-2007') then dostereo=1 else dostereo=0 
if keyword_set(nostab) then dostereo=0

if keyword_set(plot_points) then window,xs=750,ys=750

plotp=sw_paths(/plotp) ;'~/science/procedures/cme_propagation/sw_prop_movie/'

;SW Properties fall off as 1/r^alpha
alpha_b=2.5d ;between 2 and 3
alpha_rho=4.d ;!!!NOT USED: switched to using empirical model in Chen 1996 (Eqn. 17.)
alpha_t=1.d
if keyword_set(no_alpha) then begin & alpha_b=0d & alpha_rho=0d & alpha_t=0d & endif

;Physical Constants
au_km=149.6d6 ;1 AU in kilometers
vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees or solar ascending node?? 
nan=0/0.
r_sun=6.955d5 ;km
omegasun=360d/(25.2d*3600d*24d) ;in degrees/second from diff. rot. of latitudes (-10 -> +10)

constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]

if n_elements(trange) eq 1 then begin 
	time_arr=time
	ntime=n_elements(time)
	trange=minmax(time_arr)
endif else begin
	ntime=round((trange[1]-trange[0])/tbin)
	time_arr=findgen(ntime)*tbin+trange[0]
endelse

;Read in situ data SC=([r,theta,v,rho,temp,bmag,t],nblob)
if dostereo then nspacecraft=3 else nspacecraft=1
if not keyword_set(res1tore) then begin
	sc_arr0=sw_get_data(trange,/omni,const=constants_arr,custom=fcustom)
	if dostereo then begin 
		sc_arr1=sw_get_data(trange,/sta,const=constants_arr,custom=fcustom)
;!!!TEMP
;		sc_arr2=sc_arr1
		sc_arr2=sw_get_data(trange,/stb,const=constants_arr,custom=fcustom)
	endif
	save,sc_arr0,sc_arr1,sc_arr2,file=sw_paths(/sav)+'sw_run_data_load_restore.sav'
endif else restore,sw_paths(/sav)+'sw_run_data_load_restore.sav'

;Create SW propagation SW(r,theta,v,rho,temp,bmag,t_sc) array sw=[nblob,nparam=7,ntime]
nswparam=8
nspiralprop=nswparam-1 ;do -1 to exclude time column

nblob0=n_elements(sc_arr0[0,*])
sw0=fltarr(nblob0,nswparam,ntime)
if dostereo then begin
	nblob1=n_elements(sc_arr1[0,*])
	sw1=fltarr(nblob1,nswparam,ntime)
	nblob2=n_elements(sc_arr2[0,*])
	sw2=fltarr(nblob2,nswparam,ntime)
endif

;Find range to solve spiral for each space craft: 1/2*Difference between earth and spacecraft + 90 degrees
wbestearth=(where(abs(time_arr[ntime/2.]-sc_arr0[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr0[6,*]))))[0]
;wbeststa=(where(abs(time_arr[ntime/2.]-sc_arr1[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr1[6,*]))))[0]
;wbeststb=(where(abs(time_arr[ntime/2.]-sc_arr2[6,*]) eq min(abs(time_arr[ntime/2.]-sc_arr2[6,*]))))[0]
if dostereo then begin 
	stbhglon=(GET_STEREO_LONLAT( anytim(time_arr[ntime/2.],/vms), 'B', system = 'HEE', /degrees ))[1]
	stahglon=(GET_STEREO_LONLAT( anytim(time_arr[ntime/2.],/vms), 'A', system = 'HEE', /degrees ))[1]
	earthsta=stahglon;abs(sc_arr0[1,wbestearth]-sc_arr1[1,wbeststa])/2.
	earthstb=stbhglon;abs(sc_arr0[1,wbestearth]-sc_arr2[1,wbeststb])/2.
endif else begin
	earthsta=0.
	earthstb=0.
endelse

if n_elements(inthetarange) gt 0 then thetarange=inthetarange else begin
;	Calculate ranges to model over for each space craft
;	thetarange=ceil([[earthstb/2., earthsta/2.], $
;		[-earthsta/2.,180.], $
;		[-180.,earthstb/2.]])
	thetarange=ceil([[earthstb/2.-10., earthsta/2.+10.], $ ;tried to fix non-matching up limits...
		[earthsta/2.+10.,180.], $
		[-180.,earthstb/2.-10.]])
;	thetarange=sw_theta_shift(thetarange)
endelse

if n_elements(inthetarange) lt 1 then begin
	if keyword_set(full360) then thetarange=ceil([[0.,360.], [0,360.], [0,360.]])
	if n_elements(spacecraft) gt 0 then begin
		if n_elements(spacecraft) eq 2 then begin
			if strjoin(strlowcase(spacecraft)) eq ['stastb'] or strjoin(strlowcase(spacecraft)) eq ['stbsta'] then $
				thetarange=ceil([[0.,0.], [-earthsta,180.-earthsta], [-180.-earthstb,-1.*earthstb]])
			if strjoin(strlowcase(spacecraft)) eq ['staomni'] or strjoin(strlowcase(spacecraft)) eq ['omnista'] then $
				thetarange=ceil([[-180-earthsta/2., earthsta/2.], [-earthsta/2.,180.-earthsta], [0.,0.]])
			if strjoin(strlowcase(spacecraft)) eq ['stbomni'] or strjoin(strlowcase(spacecraft)) eq ['omnistb'] then $
				thetarange=ceil([[earthstb/2.-10., 180+earthsta/2.], [0,0], [-180.-earthstb,-1.*earthstb/2.]])
		endif
		if n_elements(spacecraft) eq 1 then begin
			if strlowcase(spacecraft) eq 'omni' then thetarange=ceil([[-180,180],[0,0],[0,0]])
			if strlowcase(spacecraft) eq 'sta' then thetarange=ceil([[0,0],[-180,180],[0,0]])
			if strlowcase(spacecraft) eq 'stb' then thetarange=ceil([[0,0],[0,0],[-180,180]])	
		endif
	endif
endif
;stop

;Run through each propagation time within TIMERANGE using TBIN
for i=0l,ntime-1l do begin
	spiral_arrays=transpose(fltarr(nspiralprop)) ;!!!added an element to stick bradial

;Run through each spacecraft
	for j=0,nspacecraft-1 do begin
	
;Skip this iteration if the range for the instrument is [0,0]
		if thetarange[0,j] eq 0 and thetarange[1,j] eq 0 then goto,skip_instrument 

		ex=execute('sw=sw'+strtrim(j,2)+' & nblob=nblob'+strtrim(j,2)+' & sc_arr=sc_arr'+strtrim(j,2))
		thistime=time_arr[i]
;Fill time when blob measured is at spacecraft
		sw[*,6,i]=sc_arr[6,*]
		
;Fill velocity
		sw[*,2,i]=sc_arr[2,*]
		
;Fill B field
		sw[*,5,i]=sc_arr[5,*]

;Fill Bradial field
		sw[*,7,i]=sc_arr[7,*]		
		
;Fill density
		sw[*,3,i]=sc_arr[3,*]
		
;Fill temperature
		sw[*,4,i]=sc_arr[4,*]
		
;Calculate radii
		sw[*,0,i]=sc_arr[0,*]-(sc_arr[6,*]-thistime)*sc_arr[2,*]

;Assume 1/R^alpha dependence of Solar Wind properties
sw=double(sw)
;B Field
	if not keyword_set(no_alpha) then begin
		sw[*,5,i]=sw[*,5,i]*(au_km/abs(sw[*,0,i]))^alpha_b
;Temperature
		sw[*,4,i]=sw[*,4,i]*(au_km/abs(sw[*,0,i]))^alpha_t
;Density
		sw[*,3,i]=calc_helio_props(sw[*,0,i], sc_arr[0,*], sw[*,3,i], /density, constants=constants_arr);sw[*,3,i]*(au_km/abs(sw[*,0,i]))^alpha_rho
	endif
	
;Calculate earth longitude array
		if j eq 0 or n_elements(this_earthlon) lt 1 then begin
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

;stop

;Plot radial propagation
;window,xs=650,ys=650
;plot,sw[*,0]/au_km>0,sw[*,1]*!dtor,/polar,ps=4,xrange=[-1.5,1.5],yrange=[-1.5,1.5],chars=2,/iso,xmargin=[5,1],ymargin=[5,1],xtit=anytim(thistime,/vms,/date)+' [AU]'
;setcolors,/sys,/silent
;plotsym,0,2,/fill
;plots,0,0,ps=8,color=!orange
;oplot,[1,1],[this_earthlon,this_earthlon]*!dtor,ps=8,color=!blue,/polar

;Test array bounds
;pmm,(sw[*,1])[where(sw[*,0]/au_km gt 0 and sw[*,0]/au_km lt 1.5)]
;pmm,(sc_arr[1,*])[where(sw[*,0]/au_km gt 0 and sw[*,0]/au_km lt 1.5)]
;pmm,(sw_theta_shift(sc_arr[1,*], /n180to180))[where(sw[*,0]/au_km gt 0 and sw[*,0]/au_km lt 1.5)]
;wlt1p5au=where(sw[*,0]/au_km gt 0 and sw[*,0]/au_km lt 1.5)

;!!!CAUSES ERROR - (WRONG)
;Crop array to only include the theta range for the given SC
;		winthetarng=where(sw_theta_shift(sw[*,1,i]-this_earthlon, /n180to180) ge thetarange[0,j] and sw_theta_shift(sw[*,1,i]-this_earthlon, /n180to180) le thetarange[1,j])
;		sw=sw[winthetarng,*,i]
;Update number of elements variable
;		nblob=n_elements(sw[*,1,i])
		
;Test spiral
;sw1=sw[winthetarng,*,i]
;plot,sw1[*,0]/au_km>0,sw1[*,1]*!dtor,/polar,ps=4,xrange=[-1.5,1.5],yrange=[-1.5,1.5],chars=2,/iso,xmargin=[5,1],ymargin=[5,1],xtit=anytim(thistime,/vms,/date)+' [AU]'
;wltau=where(sw1[*,0]/au_km gt 0 and sw1[*,0]/au_km lt 1.5)
;tspirals=sw_spiral(sw1[554,*,i],const=constants_arr, theta_range=thetarange[*,j], nparam=nspiralprop)
;oplot,tspirals[*,0]/au_km,tspirals[*,1]*!dtor,/polar,ps=4,color=!red
;oplot,[0,1],([thetarange[0,j],thetarange[0,j]]+this_earthlon)*!dtor,lines=2,/polar
;oplot,[0,2],([thetarange[1,j],thetarange[1,j]]+this_earthlon)*!dtor,lines=2,/polar
;plots,0,0,ps=8,color=!orange
;oplot,[1,1],[this_earthlon,this_earthlon]*!dtor,ps=8,color=!blue,/polar

;Calculate spirals for each blob
 		spirals=sw_spiral(sw[*,*,i],const=constants_arr, theta_range=thetarange[*,j], nparam=nspiralprop, earthlon=this_earthlon, no_alpha=no_alpha)
		if n_elements(spirals) gt 1. then spiral_arrays=[spiral_arrays,spirals]
		
;Get rid of meaningless data  (R < R_sun)
		wbad=where(sw[*,0,i] lt r_sun)
		if wbad[0] ne -1 then sw[wbad,0,i]=nan
		
if not keyword_set(plot_points) then goto,skip_plot_points

;Plot the axes
		if j eq 0 then begin
			setcolors,/sys,/sil,/qui
			plot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,ps=3,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso,title=anytim(thistime,/vms)
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
			plot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,/nodata,/noerase,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso
;Plot the Sun
			plotsym,0,1,/fill
			oplot,/polar,[0,0]*au_km,[0,0],ps=8,color=!orange
			
;Plot Stereo B
			stb_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'B', system = 'HCI', /degrees )
			;wbest2=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
			plotsym,0,1,/fill
			oplot,/polar,[stb_pos[0],stb_pos[0]],[stb_pos[1],stb_pos[1]]*!dtor,ps=8,color=!white
			;oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!green
			;plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest2],sc_arr[0,wbest2]],[sc_arr[1,wbest2],sc_arr[1,wbest2]]*!dtor,ps=8,color=!black
;Plot Stereo A
			sta_pos=GET_STEREO_LONLAT( anytim(thistime,/vms), 'A', system = 'HCI', /degrees )
			;wbesta=where(min(abs(sc_arr[6,*]-thistime)) eq abs(sc_arr[6,*]-thistime))
			plotsym,0,1,/fill
			oplot,/polar,[sta_pos[0],sta_pos[0]],[sta_pos[1],sta_pos[1]]*!dtor,ps=8,color=!white
			;oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!red
			;plotsym,0,2,/fill & oplot,/polar,[sc_arr[0,wbest1],sc_arr[0,wbest1]],[sc_arr[1,wbest1],sc_arr[1,wbest1]]*!dtor,ps=8,color=!black
;Plot Earth
			polrec, earth_rad[0]*au_km, this_earthlon, earthpos_x, earthpos_y, /degrees
			plotsym,0,1,/fill
			oplot,[earthpos_x,earthpos_x],[earthpos_y,earthpos_y],ps=8,color=!white
			;oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!yellow
			;plotsym,0,2,/fill & oplot,/polar,[earth_rad[0],earth_rad[0]]*au_km,[this_earthlon,this_earthlon]*!dtor,ps=8,color=!black
		endif

		oplot,/polar,sw[*,0,i],sw[*,1,i]*!dtor,ps=3		

skip_plot_points:

skip_instrument:

	endfor	

	if not keyword_set(test) then if keyword_set(plot_points) then window_capture,file=plotp+'frame_'+time2file(anytim(thistime,/vms)),/png

	sw_arrays=spiral_arrays[1:*,*]
	sw_arrays[*,1]=sw_theta_shift(sw_arrays[*,1])
;Interpolate points to create solution images
	if keyword_set(interp_maps) then begin
		;sw_arrays=spiral_arrays[1:*,*]
		;sw_arrays[*,1]=sw_theta_shift(sw_arrays[*,1])
		mapstack=sw_interp(sw_arrays,thistime,/heliocentric,plot_interp=plot_interp,save_maps=save_maps,constants=constants_arr) ;,/cylindrical
	endif
	
	if not keyword_set(test) then begin
		thisfsav=savpath+'sw_properties_points_'+time2file(anytim(thistime,/vms))+'_'+filetag+'.sav'
		save,sw_arrays,file=thisfsav,/compress
		fsavlist=[fsavlist,thisfsav]
	endif

	if keyword_set(test) then begin
		blah=''
		read,blah
	endif

endfor

if n_elements(fsavlist) gt 1 then fsavlist=fsavlist[1:*]

end