pro sw_plot_movie, intime, trange=intrange, tbin=intbin, _extra=_extra, $
	velocity=velocity,density=density,magnetic=magnetic, wsize=wsize

time=anytim(intime)
if n_elements(intrange) lt 1 then trange=[time-5.*24.*3600.,time+10.*24.*3600.] else trange=time+intrange*24.*3600.
if n_elements(intbin) lt 1 then tbin=12.*3600. else tbin=intbin*24.*3600.
if n_elements(wsize) ne 1 then wsize=500

savp='~/science/data/cme_propagation/sw_prop_save/'
plotp='~/science/data/cme_propagation/heliosphere_property_plots/'
;moviep='~/science/procedures/cme_propagation/claire_events/claire_20100302_movie/'
moviep='~/science/data/cme_propagation/mars_event/'

if n_elements(trange) eq 1 then ntime=1. else ntime=round((trange[1]-trange[0])/tbin)
if ntime eq 1. then time_arr=time else time_arr=findgen(ntime)*tbin+trange[0]

if ntime lt 2 then begin
	print,'MUST HAVE > 1 FRAMES FOR A MOVIE!'
	print,'time='+strtrim(time,2)+' ntime='+strtrim(ntime,2)+' trange='+strtrim(trange,2)+' tbin='+strtrim(tbin,2)
endif
	
;!p.position=[.07,.05,.97,.95]
!p.background = 255
!p.color = 0   
!p.charthick=1
!p.charsize = 1.4
nplot=keyword_set(velocity)+keyword_set(density)+keyword_set(magnetic)
window, xsize=nplot*wsize, ysize=wsize	

for i=0,ntime-1 do begin
	!p.multi=[nplot,nplot,1]
	print,[nplot,nplot,1]
	
	if keyword_set(velocity) then begin
		!p.color = 0
		sw_plot_points, anytim(time_arr[i],/vms), /velocityplot, $
			save_plot=save_plot, title_string='Velocity', cbarpos=[.05,.92,.3,.98], /nowindow, _extra=_extra;footpoints=footpoints, chvelocity=inchvelocity
	endif
	
	!p.multi=[nplot-1,nplot,1]
	if keyword_set(density) then begin
		!p.color = 0
		sw_plot_points, anytim(time_arr[i],/vms), /densityplot, $
			save_plot=save_plot,title_string='Density', cbarpos=[.4,.92,.65,.98], /nowindow, _extra=_extra;footpoints=footpoints, chvelocity=inchvelocity
	endif
	
	!p.multi=[nplot-2,nplot,1]
	if keyword_set(magnetic) then begin
		!p.color = 0
		sw_plot_points, anytim(time_arr[i],/vms), /magneticfieldplot, $
			save_plot=save_plot, title_string='B Field', cbarpos=[.7,.92,.95,.98], /nowindow, _extra=_extra;footpoints=footpoints, chvelocity=inchvelocity
	endif

	window_capture,file=moviep+'shillelagh_vel_dens_b_'+string(i,form='(I03)'),/png
	erase
endfor


end