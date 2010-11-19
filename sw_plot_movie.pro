pro sw_plot_movie, intime, trange=intrange, tbin=intbin

time=anytim(intime)
if n_elements(intrange) lt 1 then trange=[time-5.*24.*3600.,time+10.*24.*3600.] else trange=time+intrange*24.*3600.
if n_elements(intbin) lt 1 then tbin=12.*3600. else tbin=intbin*24.*3600.

savp='~/science/procedures/cme_propagation/sw_prop_save/'
plotp='~/science/procedures/cme_propagation/heliosphere_property_plots/'
moviep='~/science/procedures/cme_propagation/claire_events/claire_20100302_movie/'

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
window, xsize=1200, ysize=450	

for i=0,ntime-1 do begin
	!p.multi=[3,3,1]
	!p.color = 0
	sw_plot_points, anytim(time_arr[i],/vms), /velocityplot, pmulti=[3,3,1], $
		save_plot=save_plot, title_string='Velocity', cbarpos=[.05,.92,.3,.98];footpoints=footpoints, chvelocity=inchvelocity
	!p.multi=[2,3,1]
	!p.color = 0
	sw_plot_points, anytim(time_arr[i],/vms), /densityplot, pmulti=[2,3,1], $
		save_plot=save_plot,title_string='Density', cbarpos=[.4,.92,.65,.98];footpoints=footpoints, chvelocity=inchvelocity
	!p.multi=[1,3,1]
	!p.color = 0
	sw_plot_points, anytim(time_arr[i],/vms), /magneticfieldplot, pmulti=[1,3,1], $
		save_plot=save_plot, title_string='B Field', cbarpos=[.7,.92,.95,.98];footpoints=footpoints, chvelocity=inchvelocity

	window_capture,file=moviep+'shillelagh_vel_dens_b_'+string(i,form='(I03)'),/png
	erase
endfor


end