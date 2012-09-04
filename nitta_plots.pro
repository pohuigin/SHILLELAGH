pro nitta_plots

path='~/science/data/cme_propagation/nitta/'

times=['12-Jun-2010 01:00', '7-Aug-2010 18:00', '14-Aug-2010 10:00', '18-Aug-2010 05:00'];'01-Jan-2009',

;sw_propagate, times, trange=0,spacecraft='omni', /nostereo;, thetarange=ceil([[0.,360], [0,360.], [0,360.]]);, /full360
sw_propagate, times, trange=0,spacecraft=['omni','stb'];, /nostereo;, thetarange=ceil([[0.,360], [0,360.], [0,360.]]);, /full360

!p.multi=0
;window,xs=700,ys=700
;window,xs=900,ys=300

for i=0,n_elements(times)-1 do begin

;	setplotenv,file=path+time2file(times[i])+'_vel_omni_sta.eps',/ps,ys=12,xs=12
	sw_plot_points, times[i], /velocity, /planets;,/nowin;,pmulti=[3,3,1];, /no_plot_sc
;	closeplotenv
;window_capture,file=path+time2file(times[i])+'_vel_omni',/png
window_capture,file=path+time2file(times[i])+'_vel_omni_stereo_b',/png

;t=findgen(339)+10.
;tt=fltarr(339)+3d8
;oplot,tt,t/359.*2.*!pi,/polar,thick=4,color=255
;stop

;	setplotenv,file=path+time2file(times[i])+'_dens_omni_sta.eps',/ps,ys=12,xs=12
	sw_plot_points, times[i], /density, /planets;,/nowin;,pmulti=[2,3,1];,/nowin, /no_plot_sc
;	closeplotenv
;window_capture,file=path+time2file(times[i])+'_dens_omni',/png
window_capture,file=path+time2file(times[i])+'_dens_omni_stereo_b',/png

;	setplotenv,file=path+time2file(times[i])+'_brad_omni_sta.eps',/ps,ys=12,xs=12
	sw_plot_points, times[i], /radial, /planets;,/nowin;,pmulti=[1,3,1];, /no_plot_sc,/nowin
;	closeplotenv
;window_capture,file=path+time2file(times[i])+'_brad_omni',/png
window_capture,file=path+time2file(times[i])+'_brad_omni_stereo_b',/png
	
;stop

endfor










end