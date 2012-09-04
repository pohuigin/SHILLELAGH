pro run_sept22


;sw_propagate,'22-sep-2011',trange=[-3,4],tbin=1./12.,spacecr=['omni','stb'],thetar=[[-45,180],[0,0],[-180,-45]]

window
!p.multi=[0,2,1]
stop

mpath='~/science/data/cme_propagation/cme_20110922/'
fcmemov=file_search(mpath+'save/*')
nframe=n_elements(fcmemov)
for i=0,nframe-1 do begin &$
erase &$
	sw_plot_points, file=fcmemov[i],/den,/label,/plan,xran=[1.5d8,-1.5d8],yran=[1.5d8,-1.5d8],/nowin,pmulti=[1,2,1],/nocolor &$
	sw_plot_points, file=fcmemov[i],/vel,/label,/plan,xran=[1.5d8,-1.5d8],yran=[1.5d8,-1.5d8],/nowin,pmulti=[2,2,1],/nocolor &$
	window_capture,file=mpath+'plot/cmeframe_'+string(i,form='(I04)') &$
endfor

stop

;	sw_plot_points, file=fcmemov[i],/radialfieldplot,xran=[1.5d8,-1.5d8],yran=[1.5d8,-1.5d8],/label,pmulti=[3,3,1],/plan,/nowin &$


;	xyouts,.1,.15,'Density',/norm,color=255,chars=2 & xyouts,.6,.15,'Velocity',/norm,color=255,chars=2 &$

;	oplot,[0,1d8],([earthpos,earthpos]-8.)*!dtor,thick=4,color=0,lines=2,/polar &$

end