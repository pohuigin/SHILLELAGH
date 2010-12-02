pro claire_insitu

plotp='~/science/data/cme_propagation/claire_events/'

time='13-mar-2010 00:00'
tau=24.*10.*3600.

trange=[anytim(time),anytim(time)+tau]
tdata=[anytim(time),anytim(time)+tau]

sw_omni=sw_get_data(tdata, /omni);, sta=sta, stb=stb, ace=ace, wind=wind, $
	;constants=constants_arr, geoindex=geoindex
w0=where(abs(tdata[0]-sw_omni[6,*]) eq min(abs(tdata[0]-sw_omni[6,*])))
w1=where(abs(tdata[1]-sw_omni[6,*]) eq min(abs(tdata[1]-sw_omni[6,*])))
sw_omni=sw_omni[*,w0:w1]

sw_stb=sw_get_data(tdata, /stb);, sta=sta, stb=stb, ace=ace, wind=wind, $
	;constants=constants_arr, geoindex=geoindex
w0=where(abs(tdata[0]-sw_stb[6,*]) eq min(abs(tdata[0]-sw_stb[6,*])))
w1=where(abs(tdata[1]-sw_stb[6,*]) eq min(abs(tdata[1]-sw_stb[6,*])))
sw_stb=sw_stb[*,w0:w1]

sw_geo=sw_get_data(tdata, /geoindex);, sta=sta, stb=stb, ace=ace, wind=wind, $
	;constants=constants_arr, geoindex=geoindex
wg0=(where(abs(tdata[0]-sw_geo[0,*]) eq min(abs(tdata[0]-sw_geo[0,*]))))[0]
wg1=(where(abs(tdata[1]-sw_geo[0,*]) eq min(abs(tdata[1]-sw_geo[0,*]))))[0]
sw_geo=sw_geo[*,wg0:wg1]

window,xs=700,ys=700
!p.multi=[0,1,3]
!p.color=0
!p.background=255
!p.charsize=2
setcolors,/sys,/sile,/quie

utplot,sw_omni[6,*]-anytim(time),sw_omni[2,*],time,ytitle='Velocity',/xsty,tit='OMNI',yran=[250,600]
utplot,sw_omni[6,*]-anytim(time),sqrt((sw_omni[7,*])^2.+(sw_omni[8,*])^2.+(sw_omni[9,*])^2.),time,ytitle='Magnetic Field',yran=[-10,10],/xsty,/ysty,thick=2
oplot,sw_omni[6,*]-anytim(time),sw_omni[8,*],color=!red,ps=-4
oplot,sw_omni[6,*]-anytim(time),sw_omni[7,*],color=!blue,ps=-4
oplot,sw_omni[6,*]-anytim(time),sw_omni[9,*],color=!forest,ps=-4
hline,0,lines=2
utplot,sw_omni[6,*]-anytim(time),sw_omni[3,*],time,ytit='Density',/xsty

window_capture,file=plotp+'insitu_cme_omni',/png

stop

utplot,sw_stb[6,*]-anytim(time),sw_stb[2,*],time,ytitle='Velocity',/xsty,tit='STB',yran=[250,600]
utplot,sw_stb[6,*]-anytim(time),sqrt((sw_stb[7,*])^2.+(sw_stb[8,*])^2.+(sw_stb[9,*])^2.),time,ytitle='Magnetic Field',yran=[-10,10],/xsty,/ysty,thick=2
oplot,sw_stb[6,*]-anytim(time),sw_stb[8,*],color=!red,ps=-4
oplot,sw_stb[6,*]-anytim(time),sw_stb[7,*],color=!blue,ps=-4
oplot,sw_stb[6,*]-anytim(time),sw_stb[9,*],color=!forest,ps=-4
hline,0,lines=2
utplot,sw_stb[6,*]-anytim(time),sw_stb[3,*],time,ytit='Density',/xsty

window_capture,file=plotp+'insitu_cme_stb',/png

end

;------------------------------------------------------------------------------>

pro claire_plot

au_km=149.6d6

restore,'~/science/data/test/sw_properties_points_20100312_0700.sav',/verb
spirals=sw_arrays

;Set up color dot plots
loadct,5,/silent
;loadct,0,/silent
plotsym,0,1,/fill

nspiral = n_elements(spirals[*,0])
propnum=2

plot,/polar,spirals[*,0],spirals[*,1]*!dtor,ps=3,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso,/nodata

case propnum of 
	2: velrange=[100.,700.]
	3: densrange=[0.,3]
	4: temprange=[3,8]
	5: brange=[0.,3.]
endcase
case propnum of 
	2: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((spirals[m,propnum])-velrange[0])/(velrange[1]-velrange[0])*255.
	3: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>densrange[0] < densrange[1])-densrange[0])/(densrange[1]-densrange[0])*255. < 255.
	4: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>temprange[0] < temprange[1])-temprange[0])/(temprange[1]-temprange[0])*255. < 255.
	5: for m=0l,nspiral-1l do oplot,[spirals[m,0],spirals[m,0]],[spirals[m,1],spirals[m,1]]*!dtor,/polar,ps=8,color=((alog10(spirals[m,propnum])>brange[0] < brange[1])-brange[0])/(brange[1]-brange[0])*255. < 255.
endcase

plot,/polar,spirals[*,0],spirals[*,1]*!dtor,/nodata,/noerase,xr=[au_km,-au_km]*1.2,yr=[au_km,-au_km]*1.2,/iso

stop

n = n_elements(sw_array[*,0])

ni = n/1000
dn = n-ni*1000

for j = 0, ni do begin
       for k = 0, 999 do begin

       i = k+(j*1000.)
       PRINT, k, j, I
       c = alog( sw_array[*,3]/max(sw_array[*,3])*255.) - min( alog( sw_array[*,3]/max(sw_array[*,3])*255.))*5.

       oplot,/polar,sw_array[i:i,0],sw_array[i:i,1]*!dtor,ps=8,color=c[i], thick = 1
       endfor

endfor

end