pro test_helio

dates=anytim(file2time(datearr('1-jan-2008','31-dec-2008',/vms)))
earth_jd_struct = anytim2jd(dates)
earth_jd = earth_jd_struct.int + earth_jd_struct.frac
helio, earth_jd, 3, earth_rad, earth_lon, earth_lat

setplotenv,/xwin
setcolors,/sys
!p.background=255
!p.color=0

plot,(anytim(dates)-min(anytim(dates)))/3600./24.,earth_lon,xtit='Day of Year 2008',ytit='Longitude of Earth'

;oplot,(anytim(dates)-min(anytim(dates)))/3600./24.,earth_lon+180.+79.,ps=4
oplot,(anytim(dates)-min(anytim(dates)))/3600./24.,earth_lon-77.,ps=4

readcol,'~/science/procedures/cme_propagation/orbits/omni_2008.txt',omndd,omntt,omnlon,form='A,A,F',delim=' '
omntim=anytim(strmid(omndd,6,4)+'-'+strmid(omndd,3,2)+'-'+strmid(omndd,0,2)+'T'+omntt)
;omntim=anytim(strtrim(omndd,2)+'T'+strtrim(omntt,2))

oplot,(omntim-min(anytim(dates)))/24./3600.,omnlon,color=!red

vline,(anytim('21-mar-2008')-min(anytim(dates)))/3600./24.,color=!blue

legend,['helios (HAE)','helios-77 deg.','heliographic inertial (HGI)','Vernal Equinox'],psym=[-3,-4,-3,-3],color=[!black,!black,!red,!blue],/left,/top

window_capture,/png,file='~/science/procedures/cme_propagation/coordinate_system_earth_pos_comparison'

stop

window,1
plot,(anytim(dates)-min(anytim(dates)))/3600./24.,earth_lat


stop

end