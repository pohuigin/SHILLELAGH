;-------------------------------------------------------------------------->

function sw_get_data, trange, omni=omni, sta=sta, stb=stb, ace=ace, wind=wind, $
	constants=constants_arr, geoindex=geoindex

orbitp=sw_paths(/insitu)
;constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun]
if n_elements(constants_arr) lt 1 then constants_arr=sw_constants()
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
	sc_br=omn_br
	sc_bt=omn_bt
	sc_bn=omn_bn
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
	sc_br=sta_br
	sc_bt=sta_bt
	sc_bn=sta_bn
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
	sc_br=stb_br
	sc_bt=stb_bt
	sc_bn=stb_bn
endif

;TIME_AT_CENTER_OF_HOUR  1AU_IP_FLOW_PRESSURE 1AU_IP_ELECTRIC_FIELD 1AU_IP_PLASMA_BETA 1AU_IP_ALFVEN_MACH_NO. 3-H_KP*10 1-H_DST 1-H_AE 3-H_AP 1-H_AL-INDEX AU-INDEX 1-H_PC(N)-INDEX
if keyword_set(geoindex) then begin
	readcol,orbitp+'omni_index_'+strmid(time2file(trange[0],/date),0,4)+'.txt',omn_dd,omn_tt,omn_press,omn_efield,omn_beta,omn_mach,omn_kp,omn_dst,omn_ae,omn_ap,omn_al,omn_au,omn_pc,form='A,A,F,F,F,F,F,F,F,F,F,F,F',delim=' '
	omn_tim=anytim(strmid(omn_dd,6,4)+'-'+strmid(omn_dd,3,2)+'-'+strmid(omn_dd,0,2)+'T'+omn_tt)
	sc_arr=fltarr(8,n_elements(omn_tim))
	sc_arr[0,*]=omn_tim
	sc_arr[1,*]=omn_kp
	sc_arr[2,*]=omn_dst
	sc_arr[3,*]=omn_ae
	sc_arr[4,*]=omn_ap
	sc_arr[5,*]=omn_al
	sc_arr[6,*]=omn_au
	sc_arr[7,*]=omn_pc
	
	if (where(sc_arr[3,*] eq -1.*10^31. or sc_arr[6,*] eq 999.9))[0] ne -1 then sc_arr[*,where(sc_arr[3,*] eq -1.*10^31. or sc_arr[6,*] eq 999.9)]=nan
	sc_arr=sc_arr[*,where(finite(sc_arr[3,*]) eq 1)]
	return,sc_arr
endif

;make an array of S=([r,theta,v,rho,temp,bmag,t],nblob)
sc_arr=fltarr(10,n_elements(sc_vel_arr))
sc_arr[0,*]=sc_r_arr
sc_arr[1,*]=sc_hgtheta_arr
sc_arr[2,*]=sc_vel_arr
sc_arr[3,*]=sc_rho
sc_arr[4,*]=sc_temp
sc_arr[5,*]=sc_bmag
sc_arr[6,*]=sc_tim_arr
sc_arr[7,*]=sc_br
sc_arr[8,*]=sc_bt
sc_arr[9,*]=sc_bn

;Check for bad data points
if (where(sc_arr[3,*] eq -1.*10^31.))[0] ne -1 then sc_arr[*,where(sc_arr[3,*] eq -1.*10^31.)]=nan
sc_arr=sc_arr[*,where(finite(sc_arr[3,*]) eq 1)]

return,sc_arr

end

