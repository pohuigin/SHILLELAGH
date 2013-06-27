function sw_constants,no_alpha=no_alpha

paramstruct=read_param_file(!SHILLELAGH_PATH+!SHILLELAGH_PARAM, meta=parammeta)

;SW Properties fall off as 1/r^alpha
alpha_b=paramstruct.alpha_b
alpha_rho=paramstruct.alpha_rho
alpha_t=paramstruct.alpha_t
if keyword_set(no_alpha) then begin & alpha_b=0d & alpha_rho=0d & alpha_t=0d & endif

;alpha_b=2.5d ;between 2 and 3
;alpha_rho=4.d ;!!!NOT USED: switched to using empirical model in Chen 1996 (Eqn. 17.)
;alpha_t=1.d

;Physical Constants
au_km=paramstruct.au_km
vernal_equinox=paramstruct.vernal_equinox
r_sun=paramstruct.r_sun

omegasun=360d/(25.2d*3600d*24d) ;in degrees/second from diff. rot. of latitudes (-10 -> +10)
omegaearth=360d/(365d*3600d*24d) ;degrees/second
nan=0/0.

;au_km=149.6d6 ;1 AU in kilometers
;vernal_equinox = -77d ; the longitude of Capella (Aries?) in degrees or solar ascending node?? 
;r_sun=6.955d5 ;km

constants_arr=[alpha_b,alpha_rho,alpha_t,au_km,vernal_equinox,nan,r_sun,omegasun,omegaearth]

return,constants_arr

end
