;rad in Km (measured at modelled spiral point)
;ref_rad in Km (measured at spacecraft)
;ref_prop in:
;	density in n/cm^-3
;	temperature in K
;	magneticfield in G??
function calc_helio_props, inrad, inref_rad, inref_prop, constants=constants_arr, $
	density=density, temperature=temperature, magneticfield=magneticfield, velocity=velocity

;Load physical constants
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]

;Put distances in units of Rsun
rad=double(inrad)/r_sun
ref_rad=double(inref_rad)/r_sun
ref_prop=double(inref_prop)

;density in n/cm^-3. equation 17. from Chen 1996 (needs distance in Rsun)
if keyword_set(density) then begin
	rho_chen_rad=4d*(3d*rad^(-12d)+rad^(-4d))*1d8+3.5d*1d5*rad^(-2d)
	rho_chen_refrad=4d*(3d*ref_rad^(-12d)+ref_rad^(-4d))*1d8+3.5d*1d5*ref_rad^(-2d)
	retval=rho_chen_rad*(inref_prop/rho_chen_refrad)
endif

if keyword_set(temperature) then begin



endif

if keyword_set(magneticfield) then begin



endif





return, retval

end