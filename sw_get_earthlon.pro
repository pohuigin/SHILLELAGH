function sw_get_earthlon, thistime

constants_arr=sw_constants()
alpha_b=constants_arr[0]
alpha_rho=constants_arr[1]
alpha_t=constants_arr[2]
au_km=constants_arr[3]
vernal_equinox=constants_arr[4]
nan=constants_arr[5]
r_sun=constants_arr[6]
omegasun=constants_arr[7]
omegaearth=constants_arr[8]

earth_jd_struct = anytim2jd( anytim(thistime,/vms) )
earth_jd = earth_jd_struct.int + earth_jd_struct.frac
helio, earth_jd, 3, earth_rad, earth_lon, earth_lat
earth_lon = earth_lon + vernal_equinox
earth_lon=sw_theta_shift(earth_lon)
this_earthlon=earth_lon[0]; & earth_lon=earth_lon[1:*]


return,this_earthlon

end
