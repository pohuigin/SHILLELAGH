;-------------------------------------------------------------------------->
;Make sure that longitudes run from 0-360
function sw_theta_shift, inarr

outarr=inarr
if (where(inarr ge 360))[0] ne -1 then outarr[where(inarr ge 360)]=outarr[where(inarr ge 360.)]-360.
if (where(inarr lt 0))[0] ne -1 then outarr[where(inarr lt 0)]=outarr[where(inarr lt 0.)]+360.

return, outarr

end

