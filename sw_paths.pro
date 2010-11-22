;-------------------------------------------------------------------------->
;Set up the paths to the in situ data files.
function sw_paths, insitup=insitup, plotp=plotp, savep=savep, plotinterpheliop=plotinterpheliop, 
	cmeplotp=cmeplotp, cmeinputdatap=cmeinputdatap, $
	rootp=rootp, verbose=verbose

result=''

root='~/science/data/cme_propagation/'
if keyword_set(rootp) then result=root

if keyword_set(insitup) then result=root+'sw_prop_insitu/'
if keyword_set(plotp) then result=root+'sw_prop_plots/'
if keyword_set(savep) then result=root+'sw_prop_save/'
if keyword_set(plotinterpheliop) then result=root+'sw_prop_interp_plots/'
if keyword_set(cmeplotp) then result=root+'plots/'
if keyword_set(cmeinputdatap) then result=root+'orbits/'

if keyword_set(verbose) then print,'SW_PATHS: '+result

return, result

end
