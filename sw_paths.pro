;-------------------------------------------------------------------------->
;Set up the paths to the in situ data files.
;OUTPUT:
;	Returned Value = the parameter requested by the keyword switch set
;
;KEYWORDS:
;	verbose = 
;	insitup = 
;	plotp = 
;	savep = 
;	plotinterpheliop = 
;	cmeplotp = 
;	cmeinputdatap = 
;	chforecast = 
;
;OUTPUT/SETUP KEYWORDS:
;	allparam = set to empty variable, will contain structure with all params
;	rootp = if set to empty variable, will contain SHILLELAgh root path
;		if set to a string, will update environment variable
;	paramf = if set to empty variable; will contain SHILLELAgh
;		param file name (NOT FULL PATH, only local name)
;		if set to a string, will update environment variable
;
;-------------------------------------------------------------------------->

function sw_paths, allparam=allparam, rootp=rootp, paramf=paramf, meta=outmeta, verbose=verbose, $
	insitup=insitup, plotp=plotp, savep=savep, plotinterpheliop=plotinterpheliop, $
	cmeplotp=cmeplotp, cmeinputdatap=cmeinputdatap, chforecast=chforecast

if keyword_set(verbose) then verbose=1 else verbose=0
result=1

if data_type(rootp) eq 7 then defsysv,'!SHILLELAGH_PATH',rootp
if data_type(paramf) eq 7 then defsysv,'!SHILLELAGH_PARAM',paramf


defsysv,'!SHILLELAGH_PATH',exist=syspath
defsysv,'!SHILLELAGH_PARAM',exist=sysparam

if not syspath then root='~/repositories/SHILLELAGH/' $
	else root=!SHILLELAGH_PATH
if not syspath then fparam=root+'sw_param.txt' $
	else fparam=root+!SHILLELAGH_PARAM

if verbose then print,'% SW_PATHS: Pointing to: ROOT = '+!SHILLELAGH_PATH
if verbose then print,'% SW_PATHS: Pointing to: PARAM = '+!SHILLELAGH_PARAM

paramstruct=read_param_file(fparam, meta=outmeta)

result=1
if keyword_set(rootp) then result=root
if keyword_set(insitup) then result=root+paramstruct.insitup
if keyword_set(plotp) then result=root+paramstruct.plotp
if keyword_set(savep) then result=root+paramstruct.savep
if keyword_set(plotinterpheliop) then result=root+paramstruct.plotinterpheliop
if keyword_set(cmeplotp) then result=root+paramstruct.cmeplotp
if keyword_set(cmeinputdatap) then result=root+paramstruct.cmeinputdatap
if keyword_set(chforecast) then result=root+paramstruct.chforecast

if keyword_set(verbose) then print,'% SW_PATHS: '+strtrim(result,2)

;Handle output keywords
allparam=paramstruct
if data_type(rootp) eq 0 then rootp=!SHILLELAGH_PATH
if data_type(paramf) eq 0 then paramf=!SHILLELAGH_PARAM


return, result

end
