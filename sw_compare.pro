;------------------------------------------------------------------------------>
;Compare two detected features to see if they are connected spatially.
;Returns 0 or 1, and optionally a structure with more detailed information on the connection.
;
;HG1 and HG2: for a surface feature this is in heliographic [Longitude, Latitude] coordinates 
;	(HG). Earth-Sun line is at longitude 0, and solar equator is at latitude 0. For 
;	heliospheric coordinates, this is in heliographic inertial [Radius, Longitude] coordinates
;	(HGI) so that 0 longitude is at the solar rising node(??) and radius is from the Sun.
;
;TYPE: the type of comparison to be made. A string with the form, 
;	"coordtype1-coordtype2;surfaceorhelio1-surfaceorhelio2;feature1-feature2"
;	There must be two things to compare: two sets of generic coordinates, coordinates and a 
;	named feature (CME/SEP/CH), or two named features.. etc.
;Examples of TYPE:
;	Comparing a flare location to an AR boundary: "point-poly;surface-surface;gen-gen"
;	Comparing a CH boundary to a HSSW stream detected in situ: "poly-struct;surface-helio;gen-hssw"
;	Comparing a planet to a propagated CME: "point-struct;helio-helio;gen-cme"
;
;TOLERANCE: how close the features have to match. For surface this is [degrees, degrees] and for
;	heliospheric this is [AU, degrees]
;
;Output: 	1 	- The features are connected in some way.
;			0 	- The features do not seem to be connected.
;			-1	- An error occurred. 
;
;OUTSTRUCT: a structure with specific information on the nature of the matching. 
;	May include: duration of association, time of arrival, location when associated, etc.
;	An info keyword tells about in what coordinates they match, and closeness of match.
;------------------------------------------------------------------------------>

function sw_compare_point(hg1, hg2, tolerance=tolerance, location=typloc)
;;;
end

;------------------------------------------------------------------------------>

function sw_compare, hg1=hg1, hg2=hg2, poly1=poly1, poly2=poly2, cmestruct=cmestruct, sepstruct=sepstruct, hsswstruct=hsswstruct, $
	type=type, tolerance=tolerance, outstruct=outstruct

if var_type(type) ne 7 then begin
	print,'TYPE string must be supplied.'
	return,-1
endif

typarr=str_sep(type,';')
typcoord=str_sep(typarr[0],'-') ;point or polygon coordinates
typloc=str_sep(typarr[1],'-') ;Surface or Heliospheric
typfeat=str_sep(typarr[2],'-') ;generic, cme, hssw, or sep

;Case of two generic point like features-------------->

if n_elements(where(typcoord eq point)) eq 2 then begin
	outishit=sw_compare_point(hg1, hg2, tolerance=tolerance, location=typloc)

endif

;----------------------------------------------------->














ishit=outishit

return, ishit

end