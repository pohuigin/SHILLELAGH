pro run_cme_forecaster

masses=findgen(20)*.2-.1

for i=0,n_elements(masses) do begin

	print,masses[i]
	cme_forecaster,'12-dec-2008 08:52',-5.,250.,50.,cmeheight=4.,cmemass=masses[i],trange=[-1./48.,3.],tbin=1./24.,/byrne,/no_plot,out_struct=this_struct,/res1
	
	if i eq 0 then struct_arr=out_struct $
		else struct_arr=[struct_arr,struct_arr]

endfor

stop

end

note!!!! change density drop off to be 1/r or 1/r^2 times some constant negative slope?