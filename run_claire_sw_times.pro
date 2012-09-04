pro run_claire_sw_times

savpath='~/science/data/cme_propagation/claire_events/sav/'

times=['13-mar-2010 18:00:00','13-mar-2010 21:00:00','14-mar-2010 01:00:00']

sw_propagate, times, trange=0, savpath=savpath, $ ;, tbin=intbin;, no_alpha=no_alpha, , full360=full360
	spacecraft='omni', filetag='omni' ;test=test, plot_points=plot_points, 

sw_propagate, times, trange=0, savpath=savpath, $ ;, tbin=intbin;, no_alpha=no_alpha, , full360=full360
	spacecraft='sta', filetag='sta' ;test=test, plot_points=plot_points, 

sw_propagate, times, trange=0, savpath=savpath, $ ;, tbin=intbin;, no_alpha=no_alpha, , full360=full360
	spacecraft='stb', filetag='stb' ;test=test, plot_points=plot_points, 










stop

end