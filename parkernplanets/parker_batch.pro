pro parker_batch, outfile
  
  outfile = 'test'
  
  !p.charsize = 1.2

  n = 0 

  for i = 0, 30 do begin
    
    parkernplanets,anytim('1-jan-2008') + i*24.*60.*60.,/inner, vel = 700
    if ( n le 9 ) then filename = outfile + '_0' + strcompress( string( n ), /rem ) + '.jpg'
    if ( n gt 9 ) then filename = outfile + '_'  + strcompress( string( n ), /rem ) + '.jpg'
    
    print, filename
    
    n = n + 1
    
    ;x2jpeg, filename
    
  endfor
  
  ;jsmovie, outfile + '.html', findfile( outfile + '_*.jpg' )
  
end
