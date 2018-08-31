;;;;; Put MLS T data on WACCM grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

MLS_path='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
output_path='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/'

;;; setup WACCM grid to GEOS grid
;LONGITUDE (lon)
;double[144]
lon=dindgen(96L)*3.75 + 1.875
nlon=n_elements(lon) & lonstep=(lon[1:*]-lon[0:nlon-2])[nlon/2]

;LATITUDE (lat)
;double[96]
lat=dindgen(72L)*2.5 -88.75
nlat=n_elements(lat) & latstep=(lat[1:*]-lat[0:nlat-2])[nlat/2]



;ALTITUDE (lev)
;double[88]

nlev=121L    
altitude = fltarr(121L)


   		    

;;; time range
start_date=20130301
end_date=20130301
;read, 'startdate, yyyymmdd:',start_date
;read, 'enddate, yyyymmdd:', end_date

start_date = strtrim(strcompress(start_date),2)
end_date = strtrim(strcompress(end_date),2)

ndays=julday(strmid(end_date,4,2),strmid(end_date,6,2),strmid(end_date,0,4)) - $
      julday(strmid(start_date,4,2),strmid(start_date,6,2),strmid(start_date,0,4)) +1
current_date = start_date
start_julday = julday(strmid(start_date,4,2),strmid(start_date,6,2),strmid(start_date,0,4))

;;; day loop
for iday=0,ndays-1 do begin
;for iday=0,1 do begin
  s0=systime(1)


 ;prepare date variables 
 
 current_julday = start_julday + iday 
 caldat, current_julday, month, day, year 
 current_doy = current_julday-julday(1,1,year)+1 
 str_year = strcompress(string(year),/r) 
 str_date_code = str_year+'d'+strcompress(string(current_doy,format='(i3.3)'),/r) 
 print, 'Starting '+str_date_code+' !!!'
  
 ;read MLS T  file=file_search(MLS_path+'cat_mls_v3.3_'+strcompress(strtrim(string(year),2))+strcompress(string(month,format='(I02)'))$
 +strcompress(string(day,format='(I02)')) + '.sav',count=nfiles)  
 if nfiles eq 1 then restore,file  
 file=file_search(MLS_path+'tpd_mls_v3.3_'+strcompress(strtrim(string(year),2))+strcompress(string(month,format='(I02)'))$
 +strcompress(string(day,format='(I02)')) + '.sav',count=nfiles)  
 if nfiles eq 1 then restore,file  
 if nfiles eq 0 then begin
    print, '   NO DATA ON '+str_date_code
    continue
  endif
  
  
  print, '   MLS data read.'


 doy = current_doy 
    pi=3.14159265
    dtor=pi/180.
    earinc=23.5
    zangle=fltarr(n_elements(latitude))
    for ii=0L,n_elements(latitude)-1 do begin
        rlat=latitude(ii)
        rlon=longitude(ii)
        gmt=time(ii)
        sinlat=sin(rlat*dtor)
        coslat=sqrt(1.-sinlat^2.)
        sinlon=sin(rlon*dtor)
        coslon=cos(rlon*dtor)
        soya=(doy-81.25)*pi/182.5           ; day angle
        soha=2.*pi*(gmt-12.)/24.            ; hour angle
        soha=-soha
        sininc=sin(earinc*dtor)
        sindec=sininc*sin(soya)
        cosdec= sqrt(1.-sindec^2.)
        coszen=cos(soha)*coslon+sin(soha)*sinlon
        coszen=coszen*cosdec*coslat
        coszen=sindec*sinlat+coszen
        coszen=min([max([coszen,-1.]),1.])
        chi = acos(coszen)
        zangle(ii) = chi/dtor
    endfor

  
  
    ;screen T
  
  
  ;store T
  MLS_t=temperature
  ntimes = n_elements(time)
  MLS_lon=longitude
  x = where(MLS_lon lt 0., nx)
  if nx gt 0 then MLS_lon[x]=360.+MLS_lon[where(MLS_lon lt 0)]
  index=where(MLS_lon lt 0 or MLS_lon gt 360,nindex)
  if nindex ne 0 then MLS_lon[index]=!values.f_nan
  MLS_lat=latitude
  index=where(MLS_lat lt -90 or MLS_lat gt 90,nindex)
  if nindex ne 0 then MLS_lat[index]=!values.f_nan
  
  xx = where(temperature lt 0. or temperature_mask lt 0.,nxx)
  if nxx gt 0 then temperature[xx] = !values.f_nan
MLS_temperature = fltarr(n_elements(lon), n_elements(lat), n_elements(altitude))*!values.f_nan

  ;store P
  MLS_p=pressure
  ntimes = n_elements(time)
  MLS_lon=longitude
  x = where(MLS_lon lt 0., nx)
  if nx gt 0 then MLS_lon[x]=360.+MLS_lon[where(MLS_lon lt 0)]
  index=where(MLS_lon lt 0 or MLS_lon gt 360,nindex)
  if nindex ne 0 then MLS_lon[index]=!values.f_nan
  MLS_lat=latitude
  index=where(MLS_lat lt -90 or MLS_lat gt 90,nindex)
  if nindex ne 0 then MLS_lat[index]=!values.f_nan
  
  xx = where(pressure lt 0.,nxx)
  if nxx gt 0 then pressure[xx] = !values.f_nan
MLS_pressure = fltarr(n_elements(lon), n_elements(lat), n_elements(altitude))*!values.f_nan
  
  

    ;check for coordinates out of range
  index=where(MLS_lon lt 0 or MLS_lon gt 360,nindex)
  if nindex ne 0 then MLS_lon[index]=!values.f_nan
  index=where(MLS_lat lt -90 or MLS_lat gt 90,nindex)
  if nindex ne 0 then MLS_lat[index]=!values.f_nan
  index=where(finite(MLS_lon) ne 1 or finite(MLS_lat) ne 1,nindex)
  
;  ;produce sunlight flag
;  earth_radius=6367500. ;in m
;  sza_sunset=180.-180./!pi*asin(earth_radius/(earth_radius+gph))
;  sunlight=transpose(rebin(sza,n_elements(sza),nlev)) lt sza_sunset ;1 = sunlight, 0 = darkness

  ;produce ascending/descending flag
  ;lat_gradient=convol(MLS_lat,[0,1,-9],/normalize,/center,/nan,/edge_truncate)
mode = findgen(n_elements(time))*0L + 1L  
node=mode
  
  ;gridding
  ;MLS_lon=MLS_lon*(MLS_lon le 180)+(MLS_lon-360.)*(MLS_lon gt 180) ;transform from 0/360 to -180/+180
   t_grid = fltarr(nlon,nlat,nlev)*!values.f_nan
   p_grid = fltarr(nlon,nlat,nlev)*!values.f_nan
   zenith_grid = fltarr(nlon,nlat)*!values.f_nan
   tmp_grid_node=fltarr(nlon,nlat,nlev,2)*!values.f_nan
   pres_grid_node=fltarr(nlon,nlat,nlev,2)*!values.f_nan
  zen_grid_node=fltarr(nlon,nlat,2)*!values.f_nan

  
  
    for ilev=0,nlev-1 do begin
    for inode=0,1 do begin
      gindex=where(finite(reform(MLS_lon[*])) and finite(reform(MLS_lat[* ])) and finite(reform(temperature[*,ilev])) and node eq inode,ngindex)
      if ngindex lt 24 then continue ;only grid when there are at least 24 points for the entire sphere
      grid_input, MLS_lon[gindex ], MLS_lat[gindex ], reform(temperature[gindex,ilev]), $
                  prep_cartesian, prep_values, $
                  /sphere,/degrees,duplicates='Avg',epsilon=.5
                  
      prep_spherical = cv_coord(/degrees, /double, from_rect=prep_cartesian, /to_sphere)
      prep_lon = reform(prep_spherical[0,*]*(prep_spherical[0,*] gt 0)+(360+prep_spherical[0,*])*(prep_spherical[0,*] lt 0)) ;transform from -180/+180 to 0/360)            
      prep_lat = reform(prep_spherical[1,*])
      
;      prep_lon=MLS_lon[gindex]
;      prep_lat=MLS_lat[gindex]
;      prep_values=reform(t_mix[ilev,gindex])
      
      qhull, prep_lon, prep_lat, tr, /delaunay, sphere=s
      
      case inode of
         0: tmp_grid_node[*,*,ilev,0]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         1: tmp_grid_node[*,*,ilev,1]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         else: stop
      endcase
  
  ;Press    
          gindex=where(finite(reform(MLS_lon[*])) and finite(reform(MLS_lat[* ])) and finite(reform(pressure[*,ilev])) and node eq inode,ngindex)
      if ngindex lt 24 then continue ;only grid when there are at least 24 points for the entire sphere
      grid_input, MLS_lon[gindex ], MLS_lat[gindex ], reform(pressure[gindex,ilev]), $
                  prep_cartesian, prep_values, $
                  /sphere,/degrees,duplicates='Avg',epsilon=.5
                  
      prep_spherical = cv_coord(/degrees, /double, from_rect=prep_cartesian, /to_sphere)
      prep_lon = reform(prep_spherical[0,*]*(prep_spherical[0,*] gt 0)+(360+prep_spherical[0,*])*(prep_spherical[0,*] lt 0)) ;transform from -180/+180 to 0/360)            
      prep_lat = reform(prep_spherical[1,*])
      
;      prep_lon=MLS_lon[gindex]
;      prep_lat=MLS_lat[gindex]
;      prep_values=reform(t_mix[ilev,gindex])
      
      qhull, prep_lon, prep_lat, tr, /delaunay, sphere=s
      
      case inode of
         0: pres_grid_node[*,*,ilev,0]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         1: pres_grid_node[*,*,ilev,1]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         else: stop
      endcase
  
      
      
      
      
      
      if ilev eq 0L then begin
          gindex=where(finite(reform(MLS_lon[*])) and finite(reform(MLS_lat[* ])) and finite(reform(zangle[*])) and node eq inode,ngindex)
      if ngindex lt 24 then continue ;only grid when there are at least 24 points for the entire sphere
      grid_input, MLS_lon[gindex ], MLS_lat[gindex ], reform(zangle[gindex]), $
                  prep_cartesian, prep_values, $
                  /sphere,/degrees,duplicates='Avg',epsilon=.5
                  
      prep_spherical = cv_coord(/degrees, /double, from_rect=prep_cartesian, /to_sphere)
      prep_lon = reform(prep_spherical[0,*]*(prep_spherical[0,*] gt 0)+(360+prep_spherical[0,*])*(prep_spherical[0,*] lt 0)) ;transform from -180/+180 to 0/360)            
      prep_lat = reform(prep_spherical[1,*])
      
;      prep_lon=MLS_lon[gindex]
;      prep_lat=MLS_lat[gindex]
;      prep_values=reform(t_mix[ilev,gindex])
      
      qhull, prep_lon, prep_lat, tr, /delaunay, sphere=s

      case inode of
         0: zen_grid_node[*,*,0]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         1: zen_grid_node[*,*,1]=reform(griddata(prep_lon, prep_lat, prep_values, /degrees, /grid, xout=lon, yout=lat, triangles=tr, /natural_neighbor, /sphere))
         else: stop
      endcase
      endif; ilev=0      
    endfor ;end of inode
  endfor ;end of ilev
  print, '   MLS data gridded.'
        
  ;smoothing
  for ilev=0,nlev-1 do begin
    for inode=0,1 do begin
      big_array=fltarr(nlon*2,nlat*3-2) ;minus two to not double up the 90S and 90N latitudes in the middle of big_array
      big_arrayz=fltarr(nlon*2,nlat*3-2) ;minus two to not double up the 90S and 90N latitudes in the middle of big_array
      big_arrayp=fltarr(nlon*2,nlat*3-2) ;minus two to not double up the 90S and 90N latitudes in the middle of big_array
      case inode of
        0: tmp_grid=reform(tmp_grid_node[*,*,ilev,0])
         1: tmp_grid=reform(tmp_grid_node[*,*,ilev,1])
         else: stop
      endcase
      case inode of
        0: pres_grid=reform(pres_grid_node[*,*,ilev,0])
         1: pres_grid=reform(pres_grid_node[*,*,ilev,1])
         else: stop
      endcase
      if ilev eq 0L then begin
      case inode of
        0: zen_grid=reform(zen_grid_node[*,*,0])
         1: zen_grid=reform(zen_grid_node[*,*,1])
         else: stop
       endcase
      endif
      for irow=0,1 do begin
        big_array[irow*nlon:(irow+1)*nlon-1,2*nlat-2:*]=rotate(reform(tmp_grid),7) ;top row
        big_array[irow*nlon:(irow+1)*nlon-1,0:nlat-1]=rotate(reform(tmp_grid),7) ;bottom row
        big_arrayp[irow*nlon:(irow+1)*nlon-1,2*nlat-2:*]=rotate(reform(pres_grid),7) ;top row
        big_arrayp[irow*nlon:(irow+1)*nlon-1,0:nlat-1]=rotate(reform(pres_grid),7) ;bottom row
      if ilev eq 0L then        big_arrayz[irow*nlon:(irow+1)*nlon-1,2*nlat-2:*]=rotate(reform(zen_grid),7) ;top row
      if ilev eq 0L then         big_arrayz[irow*nlon:(irow+1)*nlon-1,0:nlat-1]=rotate(reform(zen_grid),7) ;bottom row
 
      endfor
      big_array[0:nlon/2-1,nlat-1:2*nlat-2]=reform(tmp_grid[nlon/2:*,*]) ;middle row, left
      big_array[nlon/2:nlon+nlon/2-1,nlat-1:2*nlat-2]=reform(tmp_grid) ;middle row, center
      big_array[nlon+nlon/2:*,nlat-1:2*nlat-2]=reform(tmp_grid[0:nlon/2-1,*]) ;middle row, right
      big_arrayp[0:nlon/2-1,nlat-1:2*nlat-2]=reform(pres_grid[nlon/2:*,*]) ;middle row, left
      big_arrayp[nlon/2:nlon+nlon/2-1,nlat-1:2*nlat-2]=reform(pres_grid) ;middle row, center
      big_arrayp[nlon+nlon/2:*,nlat-1:2*nlat-2]=reform(pres_grid[0:nlon/2-1,*]) ;middle row, right
            if ilev eq 0L then big_arrayz[0:nlon/2-1,nlat-1:2*nlat-2]=reform(zen_grid[nlon/2:*,*]) ;middle row, left
            if ilev eq 0L then big_arrayz[nlon/2:nlon+nlon/2-1,nlat-1:2*nlat-2]=reform(zen_grid) ;middle row, center
           if ilev eq 0L then big_arrayz[nlon+nlon/2:*,nlat-1:2*nlat-2]=reform(zen_grid[0:nlon/2-1,*]) ;middle row, right
      
      ;big_array properties
      np_ilat=2*nlat-2
      sp_ilat=nlat-1
      lat3=[lat,lat[1:nlat-1],lat[1:*]]
      lon2=[lon,lon]
      
      ;smoothing
      big_array_new=big_array
      big_arrayp_new=big_arrayp
      if ilev eq 0L then      big_arrayz_new=big_arrayz
      nlat_width=3
      nlat_width+=((nlat_width+1) mod 2) ;force nlat_width to be odd (increase by 1 if necessary)
      lat_width=(nlat_width-1)*latstep
      for ilat=nlat_width/2, n_elements(lat3)-1-nlat_width do begin ;convolution changes with latitude, so lat by lat convolution
        ;determine nlon width of kernel
        nlon_width=intarr(nlat_width)
        for ilat_width=-nlat_width/2,nlat_width/2 do nlon_width[ilat_width+nlat_width/2]=round(nlat_width*latstep/lonstep/cos(!dtor*(lat3[ilat+ilat_width]))) < (nlon-1)
        nlon_width=max(nlon_width)
        nlon_width+=((nlon_width+1) mod 2) ;force nlon_width to be odd (increase by 1 if necessary)
      
        ;create kernel with distance weights
        kernel_smooth=fltarr(nlon_width,nlat_width)
        for iklat=-nlat_width/2,nlat_width/2 do begin
          for iklon=-nlon_width/2,nlon_width/2 do begin
            ;fill kernel width distance weights to center of kernel (map_2points)
            p1_lon=lonstep*iklon ;arbitrary longitude since only the difference in longitude matters 
            p1_lat=lat3[ilat]+latstep*iklat
            if abs(p1_lat) gt 90 then begin
              p1_lat=(2*(p1_lat gt 0)-1)*180-p1_lat
              ;if p1_lat gt 90 then p1_lat=180-p1_lat
              ;if p1_lat lt -90 then p1_lat=-180-p1_lat
              p1_lon=(p1_lon+180) mod 360
            end
            if p1_lon gt 180 then p1_lon=360-p1_lon ;make the longitude -180|180
            if p1_lon lt -180 then p1_lon=360+p1_lon ;make the longitude -180|180
            kernel_smooth[nlon_width/2+iklon,nlat_width/2+iklat]=map_2points(0,lat3[ilat],p1_lon,p1_lat,/meters,radius=1)
          endfor ;end of ilon
        endfor ;end of iklat
        kernel_smooth=-1.*(kernel_smooth-max(kernel_smooth)) ;kernel is not normalized (will be done by convol)
        
        ;apply smoothing
        tmp=convol(reform(big_array[*,(ilat-nlat_width/2):(ilat+nlat_width/2)]),kernel_smooth,/center,/normalize,/edge_truncate,/nan)
        big_array_new[*,ilat]=reform(tmp[*,nlat_width/2])
        pres=convol(reform(big_arrayp[*,(ilat-nlat_width/2):(ilat+nlat_width/2)]),kernel_smooth,/center,/normalize,/edge_truncate,/nan)
        big_arrayp_new[*,ilat]=reform(pres[*,nlat_width/2])
      if ilev eq 0L then        zen=convol(reform(big_arrayz[*,(ilat-nlat_width/2):(ilat+nlat_width/2)]),kernel_smooth,/center,/normalize,/edge_truncate,/nan)
      if ilev eq 0L then   big_arrayz_new[*,ilat]=reform(zen[*,nlat_width/2])
        
      endfor ;end of ilat
      big_array=big_array_new
      big_arrayp=big_arrayp_new
      case inode of
         0: tmp_grid_dsc=reform(big_array[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         1: tmp_grid_asc=reform(big_array[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         else: stop
      endcase 
      case inode of
         0: pres_grid_dsc=reform(big_arrayp[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         1: pres_grid_asc=reform(big_arrayp[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         else: stop
      endcase 
            if ilev eq 0L then begin
      case inode of
         0: zen_grid_dsc=reform(big_arrayz[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         1: zen_grid_asc=reform(big_arrayz[nlon/2:nlon+nlon/2-1,sp_ilat:np_ilat])
         else: stop
      endcase 
		endif
    endfor ;end of inode loop
    t_grid[*,*,ilev]=(tmp_grid_asc)
    p_grid[*,*,ilev]=(pres_grid_asc)
      if ilev eq 0L then    zenith_grid[*,*]=(zen_grid_asc)
  endfor ;end of ilev
  print, '   Gridded MLS data smoothed.'
  

  
  
  
;geometricAltitude = t_grid*0.    
;rtd=double(180./!pi)
;dtr=1./rtd
;ks=1.931853d-3
;ecc=0.081819
;gamma45=9.80
;for jj=0L,n_elements(lat)-1L do begin
;    for kk=0L,n_elements(lev)-1L do begin
;    
;            sin2=sin( (lat(jj)*dtr)^2.0 )
;        numerator=1.0+ks*sin2
;        denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
;        gammas=gamma45*(numerator/denominator)
;        r=6378.137/(1.006803-(0.006706*sin2))
;        geometricAltitude(*,jj,kk)=(r*gpaltitude[*,jj,kk])/ ( (gammas/gamma45)*r - gpaltitude[*,jj,kk] )
;  endfor
;endfor
 


  
;for ii = 0, n_elements(lon) - 1L do begin
;for jj = 0, n_elements(lat) - 1L do begin
;MLS_temperature[ii,jj,*] = interpol( t_grid[ii,jj,*], reform(geometricAltitude[ii,jj,*]),altitude[*])
			  
    
;endfor
;endfor

  ;save daily file

temperature = t_grid

date = current_date

save, filename = output_path+'MLS_T_grid5_v3.3_'+strcompress(strtrim(string(year),2))+strcompress(string(month,format='(I02)'))$
+strcompress(string(day,format='(I02)')) + '.sav', lon, lat, t_grid, P_grid, altitude


print,'   Gridded MLS data saved.'
  
print, 'Finished '+str_date_code+' in '+strcompress(string(systime(1)-s0,format='(i4)'),/r)+'s !!!' 
endfor ;end of iday loop

print,'...done.'
end
