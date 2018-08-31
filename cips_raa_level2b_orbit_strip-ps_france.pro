; Laboratory for Atmospheric and Space Physics
; University of Colorado, Boulder, Colorado, USA
;
; FILENAME:
;   CIPS_raa_orbit_strip_plot.pro
;
; AUTHOR:
;    Jeff France
;
; DATE:    Feb 19, 2018
;
; PURPOSE:
;   This program reads in Level 2b RAA files and creates orbit strip plots and corresponding thumbnails
;  
; NOTES:
;    

;Identify orbits to loop over
startorbit = 58000
endorbit = 59779

RAA_Level2_path='/atmos/harvey/CIPS_data/Datfiles/RAA/'
outpath = '/Users/franceja/CIPS_GW_paper/Figures/'


;Read in MERRA2 file for determining land/water using surface GPHT
file='/Users/franceja/CIPS_GW_paper/MERRA2_300.inst3_3d_asm_Np.20040101.SUB.nc4'
ncfile = file_search(file,count = count);Open the *.nc file
ncid=ncdf_open(ncfile)
result=ncdf_inquire(ncid)   ;Inquire about the data
nvars=result.nvars        ;# variables in the file
;Read in the data
for ivar=0,nvars-1 do begin
		result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
		;Puts data into array called "data" and variable name into "result.name":
		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
		
	 ; Invoke some magic to check for string data
	 ; masquerading as byte data, but don't convert
	 ; byte data blindly, i.e., quality_flags is a 2-dimensional
	 ; array of byte data. 
	 ; (This is done to account for a bug in the ncdf write routine.)
		  if ( ( size( data, /n_dimensions )  EQ 1 ) && $
			( size( data, /type ) EQ 1 ) ) then $
			data = string( data )

	 ;Extract each variable from the "data" structure and name it
	;  the corresponding "name" from "result.name":
		   if Execute(result.name + ' = data') eq 0 then $
		   Print, ' "Execute" command failed -- are you in Virtual Machine mode?'            
endfor
NCDF_CLOSE, ncid

;Assume surface GPHT less than 1m and greater than -1m is water
x = where(abs(mean(phis ,dim=3)) gt .1)
continent_mask = reform(phis[*,*,0])*0.
continent_mask[x] = 1.
continent_mask = smooth(continent_mask,5,/nan)
lon_mask = lon
x = where(lon_mask lt 0.)
lon_mask[x] = lon_mask[x] + 360.
lat_mask = lat


a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill


orbitList = findgen(endorbit-startorbit)+startorbit
strOrbits = string(orbitlist, format = '(i5.5)')
;LOOP OVER ORBITS
for iorbit = 0, n_elements(orbitlist) - 1L do begin ; loop over orbits

SETPLOT='ps'
;read,'setplot',setplot
nxdim=3000*1
nydim=450*1
xsize=nxdim/100.
ysize=nydim/100.
set_plot,'PS'
device,/inch,xoff=0,yoff=0.,xsize=xsize,ysize=ysize,/encapsulate,$
	/bold,/color,bits_per_pixel=8,/helvetica,filename='test2_'+strorbits[iorbit]+'.ps',decompose =0
!p.charsize=1.4    ; test2 with 1.5
!p.thick=2
!p.charthick=5
!p.charthick=5
!y.thick=2
!x.thick=2
!p.font=0

erase
loadct,0
	print, strorbits[iorbit]
	;---------------------READ IN CIPS DATA USING READ_CIPS_NCDF-------------------------------------------
	file=RAA_Level2_path+'cips_raa_2b_orbit_'+strOrbits[iorbit]+'*v01.10_r03_alb.nc'
	ncfile = file_search(file,count = count);Open the *.nc file
	if count ne 1 then continue
	 ncid=ncdf_open(ncfile)
	 result=ncdf_inquire(ncid)   ;Inquire about the data
	 nvars=result.nvars        ;# variables in the file
	 ;Read in the data
	for ivar=0,nvars-1 do begin
			result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
			;Puts data into array called "data" and variable name into "result.name":
			ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
			
		 ; Invoke some magic to check for string data
		 ; masquerading as byte data, but don't convert
		 ; byte data blindly, i.e., quality_flags is a 2-dimensional
		 ; array of byte data. 
		 ; (This is done to account for a bug in the ncdf write routine.)
			  if ( ( size( data, /n_dimensions )  EQ 1 ) && $
				( size( data, /type ) EQ 1 ) ) then $
				data = string( data )

		 ;Extract each variable from the "data" structure and name it
		;  the corresponding "name" from "result.name":
			   if Execute(result.name + ' = data') eq 0 then $
			   Print, ' "Execute" command failed -- are you in Virtual Machine mode?'            
	endfor
	NCDF_CLOSE, ncid

	file=RAA_Level2_path+'cips_raa_2b_orbit_'+strOrbits[iorbit]+'*v01.10_r03_cat.nc'
	ncfile = file_search(file,count = count);Open the *.nc file
	if count ne 1 then continue
	ncid=ncdf_open(ncfile)
	result=ncdf_inquire(ncid)   ;Inquire about the data
	nvars=result.nvars        ;# variables in the file
	;Read in the data
	for ivar=0,nvars-1 do begin
		result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
		;Puts data into array called "data" and variable name into "result.name":
		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
		
		; Invoke some magic to check for string data
		; masquerading as byte data, but don't convert
		; byte data blindly, i.e., quality_flags is a 2-dimensional
		; array of byte data.
		; (This is done to account for a bug in the ncdf write routine.)
		if ( ( size( data, /n_dimensions )  EQ 1 ) && $
			( size( data, /type ) EQ 1 ) ) then $
			data = string( data )
			
		;Extract each variable from the "data" structure and name it
		;  the corresponding "name" from "result.name":
		if Execute(result.name + ' = data') eq 0 then $
			Print, ' "Execute" command failed -- are you in Virtual Machine mode?'
	endfor
	NCDF_CLOSE, ncid
	
	;-----------------------------------------------------------------------




	;Filter out high SZA
	x=where(ZENITH_ANGLE ge 90. or finite(RAYLEIGH_ALBEDO_ANOMALY) eq 0., nx)
	RAYLEIGH_ALBEDO_ANOMALY[x] = !values.f_nan

	;redefine plot variables and shift longitude to 0-360
	plot_RAA = RAYLEIGH_ALBEDO_ANOMALY
	plot_lon = longitude
  x = where(plot_lon lt 0.)
  plot_lon[x] = plot_lon[x] +360.
  x = where(longitude lt 0.)
  longitude[x] = longitude[x] +360.
	plot_lat = latitude

	
	;Determine orbit start time in UT hours:minutes
	x=where(finite(ut_time))                  
	UT = ut_time[x[0]] - floor(ut_time[x[0]])
	hours = strmid(ORBIT_START_TIME_UT,9,2);  STRING    = '2017/266-06:59:44'
	minutes = strmid(ORBIT_START_TIME_UT,12,2)
	
	;Set file names and plot titles 
 	plot_title = 'CIPS RAA (%) -- '+' Orbit Number = '+strOrbits[iorbit]+', Date = '+UT_DATE_ORBIT_START+' '+hours+':'+minutes+'Z'
	RAA_Plot_filename = outpath+'RAA_Level2b_v01.1_rev03_Orbit-'+strOrbits[iorbit]+'_'+UT_DATE_ORBIT_START+'.png'
	RAA_Plot_filename_thumb = outpath+'RAA_Level2b_Orbit-'+strOrbits[iorbit]+'_'+UT_DATE_ORBIT_START+'.thumb.png'

	
	; determine whether each CIPS pixel is over land or water using MERRA2 surface GPHT for plotting continents
	continent_plot = RAYLEIGH_ALBEDO_ANOMALY*0.	
	for ii = 0, n_elements(RAYLEIGH_ALBEDO_ANOMALY) - 1L do begin
		xlat = latitude[ii]
		xlon = longitude[ii]
		x =min(abs(lat_mask - xlat),latindex)
		x =min(abs(lon_mask - xlon),lonindex)
		continent_plot[ii] = continent_mask[lonindex[0],latindex[0]]
	endfor


	erase
; add white box around title text



	;Define linear colortable	
	levels1 = (findgen(22)*.3-3.)
	nlvls  = n_elements(levels1)
	col1 = (1 + indgen(nlvls)) * 250. / nlvls	; define colors
  
  ;create generic plotting arrays
	xarray = findgen(n_elements(RAYLEIGH_ALBEDO_ANOMALY[*,0]))
	yarray = findgen(n_elements(RAYLEIGH_ALBEDO_ANOMALY[0,*]))

  contour,continent_plot,xarray, yarray,/fill,level = [0,.85], c_color = [0,150],$
  position = [.01,.01,.9,.9], xstyle=5,ystyle=5, color = 0, /nodata
  x = where(finite(continent_plot) eq 0,comp = y)
  continent_plot1 =continent_plot
  continent_plot1[x] = 1.
  continent_plot1[y] = 0.
  contour,continent_plot1,xarray, yarray,/fill,level = [0], c_color = [0],/overplot

	;Determine if orbit is during PMC season and add red shading poleward of 60
	year = strmid(UT_DATE_ORBIT_START,0,4)
	month = strmid(UT_DATE_ORBIT_START,4,2)
	day = strmid(UT_DATE_ORBIT_START,6,2)
	current_doy = julday(month,day,year)-julday('01','01',year)+1
	if current_doy ge 321 or current_doy le 54 then begin
		SH_PMC = 1
	endif else SH_PMC = 0
	if current_doy ge 135 and current_doy le 241 then begin
		NH_PMC = 1
	endif else NH_PMC = 0

	loadct, 70
	if NH_PMC eq 1 then begin
		x = where(latitude ge 60.,comp = y)
		continent_plot2 = continent_plot*0. + 1.
		continent_plot2[y] = 0
		contour,continent_plot2,xarray, yarray,/fill,level = [1], c_color = [80],/overplot
	endif
	if SH_PMC eq 1 then begin
		x = where(latitude le -60.,comp = y)
		continent_plot2 = continent_plot*0. + 1.
		continent_plot2[y] = 0
		contour,continent_plot2,xarray, yarray,/fill,level = [1], c_color = [80],/overplot
	endif
	
	
  ;add continents	
	loadct,0
	continent_plot = smooth(continent_plot,5,/edge_truncate)
  contour,continent_plot,xarray, yarray,/fill,level = [.85], c_color = [150],/overplot

	;Add lat/lon contours	
	x = where(longitude le 9.)
	longitude[x] = !values.f_nan
	x=where(longitude gt 359.0)
	longitude[x] = 360.
	
	x = where(abs(latitude) gt 70.)
	xlons = longitude
	xlons[x] = !values.f_nan

  contour, latitude, xarray,yarray,levels = [-80,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,80],c_linestyle = 0, /follow, /overplot,color = 250,$
     c_charsize = 2, c_thick = 2,C_Labels=Replicate(1, 13)
  contour, latitude, xarray,yarray,levels = [-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80,90],c_linestyle = 0, /overplot,color = 250,$
       c_charsize = 2, c_thick = 2,C_Labels=Replicate(0, 19)
  contour, latitude, xarray,yarray,levels = [0],c_linestyle = 0, /follow, /overplot,color = 250, c_charsize = 2, c_thick = 5,C_Labels=[0],C_annot = ['Eq']
  contour, xlons, xarray,yarray,levels = findgen(37)*10.,c_linestyle = 0, /follow, /overplot,color = 250, c_charsize = 2, c_thick =2
  contour, longitude, xarray,yarray,levels = findgen(37)*10.,c_linestyle = 0, /overplot,color = 250, c_charsize = 0, c_thick =2,C_Labels=Replicate(0, 36)
  contour,continent_plot,xarray, yarray,/fill,level = [0,.85], c_color = [0,150],$
    position = [.01,.01,.9,.9], xstyle=5,ystyle=5, color = 0, /nodata, /noerase

	
	
	;loop over CIPS pixels and plot using psym
	loadct,60
	for i=0, n_elements(rayleigh_albedo_anomaly[*,0]) - 1L do begin
	for j=0, n_elements(rayleigh_albedo_anomaly[0,*]) - 1L do begin
		dum = min(abs(levels1 - rayleigh_albedo_anomaly[i,j]),icol)
		if finite(dum) then oplot, [xarray[i],xarray[i]],[yarray[j],yarray[j]],color = col1[icol], psym=2,symsize = .05, thick = 5
	endfor
	endfor
loadct,0
	plot,[0,0],[0,1],xrange=[0,1],yrange=[0,1],/noeras,/nodata,$
	xstyle=4,ystyle = 4,position = [.01,.9,.8,.99]
		xbox=[0,1,1,0,0]
		ybox=[0,0,1,1,0]
	polyfill,xbox,ybox,color=250
	!type=2^2+2^3+2^6     ; no y title or ticsks
	plot, [0,0],[0,0],xstyle = 4, ystyle =  4,$
	  position = [.01,.9,.9,.99], charsize = 1.125,thick = 4,/noerase ; style= 4 supresses axis
	xyouts, .5,.5,plot_title, align = .5, charsize = 1.5,color=0

	; -----------------plot the color bar-----------------------	
	level = levels1
	slab=' '+strarr(n_elements(level))

	loadct,0
	!p.title = ' '
	position = [.9,0.01, .92,.9]
	plot,[0,0],[0,1],xrange=[0,1],yrange=[0,1],/noeras,/nodata,$
	xstyle=4,ystyle = 4,position = [.9,0.01, .94,.9]
		xbox=[0,1,1,0,0]
		ybox=[0,0,1,1,0]
		polyfill,xbox,ybox,color=250
		!type=2^2+2^3+2^6     ; no y title or ticsks
		plot, [0,0],[0,0],xstyle = 4, ystyle =  4,$
		position = [.01,.01,.99,.96], charsize = 1.125,thick = 4,/noerase ; style= 4 supresses axis
	xyouts, .94,.5,'Rayleigh Albedo Anomaly (%)', orientation = 90, align = .5, charsize = 1.5

	plot,[0,0],[0,1],xrange=[0,10],yrange=[0,1],/noeras,xticks=n_elements(levels1)-1L,$
		position = position,xstyle=1,xtickname=slab, color = 250

	loadct,60

		xbox=[0,10,10,0,0]
 	x2=0
	for j=1,n_elements(col1)-1 do begin
		dx= 1./(n_elements(level)-1.)
		ybox=[x2,x2,x2+dx,x2+dx,x2]
		polyfill,xbox,ybox,color=col1[j-1]
		x2=x2+dx
	endfor
	loadct,0
	slab=strcompress(string(format='(f8.3)',level),/remove_all)
	slabcolor = fltarr(n_elements(level))*0.
	;slabcolor[0:2] = 255	
;	slabcolor[-4:-1] = 255
	x1=dx/2
	
	for i=0L,n_elements(slab)-2L do begin
		slab0=slab[i]
		flab0=float(slab[i])
		slab0=strcompress(string(format='(f5.1)',flab0),/remove_all)
		xyouts,5,x1-dx/2.,slab0,charsize=1.3,/data,color=slabcolor[i], align = .5
		x1=x1+dx
	endfor

 ; convert .ps to png and save image and thumbnail
	device, /close
	spawn,'convert test2_'+strorbits[iorbit]+'.ps '+RAA_Plot_filename
	spawn,'convert test2_'+strorbits[iorbit]+'.ps -resize 8.5% '+RAA_Plot_filename_thumb
endfor
end