;
;Purpose: create a 5 panel plot of daily zonally averaged temperatures for WACCM, MLS, and SABER for all days.
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
@/Volumes/MacD68-1/france/idl_files/stddat			; Determines the number of days since Jan 1, 1956
@/Volumes/MacD68-1/france/idl_files/kgmt			; This function computes the Julian day number (GMT) from the
								;    day, month, and year information.
@/Volumes/MacD68-1/france/idl_files/ckday			; This routine changes the Julian day from 365(6 if leap yr)
								;    to 1 and increases the year, if necessary.
@/Volumes/MacD68-1/france/idl_files/kdate			; gives back kmn,kdy information from the Julian day #.
@/Volumes/MacD68-1/france/idl_files/rd_ukmo_nc3			; reads the data from nc3 files
@/Volumes/MacD68-1/france/idl_files/date2uars			; This code returns the UARS day given (jday,year) information.
@/Volumes/MacD68-1/france/idl_files/plotPosition		; defines positions for n number of plots
@/Volumes/MacD68-1/france/idl_files/rd_GEOS5_nc3



SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.6,0.1,0.6,0.1,0.6]
yorig=[0.7,0.7,0.4,0.4,0.1,0.1]

xlen=0.25
ylen=0.25
cbaryoff=0.02
cbarydel=0.01

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

if setplot eq 'ps' then begin
  xsize=nxdim/100.
  ysize=nydim/100.
  set_plot,'ps'
  device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
         /bold,/color,bits_per_pixel=8,/helvetica,filename='idl.ps'
  !p.charsize=1.25    ; test with 1.5
  !p.thick=2
  !p.charthick=5
  !p.charthick=5
  !y.thick=2
  !x.thick=2
endif

loadct,13
mcolor=!p.color
icolmax=250
mcolor=icolmax
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
;
WACCMtemp = 0L 

spawn,'ls /Volumes/MacD68-1/france/WACCM_paper/Post_process/timeseries_SH_20*.sav',ifiles
nfile=n_elements(ifiles)

for iyear = 0L, nfile-1L do begin
restore,ifiles(iyear)
waccmtemps = tempseries

x = where(waccmtemps lt 0.,nx)
if nx gt 0L then waccmtemps[x] = !values.f_nan

days = dateseries

waccmstratopause = fltarr(n_elements(days))-99.
saberstratopause = fltarr(n_elements(days))-99.
mlsstratopause = fltarr(n_elements(days))-99.

waccmtemp = waccmtemps

waccmt = fltarr(n_elements(days), 400)

for ii = 0l, n_elements(days) - 1L do begin
;Result = POLY_FIT( X, Y, Degree [, CHISQ=variable] [, COVAR=variable] [, /DOUBLE] [, MEASURE_ERRORS=vector]
; [, SIGMA=variable] [, STATUS=variable] [, YBAND=variable] [, YERROR=variable] [, YFIT=variable] )

	st = reform(waccmtemp[ii,15:80])
	z1 = reform(altitude[15:80])
	x = where(st lt 0. or finite(st) eq 0.,nx, comp = y)
	if nx gt 0L then st[x] = !values.f_nan
	if max(y) ge 0 then begin
	x = poly_fit(z1[y],st[y],6)
	z = findgen(400)*.2 + 20
 	waccmt[ii,*] = x[0,6]*z^6 + x[0,5]*z^5 + x[0,4]*z^4 + x[0,3]*z^3+x[0,2]*z^2+x[0,1]*z+x[0,0]
	
	yz = where(finite(st))
	yzalt = max(z1[yz])
	top = where(z gt yzalt, nt)
	if nt gt 0L then waccmt[ii,top] = !values.f_nan
	x = where(z gt 25. and z le yzalt)	
	xx = where(waccmt[ii,x] eq max(waccmt[ii,x], /nan), nxx)
	if nxx ge 1L then begin
	if waccmt[ii,x[xx[0]]] gt 0. and waccmt[ii,x[xx[0]]]lt 1000. then waccmstratopause[ii] = z[x[xx[0]]]
	endif
	endif
endfor



x = where(waccmstratopause lt 21)
waccmstratopause[x] = 0

x = where(waccmt le 0., nx)
if nx gt 0L then waccmt[x] = !values.f_nan

waccmt = smooth(waccmt,3,/nan)
; ----------------PLOT CODE - ---------------------------------
n = 6 ; number of plots

erase
!type=2^2+2^3
loadct, 39



plot, [0,0],[0,0],xstyle = 4, ystyle =  4,$
position = [.01,.01,.99,.96], charsize = 1.125,thick = 4 ; style= 4 supresses axis

;Contour plot of zonal mean T
;level = 10.*findgen(11)+100.
level  = 6.*findgen(16)+210.		
levela  = 6.*findgen(21)+198.		
nlvls  = n_elements(level)

col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors

fake = waccmt
x = where(waccmt gt max(level),nx)
if nx gt 0L then fake[x] = max(level)-.1
x = where(waccmt lt min(level),nx)
if nx gt 0L then fake[x] = min(level) +.1




  contour, fake,days, z ,  levels=level, /cell_fill, c_color = col1,$
           xrange = [1,185],yrange=[20,90],ytitle='Altitude (km)', ystyle = 1,title = 'SH 20' + string(FORMAT='(I2.2)',iyear),$
           min_value = 0.,max_value = 500., xticks = 6, position = [.1,.71,.9,.93],$
		   xtickname = ['May1', 'Jun1', 'Jul1', 'Aug1', 'Sep1', 'Oct1','Nov1'], xtickv = [1, 32, 62, 93, 124, 154, 185]

  contour, waccmt,days, z ,  levels=levela, color = 0, /overplot

loadct, 0
  contour, waccmt,days, z ,  levels=[200], color =280, /follow,$
           /overplot,min_value = 0.,max_value = 500., thick = 10
  contour, waccmt,days, z ,  levels=[216], color =0, /follow,$
           /overplot,min_value = 0.,max_value = 500., thick = 10
loadct, 39

oplot, days,waccmstratopause,  psym = 8, color = 0,symsize = 1


loadct, 39


;-------------------------------------------------------------

    ; -----------------plot the color bar-----------------------

      ;print, max(meanGEOS5strats)
      ;plot the color bar
      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level)
      imax=max(level)
      slab=' '+strarr(n_elements(level))

      !p.title = ' '
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xticks=n_elements(level)-1L,$
          position = [.15,.62,.85,.65],xstyle=1,xtickname=slab,$
          xtitle='Temperature (K)', charsize =1.

      ybox=[0,10,10,0,0]

      x2=imin
      for j=1,n_elements(col1)-1 do begin
        dx= level[j] - level[j-1]
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

loadct,0
     slab=strcompress(string(format='(f7.3)',level),/remove_all)
		slabcolor = fltarr(n_elements(level))*0.
		slabcolor[0:4] = 255

      x1=min(level)+dx/2
for i=0L,n_elements(slab)-2L do begin
   slab0=slab(i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(I5.3)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=1.,/data,color=slabcolor[i],charthick=1
   x1=x1+dx
endfor
loadct, 39


  if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/WACCM_paper/Figures/timeseries/WACCM_maxtemp_timeseries_SH_20' + string(FORMAT='(I2.2)',iyear) + '.png'
	spawn,'rm idl001.ppm idl.ps'
  endif



endfor

end