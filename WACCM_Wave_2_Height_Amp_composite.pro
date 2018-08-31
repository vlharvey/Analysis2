;
;Purpose: create a 5 panel plot of daily zonally averaged temperatures for HIRDLS, MLS, and WACCM for all days.
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
nydim=300
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
  !p.charsize=1.    ; test with 1.5
  !p.thick=2
  !p.charthick=5
  !p.charthick=5
  !y.thick=2
  !x.thick=2
  !p.font=0
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
hirdlstemp = 0L 




restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'
restore, '/Volumes/MacD68-1/france/ES_paper/Post_process/Platentary_Wave_2_amplitudes.sav'

idayzero = fltarr(15L)
x = dayzeros
for iES = 0L, n_elements(dayzeros) - 1L do begin
	idayzero[ies] = where(strmatch(esdates, dates[x[ies]]) eq 1L)
	print, idayzero
endfor

waccmheight = fltarr(61, n_elements(reform(waccmheight60[0,*])))
waccmstrat = fltarr(61)
x = where(waccmheight60 le 0.,nx)
if nx gt 0L then waccmheight60[x] = !values.f_nan
x = where(waccmstratopausepole le 0.,nx)
if nx gt 0L then waccmstratopausepole[x] = !values.f_nan

for iday = 0, 61L - 1L do begin
	xday = findgen(15) * 61L + iday
	waccmstrat[iday] = exp(mean(alog(waccmstratopausepole[xday]),/nan))
	for ilev = 0, n_elements(reform(waccmheight60[0,*])) - 1L do begin
		waccmheight[iday,ilev] = mean(waccmheight60[xday,ilev],/nan)
	endfor
endfor
x = where(WACCMheight le 0.,nx)
if nx gt 0L then WACCMheight[x] = !values.f_nan
WACCMheight = smooth(WACCMheight,3,/nan,/edge)
WACCMz60= WACCMheight




; ----------------PLOT CODE - ---------------------------------
n = 6 ; number of plots

erase
!type=2^2+2^3
loadct, 39



plot, [0,0],[0,0],xstyle = 4, ystyle =  4,$
position = [.01,.01,.99,.96], charsize = 1.125,thick = 4 ; style= 4 supresses axis

;Contour plot of zonal mean T
;level = 10.*findgen(11)+100.
level  = findgen(17)*125.
levela  = findgen(25)*200.

nlvls  = n_elements(level)


col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors


xyouts, .12,.95,'WACCM 60 N GPH Wave 2 Amplitude Composite'
xyouts, .5, .01,'Wave 1 Height Amplitude at 60!Uo!N N (m)', charsize =1., alignment = .5

;WACCM
WACCMz60 = WACCMz60 ; km to meters
fake = WACCMz60
x = where(WACCMz60 gt max(level),nx)
if nx gt 0L then fake[x] = max(level)-.1
x = where(WACCMz60 lt min(level),nx)
if nx gt 0L then fake[x] = min(level) +.1




  contour, /ylog, fake,findgen(61) - 30., pgrid,  levels=level, /cell_fill, c_color = col1,$
           xrange = [-30,30],yrange=[100.,.0001],ytitle='Pressure (hPa)', $
           min_value = 0.,max_value = 100000., xticklen = -.02, position = [.1,.4,.9,.85], ystyle = 8

	   
	  
	   
	   
  contour, /ylog, WACCMz60,findgen(61) - 30., pgrid,  levels=levela, color = 0, /overplot,/follow

  oplot, findgen(61) - 30.,WACCMstrat,  psym = 8, color = 255,symsize = .5
	oplot, [0,0],[1000,0], linestyle = 0, thick = 10
loadct, 0


loadct, 0
plot, [0,0],[0,0],  psym = 8, color = 255,symsize = .01, /noerase,$
	position = [.1,.4,.9,.85], yrange = [16,96],xrange = [1,75],ystyle = 4, xstyle = 4
 AXIS, YAXIS=1, YRANGE = [16,96.], YSTYLE = 1, $
  YTITLE = 'Approximate Altitude (km)', ytickv = [16,32,48,64,80,96],yticks = 5,ytickname = ['16','32','48','64','80','96']
loadct, 39

       		; -----------------plot the color bar (POLAR PLOTS)-----------------------

  			loadct, 39
      			!type=2^2+2^3+2^6			; no y title or ticsks
      			imin=min(level)
      			imax=max(level)
      			slab=' '+strarr(n_elements(level))

      			!p.title = ' '
      			plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xticks=n_elements(level)-1L,$
          		position = [.1,.13,.9,.25],xstyle=1,xtickname=slab

	    		ybox=[0,10,10,0,0]

    			x2=imin
  			for j=1,n_elements(col1)-1 do begin
     	   			dx= level[j] - level[j-1]
     	   			xbox=[x2,x2,x2+dx,x2+dx,x2]
     	   			polyfill,xbox,ybox,color=col1(j-1)
    				x2=x2+dx
    			endfor

			loadct,0
    			slab=strcompress(string(format='(i7.3)',level),/remove_all)
			slabcolor = fltarr(n_elements(level))*0.
			slabcolor[0:2] = 255
	
			slabcolor = fltarr(n_elements(level))*0.
			slabcolor[0:2] = 255
    			x1=min(level)+dx/2 + dx

			for i=1L,n_elements(slab)-1L do begin
   				slab0=slab(i)
   				flab0=float(slab(i))
      				slab0=strcompress(string(format='(i5.2)',flab0),/remove_all)
      				xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i], orientation= 90.
   				x1=x1+dx
			endfor




  if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/WACCM_PW_2_Amp_composite.png'
	spawn,'rm idl001.ppm idl.ps'
  endif



end