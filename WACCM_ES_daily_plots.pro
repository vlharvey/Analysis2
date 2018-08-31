;-----------------------------------------------------------------------------------------------------------------------------
;FIGURE 5 FOR STRATOPAUSE CLIMATOLOGY
;		 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 11/18/2011      |
;		 -------------------------------
;
;
;
;-----------------------------------------------------------------------------------------------------------------------------
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
@/Volumes/MacD68-1/france/idl_files/rd_waccm_nc3
@/Volumes/MacD68-1/france/idl_files/frac_index


;--------------------------




px1a = .07
px1b = 0.3
px2a = .32
px2b = .55

py1a = .75
py1b = .96
py2a = .46
py2b = .66
py3a = .15
py3b = .35


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

loadct,39
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

days = 0

months = ['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
MONTH = ['04','05','06','07','08','09','10','11']
titles1 = ['a)-(60-30)','c)-(30-0)','e)0-30','g)30-60','i)','k) Sep','m) Oct','o) Nov']
titles2 = ['b)','d)','f)','h)','j)','l)','n)','p)','r)']



restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/Post_process_1_height_temps_mark.sav'

restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/Figure_4_pre_process_composit_ES_WACCM.sav'

;




restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'



dailymaxtemplat = fltarr(n_elements(dayzeros), 61L)
dailymaxtemplon = fltarr(n_elements(dayzeros), 61L)
dailymaxheightlat = fltarr(n_elements(dayzeros), 61L)
dailymaxheightlon = fltarr(n_elements(dayzeros), 61L)
ESdate = strarr(n_elements(dayzeros)*61L)
inday = 0L

	
	
	for nnday = 0, n_elements(dayzeros) - 1L do begin
    sdate = dates[dayzeros(nnday)]
    print,sdate
    
    
    
;
; determine date of day -30 
;
    imn=long(strmid(sdate,4,2))
    idy=long(strmid(sdate,6,2))
    iyr=long(strmid(sdate,0,4))
    jday0 = JULDAY(imn, idy, iyr)
    jdaym30=jday0-30L
    CALDAT, jdaym30, imnm30 , idym30 , iyrm30

    for ndays=0,60L do begin
        jday=jdaym30+ndays
        CALDAT, jday, imn , idy , iyr
        print,ndays-30L,imn,idy,iyr
	inday = inday + 1L

    sdateindex = dayzeros[nnday] + ndays - 30L
    date = dates[sdateindex] 

    waccmdailytemp = reform(waccmtemps[sdateindex,*,*])
    waccmdailyheight = reform(waccmheight[sdateindex,*,*])
    waccmdailymark = reform(waccmmark[sdateindex,*,*])
            
        ;
; plot and save plot with ES#nnday(1-15) and days=ndays-30 in the filename
;



    ;-------------create wrap around for plot------------
    marker = fltarr(n_elements(waccmdailymark[*,0L]) + 2L, n_elements(waccmdailymark[0L,*]))
    marker[0L:n_elements(waccmdailymark[*,0L]) - 1L, *] = waccmdailymark[*,*]
   marker[n_elements(waccmdailymark[*,0L])-1L:n_elements(waccmdailymark[*,0L])+1L,*] = waccmdailymark[0L:2L,*]
   lons = fltarr(n_elements(lon) + 2L)
   lons[0:n_elements(lon)-1L] = lon
   lons[n_elements(lon)-1L:n_elements(lon)+1L] = lon[0L:2L]


   WACCMt = fltarr(n_elements(waccmdailytemp[*,0L]) + 2L, n_elements(waccmdailytemp[0L,*]))
    WACCMt[0L:n_elements(waccmdailytemp[*,0L]) - 1L, *] = waccmdailytemp
   WACCMt[n_elements(waccmdailytemp[*,0L])-1L:n_elements(waccmdailytemp[*,0L])+1L,*] = waccmdailytemp[0L:2L,*]
     WACCMZ = fltarr(n_elements(waccmdailyheight[*,0L]) + 2L, n_elements(waccmdailyheight[0L,*]))
    WACCMZ[0L:n_elements(waccmdailyheight[*,0L]) - 1L, *] = waccmdailyheight
   WACCMZ[n_elements(waccmdailyheight[*,0L])-1L:n_elements(waccmdailyheight[*,0L])+1L,*] = waccmdailyheight[0L:2L,*]

   
	lats1 = fltarr(n_elements(lons),n_elements(lat)) 
	for ii = 0, n_elements(lat) - 1L do lats1[*,ii] = lat[ii]
	lons1 = fltarr(n_elements(lons),n_elements(lat)) 
	for ii = 0, n_elements(lons) - 1L do lons1[ii,*] = lons[ii]

   x = where(waccmz le 0.,nx)
   if nx gt 0L then waccmz[x] = !values.f_nan
   x = where(waccmt le 0.,nx)
   if nx gt 0L then waccmt[x] = !values.f_nan
;   WACCMz = smooth(WACCMz,3,/nan)
;   WACCMt = smooth(WACCMt,3,/nan)



  
x = where(marker lt 0. and marker gt -99.,nx)
if nx gt 0L then marker[x] = -1.

x = where(marker gt 0.,nx)
if nx gt 0L then marker[x] = 1.

x = where(marker le -99.,nx)
marker[x] = 0.

marker = smooth(marker,7,/nan)




    ;----------------------plot code---------------------      ;x = where(lons eq 0.00)





plot,[0,0],[0,0], position=[.001,.001,.5,.999], xstyle = 4, ystyle = 4,$
title = strtrim(strcompress(date),2) +' Day '+ strtrim(strcompress(string(format='(I3.2)',ndays - 30L)),2)+$
	', ES event ' + strtrim(strcompress(string(format='(I2.2)',nnday + 1L)),2)

loadct, 0
xyouts, .15,.65,'Temperature (K)'
xyouts, .7,.65,'Height (km)'

; ------------- TEMPERATURE PLOTS -----------------------

loadct, 39

    ;----------------------plot code---------------------      ;x = where(lons eq 0.00)
      ;lons[x] = 360.
      level1 = findgen(11)*3.+246.
      level1a = findgen(26)*3.+222.
            level3 = findgen(1)+.2
      nlvls  = n_elements(level1)
      col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors
x = where(lat lt 99.)
xx = where(lat le 99.)

         map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,$
          position = [px1a,py1a,px1b,py1b]


PlotArray = smooth(WACCMt,5,/edge,/nan)
x = where(WACCMt gt max(level1),nx)
if nx gt 0 then PlotArray[x] = max(level1)-1.
x = where(WACCMt lt min(level1),nx)
if nx gt 0 then PlotArray[x] = min(level1) +1.



    contour, PlotArray[*,xx], lons, lat[xx], levels=level1, /cell_fill, c_color = col1,/overplot,$
      min_value = -2., max_value =40000.,xticks = 4, /close, color = 0

    contour,WACCMt[*,xx], lons, lat[xx], levels=level1a,color=0,charsize=3,$
        min_value = -2., max_value =40000.,/overplot, xticks=5, thick=1

     map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin, $
     position = [px1a,py1a,px1b,py1b]

     oplot, findgen(361), 0.1+0.*findgen(361), color = 0

loadct, 0

	contour,marker[*,xx], lons, lat[xx], levels = [.1], /overplot, color =0.,thick = 8,$
	 position = [px1a,py1a,px1b,py1b]

	 	contour,marker[*,xx], lons, lat[xx], levels = [-.1], /overplot, color =275.,thick = 8,$
	 position = [px1a,py1a,px1b,py1b]




	maxtemp = 0.*lat                                                              
	for i = 0L, n_elements(lat) - 1L do maxtemp[i]= max(plotarray[*,i],/nan) 
	; FIND ALL LOCAL TEMPERATURE MAXIMA (2 ADJACENT POINTS ARE LESS THAN CENTER POINT)
	localtmax = 0. * lat
	for j = 1L, n_elements(lat) - 2 do if maxtemp[j-1L] lt maxtemp[j] and maxtemp[j+1L] lt maxtemp[j] then localtmax[j] = maxtemp[j]
	index = where(localtmax ne 0, nindex)
	if nindex gt 0L then  begin
		ir = where(lat eq max(lat[index]))
		ir = ir[0]
		t = reform(waccmt[*,ir:*])
		lonsReform = reform(lons1[*,ir:*])
		latReform = reform(lats1[*,ir:*])
		y = where(t eq max(t,/nan),ny)
		dailymaxtemplon[nnday, ndays] = lonsReform[y[0]]
		dailymaxtemplat[nnday, ndays] = latReform[y[0]]
		ESdate[inday] = dates[sdateindex]
	endif


;;;;;;;;
    ; -----------------plot the color bar-----------------------
loadct, 39

      ;print, max(meanGEOS5strats)
      ;plot the color bar
      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level1)
      imax=max(level1)
      slab=' '+strarr(n_elements(level1))

      !p.title = ' '
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xticks=n_elements(level1)-1L,$
      position = [px1a-.03,0.7, px1b-.01,.73],xstyle=1,xtickname=slab

      ybox=[0,10,10,0,0]

      x2=imin
      for j=1,n_elements(col1)-1 do begin
        dx= level1[j] - level1[j-1]
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

loadct,0
     slab=strcompress(string(format='(f7.3)',level1),/remove_all)
		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255

		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255
    x1=min(level1)+dx/2 + dx

for i=1L,n_elements(slab)-1L do begin
   slab0=slab(i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(I5.3)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
   x1=x1+dx
endfor

loadct, 39

;----------------------


;-----------------HEIGHT PLOTS---------------------------

    ;----------------------plot code---------------------      ;x = where(lons eq 0.00)
      ;lons[x] = 360.
      level1 = findgen(11)*4.+40.
      level1a = findgen(40)*4.+20.
;level1 = [40,42,44,46,46.5,47,47.5,48,48.5,49,49.5,50,50.5,51,52,53,54,55]

            level3 = findgen(1)+.2
                  nlvls  = n_elements(level1)
      col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors
x = where(lat lt 99.)
xx = where(lat le 99.)


         map_set,90.,-90.,0,/ortho,  /grid,/noeras,/noborder,/contin,$
          position = [px2a,py1a,px2b,py1b]
 
PlotArray = smooth(WACCMz,5,/edge,/nan)

x = where(WACCMz gt max(level1),nx)
if nx gt 0 then PlotArray[x] = max(level1)-.5
x = where(WACCMz lt min(level1),nx)
if nx gt 0 then PlotArray[x] = min(level1) +.5



    contour, PlotArray[*,xx], lons, lat[xx], levels=level1, /cell_fill, c_color = col1,/overplot,$
      min_value = -2., max_value =40000.,xticks = 4, /close, color = 0

    contour,WACCMz[*,xx], lons, lat[xx], levels=level1a,color=0,charsize=3,$
        min_value = -2., max_value =40000.,/overplot, xticks=5, thick=1

         map_set,90.,-90.,0,/ortho,  /grid,/noeras,/noborder,/contin, $
          position = [px2a,py1a,px2b,py1b]
         oplot, findgen(361), 0.1+0.*findgen(361), color = 0

loadct, 0

	contour,marker[*,xx], lons, lat[xx], levels = [.1], /overplot, color =0.,thick = 8,$
	 position = [px1a,py1a,px1b,py1b]

	 	contour,marker[*,xx], lons, lat[xx], levels = [-.1], /overplot, color =275.,thick = 8,$
	 position = [px1a,py1a,px1b,py1b]


	 

	index = where(lat gt 20.)                                                                          
	t = reform(plotarray[*,index[0]:*])
	lonsReform = reform(lons1[*,index[0]:*])
	latReform = reform(lats1[*,index[0]:*])
	y = where(t eq max(t,/nan),ny)
	dailymaxheightlon[nnday, ndays] = lonsReform[y[0]]
	dailymaxheightlat[nnday, ndays] = latReform[y[0]]

	 	 	 
    ; -----------------plot the color bar-----------------------

  	 loadct, 39
    ;print, max(meanGEOS5strats)
      ;plot the color bar
      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level1)
      imax=max(level1)
      slab=' '+strarr(n_elements(level1))

      !p.title = ' '
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xticks=n_elements(level1)-1L,$
      position = [px2a+.01,.7,px2b,.73],xstyle=1,xtickname=slab


      ybox=[0,10,10,0,0]

      x2=imin
      for j=1,n_elements(col1)-1 do begin
        dx= level1[j] - level1[j-1]
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

;;;;;;;;      
     slab=strcompress(string(format='(f7.3)',level1),/remove_all)
		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255

		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:3] = 255

    x1=min(level1)+dx/2 + dx

for i=1L,n_elements(slab)-1L do begin
   slab0=slab(i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(I5.2)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
   x1=x1+dx
endfor


;-----------------------------------------------------------------------


   if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/ES_daily/ES_event_' + strtrim(strcompress(string(format='(I2.2)',nnday + 1L)),2) +$
	'_ES_Day_'+ strtrim(strcompress(string(format='(I3.2)',ndays - 30L)),2)+'_'+strtrim(strcompress(date),2)+'.png'
	spawn,'rm idl001.ppm idl.ps'
    endif

endfor

endfor

end
