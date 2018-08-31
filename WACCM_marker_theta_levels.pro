;-----------------------------------------------------------------------------------------------------------------------------
; Reads in WACCM data and returns stratopause height and temperature
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 05/24/2010      |
;	 -------------------------------
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


;-----------------------------------------------------


px1a = .22
px1b = 0.73
px2a = .52
px2b = .95

py1a = .5
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
icolmax=255
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



;The following determines the date range

restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'
restore, '/Volumes/MacD68-1/france/ES_paper/Post_process/WACCM_ES_daily_max_T_Z.sav'

ndates = 30.*n_elements(dayzeros)
dir='/Volumes/MacD68-2/france/data/WACCM_4/Datfiles/mee00fpl_FW2.cam2.h3.dyns.'

ESday = fltarr(n_elements(dayzeros)*30L)
Date = fltarr(n_elements(dayzeros)*30L)
nday = 0L
dayofES = 1L
inday = 0L


	for ies = 0, n_elements(dayzeros) - 1L do begin
    sdate = dates[dayzeros(ies)]
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

    for iday=0,60L do begin
        jday=jdaym30+iday
        CALDAT, jday, imn , idy , iyr
        print,iday-30L,imn,idy,iyr
    sdateindex = dayzeros[ies] + iday - 30L
    date = dates[sdateindex] 

	ifile = date
		nc = 0L
		ifiles=file_search(dir+ifile+'_3D_dyn.nc3',count=nfile)
		if ifiles[0] eq '' then continue
		result=strsplit(ifiles(0),'.',/extract)
    		result2=strsplit(result(4),'_',/extract)
    		sdate=result2(0)
    		iflag=0
		;
		; read daily file
		;
    		ncfile0=ifiles(0)

       		ncid=ncdf_open(ncfile0)
    		result0=ncdf_inquire(ncid)
    		for idim=0,result0.ndims-1 do begin
        		ncdf_diminq,ncid,idim,name,dim
        		if name eq 'number_of_latitudes' then nr=dim
        		if name eq 'number_of_longitudes' then nc=dim
        		if name eq 'number_of_levels' then nth=dim
;       		print,'read ',name,' dimension ',dim
    		endfor
    		for ivar=0,result0.nvars-1 do begin
        		result=ncdf_varinq(ncid,ivar)
        		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        		if result.name eq 'latitude' then alat=data
        		if result.name eq 'longitude' then alon=data
        		if result.name eq 'theta' then th=data
        		if result.name eq 'IPV' then pv2=data
        		if result.name eq 'P' then p2=data
        		if result.name eq 'U' then u2=data
        		if result.name eq 'V' then v2=data
        		if result.name eq 'QDF' then qdf2=data
        		if result.name eq 'Q' then q2=data
        		if result.name eq 'GPH' then gph2=data
        		if result.name eq 'TTGW' then ttgw2=data
        		if result.name eq 'SF' then sf2=data
        		if result.name eq 'MARK' then marksf2=data
        		print,ivar,result.name,min(data),max(data)
		endfor
    		ncdf_close,ncid

  lons = fltarr(n_elements(alon) + 2L)
  lons[0:n_elements(alon)-1L] = alon
  lons[n_elements(alon)-1L:n_elements(alon)+1L] = alon[0L:2L]

lat = alat


    ;----------------------plot code---------------------      ;x = where(lons eq 0.00)




plot,[0,0],[0,0], position=[.001,.001,.999,.999], xstyle = 4, ystyle = 4,$
title = 'Day '+ strtrim(strcompress(string(format='(I3.2)',iday)),2)+', ES event ' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)

loadct, 39
xyouts, .4, .3, esdate[inday], color = 0L


    ;----------------------plot code---------------------      ;x = where(lons eq 0.00)


    
    
zth = where(abs(th - 5000.) eq min(abs(th - 5000.),/nan),nzth)
thlevel = reverse(alog(th[zth:*]) - min(alog(th[zth:*]),/nan))/(alog(5000.) - min(alog(th[zth:*]),/nan))  *  254.
maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(5000.)))*254.

if maxheighttheta[ies,iday] gt 5000. then begin
	zth = where(abs(th - maxheighttheta[ies,iday]) eq min(abs(th - maxheighttheta[ies,iday]),/nan),nzth)
	thlevel = reverse(alog(th[zth:*]) - min(alog(th[zth:*]),/nan))/(alog(maxheighttheta[ies,iday]) - min(alog(th[zth:*]),/nan))*254.
	maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(maxheighttheta[ies,iday])))*254.
endif		



		map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,$
		position = [px1a,py1a,px1b,py1b]


          
		xtheta = th - maxheighttheta[ies,iday]
		x = where(xtheta lt 0L, nx)
		for ii = 0L, nx - 1L do begin

			waccmdailymark = transpose(reform(marksf2[*,*,x[(nx-1L) - ii]]))

    		;-------------create wrap around for plot------------
    		marker = fltarr(n_elements(waccmdailymark[*,0L]) + 2L, n_elements(waccmdailymark[0L,*]))
    		marker[0L:n_elements(waccmdailymark[*,0L]) - 1L, *] = waccmdailymark[*,*]
   			marker[n_elements(waccmdailymark[*,0L])-1L:n_elements(waccmdailymark[*,0L])+1L,*] = waccmdailymark[0L:2L,*]


			contour,marker, lons, lat, levels = [.1], /overplot, color =thlevel[ii],thick = 4,$
	 		position = [px1a,py1a,px1b,py1b]


		endfor
		oplot, [dailymaxheightlon[iES,iday],dailymaxheightlon[iES,iday]],$
				[dailymaxheightlat[iES,iday],dailymaxheightlat[iES,iday]], psym = 8, color = maxthetalevel, symsize = 2


   ; -----------------plot the color bar-----------------------
	level1 = findgen(15)
	nlvls  = n_elements(level1)
	col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors
  	 loadct, 39
    ;print, max(meanGEOS5strats)
      ;plot the color bar
      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level1)
      imax=max(level1)
      slab=' '+strarr(n_elements(level1))

      !p.title = ' '
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[0,10],/noeras,xticks=n_elements(level1)-1L,$
      position = [px1a,.4,px1b,.45],xstyle=1,xtickname=slab, xtitle = 'Theta (K)'


      ybox=[0,10,10,0,0]

      x2=0
      for j=1,n_elements(col1)-1 do begin
        dx= 10./(n_elements(col1)-1L)
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

;;;;;;;;      
     slab=strcompress(string(format='(f10.0)',reverse(th[[(n_elements(thlevel)/nlvls)*level1 + zth[0]]])),/remove_all)
		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255

		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:3] = 255

    x1=min(level1)+dx/2 + dx

for i=1L,n_elements(slab)-1L do begin
   slab0=slab((n_elements(slab)-1L)-i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(I6.2)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
   x1=x1+dx
endfor





  if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/ES_daily/vortex_ES_event_' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)+ '_ES_Day_'+$
	strtrim(strcompress(string(format='(I3.2)',iday)),2)+'_'+esdate[inday]+'.png'
	spawn,'rm idl001.ppm idl.ps'
    endif
    
inday = inday + 1L

endfor
endfor
end


 