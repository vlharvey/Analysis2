;-----------------------------------------------------------------------------------------------------------------------------
; Reads in WACCM data and returns stratopause height and temperature
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 05/24/2010      |
;	 -------------------------------
;
@stddat			; Determines the number of days since Jan 1, 1956
@kgmt			; This function computes the Julian day number (GMT) from the
@ckday			; This routine changes the Julian day from 365(6 if leap yr)
@kdate			; gives back kmn,kdy information from the Julian day #.
@date2uars		; This code returns the UARS day given (jday,year) information.

;-----------------------------------------------------

px1a = .22
px1b = .73
px2a = .52
px2b = .95
py1a = .50
py1b = .95
py2a = .45
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

days = 0
months = ['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
MONTH = ['04','05','06','07','08','09','10','11']
titles1 = ['a)-(60-30)','c)-(30-0)','e)0-30','g)30-60','i)','k) Sep','m) Oct','o) Nov']
titles2 = ['b)','d)','f)','h)','j)','l)','n)','p)','r)']

;The following determines the date range

restore, '/Users/harvey/Harvey_etal_2014/Post_process/elevated_strat.sav'
restore, '/Users/harvey/Harvey_etal_2014/Post_process/WACCM_ES_daily_max_T_Z.sav'
ndates = 30.*n_elements(dayzeros)
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
ESday = fltarr(n_elements(dayzeros)*30L)
Date = fltarr(n_elements(dayzeros)*30L)
nday = 0L
dayofES = 1L
niday = 0L
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)

for iES = 0L, n_elements(dayzeros) - 1L do begin
    sevent=string(format='(I2.2)',ies)
    ydate = where(dates eq dayzerodates[iES])
    icount=0
    for iday = 0L, 29L do begin
        sday=string(format='(I2.2)',iday)
	ifile = dates[ydate+iday]
	print, ifile
;
; save postscript version
;
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
;        !p.font=0
         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                /bold,/color,bits_per_pixel=8,/times,filename='../Figures/vortex_phase_ES_event_'+sevent+'_ES_Day_'+sday+'.ps'
         !p.charsize=1.25
         !p.thick=2
         !p.charthick=5
         !p.charthick=5
         !y.thick=2
         !x.thick=2
      endif

	nc = 0L
	ifiles=file_search(dir+ifile+'_3D_dyn.nc3',count=nfile)
	if ifiles[0] eq '' then continue
	result=strsplit(ifiles(0),'.',/extract)
    	result2=strsplit(result(4),'_',/extract)
    	sdate=result2(0)
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
;
; wrap around point
;
        lons = fltarr(nc+1L)
        lons[0:nc-1L] = alon
        lons[nc] = alon[0L]
;
; coordinate transformation
;
        nr2=nr/2
        x2d=fltarr(nc+1,nr/2)
        y2d=fltarr(nc+1,nr/2)
        for i=0,nc do y2d(i,*)=alat(nr/2:nr-1)
        for j=0,nr/2-1 do x2d(*,j)=lons
        xcn=fltarr(nc+1,nr2)
        ycn=fltarr(nc+1,nr2)
        for j=nr2,nr-1 do begin
            ANG = (90. - alat(j)) * RADG * 0.5
            FACTOR = TAN(ANG) * FAC20
            for i=0,nc do begin
                THETA = (lons(i) - 90.) * RADG
                xcn(i,j-nr2) = FACTOR * COS(THETA)
                ycn(i,j-nr2) = FACTOR * SIN(THETA)
            endfor
        endfor
;
; theta below stratopause max height
;
        if maxheighttheta[ies,iday] le 0. or maxheighttheta[ies,iday] gt 10000. then maxheighttheta[ies,iday]= maxheighttheta[ies,iday-1]
        xtheta = th - maxheighttheta[ies,iday]
        x = where(xtheta lt 0L and th ge 500. and th le 10000., nx)
        th2=th(x)
        marksf2=marksf2(*,*,x)
        sf2=sf2(*,*,x)
        xphase=0*th2
        yphase=0*th2
;
; max theta below stratopause max
;
        thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(10000.) - min(alog(th2),/nan))  *  254.
        maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(10000.)))*254.
        if maxheighttheta[ies,iday] gt 10000. then begin
           thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(maxheighttheta[ies,iday]) - min(alog(th2),/nan))*254.
           maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(maxheighttheta[ies,iday])))*254.
        endif
;
; plot
;
        erase
	map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,position = [px1a,py1a,px1b,py1b],color=0,$
                title='ES event ' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)+', Day '+ strtrim(strcompress(string(format='(I3.2)',iday)),2)+' ('+sdate+')'
          
	for ii = 0L, nx - 1L do begin
;       erase
;       map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,position = [px1a,py1a,px1b,py1b],color=0,$
;               title='ES event ' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)+', Day '+ strtrim(strcompress(string(format='(I3.2)',iday)),2)+' ('+sdate+')'

            waccmdailymark = transpose(reform(marksf2[nr2:nr-1,*,nx-1L-ii]))
            marker=fltarr(nc+1L,nr2)
            marker[0L:nc-1L,*] = waccmdailymark
            marker[nc,*] = marker(0,*)
            waccmdailysf = transpose(reform(sf2[nr2:nr-1,*,nx-1L-ii]))
            sf=fltarr(nc+1L,nr2)
            sf[0L:nc-1L,*] = waccmdailysf
            sf[nc,*] = sf(0,*)

            contour,marker,lons,alat(nr2:nr-1),levels=[.1],/overplot,color=thlevel[ii],thick=4
            smin=min(sf)
            smax=max(sf)
            slevel=smin+((smax-smin)/39.)*findgen(40)
;           contour,sf,lons,alat(nr2:nr-1),levels=slevel,/overplot,color=0
;
; superimpose center of vortex
;
            index=where(marker gt 0.)
            if index(0) ne -1L then begin
index2=where(sf(index) eq min(sf))
;              xavg=mean(xcn(index))
;              yavg=mean(ycn(index))
;              xindex=where(abs(xcn-xavg) eq min(abs(xcn-xavg)))
;              yindex=where(abs(ycn-yavg) eq min(abs(ycn-yavg)))
               oplot,x2d(index(index2)),y2d(index(index2)),color=thlevel(ii),psym=8,symsize=2
xphase(nx-1-ii)=x2d(index(index2(0)))
yphase(nx-1-ii)=y2d(index(index2(0)))

print,x2d(index(index2)),y2d(index(index2)),th2(nx-1-ii)
            endif
	endfor
	oplot, [dailymaxheightlon[iES,iday],dailymaxheightlon[iES,iday]],$
               [dailymaxheightlat[iES,iday],dailymaxheightlat[iES,iday]], psym = 8, color = maxthetalevel, symsize = 2
        a=findgen(9)*(2*!pi/8.)
        usersym,1.5*cos(a),1.5*sin(a)
        oplot, [dailymaxheightlon[iES,iday],dailymaxheightlon[iES,iday]],$
               [dailymaxheightlat[iES,iday],dailymaxheightlat[iES,iday]], psym = 8, color = 0, symsize = 2
        a=findgen(8)*(2*!pi/8.)
        usersym,1.5*cos(a),1.5*sin(a),/fill
;
; color bar
;
        level1 = findgen(nx)
        nlvls  = n_elements(level1)
        col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors
        !type=2^2+2^3+2^6			; no y title or ticsks
        imin=min(level1)
        imax=max(level1)
        slab=' '+strarr(n_elements(level1))
        plot,[imin,imax],[0,0],yrange=[0,10],xrange=[0,10],/noeras,xticks=n_elements(level1)-1L,$
        position = [px1a,.47,px1b,.5],xstyle=1,xtickname=slab, xtitle = 'Theta (K)',color=0
        ybox=[0,10,10,0,0]
        x2=0
        for j=1,n_elements(col1)-1 do begin
            dx= 10./(n_elements(col1)-1L)
            xbox=[x2,x2,x2+dx,x2+dx,x2]
            polyfill,xbox,ybox,color=col1(j-1)
            x2=x2+dx
        endfor
        slab=strcompress(long(reverse(th2)),/remove_all)
        slabcolor = fltarr(n_elements(slab))*0.
        slabcolor[0:3] = 255
        x1=min(level1)+dx/2 + dx
        for i=1L,n_elements(slab)-1L do begin
           slab0=slab((n_elements(slab)-1L)-i)
           flab0=float(slab(i))
           slab0=strcompress(string(format='(I6.2)',flab0),/remove_all)
           xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
           x1=x1+dx
        endfor

theta=th2
index=where(xphase ne 0.)
theta=theta(index)
xphase=xphase(index)
yphase=yphase(index)

index=where(xphase ge 90.)
if index(0) ne -1L then xphase(index)=xphase(index)-360.
!type=2^2+2^3
set_viewport,.1,.45,.1,.4
if dailymaxheightlon[iES,iday] gt 90. then dailymaxheightlon[iES,iday]=dailymaxheightlon[iES,iday]-360.
if icount eq 0L then begin
   plot,xphase,theta,color=0,xtitle='Longitude',xrange=[-270.,90.],thick=5,yrange=[500.,10000.],xticks=4,/noeras
   oplot,[dailymaxheightlon[iES,iday],dailymaxheightlon[iES,iday]],$
         [maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=0,symsize=1.5
endif
if icount gt 0L then begin
   plot,xphase,theta,color=0,xtitle='Longitude',xrange=[-270.,90.],thick=5,yrange=[500.,10000.],xticks=4,/noeras
   oplot,xphase,theta,color=(iday/31.)*mcolor,thick=5
   oplot,[dailymaxheightlon[iES,iday],dailymaxheightlon[iES,iday]],$
         [maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=(iday/31.)*mcolor,symsize=1.5
endif
set_viewport,.55,.9,.1,.4
if icount eq 0L then begin
   plot,yphase,theta,color=0,xtitle='Latitude',xrange=[40.,90.],yrange=[500.,10000.],thick=5
   oplot,[dailymaxheightlat[iES,iday],dailymaxheightlat[iES,iday]],$
         [maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=0,symsize=1.5
endif
if icount gt 0L then begin
   plot,yphase,theta,color=0,xtitle='Latitude',xrange=[40.,90.],yrange=[500.,10000.],thick=5,/noeras
   oplot,yphase,theta,color=(iday/31.)*mcolor,thick=5
   oplot,[dailymaxheightlat[iES,iday],dailymaxheightlat[iES,iday]],$
         [maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=(iday/31.)*mcolor,symsize=1.5
endif

        if setplot ne 'ps' then stop
        if setplot eq 'ps' then begin
           device, /close
           spawn,'convert -trim ../Figures/vortex_phase_ES_event_'+sevent+'_ES_Day_'+sday+'.ps -rotate -90 ../Figures/vortex_phase_ES_event_'+sevent+'_ES_Day_'+sday+'.png'
           spawn,'rm -f ../Figures/vortex_phase_ES_event_'+sevent+'_ES_Day_'+sday+'.ps'
        endif
    
        icount=icount+1L
        niday = niday + 1L

endfor		; loop over days 0 to +30
endfor		; loop over ES events
end
