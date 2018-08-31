;---------------------------------------------------------------------------------------------------
; Reads in WACCM data and plots zonal and meridional phase of the Arctic vortex (tilt with altitude) 
; -30 to +30 days around composite ES event
;
;	 -------------------------------
;       |         Lynn Harvey           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 11/24/2012      |
;	 -------------------------------
;
@stddat			; Determines the number of days since Jan 1, 1956
@kgmt			; This function computes the Julian day number (GMT) from the
@ckday			; This routine changes the Julian day from 365(6 if leap yr)
@kdate			; gives back kmn,kdy information from the Julian day #.

;-----------------------------------------------------

SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.55,0.1,0.55]
yorig=[0.55,0.55,0.1,0.1]
xlen=0.35
ylen=0.35
cbaryoff=0.02
cbarydel=0.01

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
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
;
; save postscript version
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/times,filename='../Figures/waccm_vortex_phase_ES_composite.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
restore,'/Users/harvey/Desktop/Harvey_etal_2014/Post_process/elevated_strat.sav'
restore,'/Users/harvey/Desktop/Harvey_etal_2014/Post_process/WACCM_ES_daily_max_T_Z.sav'
ndates = 30.*n_elements(dayzeros)
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
ESday = fltarr(n_elements(dayzeros)*30L)
Date = fltarr(n_elements(dayzeros)*30L)
ndays = 61L
dayofES = 1L
niday = 0L
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
restore,'ES_composite_days_waccm_dm30-dp30.sav' ;,alon,alat,th,nc,nr,nth,pv2_comp,p2_comp,u2_comp,v2_comp,qdf2_comp,q2_comp,gph2_comp,ttgw2_comp,sf2_comp,mark2_comp
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
; timeseries of max stratopause height
;
sdays=string(format='(I3.2)',-30+indgen(ndays))
erase
!type=2^2+2^3
;set_viewport,.2,.9,.85,.95
;plot,findgen(ndays),maxheighttheta[0,*],color=0,xtitle='Days From ES Onset',thick=6,yrange=[500.,10000.],/noeras,ytitle='Theta (K)',$
;     xticks=ndays/10,xtickname=sdays(0:ndays-1:10),title='WACCM ES Composite',/nodata
;for i=0L,14L do oplot,findgen(ndays),maxheighttheta[i,*],color=0,psym=0
;plots,30,500
;plots,30,10000.,/continue,color=0,thick=3
;;plot,findgen(ndays),maxheighttheta[iES,*],color=0,thick=6
;for ies=0L,14L do begin 
;for iday=0,ndays-1L do begin
;    oplot,[iday,iday],[maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=(iday/(ndays+1.))*mcolor,symsize=1.5
;    a=findgen(9)*(2*!pi/8.)
;    usersym,1.5*cos(a),1.5*sin(a)
;    oplot,[iday,iday],[maxheighttheta[iES,iday],maxheighttheta[iES,iday]],psym=8,color=0,symsize=1.5
;    a=findgen(8)*(2*!pi/8.)
;    usersym,1.5*cos(a),1.5*sin(a),/fill
;endfor
;endfor

sdays=strarr(61)
    for iday=0,60L do begin
        iday0=iday-30
        if iday0 lt 0L then sday=string(format='(i3.2)',iday0)
        if iday0 ge 0L then sday=string(format='(i2.2)',iday0)
        sdays(iday)=sday
;
; read daily file
;
        pv2=reform(pv2_comp(*,*,*,iday))
        p2=reform(p2_comp(*,*,*,iday))
        u2=reform(u2_comp(*,*,*,iday))
        v2=reform(v2_comp(*,*,*,iday))
        qdf2=reform(qdf2_comp(*,*,*,iday))
        q2=reform(q2_comp(*,*,*,iday))
        gph2=reform(gph2_comp(*,*,*,iday))
        ttgw2=reform(ttgw2_comp(*,*,*,iday))
        sf2=reform(sf2_comp(*,*,*,iday))
        mark2=reform(mark2_comp(*,*,*,iday))
;
        for ies=0L,14L do begin 
            if maxheighttheta[ies,iday] le 0. or maxheighttheta[ies,iday] gt 10000. then maxheighttheta[ies,iday]= maxheighttheta[ies,iday-1]
            xtheta = th - maxheighttheta[ies,iday]
        endfor
	x = where(xtheta lt 0L and th ge 500. and th le 10000., nx)
        th2=th	;(x)
        mark2=mark2	;(*,*,x)
        sf2=sf2	;(*,*,x)
        xphase=0*th2
        yphase=0*th2
;
; max theta below stratopause max
;
;       thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(10000.) - min(alog(th2),/nan))  *  254.
;       for ies=0L,14L do begin
;       maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(10000.)))*254.
;       if maxheighttheta[ies,iday] gt 10000. then begin
;          thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(maxheighttheta[ies,iday]) - min(alog(th2),/nan))*254.
;          maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(maxheighttheta[ies,iday])))*254.
;       endif
;       endfor
;
; plot
;
        if iday le 20 then begin
           set_viewport,.225,.475,.65,.9
           map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,color=0
	   for ii = 0L, nx - 1L do begin
               waccmdailymark = transpose(reform(mark2[nr2:nr-1,*,nx-1L-ii]))
               marker=fltarr(nc+1L,nr2)
               marker[0L:nc-1L,*] = waccmdailymark
               marker[nc,*] = marker(0,*)
               waccmdailysf = transpose(reform(sf2[nr2:nr-1,*,nx-1L-ii]))
               sf=fltarr(nc+1L,nr2)
               sf[0L:nc-1L,*] = waccmdailysf
               sf[nc,*] = sf(0,*)
;              contour,marker,lons,alat(nr2:nr-1),levels=[.1],/overplot,color=thlevel[ii],thick=8
;
; superimpose center of vortex
;
               index=where(marker gt 0.)
               if index(0) ne -1L then begin
                  index2=where(sf(index) eq min(sf(index)))
                  oplot,x2d(index(index2)),y2d(index(index2)),color=(iday/(ndays+1.))*mcolor,psym=8
                  a=findgen(9)*(2*!pi/8.)
                  usersym,1.5*cos(a),1.5*sin(a)
;                 oplot,x2d(index(index2)),y2d(index(index2)),color=mcolor,psym=8
                  a=findgen(9)*(2*!pi/8.)
                  usersym,1.5*cos(a),1.5*sin(a),/fill
                  xphase(nx-1-ii)=x2d(index(index2(0)))
                  yphase(nx-1-ii)=y2d(index(index2(0)))
               endif
            endfor
        endif

        if iday ge 40 then begin
           set_viewport,.625,.875,.65,.9
           map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,color=0
           for ii = 0L, nx - 1L do begin
               waccmdailymark = transpose(reform(mark2[nr2:nr-1,*,nx-1L-ii]))
               marker=fltarr(nc+1L,nr2)
               marker[0L:nc-1L,*] = waccmdailymark
               marker[nc,*] = marker(0,*)
               waccmdailysf = transpose(reform(sf2[nr2:nr-1,*,nx-1L-ii]))
               sf=fltarr(nc+1L,nr2)
               sf[0L:nc-1L,*] = waccmdailysf
               sf[nc,*] = sf(0,*)
;              contour,marker,lons,alat(nr2:nr-1),levels=[.1],/overplot,color=thlevel[ii],thick=8
;
; superimpose center of vortex
;
               index=where(marker gt 0.)
               if index(0) ne -1L then begin
                  index2=where(sf(index) eq min(sf(index)))
                  oplot,x2d(index(index2)),y2d(index(index2)),color=(iday/(ndays+1.))*mcolor,psym=8
                  a=findgen(9)*(2*!pi/8.)
                  usersym,1.5*cos(a),1.5*sin(a)
;                 oplot,x2d(index(index2)),y2d(index(index2)),color=mcolor,psym=8
                  a=findgen(9)*(2*!pi/8.)
                  usersym,1.5*cos(a),1.5*sin(a),/fill
                  xphase(nx-1-ii)=x2d(index(index2(0)))
                  yphase(nx-1-ii)=y2d(index(index2(0)))
               endif
           endfor
        endif
;
; phase plots
;
        theta=th2
        index=where(xphase ne 0.)
        theta=theta(index)
        xphase=xphase(index)
        yphase=yphase(index)
        index=where(theta le 9000.)
        theta=theta(index)
        xphase=xphase(index)
        yphase=yphase(index)
        irot=135.
        index=where(xphase ge irot)
        if index(0) ne -1L then xphase(index)=xphase(index)-360.
;       if dailymaxheightlon[iES,iday] gt 180. then dailymaxheightlon[iES,iday]=dailymaxheightlon[iES,iday]-360.
;
; longitude phase
;
        if iday le 20L then begin
           set_viewport,.2,.5,.475,.625
           plot,xphase,theta,color=0,xtitle='Longitude',ytitle='Theta (K)',xrange=[irot-360.,irot],thick=6,yrange=[0.,10000.],xticks=3,/noeras,/nodata,title='-30 to -10'
           oplot,xphase,theta,color=(iday/(ndays+1.))*mcolor,psym=8,symsize=0.5
        endif
        if iday ge 40L then begin
           set_viewport,.6,.9,.475,.625
           plot,xphase,theta,color=0,xtitle='Longitude',xrange=[irot-360.,irot],thick=6,yrange=[0.,10000.],xticks=3,/noeras,/nodata,title='+10 to +30'
           oplot,xphase,theta,color=(iday/(ndays+1.))*mcolor,psym=8,symsize=0.5
        endif
;
; latitude phase
;
        if iday le 20L then begin
           set_viewport,.2,.5,.25,.4
           plot,yphase,theta,color=0,xtitle='Latitude',ytitle='Theta (K)',xrange=[40.,90.],yrange=[0.,10000.],thick=6,/noeras,/nodata
           oplot,yphase,theta,color=(iday/(ndays+1.))*mcolor,thick=8
        endif
        if iday ge 40L then begin
           set_viewport,.6,.9,.25,.4
           plot,yphase,theta,color=0,xtitle='Latitude',xrange=[40.,90.],yrange=[0.,10000.],thick=6,/noeras,/nodata
           oplot,yphase,theta,color=(iday/(ndays+1.))*mcolor,thick=8
        endif

        jumpday:
endfor		; loop over days 
;
; color bar
;
        level1 = findgen(20)
        nlvls  = n_elements(sdays)
        col1 = (1 + indgen(nlvls)) * icolmax / nlvls    ; define colors
        !type=2^2+2^3+2^6                       ; no y title or ticsks
        imin=min(level1)
        imax=max(level1)
        slab=['-30','-20','-10','0','10','20','30']
        set_viewport,.2,.9,.15,.17
        plot,[imin,imax],[0,0],yrange=[0,10],xrange=[0,10],/noeras,xticks=n_elements(slab)-1L,$
              xstyle=1,xtickname=slab, xtitle = 'Days From ES Onset',color=0
        ybox=[0,10,10,0,0]
        x2=0
        for j=1,n_elements(col1)-1 do begin
            dx= 10./(n_elements(col1)-1L)
            xbox=[x2,x2,x2+dx,x2+dx,x2]
            polyfill,xbox,ybox,color=col1(j-1)
            x2=x2+dx
        endfor
;       slabcolor = fltarr(n_elements(slab))*0.
;       slabcolor[0:7] = 255
;       x1=min(level1)+dx/2 + dx
;       for i=0L,20-1L do begin
;          xyouts,x1-dx/2,.76,sdays(i*3),charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
;          x1=x1+dx*3
;       endfor

        if setplot ne 'ps' then stop
        if setplot eq 'ps' then begin
           device, /close
           spawn,'convert -trim ../Figures/waccm_vortex_phase_ES_composite.ps -rotate -90 /Users/harvey/Desktop/Harvey_etal_2014/Figures/waccm_vortex_phase_ES_composite.png'
           spawn,'rm -f ../Figures/waccm_vortex_phase_ES_composite.ps'
        endif
end
