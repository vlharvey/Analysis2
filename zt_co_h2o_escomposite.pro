;
; ES composite
; scatter plot of CO vs H2O
; +/- 30 days around all ES events
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

sver='v3.3'

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=1000
nydim=700
xorig=[.1,0.1,0.55,0.55]
yorig=[.6,.2,.55,.1]
xlen=0.7
ylen=0.3
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
!noeras=1
dirm='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
mlsesdates=['20060130','20090205','20120130','20130123']
;sabesdates=['20060128','20090205','20120128','20130123']

nevent=n_elements(mlsesdates)
for ievent=0,nevent-1L do begin
    esdate0=mlsesdates(ievent)
    iyr=long(strmid(esdate0,0,4))
    imn=long(strmid(esdate0,4,2))
    idy=long(strmid(esdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
    jday0=jday-30
    jday1=jday+30
    CALDAT, jday0, lstmn ,lstdy , lstyr
    CALDAT, jday1, ledmn ,leddy , ledyr

    lstday=0L & ledday=0L
    if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
    if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
    restore,'/Users/harvey/Harvey_etal_2014/Post_process/yt_geos5_co+mark_'+yearlab+'.sav'
;
; average all events to create composite
;
    if ievent eq 0L then begin
       ytco_all=0.*ytco
       nytco_all=0.*ytco
       yth2o_all=0.*ytco
       nyth2o_all=0.*ytco
       ytmark_all=0.*ytco
       nytmark_all=0.*ytco
       ytspeed_all=0.*ytco
       nytspeed_all=0.*ytco
    endif

    index=where(ytco ne -9999.)
    if index(0) ne -1L then ytco_all(index)=ytco_all(index)+ytco(index)
    if index(0) ne -1L then nytco_all(index)=nytco_all(index)+1.0
    index=where(yth2o ne -9999.)
    if index(0) ne -1L then yth2o_all(index)=yth2o_all(index)+yth2o(index)
    if index(0) ne -1L then nyth2o_all(index)=nyth2o_all(index)+1.0
    index=where(ytmark ne -9999.)
    if index(0) ne -1L then ytmark_all(index)=ytmark_all(index)+ytmark(index)
    if index(0) ne -1L then nytmark_all(index)=nytmark_all(index)+1.0
    index=where(ytspeed ne -9999.)
    if index(0) ne -1L then ytspeed_all(index)=ytspeed_all(index)+ytspeed(index)
    if index(0) ne -1L then nytspeed_all(index)=nytspeed_all(index)+1.0
endfor
index=where(nytco_all gt 0.)
if index(0) ne -1L then ytco_all(index)=ytco_all(index)/float(nytco_all(index))
index=where(nyth2o_all gt 0.)
if index(0) ne -1L then yth2o_all(index)=yth2o_all(index)/float(nyth2o_all(index))
index=where(nytmark_all gt 0.)
if index(0) ne -1L then ytmark_all(index)=ytmark_all(index)/float(nytmark_all(index))
index=where(nytspeed_all gt 0.)
if index(0) ne -1L then ytspeed_all(index)=ytspeed_all(index)/float(nytspeed_all(index))
;
; rename for convenience
;
rlat=71.2500
print,alat
;read,'Enter Latitude ',rlat
index=where(alat eq rlat)
slat=strcompress(long(rlat),/remove_all)
ilat=index(0)
ytco=reform(ytco_all(*,ilat,*))
yth2o=reform(yth2o_all(*,ilat,*))
ytmark=reform(ytmark_all(*,ilat,*))
ytspeed=reform(ytspeed_all(*,ilat,*))
;
; save postscript version
;
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename='zt_geos5_co+mark_'+slat+'_composite.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.thick=2.0                   ;Plotted lines twice as thick
         !p.charsize=2.0
      endif

erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
plotarray=ytco
stitle='Composite MLS, Arctic Vortex, and Wind Speed'
index=where(plotarray ne -9999.)
imin=min(plotarray(index))
imax=max(plotarray(index))
iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
level=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,.2,.5,1,2,3,4,5]
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(plotarray eq -9999.)
if index(0) ne -1L then plotarray(index)=0./0.
contour,plotarray,-30.+findgen(kday),th(zindex),levels=level,color=0,c_color=col1,/noeras,charsize=2,charthick=2,$
        xrange=[-30,30],yrange=[500,4000],/fill,title=stitle,min_value=-9999.,ytitle='Theta (K)',xtitle='Days Since ES'
;contour,smooth(ytco,3,/edge_truncate),-30.+findgen(kday),th(zindex),levels=level,color=0,/noeras,/follow,/overplot,min_value=-9999.
loadct,0
;contour,smooth(ytmark,3,/edge_truncate),-30.+findgen(kday),th(zindex),/overplot,levels=[0.3,0.5,0.7],color=100,thick=10,/noeras,/follow,min_value=-9999.
contour,smooth(ytspeed,3,/edge_truncate),-30.+findgen(kday),th(zindex),/overplot,levels=30.+10*findgen(20),color=255,thick=10,/noeras,/follow,min_value=-9999.
plots,0,500
plots,0,4000,/continue,thick=5,color=150
loadct,39
omin=min(level)
omax=max(level)
xmnb=xorig(0)+xlen+0.05
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(0),yorig(0)+ylen
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(0)+ylen/5.0,slat+' N MLS CO (ppmv)',/normal,charsize=1.5,orientation=90.,color=0,charthick=1.5
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0,ytickname=string(format='(f3.1)',level(0:-1:2)),charsize=1.5,charthick=1.5
;
!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
plotarray=yth2o
stitle=' '
index=where(plotarray ne -9999.)
imin=min(plotarray(index))
imax=max(plotarray(index))
iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
index=where(plotarray eq -9999.)
if index(0) ne -1L then plotarray(index)=0./0.
contour,plotarray,-30.+findgen(kday),th(zindex),levels=level,color=0,c_color=col1,/noeras,charsize=2,charthick=2,$
        xrange=[-30,30],yrange=[500,4000],/fill,title=stitle,min_value=-9999.,ytitle='Theta (K)',xtitle='Days Since ES'
;contour,smooth(ytco,3,/edge_truncate),-30.+findgen(kday),th(zindex),levels=level,color=0,/noeras,/follow,/overplot,min_value=-9999.
loadct,0
contour,smooth(ytmark,3,/edge_truncate),-30.+findgen(kday),th(zindex),/overplot,levels=[0.3,0.5,0.7],color=0,thick=10,/noeras,/follow,min_value=-9999.
;contour,smooth(ytspeed,3,/edge_truncate),-30.+findgen(kday),th(zindex),/overplot,levels=30.+10*findgen(20),color=255,thick=10,/noeras,/follow,min_value=-9999.
plots,0,500
plots,0,4000,/continue,thick=5,color=150
loadct,39
omin=min(level)
omax=max(level)
xmnb=xorig(1)+xlen+0.05
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(1),yorig(1)+ylen
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(1)+ylen/5.0,slat+' N MLS H2O (ppmv)',/normal,charsize=1.5,orientation=90.,color=0,charthick=1.5
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0,ytickname=string(format='(f3.1)',level(0:-1:2)),charsize=1.5,charthick=1.5

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim zt_geos5_co+mark_'+slat+'_composite.ps -rotate -90 '+$
         'zt_geos5_co+mark_'+slat+'_composite.jpg'
endif

end
