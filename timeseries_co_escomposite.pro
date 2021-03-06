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
xorig=[.15]
yorig=[.2]
xlen=0.7
ylen=0.7
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
idir='/Users/harvey/Harvey_etal_2014/Post_process/'
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
    restore,idir+'yt_geos5_co+mark_'+yearlab+'.sav'
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
yindex=where(alat ge 60.)
nr=n_elements(yindex)
ytco=reform(ytco_all(*,yindex,1))
yth2o=reform(yth2o_all(*,yindex,1))
ytmark=reform(ytmark_all(*,yindex,1))
ytspeed=reform(ytspeed_all(*,yindex,1))
day3d=0.*ytspeed
lat3d=0.*ytspeed
;for k=0,1L do for j=0,nr-1 do day3d(*,j,k)=-30+findgen(kday)
;for k=0,1L do for i=0,kday-1L do lat3d(i,*,k)=alat(yindex)
for j=0,nr-1 do day3d(*,j)=-30+findgen(kday)
for i=0,kday-1L do lat3d(i,*)=alat(yindex)

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
         device,/landscape,bits=8,filename='timeseries_co_escomposite.ps'
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
plotarray2=yth2o
stitle=' '
index=where(plotarray ne -9999.)
imin=min(plotarray(index))
imax=max(plotarray(index))
iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
index=where(plotarray eq -9999.)
if index(0) ne -1L then plotarray(index)=0./0.
 plot,findgen(kday),findgen(kday),xrange=[-30,30],yrange=[1,7],color=0,/nodata,xtitle='Days since ES onset',ytitle='CO (ppmv)',charsize=2,charthick=2,/noeras
for i=0L,kday-1L do begin
    arr=reform(plotarray(i,*))
    index=where(arr ne 0.)
    plots,i-30,min(arr(index))
    plots,i-30,max(arr(index)),/continue,thick=5,color=0
    oplot,[i-30,i-30],[min(arr(index)),min(arr(index))],psym=8,color=0
    oplot,[i-30,i-30],[max(arr(index)),max(arr(index))],psym=8,color=0
endfor
restore,'yt_waccm_co.sav'
 for i=-30,30 do oplot,[i,i],[min(plotarray(i+30,*)),min(plotarray(i+30,*))],psym=8,color=mcolor*.9
 for i=-30,30 do oplot,[i,i],[max(plotarray(i+30,*)),max(plotarray(i+30,*))],psym=8,color=mcolor*.9
for i=0L,kday-1L do begin
    arr=reform(plotarray(i,*))
    index=where(arr ne 0.)
    plots,i-30,min(arr(index))
    plots,i-30,max(arr(index)),/continue,thick=5,color=mcolor*.9
endfor
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim timeseries_co_escomposite.ps -rotate -90 timeseries_co_escomposite.jpg'
endif

end
