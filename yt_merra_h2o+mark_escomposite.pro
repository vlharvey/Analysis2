;
; ES composite with MERRA
; contour H2O and the vortex edge as a function of latitude and time
; +/- 30 days around all ES events
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

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
yorig=[.25,.1,.55,.1]
xlen=0.8
ylen=0.6
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
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
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
    restore,idir+'yt_merra_co+mark_'+yearlab+'.sav'
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
ytco=ytco_all
yth2o=yth2o_all
ytmark=ytmark_all
ytspeed=ytspeed_all
;
; loop over theta
;
      for kk=0,nth2-1L do begin
          ilev=zindex(kk)
          rlev=th(ilev)
          slev=strcompress(long(rlev),/remove_all)+'K'
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
         device,/landscape,bits=8,filename='yt_merra_h2o+mark_'+slev+'_composite.ps'
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
plotarray=reform(yth2o(*,*,kk))
stitle=' '
;stitle='Composite CO, Arctic Vortex, and Wind Speed'
cbartitle=slev+' MLS H2O (ppmv)'
;if rlev le 2000. then begin
;   plotarray=reform(yth2o(*,*,kk))
;   stitle='Composite H2O, Arctic Vortex, and Wind Speed'
;   cbartitle=slev+' MLS H2O (ppmv)'
;endif
index=where(plotarray gt 0.)
imin=0.
imax=max(plotarray(index))
print,rlev,min(plotarray(index)),max(plotarray(index))
if rlev eq      4000.00 then begin
imin=0.5
imax=5.
endif ;6.16748      1.65434
if rlev eq      3000.00      then begin
imin=0.5
imax=8.
endif ;6.98830      3.40038
if rlev eq      2000.00      then begin
imin=5.
imax=8.
endif ;7.15542      5.18743
if rlev eq      1000.00      then begin
imin=4.
imax=7.
endif ; 6.86614      4.22600
if rlev eq      500.000      then begin
imin=1.
imax=6.
endif ;5.36497      1.66931

iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
index=where(plotarray eq -9999.)
if index(0) ne -1L then plotarray(index)=0./0.
contour,smooth(plotarray,3,/edge_truncate),-30.+findgen(kday),alat,levels=level,color=0,c_color=col1,/noeras,charsize=2,charthick=2,$
        xrange=[-30,30],yrange=[30,80],/fill,title=stitle,min_value=-9999.,ytitle='Latitude',xtitle='Days Since ES'
contour,smooth(plotarray,3,/edge_truncate),-30.+findgen(kday),alat,levels=level,color=0,/noeras,/follow,/overplot,min_value=-9999.
contour,smooth(ytmark(*,*,kk),5,/edge_truncate),-30.+findgen(kday),alat,/overplot,levels=[0.5,0.7,0.9],color=0,thick=15,/noeras,/follow,min_value=-9999.
loadct,0
contour,smooth(ytspeed(*,*,kk),3,/edge_truncate),-30.+findgen(kday),alat,/overplot,levels=30.+10*findgen(20),color=255,thick=10,/noeras,/follow,min_value=-9999.
;plots,0,30
;plots,0,80,/continue,thick=5,color=150
loadct,39
omin=min(level)
omax=max(level)
set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
!type=2^2+2^3+2^6
plot,[omin,omax],[0,0],yrange=[0,10],charsize=2,charthick=2,$
      xrange=[omin,omax],xtitle=cbartitle,/noeras,xstyle=1,color=0
ybox=[0,10,10,0,0]
x1=omin
dx=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim yt_merra_h2o+mark_'+slev+'_composite.ps -rotate -90 '+$
         'yt_merra_h2o+mark_'+slev+'_composite.jpg'
endif

endfor	; loop over altitude
end
