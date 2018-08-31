;
; average all 8 years
; plot CIPS summary albedo, radius, MLS Temp, MLS H2O
; VLH 10/20/2011
;
@stddat
@kgmt
@ckday
@kdate
@mkltime
@smoothit
@fillit

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.6,0.15,0.6]
yorig=[0.55,0.55,0.1,0.1]
xlen=0.3
ylen=0.3
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
cdir='/Volumes/Data/CIPS_data/Datfiles/cips_3c_'
mdir='/Volumes/Data/CIPS_data/Datfiles_MLS_DMP/cips_3c_mls_'
;
; Ask interactive questions- get starting/ending date and p surface
;
shem='n'
lstyr=7
ialb=2
read,' Enter hemisphere (n,s) ',shem
read,' Enter PMC onset year (yy) ',lstyr
;read,' Enter Albedo threshhold (1,2,5,10) ',ialb
salb=string(format='(i2.2)',ialb)+'G'
if shem eq 'n' then begin
   shem='north'
   syr=string(format='(i2.2)',lstyr)
endif
if shem eq 's' then begin
   shem='south'
   syr=string(format='(i2.2)',lstyr)+string(format='(i2.2)',lstyr+1)
endif
restore,mdir+shem+'_'+syr+'_v04.20_r05_LT.sav'
;
; restore CIPS season
;
; ALB             FLOAT     = Array[589, 70]
; ALB_STD         FLOAT     = Array[589, 70]
; DFS             LONG      = Array[589]
; DOY             INT       = Array[589]
; IWC             FLOAT     = Array[589, 70]
; IWC_STD         FLOAT     = Array[589, 70]
; LATHI           INT       = Array[70]
; LATLO           INT       = Array[70]
; LON             FLOAT     = Array[589, 70]
; LTIME           FLOAT     = Array[589, 70]
; NBIN            INT       =       70
; NREV            LONG      =          589
; NUM_CLD         INT       = Array[589, 70]
; NUM_OBS         INT       = Array[589, 70]
; RAD             FLOAT     = Array[589, 70]
; RAD_STD         FLOAT     = Array[589, 70]
; REV             INT       = Array[589]
; SZA             FLOAT     = Array[589, 70]
; UT              FLOAT     = Array[589, 70]
; YEAR            INT       = Array[589]
;
restore,cdir+shem+'_'+syr+'_v04.20_r05_'+salb+'_cld.sav'
print,'restored '+cdir+shem+'_'+syr+'_v04.20_r05_'+salb+'_cld.sav'
;if shem eq 'south' then latlo=-1.*latlo
;if shem eq 'south' then lathi=-1.*lathi
;
; convert date to DOY
;
sdate=strcompress(date,/remove_all)
nn=n_elements(sdate)
doy=fltarr(nn)
for i=0L,nn-1L do begin
    iyr=long(strmid(sdate(i),0,4))
    imn=long(strmid(sdate(i),4,2))
    idy=long(strmid(sdate(i),6,2))
    z = kgmt(imn,idy,iyr,iday)
    doy(i)=1.0*iday
endfor
year=long(strmid(sdate,0,4))
;
; loop over DOYs in the cips summary file
;
n=findgen(nrev)
n1=1+findgen(nrev)
index=where(doy(n)-doy(n1) ne 0)
days=[doy(index),doy(nrev-1)]
years=[year(index),year(nrev-1)]
nday=n_elements(days)
;
; days array has to be monotonically increasing
;
index=where(days lt 100.)
if index(0) ne -1L then days(index)=days(index)+max(days)
;
; compute DFS
;
dfs=0.*doy
if shem eq 'north' then dfs=doy-172.			; June 21
if shem eq 'south' then begin
   index=where(doy gt 100.)
   dfs(index)=doy(index)-355.				; Dec 21
   index=where(doy lt 100.)
   dfs(index)=doy(index)+max(doy)-355. 
endif
;
; postscript file
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,/helvetica,bits_per_pixel=8,filename='contour_cips_mls_4pan_'+syr+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=2
   !p.charthick=2
   !y.thick=2
   !x.thick=2
endif
;
index=where(alb eq -999.)
if index(0) ne -1L then alb(index)=0
index=where(rad eq -999.)
if index(0) ne -1L then rad(index)=0
index=where(iwc eq -999.)
if index(0) ne -1L then iwc(index)=0
index=where(mls_tp eq -999.)
if index(0) ne -1L then mls_tp(index)=0.
index=where(mls_tp eq -999.)
if index(0) ne -1L then mls_tp(index)=0.
mls_tpsm=0.*mls_tp
fillit,mls_tp,mls_tpsm
mls_tp=mls_tpsm

mls_h2osm=0.*mls_h2o
fillit,mls_h2o,mls_h2osm
mls_h2o=mls_h2osm

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nbin2=nbin/2
latlosave=latlo
latlo(nbin2:nbin-1)=latlo(nbin2)+findgen(nbin2)
index=where(latlo mod 10 eq 0,nyticks)
tmin=1.	
tmax=40.
ival=1.
if shem eq 'south' then ival=-1.
tlevel=tmin+findgen(25)
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/nlvls)*mcolor
contour,smooth(alb,3),dfs,latlo,xtitle='DFS 20'+syr,ytitle='PM    Latitude    AM',color=0,title='CIPS Albedo > '+salb,charsize=1.25,ytickv=latlo(index),$
        yticks=nyticks-1,ytickname=strcompress(long(ival*latlosave(index)),/remove_all),/fill,levels=tlevel,c_color=col1,xrange=[-51,71]
if shem eq 'north' then contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[140],color=.9*mcolor,thick=3
if shem eq 'south' then contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[150],color=.9*mcolor,thick=3
imin=min(tlevel)
imax=max(tlevel)
xmnb=xmx+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,title='Garys',charsize=1.25,charthick=1.5
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tmin=20.
tmax=40.
tlevel=tmin+2.*findgen(25)
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/nlvls)*mcolor
contour,smooth(rad,3),dfs,latlo,xtitle='DFS 20'+syr,ytitle='PM    Latitude    AM',color=0,title='CIPS Radius > '+salb,charsize=1.25,ytickv=latlo(index),$
        yticks=nyticks-1,ytickname=strcompress(long(ival*latlosave(index)),/remove_all),/fill,levels=tlevel,c_color=col1,xrange=[-51,71]
if shem eq 'north' then contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[140],color=.9*mcolor,thick=3
if shem eq 'south' then contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[150],color=.9*mcolor,thick=3
imin=min(tlevel)
imax=max(tlevel)
xmnb=xmx+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,title='nm',charsize=1.25,charthick=1.5
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tmin=120.
if shem eq 'south' then tmin=130.
tmax=40.
tlevel=tmin+2.*findgen(25)
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/nlvls)*mcolor
contour,smooth(mls_tp,3),dfs,latlo,xtitle='DFS 20'+syr,ytitle='PM    Latitude    AM',color=0,title='MLS Temperature',charsize=1.25,ytickv=latlo(index),$
        yticks=nyticks-1,ytickname=strcompress(long(ival*latlosave(index)),/remove_all),/fill,levels=tlevel,c_color=col1,xrange=[-51,71]
contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[150],color=mcolor,thick=7
contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[140],color=mcolor,thick=7
contour,smooth(mls_tp,3),dfs,latlo,/overplot,/follow,level=[130],color=mcolor,thick=7
imin=min(tlevel)
imax=max(tlevel)
xmnb=xmx+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,title='K',charsize=1.25,charthick=1.5
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tmin=1.
tmax=40.
tinc=0.25
if shem eq 'south' then tinc=0.2
tlevel=tmin+tinc*findgen(25)
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/nlvls)*mcolor
contour,smooth(mls_h2o,3),dfs,latlo,xtitle='DFS 20'+syr,ytitle='PM    Latitude    AM',color=0,title='MLS Water Vapor',charsize=1.25,ytickv=latlo(index),$
        yticks=nyticks-1,ytickname=strcompress(long(ival*latlosave(index)),/remove_all),/fill,levels=tlevel,c_color=col1,xrange=[-51,71]
;contour,smooth(mls_h2o,3),dfs,latlo,/overplot,/follow,level=[5],color=mcolor,thick=5
imin=min(tlevel)
imax=max(tlevel)
xmnb=xmx+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,title='ppmv',charsize=1.25,charthick=1.5
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim contour_cips_mls_4pan_'+syr+'.ps -rotate -90 contour_cips_mls_4pan_'+syr+'.jpg'
;  spawn,'rm -f contour_cips_mls_4pan_'+syr+'.ps'
endif
end
