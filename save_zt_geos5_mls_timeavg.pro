;
; plot and save user entire season of daily GEOS and MLS at specified latitude range 
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
yorig=[0.7,0.4,0.1]
xorig=[0.1,0.1,0.1]
xlen=0.75
ylen=0.2
cbaryoff=0.075
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
mdir='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
month=['January','February','March','April','May','June','July','August','September','October','November','December']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
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
lstmn=11L & lstdy=1L & lstyr=2010L
ledmn=3L & leddy=1L & ledyr=2011L
lstday=0L & ledday=0L
firstdate=string(FORMAT='(i4)',lstyr)
lastdate=string(FORMAT='(i4)',ledyr)
daterange=firstdate+'-'+lastdate

solday=355L     ; doy of SH summer solstice
;
; get date range
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L	; number of days
dfs=fltarr(kday)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit

      dfs(icount)=iday-solday
      if iday lt 200 then dfs(icount)=iday-solday+365.
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; read daily zonal mean data
;
      dum=findfile(mdir+'yz_geos5_mls_'+sdate+'.sav')
      if dum(0) eq '' then goto,jumpday
      restore,mdir+'yz_geos5_mls_'+sdate+'.sav'	; alat,altitude,tbar,cobar,h2obar,tbarz,ubarz,vbarz
      nz=n_elements(altitude)
      if icount eq 0L then begin
         tbaravg=fltarr(kday,nz)
         cobaravg=fltarr(kday,nz)
         h2obaravg=fltarr(kday,nz)
         tbarzavg=fltarr(kday,nz)
         ubarzavg=fltarr(kday,nz)
         vbarzavg=fltarr(kday,nz)
print,alat
rlat=60.
read,'Enter desired latitude ',rlat
index=where(rlat eq alat)
ilat=index(0)
slat=strcompress(rlat,/remove_all)

      endif
      if icount gt 0L then begin
         tbaravg(icount,*)=tbar(ilat,*)
         cobaravg(icount,*)=cobar(ilat,*)
         h2obaravg(icount,*)=h2obar(ilat,*)
         tbarzavg(icount,*)=tbarz(ilat,*)
         ubarzavg(icount,*)=ubarz(ilat,*)
         vbarzavg(icount,*)=vbarz(ilat,*)
      endif
jumpday:
icount=icount+1
goto,jump
;
; plot and save avg file
;
plotit:
tbar=tbaravg
cobar=cobaravg
h2obar=h2obaravg
tbarz=tbarzavg
ubarz=ubarzavg
vbarz=vbarzavg
save,file=mdir+'zt_geos5_mls_'+daterange+'_'+slat+'.sav',rlat,slat,dfs,altitude,tbar,cobar,h2obar,tbarz,ubarz,vbarz
quick:
restore,mdir+'zt_geos5_mls_'+daterange+'_'+slat+'.sav'

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_geos5_mls_temp_'+daterange+'_'+slat+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot temperature
;
erase
xyouts,.2,.95,'Temperature '+daterange+'  Lat='+slat,charsize=2,charthick=2,/normal,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=110.+10.*findgen(nlvls)
index=where(tbarz eq 0.)
if index(0) ne -1L then tbarz(index)=-9999.
contour,tbarz,dfs,altitude,/noeras,xrange=[-50.,70.],yrange=[0.,100.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='DFS',title='GEOS-5',charthick=2
index=where(level gt 0.)
contour,tbarz,dfs,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,dfs,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,dfs,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
      xmnb=xorig(0)+xlen+0.1
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(level)
      omax=max(level)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0,charthick=2,charsize=1.5,title='(K)'
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(nlvls)
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
index=where(tbar eq 0.)
if index(0) ne -1L then tbar(index)=-9999.
contour,tbar,dfs,altitude,/noeras,xrange=[-50.,70.],yrange=[0.,100.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='DFS',title='MLS',charthick=2
index=where(level gt 0.)
contour,tbar,dfs,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,dfs,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,dfs,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
      xmnb=xorig(1)+xlen+0.1
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(level)
      omax=max(level)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0,charthick=2,charsize=1.5,title='(K)'
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tbarz
index=where(tbarz gt 0. and tbar gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tbarz(index)-tbar(index)
level=-50.+10.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,dfs,altitude,xrange=[-50.,70.],yrange=[0.,100.],$
        xtitle='DFS',charsize=1.5,levels=level,/cell_fill,$
        title='GEOS-5 - MLS',c_color=col2,color=0,min_value=-99.,charthick=2
index=where(level gt 0.)
contour,pdiff,dfs,altitude,/overplot,levels=level(index),color=0,/follow,min_value=-99.,c_labels=0*level(index)
index=where(level lt 0.)
contour,pdiff,dfs,altitude,/overplot,levels=level(index),color=mcolor,/follow,min_value=-99.,c_labels=0*level(index)
contour,pdiff,dfs,altitude,/overplot,levels=[0],color=0,thick=3,min_value=-99.
      xmnb=xorig(2)+xlen+0.1
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(level)
      omax=max(level)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0,charthick=2,charsize=1.5,title='(K)'
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(11)
      for j=0,11-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col2(j)
          y1=y1+dy
      endfor
loadct,39

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_geos5_mls_temp_'+daterange+'_'+slat+'.ps -rotate -90 zt_geos5_mls_temp_'+daterange+'_'+slat+'.jpg'
;  spawn,'/usr/bin/rm -f zt_geos5_mls_temp_'+daterange+'_'+slat+'.ps'
endif
end
