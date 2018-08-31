;
; plot user specified time average zonal mean GEOS-5 T, and MLS T, and difference
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

sver='v2.2'
sver='v3.3'

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
xorig=[0.1,0.4,0.7]
yorig=[0.4,0.4,0.4]
xlen=0.25
ylen=0.25
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
lstmn=12L & lstdy=21L & lstyr=2010L
ledmn=1L & leddy=10L & ledyr=2011L
lstday=0L & ledday=0L
;
; get date range
;
print, ' '
print, '      GEOS-5 Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
firstdate=strcompress(long(lstyr),/remove_all)+string(FORMAT='(i2.2)',lstmn)+string(FORMAT='(i2.2)',lstdy)
lastdate=strcompress(long(ledyr),/remove_all)+string(FORMAT='(i2.2)',ledmn)+string(FORMAT='(i2.2)',leddy)
daterange=firstdate+'-'+lastdate
dum=findfile(mdir+'yz_geos5_mls_'+daterange+'.sav')
;if dum(0) ne '' then goto,quick
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
      if dum(0) eq '' then goto,jump
      restore,mdir+'yz_geos5_mls_'+sdate+'.sav'	; alat,altitude,tbar,cobar,h2obar,tbarz,ubarz,vbarz
      if icount eq 0L then begin
         tbaravg=0.*tbar
         ntbaravg=0.*tbar
         cobaravg=0.*cobar
         ncobaravg=0.*cobar
         h2obaravg=0.*h2obar
         nh2obaravg=0.*h2obar
         tbarzavg=0.*tbarz
         ntbarzavg=0.*tbarz
         ubarzavg=0.*ubarz
         nubarzavg=0.*ubarz
         vbarzavg=0.*vbarz
         nvbarzavg=0.*vbarz
         icount=1L
      endif
      if icount gt 0L then begin
         index=where(tbar gt 0.)
         if index(0) ne -1L then tbaravg(index)=tbaravg(index)+tbar(index)
         if index(0) ne -1L then ntbaravg(index)=ntbaravg(index)+1.
         index=where(cobar gt 0.)
         if index(0) ne -1L then cobaravg(index)=cobaravg(index)+cobar(index)
         if index(0) ne -1L then ncobaravg(index)=ncobaravg(index)+1.
         index=where(h2obar gt 0.)
         if index(0) ne -1L then h2obaravg(index)=h2obaravg(index)+h2obar(index)
         if index(0) ne -1L then nh2obaravg(index)=nh2obaravg(index)+1.
         index=where(tbarz gt 0.)
         if index(0) ne -1L then tbarzavg(index)=tbarzavg(index)+tbarz(index)
         if index(0) ne -1L then ntbarzavg(index)=ntbarzavg(index)+1.
         index=where(ubarz ne 0.)
         if index(0) ne -1L then ubarzavg(index)=ubarzavg(index)+ubarz(index)
         if index(0) ne -1L then nubarzavg(index)=nubarzavg(index)+1.
         index=where(vbarz ne 0.)
         if index(0) ne -1L then vbarzavg(index)=vbarzavg(index)+vbarz(index)
         if index(0) ne -1L then nvbarzavg(index)=nvbarzavg(index)+1.
      endif
goto,jump
;
; plot and save avg file
;
plotit:
index=where(ntbaravg ne 0.)
if index(0) ne -1L then tbaravg(index)=tbaravg(index)/ntbaravg(index)
index=where(ncobaravg gt 0.)
if index(0) ne -1L then cobaravg(index)=cobaravg(index)/ncobaravg(index)
index=where(nh2obaravg gt 0.)
if index(0) ne -1L then h2obaravg(index)=h2obaravg(index)/nh2obaravg(index)
index=where(ntbarzavg gt 0.)
if index(0) ne -1L then tbarzavg(index)=tbarzavg(index)/ntbarzavg(index)
index=where(nubarzavg gt 0.)
if index(0) ne -1L then ubarzavg(index)=ubarzavg(index)/nubarzavg(index)
index=where(nvbarzavg gt 0.)
if index(0) ne -1L then vbarzavg(index)=vbarzavg(index)/nvbarzavg(index)
tbar=tbaravg
cobar=cobaravg
h2obar=h2obaravg
tbarz=tbarzavg
ubarz=ubarzavg
vbarz=vbarzavg
save,file=mdir+'yz_geos5_mls_'+daterange+'.sav',alat,altitude,tbar,cobar,h2obar,tbarz,ubarz,vbarz
quick:
restore,mdir+'yz_geos5_mls_'+daterange+'.sav'

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_geos5_mls_temp_'+daterange+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot temperature
;
erase
xyouts,.3,.7,'Avg Temperature '+daterange,charsize=2,charthick=2,/normal,color=0
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
contour,tbarz,alat,altitude,/noeras,xrange=[-90.,90.],yrange=[0.,100.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',xticks=6,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title='GEOS-5',charthick=2
index=where(level gt 0.)
contour,tbarz,alat,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,alat,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,alat,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(tbar eq 0.)
if index(0) ne -1L then tbar(index)=-9999.
contour,tbar,alat,altitude,/noeras,xrange=[-90.,90.],yrange=[0.,100.],charsize=1.5,color=0,$
      xticks=6,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title='MLS',charthick=2
index=where(level gt 0.)
contour,tbar,alat,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,alat,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,alat,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
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
contour,pdiff,alat,altitude,xrange=[-90.,90.],yrange=[0.,100.],$
        xticks=6,xtitle='Latitude',charsize=1.5,levels=level,/cell_fill,$
        title='GEOS-5 - MLS',c_color=col2,color=0,min_value=-99.,charthick=2
index=where(level gt 0.)
contour,pdiff,alat,altitude,/overplot,levels=level(index),color=0,/follow,min_value=-99.,c_labels=0*level(index)
index=where(level lt 0.)
contour,pdiff,alat,altitude,/overplot,levels=level(index),color=mcolor,/follow,min_value=-99.,c_labels=0*level(index)
contour,pdiff,alat,altitude,/overplot,levels=[0],color=0,thick=3,min_value=-99.
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle='(K)',color=0,xticks=n_elements(level)/2,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dx
endfor
loadct,39

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_geos5_mls_temp_'+daterange+'.ps -rotate -90 yz_geos5_mls_temp_'+daterange+'.jpg'
;  spawn,'/usr/bin/rm -f yz_geos5_mls_temp_'+daterange+'.ps'
endif
end
