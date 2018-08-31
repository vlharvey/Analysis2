;
; daily zonal mean GEOS-5 and MLS temperature and differences
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

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
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
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
dir='/atmos/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
lstmn=1L & lstdy=1L & lstyr=2013L
ledmn=1L & leddy=31L & ledyr=2013L
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
      if ndays gt ledday then stop,' Normal termination condition '
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
print,sdate
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.
      y2d=fltarr(nc,nr)
      x2d=fltarr(nc,nr)
      for i=0,nc-1 do y2d(i,*)=alat
      for i=0,nr-1 do x2d(*,i)=alon
      t2=0.*pv2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
;
; GEOS-5 zonal mean temperature and geopotential height
;
    tbarg=fltarr(nr,nth)
    zbarg=fltarr(nr,nth)
    for k=0L,nth-1L do begin
    for j=0L,nr-1L do begin
        tpts=reform(t2(j,*,k))
        zpts=reform(z2(j,*,k))
        index=where(tpts ne 0.,nx)
        if index(0) ne -1L then begin
           tbarg(j,k)=total(tpts(index))/float(nx)
           zbarg(j,k)=total(zpts(index))/float(nx)
        endif
    endfor
    endfor
    yyyymmdd=sdate
;
; restore MLS tpd_mls_v2.2_20080614.sav file
;
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[7]
; DENTOT          FLOAT     = Array[3494, 121]
; ID              STRING    = Array[3494]
; PRESSURE        FLOAT     = Array[3494, 121]
; TEMPERATURE     FLOAT     = Array[3494, 121]
; TEMPERATURE_ERROR
; TEMPERATURE_MASK
;
    dum=findfile(mdir+'tpd_mls_v3.3_'+yyyymmdd+'.sav')
    if dum(0) eq '' then goto,jump
    restore,mdir+'tpd_mls_v3.3_'+yyyymmdd+'.sav'
    restore,mdir+'cat_mls_v3.3_'+yyyymmdd+'.sav'   ; latitude
    nlv=n_elements(altitude)
;
; compute zonal mean T, Z
;
    dy=alat(1)-alat(0)
    tbar=fltarr(nr,nlv)
    nbar=lonarr(nr,nlv)
    for ii=0L,n_elements(id)-1L do begin
        tmask_prof=reform(TEMPERATURE_MASK(ii,*))
        good=where(tmask_prof ne -99.,ngood)
        if good(0) ne -1L then begin
           ymean=latitude(ii)
           for j=0L,nr-1L do begin
               if ymean ge alat(j)-dy/2. and ymean lt alat(j)+dy/2. then begin
                  tbar(j,good)=tbar(j,good)+temperature(ii,good)
                  nbar(j,good)=nbar(j,good)+1L
               endif
           endfor
        endif
    endfor

good=where(nbar gt 0L)
if good(0) ne -1L then tbar(good)=tbar(good)/float(nbar(good))
;
; interpolate GEOS-5 temperature to height surfaces
;
tbarz=fltarr(nr,nlv)
for kk=0L,nlv-1L do begin
    zz=altitude(kk)
    for j=0L,nr-1L do begin
        zprof=reform(zbarg(j,*))
        if zz gt max(zprof) then goto,jumplev
        for k=1L,nth-1L do begin
            zup=zbarg(j,k-1) & zlw=zbarg(j,k)
            if zup ge zz and zlw le zz then begin
               zscale=(zup-zz)/(zup-zlw)
               tbarz(j,kk)=tbarg(j,k-1)+zscale*(tbarg(j,k)-tbarg(j,k-1))
;print,zlw,zz,zup,zscale
;print,tbarg(j,k),tbarz(j,kk),tbarg(j,k-1)
;stop
            endif
         endfor
      endfor
jumplev:
endfor
;
; save monthly mean gridded MLS data
;
;save,file=mdir+'geos5_mls_v2.2_'+sdate+'_temp.sav',alat,altitude,tbar,tbarz
plotit:
;restore,mdir+'geos5_mls_v2.2_'+sdate+'_temp.sav'

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_geos5_mls_temp_'+sdate+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot temperature
;
erase
xyouts,.3,.7,sdate+' Temperature',charsize=2,charthick=2,/normal,color=0
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
   spawn,'convert -trim yz_geos5_mls_temp_'+sdate+'.ps -rotate -90 yz_geos5_mls_temp_'+sdate+'.jpg'
   spawn,'/usr/bin/rm -f yz_geos5_mls_temp_'+sdate+'.ps'
endif
goto,jump
end
