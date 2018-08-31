;
; GEOS-5 version plus MLS
; compute daily average Arctic Temp and Mark
; plot altitude-time sections
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3

sver='v2.2'

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15,0.55,0.55]
yorig=[0.55,0.15,0.55,0.15]
xlen=0.3
ylen=0.3
cbaryoff=0.06
cbarydel=0.01
!noeras=1
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
dirm='/aura6/data/MLS_data/Datfiles_SOSST/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
smon=['J','J','J','J','F','F','F','F','M','M','M','M','M','A','M','J','J','A','S','O','N','D']
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
lstmn=11L & lstdy=1L & lstyr=9L
ledmn=2L & leddy=1L & ledyr=10L
lstday=0L & ledday=0L
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
minyear=lstyr
maxyear=ledyr
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
goto,quick

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=-1L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      ifile=mon(imn-1)+sdy+'_'+uyr
      rd_geos5_nc3,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
               pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      kcount=kcount+1L
      if iflag eq 1 then goto,jump
;
; read new marker field
;
      ncid=ncdf_open(dir+sdate+'_AVG.V01.nc4')
      mark2new=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,mark2new
      ncdf_close,ncid

      if icount eq 0L then begin
         mls_temp_nh_zt=fltarr(kday,nth)
         mls_temp_sh_zt=fltarr(kday,nth)
         mls_co_nh_zt=fltarr(kday,nth)
         mls_co_sh_zt=fltarr(kday,nth)
         temp_nh_zt=fltarr(kday,nth)
         temp_sh_zt=fltarr(kday,nth)
         u_nh_zt=fltarr(kday,nth)
         u_sh_zt=fltarr(kday,nth)
         sdate_all=strarr(kday)
         time_all=fltarr(kday)
         icount=1
      endif
      sdate_all(kcount)=sdate
      time_all(kcount)=0.
;
; temperature
;
      temp2=0.*p2
      for k=0L,nth-1L do temp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^.286
;
; average 60-90N
;
      x2d=fltarr(nr,nc)
      y2d=fltarr(nr,nc)
      for i=0L,nc-1 do y2d(*,i)=alat
      for j=0L,nr-1 do x2d(j,*)=alon
      for k=0L,nth-1L do begin
          temp1=reform(temp2(*,*,k))
          u1=reform(u2(*,*,k))
          nhindex=where(y2d ge 60. and temp1 ne 0.,nnh)
          shindex=where(y2d le -60. and temp1 ne 0.,nsh)
          temp_nh_zt(kcount,k)=mean(temp1(nhindex))
          temp_sh_zt(kcount,k)=mean(temp1(shindex))
          u_nh_zt(kcount,k)=mean(u1(nhindex))
          u_sh_zt(kcount,k)=mean(u1(shindex))
      endfor
;
; MLS daily polar temp interpolated to GEOS theta surfaces
;
    dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
    if dum(0) eq '' then goto,jump
    restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
    restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
    restore,dirm+'co_mls_'+sver+'_'+sdate+'.sav'              ; mix
;   restore,dirm+'h2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
;   restore,dirm+'n2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
    nz=n_elements(altitude)
    nthlev=n_elements(thlev)
    mprof=n_elements(longitude)
    mlev=n_elements(altitude)
    muttime=time
    mlat=latitude
    mlon=longitude
    bad=where(mask eq -99.)
    if bad(0) ne -1L then mix(bad)=-99.
    good=where(mix ne -99.)
    if good(0) eq -1L then goto,jump
    mco=mix*1.e6
    mtemp=temperature
    mpress=pressure
;
; eliminate bad uttimes tropics
;
    index=where(muttime gt 0. and abs(mlat) ge 60.,mprof)
    if index(0) eq -1L then goto,jump
    muttime=reform(muttime(index))
    mlat=reform(mlat(index))
    mlon=reform(mlon(index))
    mtemp=reform(mtemp(index,*))
    mpress=reform(mpress(index,*))
    mco=reform(mco(index,*))
    mtheta=mtemp*(1000./mpress)^0.286
    index=where(mtemp lt 0.)
    if index(0) ne -1L then mtheta(index)=-99.
;
; interpolate CO data to GEOS-5 theta surfaces
;
    mco_th=fltarr(mprof,nth)
    mtemp_th=fltarr(mprof,nth)
    for k=nth-1L,0L,-1L do begin
        zlev=th(k)
        for iprof=0L,mprof-1L do begin
            for kk=2L,nz-2L do begin

if mco(iprof,kk) ne -9.90000e+07 and mco(iprof,kk+1) ne -9.90000e+07 then begin
if mtheta(iprof,kk) lt zlev and mtheta(iprof,kk+1) ge zlev then begin
   zscale=(mtheta(iprof,kk+1)-zlev)/(mtheta(iprof,kk+1)-mtheta(iprof,kk))
   mco_th(iprof,k)= mco(iprof,kk+1)+zscale*(mco(iprof,kk)-mco(iprof,kk+1))
;
;;print,mtheta(iprof,kk),zlev,mtheta(iprof,kk+1),zscale
;;print,mco(iprof,kk),mco_th(iprof,k),mco(iprof,kk+1)
;;stop
endif
endif
if mtemp(iprof,kk) gt 0. and mtemp(iprof,kk+1) gt 0. then begin
if mtheta(iprof,kk) lt zlev and mtheta(iprof,kk+1) ge zlev then begin
   zscale=(mtheta(iprof,kk+1)-zlev)/(mtheta(iprof,kk+1)-mtheta(iprof,kk))
   mtemp_th(iprof,k)= mtemp(iprof,kk+1)+zscale*(mtemp(iprof,kk)-mtemp(iprof,kk+1))
endif
endif
            endfor
        endfor
;   endfor
;
; MLS polar mean temps
;
;   for k=0L,nth-1L do begin
        mtemp1=reform(mtemp_th(*,k))
        mco1=reform(mco_th(*,k))
        nhindex=where(mlat ge 60. and mtemp1 gt 0.,nnh)
        shindex=where(mlat le -60. and mtemp1 gt 0.,nsh)
        mls_co_nh_zt(kcount,k)=mean(mco1(nhindex))
        mls_co_sh_zt(kcount,k)=mean(mco1(shindex))
        mls_temp_nh_zt(kcount,k)=mean(mtemp1(nhindex))
        mls_temp_sh_zt(kcount,k)=mean(mtemp1(shindex))
    endfor
goto, jump

plotit:
;
; interpolate small gaps in time
;
for k=0,nth-1 do begin
    dlev=reform(mls_temp_nh_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    mls_temp_nh_zt(*,k)=dlev
    dlev=reform(mls_temp_sh_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump2
               endif
           endfor
jump2:
        endif
    endfor
    mls_temp_sh_zt(*,k)=dlev
endfor
;
; year date label
;
syear=strmid(sdate_all,0,4)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
; comment out if generating quick plot
;
SAVE, /VARIABLES, FILENAME = 'zt_geos5+mls_polarT_'+yearlab+'.sav'
quick:
restore,'zt_geos5+mls_polarT_'+yearlab+'.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
; smooth u
;
u_nh_zt=smooth(u_nh_zt,7)
;u_sh_zt=smooth(u_sh_zt,7)
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   mpsym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_geos5+mls_polarT_'+yearlab+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
slab=['NH','SH']
;
; plot time-altitude sections of Temperature and zonal wind speed in NH and SH
;
xyouts,.4,.95,yearlab,/normal,color=0,charsize=3,charthick=3
;
; print end date
;
maxdate=max(long(sdate_all))
smaxdate=strcompress(maxdate,/remove_all)
datelab=strmid(smaxdate,4,2)+'/'+strmid(smaxdate,6,2)
xyouts,.35,.9,'Last Day '+datelab,color=0,/normal,charsize=3,charthick=3

tlevel=190.+5.*findgen(17)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*mcolor/nlvls
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
ylab='Theta (K)'
contour,temp_nh_zt,findgen(kday),th,charsize=1.5,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=ylab,levels=tlevel,c_color=col1,/cell_fill,color=0,title='GEOS-5 '+slab(0),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
loadct,0
contour,temp_nh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
loadct,0
;contour,temp_nh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
ulevel=[-90.,-60.,-30.,-10.,10.,30.,60.,90.]
index=where(ulevel lt 0.)
contour,u_nh_zt,findgen(kday),th,/noeras,/overplot,levels=ulevel(index),color=0,thick=3,c_linestyle=5
index=where(ulevel gt 0.)
contour,u_nh_zt,findgen(kday),th,/noeras,/overplot,levels=ulevel(index),color=mcolor,thick=3
;
; MLS
;
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,mls_temp_nh_zt,findgen(kday),th,charsize=1.5,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=' ',levels=tlevel,c_color=col1,/cell_fill,color=0,title='MLS '+slab(0),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
loadct,0
contour,mls_temp_nh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,temp_sh_zt,findgen(kday),th,charsize=1.5,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=ylab,levels=tlevel,c_color=col1,/cell_fill,color=0,title='GEOS-5 '+slab(1),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
loadct,0
contour,temp_sh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
loadct,0
;contour,temp_sh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
index=where(ulevel lt 0.)
contour,u_sh_zt,findgen(kday),th,/noeras,/overplot,levels=ulevel(index),color=0,thick=3,c_linestyle=5
index=where(ulevel gt 0.)
contour,u_sh_zt,findgen(kday),th,/noeras,/overplot,levels=ulevel(index),color=mcolor,thick=3
;
; MLS SH
;
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,mls_temp_sh_zt,findgen(kday),th,charsize=1.5,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=' ',levels=tlevel,c_color=col1,/cell_fill,color=0,title='MLS '+slab(1),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
loadct,0
contour,mls_temp_sh_zt,findgen(kday),th,levels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xorig(1),xorig(3)+xlen,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='60-90 Mean Temperature (K)',charsize=2,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_geos5+mls_polarT_'+yearlab+'.ps -rotate -90 zt_geos5+mls_polarT_'+yearlab+'.jpg'
;  spawn,'/usr/bin/rm zt_geos5+mls_polarT_'+yearlab+'.ps'
endif
end
