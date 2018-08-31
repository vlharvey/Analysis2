;
; save daily zonal mean anticyclone area, mean height in anticyclones, MLS temp and CO
; VLH March 11 2010
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.20]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.07
cbarydel=0.01
!NOERAS=-1
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
ks=1.931853d-3
ecc=0.081819
gamma45=9.80
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
lstmn=1
lstdy=21
lstyr=2008
ledmn=3
leddy=1
ledyr=2008
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
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
;goto,quick

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L
gcount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      sdate_all(icount)=sdate
;
; read GEOS data
;
;     rd_geos5_nc3,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
;              pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
;     if iflag eq 1 then goto,jump
      file1=dir+sdate+'_AVG.V01.nc3'
      dum1=findfile(file1)
      if dum1(0) ne '' then begin
         ncid=ncdf_open(file1)
         print,'opening ',file1
      endif
      if dum1(0) eq '' then goto,jump
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      p2=fltarr(nr,nc,nth)
      msf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,4,p2
      ncdf_varget,ncid,5,msf2
      ncdf_close,ncid
;
; read new marker field
;
      ncid=ncdf_open(dir+sdate+'_AVG.V01.nc4')
      mark2new=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,mark2new
      ncdf_close,ncid
      index=where(mark2new lt 0.)
      if index(0) ne -1L then mark2new(index)=-1.*(mark2new(index)/mark2new(index))
;
; temperature and geopotential height
;
      t2=0.*p2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      g2=(msf2-1004.*t2)/(9.86*1000.)
;
; convert geopotential to geometric height
;
      z2=0.*g2
      for j=0L,nr-1L do begin
          sin2=sin( (alat(j)*dtr)^2.0 )
          numerator=1.0+ks*sin2
          denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
          gammas=gamma45*(numerator/denominator)
          r=6378.137/(1.006803-(0.006706*sin2))
          z2(j,*,*)=(r*g2(j,*,*))/ ( (gammas/gamma45)*r - g2(j,*,*) )
      endfor
;
; declare time altitude and area arrays on first day
;
      if gcount eq 0L then begin
         aarea_nh_yzt=fltarr(kday,nr,nth)
         varea_nh_yzt=fltarr(kday,nr,nth)
         az_nh_yzt=fltarr(kday,nr,nth)	; average height in anticyclones
         vz_nh_yzt=fltarr(kday,nr,nth)	; average height in the vortex
         x2d=fltarr(nr,nc)
         y2d=fltarr(nr,nc)
         for i=0L,nc-1 do y2d(*,i)=alat
         for j=0L,nr-1 do x2d(j,*)=alon
         area=0.*y2d
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(j,*)=dx*hy    ; area of each grid point
         endfor
         gcount=1
      endif
;
; zonal mean temperature and geopotential height. marker
;
      tbarg=fltarr(nr,nth)
      zbarg=fltarr(nr,nth)
      amarkbarg=fltarr(nr,nth)
      vmarkbarg=fltarr(nr,nth)
;
; percent area poleward of a certain latitude circle
;
aarea_prof=fltarr(nth)
varea_prof=fltarr(nth)
lat0=60.
index0=where(y2d ge lat0)
area0=total(area(index0))

      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          tpts=reform(t2(j,*,k))
          zpts=reform(z2(j,*,k))
          markpts=reform(mark2new(j,*,k))
          areapts=reform(area(j,*))
          index=where(tpts ne 0.,nx)
          if index(0) ne -1L then begin
              tbarg(j,k)=total(tpts(index))/float(nx)
              zbarg(j,k)=total(zpts(index))/float(nx)
          endif
;
; area in the vortex and anticyclones at each latitude and altitude
;
          index=where(tpts ne 0. and markpts gt 0.,nn)
          if index(0) ne -1L then begin
             vmarkbarg(j,k)=100.*float(nn)/float(nc)		; % of longitude points
             varea_nh_yzt(icount,j,k)=total(areapts(index))
             vz_nh_yzt(icount,j,k)=total(zpts(index))/float(nx)
          endif
          index=where(tpts ne 0. and markpts lt 0.,nn)
          if index(0) ne -1L then begin
             amarkbarg(j,k)=100.*float(nn)/float(nc)
             aarea_nh_yzt(icount,j,k)=total(areapts(index))
             az_nh_yzt(icount,j,k)=total(zpts(index))/float(nx)
          endif
      endfor
;
; what % of polar cap is occupied by anticyclones?
;
          mark1=reform(mark2new(*,*,k))

          index=where(y2d ge lat0 and mark1 lt 0.0,nn)
          if index(0) ne -1 then aarea_prof(k)=100.*total(area(index))/area0
          index=where(y2d ge lat0 and mark1 gt 0.0,nn)
          if index(0) ne -1 then varea_prof(k)=100.*total(area(index))/area0

      endfor
;
; check zonal means
;
erase
xyouts,.425,.925,sdate,/normal,color=0,charsize=3
set_viewport,0.2,0.5,0.55,0.85
!type=2^2+2^3+2^7       ; ticks outward
index=where(tbarg eq 0.)
if index(0) ne -1L then tbarg(index)=0./0.
tlevel=160.+5.*findgen(29)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,tbarg,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],title='Temperature',levels=tlevel,$
      charsize=1.5,color=0,ytitle='Theta (K)',/cell_fill,c_color=col1,xticks=6,xtitle='Latitude'
contour,tbarg,alat,th,levels=tlevel,color=0,/follow,/overplot,c_labels=1+0*intarr(nlvls)

set_viewport,0.2,0.5,0.15,0.45
!type=2^2+2^3+2^7       ; ticks outward
index=where(zbarg eq 0.)
if index(0) ne -1L then zbarg(index)=0./0.
tlevel=10.+2.5*findgen(29)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,zbarg,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],levels=tlevel,title='Altitude',$
      charsize=1.5,color=0,ytitle='Theta (K)',/cell_fill,c_color=col1,xticks=6,xtitle='Latitude'
contour,zbarg,alat,th,levels=tlevel,color=0,/follow,/overplot,c_labels=1+0*intarr(nlvls)

set_viewport,0.6,0.9,0.55,0.85
!type=2^2+2^3+2^7       ; ticks outward
index=where(vmarkbarg eq 0.)
if index(0) ne -1L then vmarkbarg(index)=0./0.
tlevel=[1.,5.,10.+10.*findgen(10)]
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,vmarkbarg,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],title='Vortex Area',levels=tlevel,$
      charsize=1.5,color=0,/cell_fill,c_color=col1,xticks=6,xtitle='Latitude'
contour,vmarkbarg,alat,th,levels=tlevel,color=0,/follow,/overplot,c_labels=1+0*intarr(nlvls)
svarea_prof=string(FORMAT='(f5.1)',varea_prof)
for k=0,nth-1 do xyouts,90.,th(k),svarea_prof(k),/data,color=0

set_viewport,0.6,0.9,0.15,0.45
!type=2^2+2^3+2^7       ; ticks outward
index=where(amarkbarg eq 0.)
if index(0) ne -1L then amarkbarg(index)=0./0.
contour,amarkbarg,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],levels=tlevel,title='Anticyclone Area',$
      charsize=1.5,color=0,/cell_fill,c_color=col1,xticks=6,xtitle='Latitude'
contour,amarkbarg,alat,th,levels=tlevel,color=0,/follow,/overplot,c_labels=1+0*intarr(nlvls)
saarea_prof=string(FORMAT='(f5.1)',aarea_prof)
for k=0,nth-1 do xyouts,90.,th(k),saarea_prof(k),/data,color=0


stop
;
; restore MLS CO on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v2.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'
      restore,mdir+'co_mls_v2.2_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v2.2_'+sdate+'.sav'
;
; apply mask
;
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlscomix=mix
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      mlscomix=mix
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         mlspolarco_yzt=fltarr(kday,nlv)
         mlspolartp_yzt=fltarr(kday,nlv)
         kcount=1
      endif
;
; compute polar CO and temp
;
      mlspolarco=fltarr(nlv)
      mlsncoprof=lonarr(nlv)
      mlspolartp=fltarr(nlv)
      mlsntpprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             co_prof=reform(mlscomix(ii,*))
             good=where(co_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolarco(good)=mlspolarco(good)+reform(co_prof(good))
                mlsncoprof(good)=mlsncoprof(good)+1L
             endif
             tp_prof=reform(temperature(ii,*))
             good=where(tp_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolartp(good)=mlspolartp(good)+reform(tp_prof(good))
                mlsntpprof(good)=mlsntpprof(good)+1L
             endif
          endif
      endfor
      good=where(mlsncoprof gt 0L)
      if good(0) ne -1L then mlspolarco(good)=mlspolarco(good)/float(mlsncoprof(good))
      good=where(mlsntpprof gt 0L)
      if good(0) ne -1L then mlspolartp(good)=mlspolartp(good)/float(mlsntpprof(good))
      mlspolarco_yzt(icount,*)=mlspolarco
      mlspolartp_yzt(icount,*)=mlspolartp
skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; interpolate small gaps in time
;
for k=0,nlv-1 do begin
    dlev=reform(mlspolarco_yzt(*,k))
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
    mlspolarco_yzt(*,k)=dlev
    dlev=reform(mlspolartp_yzt(*,k))
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
    mlspolartp_yzt(*,k)=dlev
endfor
;
; interpolate GEOS to MLS altitudes
;
aarea_nh_yzt_alt=fltarr(kday,nlv)
varea_nh_yzt_alt=fltarr(kday,nlv)
for kk=10L,nlv-1L do begin
    zz=altitude(kk)
    for i=0L,kday-1L do begin
        zprof=reform(az_nh_yzt(i,*))
        if zz gt max(zprof) then goto,jumplev
        for k=1L,nth-1L do begin
            zup=zprof(k-1) & zlw=zprof(k)
            if zup ge zz and zlw le zz and zup ne 0. and zlw ne 0. then begin
               zscale=(zup-zz)/(zup-zlw)
               aarea_nh_yzt_alt(i,kk)=aarea_nh_yzt(i,k-1)+zscale*(aarea_nh_yzt(i,k)-aarea_nh_yzt(i,k-1))
               varea_nh_yzt_alt(i,kk)=varea_nh_yzt(i,k-1)+zscale*(varea_nh_yzt(i,k)-varea_nh_yzt(i,k-1))

;print,zlw,zz,zup,zscale
;print,aarea_nh_yzt(i,k),aarea_nh_yzt_alt(i,kk),aarea_nh_yzt(i,k-1)
;stop
            endif
         endfor
      endfor
jumplev:
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
mlspolarco_yzt=mlspolarco_yzt*1.e6
;
; save temp, co, etc
;
;save,file='yzt_mls_temp+co+anticylone_area_'+yearlab+'.sav',mlspolarco_yzt,mlspolartp_yzt,kday,$
;     varea_nh_yzt_alt,aarea_nh_yzt_alt,altitude,sdate_all,alat
quick:
;restore,'yzt_mls_temp+co+anticylone_area_'+yearlab+'.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15' or sday eq '01',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yzt_mls_temp+co+anticylone_area_'+yearlab+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean temperature and CO + anticyclone area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolarco_yzt eq 0.)
if index(0) ne -1L then mlspolarco_yzt(index)=0./0.
mlspolarco_yzt=smooth(mlspolarco_yzt,7,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarco_yzt(index)=99.
tlevel=160.+5.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlspolartp_yzt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,100.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',/cell_fill,c_color=col1,$
      levels=tlevel,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlspolartp_yzt,1.+findgen(kday),altitude,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
;contour,mlspolarco_yzt,1.+findgen(kday),altitude,levels=[0.05,0.1,0.25,0.5,1,2,5,7.5],color=mcolor,/follow,/overplot,$
;        max_value=99.,c_labels=[1,1,1,1,1,1,1,1],thick=7

index=where(aarea_nh_yzt_alt eq 0.)
if index(0) ne -1L then aarea_nh_yzt_alt(index)=0./0.
aarea_nh_yzt_alt=smooth(aarea_nh_yzt_alt,3,/nan)
index=where(finite(aarea_nh_yzt_alt) eq 1)
if index(0) ne -1L then aarea_nh_yzt_alt(index)=aarea_nh_yzt_alt(index)/1000.	; thousands of km^2
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[100],color=0,/follow,/overplot,thick=8
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[200],color=mcolor*.2,/follow,/overplot,thick=8
;contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[250],color=mcolor*.35,/follow,/overplot,thick=8
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[300],color=mcolor*.75,/follow,/overplot,thick=8
;contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[350],color=mcolor*.875,/follow,/overplot,thick=8
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[400],color=mcolor*.85,/follow,/overplot,thick=8
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[500],color=mcolor*.95,/follow,/overplot,thick=8
contour,aarea_nh_yzt_alt,1.+findgen(kday),altitude,levels=[600],color=mcolor,/follow,/overplot,thick=8
xyouts,xmn+0.02,ymx-0.05,yearlab,/normal,color=mcolor,charsize=3,charthick=3
;
; print end date
;
maxdate=max(long(sdate_all))
smaxdate=strcompress(maxdate,/remove_all)
datelab=strmid(smaxdate,4,2)+'/'+strmid(smaxdate,6,2)
;xyouts,(xmx+xmn)/2.,ymn+.02,'Last Day '+datelab,color=0,/normal,charsize=2

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Avg Temperature > 60 N (K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim yzt_mls_temp+co+anticylone_area_'+yearlab+'.ps -rotate -90 yzt_mls_temp+co+anticylone_area_'+yearlab+'.jpg'
    endif
end
