;
; user enter desired latitude circle
; theta-time series of MLS temperature + CO
; superimpose anticyclone area from GEOS-5. Percent of area poleward of 60 N
; VLH March 12 2010
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
lat0=60.
read,'Enter latitude ',lat0
slat0=strcompress(long(lat0),/remove_all)
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
dir2='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
lstmn=1
lstdy=15
lstyr=2009
ledmn=3
leddy=1
ledyr=2009

lstmn=11
lstdy=1
lstyr=2008
ledmn=4
leddy=1
ledyr=2009

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
;     if iflag eq 1 then goto,skipday
      file1=dir+sdate+'_AVG.V01.nc3'
      dum1=findfile(file1)
      if dum1(0) ne '' then begin
         ncid=ncdf_open(file1)
         print,'opening ',file1
      endif
      if dum1(0) eq '' then begin
         file1=dir2+sdate+'_AVG.V01.nc3'
         dum1=findfile(file1)
         if dum1(0) ne '' then begin
            ncid=ncdf_open(file1)
            print,'opening ',file1
         endif
         if dum1(0) eq '' then goto,skipday
      endif
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
      file1=dir+sdate+'_AVG.V01.nc4'
      dum1=findfile(file1)
      if dum1(0) ne '' then begin
         ncid=ncdf_open(file1)
         print,'opening ',file1
      endif
      if dum1(0) eq '' then begin
         file1=dir2+sdate+'_AVG.V01.nc4'
         dum1=findfile(file1)
         if dum1(0) ne '' then begin
            ncid=ncdf_open(file1)
            print,'opening ',file1
         endif
         if dum1(0) eq '' then goto,skipday
      endif
      mark2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,mark2
      ncdf_close,ncid
      index=where(mark2 lt 0.)
      if index(0) ne -1L then mark2(index)=-1.*(mark2(index)/mark2(index))
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
; zonal mean temperature and geopotential height. marker
;
;     tbarg=fltarr(nr,nth)
;     zbarg=fltarr(nr,nth)
;     markbarg=fltarr(nr,nth)
;     for k=0L,nth-1L do begin
;     for j=0L,nr-1L do begin
;         tpts=reform(t2(j,*,k))
;         zpts=reform(z2(j,*,k))
;         markpts=reform(mark2(j,*,k))
;         index=where(tpts ne 0.,nx)
;         if index(0) ne -1L then begin
;             tbarg(j,k)=total(tpts(index))/float(nx)
;             zbarg(j,k)=total(zpts(index))/float(nx)
;             markbarg(j,k)=min(zpts(index))	; for now just set to -1 if there is an anticyclone at this lat
;         endif
;     endfor
;     endfor
;
; declare time altitude and area arrays on first day
;
      if gcount eq 0L then begin
         aarea_nh_zt=fltarr(kday,nth)
         varea_nh_zt=fltarr(kday,nth)
         az_nh_zt=fltarr(kday,nth)	; average height in anticyclones
         vz_nh_zt=fltarr(kday,nth)	; average height in the vortex
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
;
; percent area poleward of a certain latitude circle
;
index0=where(y2d ge lat0)
area0=total(area(index0))

         gcount=1
      endif
;
; loop over theta and compute area in vortex and in anticyclones. Use nc4 file for anticyclones (mark2)
; compute average height in anticyclones and in vortex for vertical interpolation to altitude levels
;
      for thlev=0,nth-1 do begin
          mark1=reform(mark2(*,*,thlev))
          z1=reform(z2(*,*,thlev))
          index=where(y2d gt lat0 and mark1 gt 0.0,nn)
          if index(0) ne -1 then begin
             varea_nh_zt(icount,thlev)=100.*total(area(index))/area0	; % of area poleward of lat0
             vz_nh_zt(icount,thlev)=total(z1(index))/float(nn)		; average height at those gridpoints
          endif
          index=where(y2d gt lat0 and mark1 lt 0.0,nn)
          if index(0) ne -1 then begin
             aarea_nh_zt(icount,thlev)=100.*total(area(index))/area0
             az_nh_zt(icount,thlev)=total(z1(index))/float(nn)
          endif
;print,icount,th(thlev),' vortex anticyclone % ',varea_nh_zt(icount,thlev),aarea_nh_zt(icount,thlev),$
;      vz_nh_zt(icount,thlev),az_nh_zt(icount,thlev)
      endfor
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
         mlspolarco_zt=fltarr(kday,nlv)
         mlspolartp_zt=fltarr(kday,nlv)
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
          if latitude(ii) ge lat0 then begin
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
      mlspolarco_zt(icount,*)=mlspolarco
      mlspolartp_zt(icount,*)=mlspolartp
skipmls:
skipday:
      icount=icount+1L
goto,jump

plotit:
;
; interpolate small gaps in time
;
for k=0,nlv-1 do begin
    dlev=reform(mlspolarco_zt(*,k))
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
    mlspolarco_zt(*,k)=dlev
    dlev=reform(mlspolartp_zt(*,k))
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
    mlspolartp_zt(*,k)=dlev
endfor
;
; interpolate GEOS to MLS altitudes
;
aarea_nh_zt_alt=fltarr(kday,nlv)
varea_nh_zt_alt=fltarr(kday,nlv)
for kk=10L,nlv-1L do begin
    zz=altitude(kk)
    for i=0L,kday-1L do begin
        zprof=reform(az_nh_zt(i,*))
;       if zz gt max(zprof) then goto,jumplev
        for k=1L,nth-1L do begin
            zup=zprof(k-1) & zlw=zprof(k)
            if zup ge zz and zlw le zz and zup ne 0. and zlw ne 0. then begin
               zscale=(zup-zz)/(zup-zlw)
               aarea_nh_zt_alt(i,kk)=aarea_nh_zt(i,k-1)+zscale*(aarea_nh_zt(i,k)-aarea_nh_zt(i,k-1))
               varea_nh_zt_alt(i,kk)=varea_nh_zt(i,k-1)+zscale*(varea_nh_zt(i,k)-varea_nh_zt(i,k-1))

;print,zlw,zz,zup,zscale
;print,aarea_nh_zt(i,k),aarea_nh_zt_alt(i,kk),aarea_nh_zt(i,k-1)
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
mlspolarco_zt=mlspolarco_zt*1.e6
;
; save temp, co, etc
;
save,file='zt_mls_temp+co+anticylone_percentarea_'+yearlab+'_'+slat0+'.sav',mlspolarco_zt,mlspolartp_zt,kday,$
     varea_nh_zt_alt,aarea_nh_zt_alt,altitude,sdate_all
quick:
restore,'zt_mls_temp+co+anticylone_percentarea_'+yearlab+'_'+slat0+'.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
xindex=where(sday eq '01',nxticks)
if minyear eq maxyear then xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_mls_temp+co+anticylone_percentarea_'+yearlab+'_'+slat0+'.ps'
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
index=where(mlspolarco_zt eq 0.)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarco_zt(index)=99.
tlevel=160.+5.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlspolartp_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,100.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',/cell_fill,c_color=col1,$
      levels=tlevel,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlspolartp_zt,1.+findgen(kday),altitude,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
;contour,mlspolarco_zt,1.+findgen(kday),altitude,levels=[0.05,0.1,0.25,0.5,1,2,5,7.5],color=mcolor,/follow,/overplot,$
;        max_value=99.,c_labels=[1,1,1,1,1,1,1,1],thick=7

index=where(aarea_nh_zt_alt eq 0.)
if index(0) ne -1L then aarea_nh_zt_alt(index)=0./0.
aarea_nh_zt_alt=smooth(aarea_nh_zt_alt,3,/nan)
loadct,0
;contour,aarea_nh_zt_alt,1.+findgen(kday),altitude,levels=[10],color=mcolor*.9,/follow,/overplot,thick=15
contour,aarea_nh_zt_alt,1.+findgen(kday),altitude,levels=[20],color=mcolor*.5,/follow,/overplot,thick=15
;contour,aarea_nh_zt_alt,1.+findgen(kday),altitude,levels=[30],color=mcolor*.3,/follow,/overplot,thick=15
contour,aarea_nh_zt_alt,1.+findgen(kday),altitude,levels=[40],color=mcolor*.3,/follow,/overplot,thick=15
contour,aarea_nh_zt_alt,1.+findgen(kday),altitude,levels=[60],color=0,/follow,/overplot,thick=15
loadct,39
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
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Average Temperature > '+slat0+'!uo!n N (K)'
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
       spawn,'convert -trim zt_mls_temp+co+anticylone_percentarea_'+yearlab+'_'+slat0+'.ps -rotate -90 zt_mls_temp+co+anticylone_percentarea_'+yearlab+'_'+slat0+'.jpg'
    endif
end
