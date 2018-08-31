;
; GEOS-5 version
;
; plot daily d(theta)/dz in altitude and time
; use nc3 geos5 data
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
mcolor=icolmax
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.7
ylen=0.3
cbaryoff=0.08
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/'
stimes=[$
'_0000.V01.',$
'_0600.V01.',$
'_1200.V01.',$
'_1800.V01.']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=10
lstdy=1
lstyr=2003
ledmn=9
leddy=30
ledyr=2008
lstday=0
ledday=0
nlv=201L
altitude=findgen(nlv)
goto,plotit
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS-5 Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-5 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'+sdate+stimes(0)+'nc3'
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+ifile,nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      t2=0.*pv2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
;
; compute temperature
;
      t2=0.*u2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; compute lapse rate dT/dz=-(g/RT)*dT/dlnp
;
;     R=287.05
;     g=9.8
;     cp=1004.
;     dtdz=0*t2
;     for i=0,nc-1 do begin
;     for j=0,nr-1 do begin
;     for l=1,nth-1 do begin
;         lm1=l-1
;         dt=t2(j,i,l)-t2(j,i,lm1)
;         pl=p2(j,i,l)
;         plm1=p2(j,i,lm1)
;         dlnp=alog(pl/plm1)
;         dtdz(j,i,l)=-1000.*(g/(t2(j,i,l)*R))*dt/dlnp
;     endfor
;     endfor
;     endfor
;
; d(theta)/dz = (theta / T) * (dT/dz + g/Cp)
;
;     dthdz=0.*dtdz
;     for k=0,nth-1 do dthdz(*,*,k)=(th(k)/t2(*,*,k)) * (dtdz(*,*,k)+(g/cp))
;
; dtheta/dz = (th1-th0) / (z1-z0)
;
      dthdz=0.*t2
      for k=0,nth-1 do begin
          lp1=k-1
          lm1=k+1
          IF K EQ 0 then LP1=0
          IF K EQ NTH-1 then LM1=NTH-1
          for i=0,nc-1 do begin
          for j=0,nr-1 do begin
              DTHDz(j,i,K)=(TH(LP1)-TH(LM1))/(z2(j,i,LP1)-z2(j,i,LM1))
          endfor
          endfor
      endfor
;
; zonal means
;
      dtbarg=fltarr(nr,nth)
      tbarg=fltarr(nr,nth)
      ubarg=fltarr(nr,nth)
      zbarg=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          dtpts=reform(dthdz(j,*,k))
          tpts=reform(t2(j,*,k))
          upts=reform(u2(j,*,k))
          zpts=reform(z2(j,*,k))
          index=where(tpts ne 0.,nx)
          if index(0) ne -1L then begin
              dtbarg(j,k)=total(dtpts(index))/float(nx)
              tbarg(j,k)=total(tpts(index))/float(nx)
              ubarg(j,k)=total(upts(index))/float(nx)
              zbarg(j,k)=total(zpts(index))/float(nx)
          endif
      endfor
      endfor
;erase
;nlvls=21
;col1=1+indgen(nlvls)*icolmax/nlvls
;level=180.+5.*findgen(nlvls)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*icolmax/nlvls
;contour,tbarg,alat,th,/noeras,xrange=[-90,90],yrange=[500.,4000.],$
;      charsize=1.5,color=0,/fill,c_color=col1,levels=level,title=sdate
;nlvls=41
;level=-100.+10.*findgen(nlvls)
;index=where(level lt 0.)
;contour,dtbarg,alat,th,color=mcolor,/overplot,c_linestyle=5,levels=level(index)
;index=where(level gt 0.)
;contour,dtbarg,alat,th,color=0,/overplot,levels=level(index),c_labels=1+0*index
;contour,dtbarg,alat,th,color=0,/overplot,levels=[0],thick=2,c_labels=1+0*index
;stop
;
; first day
;
      if icount eq 0L then begin
         sdate0=sdate
         dttz=fltarr(kday,nth,nr)
         ttz=fltarr(kday,nth,nr)
         utz=fltarr(kday,nth,nr)
         sdates=strarr(kday)
;        rlat=0.
;        print,alat
;        read,' Enter Desired Latitude ',rlat
;        slat=string(format='(f5.1)',rlat)
      endif
      sdates(icount)=sdate
;
; interpolate GEOS temperature to SABER height surfaces
;
;tbarz=fltarr(nr,nlv)
;ubarz=fltarr(nr,nlv)
;for kk=0L,nlv-1L do begin
;    zz=altitude(kk)
;    for j=0L,nr-1L do begin
;zprof=reform(zbarg(j,*))
;if zz gt max(zprof) then goto,jumplev
;        for k=1L,nth-1L do begin
;            zup=zbarg(j,k-1) & zlw=zbarg(j,k)
;            if zup ge zz and zlw le zz then begin
;               zscale=(zup-zz)/(zup-zlw)
;               tbarz(j,kk)=tbarg(j,k-1)+zscale*(tbarg(j,k)-tbarg(j,k-1))
;               ubarz(j,kk)=ubarg(j,k-1)+zscale*(ubarg(j,k)-ubarg(j,k-1))
;
;;print,zlw,zz,zup,zscale
;;print,tbarg(j,k),tbarz(j,kk),tbarg(j,k-1)
;;stop
;            endif
;         endfor
;      endfor
;jumplev:
;endfor
;
; retain zonal mean temperature each day
;
      for k=0L,nlv-1L do begin
      for j=0L,nr-1L do begin
          utz(icount,*,*)=transpose(reform(ubarg))
          ttz(icount,*,*)=transpose(reform(tbarg))
          dttz(icount,*,*)=transpose(reform(dtbarg))
      endfor
      endfor

icount=icount+1L
goto,jump
;
; dtdz should really be called dthdztz
;
saveit:
save,file='zt_geos5_dthdz.sav',dttz,utz,ttz,alat,th,nr,nth,kday,sdates
plotit:
restore,'zt_geos5_dthdz.sav'

for j=nr-1L,nr/2,-1L do begin
slat=strcompress(string(format='(f6.2)',alat(j)),/remove_all)
index=where(alat eq -1.*alat(j))
kk=index(0)
slat2=strcompress(string(format='(f6.2)',alat(kk)),/remove_all)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_geos5_dthdz_'+slat+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; date labels
;
syear=strmid(sdates,0,4)
smon=strmid(sdates,4,2)
sday=strmid(sdates,6,2)
xindex=where(sday eq '15' and smon eq '01',nxticks)
xlabs=smon(xindex)+'/'+syear(xindex)
;
; plot zonal mean zonal wind
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=16
col1=1+indgen(nlvls)*icolmax/nlvls
level=10.*findgen(nlvls)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
;dthdzlat=smooth(reform(dttz(*,*,j)),7)
tlat=smooth(reform(ttz(*,*,j)),7)
for i=0,3 do dthdzlat=smooth(reform(dttz(*,*,j)),7)
;for i=0,9 do tlat=smooth(reform(ttz(*,*,j)),7)
;for i=0,3 do ulat=smooth(reform(utz(*,*,j)),7)
contour,dthdzlat,1.+findgen(kday),th,/noeras,xrange=[1.,kday-30],yrange=[500.,4000.],$
      charsize=1.5,color=0,ytitle='Theta (K)',title='!9d!X(theta)/!9d!Xz at '+slat+' N',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
for i=0,3 do dthdzlat=smooth(reform(dttz(*,*,j)),7)
contour,dthdzlat,1.+findgen(kday),th,levels=[25,50,75,100],color=0,/follow,/overplot,c_labels=1+0*indgen(8)
;contour,tlat,1.+findgen(kday),th,levels=[280],color=mcolor*.9,/follow,/overplot,c_labels=0,thick=5
;contour,tlat,1.+findgen(kday),th,levels=[270],color=mcolor*.8,/follow,/overplot,c_labels=0,thick=2
;contour,tlat,1.+findgen(kday),th,levels=[250],color=mcolor*.75,/follow,/overplot,c_labels=0,thick=2

;contour,ulat,1.+findgen(kday),th,levels=[0],color=0,/follow,/overplot,c_labels=0,thick=2
;contour,ulat,1.+findgen(kday),th,levels=[25,50,75],color=mcolor,/follow,/overplot,c_labels=0
;contour,ulat,1.+findgen(kday),th,levels=[-75,-50,-25],color=0,/follow,/overplot,c_labels=0

;xbox1=where(sday eq '01' and smon eq '06',nxbox1)
;for i=0L,nxbox1-1 do begin
;    plots,xbox1(i),1500.
;    plots,xbox1(i),2500.,thick=4,/continue,color=mcolor
;endfor
;xbox2=where(sday eq '31' and smon eq '07',nxbox2)
;for i=0L,nxbox2-1 do begin
;    plots,xbox2(i),1500.
;    plots,xbox2(i),2500.,thick=4,/continue,color=mcolor
;endfor
xyouts,xmx+0.1,ymn-0.02,'Approximate Altitude (km)',color=0,/normal,orientation=90.
xyouts,kday+40L,500.,'20',color=0,/data
xyouts,kday+40L,1000.,'40',color=0,/data
xyouts,kday+40L,1500.,'45',color=0,/data
xyouts,kday+40L,2000.,'50',color=0,/data
xyouts,kday+40L,2500.,'55',color=0,/data
xyouts,kday+40L,3000.,'60',color=0,/data
xyouts,kday+40L,3500.,'65',color=0,/data
xyouts,kday+40L,4000.,'70',color=0,/data
;
; SH
;
syear=strmid(sdates,0,4)
smon=strmid(sdates,4,2)
sday=strmid(sdates,6,2)
xindex=where(sday eq '15' and smon eq '07',nxticks)
xlabs=smon(xindex)+'/'+syear(xindex)
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=16
col1=1+indgen(nlvls)*icolmax/nlvls
level=10.*findgen(nlvls)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
;dthdzlat=reform(dttz(*,*,kk))
tlat=smooth(reform(ttz(*,*,kk)),7)
for i=0,3 do dthdzlat=smooth(reform(dttz(*,*,kk)),7)
;for i=0,9 do tlat=smooth(reform(ttz(*,*,kk)),7)
;for i=0,3 do ulat=smooth(reform(utz(*,*,kk)),7)
contour,dthdzlat,1.+findgen(kday),th,/noeras,xrange=[1.,kday-30],yrange=[500.,4000.],$
      charsize=1.5,color=0,ytitle='Theta (K)',title='!9d!X(theta)/!9d!Xz at '+slat2+' S',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
for i=0,3 do dthdzlat=smooth(reform(dttz(*,*,kk)),7)
contour,dthdzlat,1.+findgen(kday),th,levels=[25,50,75,100],color=0,/follow,/overplot,c_labels=1+0*indgen(8)
;contour,tlat,1.+findgen(kday),th,levels=[280],color=mcolor*.9,/follow,/overplot,c_labels=0,thick=5
;contour,tlat,1.+findgen(kday),th,levels=[270],color=mcolor*.8,/follow,/overplot,c_labels=0,thick=2
;contour,tlat,1.+findgen(kday),th,levels=[250],color=mcolor*.75,/follow,/overplot,c_labels=0,thick=2

;contour,ulat,1.+findgen(kday),th,levels=[0],color=0,/follow,/overplot,c_labels=0,thick=2
;contour,ulat,1.+findgen(kday),th,levels=[25,50,75],color=mcolor,/follow,/overplot,c_labels=0
;contour,ulat,1.+findgen(kday),th,levels=[-75,-50,-25],color=0,/follow,/overplot,c_labels=0

;xbox1=where(sday eq '01' and smon eq '06',nxbox1)
;for i=0L,nxbox1-1 do begin
;    plots,xbox1(i),1500.
;    plots,xbox1(i),2500.,thick=4,/continue,color=mcolor
;endfor
;xbox2=where(sday eq '31' and smon eq '07',nxbox2)
;for i=0L,nxbox2-1 do begin
;    plots,xbox2(i),1500.
;    plots,xbox2(i),2500.,thick=4,/continue,color=mcolor
;endfor
xyouts,xmx+0.1,ymn-0.02,'Approximate Altitude (km)',color=0,/normal,orientation=90.
xyouts,kday+40L,500.,'20',color=0,/data
xyouts,kday+40L,1000.,'40',color=0,/data
xyouts,kday+40L,1500.,'45',color=0,/data
xyouts,kday+40L,2000.,'50',color=0,/data
xyouts,kday+40L,2500.,'55',color=0,/data
xyouts,kday+40L,3000.,'60',color=0,/data
xyouts,kday+40L,3500.,'65',color=0,/data
xyouts,kday+40L,4000.,'70',color=0,/data

imin=min(level)
imax=max(level)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K/km)'
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
   spawn,'convert -trim zt_geos5_dthdz_'+slat+'.ps -rotate -90 zt_geos5_dthdz_'+slat+'.jpg'
endif

endfor	; loop over latitude
end
