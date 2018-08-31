;
; GEOS-4 version
;
; store and plot daily timeseries of polar (70 degrees latitude to the pole)
; stratopause pressure as a function of time-of-year in each hemisphere.
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_dat

loadct,38
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
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
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='timeseries_polar_stratopause_pressure_geos4.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
;goto,plotit
dir='/aura7/harvey/GEOS4_data/Datfiles/'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nlg=0l
nlat=0l
nlv=0l
lstmn=3
lstdy=31
lstyr=2004
ledmn=5
leddy=15
ledyr=2007
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS-4 Version '
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
nday=365L
nyear=ledyr-lstyr+1L
nhstratpavg=-99.+0.*fltarr(nday,nyear)
nhstratpmin=-99.+0.*fltarr(nday,nyear)
nhstratpmax=-99.+0.*fltarr(nday,nyear)
shstratpavg=-99.+0.*fltarr(nday,nyear)
shstratpmin=-99.+0.*fltarr(nday,nyear)
shstratpmax=-99.+0.*fltarr(nday,nyear)
sdate_save=strarr(nday,nyear)
syear=strarr(nyear)
for n=0L,nyear-1L do syear(n)=strcompress(lstyr+n,/remove_all)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then begin
         save,file='timeseries_polar_stratopause_pressure_geos4.sav',nday,nyear,syear,sdate_save,nhstratpavg,$
              shstratpavg,nhstratpmin,shstratpmin,nhstratpmax,shstratpmax
         goto,plotit
      endif
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
; skip leap day for simplicity
;
      if imn eq 2L and idy eq 29L then goto, jump
;
; accomodate leap days
;
      leapday=0L
      if iyr mod 4 eq 0L and iday gt 59L then leapday=1L
;
;***Read GEOS-4 data
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      ifile='DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'+sdate+'_1200.V01.dat'
      rd_geos5_dat,dir+ifile,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,zp,tp,up,vp,qp
      if iflg ne 0 then goto, jump
      zp=zp/1000.
;
; compute lapse rate dT/dz=-(g/RT)*dT/dlnp
;
      R=287.05
      g=9.8
      laps=0*tp
      for i=0,nlg-1 do begin
      for j=0,nlat-1 do begin
      for l=1,nlv-1 do begin
          lm1=l-1
          dt=tp(i,j,l)-tp(i,j,lm1)
          pl=p(l)
          plm1=p(lm1)
          dlnp=alog(pl/plm1)
          laps(i,j,l)=-1000.*(g/(tp(i,j,l)*R))*dt/dlnp
      endfor
      endfor
      endfor
;
; zonal mean temperature
;
      tpzm=fltarr(nlat,nlv)
      for k=0L,nlv-1L do begin
      for j=0L,nlat-1L do begin
          tpzm(j,k)=total(tp(*,j,k))/float(nlg)
      endfor
      endfor

      stratplev=fltarr(nlg,nlat)
      stratzlev=fltarr(nlg,nlat)
;erase
;!type=2^2+2^3
;set_viewport,xorig(0),xorig(0)+xlen,yorig(0),yorig(0)+ylen
;plot,findgen(10),p,/ylog,/nodata,yrange=[100.,0.01],xrange=[-10.,10.],/noeras,color=0
;plots,0.,100.
;plots,0.,0.01,color=0,thick=2,/continue
      for j=0,nlat-1 do begin
          if abs(alat(j)) ge 70. then begin
             for i=0,nlg-1 do begin
                 laps0=reform(laps(i,j,*))
                 zprof=reform(zp(i,j,*))
                 for k=nlv-2L,11,-1 do begin
;if alat(j) gt 0. then oplot,laps0,p,color=icolmax*.3
;if alat(j) lt 0. then oplot,laps0,p,color=icolmax*.9
                     if laps0(k) eq 0. then begin
;
; require negative below and positive above
;
                        if k eq nlv-2 then $
                          if laps0(k+1) lt 0. and laps0(k-1) gt 0. then stratplev(i,j)=p(k)
                        if k lt nlv-2 then $
                          if laps0(k+2) lt 0. and laps0(k+1) lt 0. and $
                             laps0(k-2) gt 0. and laps0(k-1) gt 0. then stratplev(i,j)=p(k)

;if alat(j) gt 0. then oplot,[laps0(k),laps0(k)],[p(k),p(k)],psym=8,symsize=2,color=0
                        goto,jumpz
                     endif
                     if laps0(k) lt 0. and laps0(k-1) gt 0. then begin
                        scale=(laps0(k)-0.)/(laps0(k)-laps0(k-1))
                        if k eq nlv-2 then $
                           if laps0(k+1) lt 0. and laps0(k-1) gt 0. then $
                           stratplev(i,j)=exp( alog(p(k))+scale*(alog(p(k-1))-alog(p(k))) )

                        if k lt nlv-2 then $
                          if laps0(k+2) lt 0. and laps0(k+1) lt 0. and $
                             laps0(k-2) gt 0. and laps0(k-1) gt 0. then $
                             stratplev(i,j)=exp( alog(p(k))+scale*(alog(p(k-1))-alog(p(k))) )

;if alat(j) gt 0. then oplot,[0.,0.],[stratplev(i,j),stratplev(i,j)],psym=8,symsize=2,color=0
;wait,1.
                        goto,jumpz
                     endif
                 endfor
jumpz:
             endfor
          endif
      endfor
if max(stratplev) gt p(10) then stop
;
; polar daily average stratopause altitude
;
      y2d=fltarr(nlg,nlat)
      for i=0,nlg-1 do y2d(i,*)=alat
      index=where(y2d gt 0. and stratplev gt 0.,npts)
      if index(0) ne -1L then begin
         nhstratpavg(iday-1-leapday,iyr-lstyr)=total(stratplev(index))/float(npts)
         nhstratpmin(iday-1-leapday,iyr-lstyr)=min(stratplev(index))
         nhstratpmax(iday-1-leapday,iyr-lstyr)=max(stratplev(index))
;set_viewport,xorig(1),xorig(1)+xlen,yorig(1),yorig(1)+ylen
;level=160.+10.*findgen(17)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*icolmax/nlvls
;contour,tpzm,alat,p,levels=level,c_color=col1,/fill,/noerase,/ylog,color=0,xrange=[-90.,90.],yrange=[100.,0.01],title=sdate
;oplot,y2d(index),stratplev(index),psym=1,color=0
;oplot,[70.,70.],[nhstratpavg(iday-1-leapday,iyr-lstyr),nhstratpavg(iday-1-leapday,iyr-lstyr)],psym=8,color=250,symsize=3
;oplot,[70.,70.],[nhstratpmin(iday-1-leapday,iyr-lstyr),nhstratpmax(iday-1-leapday,iyr-lstyr)],psym=0,color=250,thick=3
      endif
      index=where(y2d lt 0. and stratplev gt 0.,npts)
      if index(0) ne -1L then begin
         shstratpavg(iday-1-leapday,iyr-lstyr)=total(stratplev(index))/float(npts)
         shstratpmin(iday-1-leapday,iyr-lstyr)=min(stratplev(index))
         shstratpmax(iday-1-leapday,iyr-lstyr)=max(stratplev(index))
;oplot,y2d(index),stratplev(index),psym=1,color=0
;oplot,[-70.,-70.],[shstratpavg(iday-1-leapday,iyr-lstyr),shstratpavg(iday-1-leapday,iyr-lstyr)],psym=8,color=250,symsize=3
;oplot,[-70.,-70.],[shstratpmin(iday-1-leapday,iyr-lstyr),shstratpmax(iday-1-leapday,iyr-lstyr)],psym=0,color=250,thick=3
      endif
;   imin=min(level)
;   imax=max(level)
;   ymnb=yorig(1)-cbaryoff
;   ymxb=ymnb+cbarydel
;   set_viewport,xorig(0),xorig(0)+xlen,ymnb,ymxb
;   !type=2^2+2^3+2^6
;   plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,xtitle='GEOS-4 Temperature'
;   ybox=[0,10,10,0,0]
;   x2=imin
;   dx=(imax-imin)/(float(nlvls)-1)
;   for j=1,nlvls-1 do begin
;       xbox=[x2,x2,x2+dx,x2+dx,x2]
;       polyfill,xbox,ybox,color=col1(j)
;       x2=x2+dx
;   endfor

      sdate_save(iday-1-leapday,iyr-lstyr)=string(FORMAT='(i4,i2.2,i2.2)',iyr,imn,idy)
      print,iyr,imn,idy,shstratpavg(iday-1-leapday,iyr-lstyr),nhstratpavg(iday-1-leapday,iyr-lstyr)
goto, jump

plotit:
restore,'timeseries_polar_stratopause_pressure_geos4.sav'	; nday,nyear,sdate_save,nhstratp,shstratp

; Autoscale if scale values for parameter/level are not defined
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,findgen(nday),findgen(10),/noeras,/nodata,title='Daily Mean Arctic Stratopause Pressure',$
     ytitle='Pressure (hPa)',xtitle='Day-of-Year',yrange=[20.,0.01],/ylog,color=0,charsize=1.5
nlvls=nyear
col1=1+indgen(nlvls)*icolmax/nlvls
yinc=(ymx-ymn)/nlvls
for n=0L,nyear-1 do begin
    nhsp=reform(nhstratpavg(*,n))
    nhspmax=reform(nhstratpmax(*,n))
    nhspmin=reform(nhstratpmin(*,n))
    good=where(nhsp ne -99.)
    oplot,good,nhsp(good),psym=8,color=col1(n)
    nhavg=nhsp(good)
    nhmax=nhspmax(good)
    nhmin=nhspmin(good)
    for nn=0,n_elements(good)-1L do begin
        plots,good(nn),nhmin(nn)
        plots,good(nn),nhmax(nn),/continue,color=col1(n)
    endfor
    bad=where(nhsp eq -99.)
    if bad(0) ne -1L then oplot,bad,0.01+0.*bad,psym=8,color=0
    xyouts,xmx+0.02,ymn+n*yinc,syear(n),charsize=1.25,color=col1(n),/normal,charthick=2
endfor

xyouts,.375,.95,'GEOS-4 Analyses',/normal,charsize=2,color=0
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,findgen(nday),findgen(10),/noeras,/nodata,title='Daily Mean Antarctic Stratopause Pressure',$
     ytitle='Pressure (hPa)',xtitle='Day-of-Year',yrange=[20.,0.01],/ylog,color=0,charsize=1.5
nlvls=nyear
col1=1+indgen(nlvls)*icolmax/nlvls
yinc=(ymx-ymn)/nlvls
for n=0L,nyear-1 do begin
    shsp=reform(shstratpavg(*,n))
    shspmax=reform(shstratpmax(*,n))
    shspmin=reform(shstratpmin(*,n))
    good=where(shsp ne -99.)
    oplot,good,shsp(good),psym=8,color=col1(n)
    shavg=shsp(good)
    shmax=shspmax(good)
    shmin=shspmin(good)
    for nn=0,n_elements(good)-1L do begin
        plots,good(nn),shmin(nn)
        plots,good(nn),shmax(nn),/continue,color=col1(n)
    endfor
    bad=where(nhsp eq -99.)
    if bad(0) ne -1L then oplot,bad,0.01+0.*bad,psym=8,color=0
    xyouts,xmx+0.02,ymn+n*yinc,syear(n),charsize=1.25,color=col1(n),/normal,charthick=2
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_polar_stratopause_pressure_geos4.ps -rotate -90 '+$
         'timeseries_polar_stratopause_pressure_geos4.jpg'
endif
end
