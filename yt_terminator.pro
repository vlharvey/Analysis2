;
; Plot SZA at noon at the GM for each latitude each day over 1 year
;
@stddat
@kgmt
@ckday
@kdate

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
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
lstmn=1
lstdy=1
lstyr=2007
ledmn=12
leddy=31
ledyr=2007
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nr=91
latitude=-90.+2.*findgen(nr)
sza_yt=fltarr(kday,nr)
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
;
; compute solar zenith angle at each latitude at the GM at noon
;
      gmt=12.
      rlon=0.
      doy=iday
      pi=3.14159265
      dtor=pi/180.
      earinc=23.5
      for ii=0L,nr-1 do begin
          rlat=latitude(ii)
          sinlat=sin(rlat*dtor)
          coslat=sqrt(1.-sinlat^2.)
          sinlon=sin(rlon*dtor)
          coslon=cos(rlon*dtor)
          soya=(doy-81.25)*pi/182.5           ; day angle
          soha=2.*pi*(gmt-12.)/24.            ; hour angle
          soha=-soha
          sininc=sin(earinc*dtor)
          sindec=sininc*sin(soya)
          cosdec= sqrt(1.-sindec^2.)
          coszen=cos(soha)*coslon+sin(soha)*sinlon
          coszen=coszen*cosdec*coslat
          coszen=sindec*sinlat+coszen
          coszen=min([max([coszen,-1.]),1.])
          chi = acos(coszen)
          sza_yt(icount,ii) = chi/dtor
      endfor
      sdate_all(icount)=sdate
      icount=icount+1L
goto,jump

plotit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
;xlabs=smon(xindex)
xlabs=smonth
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yt_terminator.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
level=5.*findgen(24)
nlvls=n_elements(level)
col1=reverse(1+indgen(nlvls)*icolmax/nlvls)
contour,sza_yt,1.+findgen(kday),latitude,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',title='Latitude of the Terminator',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.,yticks=6,/nodata
;contour,sza_yt,1.+findgen(kday),latitude,levels=level,color=0,/follow,/overplot,c_labels=1+fltarr(nlvls)
contour,sza_yt,1.+findgen(kday),latitude,levels=[90.],thick=8,color=0,/follow,/overplot,c_labels=[0]
;imin=min(level)
;imax=max(level)
;ymnb=yorig(0) -cbaryoff
;ymxb=ymnb  +cbarydel
;set_viewport,xmn,xmx,ymnb,ymxb
;!type=2^2+2^3+2^6
;plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='SZA',charsize=1.5
;ybox=[0,10,10,0,0]
;x1=imin
;dx=(imax-imin)/float(nlvls)
;for jj=0,nlvls-1 do begin
;    xbox=[x1,x1,x1+dx,x1+dx,x1]
;    polyfill,xbox,ybox,color=col1(jj)
;    x1=x1+dx
;endfor
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yt_terminator.ps -rotate -90 yt_terminator.jpg'
endif
end
