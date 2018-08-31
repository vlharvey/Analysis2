;
; latitude-time slice of PV - allow user to choose individual longitudes
; MERRA-2 version
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra2_nc3

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.20]
yorig=[0.3]
cbaryoff=0.05
cbarydel=0.01
xlen=0.6
ylen=0.4
device,decompose=0
!NOERAS=-1
nlvls=20L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016']
;dum=1979+indgen(36)
;dum=2008
;syear=strcompress(dum,/remove_all)
nyear=n_elements(syear)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
;
; user enters dynamic date range
;

lstmn=12
lstdy=15
lstyr=2014
ledmn=2
leddy=15
ledyr=2015
lstday=0
ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 79 then lstyr=lstyr+2000
if ledyr lt 79 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1979 then stop,'Year out of range '
if ledyr lt 1979 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
dayno=lonarr(kday)
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
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
      sdate_all(icount)=sdate
      dayno(icount)=iday
;
; read data
;
      dum=findfile(dir+sdate+'00.nc3')
      if dum ne '' then ncfile0=dir+sdate+'00.nc3'
      rd_merra2_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag
      if iflag ne 0L then goto,jump
      tmp2=0.*p2
      for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
;
; normalize marker
;
      index=where(mark2 lt 0.)
      if index(0) ne -1L then mark2(index)=-1.
;
; first day
;
      if icount eq 0L then begin
         rth=2000.
;        print,th
;        read,'Enter desired theta surface ',rth
         index=where(abs(rth-th) eq min(abs(rth-th)))
         ith=index(0)
         sth=strcompress(long(th(ith)),/remove_all)

         rlon=270.
;        print,alon
;        read,'Enter desired longitude ',rlon
         index=where(abs(rlon-alon) eq min(abs(rlon-alon)))
         ilon=index(0)
         slon=strcompress(alon(ilon),/remove_all)

         pv_yt=fltarr(kday,nr)
         z_yt=fltarr(kday,nr)
         mark_yt=fltarr(kday,nr)
         u_yt=fltarr(kday,nr)
         tmp_yt=fltarr(kday,nr)

         pvbar_yt=fltarr(kday,nr)
         zbar_yt=fltarr(kday,nr)
         markbar_yt=fltarr(kday,nr)
      endif
;
; save zonal means on each day
;
      pv_yt(icount,*)=reform(pv2(*,ilon,ith))
      pv_yt(icount,*)=min(pv2(*,ilon-2:ilon+2,ith),dim=2)

      z_yt(icount,*)=reform(z2(*,ilon,ith))
      mark_yt(icount,*)=reform(mark2(*,ilon,ith))
;     mark_yt(icount,*)=min(mark2(*,*,ith),dim=2)

      u_yt(icount,*)=reform(u2(*,ilon,ith))
      tmp_yt(icount,*)=reform(tmp2(*,ilon,ith))

      pvbar_yt(icount,*)=mean(reform(pv2(*,*,ith)),dim=2)
      zbar_yt(icount,*)=mean(reform(z2(*,*,ith)),dim=2)
      markbar_yt(icount,*)=mean(reform(mark2(*,*,ith)),dim=2)

      icount=icount+1L
goto,jump

plotit:

;sdate0=sdate_all(0)
;sdate1=sdate_all(n_elements(sdate_all)-1)
;save,file='yt_merra2_pv_'+sdate0+'-'+sdate1+'.sav',zbar_yt,tbar_yt,ubar_yt,$
;     vbar_yt,qbar_yt,pvbar_yt,markbar_yt,pbar_yt,sdate_all,alat,alat,th,kday
;restore,'yt_merra2_pv_20150101-20150201.sav

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yt_merra2_pv_'+sdate0+'-'+sdate1+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
syr=strmid(sdate_all,0,4)
smn=strmid(sdate_all,4,2)
sdy=strmid(sdate_all,6,2)
xindex=where(sdy eq '01' or sdy eq '15',nxticks)
xlabs=smn(xindex)+'/'+sdy(xindex)

nlvls=31
col1=1+indgen(nlvls)*mcolor/nlvls
imin=min(pv_yt)
imax=max(pv_yt)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
pv_yt=smooth(pv_yt,3,/edge_truncate)
pv_yt(*,0:1)=0./0.
pv_yt(*,-2:-1)=0./0.
contour,pv_yt,findgen(kday),alat,charsize=1.5,/noeras,yrange=[0.,90.],$
        ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=sth+' K & '+slon+' E',$
        xrange=[0.,kday-1],charthick=1.5,yticks=6,xticks=nxticks-1,xtickname=xlabs
contour,pv_yt,findgen(kday),alat,charsize=1.5,/noeras,/overplot,levels=[0.],color=0,thick=5
index=where(level gt 0.)
contour,pv_yt,findgen(kday),alat,charsize=1.5,/noeras,/overplot,levels=level(index),color=0
index=where(level lt 0.)
contour,pv_yt,findgen(kday),alat,charsize=1.5,/noeras,/overplot,levels=level(index),c_linestyle=5,color=mcolor

contour,smooth(mark_yt,3,/edge_truncate),findgen(kday),alat,levels=[0.1],color=mcolor,/follow,/noeras,/overplot,thick=15,c_labels=[0]
;contour,smooth(mark_yt,3,/edge_truncate),findgen(kday),alat,levels=[-0.1],color=mcolor*.9,/follow,/noeras,/overplot,thick=5

x2d=0.*pv_yt
for j=0L,nr-1L do x2d(*,j)=findgen(kday)
y2d=0.*pv_yt
for i=0L,kday-1L do y2d(i,*)=alat
index=where(y2d gt 10. and pv_yt lt 0.)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=0
;
yearlab=syr(0)
if syr(0) ne syr(-1) then yearlab=syr(0)+'-'+syr(-1)
xyouts,0.,-85.,yearlab,/data,color=mcolor,charsize=3,charthick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MERRA-2 PV'
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
   spawn,'convert yt_merra2_pv_'+sdate0+'-'+sdate1+'.ps -rotate -90 yt_merra2_pv_'+sdate0+'-'+sdate1+'.jpg'
endif
end
