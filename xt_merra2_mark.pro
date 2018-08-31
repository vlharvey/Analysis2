;
; Hovmoller of MERRA-2 marker - allow user to choose individual latitudes
; apply to Asian monsoon
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
xorig=[0.3]
yorig=[0.2]
cbaryoff=0.1
cbarydel=0.01
xlen=0.4
ylen=0.6
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

;goto,skipit
;
; user enters dynamic date range
;
for iyear=1980, 2016 do begin
lstmn=4
lstdy=1
lstyr=iyear
ledmn=10
leddy=31
ledyr=iyear
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
         rth=350
;        print,th
;        read,'Enter desired theta surface ',rth
         index=where(abs(rth-th) eq min(abs(rth-th)))
         ith=index(0)
         sth=strcompress(long(th(ith)),/remove_all)

         rlat=35.
;        print,alat
;        read,'Enter desired latitude ',rlat
         index=where(abs(rlat-alat) eq min(abs(rlat-alat)))
         ilat=index(0)
         slat=strcompress(alat(ilat),/remove_all)

         pv_xt=fltarr(nc,kday)
         z_xt=fltarr(nc,kday)
         mark_xt=fltarr(nc,kday)
         u_xt=fltarr(nc,kday)
         tmp_xt=fltarr(nc,kday)
      endif
;
; save longitude strip at ilat on each day (when rlat = 35 the lat band is 25-40)
;
      pv_xt(*,icount)=mean(pv2(ilat-5:ilat+3,*,ith),dim=1)
      z_xt(*,icount)=mean(z2(ilat-5:ilat+3,*,ith),dim=1)
      mark_xt(*,icount)=mean(mark2(ilat-5:ilat+3,*,ith),dim=1)
      u_xt(*,icount)=mean(u2(ilat-5:ilat+3,*,ith),dim=1)
      tmp_xt(*,icount)=mean(tmp2(ilat-5:ilat+3,*,ith),dim=1)

      icount=icount+1L
goto,jump

plotit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='xt_merra2_'+sdate0+'-'+sdate1+'_'+sth+'_'+slat+'.sav',pv_xt,z_xt,mark_xt,u_xt,tmp_xt,alon,sdate_all,kday,sdate0,sdate1,slat,sth,nc,ilat,alat
skipit:
;restore,'xt_merra2_20140401-20141031_350_29.3686.sav
;restore,'xt_merra2_20150401-20151031_350_29.3686.sav

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
   device,/landscape,bits=8,filename='xt_merra2_mark_'+sdate0+'-'+sdate1+'_'+sth+'_'+slat+'.ps'
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
pv_xt=pv_xt*1.e6
imin=min(pv_xt)
imax=5.	;max(pv_xt)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
;pv_xt=smooth(pv_xt,3,/edge_truncate)
loadct,23
contour,pv_xt,alon,findgen(kday),charsize=1.5,/noeras,xrange=[0.,360.],$
        xtitle='Longitude',levels=level,c_color=col1,/cell_fill,color=0,title=sth+' K & '+slat+' N',$
        yrange=[kday-1,0],charthick=1.5,xticks=6,yticks=nxticks-1,ytickname=xlabs,ytickv=xindex
loadct,0
xyouts,250.,15,syr(0),charsize=3,charthick=4,color=mcolor,/data
loadct,23
;contour,pv_xt,alon,findgen(kday),charsize=1.5,/noeras,/overplot,levels=[0.],color=0,thick=5
;index=where(level gt 0.)
;contour,pv_xt,alon,findgen(kday),charsize=1.5,/noeras,/overplot,levels=level(index),color=0
;index=where(level lt 0.)
;contour,pv_xt,alon,findgen(kday),charsize=1.5,/noeras,/overplot,levels=level(index),c_linestyle=5,color=mcolor

;index=where(mark_xt eq 0.)
;if index(0) ne -1L then mark_xt(index)=0./0.
;mark_xt=smooth(mark_xt,3,/edge_truncate,/Nan)
;index=where(finite(mark_xt) eq 0)
;if index(0) ne -1L then mark_xt(index)=0.
contour,mark_xt,alon,findgen(kday),levels=[0.1],color=mcolor,/follow,/noeras,/overplot,thick=15,c_labels=[0]
contour,mark_xt,alon,findgen(kday),levels=[-0.1],color=mcolor*.95,/follow,/noeras,/overplot,thick=15

x2d=0.*pv_xt
y2d=0.*pv_xt
for i=0L,nc-1L do y2d(i,*)=findgen(kday)
for j=0L,kday-1L do x2d(*,j)=alon
index=where(mark_xt lt 0.)
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
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MERRA-2 PV (PVU)'
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
   spawn,'convert xt_merra2_mark_'+sdate0+'-'+sdate1+'_'+sth+'_'+slat+'.ps -rotate -90 xt_merra2_mark_'+sdate0+'-'+sdate1+'_'+sth+'_'+slat+'.jpg'
endif

endfor	; loop over years
end
