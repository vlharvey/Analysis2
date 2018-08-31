;
; NOTE: ERA40 LONGITUDES ARE -180 TO 180
; ERA40 zonal mean temperature and zonal wind at all latitudes for 45 year record
;
; VLH 9/10/09
;
@stddat
@kgmt
@ckday
@kdate
@rd_era40_nc

; define viewport location
loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=750
nydim=750
xorig=[0.125]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
!p.background=mcolor
!NOERAS=-1
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

dir='/Volumes/earth/harvey/ERA40_data/Datfiles/era40_ua_12Z_'
;lstmn=9L & lstdy=1L & lstyr=1957L
;ledmn=8L & leddy=31L & ledyr=2002L

lstmn=2L & lstdy=23L & lstyr=1967L
ledmn=12L & leddy=27L & ledyr=1991L
;
; Get start and end dates
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1950 or lstyr gt 2002 then stop,'Year out of range '
if ledyr lt 1950 or ledyr gt 2002 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' Starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
;
; read ERA40 pressure data, i.e. /aura7/harvey/ERA40_data/Datfiles/era40_ua_12Z_19600131.nc
;
      rd_era40_nc,dir+date+'.nc',nc,nr,nl,alon,alat,press,tp,uu,vv,gp,iflg
      if iflg eq 1 then goto,jump
;
; calculate zonal mean temperature and zonal wind
;
      uzm=-9999.+0.*fltarr(nr,nl)
      tzm=-9999.+0.*fltarr(nr,nl)
      gzm=-9999.+0.*fltarr(nr,nl)
      for k=0,nl-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tp(*,j,k))/float(nc)
              uzm(j,k)=total(uu(*,j,k))/float(nc)
              gzm(j,k)=total(gp(*,j,k))/float(nc)
              if tzm(j,k) gt 400. or tzm(j,k) lt 100. then stop,'check temperature'
              if abs(uzm(j,k)) gt 200. then stop,'check zonal wind'
          endfor
      endfor
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_tbar_ubar_era40_'+date+'.ps'
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=27
col1=1+indgen(nlvls)*icolmax/nlvls
tlevel=170.+5.*findgen(nlvls)
contour,tzm,alat,press,/noeras,xrange=[-90.,90.],yrange=[100.,min(press)],/ylog,charsize=1.5,color=0,$
      ytitle='Pressure (hPa)',title=date,xticks=6,/cell_fill,c_color=col1,levels=tlevel,xtitle='Latitude'
index=where(tlevel mod 10. eq 0)
contour,tzm,alat,press,levels=tlevel(index),color=0,/follow,/overplot,c_labels=0*index,/noerase
contour,uzm,alat,press,levels=30.+10*findgen(20),color=0,/follow,/overplot,c_labels=1+0*findgen(20),/noerase,thick=4
contour,uzm,alat,press,levels=-200.+10*findgen(20),color=mcolor,/follow,/overplot,c_labels=1+0*findgen(20),/noerase,$
        thick=8,c_linestyle=5
;contour,markzm,alat,press,levels=[0.5],color=0,/follow,/overplot,thick=5,/noerase
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='ERA-40 Temperature (K)'
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
       device,/close
       spawn,'convert -trim yz_tbar_ubar_era40_'+date+'.ps -rotate -90 yz_tbar_ubar_era40_'+date+'.jpg'
       spawn,'rm -f yz_tbar_ubar_era40_'+date+'.ps'
    endif

goto,jump
end
