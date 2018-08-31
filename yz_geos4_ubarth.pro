;
; GEOS-4 version
;
; plot daily Ubar from theta data
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3

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
dir='/Users/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
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
      if ndays gt ledday then stop,' normal termination condition '
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-4 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      rd_geos5_nc3_meto,dir+sdate+'_1200.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,markold2,sf2,vp2,iflag
      if iflag ne 0 then goto, jump
;
; compute temperature
;
      t2=0.*u2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; zonal mean zonal wind
;
      tpzm=fltarr(nr,nth)
      uzm=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          tpzm(j,k)=total(t2(j,*,k))/float(nc)
          uzm(j,k)=total(u2(j,*,k))/float(nc)
      endfor
      endfor

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_geos4_ubarth_'+sdate+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
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
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=-120.+10.*findgen(nlvls)
contour,uzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=1.5,color=0,$
      ytitle='Pressure (hPa)`',title='GEOS-4 Zonal Mean Zonal Wind '+sdate,xticks=6,/fill,c_color=col1,$
      levels=level
index=where(level gt 0.)
contour,uzm,alat,th,levels=level(index),color=mcolor,/follow,/overplot
index=where(level lt 0.)
contour,uzm,alat,th,levels=level(index),color=0,/follow,/overplot
contour,uzm,alat,th,levels=[0],color=0,/follow,/overplot,thick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K/day)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot ne 'ps' then wait,.5
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_geos4_ubarth_'+sdate+'.ps -rotate -90 yz_geos4_ubarth_'+sdate+'.jpg'
endif

goto,jump
end
