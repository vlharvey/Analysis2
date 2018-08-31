;
; zonal mean temperature and zonal wind
; SD-WACCM
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=3L & lstdy=1L & lstyr=2009L 
ledmn=10L & leddy=1L & ledyr=2009L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)

lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
dir= '/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.zm.'
dir2= '/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.zm.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=long(ledday-lstday+1L)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition'
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; restore SDW on this day
; H2OBAR          FLOAT     = Array[96, 88]
; ILEV            DOUBLE    = Array[89]
; LAT             DOUBLE    = Array[96]
; LEV             DOUBLE    = Array[88]
; TBAR            FLOAT     = Array[96, 88]
; UBAR            FLOAT     = Array[96, 88]
; VBAR            FLOAT     = Array[96, 88]
; VSTAR           FLOAT     = Array[96, 89]
; WSTAR           FLOAT     = Array[96, 89]
; ZBAR            FLOAT     = Array[96, 88]
;
      if sdate lt '20120301' then dum=findfile(dir+sdate+'.sav')
      if sdate ge '20120301' then dum=findfile(dir2+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,dum
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='yz_sdwaccm_tbar_ubar_'+sdate+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot zt vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=23
col1=1+indgen(nlvls)*icolmax/nlvls
level=100.+10.*findgen(nlvls)
contour,tbar,lat,zbar/1000.,color=0,xtitle='Latitude',thick=6,yrange=[10,100],/noeras,ytitle='Approximate Altitude (km)',/fill,$
     title=sdate,levels=level,c_color=col1,charthick=2,charsize=1.5,xticks=6
contour,ubar,lat,zbar/1000.,levels=10+10*findgen(20),c_color=0,/follow,/overplot,c_labels=1+0*level,thick=3
contour,ubar,lat,zbar/1000.,levels=-200+10*findgen(20),c_color=mcolor,/follow,/overplot,c_labels=1+0*level,thick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='SD-WACCM Temperature (K)',charthick=2,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_sdwaccm_tbar_ubar_'+sdate+'.ps -rotate -90 yz_sdwaccm_tbar_ubar_'+sdate+'.jpg'
endif

      jumpday:
      icount=icount+1L
goto,jump

end
