;
; GEOS-5 version
;
; plot daily Ubar in altitude and time
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_dat_origp

loadct,39
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
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nlg=0l
nlat=0l
nlv=0l
lstmn=1
lstdy=1
lstyr=2008
ledmn=3
leddy=31
ledyr=2008
lstday=0
ledday=0
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
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'+sdate+'_1200.V01.dat'
      rd_geos5_dat_origp,dir+ifile,iflg,nlg,nlat,nlv,alon,alat,height,p3d,g3d,t3d,u3d,v3d,q3d,pv3d
      if iflg ne 0 then goto, jump
      zp=g3d/1000.
      altitude=height
;
; first day
;
      if icount eq 0L then begin
         sdate0=sdate
         utz=fltarr(kday,nlv,nlat)
         sdates=strarr(kday)
;        rlat=0.
;        print,alat
;        read,' Enter Desired Latitude ',rlat
;        slat=string(format='(f5.1)',rlat)
      endif
      sdates(icount)=sdate
      print,sdate
;
; zonal mean zonal wind
;
      for k=0L,nlv-1L do begin
      for j=0L,nlat-1L do begin
          utz(icount,k,j)=total(u3d(*,j,k))/float(nlg)
      endfor
      endfor

icount=icount+1L
goto,jump

saveit:
save,file='zt_geos5_ubarp.sav',utz,alat,nlat,kday,sdates,altitude
plotit:
restore,'zt_geos5_ubarp.sav'

for j=0L,nlat-1L do begin
slat=strcompress(string(format='(f5.1)',alat(j)),/remove_all)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_geos5_ubarp_'+slat+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
;
; date labels
;
syear=strmid(sdates,0,4)
smon=strmid(sdates,4,2)
sday=strmid(sdates,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)
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
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
level=-100.+10.*findgen(nlvls)
utzlat=reform(utz(*,*,j))
print,slat,' ',min(utzlat),max(utzlat)
contour,utzlat,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,max(altitude)],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='GEOS-5 Ubar at '+slat+' Latitude',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
index=where(level lt 0)
contour,utzlat,1.+findgen(kday),altitude,levels=level(index),color=0,/follow,/overplot,c_labels=1+0*index,c_linestyle=5
index=where(level gt 0)
contour,utzlat,1.+findgen(kday),altitude,levels=level(index),color=mcolor,/follow,/overplot,c_labels=1+0*index
contour,utzlat,1.+findgen(kday),altitude,levels=[0],color=0,/follow,/overplot,c_labels=[0],thick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
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
   spawn,'convert -trim zt_geos5_ubarp_'+slat+'.ps -rotate -90 zt_geos5_ubarp_'+slat+'.jpg'
endif

endfor	; loop over latitude
end
