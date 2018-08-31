;
; plot zonal mean diabatic heating rate
;
@stddat
@kgmt
@ckday
@kdate

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=1L & lstdy=1L & lstyr=6L 
ledmn=12L & leddy=31L & ledyr=6L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
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
yorig=[0.15]
xlen=0.8
ylen=0.8
cbaryoff=0.08
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
dir='/aura7/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
ledyr=lstyr
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
      if ndays gt ledday then stop

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
      print,date

      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      rd_geos5_nc3_meto,dir+date+'_1200.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,markold2,sf2,vp2,iflag
      if iflag eq 1 then goto,skipit

      dum1=findfile(dir+date+'_1200.V01.nc4')
      if dum1(0) ne '' then ncid=ncdf_open(dir+date+'_1200.V01.nc4')
      if dum1(0) eq '' then goto,skipit
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      mark2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,mark2
      ncdf_close,ncid

      if kcount eq 0L then begin
         qyz=fltarr(nr,nth)
         sfile=strarr(kday)
         dum=transpose(mark2(*,*,0))
         lon=0.*dum
         lat=0.*dum
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy    ; area of each grid point
         endfor
         kcount=1L
      endif
      sfile(icount)=lfile
;
; loop over theta
;
      for thlev=0,nth-1 do begin
      for j=0,nr-1 do begin
          q1=q2(j,*,thlev)
          index=where(q1 ne 0.,nn)
          if index(0) ne -1 then qyz(j,thlev)=total(q1(index))/float(nn)
      endfor
      endfor
;
; plot 
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_geos4_q_'+date+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[min(alat),max(alat),max(alat),min(alat),min(alat)],[500.,500.,5000.,5000.,500.],min_value=0.,$
      xrange=[-90.,90.],yrange=[500.,5000.],/nodata,charsize=1.5,color=0,$
      ytitle='Theta (K)',title='GEOS-4 Diabatic Heating Rate '+date,xticks=6
nlvls=41
col1=1+indgen(nlvls)*icolmax/nlvls
level=-4000.+200.*findgen(nlvls)
contour,qyz,alat,th,levels=level,/fill,/cell_fill,/overplot,c_color=col1
index=where(level gt 0.)
contour,qyz,alat,th,levels=level(index),c_color=mcolor,/follow,/overplot
index=where(level lt 0.)
contour,qyz,alat,th,levels=level(index),c_color=0,/follow,/overplot
contour,qyz,alat,th,levels=[0],c_color=0,/follow,/overplot,thick=3
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

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_geos4_q_'+date+'.ps -rotate -90 yz_geos4_q_'+date+'.jpg'
endif

skipit:
icount=icount+1L
goto,jump

end
