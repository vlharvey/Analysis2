;
; plot altitude-time of MERRA2 Ubar
; user specified latitude and date range
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
cbaryoff=0.07
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
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif

goto,quick
;
; get file listing
;
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'

lstmn=1
lstdy=1
lstyr=2007
ledmn=1
leddy=1
ledyr=2017
lstday=0
ledday=0
if lstyr lt 1979 then stop,'Year out of range '
if ledyr lt 1979 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
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
      dum=findfile(dir+sdate+'00.nc3')
      if dum ne '' then ncfile0=dir+sdate+'00.nc3'
      rd_merra2_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag
      if iflag ne 0L then goto,jumpstep
      tmp2=0.*p2
      for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
;
; monthly files
;
      if icount eq 0L then begin
         sdate0=sdate
         rlat=20.
;        print,alat
;        read,'Enter desired latitude ',rlat
         index=where(abs(rlat-alat) eq min(abs(rlat-alat)))
         ilat=index(0)
         slat=strcompress(long(alat(ilat)),/remove_all)

         zbar=fltarr(kday,nth)
         ubar=fltarr(kday,nth)
      endif
      ubar(icount,*)=reform(mean(u2(ilat,*,*),dim=2))
      zbar(icount,*)=reform(mean(z2(ilat,*,*),dim=2))
;
icount=icount+1L
jumpstep:
goto,jump

plotit:
sdate1=sdate

quick:
daterange='20070101-20170101'
;restore,'merra2_ubar_20070101-20170101_0.sav
restore,'merra2_ubar_20070101-20170101_10.sav

;daterange=sdate0+'-'+sdate1
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='zt_merra2_ubar_'+daterange+'_'+slat+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
imin=-50.
int=10
nlvls=11
col1=1+indgen(nlvls)*mcolor/nlvls
level=imin+int*findgen(nlvls)
imax=max(level)
index=where(ubar eq 0. or zbar eq 0.)
if index(0) ne -1L then ubar(index)=0./0.
if index(0) ne -1L then zbar(index)=0./0.
;ubarsm=0*ubar
;result=size(ubar)
;nl=result(2)
;for k=0L,nl-1L do begin
;ubarsm(*,k)=smooth(reform(ubar(*,k)),31,/Nan)
;endfor
;ubar=ubarsm
xindex=where(strmid(sdate_all,4,4) eq '0701',nx)
xlabs=strmid(sdate_all(xindex),0,4)
zprof=mean(zbar,dim=1)
contour,ubar,findgen(icount),zprof,levels=level,/cell_fill,c_color=col1,color=0,xrange=[0,icount-1],$
            ytitle='Altitude (km)',yrange=[15.,75.],charsize=1.5,charthick=2,xticks=nx-1,xtickname=xlabs,xtickv=xindex
;contour,ubar,findgen(icount),zbar,levels=[20],/overplot,/follow,color=0,thick=2
;contour,ubar,findgen(icount),zbar,levels=[-20],/overplot,/follow,color=mcolor,c_linestyle=5,thick=2
;contour,ubar,findgen(icount),zbar,levels=0,/overplot,/follow,color=0,thick=1

index=where(sdate_all eq '20140701')
oplot,index(0)+0*findgen(60),15+findgen(60),thick=3,color=250
oplot,index(0)+0*findgen(60),15+findgen(60),thick=3,color=255,linestyle=5
index=where(sdate_all eq '20090701')
oplot,index(0)+0*findgen(60),15+findgen(60),thick=3,color=250
oplot,index(0)+0*findgen(60),15+findgen(60),thick=3,color=255,linestyle=5

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MERRA2 Ubar at '+slat+' N (m/s)',/noeras,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;save,file='merra2_ubar_'+daterange+'_'+slat+'.sav',ubar,zbar,icount,sdate_all,slat

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim zt_merra2_ubar_'+daterange+'_'+slat+'.ps -rotate -90 '+$
                            'zt_merra2_ubar_'+daterange+'_'+slat+'.jpg'
     endif
end
