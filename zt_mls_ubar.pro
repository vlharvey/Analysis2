;
; ZT of MLS Ubar
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.20]
yorig=[0.3]
cbaryoff=0.1
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
syear=['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016']
syear=['2016']
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
;
; get file listing
;
mdir3='/atmos/aura6/data/MLS_data/Datfiles_Grid/'

for iyear=0L,nyear-1L do begin

lstmn=4
lstdy=1
lstyr=long(syear(iyear))
ledmn=5
leddy=22
ledyr=long(syear(iyear))
lstday=0
ledday=0
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
      sdate_all(icount)=sdate
      dayno(icount)=iday
      print,iday,' ',sdate
;
; restore MLS gridded GPH on this day
; CO_GRID         FLOAT     = Array[144, 96, 37]
; GP_GRID         FLOAT     = Array[144, 96, 55]
; H2O_GRID        FLOAT     = Array[144, 96, 55]
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; N2O_GRID        FLOAT     = Array[144, 96, 37]
; O3_GRID         FLOAT     = Array[144, 96, 55]
; PMLS            FLOAT     = Array[37]
; PMLS2           FLOAT     = Array[55]
; TP_GRID         FLOAT     = Array[144, 96, 55]
;
      dum=findfile(mdir3+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,jumpstep
      restore,dum(0)
      if icount eq 0L then begin
         pressure=PMLS2
         rlat=60.
;        print,lat
;        read,'Enter desired latitude ',rlat
         index=where(abs(lat-rlat) eq min(abs(lat-rlat)))
         ilat=index(0)
         slat=strcompress(long(lat(ilat)),/r)
         pressure_altitude=alog(1000./pressure)*7.
         nlv2=n_elements(pressure)

         zbar=fltarr(kday,nlv2)
         ubar=fltarr(kday,nlv2)
    endif
    umean=mean(u,dim=4)		; mean over both nodes
    zmean=mean(gph,dim=4)
    utmp=reform(umean(*,ilat,*))
    ztmp=reform(zmean(*,ilat,*))
    ubar(icount,*)=mean(utmp,dim=1)
    zbar(icount,*)=mean(ztmp,dim=1)
jumpstep:
icount=icount+1L
goto,jump

plotit:
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='zt_mls_ubar_'+syear(iyear)+'_'+slat+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=-35.
int=5
nlvls=13
col1=1+indgen(nlvls)*mcolor/nlvls
level=imin+int*findgen(nlvls)
imax=max(level)
index=where(ubar eq 0. or zbar eq 0.)
if index(0) ne -1L then ubar(index)=0./0.
if index(0) ne -1L then zbar(index)=0./0.
contour,ubar,dayno,pressure,/ylog,levels=level,/cell_fill,c_color=col1,color=0,title=syear(iyear),$
            ytitle='Pressure (hPa)',yrange=[100.,0.01],charsize=1.5,charthick=2,xtitle='DOY',xrange=[dayno(0),dayno(-1)]
contour,ubar,dayno,pressure,/ylog,levels=10.+10.*findgen(10),/overplot,/follow,color=0,thick=2
contour,ubar,dayno,pressure,/ylog,levels=-100.+10.*findgen(10),/overplot,/follow,color=mcolor,c_linestyle=5,thick=2
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Ubar at '+slat+' N (m/s)',/noeras,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;save,file='mls_ubar_'+syear(iyear)+'_'+slat+'.sav',ubar,dayno,pressure,zbar

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim zt_mls_ubar_'+syear(iyear)+'_'+slat+'.ps -rotate -90 '+$
                            'zt_mls_ubar_'+syear(iyear)+'_'+slat+'.jpg'
     endif

endfor	; loop over years
end
