;
; read MERRA2 IDL save data on pressure. Interpolate Q to MERRA pressure levels.
; average same UT for multiple days and save out "composite day" (averages at 00,06,12,18)
; Input data:  IDL> restore,'/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_press_YYYYMMDD.sav
; Input data:  IDL> restore,'/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_rad_YYYYMMDD.sav
; 
; LATITUDE_WACCM  FLOAT     = Array[96]
; LONGITUDE_WACCM FLOAT     = Array[144]
; PRESSURE        FLOAT     = Array[41]
; PSGRD           FLOAT     = Array[144, 96]
; O3GRD           FLOAT     = Array[144, 96, 41]
; QVGRD           FLOAT     = Array[144, 96, 41]
; TGRD            FLOAT     = Array[144, 96, 41]
; UGRD            FLOAT     = Array[144, 96, 41]
; VGRD            FLOAT     = Array[144, 96, 41]
; ZGRD            FLOAT     = Array[144, 96, 41]
; add QGRD            FLOAT     = Array[144, 96, 41]

; LATITUDE_WACCM  FLOAT     = Array[96]
; LEVELS_RAD      DOUBLE    = Array[42]
; LONGITUDE_WACCM FLOAT     = Array[144]
; QNEW            FLOAT     = Array[144, 96, 42]

loadct,38
device,decompose=0
mcolor=byte(!p.color)
nlvls=30L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
dirw='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_'
stime=['00','06','12','18']
ntime=n_elements(stime)
for itime=0L,ntime-1L do begin
;
; test on first 9 days in January 2015
;
ifiles=file_search(dirw+'press_2015010?'+stime(itime)+'.sav')
nfile=n_elements(ifiles)
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN

result=strsplit(ifiles(n),'_',/extract)
result2=strsplit(result(3),'.',/extract)
sdate=result2(0)
print,sdate
;dum=findfile(dirw+'theta_'+sdate+'.nc3')
;if dum(0) ne '' then goto,jumpfile
restore,dirw+'press_'+sdate+'.sav'
print,'reading '+dirw+'press_'+sdate+'.sav'
restore,dirw+'rad_'+sdate+'.sav'
alon=LONGITUDE_WACCM
alat=LATITUDE_WACCM
nc=n_elements(alon)
nr=n_elements(alat)
nl=n_elements(pressure)
;
; remove tropospheric levels for diabatic heating rate for which there is no Nv data
;
index=where(levels_rad le max(pressure))
levels_rad=levels_rad(index)
qnew=reform(qnew(*,*,index))
;
; interpoloate qnew to qgrd with same pressure levels as Nv data (fill above 0.1 with NaN)
;
qgrd=0.*tgrd
for k=0L,nl-1L do begin
    zp=pressure(k)
    if zp lt min(levels_rad) then qgrd(*,*,k)=0./0.
    if zp lt min(levels_rad) then goto,jumplev
    for kk=0L,n_elements(levels_rad)-2L do begin
        if levels_rad(kk) le zp and levels_rad(kk+1L) ge zp then begin
           zscale=(alog(zp)-alog(levels_rad(kk)))/(alog(levels_rad(kk+1L))-alog(levels_rad(kk)))
           qgrd(*,*,k)=qnew(*,*,kk)-zscale*(qnew(*,*,kk)-qnew(*,*,kk+1))
        endif
    endfor
    jumplev:
endfor
;
; average multiple days together at this time
;
if itime eq 0L and n eq 0L then begin
   o3avg=fltarr(ntime,nc,nr,nl)
   qavg=fltarr(ntime,nc,nr,nl)
   qvavg=fltarr(ntime,nc,nr,nl)
   tavg=fltarr(ntime,nc,nr,nl)
   uavg=fltarr(ntime,nc,nr,nl)
   vavg=fltarr(ntime,nc,nr,nl)
   zavg=fltarr(ntime,nc,nr,nl)
endif
o3avg(itime,*,*,*)=o3avg(itime,*,*,*)+O3GRD
qavg(itime,*,*,*)=qavg(itime,*,*,*)+QGRD
qvavg(itime,*,*,*)=qvavg(itime,*,*,*)+QVGRD
tavg(itime,*,*,*)=tavg(itime,*,*,*)+TGRD
uavg(itime,*,*,*)=uavg(itime,*,*,*)+UGRD
vavg(itime,*,*,*)=vavg(itime,*,*,*)+VGRD
zavg(itime,*,*,*)=zavg(itime,*,*,*)+ZGRD

jumpfile:
ENDFOR		; LOOP OVER files
;
;*** check***
;
;for kk=0L,nl-1L do begin
kk=10           ; 0.2 hPa
rlev=pressure(kk)
;print,pressure
;;read,'Enter theta surface ',rlev
index=where(pressure eq rlev)
ilev=index(0)
slev=string(rlev)
pp=reform(qavg(itime,*,*,ilev))/float(nfile)
pp2=reform(zavg(itime,*,*,ilev))/float(nfile)
level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
!type=2^2+2^3
erase
contour,pp,alon,alat,levels=level,/cell_fill,c_color=col1,/noeras,xrange=[0.,360.],yrange=[-90.,90.],title=stime(itime)+'  '+slev+' K'
index=where(level lt 0.)
if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=mcolor,c_linestyle=5,/noeras,/overplot
index=where(level gt 0.)
if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=0,/noeras,/overplot
level2=min(pp2)+((max(pp2)-min(pp2))/float(nlvls))*findgen(nlvls)
contour,pp2,alon,alat,levels=level2,/follow,c_color=0,/noeras,/overplot,thick=3
;endfor
;
;*** end check***
;
ENDFOR		; LOOP OVER times/day
;
; average multiple days at each UT
;
o3avg=o3avg/float(nfile)
qavg=qavg/float(nfile)
qvavg=qvavg/float(nfile)
tavg=tavg/float(nfile)
uavg=uavg/float(nfile)
vavg=vavg/float(nfile)
zavg=zavg/float(nfile)
;
; calculate LT from alon and stime
;
ltime=fltarr(ntime,nc)
for itime=0L,ntime-1L do begin
    ltime(itime,*)=stime(itime)+alon*24./360.
endfor
index=where(ltime gt 24.)
ltime(index)=ltime(index)-24.
;
; write composite-day
;
ofile=dirw+'composite-day.sav'
print,'writing ',ofile
save,file=ofile,nr,nc,nl,ntime,alon,alat,pressure,stime,ltime,uavg,vavg,zavg,tavg,qvavg,qavg,o3avg

end
