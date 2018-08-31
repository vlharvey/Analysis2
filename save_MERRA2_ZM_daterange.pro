;---------------------------------------------------------------------------------------------------
; pressure data
; save zonal mean zonal data as a function of altitude and latitude based on MERRA2 data - take 12Z (for now) for comparison to MERRA. stop at Jan 31 2016 for comparison to MERRA
; save over specified dates
;
;	 -------------------------------
;       |         Lynn Harvey           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 8/2/2016      |
;	 -------------------------------
;
@stddat			; Determines the number of days since Jan 1, 1956
@kgmt			; This function computes the Julian day number (GMT) from the
@ckday			; This routine changes the Julian day from 365(6 if leap yr)
@kdate			; gives back kmn,kdy information from the Julian day #.
@rd_merra_nc3

;-----------------------------------------------------

nxdim=750
nydim=750
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=!p.color
icolmax=255
mcolor=icolmax
device,decompose=0
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=1980
ledmn=1
leddy=31
ledyr=2016
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date ',lstmn,lstdy,lstyr
;read,' Enter ending date ',ledmn,leddy,ledyr
if lstyr lt 79 then lstyr=lstyr+2000
if ledyr lt 79 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
minyear=strcompress(lstyr,/remove_all)+string(FORMAT='(I2.2)',lstmn)+string(FORMAT='(I2.2)',lstdy)
maxyear=strcompress(ledyr,/remove_all)+string(FORMAT='(I2.2)',ledmn)+string(FORMAT='(I2.2)',leddy)
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_press_'
dir2='/atmos/harvey/MERRA2_data/Datfiles_PW/MERRA2_PW1-2_'
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;itimes=['00','06','12','18']
itimes=['12']
ntimes=n_elements(itimes)
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
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; read data
; LATITUDE_WACCM  FLOAT     = Array[96]
; LONGITUDE_WACCM FLOAT     = Array[144]
; O3GRD           FLOAT     = Array[144, 96, 41]
; PRESSURE        FLOAT     = Array[41]
; PSGRD           FLOAT     = Array[144, 96]
; QVGRD           FLOAT     = Array[144, 96, 41]
; TGRD            FLOAT     = Array[144, 96, 41]
; UGRD            FLOAT     = Array[144, 96, 41]
; VGRD            FLOAT     = Array[144, 96, 41]
; ZGRD            FLOAT     = Array[144, 96, 41]
;
; COMMENT STRING = 'Wave-1 phase in degrees, represents inflection point 90 degrees east of the trough. Wave-2 phase in degrees, represents inflection point in Eastern Hemisphere, 45 degrees east of the trough.'
; LAT             FLOAT     = Array[96]
; MERRAWAVE1HEIGHT FLOAT     = Array[96, 41]
; MERRAWAVE1PHASE FLOAT     = Array[96, 41]
; MERRAWAVE2HEIGHT FLOAT     = Array[96, 41]
; MERRAWAVE2PHASE FLOAT     = Array[96, 41]
; PRESSURE        FLOAT     = Array[41]
;
ihour=0L
;     for ihour=0,ntimes-1L do begin
          shour=itimes(ihour)
          print,sdate+shour
          dum=findfile(dir+sdate+shour+'.sav')
          if dum eq '' then stop,'Pressure data missing on '+sdate+shour
;         dum2=findfile(dir2+sdate+shour+'.sav')
;         if dum2 eq '' then stop,'PW data missing on '+sdate+shour
          restore,dum
;         restore,dum2	; cannot locate code that computes PW diagnostics. resolve later
;
; declare arrays on first days
;
          if icount eq 0 then begin
             nr=n_elements(LATITUDE_WACCM)
             nl=n_elements(pressure)
             sdates=strarr(kday)
             zbar2d=fltarr(kday,nr,nl)
             ubar2d=fltarr(kday,nr,nl)
             tbar2d=fltarr(kday,nr,nl)
;             pw1amp=fltarr(kday,nr,nl)
;             pw1phase=fltarr(kday,nr,nl)
;             pw2amp=fltarr(kday,nr,nl)
;             pw2phase=fltarr(kday,nr,nl)
          endif
          sdates(icount)=sdate		;+shour
;
; save lat/alt of Tbar and Ubar each day
;
          ubar2d(icount,*,*)=mean(ugrd,dim=1)
          tbar2d(icount,*,*)=mean(tgrd,dim=1)
          zbar2d(icount,*,*)=mean(zgrd,dim=1)
;          pw1amp(icount,*,*)=MERRAWAVE1HEIGHT
;          pw1phase(icount,*,*)=MERRAWAVE1PHASE
;          pw2amp(icount,*,*)=MERRAWAVE2HEIGHT
;          pw2phase(icount,*,*)=MERRAWAVE2PHASE

;erase
;!type=2^2+2^3
;set_viewport,.1,.9,.1,.9
;nlvls=21
;level=-100+10*findgen(nlvls)
;col1=(findgen(nlvls)/float(nlvls))*mcolor
;index=where(zbar2d eq 0.)
;if index(0) ne -1L then zbar2d(index)=0./0.
;contour,reform(ubar2d(icount,*,*)),lat,zbar2d(icount,*,*),levels=level,c_color=col1,/cell_fill,/noeras,title=sdate+shour,ytitle='Altitude (km)',xtitle='Latitude',charsize=2,charthick=2
;index=where(level lt 0.)
;contour,ubar2d(icount,*,*),lat,zbar2d(icount,*,*),levels=level(index),/foll,color=mcolor,c_linestyle=5,/noeras,/overplot
;nlvls=25
;levels=15*findgen(25)
;col1=findgen(nlvls)*255/float(nlvls)
;contour,MERRAWAVE1HEIGHT,lat,zbar2d(icount,*,*),/overplot,thick=2,levels=100+100*findgen(30),/noerase
;contour,MERRAWAVE1PHASE,lat,zbar2d(icount,*,*),/overplot,levels=levels,c_color=col1,/follow,thick=3 
;stop
           icount=icount+1L
;      endfor

goto,jump

saveit:
lat=latitude_waccm
;
; save
;
save,filename='merra2_ZM-Z_'+yearlab+'.sav',sdates,pressure,lat,zbar2d
save,filename='merra2_ZM-U_'+yearlab+'.sav',sdates,pressure,lat,ubar2d
save,filename='merra2_ZM-T_'+yearlab+'.sav',sdates,pressure,lat,tbar2d
;save,filename='merra2_ZM-PW1A_'+yearlab+'.sav',sdates,pressure,lat,pw1amp
;save,filename='merra2_ZM-PW1P_'+yearlab+'.sav',sdates,pressure,lat,pw1phase
;save,filename='merra2_ZM-PW2A_'+yearlab+'.sav',sdates,pressure,lat,pw2amp
;save,filename='merra2_ZM-PW2P_'+yearlab+'.sav',sdates,pressure,lat,pw2phase
end
