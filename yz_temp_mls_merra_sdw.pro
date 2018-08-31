;
; plot zonal mean daily temperature from MLS, MERRA, and SD-WACCM
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
xorig=[0.1,0.4,0.7]
yorig=[0.3,0.3,0.3]
cbaryoff=0.1
cbarydel=0.01
xlen=0.2
ylen=0.4
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
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
dirm='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
dirw='/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h5.'
dir='/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_'

lstmn=12
lstdy=1
lstyr=2007
ledmn=2
leddy=28
ledyr=2008
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
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
kcount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then begin
         save,filename='tbar_mls_merra_sdw_'+sdate_all(0)+'.sav',mlstt,merratt,sdwtt,sdate_all,mlat,latitude_waccm,lat,pressure,lev,pmls2
         stop,' Normal Termination Condition '
      endif
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; read gridded MLS data on pressure
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
        ifile=dir+sdate+'.sav'
        dum=findfile(ifile)
        if dum(0) ne '' then restore,ifile
        if dum(0) eq '' then TP_GRID=0.*fltarr(144,96,55)
        mlstbar=mean(TP_GRID,dim=1)
        mlat=lat
;
; read MERRA
; LATITUDE_WACCM  FLOAT     = Array[96]
; LONGITUDE_WACCM FLOAT     = Array[144]
; PRESSURE        FLOAT     = Array[41]
; PSGRD           FLOAT     = Array[144, 96]
; QVGRD           FLOAT     = Array[144, 96, 41]
; TGRD            FLOAT     = Array[144, 96, 41]
; UGRD            FLOAT     = Array[144, 96, 41]
; VGRD            FLOAT     = Array[144, 96, 41]
; ZGRD            FLOAT     = Array[144, 96, 41]
;
        ifile=dirm+sdate+'.sav'
        restore,ifile
        merratbar=mean(tgrd,dim=1)
;
; read SD-WACCM
; H2O             FLOAT     = Array[144, 96, 88]
; LAT             DOUBLE    = Array[96]
; LEV             DOUBLE    = Array[88]
; LON             DOUBLE    = Array[144]
; N2O             FLOAT     = Array[144, 96, 88]
; O3              FLOAT     = Array[144, 96, 88]
; P               FLOAT     = Array[144, 96, 88]
; QSUM            FLOAT     = Array[144, 96, 88]
; T               FLOAT     = Array[144, 96, 88]
; U               FLOAT     = Array[144, 96, 88]
; V               FLOAT     = Array[144, 96, 88]
; Z               FLOAT     = Array[144, 96, 88]
; 
        ifile=dirw+sdate+'.sav'
        restore,ifile
        sdwtbar=mean(t,dim=1)
;
; retain temps for temperature tendency calculation
;
        if icount eq 0L then begin
           mlstt=fltarr(kday,n_elements(mlat),n_elements(pmls2))
           merratt=fltarr(kday,n_elements(latitude_waccm),n_elements(pressure))
           sdwtt=fltarr(kday,n_elements(lat),n_elements(lev))
           sdate_all=strarr(kday)
        endif
        mlstt(icount,*,*)=mlstbar
        merratt(icount,*,*)=merratbar
        sdwtt(icount,*,*)=sdwtbar
        sdate_all(icount)=sdate
;
; postscript file
;
        if setplot eq 'ps' then begin
           lc=0
           xsize=nxdim/100.
           ysize=nydim/100.
           set_plot,'ps'
           device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                  /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_temp_mls_merra_sdw_'+sdate+'.ps'
           !p.charsize=1.25
           !p.thick=2
           !p.charthick=5
           !y.thick=2
           !x.thick=2
        endif
;
; plot
;
       imin=110.
       level=imin+10.*findgen(21)
       imax=max(level)
       nlvls=n_elements(level)
       col1=1+(indgen(nlvls)/float(nlvls))*mcolor

       erase
       xyouts,.4,.9,sdate,color=0,charsize=2,charthick=2,/normal
       xmn=xorig(0)
       xmx=xorig(0)+xlen
       ymn=yorig(0)
       ymx=yorig(0)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,mlstbar,mlat,pmls2,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='MLS'
       contour,mlstbar,mlat,pmls2,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow

       xmn=xorig(1)
       xmx=xorig(1)+xlen
       ymn=yorig(1)
       ymx=yorig(1)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,merratbar,latitude_waccm,pressure,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='MERRA'
       contour,merratbar,latitude_waccm,pressure,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow

       xmn=xorig(2)
       xmx=xorig(2)+xlen
       ymn=yorig(2)
       ymx=yorig(2)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,sdwtbar,lat,lev,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='SD-WACCM'
       contour,sdwtbar,lat,lev,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow
       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xorig(0),xorig(2)+xlen,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='Zonal Mean Temperature (K)
       ybox=[0,10,10,0,0]
       x2=imin
       dx=(imax-imin)/(float(nlvls)-1)
       for j=1,nlvls-1 do begin
           xbox=[x2,x2,x2+dx,x2+dx,x2]
           polyfill,xbox,ybox,color=col1(j)
           x2=x2+dx
       endfor
;
; Close PostScript file and return control to X-windows
;      if setplot ne 'ps' then stop	;wait,1
       if setplot eq 'ps' then begin
          device, /close
          spawn,'convert -trim yz_temp_mls_merra_sdw_'+sdate+'.ps -rotate -90 '+$
                              'yz_temp_mls_merra_sdw_'+sdate+'.jpg'
          spawn,'rm -f yz_temp_mls_merra_sdw_'+sdate+'.ps'
       endif

skip:
      icount=icount+1L
goto,jump

end
