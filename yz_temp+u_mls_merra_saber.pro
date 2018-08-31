;
; plot zonal mean daily temperature and zonal wind from MLS, MERRA, and SABER
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
xorig=[0.1,0.325,0.55,0.775]
yorig=[0.6,0.6,0.6,0.6]
cbaryoff=0.1
cbarydel=0.01
xlen=0.2
ylen=0.3
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
dir='/atmos/aura6/data/MLS_data/Datfiles_Grid/'
dirs='/atmos/harvey/SABER_data/Datfiles_Grid/'

lstmn=1
lstdy=16
lstyr=2008
ledmn=2
leddy=28
ledyr=2009
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
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; read SABER
; restore,'SABER_grid5_ALL_v2.0_20100701.sav
; GP              FLOAT     = Array[144, 96, 121, 2]
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; O3_127          FLOAT     = Array[144, 96, 121, 2]
; O3_96           FLOAT     = Array[144, 96, 121, 2]
; PRESSURE        FLOAT     = Array[121]
; TP              FLOAT     = Array[144, 96, 121, 2]
; U               FLOAT     = Array[144, 96, 121, 2]
; V               FLOAT     = Array[144, 96, 121, 2]
;
      ifile=dirs+'SABER_grid5_ALL_v2.0_'+sdate+'.sav'
        dum=findfile(ifile)
        if dum(0) ne '' then begin
           restore,ifile
           print,'read '+ifile
           spress=pressure
           TP_GRID=mean(TP,dim=4,/Nan)
           U_GRID=mean(U,dim=4,/Nan)
        endif
        if dum(0) eq '' then begin
           TP_GRID=0.*fltarr(144,96,121)
           U_GRID=0.*fltarr(144,96,121)
        endif
        sabertbar=mean(TP_GRID,dim=1)
        saberubar=mean(U_GRID,dim=1)
        slat=lat
;
; read gridded MLS data on pressure
; BRO             FLOAT     = Array[144, 96, 37, 2]
; CLO             FLOAT     = Array[144, 96, 37, 2]
; CO              FLOAT     = Array[144, 96, 37, 2]
; GPH             FLOAT     = Array[144, 96, 55, 2]
; H2O             FLOAT     = Array[144, 96, 55, 2]
; HCL             FLOAT     = Array[144, 96, 37, 2]
; HNO3            FLOAT     = Array[144, 96, 37, 2]
; HO2             FLOAT     = Array[144, 96, 49, 2]
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; N2O             FLOAT     = Array[144, 96, 37, 2]
; NODE            STRING    = Array[2]
; O3              FLOAT     = Array[144, 96, 55, 2]
; OH              FLOAT     = Array[144, 96, 49, 2]
; PMLS            FLOAT     = Array[37]
; PMLS2           FLOAT     = Array[55]
; PMLS3           FLOAT     = Array[49]
; T               FLOAT     = Array[144, 96, 55, 2]
; U               FLOAT     = Array[144, 96, 55, 2]
; V               FLOAT     = Array[144, 96, 55, 2]
;
        ifile=dir+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav'
        dum=findfile(ifile)
        if dum(0) ne '' then begin
           restore,ifile
           print,'read '+ifile
           TP_GRID=mean(T,dim=4,/Nan)
           U_GRID=mean(U,dim=4,/Nan)
        endif
        if dum(0) eq '' then begin
           TP_GRID=0.*fltarr(144,96,55)
           U_GRID=0.*fltarr(144,96,55)
        endif
        mlstbar=mean(TP_GRID,dim=1)
        mlsubar=mean(U_GRID,dim=1)
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
        print,'read '+ifile
        merratbar=mean(tgrd,dim=1)
        merraubar=mean(ugrd,dim=1)
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
        print,'read '+ifile
        sdwtbar=mean(t,dim=1)
        sdwubar=mean(u,dim=1)
;
; postscript file
;
        if setplot eq 'ps' then begin
           lc=0
           xsize=nxdim/100.
           ysize=nydim/100.
           set_plot,'ps'
           device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                  /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_temp+u_mls_merra_sdw_'+sdate+'.ps'
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
       xyouts,.4,.95,sdate,color=0,charsize=2,charthick=2,/normal
       xmn=xorig(0)
       xmx=xorig(0)+xlen
       ymn=yorig(0)
       ymx=yorig(0)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,mlstbar,mlat,pmls2,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='MLS'
       contour,mlstbar,mlat,pmls2,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow
index=where(abs(mlat) lt 1.)
mlsubar(index,*)=0./0.
       contour,mlsubar,mlat,pmls2,levels=[10.,20,30,40,50,60,70,80,90,100],/overplot,/noeras,c_color=0,/follow
       contour,mlsubar,mlat,pmls2,levels=[-100.,-90,-80,-70,-60,-50,-40,-30,-20,-10],/overplot,/noeras,c_color=mcolor,/follow,c_linestyle=5

       xmn=xorig(1)
       xmx=xorig(1)+xlen
       ymn=yorig(1)
       ymx=yorig(1)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,sabertbar,slat,spress,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='SABER'
       contour,sabertbar,slat,spress,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow
       contour,saberubar,slat,spress,levels=[10.,20,30,40,50,60,70,80,90,100],/overplot,/noeras,c_color=0,/follow
       contour,saberubar,slat,spress,levels=[-100.,-90,-80,-70,-60,-50,-40,-30,-20,-10],/overplot,/noeras,c_color=mcolor,/follow,c_linestyle=5

       xmn=xorig(2)
       xmx=xorig(2)+xlen
       ymn=yorig(2)
       ymx=yorig(2)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,merratbar,latitude_waccm,pressure,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='MERRA'
       contour,merratbar,latitude_waccm,pressure,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow
       contour,merraubar,latitude_waccm,pressure,levels=[10.,20,30,40,50,60,70,80,90,100],/overplot,/noeras,c_color=0,/follow
       contour,merraubar,latitude_waccm,pressure,levels=[-100.,-90,-80,-70,-60,-50,-40,-30,-20,-10],/overplot,/noeras,c_color=mcolor,/follow,c_linestyle=5

       xmn=xorig(3)
       xmx=xorig(3)+xlen
       ymn=yorig(3)
       ymx=yorig(3)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       contour,sdwtbar,lat,lev,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[1000.,1.e-5],/ylog,color=0,title='SD-WACCM'
       contour,sdwtbar,lat,lev,levels=[130.,140.,150.],/overplot,/noeras,c_color=mcolor,/follow
       contour,sdwubar,lat,lev,levels=[10.,20,30,40,50,60,70,80,90,100],/overplot,/noeras,c_color=0,/follow
       contour,sdwubar,lat,lev,levels=[-100.,-90,-80,-70,-60,-50,-40,-30,-20,-10],/overplot,/noeras,c_color=mcolor,/follow,c_linestyle=5

       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xorig(0),xorig(3)+xlen,ymnb,ymxb
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
; plot Ubar as a function of latitude in the middle mesosphere
;
       !type=2^2+2^3
       set_viewport,min(xorig),max(xorig)+xlen,0.1,0.4
       rlev=0.001
index=where(finite(spress) ne 1)
if index(0) ne -1L then spress(index)=-99.
       index=where(abs(spress-rlev) eq min(abs(spress-rlev)))
       sdat=reform(saberubar(*,index(0)))
       index=where(abs(slat) le 10)
       sdat(index)=0./0.
       plot,slat,sdat,thick=5,ytitle='Ubar at '+strcompress(rlev,/r),color=0,xtitle='Latitude',yrange=[-100,100]
       plots,slat(0),0
       plots,slat(-1),0,thick=5,linestyle=5,color=0,/continue
loadct,0
       index=where(abs(pmls2-rlev) eq min(abs(pmls2-rlev)))
       mdat=reform(mlsubar(*,index(0)))
       index=where(abs(mlat) le 10)
       mdat(index)=0./0.
       oplot,mlat,mdat,thick=5,color=200
loadct,39
       if rlev gt min(pressure) then begin
          index=where(abs(pressure-rlev) eq min(abs(pressure-rlev)))
          oplot,latitude_waccm,reform(merraubar(*,index(0))),thick=5,color=250
       endif
       index=where(abs(lev-rlev) eq min(abs(lev-rlev)))
       oplot,lat,reform(sdwubar(*,index(0))),thick=5,color=50
       xyouts,-88,90,'SABER',/data,charthick=2,color=0
loadct,0
       xyouts,-88,70,'MLS',/data,charthick=2,color=200
loadct,39
       xyouts,-88,50,'MERRA',/data,charthick=2,color=250
       xyouts,-88,30,'SDW',/data,charthick=2,color=50
;
; Close PostScript file and return control to X-windows
       if setplot ne 'ps' then stop	;wait,1
       if setplot eq 'ps' then begin
          device, /close
          spawn,'convert -trim yz_temp+u_mls_merra_sdw_'+sdate+'.ps -rotate -90 '+$
                              'yz_temp+u_mls_merra_sdw_'+sdate+'.jpg'
          spawn,'rm -f yz_temp+u_mls_merra_sdw_'+sdate+'.ps'
       endif

skip:
      icount=icount+1L
goto,jump

end
