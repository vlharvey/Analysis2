;
; save daily MLS zonal means for entire data record
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=255
device,decompose=0
;
; Ask interactive questions- get starting/ending dates
;
!type=2^2+2^3
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
lstmn=8 & lstdy=8 & lstyr=2004
ledmn=5 & leddy=15 & ledyr=2017
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
icount=0L
kcount=0L
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(i4)',iyr)
      syr1=strmid(syr,2,2)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
      sdate_all(icount)=sdate
      if icount eq 0L then sdate0=sdate
;
; restore twice daily gridded data
; UT_MEAN         INT       = Array[2]                                                                  
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; PMLS            FLOAT     = Array[37]
; PMLS2           FLOAT     = Array[55]
; PMLS3           FLOAT     = Array[49]
; BRO             FLOAT     = Array[144, 96, 37, 2]
; CLO             FLOAT     = Array[144, 96, 37, 2]
; CO              FLOAT     = Array[144, 96, 37, 2]
; HCL             FLOAT     = Array[144, 96, 37, 2]
; HNO3            FLOAT     = Array[144, 96, 37, 2]
; HO2             FLOAT     = Array[144, 96, 37, 2]
; N2O             FLOAT     = Array[144, 96, 37, 2]
; OH              FLOAT     = Array[144, 96, 49, 2]
; T               FLOAT     = Array[144, 96, 55, 2]
; U               FLOAT     = Array[144, 96, 55, 2]
; V               FLOAT     = Array[144, 96, 55, 2]
; GPH             FLOAT     = Array[144, 96, 55, 2]
; H2O             FLOAT     = Array[144, 96, 55, 2]
; O3              FLOAT     = Array[144, 96, 55, 2]
;
      dum=findfile('/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_U_V_2xdaily_v4.2_'+sdate_all(icount)+'.sav')
      if dum(0) eq '' then goto,jumpday
      restore,dum
      print,'read MLS data on '+sdate
;
; on the first day of data
;
      if kcount eq 0L then begin
         press55=PMLS2
         press37=PMLS
         nr=n_elements(lat)
         nlv=n_elements(press55)
         nlv2=n_elements(press37)
         ubar=0./0.*fltarr(nr,nlv,kday)
         vbar=0./0.*fltarr(nr,nlv,kday)
         tbar=0./0.*fltarr(nr,nlv,kday)
         zbar=0./0.*fltarr(nr,nlv,kday)
         h2obar=0./0.*fltarr(nr,nlv,kday)
         o3bar=0./0.*fltarr(nr,nlv,kday)
         cobar=0./0.*fltarr(nr,nlv2,kday)
         n2obar=0./0.*fltarr(nr,nlv2,kday)
         kcount=1
      endif
;
; daily 3d grids (average over two local times)
;
      umean=mean(u,dim=4,/Nan)
      vmean=mean(v,dim=4,/Nan)
      zmean=mean(gph,dim=4,/Nan)
      tmean=mean(t,dim=4,/Nan)
      h2omean=mean(h2o,dim=4,/Nan)
      o3mean=mean(o3,dim=4,/Nan)
      comean=mean(co,dim=4,/Nan)
      n2omean=mean(n2o,dim=4,/Nan)
;
; zonal means
;
      ubar(*,*,icount)=mean(umean,dim=1,/Nan)
      vbar(*,*,icount)=mean(vmean,dim=1,/Nan)
      tbar(*,*,icount)=mean(tmean,dim=1,/Nan)
      zbar(*,*,icount)=mean(zmean,dim=1,/Nan)
      h2obar(*,*,icount)=mean(h2omean,dim=1,/Nan)
      o3bar(*,*,icount)=mean(o3mean,dim=1,/Nan)
      cobar(*,*,icount)=mean(comean,dim=1,/Nan)
      n2obar(*,*,icount)=mean(n2omean,dim=1,/Nan)

;contour,ubar(*,*,icount),lat,reform(zbar(*,*,icount))/1000.,xrange=[-90.,90.],yrange=[10.,100.],title='Ubar '+sdate,$
;        xtitle='Latitude',ytitle='Altitude (km)',charsize=2,charthick=2,xticks=6,/nodata
;contour,ubar(*,*,icount),lat,reform(zbar(*,*,icount))/1000.,/overplot,levels=10+10*findgen(20)
;contour,ubar(*,*,icount),lat,reform(zbar(*,*,icount))/1000.,/overplot,levels=-200+10*findgen(20),c_linestyle=5
;contour,ubar(*,*,icount),lat,reform(zbar(*,*,icount))/1000.,/overplot,levels=[0],thick=3
;cobar=cobar*1.e6
;cobar0=reform(cobar(*,*,icount))
;dcobar0=0.*cobar0
;for k=0L,nlv2-1L do dcobar0(*,k)=deriv(cobar0(*,k))/mean(cobar0(*,k),/Nan)	;  normalized
;contour,cobar(*,*,icount),lat,alog(1000./press37)*7.,/overplot,levels=[0.1,0.5,1+findgen(18)],c_color=1+(findgen(20)/20.)*255.,thick=5
;contour,dcobar0,lat,alog(1000./press37)*7.,/overplot,levels=[-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2],c_color=[50,60,70,80,90,100,110,120,130],thick=5
;contour,dcobar0,lat,alog(1000./press37)*7.,/overplot,levels=0.2+0.1*findgen(9),c_color=200+5*findgen(9),thick=5
;wait,.5

      jumpday:
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
sdate1=sdate
datelab=sdate0+'-'+sdate1
ofile='mls_daily_zonal_means_from_gridded_'+datelab+'.sav'
save,filename=ofile,sdate_all,kday,nr,nlv,nlv2,lat,press55,press37,ubar,vbar,tbar,zbar,h2obar,o3bar,cobar,n2obar
end
