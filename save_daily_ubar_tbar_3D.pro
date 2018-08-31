;
; save daily UKMO Ubar in IDL save format
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nwp
;
; Ask interactive questions- get starting/ending dates
;
!type=2^2+2^3
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
lstmn=9 & lstdy=28 & lstyr=1991
ledmn=5 & leddy=3 & ledyr=2015
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate=strarr(kday)
ofile='/atmos/harvey/UKMO_data/Datfiles/ukmo_12Z_Ubar_Tbar_3D.sav'
icount=0L
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
      date=syr+smn+sdy
      sdate(icount)=date
;
; have Met Office pressure data from Sep 28th 1991 to Nov 10th 2007 in IDL save format with filenames like this - ppassm_y07_m04_d25_h12.pp.sav
; NOTE - latitude goes from North to South
;
; ALAT            FLOAT     = Array[73]
; ALON            FLOAT     = Array[96]
; P               FLOAT     = Array[22]
; T3D             FLOAT     = Array[96, 73, 22]
; U3D             FLOAT     = Array[96, 72, 22]
; V3D             FLOAT     = Array[96, 72, 22]
; WLAT            FLOAT     = Array[72]
; WLON            FLOAT     = Array[96]
; Z3D             FLOAT     = Array[96, 73, 22]
;
; have Met Office pressure data from Jan 1st 2007 to the present in netCDF format with filenames like this - ukmo-nwp-strat_gbl-std_2007010112_u-v-gph-t-w_uars.nc
;
; read save data prior to 2007
; read netcdf data beginning in 2007
;
      if iyr le 2006 then begin
         restore,'/atmos/harvey/UKMO_data/Datfiles/ppassm_y'+syr1+'_m'+smn+'_d'+sdy+'_h12.pp.sav'
         nc=n_elements(wlon)
         nr=n_elements(wlat)
         nr1=n_elements(alat)
         nlv=n_elements(p)
;
; reverse latitude
;
         u3d_sorted=0.*u3d
         t3d_sorted=0.*t3d
         index=sort(wlat)
         index2=sort(alat)
         for i=0L,nc-1L do begin
             for k=0L,nlv-1L do begin
                 u3d_sorted(i,*,k)=u3d(i,index,k)
                 t3d_sorted(i,*,k)=t3d(i,index2,k)
             endfor
         endfor
         u3d=u3d_sorted
         t3d=t3d_sorted
         wlat=wlat(index)
         alat=alat(index2)
      endif
      if iyr gt 2006 then begin
         ifile='/atmos/harvey/UKMO_data/Datfiles/ukmo-nwp-strat_gbl-std_'+date+'12_u-v-gph-t-w_uars.nc'
         iflg=0L
         rd_ukmo_nwp,ifile,nc,nr,nc1,nr1,nlv,wlon,alon,wlat,alat,p,z3d,t3d,u3d,v3d,iflg
         if iflg gt 0L then goto,jumpday
      endif
;
; remove data above 1 hPa
;
      good=where(p ge 1.,nlv)
      u3d=reform(u3d(*,*,good))
      t3d=reform(t3d(*,*,good))
      p=reform(p(good))
      if min(p) ne 1.0 then goto,jumpday
      print,'read Met Office data on '+date,nr,nlv
      if icount eq 0L then begin
         ubar=-99.+0.*fltarr(nr,nlv,kday)
         tbar=-99.+0.*fltarr(nr1,nlv,kday)
      endif
;
; zonal mean U and T
;
       for k=0L,nlv-1 do begin
          for j=0L,nr-1L do ubar(j,k,icount)=mean(u3d(*,j,k))
          for j=0L,nr1-1L do tbar(j,k,icount)=mean(t3d(*,j,k))
       endfor
       contour,ubar(*,*,icount),wlat,p,/ylog,xrange=[-90.,90.],yrange=[1000.,1.],title='Ubar '+date,$
               xtitle='Latitude',ytitle='Pressure (hPa)',charsize=2,charthick=2,xticks=6,/nodata
       contour,ubar(*,*,icount),wlat,p,/ylog,/overplot,levels=10+10*findgen(20)
       contour,ubar(*,*,icount),wlat,p,/ylog,/overplot,levels=-200+10*findgen(20),c_linestyle=5
       xyouts,-10.,70.,month(imn-1),/data,charsize=2,charthick=2
;print,min(p),nlv
      jumpday:
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
save,filename=ofile,ubar,tbar,alat,wlat,p,sdate
end
