;
; save daily UKMO Ubar at 50 hPa in IDL save format
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nwp
;
; update existing file
;
restore,'/Users/harvey/Benze/ukmo_12Z_Ubar_50hPa.sav
good=where(ubar(10,*) ne 0.)
lastday=max(sdate(good))
lstyr=long(strmid(lastday,0,4))
lstmn=long(strmid(lastday,4,2))
lstdy=long(strmid(lastday,6,2))
print,'Last Day ',lstyr,lstmn,lstdy

!type=2^2+2^3
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
;lstmn=9 & lstdy=28 & lstyr=1991
ledmn=10 & leddy=6 & ledyr=2013
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nr=n_elements(wlat)
sdate_save=strarr(n_elements(sdate)+kday)
ubar_save=fltarr(nr,n_elements(sdate)+kday)
sdate_save(0:n_elements(sdate)-1)=sdate
ubar_save(*,0:n_elements(sdate)-1)=ubar

ofile='/Users/harvey/Benze/ukmo_12Z_Ubar_50hPa.sav'
icount=n_elements(sdate)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
;iday = iday - 1
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

print,iyr,imn,idy
; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(i4)',iyr)
      syr1=strmid(syr,2,2)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      date=syr+smn+sdy
;     print,icount,' ',date
      sdate_save(icount)=date
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
         restore,'/Volumes/earth/harvey/UKMO_data/Datfiles/ppassm_y'+syr1+'_m'+smn+'_d'+sdy+'_h12.pp.sav'
         nc=n_elements(wlon)
         nr=n_elements(wlat)
         nlv=n_elements(p)
;
; reverse latitude
;
         u3d_sorted=0.*u3d
         index=sort(wlat)
         for i=0L,nc-1L do begin
             for k=0L,nlv-1L do begin
                 u3d_sorted(i,*,k)=u3d(i,index,k)
             endfor
         endfor
         u3d=u3d_sorted
         wlat=wlat(index)
      endif
      if iyr gt 2006 then begin
         ifile='/Volumes/earth/harvey/UKMO_data/Datfiles/ukmo-nwp-strat_gbl-std_'+date+'12_u-v-gph-t-w_uars.nc'
         iflg=0L
         rd_ukmo_nwp,ifile,nc,nr,nc1,nr1,nlv,wlon,alon,wlat,alat,p,z3d,t3d,u3d,v3d,iflg
         if iflg gt 0L then goto,jumpday
      endif
      print,'read Met Office data on '+date
      if icount eq 0L then ubar=-99.+0.*fltarr(nr,kday)
;
; zonal mean U at pressure level closest to 50 hPa
;
      zindex=where(abs(p-50.) eq min(abs(p-50.)))
      if zindex(0) ne -1L then begin
         for j=0L,nr-1L do ubar_save(j,icount)=total(u3d(*,j,zindex(0)))/float(nc)
      endif
      plot,wlat,ubar_save(*,icount),xrange=[-90.,90.],yrange=[-40.,80.],title=date,xtitle='Latitude',ytitle='Ubar (m/s)',charsize=2,charthick=2,xticks=6
      xyouts,-10.,70.,month(imn-1),/data,charsize=2,charthick=2
      jumpday:
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
sdate=sdate_save
ubar=reform(ubar_save(0:nr-1,*))
save,filename=ofile,ubar,wlat,sdate
end
