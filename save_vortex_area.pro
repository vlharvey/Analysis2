;
; calculate the area enclosed by the vortex at a given theta level
; expressed in equivalent latitude and plot winter time series
; with average and sigma
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nc3

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=0L & lstdy=0L & lstyr=0L & ledmn=0L
leddy=0L & ledyr=0L & lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
;
; Ask interactive questions- get starting/ending date and p surface
;
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
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
      if ndays eq ledday then begin
         save,area_ave,area_num,area_all,filename='arctic_vortex_area.sav'
         stop
      endif

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy

      ifile=mon(imn-1)+sdy+'_'+uyr
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      dum1=findfile(diru+ifile+'.nc3')
      if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
      if dum1(0) eq '' then goto,jump
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      marksf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,10,marksf2
      ncdf_close,ncid

      if icount eq 0L then begin
         icount=1L
         area_ave=fltarr(366,nth)
         area_num=lonarr(366,nth)
         area_all=fltarr(15,366,nth)
         dum=transpose(marksf2(*,*,0))
         lon=0.*dum
         lat=0.*dum
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         nrr=91
         yeq=findgen(nrr)
         latcircle=fltarr(nrr)
         hem_frac=fltarr(nrr)
         for j=0,nrr-2 do begin
             hy=re*dtr
             dx=re*cos(yeq(j)*dtr)*360.*dtr
             latcircle(j)=dx*hy	; area in each latitude circle
         endfor
         for j=0,nrr-1 do begin
             index=where(yeq ge yeq(j))
  
; fraction of the hemisphere of each latitude circle
             if index(0) ne -1 then $
                hem_frac(j)=100.*total(latcircle(index))/hem_area
             if yeq(j) eq 0. then hem_frac(j)=100.
         endfor
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy	; area of each grid point
         endfor
      endif
;
; sum area of gridpoints in Arctic vortex
;
      myr=iyr-1991L
      for thlev=0,nth-1 do begin
          mark=transpose(marksf2(*,*,thlev))
          index=where(lat gt 0. and mark eq 1.0)
          if index(0) ne -1 then begin
             a0=total(area(index))
             area_ave(iday-1,thlev)=area_ave(iday-1,thlev)+a0
             area_num(iday-1,thlev)=area_num(iday-1,thlev)+1L
             area_all(myr,iday-1,thlev)=a0
          endif
      endfor
if area_num(iday-1,0) gt 0L then print,ifile,area_ave(iday-1,0)/area_num(iday-1,0)
if area_num(iday-1,0) eq 0L then print,ifile
goto,jump
end
