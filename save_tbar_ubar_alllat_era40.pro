;
; NOTE: ERA40 LONGITUDES ARE -180 TO 180
; ERA40 zonal mean temperature and zonal wind at all latitudes for 45 year record
;
; VLH 9/10/09
;
@stddat
@kgmt
@ckday
@kdate
@rd_era40_nc

dir='/aura7/harvey/ERA40_data/Datfiles/era40_ua_12Z_'
lstmn=9L & lstdy=1L & lstyr=1957L
ledmn=8L & leddy=31L & ledyr=2002L
;
; Get start and end dates
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1950 or lstyr gt 2002 then stop,'Year out of range '
if ledyr lt 1950 or ledyr gt 2002 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd_all=lonarr(nfile)
syyyymmdd_all=strarr(nfile)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' Starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
      yyyymmdd_all(icount)=long(syr+smn+sdy)
      syyyymmdd_all(icount)=date
;
; read ERA40 pressure data, i.e. /aura7/harvey/ERA40_data/Datfiles/era40_ua_12Z_19600131.nc
;
      rd_era40_nc,dir+date+'.nc',nc,nr,nl,alon,alat,press,tp,uu,vv,gp,iflg
      if iflg eq 1 then goto,skip
;
; declare time period arrays
;
      if kcount eq 0L then begin
         tbar_all=-9999.+0.*fltarr(nfile,nr,nl)
         ubar_all=-9999.+0.*fltarr(nfile,nr,nl)
         gbar_all=-9999.+0.*fltarr(nfile,nr,nl)
         kcount=1
      endif
;
; calculate zonal mean temperature and zonal wind
;
      uzm=-9999.+0.*fltarr(nr,nl)
      tzm=-9999.+0.*fltarr(nr,nl)
      gzm=-9999.+0.*fltarr(nr,nl)
      for k=0,nl-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tp(*,j,k))/float(nc)
              uzm(j,k)=total(uu(*,j,k))/float(nc)
              gzm(j,k)=total(gp(*,j,k))/float(nc)
              if tzm(j,k) gt 400. or tzm(j,k) lt 100. then stop,'check temperature'
              if uzm(j,k) lt -200. or uzm(j,k) gt 200. then stop,'check zonal wind'
          endfor
      endfor
;
; retain all daily zonal means 
;
      tbar_all(icount,*,*)=tzm
      ubar_all(icount,*,*)=uzm
      gbar_all(icount,*,*)=gzm
;
; skip if file is missing but still increment day
;
      skip:
      icount=icount+1L
goto,jump
saveit:

syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
sday=strmid(syyyymmdd_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nday=icount
save,file='ERA40_tbar_ubar_alllat_era40_'+yearlab+'.sav',nl,nr,alat,nday,yyyymmdd_all,press,tbar_all,ubar_all,gbar_all
end
