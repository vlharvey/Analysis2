;
; check for missing MERRA data
;
@stddat
@kgmt
@ckday
@kdate

lstmn=3
lstdy=1
lstyr=1979
ledmn=2
leddy=29
ledyr=2012
lsfc=0
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;
; /Volumes/Data/MERRA_data/Datfiles/MERRA300.prod.assim.inst6_3d_ana_Nv.20111231.SUB.nc
; /Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_YYYYMMDD.sav
; /Volumes/Data/MERRA_data/Datfiles/MERRA300.prod.assim.tavg3_3d_rad_Cp.20111231.SUB.nc
;
dir='/Volumes/Data/MERRA_data/Datfiles/'
pre0='MERRA*.prod.assim.inst6_3d_ana_Nv.'
pre1='MERRA-on-WACCM_'
pre2='MERRA*.prod.assim.tavg3_3d_rad_Cp.'

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '

      syr=string(FORMAT='(i4)',iyr)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy

      dum1=file_search(dir+pre1+sdate+'.sav',count=count1)
      dum2=file_search(dir+pre2+sdate+'.SUB.nc.gz',count=count2)
      if dum1(0) eq '' then begin
               dum1=file_search(dir+pre0+sdate+'.SUB.nc.gz',count=count1)
               if dum1(0) eq '' then print,'missing analysis '+sdate
      endif
      if dum2(0) eq '' then print,'missing rad '+sdate
goto,jump
end
