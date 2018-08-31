;
; save daily 50 hPa geopotential height from ECMWF and MetO for Cora
; VLH 9/23/5
;
@stddat
@kgmt
@ckday
@kdate
@rd_ecmwf
@rd_ukmo

dir='/aura5/harvey/ECMWF_data/Datfiles/ecmwf_'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
lstmn=11 & lstdy=1 & lstyr=91
ledmn=8 & leddy=31 & ledyr=2
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 57 then lstyr=lstyr+2000
if ledyr lt 57 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1957 then stop,'Year out of range '
if ledyr lt 1957 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

; Compute initial Julian date
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
      print,imn,idy,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
;
;***Read ECMWF data
      restore,dir+smn+'_'+sdy+'_'+syr+'_12Z.sav'
;
;***Read UKMO data
      file='/aura4/harvey/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nlg,nlat,nlv,ulon,ulat,wlon,wlat,p,g3d,t3d,u3d,v3d
      if iflg ne 0 then goto, jump
      if icount eq 0L then begin 
         icount=1L
         rsfc=50.
;        print,double(press)
;        read,' Enter desired pressure surface ',rsfc
         index=where(rsfc eq press)
         asfc=alog(rsfc)
         lsfc=index(0)
         longitude_ecmwf=alon
         latitude_ecmwf=alat
         longitude_meto=ulon
         latitude_meto=ulat
      endif

; save daily 50 hPa geopotential height
      gp50hpa_ecmwf=reform(gp(*,*,lsfc))
;
; interpolate MetO to 50 hPa
;
      gp50hpa_meto=0.*reform(g3d(*,*,0))
      for l=0L,nlv-2L do begin
          lp1=l+1L
          p0=alog(p(l))
          p1=alog(p(lp1))
          if asfc le p0 and asfc ge p1 then begin
             zscale=(asfc-p0)/(p1-p0)
             gp0=reform(g3d(*,*,l))
             gp1=reform(g3d(*,*,lp1))
             gp50hpa_meto=gp0+zscale*(gp1-gp0)
;print,p0,asfc,p1,zscale
;print,gp0(0,0),gp50hpa_meto(0,0),gp1(0,0)
          endif
      endfor
      save,file='/aura3/data/ECMWF_data/Datfiles/gp50hpa_ecmwf_meto_'+smn+'_'+sdy+'_'+syr+'_12Z.sav',$
           longitude_ecmwf,latitude_ecmwf,longitude_meto,latitude_meto,gp50hpa_ecmwf,gp50hpa_meto
goto, jump
end
