;
; save daily 50 hPa geopotential height from WACCM3 for Annika Seppala
; VLH 12/14/5
;
@stddat
@kgmt
@ckday
@kdate

dir='/aura3/data/WACCM_data/Datfiles/'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
lstmn=1 & lstdy=1 & lstyr=90
ledmn=1 & leddy=1 & ledyr=4
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
;***Read WACCM3 data
;
      ifile=dir+'wa3_tnv3_'+syr+smn+sdy+'.sav'
      dum=findfile(ifile)
      if dum(0) eq '' then goto,jump
      restore,ifile
      if icount eq 0L then begin
         icount=1L
         rsfc=50.
         asfc=alog(rsfc)
         p=lev
         nlv=n_elements(p)
         longitude_wa3=alon
         latitude_wa3=alat
         ncw=n_elements(longitude_wa3)
         nrw=n_elements(latitude_wa3)
      endif
;
; interpolate WACCM3 to 50 hPa
;
      gp50hpa_wa3=0.*reform(zgrd(*,*,0))
      for l=1L,nlv-2L do begin
          lp1=l-1L
          p0=alog(p(l))
          p1=alog(p(lp1))
          if asfc le p0 and asfc ge p1 then begin
             zscale=(asfc-p0)/(p1-p0)
             gp0=reform(zgrd(*,*,l))
             gp1=reform(zgrd(*,*,lp1))
             gp50hpa_wa3=gp0+zscale*(gp1-gp0)
;print,p(l),rsfc,p(lp1)
;print,p0,asfc,p1,zscale
;print,gp0(0,0),gp50hpa_wa3(0,0),gp1(0,0)
;stop
          endif
      endfor
      ofile=dir+'gp50hpa_wa3_'+syr+smn+sdy+'.ASCII'
      print,'writing '+ofile
      openw,1,ofile
      printf,1,nrw
      printf,1,ncw
      printf,1,ncw,nrw
      printf,1,longitude_wa3
      printf,1,latitude_wa3
      printf,1,gp50hpa_wa3
      close,1
goto, jump
end
