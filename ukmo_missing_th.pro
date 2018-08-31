; IDL code to view horizontal slices (on pressure surfaces) of UKMO
; Z,T,u,v data.  Plots are displayed in a polar stereographic view.

@stddat
@kgmt
@ckday
@kdate
@date2uars

dir='/usr71/users/ukmo/ukmo_'
dum=strarr(1)
uday=0L
lstmn=0
lstdy=0
lstyr=0
ledmn=0
leddy=0
ledyr=0
lsfc=0
lstday=0
ledday=0
month=['jan','feb','mar','apr','may','jun',$
       'jul','aug','sep','oct','nov','dec']
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      UKMO Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr

if lstyr ge 91 and lstyr le 99 then lstyr=lstyr+1900
if ledyr ge 91 and ledyr le 99 then ledyr=ledyr+1900
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1991 or lstyr gt 2000 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 2000 then stop,'Year out of range '
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

      if iyr lt 2000 then iyr1=iyr-1900
      if iyr ge 2000 then iyr1=iyr-2000
      ifile=dir+month(imn-1)+'_'+(string(FORMAT='(I2.2)',idy))+'_'+$
            (string(FORMAT='(I2.2)',iyr1))+'.nc3'
      dum=findfile(ifile)
      if dum(0) eq '' then begin
         z = date2uars(imn,idy,iyr,uday)
         print,imn,idy,iyr,'  UARS DAY= ',uday
      endif
      goto, jump
end
