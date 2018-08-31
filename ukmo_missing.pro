; IDL code to view horizontal slices (on pressure surfaces) of UKMO
; Z,T,u,v data.  Plots are displayed in a polar stereographic view.

@stddat
@kgmt
@ckday
@kdate
@date2uars

dir='/mars/couk'
dir2='/usr71/users/ukmo/couk'
lstmn=0
lstdy=0
lstyr=0
ledmn=0
leddy=0
ledyr=0
lsfc=0
lstday=0
ledday=0
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

; --- Calculate UARS day from (imn,idy,iyr) information.
      z = date2uars(imn,idy,iyr,uday)

      ifile=dir+(string(FORMAT='(I4.4,A4)',fix(uday),'.dat'))
      dum=findfile(ifile)
      ifile2=dir2+(string(FORMAT='(I4.4,A4)',fix(uday),'.dat'))
      dum2=findfile(ifile2)
      if dum(0) eq '' and dum2(0) eq '' then $
         print,imn,idy,iyr,' = UARS day ',fix(uday)
      goto, jump
end
