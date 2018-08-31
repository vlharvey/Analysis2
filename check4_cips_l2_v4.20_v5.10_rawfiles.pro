;
; check to make sure that there is a corresponding v5.1 file for each v4.2 file (1 file per orbit)
; -prep step for: compare CIPS level 2 v4.2 and v5.10
; VLH 1/18/2017
;
@stddat
@kgmt
@ckday
@kdate

pth='/atmos/harvey/CIPS_data/Datfiles/Level_2/cips_sci_2_orbit_'
;
; loop over years
;
for iyear=2007,2016 do begin
syear=strcompress(long(iyear),/r)

lstmn=1
lstdy=1
lstyr=iyear
ledmn=12
leddy=31
ledyr=iyear
lstday=0
ledday=0

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
;
; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotyear
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sday=string(FORMAT='(I3.3)',iday)
      sdate=syr+smn+sdy
      if smn eq '03' or smn eq '04' or smn eq '09' or smn eq '10' then goto,jump
;
; convert YYYYMMDD to DFS
;
      jday=julday(smn,sdy,syr)
      yyyy=iyr
      dfs=jday-julday(6,21,yyyy)

      if iday lt 120L or iday gt 300L then begin
         if iday lt 120L then yyyy=iyr-1L
         dfs=jday-julday(12,21,yyyy)
      endif
;
; get nc filenames on this day (I get an error if I try to read the gzipped file)
;
      if dfs ge -40 and dfs le 80 then begin

         fnames20=file_search(pth+'*'+syr+'-'+sday+'_v04.20_r05_cld.nc',count=n20)
         fnamescat20=file_search(pth+'*'+syr+'-'+sday+'_v04.20_r05_cat.nc',count=ncat20)
         fnames51=file_search(pth+'*'+syr+'-'+sday+'_v05.10_r01_cld.nc',count=n51)
         fnamescat51=file_search(pth+'*'+syr+'-'+sday+'_v05.10_r01_cat.nc',count=ncat51)
;
; missing cat or cld files
;
         if n20 ne ncat20 then print,sdate,' ',sday,' v4.2 ncld= ',n20,'  ncat= ',ncat20
         if n51 ne ncat51 then print,sdate,' ',sday,' v5.1 ncld= ',n51,'  ncat= ',ncat51
;
; check for different number of orbit files between versions
;
         if n20 ne n51 then print,sdate,' ',sday,' v4.2= ',n20,ncat20,'  v5.1= ',n51,ncat51
      endif
goto,jump

plotyear:
endfor	; loop over years
end
