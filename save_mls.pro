;
; convert raw MLS data to daily IDL save files
;
@stddat
@kgmt
@ckday
@kdate
@date2uars

lstmn=10L & lstdy=15L & lstyr=91L 
ledmn=1L & leddy=1L & ledyr=99L
lstday=0L & ledday=0L & uday=0L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 or lstyr gt 2011 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 2011 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
dir='/aura6/data/MLS_data/Datfiles/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Calculate UARS day from (imn,idy,iyr) information.
      z = date2uars(imn,idy,iyr,uday)

      z = stddat(imn,idy,iyr,ndays)
      if ndays gt ledday then stop,' Normal termination condition '
;
; read v4 MLS ozone and temperature
;
      tpfile=dir+'MLS_L3AT_STEMP_D'+string(FORMAT='(i4.4)',uday)+'.V0004_C01_PROD'
      o3file=dir+'MLS_L3AT_SO3_205_D'+string(FORMAT='(i4.4)',uday)+'.V0004_C01_PROD'
      dum=findfile(tpfile)
      if dum eq '' then goto,jump
      l3atread,temp,tpfile
      l3atread,o3,o3file
      print,'read '+o3file
      mtime=o3.SOLTIME
      mlon=o3.LONG
      mlat=o3.LAT
;
; The standard pressure level values in millibars are given by:
; P(i) = 1000.0 * (10**(-i/6)), i=0,1,2...
;
      no3lev=o3(0).TOTPTS
      ntplev=temp(0).TOTPTS
      press=1000.*10.^(-1.*findgen(no3lev)/6.)
      p0=press(o3(0).spare)
      o3press=p0*10.^(-1.*findgen(no3lev)/6.)		; base of profile is given by o3.spare
      mo3=o3.qu*1.e6
      mo3err=o3.err*1.e6
      mtp=temp.qu
      press=1000.*10.^(-1.*findgen(ntplev)/6.)
      p0=press(temp(0).spare)
      tppress=p0*10.^(-1.*findgen(ntplev)/6.)         ; base of profile is given by temp.spare
;
; save IDL save file
; 
      mfile=dir+'mls_'+mon(imn-1)+string(FORMAT='(i2.2,a1,i4.4)',idy,'_',iyr)+'.sav'
      save,file=mfile,mtime,mlon,mlat,mo3,mo3err,o3press,mtp,tppress
print,'saved '+mfile
goto,jump
end
