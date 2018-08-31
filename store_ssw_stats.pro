;
; store ssw dates.  use standard WMO definition and explore using Elat
;
@stddat
@kgmt
@ckday
@kdate

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.10,0.55,0.10,0.55]
yorig=[0.55,0.55,0.15,0.15]
xlen=0.35
ylen=0.35
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

idir='/aura3/data/WACCM_data/Datfiles/wa3_tnv3_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1
lstdy=1
lstyr=1990
ledmn=1
leddy=1
ledyr=2004
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      UKMO Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 50 then lstyr=lstyr+2000
if ledyr lt 50 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
yyyymmdd=lonarr(kday)

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

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      syr=string(FORMAT='(i4)',iyr)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      ifile=idir+string(FORMAT='(i4,i2.2,i2.2,a4)',iyr,imn,idy,'.sav')
      dum=findfile(ifile)
      if dum(0) eq '' then goto,jump
      restore,ifile
      nc=n_elements(alon) & nr=n_elements(alat) & nth=n_elements(lev)
      print,long(syr+smn+sdy)

      yyyymmdd(icount)=long(syr+smn+sdy)
      if icount eq 0L then begin
         nhdt=fltarr(kday,nth)
         shdt=fltarr(kday,nth)
         nhu60=fltarr(kday,nth)
         shu60=fltarr(kday,nth)
      endif
;
; zonal mean T and U
;
      tzm=fltarr(nr,nth)
      uzm=fltarr(nr,nth)
      for k=0L,nth-1L do begin
          for j=0L,nr-1L do begin
              tzm(j,k)=total(tgrd(*,j,k))/float(nc)
              uzm(j,k)=total(ugrd(*,j,k))/float(nc)
          endfor
      endfor
;
; meridional temperature gradient and zonal mean wind at 60
;
      p85=where(abs(alat-85.) eq min(abs(alat-85.)))
      p65=where(abs(alat-65.) eq min(abs(alat-65.)))
      p60=where(abs(alat-60.) eq min(abs(alat-60.)))
      n85=where(abs(alat+85.) eq min(abs(alat+85.)))
      n65=where(abs(alat+65.) eq min(abs(alat+65.)))
      n60=where(abs(alat+60.) eq min(abs(alat+60.)))
      for k=0L,nth-1L do begin
          nhdt(icount,k)=tzm(p85(0),k)-tzm(p60(0),k)
          shdt(icount,k)=tzm(n85(0),k)-tzm(n60(0),k)
          nhu60(icount,k)=uzm(p65(0),k)
          shu60(icount,k)=uzm(n65(0),k)
      endfor
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
index=where(yyyymmdd gt 0L)
yyyymmdd=yyyymmdd(index)
nhdt=nhdt(index,*)
shdt=shdt(index,*)
nhu60=nhu60(index,*)
shu60=shu60(index,*)
save,file='wa3_wmo_ssw_diagnostics.sav',yyyymmdd,lev,nhdt,shdt,nhu60,shu60

end
