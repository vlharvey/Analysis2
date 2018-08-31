;
; mee01fco run
; store number of minor ssw per year.  use standard WMO definition
; (no winds for major warmings)
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
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif

idir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_Mills/mee01fco/mee01fco.vars.h3.'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1
lstdy=1
lstyr=2028
ledmn=1
leddy=31
ledyr=2050
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
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
kcount=0L

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
alon=longitude
alat=latitude
lev=pressure
      nc=n_elements(alon) & nr=n_elements(alat) & nth=n_elements(lev)
      print,long(syr+smn+sdy)

      yyyymmdd(icount)=long(syr+smn+sdy)
      if kcount eq 0L then begin
         nhdt=fltarr(kday,nth)
         shdt=fltarr(kday,nth)
         nhu60=fltarr(kday,nth)
         shu60=fltarr(kday,nth)
         kcount=1
      endif
;
; zonal mean T and U
;
tgrd=temp
      tzm=fltarr(nr,nth)
      uzm=fltarr(nr,nth)
      for k=0L,nth-1L do begin
          for j=0L,nr-1L do begin
              tzm(j,k)=total(tgrd(*,j,k))/float(nc)
;             uzm(j,k)=total(ugrd(*,j,k))/float(nc)
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
;         nhu60(icount,k)=uzm(p65(0),k)
;         shu60(icount,k)=uzm(n65(0),k)
      endfor
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
index=where(yyyymmdd gt 0L,kday)
yyyymmdd=yyyymmdd(index)
nhdt=nhdt(index,*)
shdt=shdt(index,*)
nhu60=nhu60(index,*)
shu60=shu60(index,*)
save,file='waccm_mee01fco_wmo_ssw_diagnostics.sav',yyyymmdd,lev,nhdt,shdt,nhu60,shu60

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=11
col1=1+indgen(nlvls)*icolmax/nlvls
level=-25.+5.*findgen(nlvls)
syyyymmdd=strcompress(yyyymmdd,/remove_all)
xindex=where(strmid(syyyymmdd,6,2) eq '15',nxticks)
xlab=strmid(syyyymmdd(xindex),4,2)
contour,nhdt,1.+findgen(kday),lev,/noeras,xrange=[1,kday],yrange=[100.,0.1],charsize=1.25,color=0,/ylog,$
      ytitle='Pressure (hPa)',/fill,c_color=col1,levels=level,charthick=1.5,title='Temp at 85 N minus Temp at 60 N'	;,$
;      xtickname=xlab,xticks=nxticks-1,xtickv=xindex
index=where(level lt 0.)
contour,nhdt,1.+findgen(kday),lev,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_linestyle=5
index=where(level gt 0.)
contour,nhdt,1.+findgen(kday),lev,levels=level(index),color=mcolor,/follow,/overplot,min_value=-9999.
contour,nhdt,1.+findgen(kday),lev,levels=[0.],color=0,thick=5,/follow,/overplot,min_value=-9999.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(delta K)',charsize=1.25,charthick=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
end
