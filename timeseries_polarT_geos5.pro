;
; save daily Temperature at 90 N for GEOS-5 record
; store in one IDL save file
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=600
nydim=600
xorig=[0.20]
yorig=[0.30]
xlen=0.65
ylen=0.4
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
smonth=['J','A','S','O','N','D','J','F','M','A','M','J']
goto,plotit
lstmn=1
lstdy=1
lstyr=2004
ledmn=9
leddy=1
ledyr=2008
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0.
nfile=kday
yyyymmdd=lonarr(nfile)
fdoy=fltarr(nfile)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4)',iyr)
      if strmid(syr,0,1) eq '2' then uyr=string(FORMAT='(I2.2)',iyr-2000)
      if strmid(syr,0,1) eq '1' then uyr=string(FORMAT='(I2.2)',iyr-1900)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      yyyymmdd(icount)=long(sdate)
      fdoy(icount)=float(iday)
      ifile=mon(imn-1)+sdy+'_'+uyr
;
; read GEOS
;
      rd_geos5_nc3_meto,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
               pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if icount eq 0L then begin
         ynindex=where(alat eq max(alat))
         ysindex=where(alat eq min(alat))
; 
; set to maximum temperature poleward of 45 degrees at each theta level
;
         nh_polart=fltarr(nfile,nth)
         sh_polart=fltarr(nfile,nth)
      endif
;
; Temperature=theta*(p/po)^R/cp and divide by 1000 for km
;
    t2=0.*p2
    for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; retain max temperature poleward of 45 degrees
;
    for k=0,nth-1 do begin
        nhtemp=reform(t2(ynindex(0),*,k))
        shtemp=reform(t2(ysindex(0),*,k))
        nh_polart(icount,k)=total(nhtemp)/float(nc)
        sh_polart(icount,k)=total(shtemp)/float(nc)
    endfor
    icount=icount+1.
goto,jump
;
; save file
;
saveit:
save,file='GEOS5_polarT.sav',fdoy,yyyymmdd,th,nh_polart,sh_polart
plotit:
restore,'GEOS5_polarT.sav'
;
; remove leap days and data voids
;
good=where(strmid(strcompress(YYYYMMDD,/remove_all),4,4) ne '0229' and yyyymmdd ne 0)
fdoy=fdoy(good)
yyyymmdd=yyyymmdd(good)
nh_polart=nh_polart(good,*)
sh_polart=sh_polart(good,*)
;
; plot T timeseries at rlev
;
rlev=1000.
;print,th
;read,'Enter level ',rlev
ilev=where(rlev eq th)
ilev=ilev(0)
slev=strcompress(long(rlev),/remove_all)
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='timeseries_polarT_'+slev+'K_geos5.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; extract level
;
nh_polart_lev=reform(nh_polart(*,ilev))
sh_polart_lev=reform(sh_polart(*,ilev))
;
; construct daily averages from multiple years
;
ndoy=365L
nh_polart_mean=fltarr(ndoy)
nh_polart_sigma=fltarr(ndoy)
nh_polart_max=fltarr(ndoy)
nh_polart_min=fltarr(ndoy)
sh_polart_mean=fltarr(ndoy)
sh_polart_sigma=fltarr(ndoy)
sh_polart_max=fltarr(ndoy)
sh_polart_min=fltarr(ndoy)
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
for idoy=0L,ndoy-1L do begin
;
; average does not include 2008
;
    index=where(fdoy eq float(idoy+1L) and nh_polart_lev ne 0. and syear ne '2008',nn)
    if nn ge 2L then begin
       result=moment(nh_polart_lev(index))
       nh_polart_mean(idoy)=result(0)
       nh_polart_sigma(idoy)=sqrt(result(1))
       nh_polart_max(idoy)=max(nh_polart_lev(index))
       nh_polart_min(idoy)=min(nh_polart_lev(index))
if nh_polart_mean(idoy)-nh_polart_sigma(idoy) lt nh_polart_min(idoy) then begin
   print,'mean-sigma lt min by ',min(nh_polart_lev(index))-(nh_polart_mean(idoy)-nh_polart_sigma(idoy)),idoy,result(0),sqrt(result(1)),max(nh_polart_lev(index)),min(nh_polart_lev(index))
;stop
endif
if nh_polart_mean(idoy)+nh_polart_sigma(idoy) gt nh_polart_max(idoy) then begin
   print,'mean+sigma gt max ',(nh_polart_mean(idoy)+nh_polart_sigma(idoy))-nh_polart_max(idoy),idoy,result(0),sqrt(result(1)),max(nh_polart_lev(index)),min(nh_polart_lev(index))
;stop
endif

    endif
    index=where(fdoy eq float(idoy+1L) and sh_polart_lev ne 0.,nn)
    if nn ge 2L then begin
       result=moment(sh_polart_lev(index))
       sh_polart_mean(idoy)=result(0)
       sh_polart_sigma(idoy)=sqrt(result(1))
       sh_polart_max(idoy)=max(sh_polart_lev(index))
       sh_polart_min(idoy)=min(sh_polart_lev(index))
    endif
endfor
;for idoy=1L,ndoy-2L do begin
;    if nh_polart_mean(idoy) lt 170. then begin
;       nh_polart_mean(idoy)=(nh_polart_mean(idoy-1)+nh_polart_mean(idoy+1))/2.
;       nh_polart_sigma(idoy)=(nh_polart_sigma(idoy-1)+nh_polart_sigma(idoy+1))/2.
;       nh_polart_max(idoy)=(nh_polart_max(idoy-1)+nh_polart_max(idoy+1))/2.
;       nh_polart_min(idoy)=(nh_polart_min(idoy-1)+nh_polart_min(idoy+1))/2.
;    endif
;endfor
;
; swap
;
nh_polart_mean_new=0.*nh_polart_mean
nh_polart_sigma_new=0.*nh_polart_mean
nh_polart_max_new=0.*nh_polart_mean
nh_polart_min_new=0.*nh_polart_mean
nh_polart_mean_new(0:182)=nh_polart_mean(182:ndoy-1)
nh_polart_mean_new(183:ndoy-1)=nh_polart_mean(0:181)
nh_polart_sigma_new(0:182)=nh_polart_sigma(182:ndoy-1)
nh_polart_sigma_new(183:ndoy-1)=nh_polart_sigma(0:181)
nh_polart_min_new(0:182)=nh_polart_min(182:ndoy-1)
nh_polart_min_new(183:ndoy-1)=nh_polart_min(0:181)
nh_polart_max_new(0:182)=nh_polart_max(182:ndoy-1)
nh_polart_max_new(183:ndoy-1)=nh_polart_max(0:181)
nh_polart_mean=nh_polart_mean_new
nh_polart_sigma=nh_polart_sigma_new
nh_polart_min=nh_polart_min_new
nh_polart_max=nh_polart_max_new

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=180.
imax=280.
plot,findgen(ndoy),nh_polart_mean,ytitle='Temperature',yrange=[imin,imax],$
     xtickname=smonth,xticks=n_elements(smonth)-1,charsize=1.5,color=0,$
     title='GEOS-5 90N Temperature at '+slev+' K'
for idoy=0L,ndoy-1L do begin
;    plots,idoy,nh_polart_mean(idoy)-nh_polart_sigma(idoy)
;    plots,idoy,nh_polart_mean(idoy)+nh_polart_sigma(idoy),thick=3,/continue,color=200
     plots,idoy,nh_polart_min(idoy)
     plots,idoy,nh_polart_max(idoy),thick=2,/continue,color=200
endfor
oplot,findgen(ndoy),nh_polart_max,color=0
oplot,findgen(ndoy),nh_polart_min,color=0
oplot,findgen(ndoy),nh_polart_mean,thick=4,color=0
;
; superimpose 2008 in red
;
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
index=where(syear eq '2008')
loadct,38
oplot,fdoy(index)+181,NH_POLART_LEV(index),color=250,thick=4
index=where(syear eq '2007' and fdoy ge 181)
oplot,fdoy(index)-181,NH_POLART_LEV(index),color=250,thick=4


if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_polarT_'+slev+'K_geos5.ps -rotate -90 timeseries_polarT_'+slev+'K_geos5.jpg'
endif

end
