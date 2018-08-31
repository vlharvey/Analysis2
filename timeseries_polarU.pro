;
; restore daily maximum zonal wind speed poleward of 45 degrees for MetO record
; update to current and plot
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nc3

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
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','A','S','O','N','D','J','F','M','A','M','J']
;
; restore save file
;
; contents: 'MetO_polarU.sav',fdoy,yyyymmdd,th,nh_polart,sh_polart
restore,'MetO_polarU.sav'
;
; remove leap days and data voids
;
good=where(strmid(strcompress(YYYYMMDD,/remove_all),4,4) ne '0229' and yyyymmdd ne 0)
fdoy=fdoy(good)
yyyymmdd=yyyymmdd(good)
nh_polart=nh_polart(good,*)
sh_polart=sh_polart(good,*)
;
; read MetO data since last day
;
maxdate=max(yyyymmdd)
smaxdate=strcompress(maxdate,/remove_all)
print,'current record through '+smaxdate
;
; start with maxdate
;
lstmn=long(strmid(smaxdate,4,2))
lstdy=long(strmid(smaxdate,6,2))
lstyr=long(strmid(smaxdate,0,4))
ledmn=3
leddy=9
ledyr=2009
if lstmn eq ledmn and lstdy eq leddy and lstyr eq ledyr then goto,plotit
lstday=0
ledday=0
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
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
icount=0.

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
      ifile=mon(imn-1)+sdy+'_'+uyr
      print,ifile,' ',long(sdate),float(iday)
      dum1=findfile(diru+ifile+'.nc3')
      if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
      if dum1(0) eq '' then goto,jump
      if icount eq 0 then begin
         nr=0L
         nc=0L
         nth=0L
         ncdf_diminq,ncid,0,name,nr
         ncdf_diminq,ncid,1,name,nc
         ncdf_diminq,ncid,2,name,nth
         alon=fltarr(nc)
         alat=fltarr(nr)
         th=fltarr(nth)
         ncdf_varget,ncid,0,alon
         ncdf_varget,ncid,1,alat
         ncdf_varget,ncid,2,th
         u2=fltarr(nr,nc,nth)
         v2=fltarr(nr,nc,nth)
         ynindex=where(alat gt 45.)
         ysindex=where(alat lt 45.)
; 
; set to maximum wind speed poleward of 45 degrees at each theta level
;
         kday=ledday-lstday
         nfile=kday
         yyyymmdd_new=lonarr(nfile)
         fdoy_new=fltarr(nfile)
         nh_polart_new=fltarr(nfile,nth)
         sh_polart_new=fltarr(nfile,nth)
    endif
    ncdf_varget,ncid,6,u2
    ncdf_varget,ncid,7,v2
    ncdf_close,ncid
;
; Temperature=theta*(p/po)^R/cp and divide by 1000 for km
;
    s2=0.*u2
    s2=sqrt(u2^2.+v2^2.)
;
; retain max temperature poleward of 45 degrees
;
    for k=0,nth-1 do begin
        nhtemp=reform(s2(ynindex(0),*,k))
        shtemp=reform(s2(ysindex(0),*,k))
        nh_polart_new(icount,k)=total(nhtemp)/float(nc)
        sh_polart_new(icount,k)=total(shtemp)/float(nc)
    endfor
    yyyymmdd_new(icount)=long(sdate)
    fdoy_new(icount)=float(iday)
    icount=icount+1.
goto,jump
saveit:
;
; concatenate old and new
;
fdoy=[fdoy,fdoy_new]
yyyymmdd=[yyyymmdd,yyyymmdd_new]
nh_polart=[nh_polart,nh_polart_new]
sh_polart=[sh_polart,sh_polart_new]
;
; comment out if not starting from scratch
;
;fdoy=fdoy_new
;yyyymmdd=yyyymmdd_new
;nh_polart=nh_polart_new
;sh_polart=sh_polart_new
;
; save file
;
save,file='MetO_polarU.sav',fdoy,yyyymmdd,th,nh_polart,sh_polart
plotit:
restore,'MetO_polarU.sav'
;
; remove leap days and data voids
;
good=where(strmid(strcompress(YYYYMMDD,/remove_all),4,4) ne '0229' and yyyymmdd ne 0)
if good(0) ne -1L then begin
   fdoy=fdoy(good)
   yyyymmdd=yyyymmdd(good)
   nh_polart=nh_polart(good,*)
   sh_polart=sh_polart(good,*)
endif
;
; plot speed timeseries at rlev
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
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_polarU_'+slev+'K.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
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
    index=where(fdoy eq float(idoy+1L) and nh_polart_lev ne 0. and syear ne '2009',nn)
    if nn ge 2L then begin
       result=moment(nh_polart_lev(index))
       nh_polart_mean(idoy)=result(0)
       nh_polart_sigma(idoy)=sqrt(result(1))
       nh_polart_max(idoy)=max(nh_polart_lev(index))
       nh_polart_min(idoy)=min(nh_polart_lev(index))
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
for idoy=1L,ndoy-2L do begin
    if nh_polart_mean(idoy) lt 0. then begin
       nh_polart_mean(idoy)=(nh_polart_mean(idoy-1)+nh_polart_mean(idoy+1))/2.
       nh_polart_sigma(idoy)=(nh_polart_sigma(idoy-1)+nh_polart_sigma(idoy+1))/2.
       nh_polart_max(idoy)=(nh_polart_max(idoy-1)+nh_polart_max(idoy+1))/2.
       nh_polart_min(idoy)=(nh_polart_min(idoy-1)+nh_polart_min(idoy+1))/2.
    endif
endfor
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
imin=0.
imax=100.
plot,findgen(ndoy),nh_polart_mean,ytitle=slev+' K Wind Speed (m/s)',yrange=[imin,imax],$
     xtickname=smonth,xticks=n_elements(smonth)-1,charsize=1.5,color=0,$
     title='MetO Max Wind Speed poleward of 45 N'
for idoy=0L,ndoy-1L do begin
    plots,idoy,nh_polart_mean(idoy)-nh_polart_sigma(idoy)
    plots,idoy,nh_polart_mean(idoy)+nh_polart_sigma(idoy),thick=5,/continue,color=200
endfor
oplot,findgen(ndoy),nh_polart_max,color=0
oplot,findgen(ndoy),nh_polart_min,color=0
oplot,findgen(ndoy),nh_polart_mean,thick=4,color=0
;
; superimpose 2009 in red
;
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
index=where(syear eq '2009')
loadct,38
oplot,fdoy(index)+181,NH_POLART_LEV(index),color=250,thick=5
index=where(syear eq '2008' and fdoy ge 181)
oplot,fdoy(index)-181,NH_POLART_LEV(index),color=250,thick=5
xyouts,10.,imin+5.,'2008-2009',color=250,/data,charsize=1.5,charthick=5
;
; print end date
;
maxdate=max(yyyymmdd)
smaxdate=strcompress(maxdate,/remove_all)
datelab=strmid(smaxdate,4,2)+'/'+strmid(smaxdate,6,2)
xyouts,200.,imin+5.,'Last Day '+datelab,color=250,/data,charsize=1.5,charthick=5

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_polarU_'+slev+'K.ps -rotate -90 timeseries_polarU_'+slev+'K.jpg'
endif
end
