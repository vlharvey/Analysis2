;
; save daily avg, min, max Temperature poleward of 45 degrees on each day and altitude
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
; restart from previously saved file?
;
start_scratch='n'
read,'Start from scratch? (y or n)',start_scratch
;
; if restoring previously saved file
;
if start_scratch eq 'n' then begin
;
; contents: 'MetO_polarT.sav',fdoy,yyyymmdd,th,nh_polart_mean,sh_polart_mean,$
;                      nh_polart_max,sh_polart_max,nh_polart_min,sh_polart_min
;
   restore,'Save_files/MetO_polarT.sav'
   fdoy_old=fdoy
   yyyymmdd_old=yyyymmdd
   nh_polart_mean_old=nh_polart_mean
   sh_polart_mean_old=sh_polart_mean
   nh_polart_max_old=nh_polart_max
   sh_polart_max_old=sh_polart_max
   nh_polart_min_old=nh_polart_min
   sh_polart_min_old=sh_polart_min
;
; remove leap days and data voids
;
   good=where(strmid(strcompress(YYYYMMDD,/remove_all),4,4) ne '0229' and yyyymmdd ne 0)
   fdoy=fdoy(good)
   yyyymmdd=yyyymmdd(good)
   nh_polart_mean_old=nh_polart_mean_old(good,*)
   sh_polart_mean_old=sh_polart_mean_old(good,*)
   nh_polart_max_old=nh_polart_max_old(good,*)
   sh_polart_max_old=sh_polart_max_old(good,*)
   nh_polart_min_old=nh_polart_min_old(good,*)
   sh_polart_min_old=sh_polart_min_old(good,*)
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
endif
;
; if starting from scratch
;
if start_scratch eq 'y' then begin
   lstmn=11
   lstdy=1
   lstyr=1991
endif
ledmn=2
leddy=1
ledyr=2010
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
      print,ifile,' ',long(sdate)
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
         p2=fltarr(nr,nc,nth)
;        ynindex=where(alat eq max(alat),nhy)
;        ysindex=where(alat eq min(alat),shy)
         ynindex=where(alat gt  45.,nhy)
         ysindex=where(alat lt -45.,shy)
         kday=ledday-lstday
         nfile=kday
         yyyymmdd=lonarr(nfile)
         fdoy=fltarr(nfile)
         nh_polart_mean=fltarr(nfile,nth)
         sh_polart_mean=fltarr(nfile,nth)
         nh_polart_max=fltarr(nfile,nth)
         sh_polart_max=fltarr(nfile,nth)
         nh_polart_min=fltarr(nfile,nth)
         sh_polart_min=fltarr(nfile,nth)
    endif
    ncdf_varget,ncid,4,p2
    ncdf_close,ncid
;
; Temperature=theta*(p/po)^R/cp and divide by 1000 for km
;
    t2=0.*p2
    for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; retain mean, min, max temperature at the poles
;
    for k=0,nth-1 do begin
        nhtemp=reform(t2(ynindex,*,k))
        shtemp=reform(t2(ysindex,*,k))
        index=where(nhtemp ne 0.,nhy)
        if index(0) ne -1L then begin
           nh_polart_mean(icount,k)=total(nhtemp(index))/float(nhy)
           nh_polart_max(icount,k)=max(nhtemp(index))
           nh_polart_min(icount,k)=min(nhtemp(index))
        endif
        index=where(shtemp ne 0.,shy)
        if index(0) ne -1L then begin
           sh_polart_mean(icount,k)=total(shtemp(index))/float(shy)
           sh_polart_max(icount,k)=max(shtemp(index))
           sh_polart_min(icount,k)=min(shtemp(index))
        endif
    endfor
    yyyymmdd(icount)=long(sdate)
    fdoy(icount)=float(iday)
    icount=icount+1.
goto,jump
saveit:
;
; if NOT started from scratch concatenate old and new
;
if start_scratch eq 'n' then begin
   fdoy=[fdoy_old,fdoy]
   yyyymmdd=[yyyymmdd_old,yyyymmdd]
   nh_polart_mean=[nh_polart_mean_old,nh_polart_mean]
   sh_polart_mean=[sh_polart_mean_old,sh_polart_mean]
   nh_polart_max=[nh_polart_max_old,nh_polart_max]
   sh_polart_max=[sh_polart_max_old,sh_polart_max]
   nh_polart_min=[nh_polart_min_old,nh_polart_min]
   sh_polart_min=[sh_polart_min_old,sh_polart_min]
endif
;
; save file
;
save,file='Save_files/MetO_polarT.sav',fdoy,yyyymmdd,th,nh_polart_mean,sh_polart_mean,$
     nh_polart_max,sh_polart_max,nh_polart_min,sh_polart_min
plotit:
restore,'Save_files/MetO_polarT.sav'
nth=n_elements(th)
;
; remove leap days and data voids
;
good=where(strmid(strcompress(YYYYMMDD,/remove_all),4,4) ne '0229' and yyyymmdd ne 0,kday)
if good(0) ne -1L then begin
   fdoy=fdoy(good)
   yyyymmdd=yyyymmdd(good)
   nh_polart_mean=nh_polart_mean(good,*)
   sh_polart_mean=sh_polart_mean(good,*)
   nh_polart_max=nh_polart_max(good,*)
   sh_polart_max=sh_polart_max(good,*)
   nh_polart_min=nh_polart_min(good,*)
   sh_polart_min=sh_polart_min(good,*)
endif
nh_polart_mean_orig=nh_polart_mean
sh_polart_mean_orig=sh_polart_mean
nh_polart_max_orig=nh_polart_max
sh_polart_max_orig=sh_polart_max
nh_polart_min_orig=nh_polart_min
sh_polart_min_orig=sh_polart_min
;
; plot T timeseries at rlev
;
for kk=0L,nth-1L do begin
rlev=1000.
rlev=th(kk)
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_polarT_'+slev+'K.ps'
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
nh_polart_lev_mean=reform(nh_polart_mean_orig(*,ilev))
sh_polart_lev_mean=reform(sh_polart_mean_orig(*,ilev))
nh_polart_lev_max=reform(nh_polart_max_orig(*,ilev))
sh_polart_lev_max=reform(sh_polart_max_orig(*,ilev))
nh_polart_lev_min=reform(nh_polart_min_orig(*,ilev))
sh_polart_lev_min=reform(sh_polart_min_orig(*,ilev))
;
; interpolate small gaps in time
;
dlevmean=nh_polart_lev_mean
dlevmax=nh_polart_lev_max
dlevmin=nh_polart_lev_min
for i=1,kday-1 do begin
    if dlevmean(i) eq 0. and dlevmean(i-1) ne 0. then begin
       for ii=i+1,kday-1 do begin
           naway=float(ii-i)
           if naway le 5.0 and dlevmean(ii) ne 0. then begin
              dlevmean(i)=(naway*dlevmean(i-1)+dlevmean(ii))/(naway+1.0)
              goto,jump1
           endif
       endfor
jump1:
    endif
    if dlevmax(i) eq 0. and dlevmax(i-1) ne 0. then begin
       for ii=i+1,kday-1 do begin
           naway=float(ii-i)
           if naway le 5.0 and dlevmax(ii) ne 0. then begin
              dlevmax(i)=(naway*dlevmax(i-1)+dlevmax(ii))/(naway+1.0)
              goto,jump2
           endif
       endfor
jump2:
    endif
    if dlevmin(i) eq 0. and dlevmin(i-1) ne 0. then begin
       for ii=i+1,kday-1 do begin
           naway=float(ii-i)
           if naway le 5.0 and dlevmin(ii) ne 0. then begin
              dlevmin(i)=(naway*dlevmin(i-1)+dlevmin(ii))/(naway+1.0)
              goto,jump3
           endif
       endfor
jump3:
    endif
endfor
nh_polart_lev_mean=dlevmean
nh_polart_lev_max=dlevmax
nh_polart_lev_min=dlevmin
;
; construct daily means, sigmas of means, maxs, mins from multiple years
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
    index=where(fdoy eq float(idoy+1L) and nh_polart_lev_mean gt 180. and syear ne '2009',nn)
    if nn ge 2L then begin
       result=moment(nh_polart_lev_mean(index))
       nh_polart_mean(idoy)=result(0)				; mean of the means
       nh_polart_sigma(idoy)=sqrt(result(1))			; sigma of the means
       nh_polart_max(idoy)=max(nh_polart_lev_max(index))	; max of the maxs
       nh_polart_min(idoy)=min(nh_polart_lev_min(index))	; min of the mins
if nh_polart_min(idoy) lt 180. then stop
    endif
    index=where(fdoy eq float(idoy+1L) and sh_polart_lev_mean gt 180.,nn)
    if nn ge 2L then begin
       result=moment(sh_polart_lev_mean(index))
       sh_polart_mean(idoy)=result(0)
       sh_polart_sigma(idoy)=sqrt(result(1))
       sh_polart_max(idoy)=max(sh_polart_lev_max(index))
       sh_polart_min(idoy)=min(sh_polart_lev_min(index))
    endif
endfor
for idoy=1L,ndoy-2L do begin
    if nh_polart_mean(idoy) lt 170. then begin
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
loadct,0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=180.
imax=300.
plot,findgen(ndoy),nh_polart_mean,ytitle='Temperature',yrange=[imin,imax],$
     xtickname=smonth,xticks=n_elements(smonth)-1,charsize=1.5,color=0,$
     title='MetO Mean Polar Temperature'
for idoy=0L,ndoy-1L do begin
    plots,idoy,nh_polart_mean(idoy)-nh_polart_sigma(idoy)
    plots,idoy,nh_polart_mean(idoy)+nh_polart_sigma(idoy),thick=5,/continue,color=200
endfor
oplot,findgen(ndoy),nh_polart_max,color=0
oplot,findgen(ndoy),nh_polart_min,color=0
oplot,findgen(ndoy),nh_polart_mean,thick=5,color=0
;
; superimpose 2008 in green
;
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
index=where(syear eq '2008' and fdoy lt 181)
loadct,39
;oplot,fdoy(index)+181,NH_POLART_LEV_MIN(index),color=100,thick=5
;oplot,fdoy(index)+181,NH_POLART_LEV_MAX(index),color=100,thick=5
oplot,fdoy(index)+181,NH_POLART_LEV_MEAN(index),color=100,thick=10
index=where(syear eq '2007' and fdoy ge 181)
;oplot,fdoy(index)-181,NH_POLART_LEV_MIN(index),color=100,thick=5
;oplot,fdoy(index)-181,NH_POLART_LEV_MAX(index),color=100,thick=5
oplot,fdoy(index)-181,NH_POLART_LEV_MEAN(index),color=100,thick=10
xyouts,10.,205.,'2007-2008',color=100,/data,charsize=1.5,charthick=5
;
; superimpose 2009 in red
;
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
index=where(syear eq '2009' and fdoy lt 181)
;oplot,fdoy(index)+181,NH_POLART_LEV_MIN(index),color=50,thick=5
;oplot,fdoy(index)+181,NH_POLART_LEV_MAX(index),color=250,thick=5
oplot,fdoy(index)+181,NH_POLART_LEV_MEAN(index),color=250,thick=10
index=where(syear eq '2008' and fdoy ge 181)
;oplot,fdoy(index)-181,NH_POLART_LEV_MIN(index),color=50,thick=5
;oplot,fdoy(index)-181,NH_POLART_LEV_MAX(index),color=250,thick=5
oplot,fdoy(index)-181,NH_POLART_LEV_MEAN(index),color=250,thick=10
xyouts,10.,195.,'2008-2009',color=250,/data,charsize=1.5,charthick=5
;
; superimpose 2010 in blue
;
syear=strmid(strcompress(YYYYMMDD,/remove_all),0,4)
index=where(syear eq '2010' and fdoy lt 181)
;oplot,fdoy(index)+181,NH_POLART_LEV_MIN(index),color=50,thick=5
;oplot,fdoy(index)+181,NH_POLART_LEV_MAX(index),color=50,thick=5
oplot,fdoy(index)+181,NH_POLART_LEV_MEAN(index),color=50,thick=10
index=where(syear eq '2009' and fdoy ge 181)
;oplot,fdoy(index)-181,NH_POLART_LEV_MIN(index),color=50,thick=5
;oplot,fdoy(index)-181,NH_POLART_LEV_MAX(index),color=50,thick=5
oplot,fdoy(index)-181,NH_POLART_LEV_MEAN(index),color=50,thick=10
xyouts,10.,185.,'2009-2010',color=50,/data,charsize=1.5,charthick=5

xyouts,10.,290.,slev+' K',color=0,/data,charsize=1.5,charthick=5
;
; print end date
;
maxdate=max(yyyymmdd)
smaxdate=strcompress(maxdate,/remove_all)
datelab=strmid(smaxdate,4,2)+'/'+strmid(smaxdate,6,2)
;xyouts,200.,185.,'Last Day '+datelab,color=250,/data,charsize=1.5,charthick=5

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_polarT_'+slev+'K.ps -rotate -90 timeseries_polarT_'+slev+'K.jpg'
endif

endfor	; loop over levels
end
