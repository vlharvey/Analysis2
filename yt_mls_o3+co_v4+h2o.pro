;
; version 4
; latitude-time series of MLS O3 + CO
;
; note, CO and Ozone have different number of altitude levels
;
@stddat
@kgmt
@ckday
@kdate
@readl2gp_std

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.8
ylen=0.275
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
dir='/atmos/aura6/data/MLS_data/Datfiles/'
;
; loop over years
;
for iyear=2004,2016 do begin

lstmn=1
lstdy=1
lstyr=iyear
ledmn=12
leddy=31
ledyr=iyear
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
minyear=lstyr
maxyear=ledyr
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nr=91L
latbin=-90.+2.*findgen(nr)

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L

goto,quick
;
; Compute initial Julian date
;
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
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; look for EOS-MLS data files for today
;
      spawn,'ls '+dir+'MLS-Aura_L2GP-CO_v04-*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',cofiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-O3_v04-*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',o3files
      result5=size(cofiles)
      result10=size(o3files)
;
; this logic will jump day if any one of the above products are missing
;
      if result5(0) eq 0L or result10(0) eq 0L then begin
         print,'MLS data missing on '+sdate
         goto,skipmls
      endif
      print,'MLS data complete on '+sdate
;     print,cofiles,o3files
;
; read EOS-MLS data today.  Original data is in the form of structures
;
      co=readl2gp_std(cofiles(0),swathName='CO',variableName=variableName,precisionName=precisionName)
      o3=readl2gp_std(o3files(0),swathName='O3',variableName=variableName,precisionName=precisionName)
;
; extract some catalog information from MLS structure... dimensions,x,y,t,p
;
      mprof=co.ntimes
      pmls=o3.pressure         ; 55 elements
      pmls2=co.pressure         ; 37 elements
      latitude=co.latitude
;
; build mask
;
      comls=transpose(co.L2GPVALUE)
      coprecision=transpose(co.L2GPPRECISION)
      costatus=co.STATUS
      coquality=co.QUALITY
      coconvergence=co.CONVERGENCE
      comask=0.*comls
      o3mls=transpose(o3.L2GPVALUE)
      o3precision=transpose(o3.L2GPPRECISION)
      o3status=o3.STATUS
      o3quality=o3.QUALITY
      o3convergence=o3.CONVERGENCE
      o3mask=0.*o3mls

      cobad=where(coprecision lt 0.)
      if cobad(0) ne -1L then comask(cobad)=-99.
      o3bad=where(o3precision lt 0.)
      if o3bad(0) ne -1L then o3mask(o3bad)=-99.

      cobad=where(costatus mod 2 ne 0L)
      if cobad(0) ne -1L then comask(cobad,*)=-99.
      o3bad=where(o3status mod 2 ne 0L)
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.

      cobad=where(coquality lt 1.5)                   ; do not use if coquality < 1.5
      if cobad(0) ne -1L then comask(cobad,*)=-99.
      o3bad=where(o3quality lt 1.)                    ; do not use if o3quality < 1
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.

      cobad=where(coconvergence gt 1.03)                   ; do not use if convergence > 1.03
      if cobad(0) ne -1L then comask(cobad,*)=-99.
      o3bad=where(o3convergence gt 1.03)                     ; do not use if convergence > 1.03
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
;
; apply mask
;
      index=where(comask eq -99.)
      if index(0) ne -1L then comls(index)=-99.
      index=where(o3mask eq -99.)
      if index(0) ne -1L then o3mls(index)=-99.
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         mlsco_yt=fltarr(kday,nr)
         mlso3_yt=fltarr(kday,nr)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute zonal mean 
;
;     plev=0.00464159
      plev=0.01
      slev=string(plev,format='(f5.3)')
      ilev=where(abs(pmls-plev) eq min(abs(pmls-plev)))
      ilev=ilev(0)
      ilev2=where(abs(pmls2-plev) eq min(abs(pmls2-plev)))
      ilev2=ilev2(0)
      mlscolev=reform(comls(*,ilev2))
      mlso3lev=reform(o3mls(*,ilev))
      for j=1L,nr-2L do begin
          alatm1=latbin(j)-((latbin(1)-latbin(0))/2.)
          alatp1=latbin(j)+((latbin(1)-latbin(0))/2.)
          good=where(latitude ge alatm1 and latitude lt alatp1 and mlscolev ne -99.,nn)
          if good(0) ne -1L then mlsco_yt(icount,j)=mean(mlscolev(good))
          good=where(latitude ge alatm1 and latitude lt alatp1 and mlso3lev ne -99.,nn)
          if good(0) ne -1L then mlso3_yt(icount,j)=mean(mlso3lev(good))
;print,latbin(j),mlso3_yt(icount,j),mlsco_yt(icount,j)
      endfor
;stop
skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; interpolate small gaps in time
;
;for k=0,nr-1 do begin
;    dlev=reform(mlsco_yt(*,k))
;    for i=1,kday-1 do begin
;        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
;           for ii=i+1,kday-1 do begin
;               naway=float(ii-i)
;               if naway le 5.0 and dlev(ii) ne 0. then begin
;                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
;                  goto,jump1
;               endif
;           endfor
;jump1:
;        endif
;    endfor
;    mlsco_yt(*,k)=dlev
;    dlev=reform(mlso3_yt(*,k))
;    for i=1,kday-1 do begin
;        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
;           for ii=i+1,kday-1 do begin
;               naway=float(ii-i)
;               if naway le 5.0 and dlev(ii) ne 0. then begin
;                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
;                  goto,jump2
;               endif
;           endfor
;jump2:
;        endif
;    endfor
;    mlso3_yt(*,k)=dlev
;endfor
;
; year date label
;
syear=strmid(sdate_all,0,4)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
mlsco_yt=mlsco_yt*1.e6
mlso3_yt=mlso3_yt*1.e6
;
; save
;
co=mlsco_yt
o3=mlso3_yt
sdates=sdate_all
save,file='yt_mls_o3_co_pmc_'+yearlab+'.sav',co,o3,latbin,sdates,plev,slev
quick:
restore,'yt_mls_o3_co_pmc_'+yearlab+'.sav'
mlsco_yt=co
mlso3_yt=o3
sdate_all=sdates
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)	;+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
restore,'yt_mls_temp_h2o_pmc_'+yearlab+'.sav'
mlsh2o_yt=h2o
mlstp_yt=temp
index=where(mlsh2o_yt eq 0.)
if index(0) ne -1L then mlsh2o_yt(index)=0./0.
mlsh2o_yt=smooth(mlsh2o_yt,3,/NaN)
index=where(mlstp_yt eq 0.)
if index(0) ne -1L then mlstp_yt(index)=0./0.
mlstp_yt=smooth(mlstp_yt,3,/NaN)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yt_mls_o3+co_'+yearlab+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot latitude-time 
;
erase
xyouts,0.3,0.95,yearlab+' MLS '+slev+' hPa',/normal,color=0,charsize=2
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(mlsco_yt eq 0.)
if index(0) ne -1L then mlsco_yt(index)=0./0.
mlsco_yt=smooth(mlsco_yt,3,/NaN)
index=where(mlso3_yt eq 0.)
if index(0) ne -1L then mlso3_yt(index)=0./0.
mlso3_yt=smooth(mlso3_yt,3,/NaN)
tlevel=.05*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlso3_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='Ozone',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
;contour,mlso3_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(ppmv)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=0.5*findgen(31)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlsco_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='CO',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlsco_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
;contour,mlsco_yt,1.+findgen(kday),latbin,levels=[4.],color=mcolor,/follow,/overplot,thick=5
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=[0.5,0.75,1.],color=mcolor,/follow,/overplot,thick=3
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=[2.,2.5],color=9.*mcolor,/follow,/overplot,thick=3
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(ppmv)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim yt_mls_o3+co_'+yearlab+'.ps -rotate -90 yt_mls_o3+co_'+yearlab+'.jpg'
    endif

endfor	; loop over years
end
