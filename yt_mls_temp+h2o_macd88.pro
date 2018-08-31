;
; altitude-time series of MLS temperature + H2O
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
dir='/Volumes/earth/aura6/data/MLS_data/Datfiles/'
lstmn=1
lstdy=1
lstyr=2004
ledmn=12
leddy=31
ledyr=2004
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
;goto,quick

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
      spawn,'ls '+dir+'MLS-Aura_L2GP-H2O_v03-*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',h2ofiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-Temperature_v03-*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',tpfiles
      result5=size(h2ofiles)
      result10=size(tpfiles)
;
; this logic will jump day if any one of the above products are missing
;
      if result5(0) eq 0L or result10(0) eq 0L then begin
         print,'MLS data missing on '+sdate
         goto,jump
      endif
      print,'MLS data complete on '+sdate
;     print,h2ofiles,tpfiles
;
; read EOS-MLS data today.  Original data is in the form of structures
;
      h2o=readl2gp_std(h2ofiles(0),swathName='H2O',variableName=variableName,precisionName=precisionName)
      tp=readl2gp_std(tpfiles(0),swathName='Temperature',variableName=variableName,precisionName=precisionName)
;
; extract some catalog information from MLS water structure... dimensions,x,y,t,p
;
      mprof=h2o.ntimes
      pmls=tp.pressure         ; 55 elements
      latitude=h2o.latitude
;
; build mask
;
      h2omls=transpose(h2o.L2GPVALUE)
      h2oprecision=transpose(h2o.L2GPPRECISION)
      h2ostatus=h2o.STATUS
      h2oquality=h2o.QUALITY
      h2oconvergence=h2o.CONVERGENCE
      h2omask=0.*h2omls
      tpmls=transpose(tp.L2GPVALUE)
      tpprecision=transpose(tp.L2GPPRECISION)
      tpstatus=tp.STATUS
      tpquality=tp.QUALITY
      tpconvergence=tp.CONVERGENCE
      tpmask=0.*tpmls

      h2obad=where(h2oprecision lt 0.)
      if h2obad(0) ne -1L then h2omask(h2obad)=-99.
      tpbad=where(tpprecision lt 0.)
      if tpbad(0) ne -1L then tpmask(tpbad)=-99.

      h2obad=where(h2ostatus mod 2 ne 0L)
      if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpstatus mod 2 ne 0L)
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.

      h2obad=where(h2oquality lt 1.3)                   ; do not use if h2oquality < 1.3
      if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpquality lt 0.65)                    ; do not use if tpquality < 0.65
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.

      h2obad=where(h2oconvergence gt 2.)                   ; do not use if convergence > 2.0
      if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; apply mask
;
      index=where(h2omask eq -99.)
      if index(0) ne -1L then h2omls(index)=-99.
      index=where(tpmask eq -99.)
      if index(0) ne -1L then tpmls(index)=-99.
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         mlsh2o_yt=fltarr(kday,nr)
         mlstp_yt=fltarr(kday,nr)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute zonal mean temp and H2O
;
      plev=0.00464159
      slev=string(plev,format='(f5.3)')
      ilev=where(abs(pmls-plev) eq min(abs(pmls-plev)))
      ilev=ilev(0)
      mlsh2olev=reform(h2omls(*,ilev))
      mlstplev=reform(tpmls(*,ilev))
      for j=1L,nr-2L do begin
          alatm1=latbin(j)-((latbin(1)-latbin(0))/2.)
          alatp1=latbin(j)+((latbin(1)-latbin(0))/2.)
          good=where(latitude ge alatm1 and latitude lt alatp1 and mlsh2olev ne -99.,nn)
          if good(0) ne -1L then mlsh2o_yt(icount,j)=mean(mlsh2olev(good))
          good=where(latitude ge alatm1 and latitude lt alatp1 and mlstplev ne -99.,nn)
          if good(0) ne -1L then mlstp_yt(icount,j)=mean(mlstplev(good))
;print,latbin(j),mlstp_yt(icount,j),mlsh2o_yt(icount,j)
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
;    dlev=reform(mlsh2o_yt(*,k))
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
;    mlsh2o_yt(*,k)=dlev
;    dlev=reform(mlstp_yt(*,k))
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
;    mlstp_yt(*,k)=dlev
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
mlsh2o_yt=mlsh2o_yt*1.e6
;
; save temp, h2o, etc
;
h2o=mlsh2o_yt
temp=mlstp_yt
sdates=sdate_all
save,file='yt_mls_temp_h2o_pmc_'+yearlab+'.sav',h2o,temp,latbin,sdates
quick:
restore,'yt_mls_temp_h2o_pmc_'+yearlab+'.sav'
mlsh2o_yt=h2o
mlstp_yt=temp
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
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yt_mls_temp+h2o_'+yearlab+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot latitude-time temperature and water
;
erase
xyouts,0.3,0.95,yearlab+' MLS '+slev+' hPa',/normal,color=0,charsize=2
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(mlsh2o_yt eq 0.)
if index(0) ne -1L then mlsh2o_yt(index)=0./0.
mlsh2o_yt=smooth(mlsh2o_yt,3,/NaN)
index=where(mlstp_yt eq 0.)
if index(0) ne -1L then mlstp_yt(index)=0./0.
mlstp_yt=smooth(mlstp_yt,3,/NaN)
tlevel=130.+5.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlstp_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='Temperature',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlstp_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
contour,mlstp_yt,1.+findgen(kday),latbin,levels=[160.],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
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
tlevel=0.25*findgen(24)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlsh2o_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='Water Vapor',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=[4.],color=mcolor,/follow,/overplot,thick=5
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
       spawn,'convert -trim yt_mls_temp+h2o_'+yearlab+'.ps -rotate -90 yt_mls_temp+h2o_'+yearlab+'.jpg'
    endif
end
