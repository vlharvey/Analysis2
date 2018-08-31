;
; timeseries of MLS water 1 May to 1 July
; loop over latitude bins every 5 degrees
; overplot CIPS frequency
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
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
xorig=[0.2]
yorig=[0.3]
xlen=0.7
ylen=0.5
cbaryoff=0.14
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
start_year=[2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]
start_dfs=[-27, -21, -24, -24, -26, -27, -34, -28,-33,-28,-27,-27]			; 2017 and 2018 DFS of -27 as a placeholder is the mean of the previous years
start_doy=[-27, -21, -24, -24, -26, -27, -34, -28,-33,-28,-27,-27]+172
end_date=[66, 65, 61, 61, 64, 61, 64, 80]
nyear=n_elements(start_year)
for ilat=50,80,5 do begin               ;80,5 do begin

kcount=0
for nn=start_year(0),start_year(nyear-1L) do begin

lstmn=6         ; NH
lstdy=9
lstyr=nn
ledmn=8         ; NH
leddy=12
ledyr=nn
lstday=0
ledday=0
;
; loop over years
;
for iyear=lstyr,ledyr do begin

syr=string(format='(i2.2)',iyear-2000)
slt=strcompress(ilat,/remove_all)
;goto,quick
;restore,'/Volumes/Data/CIPS_data/Pre_process/xt_cips_3c_NH'+syr+'_DSC_'+slt+'Lat_2G_freq.sav
;decfreq=XTFREQ
;restore,'/Volumes/Data/CIPS_data/Pre_process/xt_cips_3c_NH'+syr+'_ASC_'+slt+'Lat_2G_freq.sav
;ascfreq=XTFREQ
;cipsfreq=0.*decfreq
;index=where(decfreq ne 0. and ascfreq eq 0.)
;if index(0) ne -1L then cipsfreq(index)=decfreq(index)
;index=where(decfreq eq 0. and ascfreq ne 0.)
;if index(0) ne -1L then cipsfreq(index)=ascfreq(index)
;index=where(decfreq ne 0. and ascfreq ne 0.)
;if index(0) ne -1L then cipsfreq(index)=(decfreq(index)+ascfreq(index))/2.

rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

;goto,quick

z = stddat(lstmn,lstdy,iyear,lstday)
z = stddat(ledmn,leddy,iyear,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
if kcount eq 0L then begin
   sdate_all=strarr(kday,nyear)
   doy_all=strarr(kday,nyear)
endif
;
; Compute initial Julian date
;
iyr = iyear
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; longitude grid
;
dx=15.
nc=long(360./dx)+1
longrid=dx*findgen(nc)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,next
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      sdate_all(icount,nn-start_year(0))=sdate
      doy_all(icount,nn-start_year(0))=iday
      print,sdate
;
; restore MLS on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v4.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v4.2_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v4.2_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v4.2_'+sdate+'.sav'
;
; apply mask
;
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      mlsh2omix=mix
      good=where(mlsh2omix ne -99.)
      if good(0) ne -1L then mlsh2omix(good)=mlsh2omix(good)*1.e6
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
;        print,altitude
         ralt=83.
         ralt2=30.
;        read,'Enter desired altitude ',ralt
         index=where(abs(altitude-ralt) eq min(abs(altitude-ralt)))
         ilev=index(0)
         index2=where(abs(altitude-ralt2) eq min(abs(altitude-ralt2)))
         ilev2=index2(0)

         salt=strcompress(long(ralt),/remove_all)
         salt2=strcompress(long(ralt2),/remove_all)
;        rlat=65.
;        read,'Enter desired latitude',rlat
;        slat=strcompress(long(rlat),/remove_all)
         mlspolarh2o_xt=fltarr(kday,nyear)
         mlspolarh2o_xt_min=fltarr(kday,nyear)
         mlspolarh2o_xt_max=fltarr(kday,nyear)

         mlspolartp_xt=fltarr(kday,nyear)
         mlspolartp_xt_min=fltarr(kday,nyear)
         mlspolartp_xt_max=fltarr(kday,nyear)

         mlspolartp_strat=fltarr(kday,nyear)

         kcount=1
      endif
;
; compute mean temp around lat circle
;
      index=where(abs(latitude-rlat) le 2.5)
      if index(0) ne -1L then begin
         tdata=reform(temperature(index,ilev))
         tdata2=reform(temperature(index,ilev2))
         good=where(tdata gt 0.)
         if good(0) ne -1L then begin
;if min(tdata(good)) lt 130. then stop
            mlspolartp_xt(icount,nn-start_year(0))=mean(tdata(good))
            mlspolartp_xt_min(icount,nn-start_year(0))=min(tdata(good))
            mlspolartp_xt_max(icount,nn-start_year(0))=max(tdata(good))
         endif
         good=where(tdata2 gt 0.)
         if good(0) ne -1L then mlspolartp_strat(icount,nn-start_year(0))=mean(tdata2(good))

      endif
;
; mean water
;
      index=where(abs(latitude-rlat) le 2.5)
      if index(0) ne -1L then begin
         hdata=reform(mlsh2omix(index,ilev))
         good=where(hdata gt 0.)
         if good(0) ne -1L then begin
            mlspolarh2o_xt(icount,nn-start_year(0))=mean(hdata(good))
            mlspolarh2o_xt_min(icount,nn-start_year(0))=min(hdata(good))
            mlspolarh2o_xt_max(icount,nn-start_year(0))=max(hdata(good))
         endif
      endif

skipmls:
      icount=icount+1L
goto,jump

next:

endfor  ; loop over days
endfor	; loop over years

plotit:
;
; year date label
;
syear=strmid(sdate_all,0,4)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
; save temp
;
save,file='timeseries_mls_h2o_'+yearlab+'_'+slat+'_pmc_midseason.sav',mlspolarh2o_xt,kday,sdate_all,mlspolarh2o_xt_min,mlspolarh2o_xt_max,doy_all
quick:
yearlab='2007-2018'
slat=slt
salt='83'
salt2='30'
restore,'timeseries_mls_h2o_'+yearlab+'_'+slat+'_pmc_midseason.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
;xindex=where(sday eq '01' or sday eq '05' or sday eq '10' or sday eq '15' or sday eq '20' or sday eq '25',nxticks)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_mls_h2o_'+yearlab+'_'+slat+'_pmc_midseason.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean H2O
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolarh2o_xt eq 0.)	;eq 0.)
if index(0) ne -1L then mlspolarh2o_xt(index)=0./0.
index=where(mlspolarh2o_xt_min eq 0.)      ;eq 0.)
if index(0) ne -1L then mlspolarh2o_xt_min(index)=0./0.
index=where(mlspolarh2o_xt_max eq 0.)      ;eq 0.)
if index(0) ne -1L then mlspolarh2o_xt_max(index)=0./0.
imin=min(mlspolarh2o_xt,/nan)
imax=max(mlspolarh2o_xt,/nan)
plot,doy_all(*,0),smooth(mlspolarh2o_xt(*,0),1,/nan,/edge_truncate),/noeras,xrange=[min(doy_all(*,0)),max(doy_all(*,0))],yrange=[imin,imax],ytitle='MLS H2O (ppmv)',$
      charsize=1.5,color=0,xtitle='DOY',charthick=2,/nodata	;,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.,/nodata

for nn=0,nyear-1L do begin
    oplot,doy_all(*,nn),mlspolarh2o_xt(*,nn),color=mcolor*(nn/float(nyear)),thick=20	;psym=nn
    oplot,[start_doy(nn),start_doy(nn)],[0,7],color=mcolor*(nn/float(nyear)),thick=5,linestyle=5
endfor
xyouts,xmx-0.2,ymn+0.01,slat+'N/'+salt+'km',/normal,charsize=2,color=0,charthick=2
loadct,39

nlvls=nyear
col1=1+indgen(nlvls)*icolmax/nlvls

imin=long(min(syear))
imax=long(max(syear))
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1,charthick=2,xticks=nyear,xtickname=['       '+strcompress(start_year,/r),' ']
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nyear)
for jj=0,nyear-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert timeseries_mls_h2o_'+yearlab+'_'+slat+'_pmc_midseason.ps -rotate -90 timeseries_mls_h2o_'+yearlab+'_'+slat+'_pmc_midseason.png'
endif

endfor	; loop over latitudes
end
