;
; timeseries of MLS temperature 1 May to 1 July
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
mdir='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
start_year=[2007,2008,2009,2010,2011,2012,2013,2014]
start_date=[-27, -21, -24, -24, -26, -27, -34, -28]
end_date=[66, 65, 61, 61, 64, 61, 64, 80]
nyear=n_elements(start_year)
kcount=0
;goto,quick
for nn=start_year(0),start_year(nyear-1L) do begin

lstmn=5         ; NH
lstdy=10
lstyr=nn
ledmn=7         ; NH
leddy=10
ledyr=nn
lstday=0
ledday=0
;
; loop over years
;
for iyear=lstyr,ledyr do begin
for ilat=50,50 do begin		;80,5 do begin

syr=string(format='(i2.2)',iyear-2000)
slt=strcompress(ilat,/remove_all)
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
if kcount eq 0L then sdate_all=strarr(kday,nyear)
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
      dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
;
; apply mask
;
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
;        print,altitude
         ralt=83.
;        read,'Enter desired altitude ',ralt
         index=where(abs(altitude-ralt) eq min(abs(altitude-ralt)))
         ilev=index(0)
         salt=strcompress(long(ralt),/remove_all)
;        rlat=75.
;        read,'Enter desired latitude',rlat
;        slat=strcompress(long(rlat),/remove_all)
         mlspolartp_xt=fltarr(kday,nyear)
         kcount=1
      endif
;
; compute mean temp around lat circle
;
      index=where(abs(latitude-rlat) le 2.0)
      if index(0) ne -1L then mlspolartp_xt(icount,nn-start_year(0))=mean(temperature(index,ilev))
skipmls:
      icount=icount+1L
goto,jump

next:

endfor  ; loop over latitudes
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
save,file='timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_onset.sav',mlspolartp_xt,kday,sdate_all
quick:
yearlab='2007-2014'
slat='50'
salt='83'
restore,'timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_onset.sav'
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_onset.ps'
   !p.charsize=1.55
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean temperature and CO
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolartp_xt lt 110.)	;eq 0.)
if index(0) ne -1L then mlspolartp_xt(index)=0./0.
plot,130.+findgen(kday)-172.,smooth(mlspolartp_xt(*,0),3,/nan,/edge_truncate),/noeras,xrange=[-40.,20],yrange=[130.,175.],ytitle='MLS Temperature (K)',$
      charsize=2,color=0,xtitle='DFS',charthick=2	;,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.,/nodata
for nn=0,nyear-1L do oplot,130.+findgen(kday)-172.,smooth(mlspolartp_xt(*,nn),3,/nan,/edge_truncate),color=mcolor*(nn/float(nyear)),thick=20

;axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[-6,-51.],/save,charsize=2			; May 1 too June 15
;axis,yaxis=1,xaxis=1,yrange=[150.,190.],xrange=[-11,-42.],/save,charsize=2,color=0,charthick=2,xtitle='DFS'    ; May 10 to July 10
;level=[1,5,10,20,30,40,50,60,70,80,90,100]
;contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=level,color=mcolor,/overplot,/noeras,$
;        c_charsize=2,c_labels=1+0*level,thick=10
;index=where(iyr eq start_year)
;loadct,0
;plots,0.,start_date(index)
;plots,360.,start_date(index),/continue,color=0,thick=20
xyouts,xmn+0.3,ymx-0.05,slat+'N/'+salt+'km',/normal,charsize=2.5,color=0,charthick=2
loadct,39

nlvls=nyear
col1=1+indgen(nlvls)*icolmax/nlvls

imin=long(min(syear))
imax=long(max(syear))
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,charthick=2
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
   spawn,'convert timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_onset.ps -rotate -90 timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_onset.png'
endif
end
