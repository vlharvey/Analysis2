;
; altitude-time series of MLS water
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,38
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
xorig=[0.20]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
lstmn=5
lstdy=1
lstyr=2009
ledmn=7
leddy=19
ledyr=2009
lstday=0
ledday=0
;goto,quick
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
      print,sdate 
;
; restore MLS CO on this day
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
      dum=findfile(mdir+'cat_mls_v2.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v2.2_'+sdate+'.sav'
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlscomix=mix
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         mlspolarco_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute polar CO
;
      mlspolarco=fltarr(nlv)
      mlsncoprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
;         if latitude(ii) ge 40. and latitude(ii) le 50. then begin
;         if latitude(ii) ge 50. and latitude(ii) le 60. then begin
          if latitude(ii) ge 60. and latitude(ii) le 70. then begin
             co_prof=reform(mlscomix(ii,*))
             good=where(co_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolarco(good)=mlspolarco(good)+reform(co_prof(good))
                mlsncoprof(good)=mlsncoprof(good)+1L
             endif
          endif
      endfor
      good=where(mlsncoprof gt 0L)
      if good(0) ne -1L then mlspolarco(good)=mlspolarco(good)/float(mlsncoprof(good))
      mlspolarco_zt(icount,*)=mlspolarco
skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; interpolate small gaps in time
;
for k=0,nlv-1 do begin
    dlev=reform(mlspolarco_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 15.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    mlspolarco_zt(*,k)=dlev
endfor
mlspolarco_zt=mlspolarco_zt*1.e6
save,file='zt_mls_h2o.sav',mlspolarco_zt,kday,altitude,sdate_all
quick:
restore,'zt_mls_h2o.sav'
;restore,'zt_mls_h2o_70N.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_mls_h2o_'+syear(0)+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; plot Arctic mean temperature and NO2
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
level=[0.3,0.4,0.5,0.6,0.8,1.,1.5,2.,3.,4.,5.,6.,6.5,7.]
slevel=['0.3','0.4','0.5','0.6','0.8','1','1.5','2','3','4','5','6','6.5','7']
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(mlspolarco_zt eq 0.)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
;mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
index=where(finite(mlspolarco_zt) ne 1)
if index(0) ne -1L then mlspolarco_zt(index)=0.
;mlspolarco_zt(kday-1,*)=mlspolarco_zt(kday-2,*)

contour,mlspolarco_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[50.,87.],$
;     charsize=1.5,color=0,ytitle='Altitude (km)',title='MLS H!l2!nO Avg > 70 N',/cell_fill,c_color=col1,$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='MLS H!l2!nO Avg 60-70 N',/cell_fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
contour,mlspolarco_zt,1.+findgen(kday),altitude,levels=level,color=0,/follow,/overplot,c_labels=1+fltarr(nlvls)
index=where(sdate_all eq '20090713' or sdate_all eq '20090714' or sdate_all eq '20090715')
oplot,index,82.+0.*index,psym=2,color=0
;index=where(mlspolarco_zt eq 0.)
;if index(0) ne -1L then mlspolarco_zt(index)=0./0.
;mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
;index=where(finite(mlspolarco_zt) ne 1)
;if index(0) ne -1L then mlspolarco_zt(index)=0.
;mlspolarco_zt(kday-1,*)=mlspolarco_zt(kday-2,*)
;contour,mlspolarco_zt,1.+findgen(kday),altitude,levels=[0.05,0.08,0.1,0.2,0.5],color=0,/follow,/overplot,$
;        min_value=0.,c_labels=[1,1,1],thick=5
xyouts,xmn+0.02,ymn+0.02,syear(0),/normal,color=0,charsize=3,charthick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(ppmv)',xticks=nlvls-1,xtickname=slevel
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
       spawn,'convert -trim zt_mls_h2o_'+syear(0)+'.ps -rotate -90 zt_mls_h2o_'+syear(0)+'.jpg'
    endif
end
