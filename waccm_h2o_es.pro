;
; save H2O for each ES event
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=9L & lstdy=1L & lstyr=2009L 
ledmn=5L & leddy=1L & ledyr=2010L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)

lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,39
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
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
diro='/Users/harvey/Harvey_etal_2014/Post_process/'
;restore, diro+'WACCM_ES_daily_max_T_Z.sav'
dir='/Volumes/Data/WACCM/WACCM4/mee00fpl_FW2/mee00fpl_FW2.cam2.h3.Year'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
esdates=[20041218,20060130,20151221,20231223,20261210,20320204,20331220,20390226,20420104]
esdates=strcompress(esdates,/remove_all)
;result=size(MAXHEIGHTTHETA)
;nevents=result(1)
nevents=n_elements(esdates)

for iES = 0L, nevents - 1L do begin
    sevent=strtrim(strcompress(string(format='(I3.2)',ies+1)),2)
    esdate0=esdates(ies)
    iyr=long(strmid(esdate0,0,4))
    imn=long(strmid(esdate0,4,2))
    idy=long(strmid(esdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
koff=30
    jday0=jday-koff
    jday1=jday+koff
    CALDAT, jday0, lstmn ,lstdy , lstyr
    CALDAT, jday1, ledmn ,leddy , ledyr

lstday=0L & ledday=0L
if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
;goto,quick
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=long(ledday-lstday+1L)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
if smn eq '02' and sdy eq '29' then goto,jump
print,sdate
;
; read daily file
;
print,'reading ',strmid(sdate,2,2),smn+sdy
        ncfile0=dir+strmid(sdate,2,2)+'_'+smn+sdy+'_H2O.sav'
        restore,ncfile0

      if icount eq 0 then begin
         nth=n_elements(lev)
         nr=n_elements(lat)
         sfile=strarr(kday)
         h2o_zt_nc4=fltarr(kday,nr,nth)
         kcount=1L
      endif
      sfile(icount)=sdate
;
      h2o_zt_nc4(icount,*,*)=h2o
;
      jumpday:
      icount=icount+1L
goto,jump

saveit:
;
; plot altitude-time series 
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
y1='20'+string(FORMAT='(I2.2)',min(yy(index)))
y2='20'+string(FORMAT='(I2.2)',max(yy(index)))
;
; save file
;
save,file=diro+'vortex_h2o_ES_'+esdate0+'_waccm.sav',h2o_zt_nc4,lat,lev,sfile,y1,y2,koff
index=where(h2o_zt_nc4 le 0.)
if index(0) ne -1L then h2o_zt_nc4(index)=0./0.

plotit:
restore,diro+'vortex_h2o_ES_'+esdate0+'_waccm.sav'
if ies eq 0L then begin
   rlat=0.
   print,lat
   read,'Enter latitude ',rlat
   index=where(lat eq rlat)
   ilat=index(0)
   slat=strcompress(long(rlat),/remove_all)
endif
kday=n_elements(SFILE)
h2o_zt_nc4_save=h2o_zt_nc4
h2o_zt_nc4_plot=reform(h2o_zt_nc4(*,ilat,*))
h2o_zt_nc4=h2o_zt_nc4_plot
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/waccm_vortex_h2o_ES_event_'+sevent+'_'+slat+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot zt 
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
level=[0.001,0.01,0.1,.2,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8]
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(h2o_zt_nc4 eq 0.)
if index(0) ne -1L then h2o_zt_nc4(index)=0./0.
h2o_zt_nc4=1.e6*smooth(h2o_zt_nc4,5,/nan)
contour,h2o_zt_nc4,-1.*koff+findgen(kday),lev,color=0,xtitle='Days From ES Onset',thick=6,yrange=[10.,0.0001],/ylog,/noeras,ytitle='Pressure (hPa)',/cell_fill,$
     title='WACCM ES event '+sevent+' ('+esdate0+') at '+slat,levels=level,c_color=col1
contour,h2o_zt_nc4,-1.*koff+findgen(kday),lev,levels=level,c_color=mcolor,/follow,/overplot,c_labels=1+0*level,min_value=0
plots,0,500
plots,0,8000.,/continue,color=mcolor,thick=5
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='H2O (ppmv)',charthick=2,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ../Figures/waccm_vortex_h2o_ES_event_'+sevent+'_'+slat+'.ps -rotate -90 ../Figures/waccm_vortex_h2o_ES_event_'+sevent+'_'+slat+'.jpg'
endif

endfor	; loop over ES events
end
