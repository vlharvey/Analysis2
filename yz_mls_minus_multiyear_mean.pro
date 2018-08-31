;
; plot zonal mean temperature deviation from multi-year mean
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
xorig=[0.175]
yorig=[0.25]
xlen=0.7
ylen=0.65
cbaryoff=0.15
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
;start_year=[2007,2008,2009,2010,2011,2012,2013,2014,2015]
;start_date=[-27, -21, -24, -24, -26, -27, -34, -28,-42]
;end_date=[66, 65, 61, 61, 64, 61, 64, 80]
;nyear=n_elements(start_year)

restore,'MLS_YZ_T_anomaly_2004-2015.sav
;
; NH summer date ranges 3/1 to 9/10
; SH summer date ranges 9/20 to 2/28
;
lstmn=12         ; NH
lstdy=15
lstyr=2012
ledmn=2         ; NH
leddy=15
ledyr=2013
lstday=0
ledday=0
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
;
sdate0=strcompress(lstyr,/r)+strcompress(string(format='(I2.2)',lstmn))+strcompress(string(format='(I2.2)',lstdy))
sdate1=strcompress(ledyr,/r)+strcompress(string(format='(I2.2)',ledmn))+strcompress(string(format='(I2.2)',leddy))
;
; extract date range at ilat
;
dateindex=where(SDATE_ALL ge sdate0 and SDATE_ALL le sdate1,nday)
plotarray=reform(TANOM(dateindex,*,*))
;
; loop over days
;
for iday=0L,nday-1L do begin

    sdate=SDATE_ALL(dateindex(iday))
    plotarray2=reform(plotarray(iday,*,*))
    xarray=latgrid
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps'
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
level=-10+2.*findgen(11)
level=[-32,-16,-8,-4,-2,0,2,4,8,16,32]	;-10+2.*findgen(11)
nlvls=n_elements(col2)
index=where(plotarray2 eq 0.0)
if index(0) ne -1L then plotarray2(index)=0./0.
contour,plotarray2,xarray,altitude,/noeras,xrange=[-90,90],yrange=[10.,90.],ytitle='Altitude (km)',charsize=1.75,color=0,xtitle='Latitude',charthick=2,title=sdate,$
        levels=level,c_color=col2,/cell_fill,xticks=6
index=where(level gt 0.)
contour,plotarray2,xarray,altitude,/noeras,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,plotarray2,xarray,altitude,/noeras,levels=level(index),color=mcolor,/follow,/overplot
;contour,plotarray2,xarray,altitude,/noeras,levels=[0.],color=0,thick=3,/follow,/overplot
;
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.25,charthick=2,xtitle='MLS Temperature Anomaly (K)',xtickname=strcompress(long(level),/r),xticks=nlvls-1
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(jj)
x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps -rotate -90 yz_mls_temp_minus_multiyear_mean_'+sdate+'.png'
   spawn,'rm -f yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps'
endif

endfor	; loop over days

end
