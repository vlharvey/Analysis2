;
; plot altitude time section of zonal mean temperature deviation from multi-year mean
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
xorig=[0.15,0.15]
yorig=[0.65,0.15]
xlen=0.7
ylen=0.3
cbaryoff=0.1
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
print,latgrid
rlat=-70.
;read,'Enter NH latitude ',rlat
index=where(abs(rlat-latgrid) eq min(abs(rlat-latgrid)))
ilat=index(0)
slat=strcompress(long(latgrid(ilat)),/r)
rlat2=70.
;read,'Enter SH latitude ',rlat2
index=where(abs(rlat2-latgrid) eq min(abs(rlat2-latgrid)))
ilat2=index(0)
slat2=strcompress(long(latgrid(ilat2)),/r)
;
; extract date range at ilat
;
dateindex=where(SDATE_ALL ge sdate0 and SDATE_ALL le sdate1)
plotarray=reform(TANOM(dateindex,ilat,*))
plotarray2=reform(TANOM(dateindex,ilat2,*))
xarray=DFS_ALL(dateindex)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_mls_temp_minus_multiyear_mean_'+sdate0+'-'+sdate1+'_'+slat+'.ps'
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
level=[-12,-8,-4,-2,-1,0,1,2,4,8,12]	;-10+2.*findgen(11)
nlvls=n_elements(col2)
index=where(plotarray eq 0.0)
if index(0) ne -1L then plotarray(index)=0./0.
contour,plotarray,xarray,altitude,/noeras,xrange=[min(xarray),max(xarray)],yrange=[10.,90.],ytitle='Altitude (km)',charsize=1.5,color=0,xtitle='DFS',charthick=2,title=sdate0+'-'+sdate1+'         '+slat,$
        levels=level,c_color=col2,/cell_fill
index=where(level gt 0.)
contour,plotarray,xarray,altitude,/noeras,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,plotarray,xarray,altitude,/noeras,levels=level(index),color=mcolor,/follow,/overplot
;contour,plotarray,xarray,altitude,/noeras,levels=[0.],color=0,thick=3,/follow,/overplot
;
; CIPS frequencies
;
syr=strmid(strcompress(lstyr,/r),2,2)
if rlat gt 0. then pre='NH'+syr
if rlat lt 0. then pre='SH'+syr
restore,'/Volumes/Data/CIPS_data/Pre_process/Line_plots/F_V4.20_80Lat_2G_'+pre+'.sav
print,'F_V4.20_80Lat_2G_'+pre+'.sav
axis,/yax,yrange=[0,100],/save,ytitle='CIPS Frequency (%)',color=0,charsize=1.5,charthick=2
oplot,ddd,freq,color=0,thick=6

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.2,charthick=2,xtitle='MLS Temperature Anomaly (K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(jj)
x1=x1+dx
endfor


xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7       ; ticks outward
level=-35+7.*findgen(11)
nlvls=n_elements(col2)
index=where(plotarray2 eq 0.0)
if index(0) ne -1L then plotarray2(index)=0./0.
contour,plotarray2,xarray,altitude,/noeras,xrange=[min(xarray),max(xarray)],yrange=[10.,90.],ytitle='Altitude (km)',charsize=1.5,color=0,xtitle='DFS',charthick=2,title='                                   '+slat2,$
        levels=level,c_color=col2,/cell_fill
index=where(level gt 0.)
contour,plotarray2,xarray,altitude,/noeras,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,plotarray2,xarray,altitude,/noeras,levels=level(index),color=mcolor,/follow,/overplot
;contour,plotarray2,xarray,altitude,/noeras,levels=[0.],color=0,thick=3,/follow,/overplot,charthick=2
axis,/yax,yrange=[0,100],/save,ytitle='CIPS Frequency (%)',color=0,charsize=1.5,charthick=2
oplot,ddd,freq,color=0,thick=6

imin=min(level)
imax=max(level)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.2,charthick=2,xtitle='MLS Temperature Anomaly (K)'
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
   spawn,'convert zt_mls_temp_minus_multiyear_mean_'+sdate0+'-'+sdate1+'_'+slat+'.ps -rotate -90 zt_mls_temp_minus_multiyear_mean_'+sdate0+'-'+sdate1+'_'+slat+'.png'
   spawn,'rm -f zt_mls_temp_minus_multiyear_mean_'+sdate0+'-'+sdate1+'_'+slat+'.ps'
endif

end
