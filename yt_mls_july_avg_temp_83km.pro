;
; impose smoothing and adjust yrange as desired
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
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.10
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
rlats=[50,55,60,65,70,75,80]
nlat=n_elements(rlats)
yt_temp=fltarr(nyear,nlat)
;
; loop over latitude and retain July average temperature each year
for ii=0,nlat-1L do begin
    ilat=rlats(ii)
    rlat=float(ilat)
    slat=strcompress(long(rlat),/remove_all)
    yearlab='2007-2018'
    salt='83'
    restore,'timeseries_mls_temp_'+yearlab+'_'+slat+'_pmc_midseason.sav'
    bad=where(MLSPOLARTP_XT eq 0.)
    if bad(0) ne -1L then MLSPOLARTP_XT(bad)=0./0.
    for iyear=0L,nyear-1L do begin
        doy_thisyear=reform(doy_all(*,iyear))
        temp_thisyear=reform(MLSPOLARTP_XT(*,iyear))
        july=where(doy_thisyear ge 212.-5. and doy_thisyear le 212.+5.)
;       aug=where(doy_thisyear gt 212.)
        yt_temp(iyear,ii)=mean(temp_thisyear(july),/Nan)
;       yt_temp(iyear,ii)=mean(temp_thisyear(aug),/Nan)
    endfor
endfor	; loop over latitudes

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yt_mls_july_temp_'+yearlab+'.ps'
   !p.charsize=2
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
imin=140.	;min(yt_temp,/Nan)-2.
imax=175.	;max(yt_temp,/Nan)+2.
nlvls=15
col1=(indgen(nlvls)/float(nlvls))*mcolor
level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
contour,yt_temp,start_year,rlats,/noeras,xrange=[min(start_year),max(start_year)],yrange=[min(rlats),max(rlats)],ytitle='Latitude',$
      charsize=1.5,color=0,charthick=2,xticklen=-0.075,levels=level,/cell_fill
;contour,yt_temp,start_year,rlats,/noeras,levels=level,/follow,color=0,c_labels=1+0*level,/overplot
contour,yt_temp,start_year,rlats,/noeras,/overplot,levels=[145,150,155,160],color=mcolor,c_labels=[1,1,1,1],thick=3
contour,yt_temp,start_year,rlats,/noeras,/overplot,levels=[165,170,175],color=0,c_labels=[1,1,1],thick=3
loadct,39

ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1,charthick=2,xticks=nlvls/2,xtitle='MLS 26 July-5 Aug Average Temperature (K) at '+salt+' km'
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
   spawn,'convert yt_mls_july_temp_'+yearlab+'.ps -rotate -90 yt_mls_july_temp_'+yearlab+'.png'
endif

end
