;
;
; plot zt vortex area and MLS stratopause altitude
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

sver='v3.3'
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
xlen=0.75
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
dir='/Users/harvey/Harvey_etal_2014/Post_process/vortex_area_ES_'
dirm='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
;restore, '/Users/harvey/Desktop/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
;result=size(MAXHEIGHTTHETA)
;nevents=result(1)

mlsesdates=['20060130','20090205','20120130','20130123']
nevents=n_elements(mlsesdates)
for iES = 0L, nevents - 1L do begin
;
; AREA_ZT_NC4     FLOAT     = Array[121, 30]
; KOFF            INT       =       60
; SFILE           STRING    = Array[121]
; TH              FLOAT     = Array[30]
; Y1              STRING    = '2001'
; Y2              STRING    = '2031'
; ZTAVG           FLOAT     = Array[121]
; ZTMAX           FLOAT     = Array[121]
;
    esdate0=mlsesdates(ies)
    restore,dir+esdate0+'_merra.sav'

       kday=n_elements(sfile)
       nth=n_elements(th)
       area_composite=fltarr(kday,nth)
       ztavg_composite=fltarr(kday)	; height of polar cap average stratopause
if ies eq 0 then begin
       ztavg_all=fltarr(kday,nevents)	; height of polar cap average stratopause
       ztmax_all=fltarr(kday,nevents)	; height of polar cap average stratopause
endif
       ztmax_composite=fltarr(kday)	; height of maximum temperature anywhere in the polar cap
    area_composite=AREA_ZT_NC4
    ztavg_composite=ztavg
    ztmax_composite=ztmax
    ztavg_all(*,ies)=ztavg
    ztmax_all(*,ies)=ztmax
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_vortex_area_ES_event_'+esdate0+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot zt vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=18
col1=1+indgen(nlvls)*icolmax/nlvls
level=2.*findgen(nlvls)
index=where(area_composite eq 0.)
if index(0) ne -1L then area_composite(index)=0./0.
area_composite=smooth(area_composite,3,/nan)
xmax=30.
contour,area_composite,-1.*koff+findgen(kday),th,color=0,xtitle='Days From ES Onset ('+esdate0+')',thick=6,yrange=[500.,max(th)],/noeras,ytitle='Theta (K)',/cell_fill,$
     levels=level,c_color=col1,charthick=2,charsize=1.5,xrange=[-30.,xmax]
contour,area_composite,-1.*koff+findgen(kday),th,levels=level,c_color=0,/follow,/overplot,c_labels=0*level,min_value=0
contour,area_composite,-1.*koff+findgen(kday),th,levels=[40,50,60],c_color=mcolor,/follow,/overplot,c_labels=0*level,min_value=0
plots,0,500
plots,0,max(th),/continue,color=mcolor,thick=5
;oplot,-1.*koff+findgen(kday),ztavg_composite,psym=8,color=0
;for i=0L,nevents-1L do oplot,-1.*koff+findgen(kday),ztavg_all(*,i),psym=8,color=mcolor*(float(i)/5.),symsize=.5
;for i=0L,kday-1L do oplot,[-1.*koff+i,-1.*koff+i],[max(ztavg_all(i,*)),max(ztavg_all(i,*))],psym=8,color=0
;for i=0L,kday-1L do oplot,[-1.*koff+i,-1.*koff+i],[max(ztmax_all(i,*)),max(ztmax_all(i,*))],psym=8,color=0
;oplot,-1.*koff+findgen(kday),ztmax_composite,psym=8,color=mcolor*.9
xyouts,xmax+10.,1500.,'Approximate Altitude (km)',color=0,/data,orientation=90.,charsize=1.5,charthick=2
xyouts,xmax+1.,500.,'20',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,800.,'30',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,1400.,'40',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,2000.,'50',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,3000.,'60',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,4000.,'70',color=0,/data,charsize=1.5,charthick=2
xyouts,xmax+1.,5000.,'80',color=0,/data,charsize=1.5,charthick=2
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='MERRA ES Vortex Area (% NH)',charthick=2,charsize=1.5
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
   spawn,'convert -trim ../Figures/merra_vortex_area_ES_event_'+esdate0+'.ps -rotate -90 ../Figures/merra_vortex_area_ES_event_'+esdate0+'.jpg'
endif

endfor

end
