;
; save the product of the average wind speed at the vortex edge
; multiplied by the normalised gradient of PV wrt Equivalent latitude
; as a function of altitude and time
;
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=11L & lstdy=1L & lstyr=0L 
ledmn=5L & leddy=1L & ledyr=0L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.35]
xlen=0.8
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
dir='/aura3/data/WACCM_data/Datfiles/waccm_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
;
; save file
;
restore,file='waccm_vortex_strength_2003_2004.sav'
ndays=n_elements(sfile)
y1=strmid(sfile(0),0,4)
y2=strmid(sfile(ndays-1),0,4)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='waccm_pvg_'+y1+'-'+y2+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot altitude-time section of vortex strength indicator, look at prod=dpv*speed
;
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,ndays,ndays,1,1],[400.,400.,2000.,2000.,400.],min_value=0.,$
      xrange=[1,ndays],yrange=[400.,2000.],/nodata,charsize=2,$
      ytitle='Theta (K)',title='WACCM '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,6,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),4,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,350.
    plots,kindex(i)+1,400.,/continue,/data
    xyouts,kindex(i)+1,225.,xlab,/data,alignment=0.5,charsize=3
endfor
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=5.*findgen(nlvls)
index=where(prod_zt eq 0.)
if index(0) ne -1 then prod_zt(index)=-999.
contour,prod_zt,1.+findgen(ndays),th,levels=level,/fill,$
        /cell_fill,/overplot,c_color=col1,min_value=-999.
contour,prod_zt,1.+findgen(ndays),th,levels=level,c_color=0,$
        /follow,/overplot,min_value=-999.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
      xtitle='Edge Avg Wind Speed x Normalised PV Gradient',charsize=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim waccm_pvg_'+y1+'-'+y2+'.ps -rotate -90 waccm_pvg_'+y1+'-'+y2+'.jpg'
endif

end
