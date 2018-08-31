;
; average individual monthly averages to get multi-year monthly averages
;
@fillit
@smoothit

setplot='x'
read,'setplot?',setplot
loadct,38
device,decompose=0
mcolor=byte(!p.color)
mcolor=fix(mcolor)
if mcolor eq 0 then mcolor=255
icmm1=mcolor-1
icmm2=mcolor-2
!noeras=1
a=findgen(6)*(2*!pi/6.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.7
ylen=0.3
cbaryoff=0.02
cbarydel=0.02
if setplot ne 'ps' then begin
   lc=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura3/data/ECMWF_data/Datfiles/ecmwf_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nlat=18
ybin=-85.+10.*findgen(nlat)
nmonth=12L
for imonth=0L,nmonth-1L do begin
    spawn,'ls -1 '+dir+'*airmass*'+mon(imonth)+'*.sav',ifiles
;   for i=0L,0L do begin	;n_elements(ifiles)-1L do begin
    for i=0L,n_elements(ifiles)-1L do begin
        restore,ifiles(i)
print,'restored '+ifiles(i)
        if i eq 0L then begin
           HIGH_O3_AVG_ALL=0.*high_o3_avg
           HIGH_O3_NUM_ALL=0L*high_o3_num
           HIGH_O3_SIG_ALL=0.*high_o3_avg
           OUT_O3_AVG_ALL=0.*out_o3_avg
           OUT_O3_NUM_ALL=0L*out_o3_num
           OUT_O3_SIG_ALL=0.*out_o3_avg
           VORTEX_O3_AVG_ALL=0.*vortex_o3_avg
           VORTEX_O3_NUM_ALL=0L*vortex_o3_num
           VORTEX_O3_SIG_ALL=0.*vortex_o3_avg
x2d=0.*vortex_o3_avg
y2d=0.*vortex_o3_avg
for ii=0,nlat-1 do y2d(ii,*)=th
for j=0,n_elements(th)-1 do x2d(*,j)=ybin
        endif
        index=where(high_o3_avg ne -9999.)
        if index(0) ne -1L then begin
        HIGH_O3_AVG_ALL(index)=HIGH_O3_AVG_ALL(index)+high_o3_avg(index)*float(high_o3_num(index))
        HIGH_O3_NUM_ALL(index)=HIGH_O3_NUM_ALL(index)+high_o3_num(index)
        HIGH_O3_SIG_ALL(index)=HIGH_O3_SIG_ALL(index)+high_o3_sig(index)*float(high_o3_num(index))
        endif
        index=where(out_O3_avg ne 0L)
        if index(0) ne -1L then begin
        OUT_O3_AVG_ALL(index)=OUT_O3_AVG_ALL(index)+out_o3_avg(index)*float(out_o3_num(index))
        OUT_O3_NUM_ALL(index)=OUT_O3_NUM_ALL(index)+out_o3_num(index)
        OUT_O3_SIG_ALL(index)=OUT_O3_SIG_ALL(index)+out_o3_sig*float(out_o3_num(index))
        endif
        index=where(vortex_O3_avg ne 0L)
        if index(0) ne -1L then begin
        VORTEX_O3_AVG_ALL(index)=VORTEX_O3_AVG_ALL(index)+vortex_o3_avg(index)*float(vortex_o3_num(index))
        VORTEX_O3_NUM_ALL(index)=VORTEX_O3_NUM_ALL(index)+vortex_o3_num(index)
        VORTEX_O3_SIG_ALL(index)=VORTEX_O3_SIG_ALL(index)+vortex_o3_sig(index)*float(vortex_o3_num(index))
        endif
    endfor
    index=where(HIGH_O3_AVG_ALL gt 0.)
    HIGH_O3_AVG_ALL(index)=HIGH_O3_AVG_ALL(index)/HIGH_O3_NUM_ALL(index)
    HIGH_O3_SIG_ALL(index)=HIGH_O3_SIG_ALL(index)/HIGH_O3_NUM_ALL(index)
    index=where(OUT_O3_AVG_ALL gt 0.)
    OUT_O3_AVG_ALL(index)=OUT_O3_AVG_ALL(index)/OUT_O3_NUM_ALL(index)
    OUT_O3_SIG_ALL(index)=OUT_O3_SIG_ALL(index)/OUT_O3_NUM_ALL(index)
    index=where(VORTEX_O3_AVG_ALL gt 0.)
    VORTEX_O3_AVG_ALL(index)=VORTEX_O3_AVG_ALL(index)/VORTEX_O3_NUM_ALL(index)
    VORTEX_O3_SIG_ALL(index)=VORTEX_O3_SIG_ALL(index)/VORTEX_O3_NUM_ALL(index)

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/portrait,bits=8,filename='ecmwf_3d_o3_airmass_'+mon(imonth)+'.ps'
   device,/color
   device,/inch,xoff=0.05,yoff=.1,xsize=xsize,ysize=ysize
endif
erase
!type=2^2+2^3
set_viewport,.3,.8,.75,.95
nlvls=31
col1=1+mcolor*findgen(nlvls)/nlvls
level=0.5*findgen(nlvls)
index=where(out_o3_avg_all eq 0.)
if index(0) ne -1L then out_o3_avg_all(index)=-9999.
contour,out_o3_avg_all,ybin,th,xrange=[-90.,90.],yrange=[240.,2000.],/noeras,/cell_fill,$
        c_color=col1,levels=level,min_value=-9999.,title='Ambient '+mon(imonth)
contour,out_o3_avg_all,ybin,th,/overplot,levels=level,/follow,color=0,min_value=-9999.

set_viewport,.2,.5,.45,.65
index=where(high_o3_avg_all eq 0.)
if index(0) ne -1L then high_o3_avg_all(index)=-9999.
;
; fill
;
;fillit,out_o3_avg_all,out_o3_avg_allfill
;out_o3_avg_all=out_o3_avg_allfill
;
; smooth
;
;smoothit,out_o3_avg_all,out_o3_avg_allsmooth
;out_o3_avg_all=out_o3_avg_allsmooth

contour,high_o3_avg_all,ybin,th,xrange=[-90.,90.],yrange=[240.,2000.],/noeras,/cell_fill,$
        c_color=col1,levels=level,min_value=-9999.,title='Anticyclone'
contour,high_o3_avg_all,ybin,th,/overplot,levels=level,/follow,color=0,min_value=-9999.
index=where(high_o3_num_all eq 0L)
high_o3_num_all(index)=-9999L
contour,high_o3_num_all,ybin,th,/overplot,levels=[1000.,5000.,10000.,20000.,30000.,40000.,50000.,60000.],$
        /follow,color=mcolor,min_value=-9999.
index=where(high_o3_num_all le 0.)
oplot,x2d(index),y2d(index),psym=4

set_viewport,.2,.5,.15,.35
index=where(vortex_o3_avg_all eq 0.)
if index(0) ne -1L then vortex_o3_avg_all(index)=-9999.
contour,vortex_o3_avg_all,ybin,th,xrange=[-90.,90.],yrange=[240.,2000.],/noeras,/cell_fill,$
        c_color=col1,levels=level,min_value=-9999.,title='Vortex'
contour,vortex_o3_avg_all,ybin,th,/overplot,levels=level,/follow,color=0,min_value=-9999.
index=where(vortex_o3_num_all eq 0L)
vortex_o3_num_all(index)=-9999L
contour,vortex_o3_num_all,ybin,th,/overplot,levels=[1000.,50000,100000.,150000.,200000.,300000.],$
        /follow,color=mcolor,min_value=-9999.
index=where(vortex_o3_num_all le 0.)
oplot,x2d(index),y2d(index),psym=4

set_viewport,.6,.9,.45,.65
nlvls=21
col1=1+mcolor*findgen(nlvls)/nlvls
level=-2.+0.2*findgen(nlvls)
index=where(high_o3_avg_all gt 0. and out_o3_avg_all gt 0.)
dum=0.*out_o3_avg_all
dum(index)=high_o3_avg_all(index)-out_o3_avg_all(index)
index=where(dum eq 0.)
dum(index)=-9999.
contour,dum,ybin,th,xrange=[-90.,90.],yrange=[240.,2000.],/noeras,/cell_fill,$
        c_color=col1,levels=level,min_value=-9999.,title='Anticyclone-Ambient'
index=where(level lt 0.)
contour,dum,ybin,th,/overplot,/noeras,color=0,levels=level(index),c_linestyle=1,min_value=-9999.
index=where(level gt 0.)
contour,dum,ybin,th,/overplot,/noeras,color=mcolor,levels=level(index),c_linestyle=0,min_value=-9999.
contour,dum,ybin,th,/overplot,/noeras,color=0,levels=0,thick=2,min_value=-9999.
index=where(dum ne -9999.)
print,'high ',min(dum(index)),max(dum)

set_viewport,.6,.9,.15,.35
index=where(vortex_o3_avg_all gt 0. and out_o3_avg_all gt 0.)
dum=0.*out_o3_avg_all
dum(index)=vortex_o3_avg_all(index)-out_o3_avg_all(index)
index=where(dum eq 0.)
dum(index)=-9999.
contour,dum,ybin,th,xrange=[-90.,90.],yrange=[240.,2000.],/noeras,/cell_fill,$
        c_color=col1,levels=level,min_value=-9999.,title='Vortex-Ambient'
index=where(level lt 0.)
contour,dum,ybin,th,/overplot,/noeras,color=0,levels=level(index),c_linestyle=1,min_value=-9999.
index=where(level gt 0.)
contour,dum,ybin,th,/overplot,/noeras,color=mcolor,levels=level(index),c_linestyle=0,min_value=-9999.
contour,dum,ybin,th,/overplot,/noeras,color=0,levels=0,thick=2,min_value=-9999.
index=where(dum ne -9999.)
print,'vortex ',min(dum(index)),max(dum)

if setplot eq 'ps' then device,/close
if setplot ne 'ps' then stop
endfor		; loop over months
end
