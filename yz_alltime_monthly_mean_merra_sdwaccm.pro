;
; plot zonal mean multi-year monthly mean MERRA in T, U, Mark
; /Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_????????.nc3
;
@rd_merra_nc3

loadct,39
mcolor=!p.color
mcolor=byte(!p.color)
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=1000
nydim=700
xorig=[0.1,0.4,0.7]
yorig=[0.3,0.3,0.3]
xlen=0.25
ylen=0.4
cbaryoff=0.1
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif

dirm='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirw='/Volumes/cloud/data/WACCM_data/Datfiles_SD_New/f.e11.FWTREFC1SD.f19.f19.ccmi30.001.cam.h5.'
;
; read monthly means
;
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
for m=0,nmonths-1 do begin
    smon=string(format='(i2.2)',m+1L)
    kfile=0L
    print,smon

    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,font_size=9
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_sdw+merra_Ubar_markbar_'+smon+'.ps'
       !p.charsize=2
       !p.thick=2
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

; read monthly means 
    restore,dirm+mon(m)+'avg.sav'	;,nc,nr,nth,alon,alat,th,p_mean,z_mean,u_mean,mark_mean
    mlon=alon
    mlat=alat
    thm=th
    nthm=n_elements(thm)
    pm_mean=p_mean
    um_mean=u_mean
    zm_mean=z_mean
    markm_mean=mark_mean
    restore,dirw+mon(m)+'avg.sav'       ;,nc,nr,nth,alon,alat,th,p_mean,z_mean,u_mean,mark_mean
;
; plot multi-year monthly mean zonal mean
;
    loadct,39
    erase
    xyouts,.45,.8,strupcase(strmid(mon(m),0,3)),/normal,charsize=3,color=0,charthick=2
    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    plt=mean(um_mean,dim=2)
    mplt=mean(markm_mean,dim=2)
    zplt=mean(zm_mean,dim=2)
    index=where(zplt eq 0.)
;   if index(0) ne -1L then zplt(index)=0./0.
    nlvls=21
    ulevel=-100+10.*findgen(nlvls)
    col1=3+indgen(nlvls)*mcolor/nlvls
    contour,plt,alat,zplt,levels=ulevel,/fill,/cell_fill,c_color=col1,/noeras,title='MERRA',xticks=6,xtitle='Latitude',color=0,yrange=[10,80]
    index=where(ulevel gt 0.)
    contour,plt,alat,zplt,levels=ulevel(index),/follow,/overplot,c_color=0,/noeras
    index=where(ulevel lt 0.)
    contour,plt,alat,zplt,levels=ulevel(index),/follow,/overplot,c_color=mcolor,c_linestyle=5,/noeras
    contour,mplt,alat,zplt,levels=[0.1,0.9],/follow,color=0,thick=10,/noeras,/overplot

mark_mean_interp=0.*plt
u_mean_interp=0.*plt
z_mean_interp=0.*plt
x2d=0.*mplt
y2d=zplt
for j=0,nthm-1L do x2d(*,j)=alat

    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    plt=mean(u_mean,dim=2)
    mplt=mean(mark_mean,dim=2)
    zplt=mean(z_mean,dim=2)
    index=where(zplt eq 0.)
    if index(0) ne -1L then zplt(index)=0./0.
    contour,plt,alat,zplt/1000.,levels=ulevel,/fill,/cell_fill,c_color=col1,/noeras,title='SD-WACCM',xticks=6,xtitle='Latitude',ytitle='Altitude (km)',color=0,yrange=[10,80]
    index=where(ulevel gt 0.)
    contour,plt,alat,zplt/1000.,levels=ulevel(index),/follow,/overplot,c_color=0,/noeras
    index=where(ulevel lt 0.)
    contour,plt,alat,zplt/1000.,levels=ulevel(index),/follow,/overplot,c_color=mcolor,c_linestyle=5,/noeras
    contour,mplt,alat,zplt/1000.,levels=[0.1,0.9],/follow,color=0,thick=10,/noeras,/overplot

    imin=min(ulevel)
    imax=max(ulevel)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,min(xorig),xorig(1)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle='Zonal Mean Wind Speed (m/s)'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor
;
; interpolate WACCM to thm grid
;
for j=0L,nr-1L do begin
    mark_mean_interp(j,*)=interpol(reform(mplt(j,*)),th,thm)
    u_mean_interp(j,*)=interpol(reform(plt(j,*)),th,thm)
    z_mean_interp(j,*)=interpol(reform(zplt(j,*)),th,thm)
endfor
index=where(mark_mean_interp lt 0.)
if index(0) ne -1L then mark_mean_interp(index)=0.
index=where(markm_mean lt 0.)
if index(0) ne -1L then markm_mean(index)=0.

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)

    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
level=-30+6.*findgen(11)
nlvls=n_elements(level)
    plt=u_mean_interp-mean(um_mean,dim=2)		; SDW - MERRA
    mplt=(mark_mean_interp-mean(markm_mean,dim=2))	;/mark_mean_interp
    zplt=z_mean_interp
    index=where(zplt eq 0.)
    if index(0) ne -1L then zplt(index)=0./0.
    contour,plt,alat,zplt/1000.,levels=level,/fill,/cell_fill,c_color=col2,/noeras,title='SD-WACCM minus MERRA',xticks=6,xtitle='Latitude',color=0,yrange=[10,80]
    index=where(level gt 0.)
    contour,plt,alat,zplt/1000.,levels=level(index),/follow,/overplot,c_color=0,/noeras
    index=where(level lt 0.)
    contour,plt,alat,zplt/1000.,levels=level(index),/follow,/overplot,c_color=mcolor,c_linestyle=5,/noeras
loadct,39
index=where(mplt gt .05)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=4,color=mcolor*.9
index=where(mplt lt -.05)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=4,color=mcolor*.3

;   contour,mplt,alat,zplt/1000.,levels=[0.01,0.1,1,5,10],/follow,c_color=[150,180,210,230,250],thick=5,/noeras,/overplot
;   contour,mplt,alat,zplt/1000.,levels=[-10,-5,-1,0.1,0.01],/follow,c_color=[0,30,60,80,100,120],thick=5,/noeras,/overplot,c_linestyle=5

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)

    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xorig(2)+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle='dU (m/s)'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col2(j)
        x2=x2+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim yz_sdw+merra_Ubar_markbar_'+smon+'.ps -rotate -90 yz_sdw+merra_Ubar_markbar_'+smon+'.png'
    endif
endfor  ; loop over months
end
