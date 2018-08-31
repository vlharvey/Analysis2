; 
; WACCM-X January
; quick plot from saved file
;
loadct,39
restore,'zt_height_waccmx_Traject_Season_70-90N.sav
device,decompose=0
result=size(zt_co)
ksteps=result(1)

set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='zt_height_waccmx_Traject_Season_70-90N.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif

xorig=[0.075,0.075,0.075,0.55,0.55,0.55]
yorig=[0.75,0.45,0.15,0.75,0.45,0.15]-0.05
xlen=0.3
ylen=0.2
;
; zt array
;
;nhindex=where(alat le -70)
;for k=0L,nl-1L do begin
;    zt_height(icount,k)=mean(aoa2grd(*,nhindex,k))
;    zt_no(icount,k)=mean(no_conc(*,nhindex,k))
;    zt_co(icount,k)=mean(cogrd(*,nhindex,k))
;    zt_t(icount,k)=mean(tgrd(*,nhindex,k))
;    zt_u(icount,k)=mean(ugrd(*,nhindex,k))
;endfor
;
;for k=0L,nl-1L do zlev(k)=mean(zgrd(*,*,k))

zlev=lev
    erase
    xyouts,.25,.95,'WACCM-X >70N Average',/normal,charsize=2,color=0
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=[0.,100.,1000.,1.e4,1.e5,5.e5,1.e6,5.e6,1.e7,5.e7,1.e8,5.e8,1.e9]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,zt_no,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Days',title='NO',/ylog
    index=where(level gt 0)
    contour,zt_no,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(molecules/cm!u3!n)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=[0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.,1.5,2.,5.,10.,15.,20.,25.,30.,40.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,zt_co,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.001],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Days',title='CO',/ylog
    index=where(level gt 0)
    contour,zt_co,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    contour,zt_co,findgen(ksteps),zlev,/overplot,levels=5,/follow,/noeras,color=0,thick=3,c_labels=[0]

    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(ppmv)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    imin=0.
    imax=20
    imax=80
    nlvls=21
    level=imin+((imax-imin)/float(nlvls))*findgen(nlvls+1)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
;   contour,zt_height,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Days',title='AOA=Altitude',/ylog
    contour,zt_gph,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Days',title='GPH',/ylog
    index=where(level gt 0)
    contour,zt_gph,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='0=bottom, 1=top'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=120.+10.*findgen(21)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,zt_t,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Days',title='Temperature',/ylog
    index=where(level gt 0)
    contour,zt_t,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(K)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(4)
    xmx=xorig(4)+xlen
    ymn=yorig(4)
    ymx=yorig(4)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=-50.+10.*findgen(14)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,zt_u,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Days',title='Zonal Wind',/ylog
    index=where(level gt 0)
    contour,zt_u,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    index=where(level lt 0)
    contour,zt_u,findgen(ksteps),zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=mcolor,c_linestyle=5
;   contour,zt_u,findgen(ksteps),zlev,/overplot,levels=[0],/follow,c_labels=[0],/noeras,color=0,thick=5
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(m/s)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(5)
    xmx=xorig(5)+xlen
    ymn=yorig(5)
    ymx=yorig(5)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=-10.+findgen(21)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,zt_v,findgen(ksteps),zlev,xrange=[0.,ksteps-1],yrange=[100.,.01],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Days',title='Meridional Wind',/ylog
    index=where(level gt 0)
    contour,zt_v,findgen(ksteps),zlev,/overplot,levels=[5.,10.,20.,30.],/follow,c_labels=[0,0,0,0],/noeras,color=0
    contour,zt_v,findgen(ksteps),zlev,/overplot,levels=[-40.,-30.,-20.,-10.,-5.],/follow,c_labels=[0,0,0,0,0],/noeras,color=mcolor,c_linestyle=5
;   contour,zt_v,findgen(ksteps),zlev,/overplot,levels=[0],/follow,c_labels=[0],/noeras,color=0,thick=5
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(m/s)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim zt_height_waccmx_Traject_Season_70-90N.ps -rotate -90 zt_height_waccmx_Traject_Season_70-90N.jpg'
       spawn,'rm -f zt_height_waccmx_Traject_Season_70-90N.ps'
    endif
end
