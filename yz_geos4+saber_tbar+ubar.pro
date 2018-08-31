;
; daily zonal mean SABER vs GEOS-4 temperature and zonal wind
; interpolate GEOS-4 to SABER pressure levels
;
@kdate
@rd_ukmo_nc3
@rd_geos5_dat

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
xorig=[0.11,0.4,0.69,0.11,0.4,0.69]
yorig=[0.62,0.62,0.62,0.18,0.18,0.18]
xlen=0.25
ylen=0.25
cbaryoff=0.08
cbarydel=0.015
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
gdir='/aura7/harvey/GEOS4_data/Datfiles/'
spawn,'ls '+gdir+'DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.origp.*_1200.V01.sav',ncfiles
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
;
; restore GEOS data
;
    restore,file=ncfiles(ifile)        ; alonold,alatold,press,znew,unew,vnew,tnew,pvnew,qnew
    glon=alonold
    glat=alatold
    p=press
    nlg=n_elements(glon)
    nlat=n_elements(glat)
    nlv=n_elements(p)
    ubar=fltarr(nlat,nlv)
    tbar=fltarr(nlat,nlv)
    for k=0L,nlv-1L do begin
    for j=0L,nlat-1L do begin
        ubar(j,k)=total(unew(*,j,k))/float(nlg)
        tbar(j,k)=total(tnew(*,j,k))/float(nlg)
    endfor
    endfor
;
; read SABER winds
;
    result=strsplit(ncfiles(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(7),'_',/extract)
    sdate=result3(0)
    sfile='GRID_PHI_WINDS.'+sdate+'.sav'
    sdum=findfile('/aura6/data/SABER_data/Datfiles_winds/'+sfile)
    if sdum(0) eq '' then goto,jump
    restore,sdum(0)
    tdata=t3d	; avoid T3D intrinsic function
    if max(tdata) eq 0. then goto,jump
    udata=u3d
;
; compute Ubar
;
    snlat=n_elements(alat)
    snlg=n_elements(alon)
    snlv=n_elements(press)
    ubar_saber=fltarr(snlat,snlv)
    tbar_saber=fltarr(snlat,snlv)
    for k=0L,snlv-1L do begin
    for j=0L,snlat-1L do begin
        index=where(udata(*,j,k) ne 0.,ngood)
        if index(0) ne -1L then begin
           ubar_saber(j,k)=total(udata(index,j,k))/float(ngood)
           tbar_saber(j,k)=total(tdata(index,j,k))/float(ngood)
        endif
    endfor
    endfor
;
; intepolate Ubar,Tbar to SABER latitude (glat->alat) and pressure (p->press) grid
;
    ubar_saberp=fltarr(snlat,snlv)
    tbar_saberp=fltarr(snlat,snlv)
    for jj=0L,snlat-1L do begin
        yp=alat(jj)
        for kk=0L,snlv-1L do begin
            pp=alog(press(kk))
            for j=0L,nlat-2L do begin
                jp1=j+1
                xlat=glat(j)
                xlatp1=glat(jp1)
                if yp le xlat and yp ge xlatp1 then begin
                   yscale=(yp-xlat)/(xlatp1-xlat)
                   for k=0L,nlv-2L do begin
                       kp1=k+1
                       p1=alog(p(k))
                       p2=alog(p(kp1))
                       if pp ge p1 and pp le p2 then begin
                          pscale=(pp-p1)/(p2-p1)
                          q1=tbar(j,k)+yscale*(tbar(jp1,k)-tbar(j,k))
                          q2=tbar(j,kp1)+yscale*(tbar(jp1,kp1)-tbar(j,kp1))
                          tbar_saberp(jj,kk)=q1+pscale*(q2-q1)
                          q1=ubar(j,k)+yscale*(ubar(jp1,k)-ubar(j,k))
                          q2=ubar(j,kp1)+yscale*(ubar(jp1,kp1)-ubar(j,kp1))
                          ubar_saberp(jj,kk)=q1+pscale*(q2-q1)
goto,jumplev
                       endif
                   endfor
                endif
           endfor
jumplev:
        endfor
    endfor

    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_geos4+saber_'+sdate+'_tbar+ubar.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
    erase
    xyouts,.45,.95,sdate,/normal,charsize=2,color=0,charthick=3
    level=170.+5.*findgen(25)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*mcolor/nlvls
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,tbar_saberp,alat,press,/ylog,charsize=1.5,/noeras,xtitle='Latitude',yrange=[100.,0.0001],$
            ytitle='Pressure (hPa)',levels=level,c_color=col1,/cell_fill,color=0,min_value=0.,title='GEOS-4 Tbar',$
            xrange=[-90.,90.],xticks=6,charthick=1.5
    index=where(level gt 0.)
    contour,tbar_saberp,alat,press,levels=level(index),color=0,/follow,/noeras,/overplot,min_value=0.
    imin=min(level)
    imax=max(level)
    ymnb=yorig(0)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(K)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor
    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,tbar_saber,alat,press,/ylog,charsize=1.5,/noeras,xtitle='Latitude',yrange=[100.,0.0001],$
            levels=level,c_color=col1,/cell_fill,color=0,min_value=0.,title='SABER Tbar',xrange=[-90.,90.],xticks=6,$
            ytickname=[' ',' ',' ',' ',' ',' ',' '],charthick=1.5
    index=where(level gt 0.)
    contour,tbar_saber,alat,press,levels=level(index),color=0,/follow,/noeras,/overplot,min_value=0.
    imin=min(level)
    imax=max(level)
    ymnb=yorig(1)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(K)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor
;
; difference
;
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    col2=1+indgen(11)
    index=where(tbar_saberp ne 0. and tbar_saber ne 0.)
    pdiff=0.*tbar_saberp
;   pdiff(index)=100.*(tbar_saberp(index)-tbar_saber(index))/tbar_saber(index)
    pdiff(index)=(tbar_saberp(index)-tbar_saber(index))
    index=where(pdiff eq 0.)
    if index(0) ne -1L then pdiff(index)=-99.
    level=-20.+4.*findgen(11)
    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,pdiff,alat,press,/ylog,levels=level,/cell_fill,c_color=col2,color=0,min_value=-99.,$
            yrange=[100.,0.0001],xrange=[-90.,90.],xticks=6,title='G4 minus SABER',charsize=1.5,$
            ytickname=[' ',' ',' ',' ',' ',' ',' '],charthick=1.5,xtitle='Latitude'
    contour,pdiff,alat,press,/overplot,levels=level,color=0,/follow,min_value=-99.,c_labels=0*level
    imin=min(level)
    imax=max(level)
    ymnb=yorig(2)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(K)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(n_elements(level))
    for j=0,n_elements(level)-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col2(j)
        x1=x1+dx
    endfor
    loadct,38

    level=-120.+10.*findgen(25)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*mcolor/nlvls
    !type=2^2+2^3
    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    index=where(ubar_saberp eq 0.)
    if index(0) ne -1L then ubar_saberp(index)=-9999.
    contour,ubar_saberp,alat,press,/ylog,charsize=1.5,/noeras,xtitle='Latitude',yrange=[100.,0.0001],$
            ytitle='Pressure (hPa)',levels=level,c_color=col1,/cell_fill,color=0,min_value=-9999.,title='GEOS-4 Ubar',$
            xrange=[-90.,90.],xticks=6,charthick=1.5
    index=where(level gt 0.)
    contour,ubar_saberp,alat,press,levels=level(index),color=0,/follow,/noeras,/overplot,min_value=-9999.
    index=where(level lt 0.)
    contour,ubar_saberp,alat,press,levels=level(index),color=mcolor,/follow,/noeras,/overplot,min_value=-9999.,$
            c_linestyle=2
    imin=min(level)
    imax=max(level)
    ymnb=yorig(3)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(m/s)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor
    !type=2^2+2^3
    xmn=xorig(4)
    xmx=xorig(4)+xlen
    ymn=yorig(4)
    ymx=yorig(4)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    index=where(ubar_saber eq 0.)
    if index(0) ne -1L then ubar_saber(index)=-9999.
    contour,ubar_saber,alat,press,/ylog,charsize=1.5,/noeras,xtitle='Latitude',yrange=[100.,0.0001],$
            levels=level,c_color=col1,/cell_fill,color=0,min_value=-9999.,title='SABER Ubar',xrange=[-90.,90.],xticks=6,$
            ytickname=[' ',' ',' ',' ',' ',' ',' '],charthick=1.5
    index=where(level gt 0.)
    contour,ubar_saber,alat,press,levels=level(index),color=0,/follow,/noeras,/overplot,min_value=-9999.
    index=where(level lt 0.)
    contour,ubar_saber,alat,press,levels=level(index),color=mcolor,/follow,/noeras,/overplot,min_value=-9999.,$
            c_linestyle=2
    imin=min(level)
    imax=max(level)
    ymnb=yorig(4)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(m/s)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor
;
; difference
;
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    col2=1+indgen(11)
    index=where(ubar_saberp ne -9999. and ubar_saber ne -9999.)
    pdiff=0.*ubar_saberp
;   pdiff(index)=100.*(ubar_saberp(index)-ubar_saber(index))/ubar_saber(index)
    pdiff(index)=(ubar_saberp(index)-ubar_saber(index))
    index=where(ubar_saberp eq -9999. or ubar_saber eq -9999.)
    if index(0) ne -1L then pdiff(index)=-9999.
    level=-50.+10.*findgen(11)
    !type=2^2+2^3
    xmn=xorig(5)
    xmx=xorig(5)+xlen
    ymn=yorig(5)
    ymx=yorig(5)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,pdiff,alat,press,/ylog,levels=level,/cell_fill,c_color=col2,color=0,min_value=-9999.,$
            yrange=[100.,0.0001],xrange=[-90.,90.],xticks=6,title='G4 minus SABER',charsize=1.5,$
            ytickname=[' ',' ',' ',' ',' ',' ',' '],charthick=1.5,xtitle='Latitude'
    contour,pdiff,alat,press,/overplot,levels=level,color=0,/follow,min_value=-9999.,c_labels=0*level
    imin=min(level)
    imax=max(level)
    ymnb=yorig(5)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,xtitle='(m/s)',charthick=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(n_elements(level))
    for j=0,n_elements(level)-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col2(j)
        x1=x1+dx
    endfor
    loadct,38

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert yz_geos4+saber_'+sdate+'_tbar+ubar.ps -rotate -90 '+$
             ' yz_geos4+saber_'+sdate+'_tbar+ubar.jpg'
    endif
    jump:
endfor		; loop over time steps
end
