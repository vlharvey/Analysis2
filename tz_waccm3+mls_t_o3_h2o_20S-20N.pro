;
; plot altitude-time section of tropical (20S-20N) T, O3, and water vapor
;
loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
yorig=[0.7,0.7,0.7,0.4,0.4,0.4,0.1,0.1,0.1]
xorig=[0.1,0.4,0.7,0.1,0.4,0.7,0.1,0.4,0.7]
ylen=0.2
xlen=0.25
cbaryoff=0.06
cbarydel=0.01
nxdim=800
nydim=800
cbaryoff=0.04
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
mlab=['J','F','M','A','M','J','J','A','S','O','N','D']
!noeras=1
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_TNV3/wa3_tnv3_'
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
nmonth=n_elements(month)
for imonth=0,nmonth-1 do begin
    spawn,'ls '+dir+mon(imonth)+'_avg.sav',ifile
    restore,ifile
    pv2=pv_mean
    p2=p_mean
    ch4=ch4_mean
    u2=u_mean
    v2=v_mean
    no2=no2_mean
    mark2=mark_mean
    sf2=sf_mean
    h2o=h2o_mean
    o3=o3_mean
    tmp2=0.*p2
    for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
;
; WACCM zonal means
;
    h2oyth_mean_waccm=fltarr(nr,nth)
    o3yth_mean_waccm=fltarr(nr,nth)
    tempyth_mean_waccm=fltarr(nr,nth)
    for k=0L,nth-1 do begin
    for j=0L,nr-1 do begin
        h2oyth_mean_waccm(j,k)=total(h2o(j,*,k))/float(nc)
        o3yth_mean_waccm(j,k)=total(o3(j,*,k))/float(nc)
        tempyth_mean_waccm(j,k)=total(tmp2(j,*,k))/float(nc)
    endfor
    endfor
;
; restore monthly mean MLS data
;
    restore,mdir+mon(imonth)+'_avg_yth.sav'	;nr,nth,alat,th,pyth_mean,h2oyth_mean,o3yth_mean,tempyth_mean

if imonth eq 0L then begin
   h2ozt=fltarr(nmonth,nth)
   o3zt=fltarr(nmonth,nth)
   tzt=fltarr(nmonth,nth)
   h2ozt_waccm=fltarr(nmonth,nth)
   o3zt_waccm=fltarr(nmonth,nth)
   tzt_waccm=fltarr(nmonth,nth)
   yindex=where(alat ge -20. and alat le 20.,nlat)
endif
;
; average tropical profile
;
   for k=0L,nth-1 do begin
       h2ozt(imonth,k)=total(h2oyth_mean(yindex,k))/float(nlat)
       o3zt(imonth,k)=total(o3yth_mean(yindex,k))/float(nlat)
       tzt(imonth,k)=total(tempyth_mean(yindex,k))/float(nlat)
       h2ozt_waccm(imonth,k)=total(h2oyth_mean_waccm(yindex,k))/float(nlat)
       o3zt_waccm(imonth,k)=total(o3yth_mean_waccm(yindex,k))/float(nlat)
       tzt_waccm(imonth,k)=total(tempyth_mean_waccm(yindex,k))/float(nlat)
   endfor

endfor	; loop over months

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='zt_waccm3+mls_20S-20N_t_o3_h2o.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
       !p.thick=2.0                   ;Plotted lines twice as thick
       !p.charsize=2.0
    endif
;
; MLS T
;
    loadct,39
    erase
    xyouts,.45,.95,'20S-20N',/normal,charsize=2,charthick=2,color=0
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    tlevel=130.+5.*findgen(31)
    nlvls=n_elements(tlevel)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,tzt,1.+findgen(nmonth),th,levels=tlevel,c_color=col1,/noeras,/cell_fill,$
            title='MLS Temperature',color=0,min_value=0.
    contour,tzt,1.+findgen(nmonth),th,levels=tlevel,color=0,/noeras,/follow,/overplot,$
            min_value=0.
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(tlevel)
    imax=max(tlevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(K)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; WACCM T
;
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    contour,tzt_waccm,1.+findgen(nmonth),th,levels=tlevel,c_color=col1,/noeras,/cell_fill,$
            title='WACCM Temperature',color=0,min_value=0.
    contour,tzt_waccm,1.+findgen(nmonth),th,levels=tlevel,color=0,/noeras,/follow,/overplot,$
            min_value=0.
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(tlevel)
    imax=max(tlevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(K)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; T difference
;
    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    tdiff=0.*tzt
    index=where(tzt gt 0. and tzt_waccm gt 0.)
    tdiff(index)=(tzt_waccm(index)-tzt(index))
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    nlvls2=11
    dlevel=-50.+10.*findgen(nlvls2)
    col2=1+indgen(nlvls2)
    contour,tdiff,1.+findgen(nmonth),th,levels=dlevel,c_color=col2,/noeras,/cell_fill,$
            title='WACCM-MLS',color=0
    index=where(dlevel gt 0.)
    contour,tdiff,1.+findgen(nmonth),th,levels=dlevel(index),color=0,/noeras,/follow,/overplot
    index=where(dlevel lt 0.)
    contour,tdiff,1.+findgen(nmonth),th,levels=dlevel(index),color=mcolor,/noeras,/follow,/overplot,c_linestyle=1
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(dlevel)
    imax=max(dlevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(K)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls2)
    for j=0,nlvls2-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col2(j)
      x1=x1+dx
    endfor
;
; MLS ozone
;
    loadct,39
    !type=2^2+2^3
    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    o3level=0.5+0.5*findgen(24)
    nlvls=n_elements(o3level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,o3zt*1.e6,1.+findgen(nmonth),th,levels=o3level,c_color=col1,/noeras,/cell_fill,$
            title='MLS O!l3!n',color=0,min_value=0.,yrange=[min(th),2000.]
    contour,o3zt*1.e6,1.+findgen(nmonth),th,levels=o3level,color=0,/noeras,/follow,/overplot,$
            min_value=0.
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(o3level)
    imax=max(o3level)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppmv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; WACCM ozone
;
    !type=2^2+2^3
    xmn=xorig(4)
    xmx=xorig(4)+xlen
    ymn=yorig(4)
    ymx=yorig(4)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,o3zt_waccm,1.+findgen(nmonth),th,levels=o3level,c_color=col1,/noeras,/cell_fill,$
            title='WACCM O!l3!n',color=0,yrange=[min(th),2000.]
    contour,o3zt_waccm,1.+findgen(nmonth),th,levels=o3level,color=0,/noeras,/follow,/overplot
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(o3level)
    imax=max(o3level)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppmv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; ozone difference
;
    !type=2^2+2^3
    xmn=xorig(5)
    xmx=xorig(5)+xlen
    ymn=yorig(5)
    ymx=yorig(5)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    o3diff=0.*o3zt
    index=where(o3zt gt 0. and o3zt_waccm gt 0.)
    o3zt(index)=o3zt(index)*1.e6
    o3diff(index)=100.*(o3zt_waccm(index)-o3zt(index))/o3zt(index)
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    nlvls2=11
    dlevel=-30.+6.*findgen(nlvls2)
    col2=1+indgen(nlvls2)
    contour,o3diff,1.+findgen(nmonth),th,levels=dlevel,c_color=col2,/noeras,/cell_fill,$
            title='(WACCM-MLS)/MLS',color=0,yrange=[min(th),2000.]
    index=where(dlevel gt 0.)
    contour,o3diff,1.+findgen(nmonth),th,levels=dlevel(index),color=0,/noeras,/follow,/overplot
    index=where(dlevel lt 0.)
    contour,o3diff,1.+findgen(nmonth),th,levels=dlevel(index),color=mcolor,/noeras,/follow,/overplot,c_linestyle=1
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(dlevel)
    imax=max(dlevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(%)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls2)
    for j=0,nlvls2-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col2(j)
      x1=x1+dx
    endfor
;
; MLS water
;
    loadct,39
    !type=2^2+2^3
    xmn=xorig(6)
    xmx=xorig(6)+xlen
    ymn=yorig(6)
    ymx=yorig(6)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    h2olevel=0.25+0.25*findgen(20)
    nlvls=n_elements(h2olevel)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,h2ozt*1.e6,1.+findgen(nmonth),th,levels=h2olevel,c_color=col1,/noeras,/cell_fill,$
            title='MLS H!l2!nO',color=0,yrange=[400.,1000.],$
            min_value=0.
    contour,h2ozt*1.e6,1.+findgen(nmonth),th,levels=h2olevel,color=0,/noeras,/follow,/overplot,$
            min_value=0.
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(h2olevel)
    imax=max(h2olevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppmv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; WACCM water
;
    !type=2^2+2^3
    !type=2^2+2^3
    xmn=xorig(7)
    xmx=xorig(7)+xlen
    ymn=yorig(7)
    ymx=yorig(7)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,h2ozt_waccm,1.+findgen(nmonth),th,levels=h2olevel,c_color=col1,/noeras,/cell_fill,$
            title='WACCM H!l2!nO',color=0,yrange=[400.,1000.]
    contour,h2ozt_waccm,1.+findgen(nmonth),th,levels=h2olevel,color=0,/noeras,/follow,/overplot
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(h2olevel)
    imax=max(h2olevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppmv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; water difference
;
    !type=2^2+2^3
    xmn=xorig(8)
    xmx=xorig(8)+xlen
    ymn=yorig(8)
    ymx=yorig(8)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    h2odiff=0.*h2ozt
    index=where(h2ozt gt 0. and h2ozt_waccm gt 0.)
    h2ozt(index)=h2ozt(index)*1.e6
    h2odiff(index)=100.*(h2ozt_waccm(index)-h2ozt(index))/h2ozt(index)
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    nlvls2=11
    dlevel=-30.+6.*findgen(nlvls2)
    col2=1+indgen(nlvls2)
    contour,h2odiff,1.+findgen(nmonth),th,levels=dlevel,c_color=col2,/noeras,/cell_fill,$
            title='(WACCM-MLS)/MLS',color=0,yrange=[min(th),1000.]
    index=where(dlevel gt 0.)
    contour,h2odiff,1.+findgen(nmonth),th,levels=dlevel(index),color=0,/noeras,/follow,/overplot
    index=where(dlevel lt 0.)
    contour,h2odiff,1.+findgen(nmonth),th,levels=dlevel(index),color=mcolor,/noeras,/follow,/overplot,c_linestyle=1
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(dlevel)
    imax=max(dlevel)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(%)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls2)
    for j=0,nlvls2-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col2(j)
      x1=x1+dx
    endfor

    !p.charthick=1.
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim zt_waccm3+mls_20S-20N_t_o3_h2o.ps -rotate -90 zt_waccm3+mls_20S-20N_t_o3_h2o.jpg'
;      spawn,'/usr/bin/rm zt_waccm3+mls_20S-20N_t_o3_h2o.ps'
    endif

end
