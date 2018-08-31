;
; plot monthly mean zonal mean WACCM and HIRDLS NO2
;
loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
yorig=[0.4,0.4,0.4]
xorig=[0.1,0.4,0.7]
ylen=0.3
xlen=0.27
cbaryoff=0.1
cbarydel=0.01
nxdim=800
nydim=800
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
!noeras=1
ver='v2.04.19'
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_TNV3/wa3_tnv3_'
hdir='/aura6/data/HIRDLS_data/Datfiles_SOSST/'
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
    no2yth_mean_waccm=fltarr(nr,nth)
    o3yth_mean_waccm=fltarr(nr,nth)
    tempyth_mean_waccm=fltarr(nr,nth)
    for k=0L,nth-1 do begin
    for j=0L,nr-1 do begin
        no2yth_mean_waccm(j,k)=total(no2(j,*,k))/float(nc)
        o3yth_mean_waccm(j,k)=total(o3(j,*,k))/float(nc)
        tempyth_mean_waccm(j,k)=total(tmp2(j,*,k))/float(nc)
    endfor
    endfor
;
; restore monthly mean HIRDLS data
;
    restore,hdir+mon(imonth)+'_avg_yth_'+ver+'.sav'	;nr,nth,alat,th,pyth_mean,no2yth_mean,o3yth_mean,tempyth_mean

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_waccm3+hirdls_'+mon(imonth)+'_no2.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
       !p.thick=2.0                   ;Plotted lines twice as thick
       !p.charsize=2.0
    endif
;
; HIRDLS NO2
;
    loadct,39
    erase
    xyouts,.4,.75,month(imonth),/normal,charsize=3,charthick=2,color=0
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    no2level=1.+findgen(18)
    nlvls=n_elements(no2level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,no2yth_mean*1.e9,alat,th,levels=no2level,c_color=col1,/noeras,/cell_fill,$
            title='HIRDLS NO!l2!n',color=0,yrange=[500.,2000.],min_value=0.,xticks=6,$
            ytitle='Potential Temperature (K)',xtitle='Latitude',xrange=[-90.,90.]
    contour,no2yth_mean*1.e9,alat,th,levels=no2level,color=0,/noeras,/follow,/overplot,$
            min_value=0.
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(no2level)
    imax=max(no2level)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppbv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; WACCM NO2
;
    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,no2yth_mean_waccm,alat,th,levels=no2level,c_color=col1,/noeras,/cell_fill,$
            title='WACCM NO!l2!n',color=0,yrange=[500.,2000.],xticks=6,xrange=[-90.,90.],$
            ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],xtitle='Latitude'
    contour,no2yth_mean_waccm,alat,th,levels=no2level,color=0,/noeras,/follow,/overplot
    set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    imin=min(no2level)
    imax=max(no2level)
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='(ppbv)',/noeras,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
;
; NO2 difference
;
    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    no2diff=0.*no2yth_mean
    index=where(no2yth_mean gt 0. and no2yth_mean_waccm gt 0.)
    if index(0) eq -1L then goto,skipmonth
    no2yth_mean(index)=no2yth_mean(index)*1.e9
    no2diff(index)=100.*(no2yth_mean(index)-no2yth_mean_waccm(index))/no2yth_mean_waccm(index)
    restore,'c11_rb.tbl'
    tvlct,c1,c2,c3
    nlvls2=11
    dlevel=-100.+20.*findgen(nlvls2)
    col2=1+indgen(nlvls2)
    contour,no2diff,alat,th,levels=dlevel,c_color=col2,/noeras,/cell_fill,$
            title='100*(H-W)/W',color=0,yrange=[500.,2000.],xticks=6,xrange=[-90.,90.],$
            ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],xtitle='Latitude'
    index=where(dlevel gt 0.)
    contour,no2diff,alat,th,levels=dlevel(index),color=0,/noeras,/follow,/overplot
    index=where(dlevel lt 0.)
    contour,no2diff,alat,th,levels=dlevel(index),color=mcolor,/noeras,/follow,/overplot,c_linestyle=1
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

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yz_waccm3+hirdls_'+mon(imonth)+'_no2.ps -rotate -90 yz_waccm3+hirdls_'+mon(imonth)+'_no2.jpg'
;      spawn,'/usr/bin/rm yz_waccm3+hirdls_'+mon(imonth)+'_no2.ps'
    endif
    skipmonth:
endfor
end
