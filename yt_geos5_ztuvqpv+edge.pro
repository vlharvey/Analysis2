;
; compute daily zonal mean GEOS-5 Z, T, U, V, Q, PV, mark for time period
; plot latitude-time section at a user specified altitude
; save entire 3-D (y,z,t) arrays for time period
;
@kdate
@rd_geos5_nc3_meto

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
xorig=[0.15,0.15,0.15]
yorig=[0.725,0.425,0.125]
xlen=0.7
ylen=0.225
cbaryoff=0.08
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
goto,plotit
spawn,'ls '+dir+'*_AVG.V01.nc3',ncfiles
icount=0L
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
;
; construct date string
;
    result=strsplit(ncfiles(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(6),'_',/extract)
    sdate=result3(0)
    syr=strmid(sdate,0,4)
    smn=strmid(sdate,4,2)
    sdy=strmid(sdate,6,2)
    sdate=syr+smn+sdy
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
;
; remove anticyclones
;
    index=where(mark2 lt 0.)
    if index(0) ne -1L then mark2(index)=0.
;
; compute temperature and geopotential height
;
    t2=0.*pv2
    z2=0.*pv2
    for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
    z2=(msf2 - 1004.*t2)/(9.86*1000.)
;
; compute daily zonal means in Z, T, U, V, Q, PV
;
    pbar=fltarr(nr,nth)
    zbar=fltarr(nr,nth)
    tbar=fltarr(nr,nth)
    ubar=fltarr(nr,nth)
    vbar=fltarr(nr,nth)
    qbar=fltarr(nr,nth)
    pvbar=fltarr(nr,nth)
    markbar=fltarr(nr,nth)
    for k=0L,nth-1L do begin
    for j=0L,nr-1L do begin
if min(p2(j,*,k)) eq 0. then stop
        zbar(j,k)=total(z2(j,*,k))/float(nc)
        tbar(j,k)=total(t2(j,*,k))/float(nc)
        ubar(j,k)=total(u2(j,*,k))/float(nc)
        vbar(j,k)=total(v2(j,*,k))/float(nc)
        qbar(j,k)=total(q2(j,*,k))/float(nc)
        pvbar(j,k)=total(pv2(j,*,k))/float(nc)
        pbar(j,k)=total(p2(j,*,k))/float(nc)
        markbar(j,k)=max(mark2(j,*,k))
    endfor
    endfor
;   markbar=smooth(markbar,3)
;
; declare time period arrays on first day
;
      if icount eq 0L then begin
         zbar_yt=fltarr(nr,nth,nfile)
         tbar_yt=fltarr(nr,nth,nfile)
         ubar_yt=fltarr(nr,nth,nfile)
         vbar_yt=fltarr(nr,nth,nfile)
         qbar_yt=fltarr(nr,nth,nfile)
         pvbar_yt=fltarr(nr,nth,nfile)
         markbar_yt=fltarr(nr,nth,nfile)
         pbar_yt=fltarr(nr,nth,nfile)
         sdate_all=strarr(nfile)
         icount=1
      endif
      sdate_all(ifile)=sdate
;
; save zonal means on each day
;
      zbar_yt(*,*,ifile)=zbar
      tbar_yt(*,*,ifile)=tbar
      ubar_yt(*,*,ifile)=ubar
      vbar_yt(*,*,ifile)=vbar
      qbar_yt(*,*,ifile)=qbar
      pvbar_yt(*,*,ifile)=pvbar
      markbar_yt(*,*,ifile)=markbar
      pbar_yt(*,*,ifile)=pbar

    jump:
endfor          ; loop over time steps

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='yt_geos5_ztuvqpv_'+sdate0+'-'+sdate1+'.sav',zbar_yt,tbar_yt,ubar_yt,$
     vbar_yt,qbar_yt,pvbar_yt,markbar_yt,pbar_yt,sdate_all,alat,alat,th,nfile

plotit:
restore,'yt_geos5_ztuvqpv_20031001-20080902.sav
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)

    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yt_geos5_tuvqpv+edge_'+sdate0+'-'+sdate1+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; 1, 0.1, 0.01 hPa levels.  3 panel
; =p(14), p(5), and p(0)
;
    erase
    level=-120.+10.*findgen(25)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*mcolor/nlvls
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    syr=strmid(sdate_all,0,4)
    smn=strmid(sdate_all,4,2)
    sdy=strmid(sdate_all,6,2)
    xindex=where(sdy eq '15',nxticks)
    xlabs=smn(xindex)

    ilev=5L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(th(ilev))+' (K)',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=transpose(reform(markbar_yt(*,ilev,*)))
;   contour,tbar,findgen(nfile),alat,levels=[170.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
;   contour,tbar,findgen(nfile),alat,levels=[230.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
;   contour,markbar,findgen(nfile),alat,levels=[1],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=3,c_labels=[0]
loadct,39
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
;   xyouts,xmx-0.2,ymx+0.01,'170 K',/normal,color=mcolor*.2,charthick=2,charsize=2
;   xyouts,xmx-0.1,ymx+0.01,'230 K',/normal,color=mcolor*.9,charthick=2,charsize=2
;
; add vortex edge, PV contours or negative regions (save again- not zonal mean), height?, date labels
;
    for ii=0L,nxticks-1L do xyouts,xindex(ii),-105.,smon(long(xlabs(ii))-1),/data,color=0    

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=12L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(th(ilev))+' (K)',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=transpose(reform(markbar_yt(*,ilev,*)))
;   contour,tbar,findgen(nfile),alat,levels=[220.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
;   contour,tbar,findgen(nfile),alat,levels=[270.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
;   contour,markbar,findgen(nfile),alat,levels=[1],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=3,c_labels=[0]
loadct,39
;   xyouts,xmx-0.2,ymx+0.01,'220 K',/normal,color=mcolor*.2,charthick=2,charsize=2
;   xyouts,xmx-0.1,ymx+0.01,'270 K',/normal,color=mcolor*.9,charthick=2,charsize=2
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
    for ii=0L,nxticks-1L do xyouts,xindex(ii),-105.,smon(long(xlabs(ii))-1),/data,color=0

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=15L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))
    pvbar=transpose(reform(pvbar_yt(*,ilev,*)))
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(th(ilev))+' (K)',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=transpose(reform(markbar_yt(*,ilev,*)))
;   contour,tbar,findgen(nfile),alat,levels=[240.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
;   contour,tbar,findgen(nfile),alat,levels=[270.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
;   contour,markbar,findgen(nfile),alat,levels=[1],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=3,c_labels=[0]
;   contour,pvbar,findgen(nfile),alat,levels=[-0.005,0.005],color=mcolor*.2,/follow,/noeras,$
;          /overplot,min_value=0.,thick=3,c_labels=[0]
;   contour,pvbar,findgen(nfile),alat,levels=[-0.002,0.002],color=mcolor*.5,/follow,/noeras,$
;          /overplot,min_value=0.,thick=3,c_labels=[0]
;   contour,pvbar,findgen(nfile),alat,levels=[-0.001,0.001],color=mcolor*.7,/follow,/noeras,$
;          /overplot,min_value=0.,thick=3,c_labels=[0]
loadct,39
;   xyouts,xmx-0.2,ymx+0.01,'240 K',/normal,color=mcolor*.2,charthick=2,charsize=2
;   xyouts,xmx-0.1,ymx+0.01,'270 K',/normal,color=mcolor*.9,charthick=2,charsize=2
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
    for ii=0L,nxticks-1L do xyouts,xindex(ii),-105.,smon(long(xlabs(ii))-1),/data,color=0
    index=where(smn eq '07' and sdy eq '01',nyear)
    for ii=0L,nyear-1L do xyouts,index(ii),-130.,syr(index(ii)),/data,color=0,charsize=2
    xyouts,nfile-60.,-130.,syr(nfile-1),/data,color=0,charsize=2
    imin=min(level)
    imax=max(level)
    ymnb=yorig(2)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(m/s)'
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert yt_geos5_tuvqpv+edge_'+sdate0+'-'+sdate1+'.ps -rotate -90 yt_geos5_tuvqpv+edge_'+sdate0+'-'+sdate1+'.jpg'
    endif
end
