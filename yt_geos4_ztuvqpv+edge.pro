;
; compute daily zonal mean GEOS-4 Z, T, U, V, Q, PV, mark for time period
; plot latitude-time section at a user specified altitude
; save entire 3-D (y,z,t) arrays for time period
;
@kdate
@rd_ukmo_nc3
@rd_geos5_dat
@rd_geos5_nc3_meto

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
dir='/aura7/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
gdir='/aura7/harvey/GEOS4_data/Datfiles/'
;dir='/Users/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
;gdir='/Users/harvey/GEOS4_data/Datfiles/'
goto,plotit
spawn,'ls '+gdir+'DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.origp.*_1200.V01.sav',ncfiles
icount=0L
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
;
; compute daily zonal means in Z, T, U, V, Q, PV
;
    zbar=fltarr(nlat,nlv)
    tbar=fltarr(nlat,nlv)
    ubar=fltarr(nlat,nlv)
    vbar=fltarr(nlat,nlv)
    qbar=fltarr(nlat,nlv)
    pvbar=fltarr(nlat,nlv)
    for k=0L,nlv-1L do begin
    for j=0L,nlat-1L do begin
        zbar(j,k)=total(znew(*,j,k))/float(nlg)
        tbar(j,k)=total(tnew(*,j,k))/float(nlg)
        ubar(j,k)=total(unew(*,j,k))/float(nlg)
        vbar(j,k)=total(vnew(*,j,k))/float(nlg)
        qbar(j,k)=total(qnew(*,j,k))/float(nlg)
        pvbar(j,k)=total(pvnew(*,j,k))/float(nlg)
    endfor
    endfor
;
; read GEOS4 marker
;
    result=strsplit(ncfiles(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(7),'_',/extract)
    sdate=result3(0)
    rd_geos5_nc3_meto,dir+sdate+'_1200.V01.nc3',nc,nr,nth,glonth,glatth,th,$
                 pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
    if iflag ne 0 then goto, jump
;
; remove anticyclones
;
    index=where(mark2 lt 0.)
    if index(0) ne -1L then mark2(index)=0.
;
; zonal mean max mark
;
      y2d=fltarr(nr,nth)
      pbar=fltarr(nr,nth)
      markbar=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      y2d(*,k)=glatth
      for j=0L,nr-1L do begin
          index=where(p2(j,*,k) ne 0.,npt)
          if index(0) ne -1L then begin
             pbar(j,k)=total(p2(j,index,k))/float(npt)
             markbar(j,k)=total(mark2(j,index,k))/float(npt)
;            markbar(j,k)=max(mark2(j,*,k))
          endif
      endfor
      endfor
;     markbar=smooth(markbar,3)
;
; declare time period arrays on first day
;
      if icount eq 0L then begin
         zbar_yt=fltarr(nlat,nlv,nfile)
         tbar_yt=fltarr(nlat,nlv,nfile)
         ubar_yt=fltarr(nlat,nlv,nfile)
         vbar_yt=fltarr(nlat,nlv,nfile)
         qbar_yt=fltarr(nlat,nlv,nfile)
         pvbar_yt=fltarr(nlat,nlv,nfile)
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
save,file='yt_geos4_ztuvqpv_'+sdate0+'-'+sdate1+'.sav',zbar_yt,tbar_yt,ubar_yt,$
     vbar_yt,qbar_yt,pvbar_yt,markbar_yt,pbar_yt,sdate_all,glat,p,glatth,th,nfile

plotit:
restore,'yt_geos4_ztuvqpv_20040331-20070515.sav
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
       device,/landscape,bits=8,filename='yt_geos4_tuvqpv+edge_'+sdate0+'-'+sdate1+'.ps'
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

    ilev=0L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=smooth(transpose(reform(markbar_yt(*,5,*))),5)	; ~0.02 hPa or ~4000 K is top continuous level
    contour,tbar,findgen(nfile),glat,levels=[170.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[230.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=5,c_labels=[0]
loadct,38
print,min(tbar),max(tbar)
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
    xyouts,xmx-0.2,ymx+0.01,'170 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'230 K',/normal,color=mcolor*.9,charthick=2,charsize=2
;
; add vortex edge, PV contours or negative regions (save again- not zonal mean), height?, date labels
;
    for ii=0L,nxticks-1L do xyouts,xindex(ii),-105.,smon(long(xlabs(ii))-1),/data,color=0    

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=5L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=smooth(transpose(reform(markbar_yt(*,10,*))),5)  ; ~0.1 hPa or ~3000 K
    contour,tbar,findgen(nfile),glat,levels=[220.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[270.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=5,c_labels=[0]
loadct,38
    xyouts,xmx-0.2,ymx+0.01,'220 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'270 K',/normal,color=mcolor*.9,charthick=2,charsize=2
print,min(tbar),max(tbar)
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
    for ii=0L,nxticks-1L do xyouts,xindex(ii),-105.,smon(long(xlabs(ii))-1),/data,color=0

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=14L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=' '+strarr(nxticks+1)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    markbar=smooth(transpose(reform(markbar_yt(*,15,*))),5)  ; ~1 hPa or ~2000 K
    contour,tbar,findgen(nfile),glat,levels=[240.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[270.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=5,c_labels=[0]
loadct,38
    xyouts,xmx-0.2,ymx+0.01,'240 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'270 K',/normal,color=mcolor*.9,charthick=2,charsize=2
print,min(tbar),max(tbar)
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
       spawn,'convert yt_geos4_tuvqpv+edge_'+sdate0+'-'+sdate1+'.ps -rotate -90 yt_geos4_tuvqpv+edge_'+sdate0+'-'+sdate1+'.jpg'
    endif
end
