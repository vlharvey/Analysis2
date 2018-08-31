;
; compute daily zonal mean Z, T, U, V, and mark from Matt Rigby's model data
; plot latitude-time section at a user specified altitude
; save entire 3-D (y,z,t) arrays for time period
;
@kdate
@rd_ukmo_nc3
@rd_rigby_dat
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
dir='/aura2/harvey/Rigby/Datfiles/ensem_03.MetO.'
goto,plotit
spawn,'ls '+dir+'*.dat',pfiles
spawn,'ls '+dir+'*.nc3',ncfiles
icount=0L
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
;
; read pressure data
;
    sday=string(FORMAT='(I3.3)',ifile)
    rd_rigby_dat,pfiles(ifile),iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,zp,tp,up,vp
    glon=alon
    glat=alat
    nlg=n_elements(glon)
    nlat=n_elements(glat)
    nlv=n_elements(p)
;
; compute daily zonal means in Z, T, U, V
;
    zbar=fltarr(nlat,nlv)
    tbar=fltarr(nlat,nlv)
    ubar=fltarr(nlat,nlv)
    vbar=fltarr(nlat,nlv)
    for k=0L,nlv-1L do begin
    for j=0L,nlat-1L do begin
        zbar(j,k)=total(zp(*,j,k))/float(nlg)
        tbar(j,k)=total(tp(*,j,k))/float(nlg)
        ubar(j,k)=total(up(*,j,k))/float(nlg)
        vbar(j,k)=total(vp(*,j,k))/float(nlg)
    endfor
    endfor
;
; read theta data
;
    rd_geos5_nc3_meto,ncfiles(ifile),nc,nr,nth,glonth,glatth,th,$
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
      maxmarkbar=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      y2d(*,k)=glatth
      for j=0L,nr-1L do begin
          index=where(p2(j,*,k) ne 0.,npt)
          if index(0) ne -1L then begin
             pbar(j,k)=total(p2(j,index,k))/float(npt)
             markbar(j,k)=total(mark2(j,index,k))/float(npt)
             maxmarkbar(j,k)=max(mark2(j,*,k))
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
         markbar_yt=fltarr(nr,nth,nfile)
         maxmarkbar_yt=fltarr(nr,nth,nfile)
         pbar_yt=fltarr(nr,nth,nfile)
         sdate_all=strarr(nfile)
         icount=1
      endif
      sdate_all(ifile)=sday
;
; save zonal means on each day
;
      zbar_yt(*,*,ifile)=zbar
      tbar_yt(*,*,ifile)=tbar
      ubar_yt(*,*,ifile)=ubar
      vbar_yt(*,*,ifile)=vbar
      markbar_yt(*,*,ifile)=markbar
      maxmarkbar_yt(*,*,ifile)=maxmarkbar
      pbar_yt(*,*,ifile)=pbar
    jump:
endfor          ; loop over time steps

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='yt_rigby_ztuv+edge_'+sdate0+'-'+sdate1+'.sav',zbar_yt,tbar_yt,ubar_yt,$
     vbar_yt,markbar_yt,maxmarkbar_yt,pbar_yt,sdate_all,glat,p,glatth,th,nfile

plotit:
restore,'yt_rigby_ztuv+edge_000-299.sav
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
       device,/landscape,bits=8,filename='yt_rigby_ztuv+edge_'+sdate0+'-'+sdate1+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; 1, 0.1, 0.01 hPa levels.  3 panel
; =p(11), p(17), and p(22)
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
    sdy=strmid(sdate_all,1,2)
    xindex=where(sdy eq '00',nxticks)
    xlabs=sdate_all(xindex)

    ilev=22L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=xlabs
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    contour,tbar,findgen(nfile),glat,levels=[200.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[210.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    markbar=transpose(reform(markbar_yt(*,5,*)))      ; ~0.02 hPa or ~4000 K is top continuous level
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=10,c_labels=[0]
    maxmarkbar=transpose(reform(maxmarkbar_yt(*,5,*)))
    contour,maxmarkbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=8,c_labels=[0]
loadct,38
print,min(tbar),max(tbar)
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
    xyouts,xmx-0.2,ymx+0.01,'200 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'210 K',/normal,color=mcolor*.9,charthick=2,charsize=2

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=17L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=xlabs
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    contour,tbar,findgen(nfile),glat,levels=[225.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[240.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    markbar=transpose(reform(markbar_yt(*,10,*)))      ; ~0.1 hPa or ~3000 K
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=10,c_labels=[0]
    maxmarkbar=transpose(reform(maxmarkbar_yt(*,10,*)))
    contour,maxmarkbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=8,c_labels=[0]
loadct,38
    xyouts,xmx-0.2,ymx+0.01,'225 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'240 K',/normal,color=mcolor*.9,charthick=2,charsize=2
print,min(tbar),max(tbar)
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=11L
    ubar=transpose(reform(ubar_yt(*,ilev,*)))
    zbar=transpose(reform(zbar_yt(*,ilev,*)))/1000.
    tbar=smooth(transpose(reform(tbar_yt(*,ilev,*))),3)
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[0.,nfile-1],charthick=1.5,yticks=6,xticks=nxticks,xtickname=xlabs
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=[40,60,80,100,120],color=0
    contour,ubar,findgen(nfile),glat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40],color=mcolor,c_linestyle=5
    contour,tbar,findgen(nfile),glat,levels=[240.],color=mcolor*.2,/follow,/noeras,/overplot,min_value=0.,thick=5
    contour,tbar,findgen(nfile),glat,levels=[270.],color=mcolor*.9,/follow,/noeras,/overplot,min_value=0.,thick=5
loadct,0
    markbar=transpose(reform(markbar_yt(*,15,*)))      ; ~1 hPa or ~2000 K
    contour,markbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.3,/follow,/noeras,/overplot,min_value=0.,thick=10,c_labels=[0]
    maxmarkbar=transpose(reform(maxmarkbar_yt(*,15,*)))
    contour,maxmarkbar,findgen(nfile),glatth,levels=[0.5],color=mcolor*.5,/follow,/noeras,/overplot,min_value=0.,thick=8,c_labels=[0]
loadct,38
    xyouts,xmx-0.2,ymx+0.01,'240 K',/normal,color=mcolor*.2,charthick=2,charsize=2
    xyouts,xmx-0.1,ymx+0.01,'270 K',/normal,color=mcolor*.9,charthick=2,charsize=2
print,min(tbar),max(tbar)
    result=moment(zbar)
    zlab=string(long(result(0)))
    xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
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
       spawn,'convert -trim yt_rigby_ztuv+edge_'+sdate0+'-'+sdate1+'.ps -rotate -90 yt_rigby_ztuv+edge_'+sdate0+'-'+sdate1+'.jpg'
    endif
end
