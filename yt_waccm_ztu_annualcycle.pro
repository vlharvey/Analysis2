;
; compute daily zonal mean WACCM Z, T, U 
; plot latitude-time section at a user specified altitude
; save entire 3-D (y,z,t) arrays for time period
;
@kgmt

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
yorig=[0.73,0.42,0.11]
xlen=0.7
ylen=0.225
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
dir='/aura7/harvey/WACCM_data/Datfiles/Aurora/waccm319_8_smax'
goto,plotit
spawn,'ls '+dir+'*h3*.sav',ncfiles
icount=0L
gcount=0L
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
;
; read WACCM daily data
;
; DATE            LONG      =     19950101
; GHGT            FLOAT     = Array[72, 46, 66]
; LATITUDE        DOUBLE    = Array[46]
; LONGITUDE       DOUBLE    = Array[72]
; PRESSURE        DOUBLE    = Array[66]
; PSFC            FLOAT     = Array[72, 46]
; TEMP            FLOAT     = Array[72, 46, 66]
; TIME            DOUBLE    =        0.0000000
; UWIND           FLOAT     = Array[72, 46, 66]
; VWIND           FLOAT     = Array[72, 46, 66]
;
    restore,ncfiles(ifile)
    print,ncfiles(ifile)
    sdate=strcompress(date,/remove_all)
;
; rename variables for convenience
;
    alon=longitude
    alat=latitude
    press=pressure
    tdata=temp
    if max(tdata) eq 0. then goto,jump
    udata=uwind
    zdata=ghgt/1000.
;
; compute zonal mean T, U, Z
;
    nlat=n_elements(alat)
    nlg=n_elements(alon)
    nlv=n_elements(press)
    p=press
    ubar=fltarr(nlat,nlv)
    tbar=fltarr(nlat,nlv)
    zbar=fltarr(nlat,nlv)
    for k=0L,nlv-1L do begin
    for j=0L,nlat-1L do begin
        index=where(udata(*,j,k) ne 0.,ngood)
        if index(0) ne -1L then begin
           ubar(j,k)=total(udata(index,j,k))/float(ngood)
           tbar(j,k)=total(tdata(index,j,k))/float(ngood)
           zbar(j,k)=total(zdata(index,j,k))/float(ngood)
        endif
    endfor
    endfor
;
; declare time period arrays on first day
;
      if icount eq 0L then begin
         zbar_yt=fltarr(nlat,nlv,nfile)
         tbar_yt=fltarr(nlat,nlv,nfile)
         ubar_yt=fltarr(nlat,nlv,nfile)
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

    jump:
endfor          ; loop over time steps

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='yt_waccm_ztu_'+sdate0+'-'+sdate1+'.sav',zbar_yt,tbar_yt,ubar_yt,sdate_all,alat,p,nfile

plotit:
restore,'yt_waccm_ztu_19950101-20250131.sav
nday=365L
nr=n_elements(alat)
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
fdoy=fltarr(n_elements(sdate_all))
for ii=0L,n_elements(sdate_all)-1L do begin
    if syear(ii) ne '' then begin
       iyr=long(syear(ii))
       imn=long(smon(ii))
       idy=long(sday(ii))
       z = kgmt(imn,idy,iyr,iday)
       fdoy(ii)=float(iday)
       if syear(ii) gt syear(0) then fdoy(ii)=float(iday)
    endif
    if syear(ii) eq '' then if fdoy(ii-1) ne 0. then fdoy(ii)=fdoy(ii-1)+1.
endfor
;
; interpolate small gaps in time
;
index=where(ubar_yt eq 0.)
if index(0) ne -1L then begin
   ubar_yt(index)=0./0.
   tbar_yt(index)=0./0.
endif
ubar_yt=smooth(ubar_yt,3,/NaN)
tbar_yt=smooth(tbar_yt,3,/NaN)
index=where(finite(tbar_yt) ne 1L)
if index(0) ne -1L then tbar_yt(index)=-999.
index=where(ubar_yt eq 0.)
if index(0) ne -1L then ubar_yt(index)=-999.
;
; save postscript file
;
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yt_waccm_ztu_annualcycle.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

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
;
; create "annual cycle" arrays dimensioned 365 days by latitude
; at 1, 0.1, 0.01 hPa levels.  3 panel
; pressure(26),p(20),p(16) = 0.95594750, 0.12827325, 0.017898000
;
    ilev=16L
    ubar_annual_yt=fltarr(nday,nr)
    for i=0,nday-1L do begin
        iday=i+1.
        today=where(fdoy eq iday)
        ubar=reform(ubar_yt(*,ilev,today))
        for j=0,nr-1L do begin
            index=where(ubar(j,*) ne -999.,nn)
            if index(0) ne -1L then ubar_annual_yt(i,j)=total(ubar(j,index))/float(nn)
        endfor
    endfor
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[1.,365.],charthick=1.5,yticks=6,xticks=11,xtickname=smonth,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=[40,80,120],color=0,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,80,40],color=mcolor,c_linestyle=5,min_value=-999
;   index=where(zbar gt 0.)
;   result=moment(zbar(index))
;   zlab=string(long(result(0)))
;   xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
;   xyouts,xmx-0.2,ymx+0.01,'170 K',/normal,color=mcolor*.2,charthick=2,charsize=2
;   xyouts,xmx-0.1,ymx+0.01,'230 K',/normal,color=mcolor*.9,charthick=2,charsize=2
;
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=20L
    ubar_annual_yt=fltarr(nday,nr)
    for i=0,nday-1L do begin
        iday=i+1.
        today=where(fdoy eq iday)
        ubar=reform(ubar_yt(*,ilev,today))
        for j=0,nr-1L do begin
            index=where(ubar(j,*) ne -999.,nn)
            if index(0) ne -1L then ubar_annual_yt(i,j)=total(ubar(j,index))/float(nn)
        endfor
    endfor
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[min(fdoy),max(fdoy)],charthick=1.5,yticks=6,xticks=11,xtickname=smonth,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=[40,80,120],color=0,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,80,40],color=mcolor,c_linestyle=5,min_value=-999
;   index=where(zbar gt 0.)
;   result=moment(zbar(index))
;   zlab=string(long(result(0)))
;   xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    ilev=26L
    ubar_annual_yt=fltarr(nday,nr)
    for i=0,nday-1L do begin
        iday=i+1.
        today=where(fdoy eq iday)
        ubar=reform(ubar_yt(*,ilev,today))
        for j=0,nr-1L do begin
            index=where(ubar(j,*) ne -999.,nn)
            if index(0) ne -1L then ubar_annual_yt(i,j)=total(ubar(j,index))/float(nn)
        endfor
    endfor
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,yrange=[-90.,90.],$
            ytitle='Latitude',levels=level,c_color=col1,/cell_fill,color=0,title=string(p(ilev))+' hPa',$
            xrange=[min(fdoy),max(fdoy)],charthick=1.5,yticks=6,xticks=11,xtickname=smonth,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=[40,80,120],color=0,min_value=-999
    contour,ubar_annual_yt,1.+findgen(nday),alat,charsize=1.5,/noeras,/overplot,levels=-1.*[120,80,40],color=mcolor,c_linestyle=5,min_value=-999
;   index=where(zbar gt 0.)
;   result=moment(zbar(index))
;   zlab=string(long(result(0)))
;   xyouts,xmx-0.08,(ymx+ymn)/2.,zlab+' km',/normal,color=0,charthick=2,charsize=2
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
       spawn,'convert -trim yt_waccm_ztu_annualcycle.ps -rotate -90 yt_waccm_ztu_annualcycle.jpg'
    endif
end
