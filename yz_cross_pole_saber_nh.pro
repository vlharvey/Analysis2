
; GEOS-4 or GEOS-5 version, depending on date
; plot polar projection and yz cross polar section

@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto
@compvort_v2	; w/4th order differentiation across the poles
@drawvectors

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

loadct,38
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.10,0.10,0.55,0.55]
yorig=[0.55,0.10,0.55,0.10]
xlen=0.4
ylen=0.4
cbaryoff=0.06
cbarydel=0.01
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir4='/aura7/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
dir5='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS501.MetO.'

lstmn=10
lstdy=1
lstyr=2004
ledmn=12
leddy=31
ledyr=2005
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      SABER Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; read SABER winds
;
    sfile='GRID_PHI_WINDS.'+sdate+'.sav'
    sdum=findfile('/aura6/data/SABER_data/Datfiles_winds/'+sfile)
    if sdum(0) eq '' then goto,jump
    restore,sdum(0)
    tdata=t3d   ; avoid T3D intrinsic function
    if max(tdata) eq 0. then goto,jump
    udata=u3d
    vdata=v3d
    slat=alat
    slon=alon
    spress=press
    snc=n_elements(slon)
    snr=n_elements(slat)
    snl=n_elements(spress)
    ubar_saber=fltarr(snr,snl)
    tbar_saber=fltarr(snr,snl)
    for k=0L,snl-1L do begin
    for j=0L,snr-1L do begin
        index=where(udata(*,j,k) ne 0.,ngood)
        if index(0) ne -1L then begin
           ubar_saber(j,k)=total(udata(index,j,k))/float(ngood)
           tbar_saber(j,k)=total(tdata(index,j,k))/float(ngood)
        endif
    endfor
    endfor
;
; 3d theta and latitude arrays
;
;     th2=0.*u2
;     alat2=0.*u2
;     for i=0L,nc-1L do $
;         for j=0L,nr-1L do th2(j,i,*)=th
;     for i=0L,nc-1L do $
;         for k=0L,nth-1L do alat2(*,i,k)=alat

;     marknew2=smooth(marknew2,3)
;     sp2=sqrt(u2^2.+v2^2.)
;     index=where(u2 lt 0.)
;     if index(0) ne -1L then sp2(index)=-1.*sp2(index)

; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
;   t2=0.*p2
;   z2=0.*p2
;   for k=0,nth-1 do begin
;       t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
;       z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
;   endfor
;
; stratopause vortex marker based on temperature
;
;   smark2=0.*u2
;   index=where(abs(alat2) gt 30. and th2 ge 2400. and t2 ge 260.)
;   if index(0) ne -1L then smark2(index)=1.
;
; declare arrays and select surfaces
;
    if icount eq 0 then begin
;   print,spress
    p0=1.11377
    p0=0.101039
;   read,' Enter desired pressure surface ',p0
    plev=where(abs(spress-p0) eq min(abs(spress-p0)))
    plev=plev(0)
    slev=strcompress(p0,/remove_all)
;   print,slon
    rlon1=180.
;   read,' Enter longitude ',rlon1
    index1=where(slon eq rlon1)
    if index1(0) eq -1 then stop,'Bad Longitude'
    ilon1=index1(0)
    rlon2=rlon1+180.
    if rlon2 gt max(slon) then rlon2=rlon2-360.
    index2=where(slon eq rlon2)
    ilon2=index2(0)
    slon1=strcompress(string(rlon1),/remove_all)+'E'
    slon2=strcompress(string(rlon2),/remove_all)+'E'
    print,'longitudes ',rlon1,rlon2
    x=fltarr(snc+1)
    x(0:snc-1)=slon(0:snc-1)
    x(snc)=slon(0)+360.
    x2d=fltarr(snc+1,snr)
    y2d=fltarr(snc+1,snr)
    for i=0,snc do y2d(i,*)=slat
    for j=0,snr-1 do x2d(*,j)=x
    xyz=fltarr(snr,snl)
    yyz=fltarr(snr,snl)
    for i=0,snr-1 do yyz(i,*)=spress
    for j=0,snl-1 do xyz(*,j)=slat 
    endif
    u1=udata(*,*,plev)
    v1=vdata(*,*,plev)
    t1=tdata(*,*,plev)
    u=fltarr(snc+1,snr)
    u(0:snc-1,0:snr-1)=u1(0:snc-1,0:snr-1)
    u(snc,*)=u(0,*)
    v=fltarr(snc+1,snr)
    v(0:snc-1,0:snr-1)=v1(0:snc-1,0:snr-1)
    v(snc,*)=v(0,*)
    t=fltarr(snc+1,snr)
    t(0:snc-1,0:snr-1)=t1(0:snc-1,0:snr-1)
    t(snc,*)=t(0,*)
    sp=sqrt(u^2.+v^2.)
    index=where(u lt 0.)
    if index(0) ne -1L then sp(index)=-1.*sp(index)
    t=fltarr(snc+1,snr)
    t(0:snc-1,0:snr-1)=t1(0:snc-1,0:snr-1)
    t(snc,*)=t(0,*)

    uyz=fltarr(snr+1,snl)
    vyz=fltarr(snr+1,snl)
    tyz=fltarr(snr+1,snl)
    for k=0,snl-1 do begin
        uyz(0:snr/2,k)=udata(ilon1,snr/2:snr-1,k)
        dum=reform(udata(ilon2,snr/2:snr-1,k))
        uyz(snr/2+1:snr,k)=reverse(dum)
        vyz(0:snr/2,k)=vdata(ilon1,snr/2:snr-1,k)
        dum=reform(vdata(ilon2,snr/2:snr-1,k))
        vyz(snr/2+1:snr,k)=reverse(dum)
        tyz(0:snr/2,k)=tdata(ilon1,snr/2:snr-1,k)
        dum=reform(tdata(ilon2,snr/2:snr-1,k))
        tyz(snr/2+1:snr,k)=reverse(dum)
    endfor
    spyz=sqrt(uyz^2.+vyz^2.)
    index=where(uyz lt 0.)
    if index(0) ne -1L then spyz(index)=-1.*spyz(index)

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='Figures/yz_cross_pole_saber_nh_'+sdate+'_'+slev+'_'+slon1+'.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

; plot
    !noeras=1
    !type=2^2+2^3
    !p.thick=1
    erase
    ipan=0
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    mtitle=slev+' K SABER Temperature on '+sdate
    !type=2^2+2^3
    MAP_SET,90,0,-90,/stereo,/noeras,title=mtitle,color=lc
    tlevel=190.+5.*findgen(nlvls)
    contour,t,x,slat,levels=tlevel,/overplot,/cell_fill,c_color=col1
    index=where(tlevel gt 0.)
    contour,t,x,slat,levels=tlevel(index),/overplot,/follow,color=0
;   drawvectors,snc+1,snr,x,slat,u,v,3,0
    MAP_SET,90,0,-90,/stereo,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*alat,alat,thick=2,color=icolmax
    oplot,rlon2+0.*alat,alat,thick=2,color=icolmax

    ipan=1
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,tyz,findgen(snr+1),spress,/ylog,yrange=[max(spress),0.0001],levels=tlevel,/fill,/cell_fill,c_color=col1,$
            xtitle=slon1+'       Latitude      '+slon2,ytitle='Pressure (hPa)',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
    index=where(tlevel gt 0.)
    contour,tyz,findgen(snr+1),spress,levels=tlevel(index),/overplot,/follow,color=0
    imin=min(tlevel)
    imax=max(tlevel)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    ipan=2
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    mtitle=slev+' hPa SABER Isotachs on '+sdate
    !type=2^2+2^3
    MAP_SET,90,0,-90,/stereo,/noeras,title=mtitle,color=lc
    level=-100.+10.*findgen(nlvls)
    contour,sp,x,slat,levels=level,/overplot,/cell_fill,c_color=col1
    index=where(level gt 0.)
    contour,sp,x,slat,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,sp,x,slat,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
    drawvectors,snc+1,snr,x,slat,u,v,3,0
    MAP_SET,90,0,-90,/stereo,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*slat,slat,thick=2,color=icolmax
    oplot,rlon2+0.*slat,slat,thick=2,color=icolmax

    ipan=3
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,spyz,findgen(snr+1),spress,/ylog,yrange=[max(spress),0.0001],levels=level,/fill,/cell_fill,c_color=col1,$
            xtitle=slon1+'       Latitude      '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
    index=where(level gt 0.)
    contour,spyz,findgen(snr+1),spress,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,spyz,findgen(snr+1),spress,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    icount=icount+1

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim Figures/yz_cross_pole_saber_nh_'+sdate+'_'+slev+'_'+slon1+'.ps -rotate -90 '+$
                            'Figures/yz_cross_pole_saber_nh_'+sdate+'_'+slev+'_'+slon1+'.jpg'
        spawn,'/usr/bin/rm Figures/yz_cross_pole_saber_nh_'+sdate+'_'+slev+'_'+slon1+'.ps'
     endif

goto,jump
end
