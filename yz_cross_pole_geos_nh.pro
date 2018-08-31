
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
ledmn=10
leddy=31
ledyr=2007
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS Version '
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
; Read GEOS-4 data
;
      dum=findfile(dir4+sdate+'_1200.V01.nc3')
      if dum(0) ne '' then begin
         rd_geos5_nc3_meto,dir4+sdate+'_1200.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,markold2,sf2,vp2,iflag

         ncid=ncdf_open(dir4+sdate+'_1200.V01.nc4')
         mark2=fltarr(nr,nc,nth)
         ncdf_varget,ncid,3,mark2
         ncdf_close,ncid
;
; read new vortex
;
         ncid=ncdf_open(dir4+sdate+'_1200.V01.nc5')
         marknew2=fltarr(nr,nc,nth)
         ncdf_varget,ncid,3,marknew2
         ncdf_close,ncid
         tlab='GEOS-4'
      endif
;
; if GEOS-4 data is not available, then read GEOS-5
;
      if dum(0) eq '' then begin
         dum5=findfile(dir5+sdate+'_1200.V01.nc3')
         if dum5(0) ne '' then begin
            tlab='GEOS-5'
            rd_geos5_nc3_meto,dir5+sdate+'_1200.V01.nc3',nc,nr,nth,alon,alat,th,$
                     pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
            if iflag eq 1 then goto,jump
         endif
      endif
      print,sdate,' ',max(marknew2)
;
; 3d theta and latitude arrays
;
      th2=0.*u2
      alat2=0.*u2
      for i=0L,nc-1L do $
          for j=0L,nr-1L do th2(j,i,*)=th
      for i=0L,nc-1L do $
          for k=0L,nth-1L do alat2(*,i,k)=alat

      marknew2=smooth(marknew2,3)
      sp2=sqrt(u2^2.+v2^2.)
      index=where(u2 lt 0.)
      if index(0) ne -1L then sp2(index)=-1.*sp2(index)

; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
    t2=0.*p2
    z2=0.*p2
    for k=0,nth-1 do begin
        t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
        z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
    endfor
;
; stratopause vortex marker based on temperature
;
    smark2=0.*u2
    index=where(abs(alat2) gt 30. and th2 ge 2400. and t2 ge 260.)
    if index(0) ne -1L then smark2(index)=1.
;
; declare arrays and select surfaces
;
    if icount eq 0 then begin
;   print,th
    th0=3000.
;   read,' Enter desired theta surface ',th0
    thlev=where(th eq th0)
    thlev=thlev(0)
    sth=strcompress(string(fix(th0)),/remove_all)
    print,alon
    rlon1=181.875
    read,' Enter longitude ',rlon1
    index1=where(alon eq rlon1)
    if index1(0) eq -1 then stop,'Bad Longitude'
    ilon1=index1(0)
    rlon2=rlon1+180.
    if rlon2 gt max(alon) then rlon2=rlon2-360.
    index2=where(alon eq rlon2)
    ilon2=index2(0)
    slon1=strcompress(string(rlon1),/remove_all)+'E'
    slon2=strcompress(string(rlon2),/remove_all)+'E'
    print,'longitudes ',rlon1,rlon2
    x=fltarr(nc+1)
    x(0:nc-1)=alon(0:nc-1)
    x(nc)=alon(0)+360.
    x2d=fltarr(nc+1,nr)
    y2d=fltarr(nc+1,nr)
    for i=0,nc do y2d(i,*)=alat
    for j=0,nr-1 do x2d(*,j)=x
    xyz=fltarr(nr,nth)
    yyz=fltarr(nr,nth)
    for i=0,nr-1 do yyz(i,*)=th
    for j=0,nth-1 do xyz(*,j)=alat 
    endif
    u1=transpose(u2(*,*,thlev))
    v1=transpose(v2(*,*,thlev))
    pv1=transpose(pv2(*,*,thlev))
    sp1=sqrt(u1^2.+v1^2.)
    t1=transpose(t2(*,*,thlev))
    z1=transpose(z2(*,*,thlev))
    qdf1=transpose(qdf2(*,*,thlev))
    l1=transpose(mark2(*,*,thlev))
    lnew1=transpose(marknew2(*,*,thlev))
    h1=transpose(mark2(*,*,thlev))
    smark1=transpose(smark2(*,*,thlev))
;
; introduce relative vorticity
;
    zeta1=u1*0.0
    compvort_v2,u1,v1,zeta1,alon,alat,nc,nr
    zeta=fltarr(nc+1,nr)
    zeta(0:nc-1,0:nr-1)=zeta1(0:nc-1,0:nr-1)*1.e5
    zeta(nc,*)=zeta(0,*)
    qdf=fltarr(nc+1,nr)
    qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
    qdf(nc,*)=qdf(0,*)
    u=fltarr(nc+1,nr)
    u(0:nc-1,0:nr-1)=u1(0:nc-1,0:nr-1)
    u(nc,*)=u(0,*)
    v=fltarr(nc+1,nr)
    v(0:nc-1,0:nr-1)=v1(0:nc-1,0:nr-1)
    v(nc,*)=v(0,*)
    pv=fltarr(nc+1,nr)
    pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
    pv(nc,*)=pv(0,*)
    sp=fltarr(nc+1,nr)
    sp(0:nc-1,0:nr-1)=sp1(0:nc-1,0:nr-1)
    sp(nc,*)=sp(0,*)
    index=where(u lt 0.)
    if index(0) ne -1L then sp(index)=-1.*sp(index)
    t=fltarr(nc+1,nr)
    t(0:nc-1,0:nr-1)=t1(0:nc-1,0:nr-1)
    t(nc,*)=t(0,*)
    z=fltarr(nc+1,nr)
    z(0:nc-1,0:nr-1)=z1(0:nc-1,0:nr-1)
    z(nc,*)=z(0,*)
    l=fltarr(nc+1,nr)
    l(0:nc-1,0:nr-1)=l1(0:nc-1,0:nr-1)
    l(nc,*)=l(0,*)
    lnew=fltarr(nc+1,nr)
    lnew(0:nc-1,0:nr-1)=lnew1(0:nc-1,0:nr-1)
    lnew(nc,*)=lnew(0,*)
    h=fltarr(nc+1,nr)
    h(0:nc-1,0:nr-1)=h1(0:nc-1,0:nr-1)
    h(nc,*)=h(0,*)
    smark=fltarr(nc+1,nr)
    smark(0:nc-1,0:nr-1)=smark1(0:nc-1,0:nr-1)
    smark(nc,*)=smark(0,*)

    uyz=fltarr(nr,nth)
    vyz=fltarr(nr,nth)
    qyz=fltarr(nr,nth)
    spyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    zyz=fltarr(nr,nth)
    lyz=fltarr(nr,nth)
    lnewyz=fltarr(nr,nth)
    hyz=fltarr(nr,nth)
    smarkyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
        uyz(0:nr/2-1,k)=u2(nr/2:nr-1,ilon1,k)
        uyz(nr/2:nr-1,k)=reverse(u2(nr/2:nr-1,ilon2,k))
        vyz(0:nr/2-1,k)=v2(nr/2:nr-1,ilon1,k)
        vyz(nr/2:nr-1,k)=reverse(v2(nr/2:nr-1,ilon2,k))
        qyz(0:nr/2-1,k)=q2(nr/2:nr-1,ilon1,k)
        qyz(nr/2:nr-1,k)=reverse(q2(nr/2:nr-1,ilon2,k))
        spyz(0:nr/2-1,k)=sp2(nr/2:nr-1,ilon1,k)
        spyz(nr/2:nr-1,k)=reverse(sp2(nr/2:nr-1,ilon2,k))
        tyz(0:nr/2-1,k)=t2(nr/2:nr-1,ilon1,k)
        tyz(nr/2:nr-1,k)=reverse(t2(nr/2:nr-1,ilon2,k))
        zyz(0:nr/2-1,k)=z2(nr/2:nr-1,ilon1,k)
        zyz(nr/2:nr-1,k)=reverse(z2(nr/2:nr-1,ilon2,k))
        lyz(0:nr/2-1,k)=mark2(nr/2:nr-1,ilon1,k)
        lyz(nr/2:nr-1,k)=reverse(mark2(nr/2:nr-1,ilon2,k))
        lnewyz(0:nr/2-1,k)=marknew2(nr/2:nr-1,ilon1,k)
        lnewyz(nr/2:nr-1,k)=reverse(marknew2(nr/2:nr-1,ilon2,k))
        hyz(0:nr/2-1,k)=mark2(nr/2:nr-1,ilon1,k)
        hyz(nr/2:nr-1,k)=reverse(mark2(nr/2:nr-1,ilon2,k))
        smarkyz(0:nr/2-1,k)=smark2(nr/2:nr-1,ilon1,k)
        smarkyz(nr/2:nr-1,k)=reverse(smark2(nr/2:nr-1,ilon2,k))
    endfor

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='Figures/yz_cross_pole_geos_nh_'+sdate+'_'+sth+'_'+slon1+'.ps'
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
    mtitle=sth+' K '+tlab+' PV on '+sdate
    !type=2^2+2^3
    MAP_SET,90,0,-90,/stereo,/noeras,title=mtitle,color=lc
;    level=-100.+10.*findgen(nlvls)
;    contour,sp,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
;    index=where(level gt 0.)
;    contour,sp,x,alat,levels=level(index),/overplot,/follow,color=0
;    index=where(level lt 0.)
;    contour,sp,x,alat,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
    index=where(pv ne 1.12 and y2d gt 0.)
    pmin=0.
    pmax=max(pv(index))
pmax=0.3
    if icount eq 0L then pvlevel=pmin+((pmax-pmin)/nlvls)*findgen(nlvls)
    contour,pv,x,alat,levels=pvlevel,/overplot,/cell_fill,c_color=col1
    index=where(pvlevel gt 0.)
    if index(0) ne -1L then contour,pv,x,alat,levels=pvlevel(index),/overplot,/follow,color=0
    index=where(pvlevel lt 0.)
    if index(0) ne -1L then contour,pv,x,alat,levels=pvlevel(index),/overplot,/follow,color=icolmax,c_linestyle=5

    tlevel=190.+5.*findgen(nlvls)
;   for ii=0,nlvls-1L do $
;   contour,t,x,alat,levels=tlevel(ii),/overplot,/follow,c_charsize=1.5,$
;           c_labels=[1],c_color=col1(ii),thick=2
;   drawvectors,nc+1,nr,x,alat,u,v,3,0
    loadct,0
    contour,l,x,alat,levels=[0.1],color=icolmax*.5,thick=5,/overplot
    contour,lnew,x,alat,levels=[0.1],color=icolmax*.1,thick=4,/overplot
;   contour,h,x,alat,levels=[-0.1],color=icolmax*.9,thick=4,/overplot
    loadct,38
;   contour,smark,x,alat,levels=[0.1],color=icolmax*.9,thick=4,/overplot
;index=where(zeta gt 0. and qdf lt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=0,psym=8
;index=where(zeta lt 0. and qdf lt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=icolmax,psym=8
;index=where(qdf gt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=.9*icolmax,psym=8
    MAP_SET,90,0,-90,/stereo,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*alat,alat,thick=2,color=icolmax
    oplot,rlon2+0.*alat,alat,thick=2,color=icolmax

    ipan=1
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    level=-150.+15.*findgen(nlvls)
    contour,tyz,alat,th,levels=tlevel,/fill,/cell_fill,c_color=col1,$
            xtitle=slon1+'       Latitude      '+slon2,ytitle='Theta',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,title='Temp+Vortex'
    index=where(tlevel ge 250.)
    contour,tyz,alat,th,levels=tlevel(index),/overplot,/follow,color=0
    index=where(tlevel lt 250.)
    contour,tyz,alat,th,levels=tlevel(index),/overplot,/follow,color=icolmax,c_linestyle=5
;   for ii=1,nlvls-1L do $
;   contour,tyz,alat,th,levels=tlevel(ii),/follow,c_color=col1(ii),/overplot,$
;           c_labels=[1],thick=3
    loadct,0
    contour,lyz,alat,th,levels=[0.1],color=icolmax*.5,thick=8,/overplot
    contour,lnewyz,alat,th,levels=[0.1],color=icolmax*.1,thick=8,/overplot
;   contour,hyz,alat,th,levels=[-0.1],color=icolmax*.9,thick=8,/overplot
    loadct,38
    contour,smarkyz,alat,th,levels=[0.1],color=icolmax*.9,thick=8,/overplot
;   oplot,alat,th0+0.*alat,thick=3,color=icolmax
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
    mtitle=sth+' K '+tlab+' Isotachs+Vortex on '+sdate
    !type=2^2+2^3
    MAP_SET,90,0,-90,/stereo,/noeras,title=mtitle,color=lc
;   level=-100.+10.*findgen(nlvls)
    contour,sp,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    index=where(level gt 0.)
    contour,sp,x,alat,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,sp,x,alat,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
;   tlevel=190.+5.*findgen(nlvls)
;   for ii=0,nlvls-1L do $
;   contour,t,x,alat,levels=tlevel(ii),/overplot,/follow,c_charsize=1.5,$
;           c_labels=[1],c_color=col1(ii),thick=3
;   drawvectors,nc+1,nr,x,alat,u,v,3,0
    loadct,0
    contour,l,x,alat,levels=[0.1],color=icolmax*.5,thick=5,/overplot
    contour,lnew,x,alat,levels=[0.1],color=icolmax*.1,thick=4,/overplot
;   contour,h,x,alat,levels=[-0.1],color=icolmax*.9,thick=4,/overplot
    loadct,38
    contour,smark,x,alat,levels=[0.1],color=icolmax*.9,thick=4,/overplot
;index=where(zeta gt 0. and qdf lt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=0,psym=8
;index=where(zeta lt 0. and qdf lt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=icolmax,psym=8
;index=where(qdf gt 0. and y2d gt 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),color=.9*icolmax,psym=8
    MAP_SET,90,0,-90,/stereo,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*alat,alat,thick=2,color=icolmax
    oplot,rlon2+0.*alat,alat,thick=2,color=icolmax

    ipan=3
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;   level=-100.+10.*findgen(nlvls)
    contour,spyz,alat,th,levels=level,/fill,/cell_fill,c_color=col1,$
            xtitle=slon1+'       Latitude      '+slon2,ytitle='Theta',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,title='Wind Speed'
    index=where(level gt 0.)
    contour,spyz,alat,th,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,spyz,alat,th,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
;   for ii=12,nlvls-1L do $
;   contour,tyz,alat,th,levels=tlevel(ii),/follow,c_color=col1(ii),/overplot,$
;           c_labels=[1],thick=3
    loadct,0
    contour,lyz,alat,th,levels=[0.1],color=icolmax*.5,thick=8,/overplot
    contour,lnewyz,alat,th,levels=[0.1],color=icolmax*.1,thick=8,/overplot
;   contour,hyz,alat,th,levels=[-0.1],color=icolmax*.9,thick=8,/overplot
    loadct,38
;   contour,smarkyz,alat,th,levels=[0.1],color=icolmax*.9,thick=8,/overplot
;   oplot,alat,th0+0.*alat,thick=3,color=icolmax
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
        spawn,'convert -trim Figures/yz_cross_pole_geos_nh_'+sdate+'_'+sth+'_'+slon1+'.ps -rotate -90 '+$
                            'Figures/yz_cross_pole_geos_nh_'+sdate+'_'+sth+'_'+slon1+'.jpg'
        spawn,'/usr/bin/rm Figures/yz_cross_pole_geos_nh_'+sdate+'_'+sth+'_'+slon1+'.ps'
     endif

goto,jump
end
