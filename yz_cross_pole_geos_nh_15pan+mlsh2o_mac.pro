;
; Mac version
;
; plot polar projections and yz cross polar sections

@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto
@compvort_v2
@drawvectors

sver='v2.2'

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

loadct,39
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=1
!p.charsize=1
!p.charthick=2
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.07,0.255,0.44,0.625,0.81,0.07,0.255,0.44,0.625,0.81,0.07,0.255,0.44,0.625,0.81]
yorig=[0.7,0.7,0.7,0.7,0.7,0.425,0.425,0.425,0.425,0.425,0.15,0.15,0.15,0.15,0.15]
xlen=0.175
ylen=0.175
cbaryoff=0.02
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/Volumes/atmos/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/Volumes/atmos/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
dirm='/Volumes/atmos/aura6/data/MLS_data/Datfiles_SOSST/'

lstmn=12
lstdy=1
lstyr=2008
ledmn=12
leddy=10
ledyr=2008
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
;
; oct through march only
;
      if iday gt 90 and iday lt 274 then goto,jumpday

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
      dum=findfile(dir+sdate+'_AVG.V01.nc3')
      if dum(0) ne '' then begin
         rd_geos5_nc3_meto,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag

;        ncid=ncdf_open(dir+sdate+'_AVG.V01.nc4')
;        mark2=fltarr(nr,nc,nth)
;        ncdf_varget,ncid,3,mark2
;        ncdf_close,ncid
      endif
      print,sdate
;
; 3d theta and latitude arrays
;
      th2=0.*u2
      alat2=0.*u2
      for i=0L,nc-1L do $
          for j=0L,nr-1L do th2(j,i,*)=th
      for i=0L,nc-1L do $
          for k=0L,nth-1L do alat2(*,i,k)=alat

      sf2=sf2*1.e6
      sp2=sqrt(u2^2.+v2^2.)
      index=where(u2 lt 0.)
      if index(0) ne -1L then sp2(index)=-1.*sp2(index)
;
; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
;
      t2=0.*p2
      z2=0.*p2
      for k=0,nth-1 do begin
          t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
          z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
      endfor
;
; dtheta/dz = (th1-th0) / (z1-z0)
;
      dthdz2=0.*t2
      for k=0,nth-1 do begin
          lp1=k-1
          lm1=k+1
          IF K EQ 0 then LP1=0
          IF K EQ NTH-1 then LM1=NTH-1
          for i=0,nc-1 do begin
          for j=0,nr-1 do begin
              DTHDz2(j,i,K)=(TH(LP1)-TH(LM1))/(z2(j,i,LP1)-z2(j,i,LM1))
          endfor
          endfor
      endfor
;
; declare arrays and select surfaces
;
      if icount eq 0 then begin
;        print,th
         th0=2200.
         th1=1000.
;        read,' Enter desired theta surface ',th0
         thlev0=where(th eq th0)
         thlev0=thlev0(0)
         thlev1=where(th eq th1)
         thlev1=thlev1(0)
         sth0=strcompress(string(fix(th0)),/remove_all)
         sth1=strcompress(string(fix(th1)),/remove_all)
;        print,alon
         rlon1=181.875
;        read,' Enter longitude ',rlon1
         index1=where(alon eq rlon1)
         if index1(0) eq -1 then stop,'Bad Longitude'
         ilon1=index1(0)
         rlon2=rlon1+180.
         if rlon2 gt max(alon) then rlon2=rlon2-360.
         index2=where(alon eq rlon2)
         ilon2=index2(0)
         slon1=string(format='(f7.3)',rlon1)+'E'
         slon2=string(format='(f7.3)',rlon2)+'E'
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
;
; compute relative vorticity
;
    zeta2=u2*0.0
    zeta1=fltarr(nc,nr)
    for k=0L,nth-1L do begin
        u1=transpose(u2(*,*,k))
        v1=transpose(v2(*,*,k))
        compvort_v2,u1,v1,zeta1,alon,alat,nc,nr
        zeta2(*,*,k)=transpose(zeta1)*1.e5
    endfor
;
; extract 3600 K level
;
    dthdz0=transpose(dthdz2(*,*,thlev0))
    u0=transpose(u2(*,*,thlev0))
    v0=transpose(v2(*,*,thlev0))
    pv0=transpose(pv2(*,*,thlev0))
    sp0=transpose(sp2(*,*,thlev0))
    t0=transpose(t2(*,*,thlev0))
    z0=transpose(z2(*,*,thlev0))
    qdf0=transpose(qdf2(*,*,thlev0))
    mark0=smooth(transpose(mark2(*,*,thlev0)),3)
    sf0=transpose(sf2(*,*,thlev0))
    zeta0=transpose(zeta2(*,*,thlev0))

    zeta00=fltarr(nc+1,nr)
    zeta00(0:nc-1,0:nr-1)=zeta0(0:nc-1,0:nr-1)
    zeta00(nc,*)=zeta00(0,*)
    dthdz00=fltarr(nc+1,nr)
    dthdz00(0:nc-1,0:nr-1)=dthdz0(0:nc-1,0:nr-1)
    dthdz00(nc,*)=dthdz00(0,*)
    pv00=fltarr(nc+1,nr)
    pv00(0:nc-1,0:nr-1)=pv0(0:nc-1,0:nr-1)
    pv00(nc,*)=pv00(0,*)
    sp00=fltarr(nc+1,nr)
    sp00(0:nc-1,0:nr-1)=sp0(0:nc-1,0:nr-1)
    sp00(nc,*)=sp00(0,*)
    t00=fltarr(nc+1,nr)
    t00(0:nc-1,0:nr-1)=t0(0:nc-1,0:nr-1)
    t00(nc,*)=t00(0,*)
    mark00=fltarr(nc+1,nr)
    mark00(0:nc-1,0:nr-1)=mark0(0:nc-1,0:nr-1)
    mark00(nc,*)=mark00(0,*)
    sf00=fltarr(nc+1,nr)
    sf00(0:nc-1,0:nr-1)=sf0(0:nc-1,0:nr-1)
    sf00(nc,*)=sf00(0,*)
;
; extract 1000 K level
;
    dthdz1=transpose(dthdz2(*,*,thlev1))
    pv1=transpose(pv2(*,*,thlev1))
    sp1=transpose(sp2(*,*,thlev1))
    t1=transpose(t2(*,*,thlev1))
    mark1=smooth(transpose(mark2(*,*,thlev1)),3)
    sf1=transpose(sf2(*,*,thlev1))
    zeta1=transpose(zeta2(*,*,thlev1))

    zeta11=fltarr(nc+1,nr)
    zeta11(0:nc-1,0:nr-1)=zeta1(0:nc-1,0:nr-1)
    zeta11(nc,*)=zeta11(0,*)
    dthdz11=fltarr(nc+1,nr)
    dthdz11(0:nc-1,0:nr-1)=dthdz1(0:nc-1,0:nr-1)
    dthdz11(nc,*)=dthdz11(0,*)
    pv11=fltarr(nc+1,nr)
    pv11(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
    pv11(nc,*)=pv11(0,*)
    sp11=fltarr(nc+1,nr)
    sp11(0:nc-1,0:nr-1)=sp1(0:nc-1,0:nr-1)
    sp11(nc,*)=sp11(0,*)
    t11=fltarr(nc+1,nr)
    t11(0:nc-1,0:nr-1)=t1(0:nc-1,0:nr-1)
    t11(nc,*)=t11(0,*)
    mark11=fltarr(nc+1,nr)
    mark11(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
    mark11(nc,*)=mark11(0,*)
    sf11=fltarr(nc+1,nr)
    sf11(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
    sf11(nc,*)=sf11(0,*)

    dthdzyz=fltarr(nr,nth)
    zetayz=fltarr(nr,nth)
    spyz=fltarr(nr,nth)
    sfyz=fltarr(nr,nth)
    pvyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    markyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
        zetayz(0:nr/2-1,k)=zeta2(nr/2:nr-1,ilon1,k)
        zetayz(nr/2:nr-1,k)=reverse(zeta2(nr/2:nr-1,ilon2,k))
        dthdzyz(0:nr/2-1,k)=dthdz2(nr/2:nr-1,ilon1,k)
        dthdzyz(nr/2:nr-1,k)=reverse(dthdz2(nr/2:nr-1,ilon2,k))
        sfyz(0:nr/2-1,k)=sf2(nr/2:nr-1,ilon1,k)
        sfyz(nr/2:nr-1,k)=reverse(sf2(nr/2:nr-1,ilon2,k))
        spyz(0:nr/2-1,k)=sp2(nr/2:nr-1,ilon1,k)
        spyz(nr/2:nr-1,k)=reverse(sp2(nr/2:nr-1,ilon2,k))

        result=moment(pv2(nr/2:nr-1,*,k))
; convert to Lait's version for vertical slice
;       pvyz(0:nr/2-1,k)=1.e4*result(0)*((th(k)/1000.))^(-9./2.)
        pvyz(0:nr/2-1,k)=1.e4*pv2(nr/2:nr-1,ilon1,k)*((th(k)/2000.))^(-9./2.)
        pvyz(nr/2:nr-1,k)=1.e4*reverse(pv2(nr/2:nr-1,ilon1,k))*((th(k)/2000.))^(-9./2.)

;       pvyz(0:nr/2-1,k)=pv2(nr/2:nr-1,ilon1,k)
;       pvyz(nr/2:nr-1,k)=reverse(pv2(nr/2:nr-1,ilon2,k))
        tyz(0:nr/2-1,k)=t2(nr/2:nr-1,ilon1,k)
        tyz(nr/2:nr-1,k)=reverse(t2(nr/2:nr-1,ilon2,k))
        markyz(0:nr/2-1,k)=mark2(nr/2:nr-1,ilon1,k)
        markyz(nr/2:nr-1,k)=reverse(mark2(nr/2:nr-1,ilon2,k))
    endfor
;
; read MLS CO data.  read temperature and pressure to interpolate to theta surfaces
;
    dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
    if dum(0) eq '' then goto,jumpday
    restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
    restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
    restore,dirm+'co_mls_'+sver+'_'+sdate+'.sav'              ; mix
;   restore,dirm+'h2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
;   restore,dirm+'n2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
    nz=n_elements(altitude)
    nthlev=n_elements(thlev)
    mprof=n_elements(longitude)
    mlev=n_elements(altitude)
    muttime=time
    mlat=latitude
    mlon=longitude
    bad=where(mask eq -99.)
    if bad(0) ne -1L then mix(bad)=-99.
    good=where(mix ne -99.)
    if good(0) eq -1L then goto,jump
    mco=mix*1.e6
    mtemp=temperature
    mpress=pressure
;
; eliminate bad uttimes and SH
;
    index=where(muttime gt 0. and mlat gt 20.,mprof)
    if index(0) eq -1L then goto,jump
    muttime=reform(muttime(index))
    mlat=reform(mlat(index))
    mlon=reform(mlon(index))
    mtemp=reform(mtemp(index,*))
    mpress=reform(mpress(index,*))
    mco=reform(mco(index,*))
    mtheta=mtemp*(1000./mpress)^0.286
    index=where(mtemp lt 0.)
    if index(0) ne -1L then mtheta(index)=-99.
;
; interpolate CO data to GEOS-5 theta surfaces
;
    mco_th=fltarr(mprof,nth)
    for k=nth-1L,0L,-1L do begin
        zlev=th(k)
        for iprof=0L,mprof-1L do begin
            for kk=2L,nz-2L do begin

if mco(iprof,kk) ne -9.90000e+07 and mco(iprof,kk+1) ne -9.90000e+07 then begin
if mtheta(iprof,kk) lt zlev and mtheta(iprof,kk+1) ge zlev then begin
   zscale=(mtheta(iprof,kk+1)-zlev)/(mtheta(iprof,kk+1)-mtheta(iprof,kk))
   mco_th(iprof,k)= mco(iprof,kk+1)+zscale*(mco(iprof,kk)-mco(iprof,kk+1))

;print,mtheta(iprof,kk),zlev,mtheta(iprof,kk+1),zscale
;print,mco(iprof,kk),mco_th(iprof,k),mco(iprof,kk+1)
;stop
endif
endif
            endfor
        endfor
    endfor
;
; extract swath closest to rlon1 and rlon2
;
    index=where(mlat ge 62. and mlat le 66.)
    mlon60=mlon(index)
    muttime60=muttime(index)
    index=where(abs(mlon60-rlon1) eq min(abs(mlon60-rlon1)))
    mtime1=muttime60(index(0))
    mindex=where(abs(muttime-mtime1) lt 1.)
;
; swap order of swath if it goes from right to left
;
    mlon_swath=mlon(mindex)
    if mlon_swath(0) gt 270. then mindex=reverse(mindex)
    mco_swath=reform(mco_th(mindex,*))
    mlon_swath=mlon(mindex)
    mlat_swath=mlat(mindex)
    muttime_swath=muttime(mindex)

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='yz_cross_pole_geos_nh_'+sdate+'_'+slon1+'_15pan+mlsh2o.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

; plot
    erase
    xyouts,.4,.95,'GEOS-5 '+sdate,/normal,color=0,charsize=3
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(t00)
    imax=max(t00)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Temperature',color=0
    contour,t00,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,t00,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp00,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark00,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*alat,alat,color=icolmax
    oplot,rlon2+0.*alat,alat,color=icolmax
    oplot,mlon(mindex),mlat(mindex),psym=2,symsize=0.5
    xyouts,xmn-0.01,ymn+0.04,sth0+' K',/normal,color=0,charsize=1.5,orientation=90.
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(dthdz00)
    imax=max(dthdz00)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Static Stability',color=0
    contour,dthdz00,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,dthdz00,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp00,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark00,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(zeta00)
    imax=max(zeta00)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Relative Vorticity',color=0
    contour,zeta00,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,zeta00,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp00,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark00,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(pv00)
    imax=max(pv00)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='PV',color=0
    contour,pv00,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,pv00,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp00,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark00,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(4)
    xmx=xorig(4)+xlen
    ymn=yorig(4)
    ymx=yorig(4)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(sf00)
    imax=max(sf00)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Stream Function',color=0
    contour,sf00,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,sf00,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp00,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark00,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor
;
; 1000 K
;
    !type=2^2+2^3
    xmn=xorig(5)
    xmx=xorig(5)+xlen
    ymn=yorig(5)
    ymx=yorig(5)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(t11)
    imax=max(t11)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Temperature',color=0
    contour,t11,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,t11,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp11,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark11,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    xyouts,xmn-0.01,ymn+0.04,sth1+' K',/normal,color=0,charsize=1.5,orientation=90.
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(6)
    xmx=xorig(6)+xlen
    ymn=yorig(6)
    ymx=yorig(6)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(dthdz11)
    imax=max(dthdz11)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Static Stability',color=0
    contour,dthdz11,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,dthdz11,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp11,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark11,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(7)
    xmx=xorig(7)+xlen
    ymn=yorig(7)
    ymx=yorig(7)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(zeta11)
    imax=max(zeta11)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Relative Vorticity',color=0
    contour,zeta11,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,zeta11,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp11,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark11,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(8)
    xmx=xorig(8)+xlen
    ymn=yorig(8)
    ymx=yorig(8)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(pv11)
    imax=max(pv11)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='PV',color=0
    contour,pv11,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,pv11,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp11,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark11,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(9)
    xmx=xorig(9)+xlen
    ymn=yorig(9)
    ymx=yorig(9)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=min(sf11)
    imax=max(sf11)
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    MAP_SET,90,0,-90,/ortho,/noeras,title='Stream Function',color=0
    contour,sf11,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    contour,sf11,x,alat,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,sp11,x,alat,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,mark11,x,alat,levels=[0.1],color=0,thick=5,/overplot
    MAP_SET,90,0,-90,/ortho,/grid,/contin,/noeras,/noborder,color=0
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(10)
    xmx=xorig(10)+xlen
    ymn=yorig(10)
    ymx=yorig(10)+ylen
    set_viewport,xmn,xmx,ymn,ymx
index=where(yyz lt 4000)
    imin=min(tyz(index))
    imax=max(tyz(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,tyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='Temperature',$
            xtitle=slon1+'       '+slon2,ytitle='Theta',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[500.,4000.]
    contour,tyz,alat,th,levels=level,/overplot,/follow,color=0
    contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(11)
    xmx=xorig(11)+xlen
    ymn=yorig(11)
    ymx=yorig(11)+ylen
    set_viewport,xmn,xmx,ymn,ymx
index=where(yyz lt 4000)
    imin=min(dthdzyz(index))
    imax=max(dthdzyz(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,dthdzyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='Static Stability',$
            xtitle=slon1+'       '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[500.,4000.]
    contour,dthdzyz,alat,th,levels=level,/overplot,/follow,color=0
    contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(12)
    xmx=xorig(12)+xlen
    ymn=yorig(12)
    ymx=yorig(12)+ylen
    set_viewport,xmn,xmx,ymn,ymx
index=where(yyz lt 4000)
    imin=min(zetayz(index))
    imax=max(zetayz(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,zetayz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='Relative Vorticity',$
            xtitle=slon1+'       '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[500.,4000.]
    contour,zetayz,alat,th,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(13)
    xmx=xorig(13)+xlen
    ymn=yorig(13)
    ymx=yorig(13)+ylen
    set_viewport,xmn,xmx,ymn,ymx
index=where(yyz lt 4000)
    imin=min(pvyz(index))
    imax=max(pvyz(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,pvyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='Laits PV',$
            xtitle=slon1+'       '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[500.,4000.]
    contour,pvyz,alat,th,levels=level(0:nlvls-1:3),/overplot,/follow,color=0
    contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    !type=2^2+2^3
    xmn=xorig(14)
    xmx=xorig(14)+xlen
    ymn=yorig(14)
    ymx=yorig(14)+ylen
    set_viewport,xmn,xmx,ymn,ymx
index=where(yyz lt 4000)
    imin=min(sfyz(index))
    imax=max(sfyz(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,sfyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='SF+H2O',$
            xtitle=slon1+'       '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[500.,4000.]
    contour,sfyz,alat,th,levels=level,/overplot,/follow,color=0
;   contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
    contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot

    axis,xrange=[0.,n_elements(mlon_swath)-1L],yrange=[500.,4000.],/save,/data
    mco_swath=smooth(mco_swath,9)
;
; CO
;
;   contour,mco_swath,findgen(n_elements(mlon_swath)),th,levels=0.05*findgen(30),/follow,c_color=mcolor,/overplot
;
; water
;
    contour,mco_swath,findgen(n_elements(mlon_swath)),th,levels=2.+0.5*findgen(30),/follow,c_color=mcolor,/overplot,thick=5
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
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
        spawn,'convert -trim yz_cross_pole_geos_nh_'+sdate+'_'+slon1+'_15pan+mlsh2o.ps -rotate -90 '+$
                            'yz_cross_pole_geos_nh_'+sdate+'_'+slon1+'_15pan+mlsh2o.jpg'
        spawn,'/usr/bin/rm yz_cross_pole_geos_nh_'+sdate+'_'+slon1+'_15pan+mlsh2o.ps'
     endif
     jumpday:
goto,jump
end
