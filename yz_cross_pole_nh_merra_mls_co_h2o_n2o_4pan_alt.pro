;
; hemispheric anomalies
; add MLS water and temp
; add N2O
;
; MERRA version
; plot polar projections and yz cross polar sections
;
; plot on altitude

@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

sver='v3.3'

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
!NOERAS=-1
!P.FONT=1
!p.charsize=1
!p.charthick=2
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.55,0.1,0.55]
yorig=[0.55,0.55,0.15,0.15]
xlen=0.325
ylen=0.325
cbaryoff=0.07
cbarydel=0.02
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirm2='/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_'

lstmn=1
lstdy=1
lstyr=2013
ledmn=2
leddy=28
ledyr=2013
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      MERRA Version '
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
; Read MERRA
;
      dum=findfile(dir+sdate+'.nc3')
      if dum(0) ne '' then begin
        rd_merra_nc3,dum(0),nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
;        ncid=ncdf_open(dir+sdate+'.nc4')
;        mark2=fltarr(nr,nc,nth)
;        ncdf_varget,ncid,3,mark2
;        ncdf_close,ncid
      endif
      index=where(mark2 lt 0.)
      if index(0) ne -1L then mark2(index)=-1.*mark2(index)/mark2(index)
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
;     if index(0) ne -1L then sp2(index)=-1.*sp2(index)
;
; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
;
      t2=0.*p2
      for k=0,nth-1 do begin
          t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
      endfor
;
; declare arrays and select surfaces
;
      if icount eq 0 then begin
;        print,th
         th0=3000.
         th1=1000.
;        read,' Enter desired theta surface ',th0
         thlev0=where(th eq th0)
         thlev0=thlev0(0)
         thlev1=where(th eq th1)
         thlev1=thlev1(0)
         sth0=strcompress(string(fix(th0)),/remove_all)
         sth1=strcompress(string(fix(th1)),/remove_all)
;        print,alon
         rlon1=150.
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
; restore MLS gridded data
;
    restore,dirm2+sver+'_'+sdate+'.sav'
    restore,'/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_U_V_'+sver+'_'+sdate+'.sav'
    index=where(finite(u) eq 1 and finite(v) eq 1)
    mlsspeed=0.*u
    if index(0) ne -1L then mlsspeed(index)=sqrt(u(index)^2.+v(index)^2.)
;
; make hemispheric anomalies
;
    h2o_anom=0.*h2o_grid
    co_anom=0.*co_grid
    n2o_anom=0.*co_grid
    for k=0,n_elements(pmls2)-1L do begin
        h2olev=reform(h2o_grid(*,nr/2:nr-1,k))
        index=where(finite(h2olev) eq 1)
        if index(0) ne -1L then begin
           h2oavg=mean(h2olev(index))
           h2o_anom(*,*,k)=h2o_grid(*,*,k)	;-h2oavg
;print,pmls2(k),h2oavg
        endif
    endfor
    for k=0,n_elements(pmls)-1L do begin
        colev=reform(co_grid(*,nr/2:nr-1,k))
        index=where(finite(colev) eq 1)
        if index(0) ne -1L then begin
           coavg=mean(colev(index))
           co_anom(*,*,k)=co_grid(*,*,k)	;-coavg
;print,pmls(k),coavg
        endif
    endfor
    for k=0,n_elements(pmls)-1L do begin
        n2olev=reform(n2o_grid(*,nr/2:nr-1,k))
        index=where(finite(n2olev) eq 1)
        if index(0) ne -1L then begin
           n2oavg=mean(n2olev(index))
           n2o_anom(*,*,k)=n2o_grid(*,*,k)	;-n2oavg
;print,pmls(k),n2oavg
        endif
    endfor
;
; cross polar sections
;
    nz=n_elements(pmls2)
    comlsyz=fltarr(nr,n_elements(pmls))
    n2omlsyz=fltarr(nr,n_elements(pmls))
    h2omlsyz=fltarr(nr,nz)
    tmlsyz=fltarr(nr,nz)
    spmlsyz=fltarr(nr,nz)
    gpmlsyz=fltarr(nr,nz)
    for k=0,nz-1L do begin
        tmlsyz(0:nr/2-1,k)=tp_grid(ilon1,nr/2:nr-1,k)
        tmlsyz(nr/2:nr-1,k)=reverse(tp_grid(ilon2,nr/2:nr-1,k))
        tmlsyz(nr/2:nr-1,k)=reverse(tmlsyz(nr/2:nr-1,k))

        spmlsyz(0:nr/2-1,k)=mlsspeed(ilon1,nr/2:nr-1,k)
        spmlsyz(nr/2:nr-1,k)=reverse(mlsspeed(ilon2,nr/2:nr-1,k))
        spmlsyz(nr/2:nr-1,k)=reverse(spmlsyz(nr/2:nr-1,k))

        h2omlsyz(0:nr/2-1,k)=h2o_anom(ilon1,nr/2:nr-1,k)
        h2omlsyz(nr/2:nr-1,k)=reverse(h2o_anom(ilon2,nr/2:nr-1,k))
        h2omlsyz(nr/2:nr-1,k)=reverse(h2omlsyz(nr/2:nr-1,k))

        gpmlsyz(0:nr/2-1,k)=gp_grid(ilon1,nr/2:nr-1,k)/1000.
        gpmlsyz(nr/2:nr-1,k)=reverse(gp_grid(ilon2,nr/2:nr-1,k))/1000.
    endfor
    for k=0,n_elements(pmls)-1L do begin
        comlsyz(0:nr/2-1,k)=co_anom(ilon1,nr/2:nr-1,k)
        comlsyz(nr/2:nr-1,k)=reverse(co_anom(ilon2,nr/2:nr-1,k))
        comlsyz(nr/2:nr-1,k)=reverse(comlsyz(nr/2:nr-1,k))

        n2omlsyz(0:nr/2-1,k)=n2o_anom(ilon1,nr/2:nr-1,k)
        n2omlsyz(nr/2:nr-1,k)=reverse(n2o_anom(ilon2,nr/2:nr-1,k))
        n2omlsyz(nr/2:nr-1,k)=reverse(n2omlsyz(nr/2:nr-1,k))
    endfor

    spyz=fltarr(nr,nth)
    sfyz=fltarr(nr,nth)
    pvyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    zyz=fltarr(nr,nth)
    markyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
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
        zyz(0:nr/2-1,k)=z2(nr/2:nr-1,ilon1,k)
        zyz(nr/2:nr-1,k)=reverse(z2(nr/2:nr-1,ilon2,k))
        tyz(0:nr/2-1,k)=t2(nr/2:nr-1,ilon1,k)
        tyz(nr/2:nr-1,k)=reverse(t2(nr/2:nr-1,ilon2,k))
        markyz(0:nr/2-1,k)=mark2(nr/2:nr-1,ilon1,k)
        markyz(nr/2:nr-1,k)=reverse(mark2(nr/2:nr-1,ilon2,k))
    endfor
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='yz_cross_pole_nh_merra_mls_temp_'+sdate+'_'+slon1+'_4pan_alt.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot
;
    erase
    xyouts,.4,.95,sdate,/normal,color=0,charsize=3
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    imin=160.
    imax=280.
    int=7.5
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
    level=imin+int*findgen(nlvls)
;   tyz=smooth(tyz,9,/edge_truncate)
    contour,tmlsyz,alat,gpmlsyz,levels=level,/cell_fill,c_color=col1,color=0,title='MLS Temp',$
            ytitle='Altitude (km)',xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[15,100],charsize=1.5
    contour,smooth(spmlsyz,3,/edge_truncate,/Nan),alat,gpmlsyz,levels=50.+20.*findgen(10),/overplot,/follow,color=mcolor,thick=2,c_labels=['']
print,'MLS ',max(spmlsyz(10:85,10:45))
    contour,tmlsyz,alat,gpmlsyz,levels=level,/overplot,/follow,color=0
    contour,markyz,alat,zyz,levels=[0.1],color=0,thick=5,/overplot
      imin=min(level)
      imax=max(level)
      xmnb=xorig(0)+xlen+cbaryoff
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,yorig(0)+cbarydel,yorig(0)+ylen-cbarydel
      !type=2^2+2^3+2^5+2^7
      plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K',color=0,xticks=4,charsize=1.5
      xbox=[0,10,10,0,0]
      y1=imin
      dy=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor

    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;
; common levels where zflag=1
;
    zflag=fltarr(n_elements(pmls2))
    for k=0L,n_elements(pmls2)-1L do begin
        index=where(pmls2(k) eq pmls)
        if index(0) ne -1L then zflag(k)=1.
    endfor
    kindex=where(zflag eq 1.)
    level=-1+0.1*findgen(21)
level=[0.001,0.01+0.01*findgen(20)]
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
    contour,n2omlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level,/cell_fill,c_color=col1,color=0,title='MLS Gridded N2O',$
            xtitle=slon1+'       '+slon2,ytitle='Altitude (km)',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[15,100],charsize=1.5
    index=where(level gt 0)
    contour,n2omlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0)
    contour,n2omlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level(index),/overplot,/follow,color=mcolor,c_linestyle=5
    contour,n2omlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=[0],/overplot,/follow,color=0,thick=5
    contour,smooth(markyz,3,/edge_truncate),alat,zyz,levels=[0.1],color=mcolor,thick=5,/overplot
      imin=min(level)
      imax=max(level)
      xmnb=xorig(1)+xlen+cbaryoff
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,yorig(1)+cbarydel,yorig(1)+ylen-cbarydel
      !type=2^2+2^3+2^5+2^7
      plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K',color=0,xticks=4,charsize=1.5
      xbox=[0,10,10,0,0]
      y1=imin
      dy=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor

    !type=2^2+2^3
    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;
; common levels where zflag=1
;
    zflag=fltarr(n_elements(pmls2))
    for k=0L,n_elements(pmls2)-1L do begin
        index=where(pmls2(k) eq pmls)
        if index(0) ne -1L then zflag(k)=1.
    endfor
    kindex=where(zflag eq 1.)
    level=-5+0.5*findgen(21)
level=[0.001,0.01,0.1+0.1+findgen(19)]

nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
    contour,comlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level,/cell_fill,c_color=col1,color=0,title='MLS Gridded CO',$
            xtitle=slon1+'       '+slon2,ytitle='Altitude (km)',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[15,100],charsize=1.5
    index=where(level gt 0)
    contour,comlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0)
    contour,comlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=level(index),/overplot,/follow,color=mcolor,c_linestyle=5
    contour,comlsyz*1.e6,alat,gpmlsyz(*,kindex),levels=[0],/overplot,/follow,color=0,thick=5
    contour,smooth(markyz,3,/edge_truncate),alat,zyz,levels=[0.1],color=mcolor,thick=5,/overplot
    imin=min(level)
    imax=max(level)
    xmnb=xorig(2)+xlen+cbaryoff
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,yorig(2)+cbarydel,yorig(2)+ylen-cbarydel
    !type=2^2+2^3+2^5+2^7
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='ppmv',color=0,xticks=4,charsize=1.5
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    !type=2^2+2^3
    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;
; set up plot bounds with no data
;
;nlvls=19
;col1=1+indgen(nlvls)*icolmax/nlvls
;    imin=160.
;    imax=280.
;    int=7.5
;    level=imin+int*findgen(nlvls)
;    contour,mtemp_swath_alt,findgen(n_elements(mlon_swath)),altitude,levels=level,/cell_fill,c_color=col1,color=0,title='MLS Profile Temp',$
;            xtitle=slon1+'       '+slon2,$	;ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
;            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[15,80],/nodata,charsize=1.5
;    contour,mtemp_swath_alt,findgen(n_elements(mlon_swath)),altitude,levels=level,/cell_fill,c_color=col1,/overplot
;    contour,mtemp_swath_alt,findgen(n_elements(mlon_swath)),altitude,levels=level,/follow,color=0,/overplot
    level=-2.5+0.5*findgen(11)
level=0.1+0.4*findgen(21)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
    contour,h2omlsyz*1.e6,alat,gpmlsyz,levels=level,/cell_fill,c_color=col1,color=0,title='MLS Gridded H2O',$
;           xtitle=slon1+'       '+slon2,ytitle='Altitude (km)',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[15,100],charsize=1.5
    index=where(level gt 0)
    contour,h2omlsyz*1.e6,alat,gpmlsyz,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0)
    contour,h2omlsyz*1.e6,alat,gpmlsyz,levels=level(index),/overplot,/follow,color=mcolor,c_linestyle=5
    contour,h2omlsyz*1.e6,alat,gpmlsyz,levels=[0],/overplot,/follow,color=0,thick=5
    contour,smooth(markyz,3,/edge_truncate),alat,zyz,levels=[0.1],color=mcolor,thick=5,/overplot

      imin=min(level)
      imax=max(level)
      xmnb=xorig(3)+xlen+cbaryoff
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,yorig(3)+cbarydel,yorig(3)+ylen-cbarydel
      !type=2^2+2^3+2^5+2^7
      plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='ppmv',color=0,xticks=4,charsize=1.5
      xbox=[0,10,10,0,0]
      y1=imin
      dy=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor

    icount=icount+1

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim yz_cross_pole_nh_merra_mls_temp_'+sdate+'_'+slon1+'_4pan_alt.ps -rotate -90 '+$
                            'yz_cross_pole_nh_merra_mls_temp_'+sdate+'_'+slon1+'_4pan_alt.jpg'
        spawn,'/usr/bin/rm yz_cross_pole_nh_merra_mls_temp_'+sdate+'_'+slon1+'_4pan_alt.ps'
     endif
     jumpday:
goto,jump
end
