;
; add MLS water and temp
;
; MERRA version
; plot polar projections and yz cross polar sections

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
dirm='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'

lstmn=1
lstdy=1
lstyr=2006
ledmn=3
leddy=1
ledyr=2006
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
      for k=0,nth-1 do begin
          t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
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
         rlon1=180.
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
    dthdzyz=fltarr(nr,nth)
    spyz=fltarr(nr,nth)
    sfyz=fltarr(nr,nth)
    pvyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    markyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
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
    restore,dirm+'h2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
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
; dtheta/dz = (th1-th0) / (z1-z0)
;
    mdthdz=0.*mtheta
    for i=0,mprof-1 do begin
    for k=0,mlev-1 do begin
        lm1=k-1
        lp1=k+1
        if k eq 0 then lm1=0
        if k eq mlev-1 then lp1=mlev-1
        mdthdz(i,k)=(mtheta(i,lp1)-mtheta(i,lm1))/(altitude(lp1)-altitude(lm1))
    endfor
    endfor
    index=where(mtemp lt 0.)
    if index(0) ne -1L then mdthdz(index)=-99.
;
; interpolate CO data to MERRA theta surfaces
;
    mco_th=fltarr(mprof,nth)
    mtemp_th=fltarr(mprof,nth)
    mdthdz_th=fltarr(mprof,nth)
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
if mtemp(iprof,kk) gt 0. and mtemp(iprof,kk+1) gt 0. then begin
if mtheta(iprof,kk) lt zlev and mtheta(iprof,kk+1) ge zlev then begin
   zscale=(mtheta(iprof,kk+1)-zlev)/(mtheta(iprof,kk+1)-mtheta(iprof,kk))
   mtemp_th(iprof,k)= mtemp(iprof,kk+1)+zscale*(mtemp(iprof,kk)-mtemp(iprof,kk+1))
   mdthdz_th(iprof,k)= mdthdz(iprof,kk+1)+zscale*(mdthdz(iprof,kk)-mdthdz(iprof,kk+1))
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
    mtemp_swath=reform(mtemp_th(mindex,*))
    mdthdz_swath=reform(mdthdz_th(mindex,*))
    mlon_swath=mlon(mindex)
    mlat_swath=mlat(mindex)
    muttime_swath=muttime(mindex)

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='yz_cross_pole_nh_merra_mls_temp_dthdz_'+sdate+'_'+slon1+'_4pan.ps'
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
    imin=180.
    imax=280.
    int=5.
    level=imin+int*findgen(nlvls)
    tyz=smooth(tyz,9,/edge_truncate)
    contour,tyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='MERRA Temp',$
            ytitle='Theta',xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[min(th),max(th)],charsize=1.5
    contour,tyz,alat,th,levels=level,/overplot,/follow,color=0
;   contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
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
    contour,dthdzyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='MLS Temp',$
            ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[min(th),max(th)],/nodata,charsize=1.5
    axis,xrange=[0.,n_elements(mlon_swath)-1L],yrange=[min(th),max(th)],/save,/data
    mtemp_swath=smooth(mtemp_swath,9,/edge_truncate)
    contour,mtemp_swath,findgen(n_elements(mlon_swath)),th,levels=level,/fill,c_color=col1,/overplot
    contour,mtemp_swath,findgen(n_elements(mlon_swath)),th,levels=level,/follow,color=0,/overplot

;   contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
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
    imin=10.
    imax=195.
    int=15.
    level=imin+int*findgen(nlvls)
    dthdzyz=smooth(dthdzyz,9,/edge_truncate)
    contour,dthdzyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='MERRA dth/dz',$
            xtitle=slon1+'       '+slon2,ytitle='Theta',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[min(th),max(th)],charsize=1.5
    contour,dthdzyz,alat,th,levels=level,/overplot,/follow,color=0
;   contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;   contour,markyz,alat,th,levels=[0.1],color=0,thick=5,/overplot
      imin=min(level)
      imax=max(level)
      xmnb=xorig(2)+xlen+cbaryoff
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,yorig(2)+cbarydel,yorig(2)+ylen-cbarydel
      !type=2^2+2^3+2^5+2^7
      plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K/km',color=0,xticks=4,charsize=1.5
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
    contour,pvyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,title='MLS dth/dz',$
            xtitle=slon1+'       '+slon2,ytickname=[' ',' ',' ',' ',' ',' ',' ',' '],$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[min(th),max(th)],/nodata,charsize=1.5
    axis,xrange=[0.,n_elements(mlon_swath)-1L],yrange=[min(th),max(th)],/save,/data
    mdthdz_swath=smooth(mdthdz_swath,9,/edge_truncate)
    contour,mdthdz_swath,findgen(n_elements(mlon_swath)),th,levels=level,/cell_fill,c_color=col1,/overplot
    contour,mdthdz_swath,findgen(n_elements(mlon_swath)),th,levels=level,/follow,color=0,/overplot
;   contour,spyz,alat,th,levels=[75.],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
;loadct,0
;    axis,xrange=[0.,n_elements(mlon_swath)-1L],yrange=[min(th),max(th)],/save,/data
;    mco_swath=smooth(mco_swath,9)
;    contour,mco_swath,findgen(n_elements(mlon_swath)),th,levels=[5.,5.5,6.],/follow,c_color=0,/overplot,thick=5
;loadct,39
      imin=min(level)
      imax=max(level)
      xmnb=xorig(3)+xlen+cbaryoff
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,yorig(3)+cbarydel,yorig(3)+ylen-cbarydel
      !type=2^2+2^3+2^5+2^7
      plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K/km',color=0,xticks=4,charsize=1.5
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
        spawn,'convert -trim yz_cross_pole_nh_merra_mls_temp_dthdz_'+sdate+'_'+slon1+'_4pan.ps -rotate -90 '+$
                            'yz_cross_pole_nh_merra_mls_temp_dthdz_'+sdate+'_'+slon1+'_4pan.jpg'
        spawn,'/usr/bin/rm yz_cross_pole_nh_merra_mls_temp_dthdz_'+sdate+'_'+slon1+'_4pan.ps'
     endif
     jumpday:
goto,jump
end
