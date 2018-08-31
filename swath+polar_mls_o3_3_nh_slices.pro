;
; plot 3 longitudinal slices in the NH each day
; plot MLS O3 in polar and vertical cross sections of high resolution swaths
;
@stddat
@kgmt
@ckday
@kdate
@aura2date
@rd_ukmo_nc3

mno=[31,28,31,30,31,30,31,31,30,31,30,31]
dlat=2.0
nlat=long((180./dlat) +1.0)
latbin=-90.+dlat*findgen(nlat)
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
icmm1=icolmax-1
icmm2=icolmax-2
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.1,0.1,0.1]
yorig=[0.65,0.35,0.05]
ylen=0.3
xlen=0.3
cbaryoff=0.015
cbarydel=0.01
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=0
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
dirm='/aura6/data/MLS_data/Datfiles/'

lstmn=10 & lstdy=21 & lstyr=2005 & lstday=0
ledmn=5 & leddy=31 & ledyr=2006 & ledday=0
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '
;
; read UKMO data
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
print,syr+smn+sdy
;
; restore MLS data
;
      dum=findfile(dirm+'MLS-Aura_L2GP_'+syr+smn+sdy+'.sav')
      if dum(0) eq '' then goto,jump
      restore,dirm+'MLS-Aura_L2GP_'+syr+smn+sdy+'.sav'
print,'opened MLS-Aura_L2GP_'+syr+smn+sdy+'.sav'
      mpress=o3.pressure
      mlev=n_elements(mpress)
      mprof=o3.ntimes
      mtime=o3.TIME
      mlat=o3.latitude
      mlon=o3.longitude
      index=where(mlon lt 0.)
      mlon(index)=mlon(index)+360.
      mtemp=tp.L2GPVALUE
      mo3=o3.L2GPVALUE*1.e6
      mo3prec=o3.L2GPPRECISION
      mo3stat=o3.STATUS
      mo3qual=o3.QUALITY
      mo3mask=0.*mo3
;
; remove bad times
;
      index=where(mtime gt 0.,mprof)
      mtime=mtime(index)
      mlat=mlat(index)
      mlon=mlon(index)
      mtemp=mtemp(*,index)
      mo3=mo3(*,index)
      mo3prec=mo3prec(*,index)
      mo3stat=mo3stat(*,index)
      mo3qual=mo3qual(*,index)
      mo3mask=mo3mask(*,index)
;
; use quality, status, and precision flags to remove suspect data
;
      o3bad=where(mo3prec lt 0.)
      if o3bad(0) ne -1L then mo3mask(o3bad)=-99.
      o3bad=where(mo3stat mod 2 ne 0L)		; o3status=0 is good, all odd values are bad
      if o3bad(0) ne -1L then mo3mask(*,o3bad)=-99.
      o3bad=where(mo3qual lt 0.1)			; do not use if quality < 0.1
      if o3bad(0) ne -1L then mo3mask(*,o3bad)=-99.
;
; convert elapsed seconds to dates (yyyymmddhh)
;
      aura2date,mdate,mtime
;
; time is elapsed seconds since midnight 1 Jan 1993 to hours today
; convert to daily UT time (0-24 hours)
;
      muttime=mtime
      istime=1993010100L
      ehr=muttime/60./60.       ; convert time from seconds to hours
      hh2=0.d*muttime
      for n=0L,mprof-1L do begin
          yy1=istime/1000000
          if yy1 mod 4 eq 0 then mno(1)=29L
          if yy1 mod 4 ne 0 then mno(1)=28L
          mm1=istime/10000L-yy1*100L
          dd1=istime/100L-yy1*10000L-mm1*100L
          dd2=dd1+long(ehr(n))/24L
          hh1=istime-yy1*1000000L-mm1*10000L-dd1*100L
          yy2=yy1 & mm2=mm1
          while dd2 gt mno(mm2-1) do begin
                dd2=dd2-mno(mm2-1)
                mm2=mm2+1L
                if mm2 gt 12L then begin
                   mm2=mm2-12L
                   yy2=yy2+1L
                   if yy2 mod 4 eq 0 then mno(1)=29
                   if yy2 mod 4 ne 0 then mno(1)=28
                endif
          endwhile
          hh2(n)=ehr(n) mod 24
          if hh2(n) ge 24. then begin
             hh2(n)=hh2(n)-24.
             dd2=dd2+1L
             if dd2 gt mno(mm2-1L) then begin
                dd2=dd2-mno(mm2-1L)
                mm2=mm2+1L
                if mm2 gt 12L then begin
                   mm2=mm2-12L
                   yy2=yy2+1L
                endif
             endif
          endif
      endfor
      muttime=hh2
;
; calculate potential temperature
;
      mtheta=0.*mtemp
      for i=0L,mlev-1L do mtheta(i,*)=mtemp(i,*)*(1000./mpress(i))^0.286
      index=where(mtemp lt 0.)
      if index(0) ne -1L then mtheta(index)=-999.
      mpress2=0.*mo3
      mtime2=0.*mo3
      mlat2=0.*mo3
      mlon2=0.*mo3
      for i=0L,mprof-1L do mpress2(*,i)=mpress
      for i=0L,mlev-1L do begin
          mtime2(i,*)=muttime
          mlat2(i,*)=mlat
          mlon2(i,*)=mlon
      endfor
;
; read UKMO data
;
      sdate=strcompress(string(mdate(mprof-1)),/remove_all)
      sdate=strmid(sdate,0,8)
      syr=strmid(sdate,0,4)
      uyr=strmid(syr,2,2)
      smn=strmid(sdate,4,2)
      imn=fix(smn)
      sdy=strmid(sdate,6,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      print,ifile
      rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
      xz2d=fltarr(nc,nth)
      yz2d=fltarr(nc,nth)
      for k=0,nth-1 do xz2d(*,k)=alon
      for j=0,nc-1 do yz2d(j,*)=th
      x=fltarr(nc+1)
      x(0:nc-1)=alon
      x(nc)=alon(0)+360.
      x2d=fltarr(nc+1,nr)
      y2d=fltarr(nc+1,nr)
      for i=0,nc do y2d(i,*)=alat
      for j=0,nr-1 do x2d(*,j)=x
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='swath+polar_mls_o3_3_nh_slices_'+sdate+'.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot MLS
;
    erase
    itheta=1000.
    stheta=strcompress(string(fix(itheta)),/remove_all)
    index=where(mtheta gt itheta-80. and mtheta le itheta+80. and mlat2 gt 0.,npt)
    th0=string(FORMAT='(I4)',itheta)
    xdata=mlon2(index) & ydata=mlat2(index)
    o3data=mo3(index)
    o3mask=mo3mask(index)
    tdata=mtemp(index)
    index=where(xdata lt 0.)
    if index(0) ne -1 then xdata(index)=xdata(index)+360.
    omin=-2.0
    omax=max(o3data)
    nlvls=30
    level=omin+((omax-omin)/nlvls)*findgen(nlvls)
    col1=1+indgen(nlvls)*mcolor/nlvls
    xyouts,.4,.97,sdate+' '+stheta+' K',charsize=2,/normal,color=0
;
; MetO at itheta
;
    index=where(itheta eq th)
    ktheta=index(0)
    pv1=transpose(pv2(*,*,ktheta))
    p1=transpose(p2(*,*,ktheta))
    sf1=transpose(sf2(*,*,ktheta))
    mark1=transpose(mark2(*,*,ktheta))
    pv=fltarr(nc+1,nr)
    pv(0:nc-1,0:nr-1)=pv1
    pv(nc,*)=pv(0,*)
    p=fltarr(nc+1,nr)
    p(0:nc-1,0:nr-1)=p1
    p(nc,*)=p(0,*)
    sf=fltarr(nc+1,nr)
    sf(0:nc-1,0:nr-1)=sf1
    sf(nc,*)=sf(0,*)
    mark=fltarr(nc+1,nr)
    mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
    mark(nc,*)=mark(0,*)

for ipan=0L,n_elements(xorig)-1L do begin
    !type=2^2+2^3
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    map_set,90,0,-180,/stereo,/contin,/grid,/noeras,charsize=1.5,/noborder,color=0
    contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
    index=where(mark gt 0. and y2d gt 0.)
    if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
    for i=0L,n_elements(o3data)-1L do begin
        if o3mask(i) ne -99. then $
           oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                 color=((o3data(i)-omin)/(omax-omin))*mcolor,symsize=1.5
    endfor
    oplot,findgen(361),0.1+0.*findgen(361),psym=8,symsize=0.2,color=0
    loadct,0
    contour,mark,x,alat,levels=[.1],c_color=mcolor*.7,/overplot,/follow,c_labels=0,/noeras,thick=5
    loadct,38
    contour,mark,x,alat,levels=[-.1],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=5
;
; superimpose MLS ozone
;
    index=where(mlat gt 80.,npt)
    tswath=muttime(index)
    flag=0.*tswath
    nswath=50L
    tsave=fltarr(nswath)
    icount=0L
    for i=0L,npt-1L do begin
        index=where(abs(tswath(i)-tswath) lt 1. and flag eq 0.)
        if index(0) ne -1L then begin
           flag(index,0)=1.0
           kindex=where(abs(muttime-tswath(index(0))) le 0.5 and mlat gt 0.)
;          oplot,mlon(kindex),mlat(kindex),psym=8,symsize=1.25,color=(float(i+1)/float(npt))*mcolor
           stime=string(FORMAT='(F4.1)',tswath(index(0)))
           tsave(icount)=tswath(index(0))
           icount=icount+1L
           xtmp=mlon(kindex)
           ytmp=mlat(kindex)
           index=where(ytmp eq min(ytmp))
;          xyouts,xtmp(index(0)),-5.,stime,/data,charsize=1.5,alignment=0.5,color=0
        endif
    endfor
    tsave=tsave(0:icount-1L)
    tplot=19.9791
;print,tsave
;read,'Central Swath times ',tplot
    if ipan eq 0L then tplot=tsave(0L)
    if ipan eq 1L then tplot=tsave(icount-1L)
    if ipan eq 2L then tplot=tsave(icount/3L)
    kindex=where(abs(muttime-tplot) le 0.5 and mlat gt 0.,mprof)
    for kk=0,mprof-1L,4L do $
        oplot,[mlon(kindex(kk)),mlon(kindex(kk))],[mlat(kindex(kk)),mlat(kindex(kk))],$
               psym=4,symsize=1.75,color=lc
    stime=string(FORMAT='(F4.1)',tplot)
;
; extract MLS swath
;
    o3swath=transpose(smooth(mo3(*,kindex),3,/edge_truncate))
    o3mask=transpose(mo3mask(*,kindex))
    index=where(o3swath lt 0.)
    if index(0) ne -1L then o3swath(index)=0.1
    o3pr=0.*o3swath
    index=where(abs(o3swath) gt 100. or o3mask eq -99.)
    if index(0) ne -1L then o3swath(index)=0.
    result=size(o3swath)
    nz=result(2)
    for k=0L,nz-1L do begin
        result=moment(o3swath(*,k))
        o3pr(*,k)=o3swath(*,k)-result(0)
    endfor
    thswath=transpose(mtheta(*,kindex))
    lonswath=transpose(mlon2(*,kindex))
    latswath=transpose(mlat2(*,kindex))
    markswath=0.*latswath
    xswath=0.*thswath
    ylabels=string(format='(f4.1)',mlat(kindex))
    xlabels=string(format='(f5.1)',mlon(kindex))
    for i=0L,mlev-1L do xswath(*,i)=findgen(mprof)
;
; interpolate MetO marker to MLS swath
;
    for ii=0L,mprof-1L do begin
    for kk=0L,mlev-1L do begin
        slon=lonswath(ii,kk)
        slat=latswath(ii,kk)
        slev=thswath(ii,kk)

        if slon lt alon(0) then slon=slon+360.
        for i=0L,nc-1L do begin
            ip1=i+1
            if i eq nc-1L then ip1=0L
            xlon=alon(i)
            xlonp1=alon(ip1)
            if i eq nc-1L then xlonp1=360.+alon(ip1)
            if slon ge xlon and slon le xlonp1 then begin
               xscale=(slon-xlon)/(xlonp1-xlon)
               goto,jumpx
            endif
        endfor
jumpx:
        for j=0L,nr-2L do begin
            jp1=j+1
            xlat=alat(j)
            xlatp1=alat(jp1)
            if slat ge xlat and slat le xlatp1 then begin
                yscale=(slat-xlat)/(xlatp1-xlat)
                goto,jumpy
            endif
        endfor
jumpy:
        for k=1L,nth-1L do begin
            kp1=k-1             ; UKMO data is "top down"
            uz=th(k)
            uzp1=th(kp1)
            if slev ge uz and slev le uzp1 then begin
               zscale=(slev-uz)/(uzp1-uz)
               pj1=mark2(j,i,k)+xscale*(mark2(j,ip1,k)-mark2(j,i,k))
               pjp1=mark2(jp1,i,k)+xscale*(mark2(jp1,ip1,k)-mark2(jp1,i,k))
               pj2=mark2(j,i,kp1)+xscale*(mark2(j,ip1,kp1)-mark2(j,i,kp1))
               pjp2=mark2(jp1,i,kp1)+xscale*(mark2(jp1,ip1,kp1)-mark2(jp1,i,kp1))
               p1=pj1+yscale*(pjp1-pj1)
               p2=pj2+yscale*(pjp2-pj2)
               markswath(ii,kk)=p1+zscale*(p2-p1)
               goto,jumpz
            endif
        endfor
jumpz:
    endfor
    endfor
;
    loadct,38    
    !type=2^2+2^3
    omin=min(o3swath)
    omax=max(o3swath)
    omin=0.
    omax=10.
    nlvls=21
    level=omin+((omax-omin)/nlvls)*findgen(nlvls)
    col1=1+indgen(nlvls)*mcolor/nlvls

    xmn=xorig(ipan)+xlen+0.12
    xmx=xorig(ipan)+xlen+0.12+xlen+0.12
    ymn=yorig(ipan)+0.05
    ymx=yorig(ipan)+ylen-0.025
    set_viewport,xmn,xmx,ymn,ymx
    contour,o3swath,xswath,thswath,levels=level,/cell_fill,title='MLS at '+stime+' UT',$
            c_color=col1,min_value=-999.,xticks=1,xrange=[0.,mprof-1L],xtickname=[' ',' '],$
            yrange=[500.,1600.],ytitle='Theta',charsize=1.25,color=0
    index=where(level gt 0.)
    contour,o3swath,xswath,thswath,levels=level(index),/follow,color=mcolor,min_value=-999.,/overplot,c_labels=0*level
    contour,o3swath,xswath,thswath,levels=[0.0],/follow,color=0,thick=3,min_value=-999.,/overplot,c_labels=0*level
    contour,markswath,xswath,thswath,levels=[-0.1],/follow,color=0,thick=5,min_value=-999.,/overplot,c_labels=0*level
    loadct,0
    contour,markswath,xswath,thswath,levels=[0.1],/follow,color=mcolor*.7,thick=5,min_value=-999.,/overplot,c_labels=0*level
    loadct,38
    for i=0L,mprof-1L,20 do begin
        plots,i,500
        plots,i,450,/data,/continue,color=0
        xyouts,i,400.,xlabels(i),/data,charsize=1.25,alignment=0.5,color=0
        xyouts,i,300.,ylabels(i),/data,charsize=1.25,alignment=0.5,color=0
    endfor
    i=mprof-1L
    plots,i,500
    plots,i,450,/data,/continue,color=0
    xyouts,i,400.,xlabels(i),/data,charsize=1.25,alignment=0.5,color=0
    xyouts,i,300.,ylabels(i),/data,charsize=1.25,alignment=0.5,color=0
    imin=omin
    imax=omax
    xmnb=xmx+.03
    xmxb=xmnb+.01
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='O3 ppmv',color=0
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

endfor
;
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim swath+polar_mls_o3_3_nh_slices_'+sdate+'.ps -rotate -90 '+$
                     'swath+polar_mls_o3_3_nh_slices_'+sdate+'.jpg'
    endif
goto,jump
END
