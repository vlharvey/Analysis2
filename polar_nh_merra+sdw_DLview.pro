;
; compare SDW and MERRA Arctic vortex at a given level
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3
@rd_sdwaccm4_nc3

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
cbaryoff=0.065
cbarydel=0.02
lstmn=2
lstdy=1
lstyr=2005
ledmn=5
leddy=1
ledyr=2005
lstday=0
ledday=0
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      MERRA Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirw='/Volumes/cloud/data/WACCM_data/Datfiles_SD/'
stime=['00Z','06Z','12Z','18Z']
ntime=n_elements(stime)

; Compute initial Julian date
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
      if ndays gt ledday then stop,' Normal termination condition '
      sdate=string(FORMAT='(i4.4,i2.2,i2.2)',iyr,imn,idy)

      for itime=0L,0L do begin

      ifile=string(FORMAT='(i4.4,i2.2,i2.2)',iyr,imn,idy)+'.nc3'
      rd_merra_nc3,dir+ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
      ifile=dirw+'f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h5.'+sdate+'.nc3'
      rd_sdwaccm4_nc3,ifile,ncsdw,nrsdw,nthsdw,alonsdw,alatsdw,thsdw,$
         pv2sdw,p2sdw,g2sdw,u2sdw,v2sdw,q2sdw,qdf2sdw,mark2sdw,sf2sdw,h2o2sdw,n2o2sdw,o32sdw,iflg
      if iflag ne 0L then goto,jump
tmp2=0.*p2
for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
tmp2sdw=0.*p2sdw
for k=0L,nthsdw-1L do tmp2sdw(*,*,k)=thsdw(k)*(p2sdw(*,*,k)/1000.)^0.286

sp2=sqrt(u2^2.+v2^2.)
sp2sdw=sqrt(u2sdw^2.+v2sdw^2.)
      if iflag eq 1 then goto,jump
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.

; select theta levels to plot
    if icount eq 0L then begin
       rtheta=1000.
;      read,'Enter theta level ',rtheta
       zindex=where(th eq rtheta,nth2)
       thlevs=reverse(strcompress(string(fix(th(zindex))))+' K')
       thlw=min(th(zindex))
       thup=max(th(zindex))
       th2=reverse(th(zindex))
       nr2=nr/2
       x2d=fltarr(nc+1,nr2)
       y2d=fltarr(nc+1,nr2)
       for i=0,nc do y2d(i,*)=alat(nr/2:nr-1)
       for j=0,nr2-1 do x2d(*,j)=x
       dy=alat(1)-alat(0)
       icount=1L
    endif

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename=sdate+'_merra+sdw_arctic.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
       !p.thick=2.0                   ;Plotted lines twice as thick
       !p.charsize=2.0
    endif

; coordinate transformation
    xcn=fltarr(nc+1,nr2)
    ycn=fltarr(nc+1,nr2)
    for j=nr/2,nr-1 do begin
        ANG = (90. - alat(j)) * RADG * 0.5
        FACTOR = TAN(ANG) * FAC20
        for i=0,nc do begin
            THETA = (x(i) - 90.) * RADG
            xcn(i,j-nr/2) = FACTOR * COS(THETA)
            ycn(i,j-nr/2) = FACTOR * SIN(THETA)
        endfor
    endfor
    xcs=fltarr(nc+1,nr2)
    ycs=fltarr(nc+1,nr2)
    for j=0,nr/2-1 do begin
        ANG = (90. + alat(j)) * RADG * 0.5
        FACTOR = TAN(ANG) * FAC20
        for i=0,nc do begin
            THETA = (x(i) - 90.) * RADG
            xcs(i,j) = FACTOR * COS(THETA)
            ycs(i,j) = -1.0 * FACTOR * SIN(THETA)
        endfor
    endfor

    erase
    !psym=0
    set_viewport,.1,.9,.1,.9
    !type=2^5+2^6
    dum=fltarr(nc+1,nr2)
    col1=fltarr(nth2)
    kk=0L
        index=where(th eq th2(kk))
        lev=index(0)
        mark1=transpose(mark2(*,*,lev))
        sf1=transpose(sf2(*,*,lev))
        sf=fltarr(nc+1,nr2)
        sf(0:nc-1,0:nr2-1)=sf1(0:nc-1,nr/2:nr-1)    ; NH
        sf(nc,*)=sf(0,*)

        pv1=transpose(pv2(*,*,lev))
        p1=transpose(p2(*,*,lev))
        sp1=transpose(sp2(*,*,lev))
        mpv1=pv1*((th(lev)/300.))^(-9./2.)
;
; extract SDW vortex
;
        index=where(thsdw eq th2(kk))
        if index(0) ne -1L then mark1sdw=transpose(mark2sdw(*,*,index(0)))
        if index(0) eq -1L then index=where(abs(thsdw-th2(kk)) eq min(abs(thsdw-th2(kk))))
print,th2(kk),index(0),thsdw(index(0))
mark1sdw=transpose(mark2sdw(*,*,index(0)))
sp1sdw=transpose(sp2sdw(*,*,index(0)))
        temp1sdw=transpose(tmp2sdw(*,*,index(0)))
        tempsdw=fltarr(nc+1,nr2)
        tempsdw(0:nc-1,0:nr2-1)=temp1sdw(0:nc-1,nr/2:nr-1)    ; NH
        tempsdw(nc,*)=tempsdw(0,*)
        sf1sdw=transpose(sf2sdw(*,*,index(0)))
        sfsdw=fltarr(nc+1,nr2)
        sfsdw(0:nc-1,0:nr2-1)=sf1sdw(0:nc-1,nr/2:nr-1)    ; NH
        sfsdw(nc,*)=sfsdw(0,*)

; temperature
        temp1=th(lev)*(p1/1000.)^.286
print,th(lev),min(temp1),max(temp1)
        temp=fltarr(nc+1,nr2)
        temp(0:nc-1,0:nr2-1)=temp1(0:nc-1,nr/2:nr-1)    ; NH
        temp(nc,*)=temp(0,*)
sp=fltarr(nc+1,nr2)
sp(0:nc-1,0:nr2-1)=sp1(0:nc-1,nr/2:nr-1)    ; NH
sp(nc,*)=sp(0,*)
;       index=where(y2d lt 30. or temp eq 0.)
        index=where(temp eq 0.)
        if index(0) ne -1 then temp(index)=1.e15

; pressure of theta surface
        index=where(p1 ne 0.)
        if n_elements(index) eq 1L then goto,jumplev
	result=moment(p1(index))
	avgz=result(0)
        savgz=strcompress(string(FORMAT='(F7.3)',avgz))

; draw latitude circles
        if kk eq 0 then begin
        !psym=0
        lon=findgen(361)
        lonp=0.*lon
        latp=0.*lon
        for k=0,0 do begin
            if k eq 0 then lat=0.*fltarr(361)
            if k eq 1 then lat=30.+0.*fltarr(361)
            if k eq 2 then lat=60.+0.*fltarr(361)
            for j=0,360 do begin
                ANG = (90. - lat(j)) * RADG * 0.5
                FACTOR = TAN(ANG) * FAC20
                THETA = (lon(j) - 90.) * RADG
                lonp(j) = FACTOR * COS(THETA)
                latp(j) = FACTOR * SIN(THETA)
            endfor
            oplot,lonp,latp,color=0,thick=2
        endfor
        MAP_SET,90,0,0,/stereo,/contin,/grid,/noborder,/noeras,londel=90.,$
            label=1,lonlab=1,charsize=2,latdel=180.,color=0
;
; fill continents grey
;
        loadct,0
        map_continents,mlinethick=2,color=mcolor*.4,/fill_continents
;
; superimpose stream function
;
        dum(0:nc-1,0:nr2-1)=sf1(0:nc-1,nr/2:nr-1)    ; NH
        dum(nc,*)=dum(0,*)
        smin=min(dum)
        smax=max(dum)
        sint=(smax-smin)/15.
        sflevel=smin+sint*findgen(15)
        contour,dum,xcn,ycn,levels=sflevel,color=0,c_labels=0+0.*sflevel,thick=2,title=sdate,charsize=2,charthick=2
        loadct,39
        contour,sfsdw,xcn,ycn,levels=sflevel,color=100,c_labels=0+0.*sflevel,thick=2,title=sdate,charsize=2,charthick=2
        endif
 
        nz2=(kk+1.)*(1./(nth2+1.))
        col1(kk)=nz2*icolmax
        dum=fltarr(nc+1,nr2)
        dum(0:nc-1,0:nr2-1)=mark1(0:nc-1,nr/2:nr-1)    ; NH
        dum(nc,*)=dum(0,*)
        dumsdw=fltarr(nc+1,nr2)
        dumsdw(0:nc-1,0:nr2-1)=mark1sdw(0:nc-1,nr/2:nr-1)    ; NH
        dumsdw(nc,*)=dumsdw(0,*)
        spsdw=fltarr(nc+1,nr2)
        spsdw(0:nc-1,0:nr2-1)=sp1sdw(0:nc-1,nr/2:nr-1)    ; NH
        spsdw(nc,*)=spsdw(0,*)
print,max(sp),max(spsdw),max(spsdw)-max(sp)
print,max(temp),max(tempsdw),max(tempsdw)-max(temp)
;
; sub-vortex modification
;
        if th(lev) le 300. then begin
           lindex=where(dum gt 0.0,nl)
           mpv=fltarr(nc+1,nr2)
           mpv(0:nc-1,0:nr2-1)=mpv1(0:nc-1,nr/2:nr-1)
           mpv(nc,*)=mpv(0,*)
           if lindex(0) eq -1 then begin
              index=where(mpv ge 0.0004 and y2d ge 55.)
              if index(0) ne -1 then dum(index)=1.
           endif
           if lindex(0) ne -1 then begin
              if min(y2d(lindex)) le 55. then begin
                 index=where(mpv ge 0.0004 and y2d ge 55.)
                 if index(0) ne -1 then dum(index)=1.
                 index=where(mpv lt 0.0004)
                 if index(0) ne -1 then dum(index)=0.
              endif
           endif
        endif

        lindex=where(dum gt 0.0,nl)
        imin=170.
        imax=310.
	if lindex(0) ne -1 then begin
;           for ii=0,nl-1 do begin
;               if temp(lindex(ii)) ne 1.e15 then $
;               oplot,[xcn(lindex(ii)),xcn(lindex(ii))],$
;                     [ycn(lindex(ii)),ycn(lindex(ii))],$
;                     psym=8,symsize=2,$
;                     color=((temp(lindex(ii))-imin)/(imax-imin))*icolmax
;;              if temp(lindex(ii)) eq 1.e15 then $
;;              oplot,[xcn(lindex(ii)),xcn(lindex(ii))],$
;;                    [ycn(lindex(ii)),ycn(lindex(ii))],$
;;                    psym=8,symsize=0.5,color=0
;           endfor
;          contour,temp,xcn,ycn,levels=[180.],color=10,$
;                  thick=3,max_value=1.e15
print,'MERRA y0 ',min(y2d(lindex))
           contour,dum,xcn,ycn,levels=[0.1],color=0,$
                   c_labels=0,thick=7,max_value=1.e15
loadct,0
           contour,sp,xcn,ycn,levels=[80],color=150,$
                   c_labels=0,thick=7,max_value=1.e15
           contour,sp,xcn,ycn,levels=[60],color=100,$
                   c_labels=0,thick=7,max_value=1.e15
loadct,39
        endif
;
; SDW
;
        lindex=where(dumsdw gt 0.0,nl)
        imin=170.
        imax=310.
        if lindex(0) ne -1 then begin
;           for ii=0,nl-1 do begin
;               if temp(lindex(ii)) ne 1.e15 then $
;               oplot,[xcn(lindex(ii)),xcn(lindex(ii))],$
;                     [ycn(lindex(ii)),ycn(lindex(ii))],$
;                     psym=8,symsize=2,$
;                     color=((temp(lindex(ii))-imin)/(imax-imin))*icolmax
;;              if temp(lindex(ii)) eq 1.e15 then $
;;              oplot,[xcn(lindex(ii)),xcn(lindex(ii))],$
;;                    [ycn(lindex(ii)),ycn(lindex(ii))],$
;;                    psym=8,symsize=0.5,color=0
;           endfor
;          contour,temp,xcn,ycn,levels=[180.],color=10,$
;                  thick=3,max_value=1.e15
print,'SDW y0 ',min(y2d(lindex))
           contour,dumsdw,xcn,ycn,levels=[0.1],color=150,$
                   c_labels=0,thick=5,max_value=1.e15
           contour,spsdw,xcn,ycn,levels=[80],color=250,$
                   c_labels=0,thick=5,max_value=1.e15
           contour,spsdw,xcn,ycn,levels=[60],color=210,$
                   c_labels=0,thick=5,max_value=1.e15
        endif
;
; anticyclones
;
        lindex=where(dum lt -100000.0,nl)
        if lindex(0) ne -1 then begin
;          oplot,xcn(lindex),ycn(lindex),psym=8,symsize=2,color=0
loadct,0
;          contour,dum,xcn,ycn,levels=[-0.1],color=mcolor*.3,$
;                  c_labels=0,thick=3
           nhigh=abs(min(dum(lindex)))
sdum=0.*dum
        sdum(0:nc-1,0:nr2-1)=sf1(0:nc-1,nr/2:nr-1)    ; NH
        sdum(nc,*)=sdum(0,*)
dx=x2d(1,0)-x2d(0,0)
           for ihigh=0,nhigh-1 do begin
               index=where(dum eq -1.0*(ihigh+1))
               if index(0) eq -1 then goto,jump1
               if min(y2d(index)) le 13.7500 then goto,jump1
;              if min(y2d(index)) le 2.5 then goto,jump1
               sedge=min(sdum(index))     ; value of SF to contour
               tmp=sdum
               xmax=max(x2d(index))+1.0*dx      ; isolate region
               xmin=min(x2d(index))-1.0*dx
               ymax=max(y2d(index))+2.0*dy
               ymin=min(y2d(index))-2.0*dy
               if xmin lt x(0) and xmax gt x(nc) then begin     ; GM
                  index=where(x2d gt 180. and dum eq -1.0*(ihigh+1))
                  xmax2=min(x2d(index))-2.0*dx
                  index=where(x2d lt 180. and dum eq -1.0*(ihigh+1))
                  xmin2=max(x2d(index))+2.0*dx
                  index=where((x2d lt xmax2 and x2d gt xmin2) or (y2d lt ymin or y2d gt ymax))
               endif
               if xmin gt x(0) or xmax lt x(nc) then $
                  index=where(x2d lt xmin or x2d gt xmax or y2d lt ymin or y2d gt ymax)
               if index(0) ne -1 then tmp(index)=-9999.
;              index=where(tmp ne -9999. and y2d gt 13.7500 and dum eq -1.0*(ihigh+1))
               index=where(tmp ne -9999. and y2d gt 0. and dum eq -1.0*(ihigh+1))
               oplot,xcn(index),ycn(index),psym=8,color=0,symsize=0.5
               contour,tmp,xcn,ycn,levels=[sedge],color=icolmax*.7,$
                 c_linestyle=0,/overplot,min_value=-9999.,thick=3
               jump1:
           endfor               ; loop over anticyclones

loadct,39
        endif
jumplev:
        xyouts,.83,.8,savgz,color=0,/normal,charsize=2,charthick=2
        xyouts,.08,.8,thlevs(kk),color=0,/normal,charsize=2,charthick=2

    !psym=0
;    set_viewport,.2,.78,.14-cbaryoff,.14-cbaryoff+cbarydel
;    !type=2^2+2^3+2^6
;    iint=(imax-imin)/12.
;    level=imin+iint*findgen(13)
;    plot,[imin,imax],[0,0],yrange=[0,10],$
;          xrange=[imin,imax],xtitle='MERRA Temperature',/noeras,$
;          xtickname=strcompress(string(fix(level)),/remove_all),$
;          xstyle=1,xticks=12,charsize=1.5,color=0,charthick=2
;    ybox=[0,10,10,0,0]
;    x1=imin
;    dx=(imax-imin)/float(nth2)
;    for j=0,nth2-1 do begin
;      xbox=[x1,x1,x1+dx,x1+dx,x1]
;      polyfill,xbox,ybox,color=col1(j)
;      x1=x1+dx
;    endfor
;    !p.charthick=1.
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim '+sdate+'_merra+sdw_arctic.ps -rotate -90 '+sdate+'_merra+sdw_arctic.jpg'
       spawn,'rm -f '+sdate+'_merra+sdw_arctic.ps'
    endif

    endfor

goto, jump

end
