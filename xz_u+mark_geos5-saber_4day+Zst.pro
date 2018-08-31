;
; Superimpose height of the stratopause from Jeff
;
; GEOS-5 and SABER longitude-altitude sections on 4 different days
; of zonal winds and differences
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

ipan=0
npp=1
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
a=findgen(4)*(2*!pi/4.)
usersym,cos(a),sin(a),/fill
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.1825,0.4,0.6175,0.1825,0.4,0.6175,0.1825,0.4,0.6175,0.1825,0.4,0.6175]+0.02
yorig=[0.73,0.73,0.73,0.52,0.52,0.52,0.31,0.31,0.31,0.1,0.1,0.1]+0.02
xlen=.2
ylen=.2
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/'
sdir='/aura6/data/SABER_data/Datfiles_winds/'
hdir='/aura6/data/HIRDLS_data/Datfiles_SOSST/'
stimes=[$
'_AVG.V01.']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1
lstdy=24
lstyr=2008
ledmn=2
leddy=23
ledyr=2008
lstday=0
ledday=0
nz=121
altitude=findgen(nz)
nlat=35L
latbin=-85+5.*findgen(nlat)
dy=latbin(1)-latbin(0)
nlon=12L
lonbin=15.+30.*findgen(nlon)
dx=lonbin(1)-lonbin(0)
rlat=50.
;print,latbin
;read,' Enter latitude ',rlat
index=where(rlat eq latbin)
ilat=index(0)
slat=strcompress(long(rlat),/remove_all)
;
; Ask interactive questions- get starting/ending date and p surface
;
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
kday=ledday-lstday+1L
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
restore,'stratopause_height_for_SABER.sav
restore,'stratopause_height_for_MLS.sav
restore,'stratopause_height_for_GEOS5.sav

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then begin
         if setplot ne 'ps' then stop
         if setplot eq 'ps' then begin
            device, /close
            spawn,'convert -trim xz_u+mark_geos5-saber_'+slat+'_4days.ps '+$
                  ' -rotate -90  xz_u+mark_geos5-saber_'+slat+'_4days.jpg'
            stop
         endif
      endif
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-5 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      if sdate ne '20080124' and sdate ne '20080206' and sdate ne '20080215' and sdate ne '20080223' then goto,jump
      print,sdate
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'+sdate+stimes(0)+'nc3'
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+ifile,nc,nr,nth,glon,glat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      x=fltarr(nc+1)
      x(0:nc-1)=glon(0:nc-1)
      x(nc)=glon(0)+360.
      t2=0.*pv2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
;
; restore SABER GRID_PHI_WINDS.20080305.sav file
;
; ALAT            FLOAT     = Array[35]
; ALON            FLOAT     = Array[12]
; PRESS           FLOAT     = Array[120]
; T3D             FLOAT     = Array[12, 35, 120]
; U3D             FLOAT     = Array[12, 35, 120]
; V3D             FLOAT     = Array[12, 35, 120]
; Z3D             FLOAT     = Array[12, 35, 120]
;
      dum=findfile(sdir+'GRID_PHI_WINDS.'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,sdir+'GRID_PHI_WINDS.'+sdate+'.sav'
      tdata=t3d   ; avoid T3D intrinsic function
      udata=u3d
;     if max(tdata) eq 0. then goto,jump
      zdata=z3d/1000.
      ncs=n_elements(alon)
      nrs=n_elements(alat)
      nlv=n_elements(press)
; 
; interpolate SABER to altitude grid at ilat
;
      tdataz=fltarr(ncs,nz)
      udataz=fltarr(ncs,nz)
      for ii=0L,ncs-1 do begin
          zprof=reform(zdata(ii,ilat,*))
          for k=20L,60L do begin
              zz=altitude(k)
              for kk=1L,nlv-1L do begin
                  zup=zprof(kk) & zlw=zprof(kk-1)       ; profiles are bottom up
                  if zup ge zz and zlw le zz then begin
                     zscale=(zup-zz)/(zup-zlw)
                     tdataz(ii,k)=tdata(ii,ilat,kk-1)+zscale*(tdata(ii,ilat,kk)-tdata(ii,ilat,kk-1))
                     udataz(ii,k)=udata(ii,ilat,kk-1)+zscale*(udata(ii,ilat,kk)-udata(ii,ilat,kk-1))
                  endif
              endfor
          endfor
      endfor
print,'interpolated SABER to altitude'
;
; interpolate GEOS temperature and zonal wind to altitude surfaces
;
      t2z=fltarr(nr,nc,nz)
      u2z=fltarr(nr,nc,nz)
      mark2z=fltarr(nr,nc,nz)
      for kk=20L,60 do begin
      zz=altitude(kk)
      if max(z2) lt zz then goto,jumplev
      for j=0L,nr-1L do begin
      for i=0L,nc-1L do begin
          for k=1L,nth-1L do begin
              zup=z2(j,i,k-1) & zlw=z2(j,i,k)
              if zup ne 0. and zlw ne 0. then begin
              if zup ge zz and zlw le zz then begin
                 zscale=(zup-zz)/(zup-zlw)
                 t2z(j,i,kk)=t2(j,i,k-1)+zscale*(t2(j,i,k)-t2(j,i,k-1))
                 u2z(j,i,kk)=u2(j,i,k-1)+zscale*(u2(j,i,k)-u2(j,i,k-1))
                 mark2z(j,i,kk)=mark2(j,i,k-1)+zscale*(mark2(j,i,k)-mark2(j,i,k-1))
;print,zlw,zz,zup,zscale
;print,t2(j,i,k),t2z(j,i,kk),t2(j,i,k-1)
;stop
              endif
              endif
          endfor
      endfor
      endfor
      jumplev:
      endfor
print,'interpolated GEOS to altitude'
;
; bin GEOS temperature and U by lonbin at ilat
;
      t2zbin=fltarr(nlon,nz)
      u2zbin=fltarr(nlon,nz)
      mark2zbin=fltarr(nlon,nz)
      nbin=lonarr(nlon,nz)
      for k=20L,60 do begin	;nlv-1L do begin
          for i=0L,nlon-1L do begin
              for jj=0L,nr-1L do begin
                  if glat(jj) ge latbin(ilat)-dy/2. and glat(jj) lt latbin(ilat)+dy/2. then begin
                  for ii=0L,nc-1L do begin
                      if glon(ii) ge lonbin(i)-dx/2. and glon(ii) lt lonbin(i)+dx/2. then begin
                         t2zbin(i,k)=t2zbin(i,k)+t2z(jj,ii,k)
                         u2zbin(i,k)=u2zbin(i,k)+u2z(jj,ii,k)
                         mark2zbin(i,k)=mark2zbin(i,k)+mark2z(jj,ii,k)
                         nbin(i,k)=nbin(i,k)+1L
                      endif
                  endfor
                  endif
              endfor
          endfor
      endfor
print,'binned GEOS to SABER lon grid'
      index=where(nbin gt 0L)
      if index(0) ne -1L then t2zbin(index)=t2zbin(index)/float(nbin(index))
      if index(0) ne -1L then u2zbin(index)=u2zbin(index)/float(nbin(index))
      if index(0) ne -1L then mark2zbin(index)=mark2zbin(index)/float(nbin(index))

      if icount eq 0 and setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
;        !p.font=0
         device,font_size=9
         device,/landscape,bits=8,$
                 filename='xz_u+mark_geos5-saber_'+slat+'_4days.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         icount=1
      endif
;
; interpolate stratopause height to rlat
;
yindex=where(geos5lat gt rlat)
ilatg=yindex(0)
gheight=(reform(GEOS5STRATHEIGHT(iday,*,ilatg-1))+reform(GEOS5STRATHEIGHT(iday,*,ilatg)))/2.
gheight=smooth(gheight,3)
stop
for j=1L,n_elements(saberlat)-1L do begin
    if saberlat(j-1) lt rlat and saberlat(j) ge rlat then begin
       yscale=(saberlat(j)-rlat)/(saberlat(j)-saberlat(j-1))
       mh0=reform(SABERSTRATHEIGHT(iday,*,j-1))
       mh1=reform(SABERSTRATHEIGHT(iday,*,j))
       sheight=mh1-yscale*(mh1-mh0)
    endif
endfor
;sheight=smooth(sheight,3)
for j=1L,n_elements(mlslat)-1L do begin
    if mlslat(j-1) lt rlat and mlslat(j) ge rlat then begin
       yscale=(mlslat(j)-rlat)/(mlslat(j)-mlslat(j-1))
       mh0=reform(MLSSTRATHEIGHT(iday,*,j-1))
       mh1=reform(MLSSTRATHEIGHT(iday,*,j))
       mheight=mh1-yscale*(mh1-mh0)
    endif
endfor
;mheight=smooth(mheight,3)
sheight=mheight
;
; day 1
;
if imn eq 1 and idy eq 24 then begin
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
xyouts,xmn-0.07,ymn+0.05,sdate,charsize=1.5,charthick=2,/normal,orientation=90
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
contour,t2zbin,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],$
        ytitle='Altitude (km)',xtickname=[' ',' ',' ',' ',' ',' '],title='GEOS-5'
contour,t2zbin,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
loadct,39

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
contour,tdataz,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' '],title='SABER'
contour,tdataz,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39

; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tdataz
index=where(tdataz gt 0. and t2zbin gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=t2zbin(index)-tdataz(index)
level=-20.+4.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,lonbin,altitude,levels=level,/cell_fill,c_color=col2,min_value=-99.,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' '],title='GEOS5 - SABER'
contour,pdiff,lonbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,lonbin,altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
endif
;
; day 2
;
if imn eq 2 and idy eq 6 then begin
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
xyouts,xorig(0)-0.07,ymn+0.05,sdate,charsize=1.5,charthick=2,/normal,orientation=90
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
contour,t2zbin,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],$
        ytitle='Altitude (km)',xtickname=[' ',' ',' ',' ',' ',' ']
contour,t2zbin,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
loadct,39

xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
contour,tdataz,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' ']
contour,tdataz,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
;
; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tdataz
index=where(tdataz gt 0. and t2zbin gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=t2zbin(index)-tdataz(index)
level=-20.+4.*findgen(11)
!type=2^2+2^3
xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,lonbin,altitude,levels=level,/cell_fill,c_color=col2,min_value=-99.,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' ']
contour,pdiff,lonbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,lonbin,altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
endif
;
; day 3
;
if imn eq 2 and idy eq 15 then begin
xmn=xorig(6)
xmx=xorig(6)+xlen
ymn=yorig(6)
ymx=yorig(6)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
xyouts,xorig(0)-0.07,ymn+0.05,sdate,charsize=1.5,charthick=2,/normal,orientation=90
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
contour,t2zbin,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],ytitle='Altitude (km)',$
        xtickname=[' ',' ',' ',' ',' ',' ']
contour,t2zbin,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
loadct,39

xmn=xorig(7)
xmx=xorig(7)+xlen
ymn=yorig(7)
ymx=yorig(7)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
contour,tdataz,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' ']
contour,tdataz,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
;
; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tdataz
index=where(tdataz gt 0. and t2zbin gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=t2zbin(index)-tdataz(index)
level=-20.+4.*findgen(11)
!type=2^2+2^3
xmn=xorig(8)
xmx=xorig(8)+xlen
ymn=yorig(8)
ymx=yorig(8)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,lonbin,altitude,levels=level,/cell_fill,c_color=col2,min_value=-99.,yrange=[20,60],$
        xtickname=[' ',' ',' ',' ',' ',' '],ytickname=[' ',' ',' ',' ',' ',' ']
contour,pdiff,lonbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,lonbin,altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
endif
;
; day 4
;
if imn eq 2 and idy eq 23 then begin
xmn=xorig(9)
xmx=xorig(9)+xlen
ymn=yorig(9)
ymx=yorig(9)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
xyouts,xorig(0)-0.07,ymn+0.05,sdate,charsize=1.5,charthick=2,/normal,orientation=90
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
contour,t2zbin,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],ytitle='Altitude (km)',$
        xtitle='Longitude'
contour,t2zbin,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
loadct,39
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle=slat+' N Temperature (K)'
ybox=[0,10,10,0,0]
x1=imin
dxx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dxx,x1+dxx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dxx
endfor
xmn=xorig(10)
xmx=xorig(10)+xlen
ymn=yorig(10)
ymx=yorig(10)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
contour,tdataz,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,60],xtitle='Longitude',$
        ytickname=[' ',' ',' ',' ',' ',' ']
contour,tdataz,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
loadct,39
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle=slat+' N Temperature (K)'
ybox=[0,10,10,0,0]
x1=imin
dxx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dxx,x1+dxx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dxx
endfor
;
; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tdataz
index=where(tdataz gt 0. and t2zbin gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=t2zbin(index)-tdataz(index)
level=-20.+4.*findgen(11)
!type=2^2+2^3
xmn=xorig(11)
xmx=xorig(11)+xlen
ymn=yorig(11)
ymx=yorig(11)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,lonbin,altitude,levels=level,/cell_fill,c_color=col2,min_value=-99.,yrange=[20,60],$
        xtitle='Longitude',ytickname=[' ',' ',' ',' ',' ',' ']
contour,pdiff,lonbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,lonbin,altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
contour,mark2zbin,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
contour,mark2zbin,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=mcolor,c_labels=0*level,thick=5
loadct,0
oplot,geos5lon,gheight,psym=8,color=100
index=where(finite(sheight) eq 1)
oplot,saberlon(index),sheight(index),color=150,psym=8,symsize=1.5
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle=slat+' N Differences (K)',xticks=n_elements(level)/2
ybox=[0,10,10,0,0]
x2=imin
dxx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dxx,x2+dxx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dxx
endfor
loadct,39
endif

goto,jump
end
