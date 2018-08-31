;
; plot daily zonal mean PV and dPV/dy in GEOS-5
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.15,0.55,0.15,0.55]
yorig=[0.6,0.6,0.15,0.15]
xlen=0.3
ylen=0.25
cbaryoff=0.065
cbarydel=0.02
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
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
gdir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
stimes=[$
'_AVG.V01.']
lstmn=1L & lstdy=1L & lstyr=2007L
ledmn=9L & leddy=30L & ledyr=2007L
lstday=0L & ledday=0L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
if lstyr lt 2006 then stop,'Year out of range '
if ledyr lt 2006 then stop,'Year out of range '
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
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,gdir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
print,'Read GEOS'
;
; compute temperature
;
      t2=0.*u2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; zonal means
;
      pvzm=fltarr(nr,nth)
      dpvdyzm=fltarr(nr,nth)
      tzm=fltarr(nr,nth)
      uzm=fltarr(nr,nth)
      sfzm=fltarr(nr,nth)
      markzm=fltarr(nr,nth)
      for j=0L,nr-1L do begin
          for k=0L,nth-1L do begin
              pvzm(j,k)=total(pv2(j,*,k))/float(nc)
              tzm(j,k)=total(t2(j,*,k))/float(nc)
              sfzm(j,k)=total(sf2(j,*,k))/float(nc)
              uzm(j,k)=total(u2(j,*,k))/float(nc)
              markzm(j,k)=total(mark2(j,*,k))/float(nc)
          endfor
      endfor
;
; dPV/dy
;
      for k=0L,nth-1 do begin
          for j=0L,nr-1L do begin
              jp1=j+1
              jm1=j-1
              if j eq 0L then jm1=0L
              if j eq nr-1L then jp1=nr-1L
              dpvdyzm(j,k)=pvzm(jp1,k)-pvzm(jm1,k)
          endfor
      endfor
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_geos5_pvbar_dpvdy_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
       !p.thick=2.0
       !p.charsize=2.0
    endif
erase
loadct,39
xyouts,.375,.95,'GEOS-5 '+sdate,/normal,charsize=2,charthick=2,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=[-1.,-0.5,-0.1,-0.05,-0.01,-0.005,-0.001,-0.0005,-0.0001,-0.00005,-0.00001,0.,$
        0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.]
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,pvzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=1.5,color=0,$
      ytitle='Theta (K)',title='Zonal Mean PV',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel gt 0.)
contour,pvzm,alat,th,levels=tlevel(index),color=mcolor,/follow,/overplot
index=where(tlevel lt 0.)
contour,pvzm,alat,th,levels=tlevel(index),color=0,/follow,/overplot
contour,markzm,alat,th,levels=[0.5],color=0,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Km!u2!nkg!u-1!ns!u-1!n'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
contour,dpvdyzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=1.5,color=0,$
        ytitle='Theta (K)',title='Zonal Mean dPV/dy',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel gt 0.)
contour,dpvdyzm,alat,th,levels=tlevel(index),color=mcolor,/follow,/overplot
index=where(tlevel lt 0.)
contour,dpvdyzm,alat,th,levels=tlevel(index),color=0,/follow,/overplot
contour,markzm,alat,th,levels=[0.5],color=0,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Km!nkg!u-1!ns!u-1!n'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=-120.+10.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,uzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=1.5,color=0,$
      ytitle='Theta (K)',title='Zonal Mean Zonal Wind',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel gt 0.)
contour,uzm,alat,th,levels=tlevel(index),color=mcolor,/follow,/overplot
index=where(tlevel lt 0.)
contour,uzm,alat,th,levels=tlevel(index),color=0,/follow,/overplot
contour,markzm,alat,th,levels=[0.5],color=0,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(2) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='ms!u-1!n'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=170.+5.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,tzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=1.5,color=0,$
      ytitle='Theta (K)',title='Zonal Mean Temperature',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel gt 0.)
contour,tzm,alat,th,levels=tlevel(index),color=mcolor,/follow,/overplot
contour,markzm,alat,th,levels=[0.5],color=0,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(3) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='K'
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
       device,/close
       spawn,'convert -trim yz_geos5_pvbar_dpvdy_'+sdate+'.ps -rotate -90 yz_geos5_pvbar_dpvdy_'+sdate+'.jpg'
    endif
goto, jump
end
