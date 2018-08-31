;
; plot daily longitude-altitude sections of temperature and zonal wind 
;
@stddat
@kgmt
@ckday
@kdate
@rd_nogaps_nc3
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
xorig=[0.1,0.4,0.7,0.1,0.4,0.7]
yorig=[0.6,0.6,0.6,0.2,0.2,0.2]
xlen=0.25
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
dir='/aura7/harvey/NOGAPS_Alpha/Datfiles/NOGAPSA_'
gdir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
stimes=[$
'_AVG.V01.']
lstmn=12L & lstdy=10L & lstyr=2008L
ledmn=2L & leddy=28L & ledyr=2009L
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
icount=0L
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
; read NOGAPS data
;
      ifile=dir+sdate+'12_MetO_sabmls_aim9c.nc3'
      rd_nogaps_nc3,ifile,nc,nr,nthn,alon,alat,thn,pv2n,p2n,msf2n,u2n,v2n,q2n,$
         qdf2n,mark2n,sf2n,vp2n,o32n,h2o2n,iflag
      if iflag eq 1 then goto,jump
      if icount eq 0L then begin
         rlat=60.
         print,alat
         read,'Enter latitude ',rlat
         index=where(alat eq rlat)
         ilat=index(0)
         slat=strcompress(rlat,/remove_all)
         icount=1
      endif
print,'Read NOGAPS'
      t2n=0.*p2n
      for lev=0L,nthn-1L do t2n(*,*,lev)=thn(lev)*(p2n(*,*,lev)/1000.)^.286
      z2n=(msf2n-1004.*t2n)/9.86/1000.
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.
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
; longitude altitude arrays
;
      tzm=fltarr(nc,nthn)
      uzm=fltarr(nc,nthn)
      markzm=fltarr(nc,nthn)
      tzmn=fltarr(nc,nthn)
      uzmn=fltarr(nc,nthn)
      markzmn=fltarr(nc,nthn)
      for k=0L,nthn-1L do begin
          tzmn(*,k)=t2n(ilat,*,k)
          uzmn(*,k)=u2n(ilat,*,k)
          markzmn(*,k)=mark2n(ilat,*,k)
          for kk=1L,nth-1L do begin
              if thn(k) eq th(kk) then begin
                 tzm(*,k)=t2(ilat,*,kk)
                 uzm(*,k)=u2(ilat,*,kk)
                 markzm(*,k)=mark2(ilat,*,kk)
              endif
          endfor
      endfor

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='xz_u_geos5_nogaps_'+sdate+'_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
       !p.thick=2.0
       !p.charsize=2.0
    endif
erase
loadct,39
xyouts,.35,.9,sdate+'  Lat= '+slat,/normal,charsize=2,charthick=2,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=26
col1=1+indgen(nlvls)*icolmax/nlvls
tlevel=170.+5.*findgen(nlvls)
index=where(tzm eq 0.)
if index(0) ne -1L then tzm(index)=0./0.
contour,tzm,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
      ytitle='Theta (K)',title='GEOS-5 Temp',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel mod 10. eq 0)
contour,tzm,alon,thn,levels=tlevel(index),color=0,/follow,/overplot,c_labels=1+0*index
contour,markzm,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzm,alon,thn,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
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
contour,tzmn,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
        yticks=1,ytickname=[' ',' '],title='NOGAPS Temp',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel mod 10. eq 0)
contour,tzmn,alon,thn,levels=tlevel(index),color=0,/follow,/overplot,c_labels=1+0*index
contour,markzmn,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzmn,alon,thn,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
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
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
nlvls=n_elements(col2)
diff=tzm-tzmn
index=where(finite(tzm) eq 0)
if index(0) ne -1L then diff(index)=0./0.
level=-20.+4.*findgen(11)
contour,diff,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
        yticks=1,ytickname=[' ',' '],title='GEOS-5 - NOGAPS',xticks=6,/cell_fill,c_color=col2,levels=level
index=where(level gt 0.)
contour,diff,alon,thn,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,diff,alon,thn,levels=level(index),color=mcolor,/follow,/overplot,c_linestyle=5
;contour,diff,alon,thn,levels=[0],color=0,/follow,/overplot,thick=3
contour,markzm,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
loadct,0
contour,markzm,alon,thn,levels=[-0.5],color=150,/follow,/overplot,thick=5
loadct,39
imin=min(level)
imax=max(level)
ymnb=yorig(2) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(j)
x1=x1+dx
endfor

loadct,39
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
slevel=-120.+10.*findgen(nlvls)
index=where(uzm eq 0.)
if index(0) ne -1L then uzm(index)=0./0.
contour,uzm,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
        ytitle='Theta (K)',title='GEOS-5 U',xticks=6,/cell_fill,c_color=col1,levels=slevel
index=where(slevel gt 0.)
contour,uzm,alon,thn,levels=slevel(index),color=0,/follow,/overplot
index=where(slevel lt 0.)
contour,uzm,alon,thn,levels=slevel(index),color=mcolor,/follow,/overplot,c_linestyle=5
;contour,uzm,alon,thn,levels=[0],color=0,/follow,/overplot,thick=3
contour,markzm,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzm,alon,thn,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
imin=min(slevel)
imax=max(slevel)
ymnb=yorig(3) -cbaryoff
ymxb=ymnb  +cbarydel
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

xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
contour,uzmn,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
        yticks=1,ytickname=[' ',' '],title='NOGAPS U',xticks=6,/cell_fill,c_color=col1,levels=slevel
index=where(slevel gt 0.)
contour,uzmn,alon,thn,levels=slevel(index),color=0,/follow,/overplot
index=where(slevel lt 0.)
contour,uzmn,alon,thn,levels=slevel(index),color=mcolor,/follow,/overplot,c_linestyle=5
;contour,uzmn,alon,thn,levels=[0],color=0,/follow,/overplot,thick=3
contour,markzmn,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzmn,alon,thn,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
imin=min(slevel)
imax=max(slevel)
ymnb=yorig(3) -cbaryoff
ymxb=ymnb  +cbarydel
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

xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
nlvls=n_elements(col2)
diff=uzm-uzmn
index=where(finite(uzm) eq 0)
if index(0) ne -1L then diff(index)=0./0.
level=-40.+8.*findgen(11)
contour,diff,alon,thn,/noeras,xrange=[0.,360.],yrange=[min(thn),max(thn)],charsize=1.5,color=0,$
        yticks=1,ytickname=[' ',' '],title='GEOS-5 - NOGAPS',xticks=6,/cell_fill,c_color=col2,levels=level
index=where(level gt 0.)
contour,diff,alon,thn,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,diff,alon,thn,levels=level(index),color=mcolor,/follow,/overplot,c_linestyle=5
;contour,diff,alon,thn,levels=[0],color=0,/follow,/overplot,thick=3
contour,markzm,alon,thn,levels=[0.5],color=0,/follow,/overplot,thick=5
loadct,0
contour,markzm,alon,thn,levels=[-0.5],color=150,/follow,/overplot,thick=5
loadct,39
imin=min(level)
imax=max(level)
ymnb=yorig(5) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(m/s)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(j)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim xz_u_geos5_nogaps_'+sdate+'_'+slat+'.ps -rotate -90 xz_u_geos5_nogaps_'+sdate+'_'+slat+'.jpg'
       spawn,'/usr/bin/rm xz_u_geos5_nogaps_'+sdate+'_'+slat+'.ps'
    endif
goto, jump
end
