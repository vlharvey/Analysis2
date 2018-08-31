;
; plot daily longitude-altitude sections of eddy height from MERRA
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.15,0.55]
yorig=[0.35,0.35]
xlen=0.3
ylen=0.3
cbaryoff=0.08
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
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
lstmn=1L & lstdy=1L & lstyr=1999L
ledmn=2L & leddy=28L & ledyr=1999L
lstday=0L & ledday=0L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 1979 then stop,'Year out of range '
if ledyr lt 1979 then stop,'Year out of range '
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
; read data
;
        dum=findfile(dir+sdate+'.nc3')
        if dum ne '' then ncfile0=dir+sdate+'.nc3'
        rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        if iflag ne 0L then goto,jump
        tmp2=0.*p2
        for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
        index=where(mark2 lt 0.)
        if index(0) ne -1L then mark2(index)=-1.*mark2(index)/mark2(index)

      if icount eq 0L then begin
         rlat=60.
         print,alat
         read,'Enter latitude ',rlat
         index=where(abs(alat-rlat) eq min(abs(alat-rlat)))
         ilat=index(0)
         slat=strcompress(rlat,/remove_all)
         icount=1
      endif
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.
;
; longitude altitude arrays
;
      zprime=fltarr(nc,nth)
      tzm=fltarr(nc,nth)
      uzm=fltarr(nc,nth)
      markzm=fltarr(nc,nth)
      xzz=fltarr(nc,nth)
      thxz=fltarr(nc,nth)
      for k=1L,nth-1L do begin
          tzm(*,k)=tmp2(ilat,*,k)
          uzm(*,k)=u2(ilat,*,k)
          zprime(*,k)=z2(ilat,*,k)-mean(z2(ilat,*,k))
          xzz(*,k)=z2(ilat,*,k)
          markzm(*,k)=mark2(ilat,*,k)
          thxz(*,k)=th(k)
      endfor

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='xz_zprime_merra_'+sdate+'_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
       !p.thick=2.0
       !p.charsize=2.0
    endif
erase
loadct,39
xyouts,.35,.8,sdate+'  Lat= '+slat,/normal,charsize=2,charthick=2,color=0
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
index=where(xzz eq 0.)
if index(0) ne -1L then xzz(index)=0./0.
contour,tzm,alon,xzz,/noeras,xrange=[0.,360.],yrange=[30.,80.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',xticks=6,/cell_fill,c_color=col1,levels=tlevel,xtitle='Longitude',c_charthick=2,c_charsize=1.5
index=where(tlevel mod 10. eq 0)
contour,tzm,alon,xzz,levels=tlevel(index),color=0,/follow,/overplot,c_labels=1+0*index
contour,markzm,alon,xzz,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzm,alon,xzz,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
contour,thxz,alon,xzz,levels=[1500.,3000.],color=mcolor*.9,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Temperature (K)'
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
nlvls=31
col1=1+indgen(nlvls)*icolmax/nlvls
tlevel=-6000.+400.*findgen(nlvls)
zprime=zprime*1000.	; meters
index=where(zprime eq 0.)
if index(0) ne -1L then zprime(index)=0./0.
index=where(markzm eq 0.)
;if index(0) ne -1L then markzm(index)=0./0.
index=where(xzz eq 0.)
if index(0) ne -1L then xzz(index)=0./0.
contour,zprime,alon,xzz,/noeras,xrange=[0.,360.],yrange=[30.,80.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',xticks=6,/cell_fill,c_color=col1,levels=tlevel,xtitle='Longitude',c_charthick=2,c_charsize=1.5
index=where(tlevel lt 0)
contour,zprime,alon,xzz,levels=tlevel(index),color=mcolor,/follow,/overplot,c_labels=1+0*index,c_linestyle=5
index=where(tlevel gt 0)
contour,zprime,alon,xzz,levels=tlevel(index),color=0,/follow,/overplot,c_labels=1+0*index
contour,markzm,alon,xzz,levels=[0.5],color=0,/follow,/overplot,thick=5
contour,markzm,alon,xzz,levels=[-0.5],color=mcolor,/follow,/overplot,thick=5
contour,thxz,alon,xzz,levels=[1500.,3000.],color=mcolor*.9,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Eddy Height (m)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

    if setplot ne 'ps' then wait,1
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim xz_zprime_merra_'+sdate+'_'+slat+'.ps -rotate -90 xz_zprime_merra_'+sdate+'_'+slat+'.jpg'
;      spawn,'rm -f xz_zprime_merra_'+sdate+'_'+slat+'.ps'
    endif
goto, jump
end
