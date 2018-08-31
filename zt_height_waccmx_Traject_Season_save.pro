; 
; WACCM-X January
; save zt arrays. 
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
!NOERAS=-1
device,decompose=0
nxdim=700
nydim=700
xorig=[0.075,0.075,0.075,0.55,0.55,0.55]
yorig=[0.75,0.45,0.15,0.75,0.45,0.15]-0.05
xlen=0.3
ylen=0.2
cbaryoff=0.03
cbarydel=0.02
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=2002
ledmn=1
leddy=31
ledyr=2002
lstday=0
ledday=0
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1

RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
;dir='/Volumes/data/WACCM/Traject_Season.cam2.h1.'
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_Liu/WAX3548T08CO_2x.cam2.h1.'

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then begin
         save,file='zt_height_waccmx_Traject_Season_70-90N.sav',/all
         stop,'Normal Termination Condition'
      endif
      syr=string(FORMAT='(i4)',iyr)
;     syr='0001'
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
;
; read WACCM data
;
      spawn,'ls '+dir+syr+'-'+smn+'-'+sdy+'-*.nc',ncfiles
      nstep=n_elements(ncfiles)
      ksteps=kday*nstep+1L
      for ifile=0L,nstep-1L do begin
result=strsplit(ncfiles(ifile),'-',/extract)
result2=strsplit(result(3),'.',/extract)
sec=result2(0)
shrs=string(FORMAT='(f5.1)',float(sec)/3600.+((idy-1)*24.))
          sdate=syr+'-'+smn+'-'+sdy+'-'+sec
          ncfile=ncfiles(ifile)
          ncid=ncdf_open(ncfile)
          result0=ncdf_inquire(ncid)
          for idim=0,result0.ndims-1 do begin
              ncdf_diminq,ncid,idim,name,dim
              if name eq 'lon' then nc=dim
              if name eq 'lat' then nr=dim
              if name eq 'lev' then nl=dim
              if name eq 'time' then nt=dim
;             print,'read ',name,' dimension ',dim
          endfor
;
; loop over variables
;
          for ivar=0,result0.nvars-1 do begin
              result=ncdf_varinq(ncid,ivar)

              ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
              if result.name eq 'P0' then p0=data
              if result.name eq 'lat' then alat=data
              if result.name eq 'lon' then alon=data
              if result.name eq 'lev' then lev=data
              if result.name eq 'ilev' then ilev=data
              if result.name eq 'time' then time=data
              if result.name eq 'datesec' then time=data
              if result.name eq 'hyai' then hyai=data
              if result.name eq 'hybi' then hybi=data
              if result.name eq 'hyam' then hyam=data
              if result.name eq 'hybm' then hybm=data
              if result.name eq 'date' then date=data
              if result.name eq 'PS' then psfc=data	;/100.
              if result.name eq 'T' then tgrd=data
              if result.name eq 'U' then ugrd=data
              if result.name eq 'V' then vgrd=data
              if result.name eq 'AOA1' then aoa1grd=data
              if result.name eq 'AOA2' then aoa2grd=data
              if result.name eq 'OMEGA' then omegagrd=data	; Pa/s
              if result.name eq 'MSKtem' then msktemgrd=data	; TEM mask (unitless)
              if result.name eq 'e' then egrd=data
;             if result.name eq 'NOY' then noygrd=data
              if result.name eq 'NO' then nogrd=data
              if result.name eq 'CO' then cogrd=data*1.e6
;             if result.name eq 'QRL_TOT' then qrtgrd=data
;             if result.name eq 'QRS_TOT' then qrsgrd=data
;             if result.name eq 'O3' then  o3grd=data
              if result.name eq 'Z3' then  zgrd=data/1000.
              print,ivar,result.name,min(data),max(data)
          endfor
          ncdf_close,ncid
;
;============================================================
; Calculate Pressure : pgrd(i,j,k) = A(k)*PO + B(k)*PS(i,j)
;============================================================
          pgrd        = fltarr(nc,nr,nl,nt)
          Pzero       = P0      ;/100.
          FOR ilon = 0, nc-1 DO $
              FOR ilat = 0, nr-1 DO $
                  FOR ialt = 0, nl-1 DO $
                      pgrd(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat)
;
; compute atmospheric density
;
; p=rho R T -> rho=P/RT where R=287 J/K kg. Pressure in Pascals.
;
rho=pgrd/(tgrd*287.)
;
; convert from Pa/s to m/s
;
wgrd=0.*omegagrd
wgrd=-1.*omegagrd/(rho*9.81)
;wgrd=wgrd/1000.		; km/s
;wgrd=wgrd*86400.	; km/day
;
; to convert species from (NO molecules/air molecules) to NO molecules/cm3
; assume the molecular weight of one molecule of air is 29 grams (weight of O is 16, weight of N is 14, atm is 80% N2 and 20% O2)
; (mol NO/mol air) * (1 molecule air/29 grams) * (1000 g air/1 kg air) * (AIR DENSITY/m^3 air) * (Avagadros #/1 mole NO) = (molecules NO/m^3 air)
; Avagadros # = 6.022e23
;
no_conc=nogrd * (1./29.) * (1000./1.) * rho * 6.022e23
no_conc=no_conc/1.e6                                 ; divide by 1.e6 for m-3 to cm-3
;
stime=strcompress(long(time),/remove_all)
x=fltarr(nc+1)
x(0:nc-1)=alon(0:nc-1)
x(nc)=alon(0)+360.

if icount eq 0L then begin
   z0grd=aoa2grd
   zt_height=fltarr(ksteps,nl)	; time-altitude array of AOA "height" poleward of 70 N
   zt_gph=fltarr(ksteps,nl)
   zt_no=fltarr(ksteps,nl)
   zt_co=fltarr(ksteps,nl)
   zt_t=fltarr(ksteps,nl)
   zt_u=fltarr(ksteps,nl)
   zt_v=fltarr(ksteps,nl)
   zt_speed=fltarr(ksteps,nl)
   zt_omega=fltarr(ksteps,nl)
endif
if icount gt 0 then daoa2grd=aoa2grd-z0grd
speedgrd=sqrt(ugrd^2+vgrd^2)
;
; zt array
;
nhindex=where(alat ge 70)
for k=0L,nl-1L do begin
    zt_height(icount,k)=mean(aoa2grd(*,nhindex,k))
    zt_gph(icount,k)=mean(zgrd(*,nhindex,k))
    zt_no(icount,k)=mean(no_conc(*,nhindex,k))
    zt_co(icount,k)=mean(cogrd(*,nhindex,k))
    zt_t(icount,k)=mean(tgrd(*,nhindex,k))
    zt_u(icount,k)=mean(ugrd(*,nhindex,k))
    zt_v(icount,k)=mean(vgrd(*,nhindex,k))
    zt_speed(icount,k)=mean(speedgrd(*,nhindex,k))
    zt_omega(icount,k)=mean(wgrd(*,nhindex,k))
endfor
;
; zonal mean arrays
;
yz_height=fltarr(nr,nl)
yz_gph=fltarr(nr,nl)
yz_no=fltarr(nr,nl)
yz_co=fltarr(nr,nl)
yz_t=fltarr(nr,nl)
yz_u=fltarr(nr,nl)
yz_v=fltarr(nr,nl)
yz_speed=fltarr(nr,nl)
yz_w=fltarr(nr,nl)
yz_tem=fltarr(nr,nl)
for k=0L,nl-1L do begin
    for j=0L,nr-1L do begin
        yz_height(j,k)=mean(aoa2grd(*,j,k))
        yz_gph(j,k)=mean(zgrd(*,j,k))
        yz_no(j,k)=mean(no_conc(*,j,k))
        yz_co(j,k)=mean(cogrd(*,j,k))
        yz_t(j,k)=mean(tgrd(*,j,k))
        yz_u(j,k)=mean(ugrd(*,j,k))
        yz_v(j,k)=mean(vgrd(*,j,k))
        yz_speed(j,k)=mean(speedgrd(*,j,k))
        yz_w(j,k)=mean(wgrd(*,j,k))
;       yz_tem(j,k)=mean(msktemgrd(*,j,k))
    endfor
endfor
if icount eq 0 then yz_gph0=yz_gph
yz_gph=yz_gph0-yz_gph
;
; mean global altitude of pressure surfaces
;
zlev=fltarr(nl)
for k=0L,nl-1L do zlev(k)=mean(zgrd(*,*,k))

zlev=lev

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_waccmx_Traject_Season_6pan_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif

    erase
    xyouts,.25,.95,'WACCM-X '+sdate,/normal,charsize=2,color=0
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=[0.,100.,1000.,1.e4,1.e5,5.e5,1.e6,5.e6,1.e7,5.e7,1.e8,5.e8,1.e9]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,yz_no,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Latitude',title='NO',/ylog
    index=where(level gt 0)
    contour,yz_no,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(molecules/cm!u3!n)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=[0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.,1.5,2.,5.,10.,15.,20.,25.,30.,40.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,yz_co,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Latitude',title='CO',/ylog
    index=where(level gt 0)
    contour,yz_co,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(ppmv)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    imin=0.
    imax=1
    imax=170
imin=-20.
imax=20.
    nlvls=21
    level=imin+((imax-imin)/float(nlvls))*findgen(nlvls+1)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
;   contour,yz_height,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Latitude',title='dZ (Zo-Z)',/ylog
    contour,yz_gph,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',xtitle='Latitude',title='dZ',/ylog
    index=where(level gt 0)
    contour,yz_gph,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    index=where(level lt 0)
    contour,yz_gph,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=mcolor,c_linestyle=5
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(km)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=120.+10.*findgen(21)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,yz_t,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Latitude',title='Temperature',/ylog
    index=where(level gt 0)
    contour,yz_t,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(K)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(4)
    xmx=xorig(4)+xlen
    ymn=yorig(4)
    ymx=yorig(4)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=-200.+20.*findgen(14)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,yz_u,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Latitude',title='Zonal Wind',/ylog
    index=where(level gt 0)
    contour,yz_u,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    index=where(level lt 0)
    contour,yz_u,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=mcolor,c_linestyle=5
;   contour,yz_u,alat,zlev,/overplot,levels=[0],/follow,c_labels=[0],/noeras,color=0,thick=5
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(m/s)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=xorig(5)
    xmx=xorig(5)+xlen
    ymn=yorig(5)
    ymx=yorig(5)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=-1.+0.1*findgen(21)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,yz_w,alat,zlev,xrange=[-90.,90.],yrange=[max(zlev),min(zlev)],levels=level,c_color=col1,/cell_fill,/noeras,color=0,xtitle='Latitude',title='Omega',/ylog
    index=where(level gt 0)
    contour,yz_w,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=0
    index=where(level lt 0)
    contour,yz_w,alat,zlev,/overplot,levels=level(index),/follow,c_labels=0*level(index),/noeras,color=mcolor,c_linestyle=5
;   contour,yz_w,alat,zlev,/overplot,levels=[0],/follow,c_labels=[0],/noeras,color=0,thick=5
    !type=2^2+2^3+2^5
    xmnb=xmx+0.075
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yrange=[min(level),max(level)],title='(m/s)'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yz_waccmx_Traject_Season_6pan_'+sdate+'.ps -rotate -90 yz_waccmx_Traject_Season_6pan_'+sdate+'.jpg'
       spawn,'rm -f yz_waccmx_Traject_Season_6pan_'+sdate+'.ps'
    endif

icount=icount+1
endfor	; loop over time steps
goto, jump

end
