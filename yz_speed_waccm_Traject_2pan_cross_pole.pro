; 
; CO
; cross polar plots and polar plots
; WACCM files from with AOA1 set to latitude and AOA2 set to altitude
; polar plots each timestep
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
!NOERAS=-1
device,decompose=0
nxdim=700
nydim=700
xorig=[0.075,0.525,0.075,0.525]
yorig=[0.6,0.6,0.15,0.15]
xlen=0.3
ylen=0.3
cbaryoff=0.03
cbarydel=0.02
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=2002
ledmn=3
leddy=29
ledyr=2002
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (monl, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (monl, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

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
monl=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
dir='/Volumes/data/WACCM/Traject_Season.cam2.h1.'

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
      if ndays gt ledday then stop,' Normal termination condition '
;     syr=string(FORMAT='(i4)',iyr)
      syr='0001'
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
;
; read WACCM data
;
      spawn,'ls '+dir+syr+'-'+smn+'-'+sdy+'-*.nc',ncfiles
      for ifile=0L,n_elements(ncfiles)-1L do begin
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
              if result.name eq 'CH4' or result.name eq 'H' or result.name eq 'H2O' or result.name eq 'NO2' or $
                 result.name eq 'NOX' or result.name eq 'NOY' or result.name eq 'QRS_TOT' or result.name eq 'QRL_TOT' or $
                 result.name eq 'O3' or result.name eq 'VTH2d' or result.name eq 'WTH2d' or result.name eq 'UV2d' or result.name eq 'UW2d' or $
                 result.name eq 'OH' or result.name eq 'HO2' or result.name eq 'MSKtem' or result.name eq 'N2O' then goto,jumpvar
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
              if result.name eq 'e' then egrd=data
;             if result.name eq 'NOY' then noygrd=data
              if result.name eq 'NO' then nogrd=data
              if result.name eq 'CO' then cogrd=data*1.e6
;             if result.name eq 'QRL_TOT' then qrtgrd=data
;             if result.name eq 'QRS_TOT' then qrsgrd=data
;             if result.name eq 'O3' then  o3grd=data
              if result.name eq 'Z3' then  zgrd=data/1000.

              print,ivar,result.name,min(data),max(data)
jumpvar:
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

    q1=reform(no_conc(*,*,0))      ; NO at top level
    q=0.*fltarr(nc+1,nr)
    q(0:nc-1,0:nr-1)=q1(0:nc-1,0:nr-1)
    q(nc,*)=q(0,*)

if icount eq 0L then begin
   nomin=fltarr(nl)
   nomax=fltarr(nl)
   emin=fltarr(nl)
   emax=fltarr(nl)
   height0grd=aoa2grd
   z0grd=zgrd
endif
if icount gt 0 then aoa2grd=aoa2grd-height0grd

speedgrd=sqrt(ugrd^2+vgrd^2)

;
; zonal mean AOA2=Altitude
;
;aoa2yz=0.*fltarr(nr,nl)
;for k=0L,nl-1L do begin
;for j=0L,nr-1L do begin
;    aoa2yz(j,k)=mean(aoa2grd(*,j,k))
;endfor
;endfor
;
; mean global altitude of pressure surfaces
;
zlev=fltarr(nl)
for k=0L,nl-1L do zlev(k)=mean(zgrd(*,*,k))
;
; user enters desired longitude
;
if icount eq 0L then begin
;   print,alon
    rlon1=180.
;   read,' Enter longitude ',rlon1
    index1=where(alon eq rlon1)
    if index1(0) eq -1 then stop,'Bad Longitude'
    ilon1=index1(0)
    rlon2=rlon1+180.
    if rlon2 gt max(alon) then rlon2=rlon2-360.
    index2=where(alon eq rlon2)
    ilon2=index2(0)
    slon1=' '+strcompress(string(rlon1),/remove_all)+' E'
    slon2=' '+strcompress(string(rlon2),/remove_all)+' E'
    print,'longitudes ',rlon1,rlon2
    x=fltarr(nc+1)
    x(0:nc-1)=alon(0:nc-1)
    x(nc)=alon(0)+360.
    x2d=fltarr(nc+1,nr)
    y2d=fltarr(nc+1,nr)
    for i=0,nc do y2d(i,*)=alat
    for j=0,nr-1 do x2d(*,j)=x
    xyz=fltarr(nr,nl)
    yyz=fltarr(nr,nl)
    for i=0,nr-1 do yyz(i,*)=zlev
    for j=0,nl-1 do xyz(*,j)=alat
endif
;
; cross pole altitude arrays
;
uyz=fltarr(nr,nl)
tyz=fltarr(nr,nl)
noyz=fltarr(nr,nl)
coyz=fltarr(nr,nl)
speedyz=fltarr(nr,nl)
aoa1yz=fltarr(nr,nl)
aoa2yz=fltarr(nr,nl)
for k=0,nl-1 do begin
    uyz(0:nr/2-1,k)=ugrd(ilon1,nr/2:nr-1,k)
    uyz(nr/2:nr-1,k)=ugrd(ilon2,nr/2:nr-1,k)
    uyz(nr/2:nr-1,k)=reverse(uyz(nr/2:nr-1,k))

    tyz(0:nr/2-1,k)=tgrd(ilon1,nr/2:nr-1,k)
    tyz(nr/2:nr-1,k)=tgrd(ilon2,nr/2:nr-1,k)
    tyz(nr/2:nr-1,k)=reverse(tyz(nr/2:nr-1,k))

    noyz(0:nr/2-1,k)=no_conc(ilon1,nr/2:nr-1,k)
    noyz(nr/2:nr-1,k)=no_conc(ilon2,nr/2:nr-1,k)
    noyz(nr/2:nr-1,k)=reverse(noyz(nr/2:nr-1,k))

    coyz(0:nr/2-1,k)=cogrd(ilon1,nr/2:nr-1,k)
    coyz(nr/2:nr-1,k)=cogrd(ilon2,nr/2:nr-1,k)
    coyz(nr/2:nr-1,k)=reverse(coyz(nr/2:nr-1,k))

    speedyz(0:nr/2-1,k)=speedgrd(ilon1,nr/2:nr-1,k)
    speedyz(nr/2:nr-1,k)=speedgrd(ilon2,nr/2:nr-1,k)
    speedyz(nr/2:nr-1,k)=reverse(speedyz(nr/2:nr-1,k))

    aoa1yz(0:nr/2-1,k)=aoa1grd(ilon1,nr/2:nr-1,k)
    aoa1yz(nr/2:nr-1,k)=aoa1grd(ilon2,nr/2:nr-1,k)
    aoa1yz(nr/2:nr-1,k)=reverse(aoa1yz(nr/2:nr-1,k))

    aoa2yz(0:nr/2-1,k)=aoa2grd(ilon1,nr/2:nr-1,k)
    aoa2yz(nr/2:nr-1,k)=aoa2grd(ilon2,nr/2:nr-1,k)
    aoa2yz(nr/2:nr-1,k)=reverse(aoa2yz(nr/2:nr-1,k))
endfor

;for ilev=22L,22L do begin	;nl-1L do begin
for ilev=9L,9L do begin	;nl-1L do begin
    slev=string(FORMAT='(f13.6)',lev(ilev))+'hPa'
    slev0=string(format='(i2.2)',ilev)
print,sdate,' ',slev

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_speed_'+sdate+'_'+slev0+'_waccm_Traject_2pan_cross_pole.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif

    zavg=mean(zgrd(*,nr/2:nr-1,ilev))
    szavg=string(FORMAT='(i3)',zavg)+'km'
    no1=reform(no_conc(*,*,ilev))
    no=0.*fltarr(nc+1,nr)
    no(0:nc-1,0:nr-1)=no1(0:nc-1,0:nr-1)
    no(nc,*)=no(0,*)
    co1=reform(cogrd(*,*,ilev))
    co=0.*fltarr(nc+1,nr)
    co(0:nc-1,0:nr-1)=co1(0:nc-1,0:nr-1)
    co(nc,*)=co(0,*)
    e1=reform(egrd(*,*,ilev))
    e=0.*fltarr(nc+1,nr)
    e(0:nc-1,0:nr-1)=e1(0:nc-1,0:nr-1)
    e(nc,*)=e(0,*)
    u1=reform(ugrd(*,*,ilev))
    u=0.*fltarr(nc+1,nr)
    u(0:nc-1,0:nr-1)=u1(0:nc-1,0:nr-1)
    u(nc,*)=u(0,*)
    v1=reform(vgrd(*,*,ilev))
    v=0.*fltarr(nc+1,nr)
    v(0:nc-1,0:nr-1)=v1(0:nc-1,0:nr-1)
    v(nc,*)=v(0,*)
    aoa1a=reform(aoa1grd(*,*,ilev))
    aoa2a=reform(aoa2grd(*,*,ilev))
    aoa1=0.*fltarr(nc+1,nr)
    aoa1(0:nc-1,0:nr-1)=aoa1a(0:nc-1,0:nr-1)
    aoa1(nc,*)=aoa1(0,*)
    aoa2=0.*fltarr(nc+1,nr)
    aoa2(0:nc-1,0:nr-1)=aoa2a(0:nc-1,0:nr-1)
    aoa2(nc,*)=aoa2(0,*)

    z1=reform(zgrd(*,*,ilev))
    z=0.*fltarr(nc+1,nr)
    z(0:nc-1,0:nr-1)=z1(0:nc-1,0:nr-1)
    z(nc,*)=z(0,*)
    speed=sqrt(u^2+v^2)
print,'Max Speed at ',lev(ilev),zavg,max(speed)

    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    xyouts,.15,.95,'WACCM-4 '+sdate+'     '+slev+' ('+szavg+')',/normal,color=0,charsize=1.5
    MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,title='Wind Speed',charsize=1.5,color=0,/noborder
    nlvls=25
    colevel=5.+5.*findgen(nlvls)
nlvls=n_elements(colevel)
    col1=1+indgen(nlvls)*icolmax/nlvls
imin=min(colevel)
imax=max(colevel)
    contour,speed,x,alat,/overplot,levels=colevel,c_color=col1,/cell_fill,/noeras,title='Wind Speed'
    contour,speed,x,alat,/overplot,levels=colevel,/follow,c_labels=1+0*colevel,/noeras,color=0
    contour,speed,x,alat,/overplot,levels=[50.],/follow,c_labels=[0],/noeras,color=mcolor,thick=5
;   contour,z,x,alat,/overplot,levels=10+findgen(150),/follow,c_labels=0*findgen(150),/noeras,color=mcolor
    oplot,rlon1+0*findgen(90),findgen(91),psym=2,color=mcolor*.9,symsize=0.75
    oplot,rlon2+0*findgen(90),findgen(91),psym=2,color=mcolor*.9,symsize=0.75
    MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,color=mcolor
;   contour,q,x,alat,/overplot,levels=1.8e-15,/follow,c_labels=0,/noeras,color=mcolor,thick=5
;   drawvectors,nc+1,nr,x,alat,u,v,20,1
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.05,xmx-0.05,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(m/s)',charsize=1.,color=0
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
    contour,speedyz,alat,lev,yrange=[max(lev),min(lev)],levels=colevel,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',title='Wind Speed',/ylog,$
            xtitle=slon1+'       Latitude      '+slon2,xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
    contour,speedyz,alat,lev,/overplot,levels=10.+10.*findgen(nlvls),/follow,c_labels=0*colevel,/noeras,color=0
    contour,speedyz,alat,lev,/overplot,levels=[50.],/follow,c_labels=[0],/noeras,color=mcolor,thick=5
    oplot,alat,lev(ilev)+0.*alat,psym=2,color=mcolor*.9,symsize=0.75
    ymnb=ymn-cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.05,xmx-0.05,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],charsize=1.,color=0,xtitle='(m/s)'
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
;   MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,title='CO',charsize=1.5,color=0,/noborder
;colevel=[0.001,0.005,0.01,0.05,0.1,0.25,0.5,0.75,1.+.25*findgen(10)]
;colevel=[0.001,0.005,0.01,0.05,0.1,0.25,0.5,0.75,1.+.5*findgen(10)]
;colevel=1.+findgen(15)
;imin=min(colevel)
;imax=max(colevel)
;    nlvls=n_elements(colevel)
;    col1=1+indgen(nlvls)*icolmax/nlvls
;    contour,co,x,alat,/overplot,levels=colevel,c_color=col1,/cell_fill,/noeras
;index=where(colevel gt 0.)
;if index(0) ne -1L then contour,co,x,alat,/overplot,levels=colevel(index),/follow,c_labels=0*colevel(index),/noeras,color=0
;    MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,color=mcolor

;colevel=0.05*findgen(21)
;imin=min(colevel)
;imax=max(colevel)
imin=min(aoa2yz)
imax=max(aoa2yz)
colevel=imin+((imax-imin)/20.)*findgen(21)
nlvls=n_elements(colevel)
col1=1+indgen(nlvls)*icolmax/nlvls
    contour,aoa2yz,alat,lev,yrange=[max(lev),min(lev)],levels=colevel,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',title='AOA2=Altitude-Z0',/ylog,$
            xtitle=slon1+'       Latitude      '+slon2,xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
index=where(colevel gt 0.)
    contour,aoa2yz,alat,lev,/overplot,levels=colevel(index),/follow,c_labels=0*colevel(index),/noeras,color=0
index=where(colevel lt 0.)
    contour,aoa2yz,alat,lev,/overplot,levels=colevel(index),/follow,c_labels=0*colevel(index),/noeras,color=mcolor,c_linestyle=5

    ymnb=ymn-cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.05,xmx-0.05,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],charsize=1.,color=0	;,xtitle='(ppmv)'
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
colevel=[0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.5,0.75,1.,1.25,1.5,1.75,2.+2.*findgen(20)]
imin=min(colevel)
imax=max(colevel)
    nlvls=n_elements(colevel)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,coyz,alat,lev,yrange=[max(lev),min(lev)],levels=colevel,c_color=col1,/cell_fill,/noeras,color=0,ytitle='Pressure (hPa)',title='CO',/ylog,$
            xtitle=slon1+'       Latitude      '+slon2,xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
index=where(colevel gt 0)
    contour,coyz,alat,lev,/overplot,levels=colevel(index),/follow,c_labels=0*colevel(index),/noeras,color=0
    ymnb=ymn-cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.05,xmx-0.05,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],charsize=1.,color=0,xtitle='(ppmv)'
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
       spawn,'convert -trim yz_speed_'+sdate+'_'+slev0+'_waccm_Traject_2pan_cross_pole.ps -rotate -90 yz_speed_'+sdate+'_'+slev0+'_waccm_Traject_2pan_cross_pole.jpg'
       spawn,'rm -f yz_speed_'+sdate+'_'+slev0+'_waccm_Traject_2pan_cross_pole.ps'
    endif
jumplev:
endfor	; loop over pressure
icount=1
endfor	; loop over time steps
goto, jump

end
