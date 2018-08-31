; 
; daily average WACCM zonal means in support of CHARM.
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
xorig=[0.1,0.4,0.7,0.1,0.4,0.7]
yorig=[0.6,0.6,0.6,0.15,0.15,0.15]
xlen=0.2
ylen=0.3
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=2001
ledmn=1
leddy=7
ledyr=2001
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
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
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_CHARM/wa4_charm_30min.cam2.h2.0001-'

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
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
;
; read WACCM data
;
      spawn,'ls '+dir+smn+'-'+sdy+'*.nc',ncfiles
      nsteps=n_elements(ncfiles)
      dlon=360./float(nsteps)
      for istep=0L,nsteps-1L do begin
          isec=long(60.*30.*istep)
          ssec=string(FORMAT='(i5.5)',isec)
          sdate=smn+'-'+sdy
          ncfile=ncfiles(istep)
          ncid=ncdf_open(ncfile)
          result0=ncdf_inquire(ncid)
          for idim=0,result0.ndims-1 do begin
              ncdf_diminq,ncid,idim,name,dim
              if name eq 'lon' then nc=dim
              if name eq 'lat' then nr=dim
              if name eq 'lev' then nl=dim
              if name eq 'time' then nt=dim
              print,'read ',name,' dimension ',dim
          endfor
;
; loop over variables
;
          for ivar=0,result0.nvars-1 do begin
              result=ncdf_varinq(ncid,ivar)
;             if result.name ne 'T' and result.name ne 'U' and result.name ne 'V' then $

              ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
              if result.name eq 'P0' then p0=data
              if result.name eq 'lat' then alat=data
              if result.name eq 'lon' then alon=data
              if result.name eq 'lev' then lev=data
              if result.name eq 'ilev' then ilev=data
              if result.name eq 'time' then time=data
              if result.name eq 'hyai' then hyai=data
              if result.name eq 'hybi' then hybi=data
              if result.name eq 'hyam' then hyam=data
              if result.name eq 'hybm' then hybm=data
              if result.name eq 'date' then date=data
              if result.name eq 'PS' then psfc=data/100.
              if result.name eq 'T' then tgrd=data
              if result.name eq 'U' then ugrd=data
              if result.name eq 'V' then vgrd=data
              if result.name eq 'CO' then cogrd=data*1.e6
              if result.name eq 'CH4' then ch4grd=data*1.e6
              if result.name eq 'NOY' then noygrd=data*1.e6
              if result.name eq 'QRL_TOT' then qrtgrd=data*86400.
              if result.name eq 'QRS_TOT' then qrsgrd=data*86400.
              if result.name eq 'O3' then  o3grd=data*1.e6
              if result.name eq 'Z3' then  zgrd=data/1000.

              print,ivar,result.name,min(data),max(data)
          endfor
          ncdf_close,ncid
;
;============================================================
; Calculate Pressure : pgrd(i,j,k) = A(k)*PO + B(k)*PS(i,j)
;============================================================
;         pgrd        = fltarr(nc,nr,nl)
;         Pzero       = P0/100.
;         FOR ilon = 0, nc-1 DO $
;             FOR ilat = 0, nr-1 DO $
;                 FOR ialt = 0, nl-1 DO $
;                     pgrd(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat,itime)
;
; compute zonal means
;
          ubar=fltarr(nr,nl)
          vbar=fltarr(nr,nl)
          tbar=fltarr(nr,nl)
          zbar=fltarr(nr,nl)
          qbar=fltarr(nr,nl)
          cobar=fltarr(nr,nl)
          ch4bar=fltarr(nr,nl)
          o3bar=fltarr(nr,nl)
          noybar=fltarr(nr,nl)
          noyeddy=fltarr(nr,nl)
          zmean=fltarr(nl)
          for k=0L,nl-1L do begin
              for j=0L,nr-1L do begin
                  ubar(j,k)=mean(ugrd(*,j,k))
                  vbar(j,k)=mean(vgrd(*,j,k))
                  tbar(j,k)=mean(tgrd(*,j,k))
                  zbar(j,k)=mean(zgrd(*,j,k))
                  cobar(j,k)=mean(cogrd(*,j,k))
                  ch4bar(j,k)=mean(ch4grd(*,j,k))
                  o3bar(j,k)=mean(o3grd(*,j,k))
                  noybar(j,k)=mean(noygrd(*,j,k))
                  noyeddy(j,k)=stdev(noygrd(*,j,k))/mean(noygrd(*,j,k))
                  qbar(j,k)=mean(qrtgrd(*,j,k))+mean(qrsgrd(*,j,k))
                  zmean(k)=mean(zgrd(*,*,k))
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
       device,/landscape,bits=8,filename='yz_charm_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif

    szmean=string(FORMAT='(i3)',zmean)+'km'
    erase
    xyouts,.35,.95,'WACCM '+sdate,/normal,color=0,charsize=2
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    imin=120.
    imax=300.
    nlvls=20
    level=imin+((imax-imin)/float(nlvls))*findgen(nlvls+1)
    level=[level,400.,500.,600.,700.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
lev=zmean
    contour,tbar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='Tbar',charsize=1.5,xticks=6
    contour,tbar,alat,lev,levels=level,/overplot,color=0,/noeras
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(K)',charsize=1,color=0
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
    imin=-120.
    imax=120.
    nlvls=12
    level=imin+((imax-imin)/float(nlvls))*findgen(nlvls+1)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,ubar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='Ubar',charsize=1.5,xticks=6
    index=where(level gt 0)
    contour,ubar,alat,lev,levels=level(index),/overplot,color=0,/noeras
    index=where(level lt 0)
    contour,ubar,alat,lev,levels=level(index),/overplot,color=mcolor,c_linestyle=5,/noeras
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(m/s)',charsize=1,color=0
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
    imin=-30.
    imax=30.
    nlvls=20
    level=imin+((imax-imin)/float(nlvls))*findgen(nlvls+1)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,vbar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='Vbar',charsize=1.5,xticks=6
    index=where(level gt 0)
    contour,vbar,alat,lev,levels=level(index),/overplot,color=0,/noeras
    index=where(level lt 0)
    contour,vbar,alat,lev,levels=level(index),/overplot,color=mcolor,c_linestyle=5,/noeras
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(m/s)',charsize=1,color=0
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
    level=[0.001,0.005,0.01,0.05,0.1,0.5,1.,2.,3.,4.,5.,10.,15.,20.,30.,40.,50.,75.,100.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,cobar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='CO bar',charsize=1.5,xticks=6
    contour,cobar,alat,lev,levels=level,/overplot,color=0,/noeras
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(ppmv)',charsize=1,color=0
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
    level=[0.0,1.,5.,10.,20.,30.,40.,50.,100.,150.,200.,300.,400.,500.,600.,800.,1000.,1200.,1400.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,noybar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='NOy bar',charsize=1.5,xticks=6
    contour,noybar,alat,lev,levels=level,/overplot,color=0,/noeras
    contour,noyeddy,alat,lev,levels=.2*findgen(30),/overplot,color=mcolor,/noeras,thick=2
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(ppmv)',charsize=1,color=0
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
    level=findgen(12)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,o3bar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,$
              yrange=[min(lev),max(lev)],xrange=[-90.,90.],title='O3 bar',charsize=1.5,xticks=6
    contour,o3bar,alat,lev,levels=level,/overplot,color=0,/noeras
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(ppmv)',charsize=1,color=0
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
       spawn,'convert -trim yz_charm_'+sdate+'.ps -rotate -90 yz_charm_'+sdate+'.jpg'
       spawn,'rm -f yz_charm_'+sdate+'.ps'
    endif

jumpstep:
noyold=noygrd

    endfor	; loop over daily timesteps
goto, jump

end
