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
xorig=[0.15]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.1
cbarydel=0.01
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=2001
ledmn=1
leddy=1
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
dir='/Volumes/data/WACCM/TEM.cam2.h3.0001-01-01-00000.nc'

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
;     spawn,'ls '+dir+smn+'-'+sdy+'*.nc',ncfiles
;     nsteps=n_elements(ncfiles)
;     dlon=360./float(nsteps)
;     for istep=0L,nsteps-1L do begin
;         isec=long(60.*30.*istep)
;         ssec=string(FORMAT='(i5.5)',isec)
;         sdate=smn+'-'+sdy
          ncfile=dir
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
              if result.name eq 'P0' or result.name eq 'lat' or result.name eq 'lon' or result.name eq 'lev' or $
                 result.name eq 'date' or result.name eq 'T' then $
              ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
              if result.name eq 'P0' then p0=data
              if result.name eq 'lat' then alat=data
              if result.name eq 'lon' then alon=data
              if result.name eq 'lev' then lev=data
;             if result.name eq 'ilev' then ilev=data
;             if result.name eq 'time' then time=data
;             if result.name eq 'hyai' then hyai=data
;             if result.name eq 'hybi' then hybi=data
;             if result.name eq 'hyam' then hyam=data
;             if result.name eq 'hybm' then hybm=data
              if result.name eq 'date' then date=data
;             if result.name eq 'PS' then psfc=data/100.
              if result.name eq 'T' then t4d=data
;             if result.name eq 'U' then ugrd=data
;             if result.name eq 'V' then vgrd=data
;             if result.name eq 'CO' then cogrd=data*1.e6
;             if result.name eq 'CH4' then ch4grd=data*1.e6
;             if result.name eq 'NOY' then noygrd=data*1.e6
;             if result.name eq 'QRL_TOT' then qrtgrd=data*86400.
;             if result.name eq 'QRS_TOT' then qrsgrd=data*86400.
;             if result.name eq 'O3' then  o3grd=data*1.e6
;             if result.name eq 'Z3' then  zgrd=data/1000.

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
; loop over timesteps
;
          for n=0L,n_elements(date)-1L do begin
sdate=strcompress(date(n),/remove_all)
imon=long(strmid(sdate,2,1))
sdate1=strmid(sdate,3,2)+' '+month(imon-1)
;
; compute zonal means
;
;         ubar=fltarr(nr,nl)
;         vbar=fltarr(nr,nl)
          tbar=fltarr(nr,nl)
;         zbar=fltarr(nr,nl)
;         qbar=fltarr(nr,nl)
;         cobar=fltarr(nr,nl)
;         ch4bar=fltarr(nr,nl)
;         o3bar=fltarr(nr,nl)
;         noybar=fltarr(nr,nl)
;         noyeddy=fltarr(nr,nl)
;         zmean=fltarr(nl)
          for k=0L,nl-1L do begin
              for j=0L,nr-1L do begin
;                 ubar(j,k)=mean(ugrd(*,j,k))
;                 vbar(j,k)=mean(vgrd(*,j,k))
                  tbar(j,k)=mean(t4d(*,j,k,n))
;                 zbar(j,k)=mean(zgrd(*,j,k))
;                 cobar(j,k)=mean(cogrd(*,j,k))
;                 ch4bar(j,k)=mean(ch4grd(*,j,k))
;                 o3bar(j,k)=mean(o3grd(*,j,k))
;                 noybar(j,k)=mean(noygrd(*,j,k))
;                 noyeddy(j,k)=stdev(noygrd(*,j,k))/mean(noygrd(*,j,k))
;                 qbar(j,k)=mean(qrtgrd(*,j,k))+mean(qrsgrd(*,j,k))
;                 zmean(k)=mean(zgrd(*,*,k))
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
       device,/landscape,bits=8,filename='yz_charm_temp_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif

;   szmean=string(FORMAT='(i3)',zmean)+'km'
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    imin=120.
    imax=300.
    nlvls=16
    level=imin+10.*findgen(nlvls+1)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,tbar,alat,lev,levels=level,/cell_fill,c_color=col1,/noeras,color=0,xtitle='Latitude',ytitle='Pressure (hPa)',$
            yrange=[max(lev),1.e-4],xrange=[-90.,90.],title='WACCM4 Tbar  '+sdate1,charsize=1.5,xticks=6,/ylog
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

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yz_charm_temp_'+sdate+'.ps -rotate -90 yz_charm_temp_'+sdate+'.jpg'
       spawn,'rm -f yz_charm_temp_'+sdate+'.ps'
    endif

jumpstep:
    endfor	; loop over days
goto, jump

end
