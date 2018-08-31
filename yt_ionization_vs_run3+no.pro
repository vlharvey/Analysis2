; 
; latitude-time of e_XFRC and NO
; VLH 1/3/11
;
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
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.8
ylen=0.25
cbaryoff=0.075
cbarydel=0.01
!NOERAS=-1
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; read Xiaohua Xiaohua_Fang_e.nc file that contains all dates in one file
;
ncfile='/Users/harvey/Desktop/CHARM/Xiaohua_Fang_e_run3_smin.nc'
ncid=ncdf_open(ncfile)
result0=ncdf_inquire(ncid)
for idim=0,result0.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'lon' then ncx=dim
    if name eq 'lat' then nrx=dim
    if name eq 'altitude' then nlx=dim
    if name eq 'time' then ntx=dim
    print,'read Xiaohua ',name,' dimension ',dim
endfor
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'lat' then xlat=data
    if result.name eq 'lon' then xlon=data
    if result.name eq 'altitude' then xaltitude=data
    if result.name eq 'date' then xdate=data
    if result.name eq 'datesec' then xdatesec=data
    if result.name eq 'e' then xegrd=data
    print,ivar,result.name,min(data),max(data)
endfor
ncdf_close,ncid
;
dir='/Volumes/data/WACCM/no_aur_run3_smin.cam2.h1.'
rtd=double(180./!pi)
dtr=1./rtd
ks=1.931853d-3
ecc=0.081819
gamma45=9.80
;
; loop over WACCM date+datesec files
;
for n=0,ntx-1L do begin
    xegrd3d=reform(xegrd(*,*,*,n))
;
; zonal mean
;
    xezm=fltarr(nrx,nlx)
    for j=0,nrx-1L do $
        for k=0,nlx-1L do $
            xezm(j,k)=mean(xegrd3d(*,j,k))
;
; extract date
;
    date0=strcompress(xdate(n),/remove_all)
    datesec0=xdatesec(n)
    while datesec0 gt 84600 do datesec0=datesec0-84600
    datesec0=string(FORMAT='(i5.5)',datesec0)
    syr=strmid(date0,0,4)
    smn=strmid(date0,4,2)
    sdy=strmid(date0,6,2)
    print,date0,' ',datesec0
;
; read WACCM data
;
    ifile=dir+syr+'-'+smn+'-'+sdy+'-'+datesec0+'.nc'
    sdum=findfile(ifile)
;   if sdum eq '' then stop
    if sdum eq '' then goto,jumpstep
    sdate=syr+'-'+smn+'-'+sdy+'-'+datesec0
    print,'Found WACCM file on '+sdate
    ncid=ncdf_open(ifile)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'lon' then nc=dim
        if name eq 'lat' then nr=dim
        if name eq 'lev' then nl=dim
        if name eq 'time' then nt=dim
;       print,'read ',name,' dimension ',dim
    endfor
;
; loop over variables
;
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)

        if result.name eq 'P0' or result.name eq 'lat' or result.name eq 'lon' or result.name eq 'lev' or result.name eq 'hyai' or $
           result.name eq 'hybi' or result.name eq 'hyam' or result.name eq 'hybm' or result.name eq 'date' or $
           result.name eq 'PS' or result.name eq 'e_XFRC' or result.name eq 'Z3' or result.name eq 'NO' or result.name eq 'T' then begin
           ncdf_varget,ncid,ncdf_varid(ncid,result.name),data

           if result.name eq 'P0' then p0=data
           if result.name eq 'lat' then lat=data
           if result.name eq 'lon' then lon=data
           if result.name eq 'lev' then lev=data
           if result.name eq 'hyai' then hyai=data
           if result.name eq 'hybi' then hybi=data
           if result.name eq 'hyam' then hyam=data
           if result.name eq 'hybm' then hybm=data
           if result.name eq 'date' then date=data
           if result.name eq 'PS' then psfc=data	;/100.
           if result.name eq 'e_XFRC' then exfrcgrd=data
           if result.name eq 'Z3' then ggrd=data/1000.
           if result.name eq 'NO' then nogrd=data
           if result.name eq 'T' then tgrd=data

;          print,ivar,result.name,min(data),max(data)
        endif
    endfor
    ncdf_close,ncid
print,'Input ',min(xegrd3d),max(xegrd3d)
print,'WACCM ',min(exfrcgrd),max(exfrcgrd)
;
; convert geopotential to geometric height
;  
    zgrd=0.*ggrd
    for k=0L,nl-1L do begin
        for j=0L,nr-1L do begin
            sin2=sin( (lat(j)*dtr)^2.0 )
            numerator=1.0+ks*sin2
            denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
            gammas=gamma45*(numerator/denominator)
            r=6378.137/(1.006803-(0.006706*sin2))
            zgrd(*,j,k)=(r*ggrd(*,j,k))/ ( (gammas/gamma45)*r - ggrd(*,j,k) )
        endfor
    endfor
;
; Calculate Pressure : pgrd(i,j,k) = A(k)*PO + B(k)*PS(i,j)
;
    pgrd        = fltarr(nc,nr,nl)
    Pzero       = P0      ;/100.
    FOR ilon = 0, nc-1 DO $
        FOR ilat = 0, nr-1 DO $
            FOR ialt = 0, nl-1 DO $
                pgrd(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat)
;
; convert NO to molecules/cm3
;
    rho=pgrd/(tgrd*287.)
    no_conc=nogrd * (1./29.) * (1000./1.) * rho * 6.022e23
    no_conc=no_conc/1.e6                                 ; divide by 1.e6 for m-3 to cm-3
    nogrd=no_conc
;
; WACCM zonal means
;
    exfrczm=fltarr(nr,nl)
    nozm=fltarr(nr,nl)
    zmean=fltarr(nl)
    for j=0L,nr-1L do begin
        for k=0,nl-1L do begin
            exfrczm(j,k)=mean(exfrcgrd(*,j,k))
            nozm(j,k)=mean(nogrd(*,j,k))
        endfor
    endfor
    for k=0,nl-1L do zmean(k)=mean(zgrd(*,*,k))
    if n eq 0L then begin
       noyt=fltarr(ntx,nr)
       exfrcyt=fltarr(ntx,nr)
       sdate_all=strarr(ntx)
       rlev=0.
       print,zmean
       read,'Enter desired altitude ',rlev
       index=where(abs(rlev-zmean) eq min(abs(rlev-zmean)))
       ilev=index(0)
       slev=strcompress(long(zmean(ilev)),/remove_all)
       dum=findfile('yt_ionization_vs_run3+no_'+slev+'_km.sav')
       if dum(0) ne '' then goto,plotit
    endif
    noyt(n,*)=nozm(*,ilev)
    exfrcyt(n,*)=exfrczm(*,ilev)
    sdate_all(n)=sdate

    jumpstep:
endfor          ; loop over time steps
save,file='yt_ionization_vs_run3+no_'+slev+'_km.sav',ntx,lat,noyt,exfrcyt,sdate_all,slev
plotit:
restore,'yt_ionization_vs_run3+no_'+slev+'_km.sav'
;
; plot
;
; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yt_ionization_vs_run3+no_'+slev+'_km.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif
;
; set zeros to NaNs for plotting
;
index=where(noyt eq 0.)
if index(0) ne -1L then noyt(index)=0./0.
index=where(exfrcyt eq 0.)
if index(0) ne -1L then exfrcyt(index)=0./0.
noyt=smooth(noyt,3,/NaN)
exfrcyt=smooth(exfrcyt,3,/NaN)

index=where(strmid(sdate_all,11,5) eq '00000',nxticks)
xlab=strmid(sdate_all(index),5,5)
    erase
    xyouts,.4,.925,'Run 3            '+slev+' km',/normal,color=0,charsize=2.5
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    nlvls=21
    level=200.*findgen(nlvls)
    imin=min(level)
    imax=max(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,exfrcyt,1+findgen(ntx),lat,levels=level,/cell_fill,c_color=col1,/noeras,color=0,ytitle='Latitude',$
            xrange=[1,ntx],yrange=[-90.,90.],title='e_XFRC',charsize=1.5,yticks=6,xticks=nxticks-1,xtickname=xlab
;   contour,exfrcyt,1+findgen(ntx),lat,levels=level,/overplot,color=0,/noeras
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,min(xorig),xorig(1)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='(ions/cm3/s)',charsize=1.5,color=0
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
nomin=2.e7
nomax=3.e8
    nlvls=21
    nolevel=nomin+((nomax-nomin)/float(nlvls))*findgen(nlvls)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,noyt,1+findgen(ntx),lat,levels=nolevel,/cell_fill,c_color=col1,/noeras,color=0,ytitle='Latitude',$
            xrange=[1.,ntx],yrange=[-90.,90.],title='NO',charsize=1.5,yticks=6,xticks=nxticks-1,xtickname=xlab
;   contour,noyt,1+findgen(ntx),lat,levels=nolevel,/overplot,color=0,/noeras

    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,min(xorig),xorig(1)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[nomin,nomax],[0,0],yrange=[0,10],xrange=[nomin,nomax],xtitle='(molec/cm3)',charsize=1.5,color=0
    ybox=[0,10,10,0,0]
    x1=nomin
    dx=(nomax-nomin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yt_ionization_vs_run3+no_'+slev+'_km.ps -rotate -90 yt_ionization_vs_run3+no_'+slev+'_km.jpg'
       spawn,'rm -f yt_ionization_vs_run3+no_'+slev+'_km.ps'
    endif
end
