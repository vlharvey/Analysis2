;
; noaurfco run
; store T85-T60 and U65 for all levs and all days
; VLH 20090614
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
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15,0.60]
yorig=[0.275,0.275]
xlen=0.325
ylen=0.5
cbaryoff=0.10
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif

idir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_Mills/noaurfco/noaurfco.vars.h3.'
;idir='/storage/ptmp/mills/waccm/epp/run/noaurfco/noaurfco.vars.h3.'
;
; spawn directory contents
;
spawn,'ls '+idir+'*nc',ifiles
nfile=n_elements(ifiles)
;
; loop over netcdf files
;
for ifile=0L,nfile-1L do begin
;
; open and read netcdf data
;
ncfile=ifiles(ifile)
print,ncfile
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
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    if result.name eq 'lat' or result.name eq 'lon' or result.name eq 'lev' or $
       result.name eq 'date' or result.name eq 'T' or result.name eq 'U' then $
       ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
;   if result.name eq 'P0' then P0=float(data)/100.             ; Reference pressure (Pa)
;   if result.name eq 'hyam' then hyam=float(data)
;   if result.name eq 'hybm' then hybm=float(data)
    if result.name eq 'lat' then alat=float(data)
    if result.name eq 'lon' then alon=float(data)
    if result.name eq 'lev' then lev=float(data)
    if result.name eq 'date' then yyyymmdd=data			; Date (YYYYMMDD)
;   if result.name eq 'PS' then psfc=float(data)/100.           ; Surface pressure (Pa)
    if result.name eq 'T' then temp_4d=float(data)		; Temperature
    if result.name eq 'U' then uwind_4d=float(data)		; Zonal wind
    print,'read ',result.name,' variable'
endfor  ; loop over variables
ncdf_close,ncid
;
; extract ncfile year from original filename
;
result=strsplit(ncfile,'.',/extract)
syear=result(3)
;
; quick: extract from native latitude grid T86-T62 and U66 (should be 60, 65, and 85)
;
n85=where(abs(alat-85.) eq min(abs(alat-85.)))
n65=where(abs(alat-65.) eq min(abs(alat-65.)))
n65=n65(0)
n60=where(abs(alat-60.) eq min(abs(alat-60.)))
n60=n60(1)	; 62N
print,'NH lats ',alat(n60),alat(n65),alat(n85)
s85=where(abs(alat+85.) eq min(abs(alat+85.)))
s65=where(abs(alat+65.) eq min(abs(alat+65.)))
s65=s65(0)
s60=where(abs(alat+60.) eq min(abs(alat+60.)))
s60=s60(0)      ; 62S
print,'SH lats ',alat(s60),alat(s65),alat(s85)
;
; declare arrays
;
nhdt=fltarr(nt,nl)
shdt=fltarr(nt,nl)
nhu60=fltarr(nt,nl)
shu60=fltarr(nt,nl)
;
; loop over days
;
for itime=0L,nt-1L do begin
    print,yyyymmdd(itime)
;
; extract 3-D grid from 4-d arrays
;
    tgrd=reform(temp_4d(*,*,*,itime))
    ugrd=reform(uwind_4d(*,*,*,itime))
;
; zonal mean T, U
;
    tzm=fltarr(nr,nl)
    uzm=fltarr(nr,nl)
    for k=0L,nl-1L do begin
        for j=0L,nr-1L do begin
            tzm(j,k)=total(tgrd(*,j,k))/float(nc)
            uzm(j,k)=total(ugrd(*,j,k))/float(nc)
        endfor
    endfor
;
; meridional temperature gradient and zonal mean wind at 60: T85-T60 and U65
;
    for k=0L,nl-1L do begin
        nhdt(itime,k)=tzm(n85,k)-tzm(n60,k)
        shdt(itime,k)=tzm(s85,k)-tzm(s60,k)
        nhu60(itime,k)=uzm(n65,k)
        shu60(itime,k)=uzm(s65,k)
    endfor
;
; plot daily zonal mean T and U (tzm and uzm)
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=11
col1=1+indgen(nlvls)*icolmax/nlvls
level=140.+15.*findgen(nlvls)
syyyymmdd=strcompress(yyyymmdd(itime),/remove_all)
if setplot eq 'ps' then begin
   set_plot,'ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='store_ssw_stats_noaurfco_mills_'+syyyymmdd+'.ps'
endif
xyouts,.4,.825,syyyymmdd,/normal,color=0,charsize=2
contour,tzm,alat,lev,/noeras,xrange=[-90,90],yrange=[100.,0.0001],color=0,/ylog,$
       ytitle='Pressure (hPa)',xtitle='Latitude',/fill,c_color=col1,levels=level,title='Temperature'
index=where(level lt 180.)
contour,tzm,alat,lev,levels=level(index),color=mcolor,/follow,/overplot,min_value=-9999.,c_linestyle=5,thick=3
index=where(level gt 250.)
contour,tzm,alat,lev,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,thick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1
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
nlvls=11
col1=1+indgen(nlvls)*icolmax/nlvls
level=-100.+20.*findgen(nlvls)
contour,uzm,alat,lev,/noeras,xrange=[-90,90],yrange=[100.,0.0001],color=0,/ylog,$
        xtitle='Latitude',/fill,c_color=col1,levels=level,title='Zonal Wind'
index=where(level lt 0.)
contour,uzm,alat,lev,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_linestyle=5
index=where(level gt 0.)
contour,uzm,alat,lev,levels=level(index),color=mcolor,/follow,/overplot,min_value=-9999.
imin=min(level)
imax=max(level)
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
;
; Close PostScript file and return control to X-windows
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim store_ssw_stats_noaurfco_mills_'+syyyymmdd+'.ps -rotate -90 '+$
                            'store_ssw_stats_noaurfco_mills_'+syyyymmdd+'.jpg'
   spawn,'/usr/bin/rm store_ssw_stats_noaurfco_mills_'+syyyymmdd+'.ps'
endif
endfor	; loop over days
;
; save yearly file
;
save,file='Datfiles/waccm_ssw_stats_noaurfco_mills_'+syear+'.sav',yyyymmdd,lev,nhdt,shdt,nhu60,shu60
endfor	; loop over files
end
