;
; plot vortex minimum temperature as a function of altitude and time
;
@rd_ukmo

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.35]
xlen=0.8
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
dir='/aura3/data/UKMO_data/Datfiles/ppassm_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
ifiles=[$
'ukmo_nh_files_91-92.fil',$
'ukmo_nh_files_92-93.fil',$
'ukmo_nh_files_93-94.fil',$
'ukmo_nh_files_94-95.fil',$
'ukmo_nh_files_95-96.fil',$
'ukmo_nh_files_96-97.fil']
;'ukmo_nh_files_97-98.fil',$
;'ukmo_nh_files_98-99.fil',$
;'ukmo_nh_files_99-00.fil',$
;'ukmo_nh_files_00-01.fil',$
;'ukmo_nh_files_01-02.fil',$
;'ukmo_nh_files_02-03.fil',$
;'ukmo_nh_files_03-04.fil',$
;'ukmo_nh_files_04-05.fil']
nyear=n_elements(ifiles)
for iyear=0L,nyear-1L do begin
ifile=''
close,1
openr,1,ifiles(iyear)
ndays=0L
readf,1,ndays
sfile=strarr(ndays)
for iday=0L,ndays-1L do begin
    readf,1,ifile
    print,ifile
    sfile(iday)=ifile
    rd_ukmo,dir+ifile+'_h12.pp.dat',iflag,nc,nr,nlv,alon,alat,wlon,wlat,p,$
            zp,tp,up,vp
    if iflag gt 0L then goto,jumpday
    press=p
    if iday eq 0L then qzt=fltarr(ndays,nlv)
    x2d=fltarr(nc,nr)
    y2d=fltarr(nc,nr)
    for i=0,nc-1 do y2d(i,*)=alat
    for j=0,nr-1 do x2d(*,j)=alon

    for thlev=0L,nlv-1L do begin
        tp1=reform(tp(*,*,thlev))
        index=where(y2d gt 60.)
        qzt(iday,thlev)=min(tp1(index))
    endfor
    jumpday:
endfor		; loop over days
yy=strmid(sfile(0),1,2)
if long(yy) lt 90L then y1='20'+yy
if long(yy) gt 90L then y1='19'+yy
yy=strmid(sfile(ndays-1),1,2)
if long(yy) lt 90L then y2='20'+yy
if long(yy) gt 90L then y2='19'+yy

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='ukmo_tmin_'+y1+'-'+y2+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,ndays,ndays,1,1],[200.,200.,0.3,0.3,200.],/ylog,min_value=0.,$
      xrange=[1,ndays],yrange=[200.,0.3],/nodata,charsize=2,$
      ytitle='Pressure (hPa)',title='MetO Analyses '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,9,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),5,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,200.
    plots,kindex(i)+1,250.,/continue,/data
    xyouts,kindex(i)+1,450.,xlab,/data,alignment=0.5,charsize=3
endfor
nlvls=21
level=170.+5.*findgen(nlvls)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,qzt,1.+findgen(ndays),press,levels=level,/fill,/ylog,$
        /cell_fill,/overplot,c_color=col1,min_value=0.
contour,qzt,1.+findgen(ndays),press,levels=level,c_color=0,/ylog,$
        /follow,/overplot,min_value=0.
index=where(level lt 195.)
contour,qzt,1.+findgen(ndays),press,levels=level(index),c_color=icolmax,/ylog,$
        /follow,/overplot,min_value=0.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
     xtitle='Minimum Vortex Temperature (K)',charsize=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ukmo_tmin_'+y1+'-'+y2+'.ps -rotate -90 ukmo_tmin_'+y1+'-'+y2+'.jpg'
endif
endfor
end
