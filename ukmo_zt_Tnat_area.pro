;
; plot the area enclosed by Tnat isopleth in 2004-2005 from
; Tnat values at theta levels and on dates in Datfiles/Tnat_SC.dat
;
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
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
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
;
; restore C. Singleton Tnat save file
;
restore,'/aura2/harvey/Singleton/Datfiles/Tnat_SC.dat
caldat,jday,m,d,y	; convert Julian day to YYYY,MM,DD
;
; loop over days
;
kday=n_elements(y)
for iday=0L,kday-1L do begin
    idy=d(iday)
    imn=m(iday)
    iyr=y(iday)
    if iyr ge 2000 then iyr1=iyr-2000
    if iyr lt 2000 then iyr1=iyr-1900
    uyr=string(FORMAT='(I2.2)',iyr1)
    syr=string(FORMAT='(I4.4)',iyr)
    smn=string(FORMAT='(I2.2)',imn)
    sdy=string(FORMAT='(I2.2)',idy)

    ifile=mon(imn-1)+sdy+'_'+uyr
    lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
    ncid=ncdf_open(diru+ifile+'.nc3')
    print,ifile
    ncdf_diminq,ncid,0,name,nr
    ncdf_diminq,ncid,1,name,nc
    ncdf_diminq,ncid,2,name,nth
    alon=fltarr(nc)
    alat=fltarr(nr)
    th=fltarr(nth)
    p2=fltarr(nr,nc,nth)
    ncdf_varget,ncid,0,alon
    ncdf_varget,ncid,1,alat
    ncdf_varget,ncid,2,th
    ncdf_varget,ncid,4,p2
    ncdf_close,ncid
;
; MetO temperature
;
    t2=0.*p2
    for k=0,nth-1 do t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
;
; do all of this only on the first day
;
    if iday eq 0L then begin
;
; interpolate AVGTNAT(124, 42) from zo levels to MetO theta levels (th)
;
       avgtnat_meto=fltarr(kday,nth)
       for i=0L,kday-1L do begin
       for kk=0L,nth-1L do begin
           th0=th(kk)
           for k=1L,n_elements(zo)-1L do begin
               thup=zo(k-1) & thlw=zo(k)
               if thup ge th0 and thlw le th0 and $
                  AVGTNAT(i,k-1) ne -99. and AVGTNAT(i,k) ne -99. then begin
                  zscale=(thup-th0)/(thup-thlw)
                  avgtnat_meto(i,kk)=AVGTNAT(i,k-1)+zscale*(AVGTNAT(i,k)-AVGTNAT(i,k-1))
               endif
           endfor
       endfor
       endfor
;
; calculate area 
;
       area_ave=fltarr(kday,nth)
       sfile=strarr(kday)
       dum=transpose(t2(*,*,0))
       lon=0.*dum
       lat=0.*dum
       for i=0,nc-1 do lat(i,*)=alat
       for j=0,nr-1 do lon(*,j)=alon
       area=0.*lat
       nrr=91
       yeq=findgen(nrr)
       latcircle=fltarr(nrr)
       latsum=fltarr(nrr)
       hem_frac=fltarr(nrr)
       for j=0,nrr-2 do begin
           hy=re*dtr
           dx=re*cos(yeq(j)*dtr)*360.*dtr
           latcircle(j)=dx*hy	; area in each latitude circle
       endfor
       for j=0L,nrr-1 do latsum(j)=total(latcircle(j:nrr-1))
       for j=0,nrr-1 do begin
           index=where(yeq ge yeq(j))
; fraction of the hemisphere in each latitude circle
           if index(0) ne -1 then $
              hem_frac(j)=100.*total(latcircle(index))/hem_area
           if yeq(j) eq 0. then hem_frac(j)=100.
       endfor
       deltax=alon(1)-alon(0)
       deltay=alat(1)-alat(0)
       for j=0,nr-1 do begin
           hy=re*deltay*dtr
           dx=re*cos(alat(j)*dtr)*deltax*dtr
           area(*,j)=dx*hy	; area of each grid point
       endfor
    endif	; if first day
    sfile(iday)=lfile
;
; sum area of gridpoints within Tnat isopleth
;
    myr=iyr-1991L
    for thlev=0,nth-1 do begin
        temp=transpose(t2(*,*,thlev))
        index=where(lat gt 0. and temp le avgtnat_meto(iday,thlev))
        if index(0) ne -1 then begin
           a0=total(area(index))
           area_ave(iday,thlev)=a0/1.e6		; millions of sqare km
        endif
    endfor

endfor		; loop over days
;
; plot altitude-time series of Arctic vortex area
;
plotit:
yy=strmid(sfile(0),6,2)
if long(yy) lt 90L then y1='20'+yy
if long(yy) gt 90L then y1='19'+yy
yy=strmid(sfile(kday-1),6,2)
if long(yy) lt 90L then y2='20'+yy
if long(yy) gt 90L then y2='19'+yy

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='ukmo_zt_tnat_area_'+y1+'-'+y2+'.ps'
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
plot,[1,kday,kday,1,1],[300.,300.,700.,700.,300.],min_value=0.,$
      xrange=[1,kday],yrange=[300.,700.],/nodata,charsize=2,$
      ytitle='Theta (K)',title='MetO Analyses '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,3,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),0,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,270.
    plots,kindex(i)+1,300.,/continue,/data
    xyouts,kindex(i)+1,250.,xlab,/data,alignment=0.5,charsize=3
endfor
nlvls=21
level=1.+4.*findgen(nlvls)
col1=1+indgen(nlvls)*icolmax/nlvls
area_ave=smooth(area_ave,3,/edge_truncate)
index=where(abs(area_ave) lt 1.)
if index(0) ne -1 then area_ave(index)=-9999.
contour,area_ave,1.+findgen(kday),th,levels=level,/fill,$
        /cell_fill,/overplot,c_color=col1,min_value=0.
contour,area_ave,1.+findgen(kday),th,levels=level,c_color=0,$
        /follow,/overplot,min_value=0.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
     xtitle='N.H. Area where T<Tnat (millions of sq. km)',charsize=2
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
   spawn,'convert -trim ukmo_zt_tnat_area_'+y1+'-'+y2+'.ps -rotate -90 ukmo_zt_tnat_area_'+y1+'-'+y2+'.jpg'
endif

end
