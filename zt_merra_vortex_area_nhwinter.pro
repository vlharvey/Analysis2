;
; save vortex area as a function of altitude and time
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=9L & lstdy=1L & lstyr=2009L 
ledmn=5L & leddy=1L & ledyr=2010L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)

lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
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
xorig=[0.15]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
goto,plotit

dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
resultf=file_search(dir+'*.nc3')
ndays=n_elements(resultf)
icount=0
for iday =0L, ndays-1L do begin
    result=strsplit(resultf(iday),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(0),'_',/extract)
    sdate=result3(2)
;
; read daily file
;
        ncfile0=resultf(iday)
        print,ncfile0
        ncid=ncdf_open(ncfile0)
ncdf_diminq,ncid,0,name,nr
ncdf_diminq,ncid,1,name,nc
ncdf_diminq,ncid,2,name,nl
alon=fltarr(nc)
alat=fltarr(nr)
thlev=fltarr(nl)
mark=fltarr(nr,nc,nl)
ncdf_varget,ncid,0,alon
ncdf_varget,ncid,1,alat
ncdf_varget,ncid,2,th
ncdf_varget,ncid,8,mark2
ncdf_close,ncid

      if iday eq 0L then begin
         nth=nl
         area_zt_nc4=fltarr(ndays,nth)
         sfile=strarr(ndays)
         dum=transpose(mark2(*,*,0))
         lon=0.*dum
         lat=0.*dum
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy    ; area of each grid point
         endfor
         kcount=1L
      endif
      sfile(icount)=sdate
;
; loop over theta
;
      for thlev=0,nth-1 do begin
          mark1=transpose(mark2(*,*,thlev))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=100.*total(area(index))/hem_area
      endfor
      jumpday:
      icount=icount+1L
endfor          ; loop over days

saveit:
;
; plot altitude-time series of Arctic vortex area
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
y1='20'+string(FORMAT='(I2.2)',min(yy(index)))
y2='20'+string(FORMAT='(I2.2)',max(yy(index)))
;
; save file
;
save,file='/Users/harvey/Harvey_etal_2014/Post_process/vortex_area_nh_merra.sav',area_zt_nc4,th,sfile,y1,y2
index=where(area_zt_nc4 le 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.

plotit:
restore,'/Users/harvey/Harvey_etal_2014/Post_process/vortex_area_nh_merra.sav'
kday=n_elements(SFILE)
;
; plot zt vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=18
col1=1+indgen(nlvls)*icolmax/nlvls
level=2.*findgen(nlvls)
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
print,th
rlev=2000.
;read,'Enter theta ',rlev
index=where(rlev eq th)
ilev=index(0)
slev=strcompress(long(rlev),/remove_all)
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_vortex_area_nh_'+slev+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

plot,findgen(366),1.+findgen(366),color=0,xtitle='DOY',thick=6,yrange=[0.,40.],/noeras,ytitle='Vortex Area (% of NH)',title='MERRA '+slev,/nodata
area2d=fltarr(35,366)
for iyear=1979,2013 do begin
    iyr=long(strmid(sfile,0,4))
    index=where(iyr eq iyear,npts)
    areaplot=reform(area_zt_nc4(index,ilev))
    area2d(iyear-1979,0:npts-1)=areaplot
    oplot,1+findgen(npts),areaplot,color=((float(iyear)-1979.)/(2014.-1979.))*mcolor,thick=5
endfor
areamedian=median(area2d,dim=1)
npts=n_elements(areamedian)
oplot,1+findgen(npts),areamedian,color=0,thick=10
;
; color bar
;
imin=1979
imax=2013
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charthick=2,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ../Figures/merra_vortex_area_nh_'+slev+'.ps -rotate -90 ../Figures/merra_vortex_area_nh_'+slev+'.jpg'
endif

;endfor	; loop 
end
