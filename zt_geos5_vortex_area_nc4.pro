;
; save GEOS-5 vortex area as a function of altitude and time
;
@stddat
@kgmt
@ckday
@kdate

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=9L & lstdy=1L & lstyr=2009L 
ledmn=5L & leddy=1L & ledyr=2010L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)
savefile='vortex_area_'+y1+'_'+y2+'_geos5.sav'

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
cbaryoff=0.06
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
;goto,plotit
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto, saveit

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate

      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      dum1=findfile(dir+sdate+stimes(0)+'nc4')
      if dum1(0) ne '' then ncid=ncdf_open(dir+sdate+stimes(0)+'nc4')
      if dum1(0) eq '' then goto,skipit
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      mark2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,mark2
      ncdf_close,ncid
      ncid=ncdf_open(dir+sdate+stimes(0)+'nc5')
      marknew2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,marknew2
      ncdf_close,ncid

      if kcount eq 0L then begin
         area_zt_nc4=fltarr(kday,nth)
         area_zt_nc5=fltarr(kday,nth)
         harea_zt=fltarr(kday,nth)
         sfile=strarr(kday)
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
      sfile(icount)=lfile
;
; loop over theta
;
      for thlev=0,nth-1 do begin
          mark1=transpose(mark2(*,*,thlev))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=100.*total(area(index))/hem_area
          index=where(lat gt 0. and mark1 lt 0.0,nn)
          if index(0) ne -1 then harea_zt(icount,thlev)=100.*total(area(index))/hem_area
          mark1=transpose(marknew2(*,*,thlev))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1 then area_zt_nc5(icount,thlev)=100.*total(area(index))/hem_area
      endfor

skipit:
icount=icount+1L
goto,jump

saveit:
;
; plot altitude-time series of Arctic vortex area
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
if long(min(yy(index))) lt 90L then y1='20'+string(FORMAT='(I2.2)',min(yy(index)))
if long(min(yy(index))) gt 90L then y1='19'+string(FORMAT='(I2.2)',min(yy(index)))
if long(max(yy(index))) lt 90L then y2='20'+string(FORMAT='(I2.2)',max(yy(index)))
if long(max(yy(index))) gt 90L then y2='19'+string(FORMAT='(I2.2)',max(yy(index)))
;
; save file
;
save,file='vortex_area_'+y1+'_'+y2+'_geos5.sav',area_zt_nc4,area_zt_nc5,harea_zt,th,sfile,y1,y2
index=where(area_zt_nc4 le 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
index=where(area_zt_nc5 le 0.)
if index(0) ne -1L then area_zt_nc5(index)=0./0.

plotit:
restore,savefile
kday=n_elements(SFILE)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_geos5_vortex_area_'+y1+'-'+y2+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,kday,kday,1,1],[400.,400.,4000.,4000.,400.],min_value=0.,$
      xrange=[1,kday],yrange=[400.,4000.],/nodata,charsize=1.5,color=0,$
      ytitle='Theta (K)',xtickname=[' ',' '],xticks=1,title=y1+'-'+y2,charthick=2
kindex=where(strmid(sfile,3,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),0,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,320.
    plots,kindex(i)+1,400.,/continue,/data,color=0
    xyouts,kindex(i)+1,180.,xlab,/data,alignment=0.5,charsize=1.5,color=0,charthick=2
endfor
nlvls=18
col1=1+indgen(nlvls)*icolmax/nlvls
level=2.*findgen(nlvls)
level2=[5.,10.,15.,20.]
;area_zt=smooth(area_zt_nc4,3,/edge_truncate,/nan)
area_zt=area_zt_nc4
harea_zt=smooth(harea_zt,5,/edge_truncate)
index=where(area_zt eq 0.)
if index(0) ne -1L then area_zt(index)=0./0.
index=where(harea_zt eq 0.)
if index(0) ne -1L then harea_zt(index)=0./0.
contour,area_zt,1.+findgen(kday),th,levels=level,/fill,/cell_fill,/overplot,c_color=col1,min_value=0
;contour,area_zt,1.+findgen(kday),th,levels=level,c_color=0,/follow,/overplot,c_labels=0*level,min_value=0
;contour,harea_zt,1.+findgen(kday),th,levels=level2,c_color=mcolor,/follow,/overplot,thick=10,min_value=0
contour,harea_zt,1.+findgen(kday),th,levels=[10.],c_color=mcolor*.75,/follow,/overplot,thick=10,min_value=0
contour,harea_zt,1.+findgen(kday),th,levels=[12.5],c_color=mcolor*.8,/follow,/overplot,thick=10,min_value=0
contour,harea_zt,1.+findgen(kday),th,levels=[15],c_color=mcolor*.85,/follow,/overplot,thick=10,min_value=0
contour,harea_zt,1.+findgen(kday),th,levels=[17.5],c_color=mcolor*.9,/follow,/overplot,thick=10,min_value=0
contour,harea_zt,1.+findgen(kday),th,levels=[20],c_color=mcolor*.95,/follow,/overplot,thick=10,min_value=0
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='GEOS-5 Vortex Area (% NH). Anticyclone Area color contours.',charthick=2,charsize=1.5
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
   spawn,'convert -trim zt_geos5_vortex_area_'+y1+'-'+y2+'.ps -rotate -90 zt_geos5_vortex_area_'+y1+'-'+y2+'.jpg'
endif
end
