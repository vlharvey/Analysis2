;
; plot timeseries of MERRA anticyclone area at a given altitude
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.20]
yorig=[0.3]
cbaryoff=0.1
cbarydel=0.01
xlen=0.6
ylen=0.4
device,decompose=0
!NOERAS=-1
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
;syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
dum=1979+indgen(36)
syear=strcompress(dum,/remove_all)
nyear=n_elements(syear)
nlvls=nyear
col1=1+(indgen(nlvls)/float(nlvls))*mcolor

smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; get file listing
;
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'

for iyear=0L,nyear-1L do begin

slat='80'
restore,'merra_ubar+area_'+syear(iyear)+'_'+slat+'.sav'	;,ubar,dayno,pressure,altitude,area1,area1h,polarh2o,polarpv
area1=smooth(area1,7,/edge_truncate)
area1h=smooth(area1h,7,/edge_truncate)
print,iyear,' ',syear(iyear)
if iyear eq 0L then begin
print,reform(ALTitudE(0,*))
ralt=27.
read,'Enter desired altitude ',ralt
height=reform(ALTitudE(0,*))
index=where(abs(height-ralt) eq min(abs(height-ralt)))
ialt=index(0)
salt=strcompress(height(index(0)),/remove_all)+'km'
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='timeseries_merra_anticarea_'+'_'+salt+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=0.
imax=40.
ubarlev=reform(area1h(*,ialt))
plot,dayno,ubarlev,color=col1(iyear),ytitle='MERRA Anticyclone Area at '+slat+' N',yrange=[imin,imax],charsize=1.5,thick=1,/nodata	;,xrange=[60.,180]
ubarlev_all=fltarr(122,nyear)
endif
ubarlev=reform(area1h(*,ialt))
ubarlev_all(*,iyear)=ubarlev
;oplot,dayno,ubarlev,color=col1(iyear),thick=3
if syear(iyear) eq '1994' then begin
;   oplot,dayno,ubarlev,color=col1(iyear),thick=10
   ubarlev94=ubarlev
endif
if syear(iyear) eq '1997' then begin
;   oplot,dayno,ubarlev,color=col1(iyear),thick=10
   ubarlev97=ubarlev
endif
if syear(iyear) eq '2008' then begin
;   oplot,dayno,ubarlev,color=col1(iyear),thick=10
   ubarlev08=ubarlev
endif
if syear(iyear) eq '2011' then begin
;   oplot,dayno,ubarlev,color=col1(iyear),thick=10
   ubarlev11=ubarlev
endif
if syear(iyear) eq '2013' then begin
;   oplot,dayno,ubarlev,color=col1(iyear),thick=10
   ubarlev13=ubarlev
endif

endfor
loadct,0
ubarmin=0.*ubarlev
ubarmax=0.*ubarlev
for i=60L,181L do begin
    oplot,[i,i],[min(ubarlev_all(i-60,*)),max(ubarlev_all(i-60,*))],thick=5,color=200
;print,i,min(ubarlev_all(i-60,*)),max(ubarlev_all(i-60,*))
    ubarmin(i-60)=min(ubarlev_all(i-60,*))
    ubarmax(i-60)=max(ubarlev_all(i-60,*))
endfor
oplot,dayno,ubarmin,color=0,thick=3
oplot,dayno,ubarmax,color=0,thick=3
loadct,39
dum=mean(ubarlev_all,dim=2)
oplot,dayno,dum,color=0,thick=10
oplot,dayno,ubarlev94,color=col1(15),thick=5
oplot,dayno,ubarlev97,color=col1(18),thick=5
oplot,dayno,ubarlev08,color=col1(10),thick=5
oplot,dayno,ubarlev11,color=col1(29),thick=5
oplot,dayno,ubarlev13,color=col1(34),thick=5
xyouts,130.,imax-5,salt,charsize=2,color=0,/data,charthick=2
xyouts,70.,imax-5,'1994',charsize=2,color=col1(15),/data,charthick=2
xyouts,70.,imax-9,'1997',charsize=2,color=col1(18),/data,charthick=2
xyouts,70.,imax-13,'2008',charsize=2,color=col1(10),/data,charthick=2
xyouts,70.,imax-17,'2011',charsize=2,color=col1(29),/data,charthick=2
xyouts,70.,imax-21,'2013',charsize=2,color=col1(34),/data,charthick=2
plots,134,0.
plots,134,-3,/continue,color=col1(34),thick=3
plots,145,0.
plots,145,-3,/continue,color=col1(10),thick=3
plots,140,0.
plots,140,-3,/continue,color=col1(29),thick=3

;imin=long(min(syear))
;imax=long(max(syear))
;ymnb=yorig(0) -cbaryoff
;ymxb=ymnb  +cbarydel
;set_viewport,xmn,xmx,ymnb,ymxb
;!type=2^2+2^3+2^6
;plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,/noeras
;ybox=[0,10,10,0,0]
;x1=imin
;dx=(imax-imin)/float(nlvls)
;for j=0,nlvls-1 do begin
;xbox=[x1,x1,x1+dx,x1+dx,x1]
;polyfill,xbox,ybox,color=col1(j)
;x1=x1+dx
;endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim timeseries_merra_anticarea_'+salt+'.ps -rotate -90 '+$
                            'timeseries_merra_anticarea_'+salt+'.jpg'
     endif

end
