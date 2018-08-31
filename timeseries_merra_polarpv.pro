;
; plot timeseries of MERRA polar PV at a given altitude
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
ralt=40.
;read,'Enter desired altitude ',ralt
height=reform(ALTitudE(0,*))
index=where(abs(height-ralt) eq min(abs(height-ralt)))
ialt=index(0)
salt=strcompress(height(index(0)),/remove_all)+'km'
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='timeseries_merra_polarh3o_'+'_'+salt+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

ubarlev_all=fltarr(122,nyear)
endif
ubarlev=reform(polarpv(*,ialt))*1.e6
ubarlev_all(*,iyear)=ubarlev

if syear(iyear) eq '1994' then begin
   ubarlev94=ubarlev
endif
if syear(iyear) eq '1997' then begin
   ubarlev97=ubarlev
endif
if syear(iyear) eq '2005' then ubarlev05=ubarlev
if syear(iyear) eq '2006' then ubarlev06=ubarlev
if syear(iyear) eq '2007' then ubarlev07=ubarlev
if syear(iyear) eq '2008' then ubarlev08=ubarlev
if syear(iyear) eq '2009' then ubarlev09=ubarlev
if syear(iyear) eq '2010' then ubarlev10=ubarlev
if syear(iyear) eq '2011' then ubarlev11=ubarlev
if syear(iyear) eq '2012' then ubarlev12=ubarlev
if syear(iyear) eq '2013' then ubarlev13=ubarlev
if syear(iyear) eq '2014' then ubarlev14=ubarlev

endfor

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(ubarlev_all ne 0.)
imin=min(ubarlev_all(index))
imax=max(ubarlev_all(index))
plot,dayno,ubarlev,color=0,ytitle='Mean PV in Anticyclone',yrange=[imin,imax],charsize=1.5,thick=1,/nodata    ;,xrange=[60.,180]

loadct,0
dum=fltarr(122)
ubarmin=0.*ubarlev
ubarmax=0.*ubarlev
for ii=0L,121L do begin
    index=where(ubarlev_all(ii,*) ne 0.)
    if index(0) ne -1L then begin
       dum(ii)=mean(ubarlev_all(ii,index))
       ubarmin(ii)=min(ubarlev_all(ii,index))
       ubarmax(ii)=max(ubarlev_all(ii,index))
       oplot,[dayno(ii),dayno(ii)],[ubarmin(ii),ubarmax(ii)],thick=5,color=200
    endif
print,ii+60
endfor
;dum=mean(ubarlev_all,dim=2)
;oplot,dayno,ubarmin,color=0,thick=3
;oplot,dayno,ubarmax,color=0,thick=3
loadct,39
oplot,dayno,dum,color=0,thick=10
;
; put individual years atop
;
for iyear=0L,nyear-1L do begin
    ubarlev=ubarlev_all(*,iyear)
    index=where(ubarlev eq 0.)
    if index(0) ne -1L then ubarlev(index)=0./0.
    if syear(iyear) eq '1982' or syear(iyear) eq '1994' or syear(iyear) eq '1997' or syear(iyear) eq '2000' or syear(iyear) eq '2002' or syear(iyear) eq '2005' or syear(iyear) eq '2007' or syear(iyear) eq '2011' or syear(iyear) eq '2013' then begin
       oplot,dayno,smooth(ubarlev,3,/Nan,/edge_truncate),color=col1(iyear),thick=4,linestyle=5
    endif else begin oplot,dayno,smooth(ubarlev,3,/Nan,/edge_truncate),color=col1(iyear),thick=1
    endelse
endfor

xyouts,70.,imax-0.2,salt,charsize=2,color=0,/data,charthick=2

imin=long(min(syear))
imax=long(max(syear))
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,/noeras
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim timeseries_merra_polarpv_'+salt+'.ps -rotate -90 '+$
                            'timeseries_merra_polarpv_'+salt+'.jpg'
     endif

end
