;
; save timeseries of qbar as a function of latitude and altitude
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=800
nydim=800
cbaryoff=0.065
cbarydel=0.02
lstmn=1
lstdy=1
lstyr=1979
ledmn=2
leddy=28
ledyr=2015
lstday=0
ledday=0
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      MERRA Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

restore,'MERRA_qbar_timeseries.sav'	;,latitude_waccm,levels_rad,qbar_timeseries,date_timeseries,doy_timeseries
year=long(strmid(strcompress(date_timeseries,/r),0,4))
index=where(year ne 0L,ndays)
year=year(index)
doy_timeseries=doy_timeseries(index)
qbar_timeseries=qbar_timeseries(index,*,*)

   erase
   !type=2^2+2^3
   set_viewport,.2,.8,.5,.9
   imin=min(year)
   imax=max(year)
   plot,doy_timeseries,qbar_timeseries(*,0,0),xrange=[0,365],yrange=[-60,20],color=0,/noeras,psym=8,ytitle='dtheta/dt (K/day)'		; South Pole at 0.1 hPa
   for ii=0,n_elements(year)-1L do oplot,[doy_timeseries(ii),doy_timeseries(ii)],[qbar_timeseries(ii,nr_waccm-1,0),qbar_timeseries(ii,nr_waccm-1,0)],$
       color=(float(year(ii)-imin)/float(imax-imin))*255.,psym=8						; North Pole at 0.1 hPa

   set_viewport,.2,.8,.1,.4
   index=where(year ne 2014L)
   plot,doy_timeseries(index)-172,qbar_timeseries(index,nr_waccm-1,0),psym=8,xtitle='DFS 2014',color=0,/noeras,xrange=[-110,110],/nodata,ytitle='dtheta/dt (K/day)'
   loadct,0
   oplot,doy_timeseries(index)-172,qbar_timeseries(index,nr_waccm-1,0),psym=8,color=200
   index=where(year eq 2014L)
   oplot,doy_timeseries(index)-172,qbar_timeseries(index,nr_waccm-1,0),thick=5,color=0

end
