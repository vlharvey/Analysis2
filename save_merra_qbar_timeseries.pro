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

dir='/Volumes/Data/MERRA_data/Datfiles/'

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      sdate=string(FORMAT='(i4.4,i2.2,i2.2)',iyr,imn,idy)
      print,sdate
      result=file_search(dir+'MERRA-on-WACCM_rad_'+sdate+'.sav')
      if result(0) eq '' then print,'missing data on '+sdate
      if result(0) eq '' then goto,jump
      restore,dir+'MERRA-on-WACCM_rad_'+sdate+'.sav'	;,longitude_waccm,latitude_waccm,levels_rad,qnew
;
; first day
;
      if icount eq 0L then begin
         nr_waccm=n_elements(latitude_waccm)
         nl_rad=n_elements(levels_rad)
         qbar_timeseries=fltarr(kday,nr_waccm,nl_rad)
         date_timeseries=lonarr(kday)
         doy_timeseries=lonarr(kday)
      endif
      qbar_timeseries(icount,*,*)=mean(qnew,dim=1)
      date_timeseries(icount)=long(sdate)
      doy_timeseries(icount)=iday
;
; test
;
;  erase
;  !type=2^2+2^3
;  set_viewport,.2,.8,.55,.9
;  qbar_new=mean(qnew,dim=1)
;  contour,qbar_new,latitude_waccm,levels_rad,/ylog,xrange=[-90,90],yrange=[1000.,0.01],levels=-20.+findgen(20),thick=3,$
;          ytitle='Pressure (hPa)',xtitle='Latitude',xticks=6,c_linestyle=5,color=0,title=sdate
;  contour,qbar_new,latitude_waccm,levels_rad,/overplot,levels=1.+findgen(19),thick=3,color=250

;  if icount gt 0L then begin
;  set_viewport,.2,.8,.1,.45
;  plot,doy_timeseries(0:icount),qbar_timeseries(0:icount,0,0),xrange=[0,365],yrange=[-50,30],color=0,/noeras,psym=8		; South Pole at 0.1 hPa
;  oplot,doy_timeseries(0:icount),qbar_timeseries(0:icount,nr_waccm-1,0),color=250,psym=8						; North Pole at 0.1 hPa
;  endif
; wait,.1

   icount=icount+1L
goto, jump
saveit:
;
; IDL save file of qbar timeseries
;
save,file='MERRA_qbar_timeseries.sav',latitude_waccm,levels_rad,qbar_timeseries,date_timeseries,doy_timeseries
end
