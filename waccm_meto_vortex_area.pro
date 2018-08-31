;
; reads in WACCM .sav and MetO and calculates vortex area
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=700
nydim=700
xorig=[0.1]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.03
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
lstmn=11
lstdy=1
lstyr=2003
ledmn=5
leddy=1
ledyr=2004
lstday=0
ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dir='/aura3/data/WACCM_data/Datfiles/PV_ELAT_S1_'
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd

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
      if ndays gt ledday then goto,plotit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=strcompress(string(iyr),/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      date=syr+smn+sdy
      ifile=date+'.sav'
      restore,dir+ifile
      pv2=pv
      elat2=elat
      if icount eq 0 then begin
         rpress=2.0793250
;        print,lev
;        read,'Enter pressure level ',rpress
         index=where(long(rpress*10000.) eq long(lev*10000.))
         if index(0) eq -1 then stop,'Invalid pressure level '
         ipress=index(0)
         spress=strcompress(string(rpress),/remove_all)
         waccm_area=fltarr(ledday-lstday+1L)
         ukmo_area=fltarr(ledday-lstday+1L)
         dates=strarr(ledday-lstday+1L)
      endif
      dates(icount)=date
      elat1=reform(elat2(*,*,ipress))
      pv1=reform(pv2(*,*,ipress))
      nc=n_elements(alon)
      nr=n_elements(alat)
      pv1(*,nr-1)=pv1(*,nr-2)
      elat1(*,nr-1)=elat1(*,nr-2)
      pv1=smooth(pv1,3,/edge_truncate)
      elat1=smooth(elat1,3,/edge_truncate)
      lon=0.*pv1
      lat=0.*pv1
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
      index=where(lat ge 40. and elat1 ge 60.)
      if index(0) ne -1 then waccm_area(icount)=total(area(index))/hem_area
      print,ifile

      ifile=mon(imn-1)+sdy+'_'+uyr
      dum1=findfile(diru+ifile+'.nc3')
      if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      pv2=fltarr(nr,nc,nth)
      marksf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      if icount eq 0L then begin
         rtheta=0.
         print,th
         read,'Enter theta level ',rtheta
         index=where(long(rtheta) eq long(th))
         if index(0) eq -1 then stop,'Invalid theta level '
         thlev=index(0)
         stheta=strcompress(string(rtheta),/remove_all)
      endif
      ncdf_varget,ncid,3,pv2
      ncdf_varget,ncid,10,marksf2
      ncdf_close,ncid
      pv1=transpose(pv2(*,*,thlev))
      elat1=calcelat2d(pv1,alon,alat)
      pv1(*,nr-1)=pv1(*,nr-2)
      elat1(*,nr-1)=elat1(*,nr-2)
      pv1=smooth(pv1,3,/edge_truncate)
      elat1=smooth(elat1,3,/edge_truncate)
      mark1=transpose(marksf2(*,*,thlev))
      lon=0.*mark1
      lat=0.*mark1
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

      index=where(lat gt 0. and mark1 eq 1.0)
;     index=where(lat gt 0. and elat1 gt 60.)
      if index(0) ne -1 then ukmo_area(icount)=total(area(index))/hem_area
;print,'MetO ',total(area(index))/hem_area
;erase
;contour,elat1,alon,alat,nlevels=20
;contour,elat1,alon,alat,level=[60.],/overplot,thick=3
;oplot,lon(index),lat(index),psym=3
;stop
         
      icount=icount+1L
goto,jump

plotit:
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='waccm_meto_vortex_area_'+spress+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
xindex=where(strmid(dates,6,2) eq '01' or strmid(dates,6,2) eq '15',nxticks)
nday=n_elements(waccm_area)
plot,findgen(nday),ukmo_area,thick=3,xticks=nxticks-1,xtickv=xindex,$
     xtickname=strmid(dates(xindex),4,4),xrange=[0,nday-1],ytitle='Hemispheric Fraction',$
     title='Arctic Vortex Area',yrange=[0.0,.4]
oplot,findgen(nday),waccm_area,thick=3,color=mcolor*.9
;for i=0L,nxticks-1L do xyouts,xindex(i),ymn-0.2
if setplot eq 'ps' then device, /close
end
