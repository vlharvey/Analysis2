;
; scatter plot of strongest jet vortex area vs NASH
; color points by theta
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto
@calcelat2d

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=1L & lstdy=1L & lstyr=4L 
ledmn=5L & leddy=1L & ledyr=4L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
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
yorig=[0.15]
xlen=0.75
ylen=0.75
cbaryoff=0.08
cbarydel=0.01
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
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
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
goto,plotit
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
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump

      ncid=ncdf_open(dir+sdate+stimes(0)+'nc5')
      marknew2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,marknew2
      ncdf_close,ncid

      if kcount eq 0L then begin
         area_zt_nc4=fltarr(kday,nth)
         area_zt_nc5=fltarr(kday,nth)
         area_zt_nash=fltarr(kday,nth)
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
          index=where(lat gt 0. and mark1 gt 0.0,nn)
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=100.*total(area(index))/hem_area
          index=where(lat gt 0. and mark1 lt 0.0,nn)
          if index(0) ne -1 then harea_zt(icount,thlev)=100.*total(area(index))/hem_area
          mark1=transpose(marknew2(*,*,thlev))
          index=where(lat gt 0. and mark1 gt 0.0,nn)
          if index(0) ne -1 then area_zt_nc5(icount,thlev)=100.*total(area(index))/hem_area
;
; area of Nash vortex
;
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          speed1=sqrt(u1^2+v1^2)
          pv1=transpose(pv2(*,*,thlev))
          elat1=calcelat2d(pv1,alon,alat)
;
; integrate wind speed and PV in Elat bins
;
          nbins=37
          dy=2.5
          latmin=0.
          latmax=90.
          elatbin=latmin+dy*findgen(nbins)
          speedbin=-999.+0.*fltarr(nbins)                               ; average windspeed per elat bin
          pvbin=0.*elatbin
          for n=0,nbins-2 do begin
              t=where(pv1 lt 1000. and lat ge latmin and elat1 ge elatbin(n) and elat1 lt elatbin(n+1),it)
              if it gt 2 then begin
                  result=moment(pv1(t))
                  pvbin(n)=result(0)
                  if min(lat(t))-latmin le dy then begin ; make sure bins are resolved (do not intersect latmin)
                     speedbin(n)=-999.
                     goto,jumpnhbin
                  endif
                  speedbin(n)=total(speed1(t))/float(it)
              endif
              jumpnhbin:
          endfor                                                        ; loop over Elat bins
          s=where(lat ge latmin and elat1 ge elatbin(nbins-1),is)
          if is gt 2 then begin
             result=moment(pv1(s))
             pvbin(n)=result(0)
             if min(lat(s))-latmin gt dy then speedbin(nbins-1)=total(speed1(s))/float(is)
          endif
;
; compute PV gradient wrt Equivalent latitude
;
          dpvbin=0.*pvbin
          for i=0,nbins-2L do dpvbin(i)=pvbin(i+1)-pvbin(i)
          dpvbin(nbins-1)=pvbin(nbins-1)-pvbin(nbins-2)
;
; impose Nash filter poleward of 80deg (and add new one Equatorward of lat0)
;
          lat0=70.
          index=where(elatbin ge lat0)                                  ; filter down poleward of 80deg
          speedbin(index)=speedbin(index)*(90.-elatbin(index))/30.
          dpvbin(index)=dpvbin(index)*(90.-elatbin(index))/30.
          lat0=25.
          if th(thlev) lt 600. then lat0=45.
          index=where(elatbin le lat0)                                  ; filter down equatorward of lat0
          speedbin(index)=speedbin(index)*(elatbin(index))/(2.*lat0)
          dpvbin(index)=dpvbin(index)*(elatbin(index))/(2.*lat0)
          dpvbin=dpvbin/max(dpvbin)                                     ; normalise
;
; vortex edge is where dPV/dElat multiplied by the wind speed integrated in Elat bins is maximum
; and integrated wind speed must be greater than 15.2 m/s
;
          prod=dpvbin*speedbin
          index=where(prod eq max(prod))
          if index(0) ne -1L then edgepv=pvbin(index)
;print,th(thlev),edgepv(0),speedbin(index(0))
          if speedbin(index(0)) lt 15.2 then goto,skiplev

          index=where(lat gt 0. and pv1 ge edgepv(0),nn)
          if index(0) ne -1 then area_zt_nash(icount,thlev)=100.*total(area(index))/hem_area
skiplev:
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
y1='20'+string(format='(i2.2)',long(min(yy(index))))
y2='20'+string(format='(i2.2)',long(max(yy(index))))
;
; save file
;
save,file='vortex_area_'+y1+'_'+y2+'_geos5_vs_nash.sav',area_zt_nc4,area_zt_nc5,area_zt_nash,$
     harea_zt,th,sfile,y1,y2

plotit:
restore,'vortex_area_2006_2007_geos5_vs_nash.sav'
kday=n_elements(SFILE)
nth=n_elements(th)
yy=strmid(sfile,6,2)
index=where(yy ne '')
y1='20'+string(format='(i2.2)',long(min(yy(index))))
y2='20'+string(format='(i2.2)',long(max(yy(index))))

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='scatter_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
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
!type=2^2+2^3
imin=500.
imax=5000.
th2d=0.*area_zt_nash
for i=0,kday-1 do th2d(i,*)=th
plot,AREA_ZT_NC4,AREA_ZT_NASH,xrange=[0.,60.],yrange=[0.,60.],title='DJF '+y1+'-'+y2,$
     xtitle='Strongest Jet Vortex Area (% of Hemisphere)',charsize=1.5,$
     ytitle='Nash Vortex Area',color=0,/nodata
mm2d=strarr(kday,nth)
for k=0L,nth-1L do mm2d(*,k)=strmid(sfile,0,2)
index=where(area_zt_nash ne 0. and area_zt_nc4 ne 0. and th2d gt imin and th2d lt imax and $
;          mm2d eq '01')
           (mm2d eq '12' or mm2d eq '01' or mm2d eq '02'))
oplot,findgen(100),findgen(100),color=0
for i=0,n_elements(index)-1L do oplot,[AREA_ZT_NC4(index(i)),AREA_ZT_NC4(index(i))],$
   [AREA_ZT_NASH(index(i)),AREA_ZT_NASH(index(i))],color=((th2d(index(i))-400.)/(5000.-400.))*mcolor
for i=0,n_elements(index)-1L do oplot,[AREA_ZT_NC4(index(i)),AREA_ZT_NC4(index(i))],$
    [AREA_ZT_NASH(index(i)),AREA_ZT_NASH(index(i))],color=((th2d(index(i))-imin)/(imax-imin))*mcolor,psym=8,symsize=2

ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Theta (K)',charsize=1.5
ybox=[0,10,10,0,0]
nlvls=20L
col1=1+indgen(nlvls)*mcolor/nlvls
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
   spawn,'convert -trim scatter_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.ps -rotate -90 '+$
         'scatter_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.jpg'
endif
end
