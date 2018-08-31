;
; latitude-altitude sections of UKMO pressure data
;
@stddat
@kgmt
@ckday
@kdate
@date2uars
@rd_ukmo
@drawvectors

month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
nlg=0l
nlat=0l
nlv=0l
lstmn=1
lstdy=1
lstyr=1992
ledmn=2
leddy=12
ledyr=2010
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      UKMO Version '
;print, ' '
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

rlon=999.0
;read,' Enter desired longitude (0-360 by 3.75, 999 for zonal mean)  ',rlon
ilon=0
ilon=rlon/3.75

; define viewport location 
!NOERAS=-1
setplot='ps'
read,'setplot= ',setplot
loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
icmm1=icolmax-1
icmm2=icolmax-2
!p.background=mcolor
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=750
nydim=750
xorig=[0.125]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01

if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy

;***Read UKMO data
      file='/aura7/harvey/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
              zp,tp,up,vp
      if iflg ne 0 then goto, jump

; Declare plotting arrays
      tp2=fltarr(nlat,nlv)
      u2=fltarr(nlat-1,nlv)
      v2=fltarr(nlat-1,nlv)
      x=fltarr(nlg+1)
      x(0:nlg-1)=alon(0:nlg-1)
      x(nlg)=x(0)+360.

; Zonal means
      if rlon eq 999. then begin      
      for il=0,nlg-1 do begin
          tp2(0:nlat-1,0:nlv-1)=tp2(0:nlat-1,0:nlv-1)+ $
                             tp(il,0:nlat-1,0:nlv-1)
           u2(0:nlat-2,0:nlv-1)=u2(0:nlat-2,0:nlv-1)+ $
                             up(il,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=v2(0:nlat-2,0:nlv-1)+ $
                             vp(il,0:nlat-2,0:nlv-1)
      endfor
      tp2=tp2/nlg
      u2=u2/nlg
      v2=v2/nlg
      endif

; Individual longitudes
      if rlon ne 999. then begin
          tp2(0:nlat-1,0:nlv-1)=tp(ilon,0:nlat-1,0:nlv-1)
           u2(0:nlat-2,0:nlv-1)=up(ilon,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=vp(ilon,0:nlat-2,0:nlv-1)
      endif

if setplot eq 'ps' then begin
   set_plot,'ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_tbar_ubar_meto_'+date+'.ps'
endif

nlvls=27
col1=1+indgen(nlvls)*icolmax/nlvls
tlevel=170.+5.*findgen(nlvls)

      erase
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      tlon=strcompress(string(FORMAT='(F7.2,A5)',rlon,' deg '))
;     mtitle=tlon+' on '+date
      mtitle=date
      !type=2^2+2^3
      contour,tp2,alat,p,yrange=[max(p),min(p)],xrange=[-90,90],/ylog,levels=tlevel,max_val=999.,color=0,$
              /fill,c_color=col1,xtitle='Latitude',ytitle='Pressure (hPa)',xticks=6,title=mtitle
      contour,tp2,alat,p,levels=tlevel,max_val=999.,/follow,/overplot,color=0

contour,u2,wlat,p,levels=30.+10*findgen(20),color=0,/follow,/overplot,c_labels=1+0*findgen(20),/noerase,thick=4
contour,u2,wlat,p,levels=-200.+10*findgen(20),color=mcolor,/follow,/overplot,c_labels=1+0*findgen(20),/noerase,$
        thick=8,c_linestyle=5

; Draw color bar
      imin=min(tlevel)
      imax=max(tlevel)
      ymnb=yorig -cbaryoff
      ymxb=ymnb  +cbarydel
      set_viewport,xmn,xmx,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MetO Zonal Mean Temperature (K)'
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(icmm1)
      for j=1,icmm1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=j
          x1=x1+dx
      endfor

; Close PostScript file and return control to X-windows
      if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yz_tbar_ubar_meto_'+date+'.ps -rotate -90 yz_tbar_ubar_meto_'+date+'.jpg'
       spawn,'/usr/bin/rm yz_tbar_ubar_meto_'+date+'.ps'
    endif

goto, jump

end
