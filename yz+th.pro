;
; latitude-altitude sections of UKMO pressure data
; superimpose theta contours
;
@stddat
@kgmt
@ckday
@kdate
@date2uars
@rd_ukmo
@drawvectors

runtitle='UKMO Analyses of '
title2=['Geopotential Height',$
        'Temperature',$
        'Zonal Wind',$
        'Meridional Wind']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
setplot='ps'
colbw='col'
nlg=0l
nlat=0l
nlv=0l
lstmn=0
lstdy=0
lstyr=0
ledmn=0
leddy=0
ledyr=0
lstday=0
ledday=0
uday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      UKMO Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

rlon=0.
ilon=0
read,' Enter desired longitude (0-95, >95 for zonal mean) ',ilon
if ilon gt 95 then rlon=999.

; define viewport location 
nxdim=750
nydim=750
xorig=[0.1]
yorig=[0.15]
xlen=0.8
ylen=0.8
cbaryoff=0.08
cbarydel=0.01

if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz+th.ps'
   if colbw ne 'bw' and colbw ne 'gs' then device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif

; set color table
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
!noeras=1

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

;***Read UKMO data
      file='/aura4/harvey/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr+'_m'+imn+'_d'+idy+'_h12.pp.dat')
      print,file
      rd_ukmo,file,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
              zp,tp,up,vp
      if iflg ne 0 then goto, jump

; Declare plotting arrays
      th2=fltarr(nlat,nlv)
      z2=fltarr(nlat,nlv)
      t2=fltarr(nlat,nlv)
      u2=fltarr(nlat-1,nlv)
      v2=fltarr(nlat-1,nlv)
      x=fltarr(nlg+1)
      x(0:nlg-1)=alon(0:nlg-1)
      x(nlg)=x(0)+360.

; Zonal means
      if rlon eq 999. then begin      
      for il=0,nlg-1 do begin
           z2(0:nlat-1,0:nlv-1)=z2(0:nlat-1,0:nlv-1)+ $
                             zp(il,0:nlat-1,0:nlv-1)/1000.
           t2(0:nlat-1,0:nlv-1)=t2(0:nlat-1,0:nlv-1)+ $
                             tp(il,0:nlat-1,0:nlv-1)
           u2(0:nlat-2,0:nlv-1)=u2(0:nlat-2,0:nlv-1)+ $
                             up(il,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=v2(0:nlat-2,0:nlv-1)+ $
                             vp(il,0:nlat-2,0:nlv-1)

; Convert temperature to potential temperature
           for k=0,nlv-1 do begin
               th2(0:nlat-2,k)=th2(0:nlat-2,k)+$
                  tp(il,0:nlat-2,k)*(1000./p(k))^.286
           endfor
      endfor
      z2=z2/nlg
      t2=t2/nlg
      u2=u2/nlg
      v2=v2/nlg
      th2=th2/nlg
      endif

; Individual longitudes
      if rlon ne 999. then begin
           z2(0:nlat-1,0:nlv-1)=zp(ilon,0:nlat-1,0:nlv-1)/1000.
           t2(0:nlat-1,0:nlv-1)=tp(ilon,0:nlat-1,0:nlv-1)
           u2(0:nlat-2,0:nlv-1)=up(ilon,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=vp(ilon,0:nlat-2,0:nlv-1)

; Convert temperature to potential temperature
      for k=0,nlv-1 do begin
          th2(*,k)=tp(ilon,*,k)*(1000./p(k))^.286
      endfor
      endif

; Find data range for autoscaling contours
      minval=min(u2)
      maxval=max(u2)

; Autoscale if scale values for parameter/level are not defined
      nlvls=15
;     level=minval+((maxval-minval)/nlvls)*findgen(nlvls)
;     level=-75.+10.*findgen(nlvls)
      level=180.+10.*findgen(nlvls)
      col1=1+indgen(nlvls)*icolmax/nlvls

      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx

      if rlon ne 999. then tlon=strcompress(string(FORMAT='(F7.2,A5)',alon(ilon),' deg '))
      if rlon eq 999. then tlon='Zonal Mean '
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
      mtitle=runtitle+title2(2)+' at '+tlon+' on '+date

; UKMO data starts at the NP
      !type=2^2+2^3
      plot,[-90,90,90,-90,-90],[p(0),p(0),p(nlv-1),p(nlv-1),p(0)],/ylog,$
               xrange=[-90,90] ,yrange=[p(0),p(nlv-1)],$
               xtitle='Latitude',ytitle='Pressure (mb)',xticks=6,$
               title=mtitle
      !type=96
      contour,t2,alat,alog10(p),levels=level,xrange=[-90,90],/fill,$
              /cell_fill,/overplot,c_color=col1

      label=1+0*level
      contour,t2,alat,alog10(p),levels=level,c_charsize=1.0,$
              /follow,/overplot,c_linestyle=(level lt 0.0),c_color=0

; superimpose theta contours
      contour,th2,alat,alog10(p),/overplot,thick=3,color=0,$
             levels=[320.,360.,600.,900.,1600.],c_charsize=0.90,c_labels=1+0*level

; Draw color bar
      imin=min(level)
      imax=max(level)
      ymnb=yorig -cbaryoff
      ymxb=ymnb  +cbarydel
      set_viewport,xmn,xmx,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax]
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(icmm1)
      for j=1,icmm1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=j
          x1=x1+dx
      endfor

; Close PostScript file and return control to X-windows
      if setplot eq 'ps' then begin
         device, /close
         set_plot, 'x'
         !p.font=0
         !p.thick=1.0
      endif

      stop
      erase

goto, jump

end
