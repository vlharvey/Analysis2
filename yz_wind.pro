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
lstyr=2002
ledmn=1
leddy=1
ledyr=2002
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
   device,/landscape,bits=8,filename='yz.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif

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

;***Read UKMO data
      file='/aura7/harvey/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
              zp,tp,up,vp
      if iflg ne 0 then goto, jump

; Declare plotting arrays
      dat=fltarr(nlat,nlv)
      u2=fltarr(nlat-1,nlv)
      v2=fltarr(nlat-1,nlv)
      x=fltarr(nlg+1)
      x(0:nlg-1)=alon(0:nlg-1)
      x(nlg)=x(0)+360.

; Zonal means
      if rlon eq 999. then begin      
      for il=0,nlg-1 do begin
          dat(0:nlat-1,0:nlv-1)=dat(0:nlat-1,0:nlv-1)+ $
                              zp(il,0:nlat-1,0:nlv-1)/1000.
           u2(0:nlat-2,0:nlv-1)=u2(0:nlat-2,0:nlv-1)+ $
                             up(il,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=v2(0:nlat-2,0:nlv-1)+ $
                             vp(il,0:nlat-2,0:nlv-1)
      endfor
      dat=dat/nlg
      u2=u2/nlg
      v2=v2/nlg
      endif

; Individual longitudes
      if rlon ne 999. then begin
          dat(0:nlat-1,0:nlv-1)=zp(ilon,0:nlat-1,0:nlv-1)/1000.
           u2(0:nlat-2,0:nlv-1)=up(ilon,0:nlat-2,0:nlv-1)
           v2(0:nlat-2,0:nlv-1)=vp(ilon,0:nlat-2,0:nlv-1)
      endif

; Convert temperature to potential temperature
;     for k=0,nlv-1 do begin
;         dat(*,k)=dat(*,k)*(1000./p(k))^.286
;     endfor

; Find data range for autoscaling contours
;     minval=min(dat)
;     maxval=max(dat)
      minval=min(u2)
      maxval=max(u2)
print,file,max(u2)
      minval=-100.
      maxval=100.

; Autoscale if scale values for parameter/level are not defined
      nlvls=21
      level=minval+10.*findgen(nlvls)
      col1=1+indgen(nlvls)*icolmax/nlvls
      erase
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx

      tlon=strcompress(string(FORMAT='(F7.2,A5)',rlon,' deg '))
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
;      mtitle=tlon+' on '+date
      mtitle=date

; UKMO data starts at the NP
;     plot_io,[alat(nlat-1),alat(0),alat(0),alat(nlat-1),alat(nlat-1)],$
;             [p(0),p(0),p(nlv-1),p(nlv-1),p(0)],$
;              xrange=[alat(nlat-1),alat(0)],yrange=[p(0),p(nlv-1)],$
;              xtitle='!6Latitude',ytitle='!6Pressure (mb)',xticks=6,$
;              title=mtitle
;     !type=96
;     contour,dat,alat,alog10(p),levels=level,max_val=999.,$
;             xrange=[alat(nlat-1),alat(0)],yrange=[3.,1.],/fill,$
;             /cell_fill,/overplot,c_color=col1,xstyle=4,ystyle=4
;     label=1+0*level
;     contour,dat,alat,alog10(p),levels=level,c_charsize=1.0,$
;             c_labels=1+0*level,max_val=999.,/follow,$
;             xrange=[alat(nlat-1),alat(0)],yrange=[3.,1.],$
;             /overplot,c_linestyle=(level lt 0.0),$
;             xstyle=4,ystyle=4,c_color=0
;     plot_io,[wlat(nlat-2),wlat(0),wlat(0),wlat(nlat-2),wlat(nlat-2)],$
;             [p(0),p(0),p(nlv-1),p(nlv-1),p(0)],$
      !type=2^2+2^3
      plot,[-90,90,90,-90,-90],[p(0),p(0),p(nlv-1),p(nlv-1),p(0)],/ylog,$
               xrange=[-90,90] ,yrange=[p(0),p(nlv-1)],color=0,$
               xtitle='Latitude',ytitle='Pressure (mb)',xticks=6,$
               title=mtitle
      !type=96
      contour,u2,wlat,alog10(p),levels=level,max_val=999.,$
              xrange=[-90,90],$
              yrange=[alog10(p(0)),alog10(p(nlv-1))],/fill,$
              /cell_fill,/overplot,c_color=col1

      label=1+0*level
      contour,u2,wlat,alog10(p),levels=level,c_charsize=1.0,$
              c_labels=1+0*level,max_val=999.,/follow,$
              xrange=[wlat(nlat-2),wlat(0)],$
              yrange=[alog10(p(0)),alog10(p(nlv-1))],$
              /overplot,c_linestyle=(level lt 0.0),c_color=0

; Draw color bar
      imin=min(level)
      imax=max(level)
      ymnb=yorig -cbaryoff
      ymxb=ymnb  +cbarydel
      set_viewport,xmn,xmx,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MetO Zonal Mean Zonal Wind (m/s)'
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
         device, /close
         set_plot, 'x'
         !p.font=0
         !p.thick=1.0
      endif

goto, jump

end
