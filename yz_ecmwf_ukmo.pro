;
; plot zonal mean temperature for ECMWF, UKMO and differences

@stddat
@kgmt
@ckday
@kdate
@rd_ecmwf
@rd_ukmo

device,decompose=0
dir='/aura5/harvey/ECMWF_data/Datfiles/ecmwf_'
title2=['Potential Vorticity ',$
        'Geopotential Height ',$
        'Temperature ',$
        'Zonal Wind ',$
        'Meridional Wind ',$
        'Vertical Wind ',$
        'Specific Humidity ',$
        'Ozone Mass Mixing Ratio ']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
setplot='x'
read,'setplot ',setplot
lstmn=11
lstdy=1
lstyr=91
ledmn=11
leddy=30
ledyr=91
lsfc=0
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      ECMWF Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 57 then lstyr=lstyr+2000
if ledyr lt 57 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1957 then stop,'Year out of range '
if ledyr lt 1957 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

; define viewport location 
nxdim=750
nydim=750
xorig=[0.1,0.1,0.1]
yorig=[0.75,0.45,0.15]
xlen=0.3
ylen=0.2
cbaryoff=0.05
cbarydel=0.01
!p.thick=1
!p.charsize=1.0
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
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
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
      print,imn,idy,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '

      print,imn,idy,iyr
      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
;
;***Read ECMWF data
      file=dir+smn+'_'+sdy+'_'+syr+'_12Z.dat'
      rd_ecmwf,file,iflg,nc,nr,nl,alon,alat,press,pv,gp,tp,uu,vv,ww,sh,oz
;
;***Read UKMO data
      file='/aura3/data/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nlg,nlat,nlv,ulon,ulat,wlon,wlat,p,g3d,t3d,u3d,v3d

      if iflg ne 0 then goto, jump
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename=$
             string(FORMAT='(a5,i2.2,a1,i2.2,a1,i4)',$
             'ecmwf+ukmo_',imn,'_',idy,'_',iyr)+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.charsize=1.0
      endif

; Declare plotting arrays
      edat=fltarr(nr,nl)
      for k=0,nl-1 do begin
      for j=0,nr-1 do begin
          edat(j,k)=total(tp(*,j,k))/float(nc)
      endfor
      endfor

      udat=fltarr(nlat,nlv)
      for k=0,nlv-1 do begin
      for j=0,nlat-1 do begin
          udat(j,k)=total(t3d(*,j,k))/float(nlg)
      endfor
      endfor
;
; interpolate EC to UKMO pressures (lats are already the same)
;
      lpress=alog(press)
      lp=alog(p)
      edatu=fltarr(nlat,nlv)
      for j=0,nlat-1 do begin
      for k=0,nlv-1 do begin
          for kk=0,nl-2 do begin
              kp1=kk+1
              if lpress(kk) ge lp(k) and lpress(kp1) le lp(k) then begin
              scale=(lpress(kk)-lp(k))/(lpress(kk)-lpress(kp1))
              edatu(j,k)=edat(j,kk)+scale*(edat(j,kp1)-edat(j,kk))
;print,press(kk),p(k),press(kp1),scale
;print,edat(j,kk),edatu(j,k),edat(j,kp1)
              endif
          endfor
      endfor
      endfor
;
; reverse ukmo latitudes
      for k=0,nlv-1 do udat(*,k)=reverse(udat(*,k))
      ulat=reverse(ulat)

; Find data range for autoscaling contours
      erase
      minval=min(edat)
      maxval=max(edat)
      nlvs=20
      cint=(maxval-minval)/float(nlvs)

; Autoscale if scale values for parameter/level are not defined
      level=minval+cint*findgen(nlvs)
      col1=1+indgen(nlvs)*icolmax/nlvs
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
      xyouts,.25,.95,title2(2)+'on '+date,/normal,charsize=2.5

; Mercator projection
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
      contour,edat,alat,press,levels=level,/cell_fill,c_color=col1,/noeras,$
              /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='ECMWF',charsize=1.5,$
              xticks=6
      contour,edat,alat,press,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0,/noeras
      imin=min(level)
      imax=max(level)
      ymnb=yorig(0) -cbaryoff
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

; Mercator projection
      xmn=xorig(1)
      xmx=xorig(1)+xlen
      ymn=yorig(1)
      ymx=yorig(1)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
      contour,udat,ulat,p,levels=level,/cell_fill,c_color=col1,/noeras,$
              /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='MetO',charsize=1.5,$
              xticks=6
      contour,udat,ulat,p,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0
      imin=min(level)
      imax=max(level)
      ymnb=yorig(1) -cbaryoff
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

      xmn=xorig(2)
      xmx=xorig(2)+xlen
      ymn=yorig(2)
      ymx=yorig(2)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
      ddat=edatu
     index=where(edatu gt 0.)
     ddat(index)=edatu(index)-udat(index)
level=-10.+findgen(20)
      contour,ddat,ulat,p,levels=level,/cell_fill,c_color=col1,/noeras,$
              /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='EC-MetO',charsize=1.5,$
              xticks=6
      contour,ddat,ulat,p,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0
      imin=min(level)
      imax=max(level)
      ymnb=yorig(2) -cbaryoff
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

      icount=icount+1L
stop
goto, jump

end
