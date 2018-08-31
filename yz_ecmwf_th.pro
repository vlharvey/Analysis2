;
; latitude-altitude slices of ECMWF water

@stddat
@kgmt
@ckday
@kdate
@rd_ecmwf_nc

loadct,38
mcolor=byte(!p.color)
device,decompose=0
dir='/aura3/data/ECMWF_data/Datfiles/ecmwf_'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
setplot='x'
read,'setplot ',setplot
lstmn=11
lstdy=1
lstyr=91
ledmn=11
leddy=1
ledyr=91
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      ECMWF Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
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
xorig=[0.1,0.4,0.7,0.15,0.55]
yorig=[0.6,0.6,0.6,0.15,0.15]
xlen=0.25
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
      file=dir+smn+'_'+sdy+'_'+syr+'_12Z.nc'
      rd_ecmwf_nc,file,nc,nr,nl,alon,alat,th,pv2,p2,msf2,u2,v2,q2,qdf2,sh2,o3,iflg
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
             'yz_ecmwf_',imn,'_',idy,'_',iyr)+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.charsize=1.0
      endif

; Declare plotting arrays
      mpv2=0.*pv2
      for k=0L,nl-1L do mpv2(*,*,k)=pv2(*,*,k)*((th(k)/300.))^(-9./2.)
      index=where(pv2 eq 1.e12)
      if index(0) ne -1 then mpv2(index)=1.e12
      edat=reform(mpv2(*,0,*),nr,nl)
      pvdat=reform(pv2(*,0,*),nr,nl)
      water=reform(sh2(*,0,*),nr,nl)*1.e6
      ozone=reform(o3(*,0,*),nr,nl)*1.e6

; Find data range for autoscaling contours
      erase
index=where(edat eq 1.00000e+12 or edat eq 0.)
if index(0) ne -1 then edat(index)=-9999.
      index=where(edat ne -9999.)
      minval=min(edat(index))
      maxval=max(edat(index))
      nlvs=30
      cint=(maxval-minval)/float(nlvs)

; Autoscale if scale values for parameter/level are not defined
      level=minval+cint*findgen(nlvs)
      col1=1+indgen(nlvs)*icolmax/nlvs
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))

; Mercator projection
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
vedge=5.e-6
      xyouts,.3,.9,'ECMWF '+date,/normal,charsize=2
      contour,edat,alat,th,levels=level,xrange=[-90.,90.],yrange=[min(th),600.],$
              /cell_fill,c_color=col1,/noeras,title='MPV',min_value=-9999.
      contour,edat,alat,th,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0,/noeras,min_value=-9999.
      contour,edat,alat,th,levels=vedge,/overplot,c_color=mcolor,thick=5

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

      xmn=xorig(1)
      xmx=xorig(1)+xlen
      ymn=yorig(1)
      ymx=yorig(1)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
      level=0.5*findgen(30)
      contour,ozone,alat,th,levels=level,xrange=[-90.,90.],yrange=[min(th),600.],$
              /cell_fill,c_color=col1,/noeras,title='Ozone',min_value=0.
      contour,ozone,alat,th,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0,/noeras,min_value=0.
      contour,ozone,alat,th,levels=0.5,/overplot,c_color=mcolor,thick=3
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
      level=1.+0.1*findgen(30)
      contour,water,alat,th,levels=level,xrange=[-90.,90.],yrange=[min(th),600.],$
              /cell_fill,c_color=col1,/noeras,title='Sp. Humidity',min_value=0.
      contour,water,alat,th,levels=level,/overplot,c_color=[0],$
              c_labels=0*level,c_linestyle=level lt 0,/noeras,min_value=0.
      contour,water,alat,th,levels=5.0,/overplot,c_color=mcolor,thick=3
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

      !type=2^2+2^3
      xmn=xorig(3)
      xmx=xorig(3)+0.3
      ymn=yorig(3)
      ymx=yorig(3)+0.3
      set_viewport,xmn,xmx,ymn,ymx
      omin=min(ozone)
      omax=max(ozone)
      index=where(pvdat gt 0. and ozone gt 0.) 
      plot,ozone(index),pvdat(index),psym=3,yrange=[1.e-8,0.1],/ylog,/xlog,xrange=[1.e-2,omax],$
           ytitle='Potential Vorticity',xtitle='Ozone',title='PV/O3 Scatterplot'
      index=where(edat ge vedge and pvdat gt 0. and ozone gt 0.)
      oplot,ozone(index),pvdat(index),psym=1,color=mcolor*.3

      xmn=xorig(4)
      xmx=xorig(4)+0.3
      ymn=yorig(4)
      ymx=yorig(4)+0.3
      set_viewport,xmn,xmx,ymn,ymx
      omin=min(water)
      omax=max(water)
      omax=5.
      index=where(pvdat gt 0. and water gt 0.)
      plot,water(index),pvdat(index),psym=3,yrange=[1.e-8,0.1],/ylog,xrange=[omin,10.],$
           ytitle='Potential Vorticity',xtitle='Specific Humidity',title='PV/H2O Scatterplot'
      index=where(edat ge vedge and pvdat gt 0. and water gt 0.)
      oplot,water(index),pvdat(index),psym=1,color=mcolor*.3

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
