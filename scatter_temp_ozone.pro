;
; read ERA40 data and plot ozonet versus temperature
;
@read_ecmwf
@stddat
@kgmt
@ckday
@kdate

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
ncfile='/aura3/data/ECMWF_data/Datfiles/ecmwf_feb-may_1984.nc'
ncfile='/aura3/data/ECMWF_data/Datfiles/ecmwf_jan_1984.nc'
read_ecmwf,ncfile,nc,nr,nl,nt,alon,alat,press,date,$
     z2,t2,u2,v2,q2,omega2,o32,varname,varunit
;
; initial date
;
result=strsplit(varunit(3),/extract)
idate=result(2)
result=strsplit(idate,'-',/extract)
iyr=long(result(0))
imn=long(result(1))
idy=long(result(2))
;
; loop over days
;
z = kgmt(imn,idy,iyr,iday)
for n=0,nt-1L do begin
    kdate,float(iday),iyr,imn,idy
    ckday,iday,iyr
    print,iyr,imn,idy
    syr=string(FORMAT='(i4)',iyr)
    smn=string(FORMAT='(i2.2)',imn)
    sdy=string(FORMAT='(i2.2)',idy)
    sdate=syr+smn+sdy

    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='scatter_temp_ozone_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

    erase
    !type=2^2+2^3
    xmn=0.05 & xmx=0.45 & ymn=0.30 & ymx=0.70
    set_viewport,xmn,xmx,ymn,ymx
    xyouts,.35,.8,sdate,charsize=3,/normal
    ozone=o32(*,nr/2:nr-1,*,n)*1.e6
    temp=t2(*,nr/2:nr-1,*,n)
    height=z2(*,nr/2:nr-1,*,n)/1000.
    plot,ozone,temp,psym=3,yrange=[170.,320.],xrange=[0.,17.],title='NH',$
         xtitle='Ozone (ppmv)',ytitle='Temperature'
    index=where(height lt 5)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.05
    index=where(height ge 5. and height lt 10.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.15
    index=where(height ge 10. and height lt 15.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.20
    index=where(height ge 15. and height lt 20.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.25
    index=where(height ge 20. and height lt 25.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.35
    index=where(height ge 25. and height lt 30.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.45
    index=where(height ge 30. and height lt 35.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.55
    index=where(height ge 35. and height lt 40.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.65
    index=where(height ge 40. and height lt 45.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.75
    index=where(height ge 45. and height lt 50.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.85
    index=where(height ge 50.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.95
;
; SH
;
    !type=2^2+2^3
    xmn=0.55 & xmx=0.95 & ymn=0.30 & ymx=0.70
    set_viewport,xmn,xmx,ymn,ymx
    ozone=o32(*,0:nr/2-1,*,n)*1.e6
    temp=t2(*,0:nr/2-1,*,n)
    height=z2(*,0:nr/2-1,*,n)/1000.
    plot,ozone,temp,psym=3,yrange=[170.,320.],xrange=[0.,17.],title='SH',$
         xtitle='Ozone (ppmv)',ytitle='Temperature'
    index=where(height lt 5)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.05
    index=where(height ge 5. and height lt 10.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.15
    index=where(height ge 10. and height lt 15.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.20
    index=where(height ge 15. and height lt 20.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.25
    index=where(height ge 20. and height lt 25.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.35
    index=where(height ge 25. and height lt 30.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.45
    index=where(height ge 30. and height lt 35.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.55
    index=where(height ge 35. and height lt 40.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.65
    index=where(height ge 40. and height lt 45.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.75
    index=where(height ge 45. and height lt 50.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.85
    index=where(height ge 50.)
    if index(0) ne -1 then oplot,ozone(index),temp(index),psym=3,color=icolmax*.95
    imin=0. & imax=max(height)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,0.05,0.95,ymnb,ymxb
    nlvls=30
    col1=1+indgen(nlvls)*icolmax/nlvls
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],charsize=1.5,xrange=[imin,imax],$
          xtitle='Geopotential Height (km)'
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor

    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert scatter_temp_ozone_'+sdate+'.ps -rotate -90 '+$
             'scatter_temp_ozone_'+sdate+'.jpg'
       spawn,'/usr/bin/rm scatter_temp_ozone_'+sdate+'.ps'
    endif
    iday=iday+1L
stop
endfor
end
