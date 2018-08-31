;
; plot zonal mean temperature for ECMWF, UKMO and differences

@stddat
@kgmt
@ckday
@kdate
@rd_ecmwf
@rd_ukmo

loadct,38
mcolor=byte(!p.color)
device,decompose=0
dir='/aura5/harvey/ECMWF_data/Datfiles/ecmwf_'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
setplot='x'
read,'setplot ',setplot
lstmn=11
lstdy=1
lstyr=91
ledmn=10
leddy=31
ledyr=92
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
xorig=[0.1,0.1,0.1,0.5,0.5,0.5]
yorig=[0.75,0.425,0.1,0.75,0.425,0.1]
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
      if icount eq 0L and idy ne 1L then stop,'Begin on first day'

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
;
;***Read ECMWF data
      file=dir+smn+'_'+sdy+'_'+syr+'_12Z.dat'
      rd_ecmwf,file,iflg,nc,nr,nl,alon,alat,press,pv,gp,tp,uu,vv,ww,sh,oz
      ss=sqrt(uu^2.0+vv^2.0)
;
;***Read UKMO data
      file='/aura3/data/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nlg,nlat,nlv,ulon,ulat,wlon,wlat,p,g3d,t3d,u3d,v3d
      if iflg ne 0 then goto, jump
      s3d=sqrt(u3d^2.0+v3d^2.0)
;
; declare mean arrays on the first day of each month
;
      if idy eq 1L then begin
         pvmean=pv
         gpmean=gp
         tpmean=tp
         ssmean=ss
         shmean=sh
         ozmean=oz
         g3dmean=g3d
         t3dmean=t3d
         s3dmean=s3d
      endif
;
; average other days in month
;
      if idy gt 1L then begin
         pvmean=pvmean+pv
         gpmean=gpmean+gp
         tpmean=tpmean+tp
         ssmean=ssmean+ss
         shmean=shmean+sh
         ozmean=ozmean+oz
         g3dmean=g3dmean+g3d
         t3dmean=t3dmean+t3d
         s3dmean=s3dmean+s3d
      endif 
;
; Plot on last day of month
;
      if idy eq mday(imn-1) then begin
         pvmean=pvmean/float(mday(imn-1))
         gpmean=gpmean/float(mday(imn-1))
         tpmean=tpmean/float(mday(imn-1))
         ssmean=ssmean/float(mday(imn-1))
         shmean=shmean/float(mday(imn-1))
         ozmean=ozmean/float(mday(imn-1))
         g3dmean=g3dmean/float(mday(imn-1))
         t3dmean=t3dmean/float(mday(imn-1))
         s3dmean=s3dmean/float(mday(imn-1))
;
; initialize zonal mean arrays
;
         etemp=fltarr(nr,nl)
         espeed=fltarr(nr,nl)
         for k=0,nl-1 do begin
         for j=0,nr-1 do begin
             etemp(j,k)=total(tpmean(*,j,k))/float(nc)
             espeed(j,k)=total(ssmean(*,j,k))/float(nc)
         endfor
         endfor
         utemp=fltarr(nlat,nlv)
         uspeed=fltarr(nlat-1,nlv)
         for k=0,nlv-1 do begin
         for j=0,nlat-1 do begin
             utemp(j,k)=total(t3dmean(*,j,k))/float(nlg)
             if j lt nlat-1 then uspeed(j,k)=total(s3dmean(*,j,k))/float(nlg)
         endfor
         endfor
;
; interpolate EC to UKMO pressures (lats are already the same)
;
         lpress=alog(press)
         lp=alog(p)
         etempu=fltarr(nlat,nlv)
         for j=0,nlat-1 do begin
         for k=0,nlv-1 do begin
             for kk=0,nl-2 do begin
                 kp1=kk+1
                 if lpress(kk) ge lp(k) and lpress(kp1) le lp(k) then begin
                 scale=(lpress(kk)-lp(k))/(lpress(kk)-lpress(kp1))
                 etempu(j,k)=etemp(j,kk)+scale*(etemp(j,kp1)-etemp(j,kk))
                 endif
             endfor
         endfor
         endfor
;
; reverse ukmo latitudes
         for k=0,nlv-1 do utemp(*,k)=reverse(utemp(*,k))
         for k=0,nlv-1 do uspeed(*,k)=reverse(uspeed(*,k))
         ulat=reverse(ulat)
         wlat=reverse(wlat)

         if setplot eq 'ps' then begin
            set_plot,'ps'
            xsize=nxdim/100.
            ysize=nydim/100.
            !p.font=0
            device,font_size=9
            device,/landscape,bits=8,filename=$
                   string(FORMAT='(a5,i2.2,a1,i4)','ecmwf+ukmo_',imn,'_',iyr)+'.ps'
            device,/color
            device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                   xsize=xsize,ysize=ysize
         endif

; Find data range for autoscaling contours
         erase
         minval=180.
         nlvs=30
         cint=5.0
         level=minval+cint*findgen(nlvs)
         col1=1+indgen(nlvs)*icolmax/nlvs
         date=strcompress(string(FORMAT='(A3,A1,I4)',month(imn-1),' ',iyr))
         xyouts,.55,.25,date,/normal,charsize=2.5
         xmn=xorig(0)
         xmx=xorig(0)+xlen
         ymn=yorig(0)
         ymx=yorig(0)+ylen
         set_viewport,xmn,xmx,ymn,ymx
         !type=2^2+2^3
         contour,etemp,alat,press,levels=level,/cell_fill,c_color=col1,/noeras,$
                 /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='ERA Temperature',$
                 charsize=1.5,xticks=6
         contour,etemp,alat,press,levels=level,/overplot,c_color=[0],$
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
         xmn=xorig(1)
         xmx=xorig(1)+xlen
         ymn=yorig(1)
         ymx=yorig(1)+ylen
         set_viewport,xmn,xmx,ymn,ymx
         !type=2^2+2^3
         contour,utemp,ulat,p,levels=level,/cell_fill,c_color=col1,/noeras,$
                 /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='MetO Temperature',$
                 charsize=1.5,xticks=6
         contour,utemp,ulat,p,levels=level,/overplot,c_color=[0],$
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
         dtemp=etempu
         index=where(etempu gt 0.)
         dtemp(index)=etempu(index)-utemp(index)
         level=-10.+findgen(20)
         contour,dtemp,ulat,p,levels=level,/cell_fill,c_color=col1,/noeras,$
                 /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='EC-MetO Tp',$
                 charsize=1.5,xticks=6
         contour,dtemp,ulat,p,levels=level,/overplot,c_color=[0],$
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

         minval=0.
         nlvs=30
         cint=4.0
         level=minval+cint*findgen(nlvs)
         col1=1+indgen(nlvs)*icolmax/nlvs
         xmn=xorig(3)
         xmx=xorig(3)+xlen
         ymn=yorig(3)
         ymx=yorig(3)+ylen
         set_viewport,xmn,xmx,ymn,ymx
         !type=2^2+2^3
         contour,espeed,alat,press,levels=level,/cell_fill,c_color=col1,/noeras,$
                 /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='ERA Isotachs',$
                 charsize=1.5,xticks=6
         contour,espeed,alat,press,levels=level,/overplot,c_color=[0],$
                 c_labels=0*level,c_linestyle=level lt 0,/noeras
         imin=min(level)
         imax=max(level)
         ymnb=yorig(3) -cbaryoff
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
         xmn=xorig(4)
         xmx=xorig(4)+xlen
         ymn=yorig(4)
         ymx=yorig(4)+ylen
         set_viewport,xmn,xmx,ymn,ymx
         !type=2^2+2^3
         contour,uspeed,wlat,p,levels=level,/cell_fill,c_color=col1,/noeras,$
                 /ylog,yrange=[1000.,1.],xrange=[-90.,90.],title='MetO Isotachs',$
                 charsize=1.5,xticks=6
         contour,uspeed,wlat,p,levels=level,/overplot,c_color=[0],$
                 c_labels=0*level,c_linestyle=level lt 0
         imin=min(level)
         imax=max(level)
         ymnb=yorig(4) -cbaryoff
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
         if setplot eq 'ps' then device, /close
      endif
      icount=icount+1L
goto, jump
end
