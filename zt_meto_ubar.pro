;
; compute daily Arctic MetO T, Mark at 4 IPY LIDAR sites
; plot altitude-time sections
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nwp

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
dir='/aura3/data/UKMO_data/Datfiles/ukmo_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1L & lstdy=1L & lstyr=7L
ledmn=12L & leddy=31L & ledyr=7L
lstday=0L & ledday=0L
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
kday=long(ledday-lstday+1L)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=-1L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
      kcount=kcount+1L
      ifile=mon(imn-1)+sdy+'_'+uyr
      ifile='/aura7/harvey/UKMO_data/Datfiles/ukmo-nwp-strat_gbl-std_'+sdate+'12_u-v-gph-t-w_uars.nc'
;
; read MetO data
;
      rd_ukmo_nwp,ifile,nc,nr,nc1,nr1,nlv,wlon,alon,wlat,alat,p,z3d,t3d,u3d,v3d,iflg

      if icount eq 0L then begin
         u_bar_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
      endif
      sdate_all(kcount)=sdate
;
; compute daily zonal mean
;
;print,wlat
rlat=61.2500
;read,'Enter latitude',rlat
index=where(wlat eq rlat)
ilat=index(0)
slat=strcompress(wlat(ilat),/remove_all)
      for k=0L,nlv-1L do begin
          u_bar_zt(icount,k)=mean(u3d(*,ilat,k))
      ENDFOR
      icount=icount+1
goto, jump

plotit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_ubar_'+sdate0+'-'+sdate1+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; extract date labels
;
    syear=strmid(sdate_all,0,4)
    smon=strmid(sdate_all,4,2)
    sday=strmid(sdate_all,6,2)
    xindex=where(sday eq '15',nxticks)
    xlabs=smon(xindex)
;
; plot time-altitude sections of zonal mean wind speed
;
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    nlvls=21
    col1=1+indgen(nlvls)*icolmax/nlvls
    level=-100.+10.*findgen(nlvls)
    contour,u_bar_zt,1.+findgen(kday),p,/noeras,xrange=[1.,kday],yrange=[max(p),min(p)],/ylog,charsize=1.5,color=0,$
            ytitle='Pressure (hPa)',xtitle=syear(0),title='MetO Ubar at '+slat+' Latitude',/fill,c_color=col1,levels=level,$
            xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
    index=where(level lt 0)
    contour,u_bar_zt,1.+findgen(kday),p,levels=level(index),color=0,/follow,/overplot,c_labels=0*index,c_linestyle=5
    index=where(level gt 0)
    contour,u_bar_zt,1.+findgen(kday),p,levels=level(index),color=mcolor,/follow,/overplot,c_labels=0*index
    contour,u_bar_zt,1.+findgen(kday),p,levels=[0],color=0,/follow,/overplot,c_labels=[0],thick=3
    imin=min(level)
    imax=max(level)
    ymnb=yorig(0)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xorig(0),xorig(0)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MetO Ubar (m/s)',charsize=2,charthick=2
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim zt_ubar_'+sdate0+'-'+sdate1+'.ps -rotate -90 zt_ubar_'+sdate0+'-'+sdate1+'.jpg'
;      spawn,'/usr/bin/rm zt_ubar_'+sdate0+'-'+sdate1+'.ps'
    endif
end
