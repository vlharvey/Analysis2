;
; compute daily Arctic GEOS-5 T, Mark at 4 IPY LIDAR sites
; plot altitude-time sections
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
device,decompose=0
!p.background=icolmax
setplot='ps'
;read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.6,0.15,0.6]
yorig=[0.55,0.55,0.15,0.15]
xlen=0.3
ylen=0.3
cbaryoff=0.06
cbarydel=0.01
stimes=[$
'_0000.V01.',$
'_0600.V01.',$
'_1200.V01.',$
'_1800.V01.']
slabs=['00Z','06Z','12Z','18Z']
ntimes=n_elements(stimes)
!noeras=1
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
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
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=10L & lstdy=15L & lstyr=7L
ledmn=3L & leddy=12L & ledyr=8L
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
kday=(ledday-lstday+1L)*4L
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
      date=syr+smn+sdy
      print,date
      ifile=mon(imn-1)+sdy+'_'+uyr
;
; loop over daily output times
;
      for itime=0L,ntimes-1L do begin

      rd_geos5_nc3,dir+date+stimes(itime)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      kcount=kcount+1L
      if iflag eq 1 then goto,jump
;
; read new marker field
;
      ncid=ncdf_open(dir+date+stimes(itime)+'nc4')
      mark2new=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,mark2new
      ncdf_close,ncid

      if icount eq 0L then begin
         temp_pf_zt=fltarr(kday,nth)
         temp_al_zt=fltarr(kday,nth)
         temp_eu_zt=fltarr(kday,nth)
         temp_so_zt=fltarr(kday,nth)
         u_pf_zt=fltarr(kday,nth)
         u_al_zt=fltarr(kday,nth)
         u_eu_zt=fltarr(kday,nth)
         u_so_zt=fltarr(kday,nth)
         date_all=strarr(kday)
         icount=1
      endif
      date_all(kcount)=date
;
; temperature
      temp2=0.*p2
      for k=0L,nth-1L do temp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^.286
;
; interpolate temperature to Poker Flat (65N, 147W)
;
      xpt=[-147.,16.,-86.,-51.] & ypt=[65.,69.,80.,67.]
      npt=n_elements(xpt)
      FOR IPT=0l,NPT-1l DO BEGIN

      if xpt(ipt) lt alon(0) then xpt(ipt)=xpt(ipt)+360.
      for i=0L,nc-1L do begin
          ip1=i+1
          if i eq nc-1L then ip1=0L
          xlon=alon(i)
          xlonp1=alon(ip1)
          if i eq nc-1L then xlonp1=360.+alon(ip1)
          if xpt(ipt) ge xlon and xpt(ipt) le xlonp1 then begin
             xscale=(xpt(ipt)-xlon)/(xlonp1-xlon)
             goto,jumpx
          endif
      endfor
jumpx:
      for j=0L,nr-2L do begin
          jp1=j+1
          xlat=alat(j)
          xlatp1=alat(jp1)
          if ypt(ipt) ge xlat and ypt(ipt) le xlatp1 then begin
              yscale=(ypt(ipt)-xlat)/(xlatp1-xlat)
              goto,jumpy
          endif
      endfor
jumpy:
      pj1=temp2(j,i,*)+xscale*(temp2(j,ip1,*)-temp2(j,i,*))
      pjp1=temp2(jp1,i,*)+xscale*(temp2(jp1,ip1,*)-temp2(jp1,i,*))
      if ipt eq 0L then temp_pf_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 1L then temp_al_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 2L then temp_eu_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 3L then temp_so_zt(kcount,*)=pj1+yscale*(pjp1-pj1)

      pj1=mark2(j,i,*)+xscale*(mark2(j,ip1,*)-mark2(j,i,*))
      pjp1=mark2(jp1,i,*)+xscale*(mark2(jp1,ip1,*)-mark2(jp1,i,*))
      if ipt eq 0L then u_pf_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 1L then u_al_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 2L then u_eu_zt(kcount,*)=pj1+yscale*(pjp1-pj1)
      if ipt eq 3L then u_so_zt(kcount,*)=pj1+yscale*(pjp1-pj1)

      ENDFOR	; loop over 4 lidar stations
      ENDFOR	; loop over 4 daily times
goto, jump

plotit:
date0=date_all(0)
date1=date_all(n_elements(date_all)-1)
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='ZT_geos5.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
erase
slat=['65N','69N','80N','67N']
slon=['147W','16E','86W','51W']
slab=['Poker Flat','Alomar','Eureka','Sondrestrom']
;
; interpolate across "small" gaps in time
;
for ipt=0L,3L do begin
if ipt eq 0L then temp_zt=temp_pf_zt
if ipt eq 1L then temp_zt=temp_al_zt
if ipt eq 2L then temp_zt=temp_eu_zt
if ipt eq 3L then temp_zt=temp_so_zt
if ipt eq 0L then u_zt=u_pf_zt
if ipt eq 1L then u_zt=u_al_zt
if ipt eq 2L then u_zt=u_eu_zt
if ipt eq 3L then u_zt=u_so_zt
for k=0,nth-1 do begin
    dlev=reform(temp_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 15.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    temp_zt(*,k)=dlev

    dlev=reform(u_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 15.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump2
               endif
           endfor
jump2:
        endif
    endfor
    u_zt(*,k)=dlev
endfor
;
; plot time-altitude sections of Temperature and zonal wind speed at 4 LIDAR locations
;
    level=175.+5.*findgen(24)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*mcolor/nlvls
    !type=2^2+2^3
    xmn=xorig(ipt)
    xmx=xorig(ipt)+xlen
    ymn=yorig(ipt)
    ymx=yorig(ipt)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    syr=strmid(date_all,0,4)
    smn=strmid(date_all,4,2)
    sdy=strmid(date_all,6,2)
    xindex=where(sdy eq '15',nxticks)
    xlabs=smn(xindex)
    ylab=' '
    if ipt eq 0L or ipt eq 2L then ylab='Theta (K)'
    contour,temp_zt,findgen(kday),th,charsize=2,/noeras,yrange=[300.,max(th)],charthick=2,$
            ytitle=ylab,levels=level,c_color=col1,/cell_fill,color=0,title=slab(ipt)+' '+'('+slon(ipt)+','+slat(ipt)+')',$
            xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
loadct,0
;   contour,temp_zt,findgen(kday),th,levels=level(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
    contour,u_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=[0.1],color=0,thick=4
    contour,u_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=[-0.1],color=mcolor,c_linestyle=5,thick=4
x2d=0.*u_zt
y2d=0.*u_zt
for i=0L,kday-1 do y2d(i,*)=th
for k=0L,nth-1 do x2d(*,k)=findgen(kday)
index=where(u_zt gt 0.)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=0,symsize=0.25
index=where(u_zt lt 0.)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=mcolor,symsize=0.25

    for ii=0L,nxticks-1L do xyouts,xindex(ii),-100.,smon(long(xlabs(ii))-1),/data,color=0,charsize=1.5,charthick=2,alignment=0.5
    for ii=0L,nxticks-1L do begin
        plots,xindex(ii),300.
        plots,xindex(ii),0.,/continue,color=0,thick=2,/data
    endfor

    endfor
    xyouts,.5,.925,'2008-2009',/normal,color=0,charsize=3,charthick=2,alignment=0.5

    imin=min(level)
    imax=max(level)
    ymnb=yorig(3)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xorig(0),xorig(1)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='GEOS-5.2 Temperature (K)',charsize=2,charthick=2
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
       spawn,'convert -trim ZT_geos5.ps -rotate -90 ZT_geos5.png'
       spawn,'/usr/bin/rm ZT_geos5.ps'
    endif
end
