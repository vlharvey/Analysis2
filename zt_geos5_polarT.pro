;
; GEOS-5 version
; compute daily average Arctic Temp and Mark
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
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.55,0.15]
xlen=0.6
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
smon=['J','J','J','J','F','F','F','F','M','M','M','M','M','A','M','J','J','A','S','O','N','D']
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
goto,plotit
lstmn=1L & lstdy=2L & lstyr=9L
ledmn=3L & leddy=19L & ledyr=9L
lstday=0L & ledday=0L
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
      FOR ITIME=0l,NTIMES-1l DO BEGIN

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
             temp_nh_zt=fltarr(kday,nth)
             temp_sh_zt=fltarr(kday,nth)
             u_nh_zt=fltarr(kday,nth)
             u_sh_zt=fltarr(kday,nth)
             date_all=strarr(kday)
             time_all=fltarr(kday)
             icount=1
          endif
          date_all(kcount)=date
          time_all(kcount)=itime*6.
;
; temperature
;
          temp2=0.*p2
          for k=0L,nth-1L do temp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^.286
;
; average 60-90N
;
          x2d=fltarr(nr,nc)
          y2d=fltarr(nr,nc)
          for i=0L,nc-1 do y2d(*,i)=alat
          for j=0L,nr-1 do x2d(j,*)=alon
          for k=0L,nth-1L do begin
              temp1=reform(temp2(*,*,k))
              u1=reform(u2(*,*,k))
              nhindex=where(y2d ge 60. and temp1 ne 0.,nnh)
              shindex=where(y2d le -60. and temp1 ne 0.,nsh)
              temp_nh_zt(kcount,k)=mean(temp1(nhindex))
              temp_sh_zt(kcount,k)=mean(temp1(shindex))
              u_nh_zt(kcount,k)=mean(u1(nhindex))
              u_sh_zt(kcount,k)=mean(u1(shindex))
          endfor
      ENDFOR	; loop over 4 daily times
goto, jump

plotit:
;
; comment out if generating quick plot
;
;SAVE, /VARIABLES, FILENAME = 'zt_geos5_polarT.sav'
;
restore,'zt_geos5_polarT.sav'
;
; smooth u
;
u_nh_zt=smooth(u_nh_zt,3)
u_sh_zt=smooth(u_sh_zt,3)
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
   device,/landscape,bits=8,filename='zt_geos5_polarT.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
erase
slab=['NH Mean >60N','SH Mean <60S']
;
; plot time-altitude sections of Temperature and zonal wind speed in NH and SH
;
tlevel=175.+5.*findgen(24)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*mcolor/nlvls
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
syr=strmid(date_all,0,4)
smn=strmid(date_all,4,2)
sdy=strmid(date_all,6,2)
xindex=where(sdy eq '15' and time_all eq 12. ,nxticks)
xlabs=smn(xindex)+'/'+strmid(syr(xindex),2,2)
ylab='Theta (K)'
contour,temp_nh_zt,findgen(kday),th,charsize=2,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=ylab,levels=tlevel,c_color=col1,/cell_fill,color=0,title=slab(0),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.5,charthick=2,alignment=0.5
loadct,0
;contour,temp_nh_zt,findgen(kday),th,tlevels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
ulevel=-100.+10.*findgen(21)
index=where(ulevel lt 0.)
contour,u_nh_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=ulevel(index),color=0,thick=4,c_linestyle=5
index=where(ulevel gt 0.)
contour,u_nh_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=ulevel(index),color=mcolor,thick=4

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,temp_sh_zt,findgen(kday),th,charsize=2,/noeras,yrange=[300.,max(th)],charthick=2,$
        ytitle=ylab,levels=tlevel,c_color=col1,/cell_fill,color=0,title=slab(1),$
        xrange=[0.,kday-1],xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
for ii=0L,nxticks-1L do begin
    plots,xindex(ii),300.
    plots,xindex(ii),0.,/continue,color=0,thick=2,/data
endfor
for ii=0L,nxticks-1L do xyouts,xindex(ii),-120.,xlabs(ii),/data,color=0,charsize=1.5,charthick=2,alignment=0.5
loadct,0
;contour,temp_sh_zt,findgen(kday),th,tlevels=tlevel(0:*:2),color=150,/follow,/noeras,/overplot,min_value=0.
loadct,38
ulevel=-100.+10.*findgen(21)
index=where(ulevel lt 0.)
contour,u_sh_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=ulevel(index),color=0,thick=4,c_linestyle=5
index=where(ulevel gt 0.)
contour,u_sh_zt,findgen(kday),th,charsize=1.5,/noeras,/overplot,levels=ulevel(index),color=mcolor,thick=4

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xorig(1),xorig(1)+xlen,ymnb,ymxb
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
       spawn,'convert -trim zt_geos5_polarT.ps -rotate -90 zt_geos5_polarT.jpg'
       spawn,'/usr/bin/rm zt_geos5_polarT.ps'
endif
end
