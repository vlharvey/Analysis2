;
; MERRA version.  
;
; zonal mean temperature and zonal wind from theta data
; store 2-D arrays (day vs. theta) of dT/dy and Ubar 
; for both hemispheres and entire data record in one IDL save file
;
@stddat
@kgmt
@ckday
@kdate

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=700
nydim=700
xorig=[.15,.15]
yorig=[.55,.15]
xlen=0.7
ylen=0.35
cbaryoff=0.075
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
!noeras=1
;goto,plotit
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1L & lstdy=1L & lstyr=1979L
ledmn=2L & leddy=28L & ledyr=2015L
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd=lonarr(nfile)
syyyymmdd=strarr(nfile)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto, saveit
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
      yyyymmdd(icount)=long(syr+smn+sdy)
      syyyymmdd(icount)=sdate
;
; read data
;
      restore,dir+sdate+'.sav'		
      alat=LATITUDE_WACCM
      nth=n_elements(pressure)
      nr=n_elements(LATITUDE_WACCM)
      nc=n_elements(LONGITUDE_WACCM)
;
      if icount eq 0L then begin
;
; WMO minor warming is T85 - T60 > 0. major is U65 < 0
; 
y60n=where(abs(alat-60) le 1.)
y60s=where(abs(alat+60) le 1.)
y65n=where(abs(alat-65) le 1.)
y65s=where(abs(alat+65) le 1.)
y85n=where(abs(alat-85) le 1.)
y85s=where(abs(alat+85) le 1.)
; 
; dTdy and Ubar in NH and SH
;
         nh_dtdy=-9999.+0.*fltarr(nfile,nth)
         sh_dtdy=-9999.+0.*fltarr(nfile,nth)
         nh_ubar=-9999.+0.*fltarr(nfile,nth)
         sh_ubar=-9999.+0.*fltarr(nfile,nth)
      endif
;
; calculate zonal mean temperature and zonal wind
;
      uzm=fltarr(nr,nth)
      tzm=fltarr(nr,nth)
      for k=0,nth-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tgrd(*,j,k))/float(nc)
              uzm(j,k)=total(ugrd(*,j,k))/float(nc)
          endfor
t60n=mean(tzm(y60n,k))
t60s=mean(tzm(y60s,k))
t85n=mean(tzm(y85n,k))
t85s=mean(tzm(y85s,k))
          if min(tzm(y60n,k)) ne 0. and min(tzm(y85n,k)) ne 0. then nh_dtdy(icount,k)=t85n-t60n
          if min(tzm(y60s,k)) ne 0. and min(tzm(y85s,k)) ne 0. then sh_dtdy(icount,k)=t85s-t60s
if abs(nh_dtdy(icount,k)) gt 100. then stop
if abs(sh_dtdy(icount,k)) gt 100. then stop

u65n=mean(uzm(y65n,k))
u65s=mean(uzm(y65s,k))
          if min(tzm(y60n,k)) ne 0. and min(tzm(y85n,k)) ne 0. then nh_ubar(icount,k)=u65n
          if min(tzm(y60s,k)) ne 0. and min(tzm(y85s,k)) ne 0. then sh_ubar(icount,k)=u65s
      endfor
;stop
      skip:
      icount=icount+1L
goto,jump
;
; save file
;
saveit:
icount=n_elements(yyyymmdd)
comment='dT/dy=Tbar90-Tbar60 and Ubar60 is zonal mean zonal wind at 60N/S'
save,file='Save_files/MERRA_dTdy_Ubar_SSW_Climo.sav',yyyymmdd,syyyymmdd,pressure,$
     nh_dtdy,sh_dtdy,nh_ubar,sh_ubar,comment,icount
;
; plot
;
plotit:
restore,'Save_files/MERRA_dTdy_Ubar_SSW_Climo.sav'
;
; save postscript version
;
sdate0=syyyymmdd(0)
sdate1=syyyymmdd(icount-1)
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='merra_dtdy_ubar_ssw_'+sdate0+'-'+sdate1+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.
   !p.charthick=2.
   !p.charsize=1.5
endif
syr=strmid(syyyymmdd,2,2)
smn=strmid(syyyymmdd,4,2)
sdy=strmid(syyyymmdd,6,2)
xindex=where(sdy eq '01' and smn eq '01',nxticks)
xlabs=syr(xindex)

erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=-50.+5.*findgen(21)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*mcolor/nlvls
!type=2^2+2^3
contour,nh_dtdy,findgen(icount),pressure,yrange=[max(pressure),min(pressure)],/ylog,/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,ytitle='Theta (K)',title='MERRA Tbar 90N-60N',$
        xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
;index=where(level lt 0.)
;contour,nh_dtdy,findgen(icount),pressure,/noeras,levels=level(index),/follow,color=0,$
;        c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_dtdy,findgen(icount),pressure,/noeras,levels=level(index),/follow,color=mcolor,/overplot
contour,nh_dtdy,findgen(icount),pressure,/noeras,levels=[0],/follow,color=0,/overplot
set_viewport,xmx+cbaryoff,xmx+cbaryoff+cbarydel,ymn,ymx
!type=2^2+2^3+2^5
omin=min(level)
omax=max(level)
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],title='(K/deg)',color=0,charsize=1.5
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=-120.+10.*findgen(25)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*mcolor/nlvls
!type=2^2+2^3
contour,nh_ubar,findgen(icount),pressure,yrange=[max(pressure),min(pressure)],/ylog,/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,ytitle='Theta (K)',title='MERRA Ubar at 60N',$
        xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
;index=where(level lt 0.)
;contour,nh_ubar,findgen(icount),pressure,/noeras,levels=level(index),$
;        /follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_ubar,findgen(icount),pressure,/noeras,levels=level(index),/follow,$
;        color=mcolor,/overplot
contour,nh_ubar,findgen(icount),pressure,/noeras,levels=[0],/follow,color=0,/overplot
xindex=where(sdy eq '15' and smn eq '07',nxticks)
xyrs=syr(xindex)
for i=0L,nxticks-1L do xyouts,xindex(i),-200.,xyrs(i),/data,color=0,charsize=2,alignment=0.5
set_viewport,xmx+cbaryoff,xmx+cbaryoff+cbarydel,ymn,ymx
!type=2^2+2^3+2^5
omin=min(level)
omax=max(level)
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],title='(m/s)',color=0,charsize=1.5
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

print,min(nh_dtdy),max(nh_dtdy),min(nh_ubar),max(nh_ubar)

if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim merra_dtdy_ubar_ssw_'+sdate0+'-'+sdate1+'.ps -rotate -90 '+$
         'merra_dtdy_ubar_ssw_'+sdate0+'-'+sdate1+'.jpg'
endif

end
