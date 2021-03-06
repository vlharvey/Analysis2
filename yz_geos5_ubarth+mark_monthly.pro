;
; monthly mean zonal mean
; GEOS-4 version
;
; plot daily Ubar from theta data
; superimpose zonal mean marker
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,38
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.1
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
month=['January','February','March','April','May','June','July','August','September','October','November','December']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
years=[2005,2006,2007,2008]
nyear=n_elements(years)
lstdy=1
for iyear=0L,nyear-1L do begin
lstyr=years(iyear)
syear=strcompress(lstyr,/remove_all)
ledyr=lstyr

for ii=0L,11L do begin
lstmn=ii+1
ledmn=lstmn
leddy=mday(ii)
if lstyr mod 4 eq 0 and ii eq 1 then leddy=leddy+1	; leap day

lstday=0
ledday=0
icount=0L
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS-4 Version '
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
nday=ledday-lstday+1L
;
; Compute initial Julian date
;
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
      if ndays gt ledday then goto,plotit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-4 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      rd_geos5_nc3_meto,dir+sdate+'_AVG.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag ne 0 then goto, jump
index=where(mark2 lt 0.)
if index(0) ne -1L then mark2(index)=0.	; zero out anticyclones
;
; compute temperature
;
      t2=0.*u2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; zonal mean zonal wind
;
      if icount eq 0L then begin
         tzm=fltarr(nr,nth)
         tzm_all=fltarr(nday,nr,nth)
         uzm=fltarr(nr,nth)
         uzm_all=fltarr(nday,nr,nth)
         markzm=fltarr(nr,nth)
         markzm_all=fltarr(nday,nr,nth)
         sdate0=sdate
      endif

      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          index=where(t2(j,*,k) ne 0.,np)
          if index(0) ne -1L then begin
          tzm(j,k)=tzm(j,k)+total(t2(j,index,k))/float(np)
          tzm_all(icount,j,k)=total(t2(j,index,k))/float(np)
          uzm(j,k)=uzm(j,k)+total(u2(j,index,k))/float(np)
          uzm_all(icount,j,k)=total(u2(j,index,k))/float(np)
;         markzm(j,k)=markzm(j,k)+max(mark2(j,*,k))
;         markzm_all(icount,j,k)=max(mark2(j,*,k))
          markzm(j,k)=markzm(j,k)+total(mark2(j,index,k))/float(np)
          markzm_all(icount,j,k)=total(mark2(j,index,k))/float(np)
          endif
      endfor
      endfor
      icount=icount+1L
goto,jump
plotit:
sdate1=sdate
tzm=tzm/float(icount)
uzm=uzm/float(icount)
markzm=markzm/float(icount)
;
; compute sigma
;
tsig=0.*tzm
usig=0.*tzm
marksig=0.*tzm
for k=0L,nth-1L do begin
for j=0L,nr-1L do begin
    dum=reform(tzm_all(*,j,k))
    index=where(dum ne 0.)
    if n_elements(index) gt 2L then begin
       result=moment(dum(index))
       tsig(j,k)=sqrt(result(1))
    endif
    dum=reform(uzm_all(*,j,k))
    index=where(dum ne 0.)
    if n_elements(index) gt 2L then begin
       result=moment(dum(index))
       usig(j,k)=sqrt(result(1))
    endif
    dum=markzm_all(*,j,k)
    index=where(dum ne 0.)
    if n_elements(index) gt 2L then begin
       result=moment(dum(index))
       marksig(j,k)=sqrt(result(1))
    endif
endfor
endfor

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_geos5_ubarth+mark_'+sdate0+'-'+sdate1+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot zonal mean zonal wind
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
index=where(uzm eq 0.)
if index(0) ne -1L then uzm(index)=-9999.
contour,uzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[min(th),max(th)],charsize=2,color=0,$
      ytitle='Potential Temperature (K)',xticks=6,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title=month(lstmn-1)+' '+syear,charthick=2
index=where(level gt 0.)
contour,uzm,alat,th,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,thick=2
index=where(level lt 0.)
contour,uzm,alat,th,levels=level(index),color=mcolor,/follow,/overplot,min_value=-9999.,c_linestyle=5,thick=2
contour,markzm,alat,th,levels=[0.5],color=0,/follow,/overplot,thick=10,c_labels=1+0*findgen(13)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(m/s)',charsize=1.5,charthick=2
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
   spawn,'convert -trim yz_geos5_ubarth+mark_'+sdate0+'-'+sdate1+'.ps -rotate -90 yz_geos5_ubarth+mark_'+sdate0+'-'+sdate1+'.jpg'
endif
endfor	; loop over months
endfor	; loop over years
end
