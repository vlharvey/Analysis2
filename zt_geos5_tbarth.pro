;
; GEOS-5 version
;
; plot daily Tbar in altitude and time
; use nc3 geos5 data
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/'
stimes=[$
'_AVG.V01.']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=2
lstdy=1
lstyr=2009
ledmn=2
leddy=28
ledyr=2009
lstday=0
ledday=0
nlv=201L
altitude=findgen(nlv)
;goto,plotit
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS-5 Version '
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
kday=ledday-lstday+1L
;
; Compute initial Julian date
;
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

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-5 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'+sdate+stimes(0)+'nc3'
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+ifile,nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      t2=0.*pv2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
;
; zonal mean temperature and geopotential height
;
      tbarg=fltarr(nr,nth)
      ubarg=fltarr(nr,nth)
      zbarg=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          tpts=reform(t2(j,*,k))
          upts=reform(u2(j,*,k))
          zpts=reform(z2(j,*,k))
          index=where(tpts ne 0.,nx)
          if index(0) ne -1L then begin
              tbarg(j,k)=total(tpts(index))/float(nx)
              ubarg(j,k)=total(upts(index))/float(nx)
              zbarg(j,k)=total(zpts(index))/float(nx)
          endif
      endfor
      endfor
;
; first day
;
      if icount eq 0L then begin
         sdate0=sdate
         ttz=fltarr(kday,nlv,nr)
         utz=fltarr(kday,nlv,nr)
         sdates=strarr(kday)
;        rlat=0.
;        print,alat
;        read,' Enter Desired Latitude ',rlat
;        slat=string(format='(f5.1)',rlat)
      endif
      sdates(icount)=sdate
;
; interpolate GEOS temperature to SABER height surfaces
;
tbarz=fltarr(nr,nlv)
ubarz=fltarr(nr,nlv)
for kk=0L,nlv-1L do begin
    zz=altitude(kk)
    for j=0L,nr-1L do begin
zprof=reform(zbarg(j,*))
if zz gt max(zprof) then goto,jumplev
        for k=1L,nth-1L do begin
            zup=zbarg(j,k-1) & zlw=zbarg(j,k)
            if zup ge zz and zlw le zz then begin
               zscale=(zup-zz)/(zup-zlw)
               tbarz(j,kk)=tbarg(j,k-1)+zscale*(tbarg(j,k)-tbarg(j,k-1))
               ubarz(j,kk)=ubarg(j,k-1)+zscale*(ubarg(j,k)-ubarg(j,k-1))

;print,zlw,zz,zup,zscale
;print,tbarg(j,k),tbarz(j,kk),tbarg(j,k-1)
;stop
            endif
         endfor
      endfor
jumplev:
endfor

;
; retain zonal mean temperature each day
;
      for k=0L,nlv-1L do begin
      for j=0L,nr-1L do begin
          utz(icount,*,*)=transpose(reform(ubarz))
          ttz(icount,*,*)=transpose(reform(tbarz))
      endfor
      endfor

icount=icount+1L
goto,jump

saveit:
save,file='zt_geos5_tbarth.sav',utz,ttz,alat,nr,kday,sdates,altitude
plotit:
restore,'zt_geos5_tbarth.sav'
altitude=th

for j=0L,nr-1L do begin
slat=strcompress(string(format='(f5.1)',alat(j)),/remove_all)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_geos5_tbarth_'+slat+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
;
; date labels
;
syear=strmid(sdates,0,4)
smon=strmid(sdates,4,2)
sday=strmid(sdates,6,2)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+strmid(syear(xindex),2,2)
;
; plot zonal mean zonal wind
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=160.+5.*findgen(nlvls)
utzlat=smooth(reform(ttz(*,*,j)),7)
print,slat,' ',min(utzlat),max(utzlat)
contour,utzlat,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[400.,4000.],$	;yrange=[20.,68.],$
      charsize=1.25,color=0,ytitle='Theta (K)',title='GEOS-5 Tbar at '+slat+' Latitude',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
contour,utzlat,1.+findgen(kday),altitude,levels=[220,240,260,280,300],color=0,/follow,/overplot
contour,utzlat,1.+findgen(kday),altitude,levels=[160,180,200],color=mcolor,/follow,/overplot
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_geos5_tbarth_'+slat+'.ps -rotate -90 zt_geos5_tbarth_'+slat+'.jpg'
endif

endfor	; loop over latitude
end
