;
; plot the product of vortex area and PV as a function of altitude and time
; see http://www.aero.jussieu.fr/~sparc/SPARC2000_new/PosterSess2/SessionP2_5/Braathen/
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nc3

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=11L & lstdy=1L & lstyr=0L 
ledmn=4L & leddy=30L & ledyr=0L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.35]
xlen=0.8
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura3/data/UKMO_data/Datfiles/ppassm_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

for lstyr=1991,2005 do begin
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
ledyr=lstyr+1L
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
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
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
      if ndays gt ledday then goto, plotit

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy

      ifile=mon(imn-1)+sdy+'_'+uyr
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      if icount gt 0L then sfile(icount)=lfile
      dum1=findfile(diru+ifile+'.nc3')
      if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
      if dum1(0) eq '' then goto,skipit
      print,ifile
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      pv2=fltarr(nr,nc,nth)
      marksf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,pv2
      ncdf_varget,ncid,10,marksf2
      ncdf_close,ncid
      pv2=pv2/100.

      if icount eq 0L then begin
         pv_zt=fltarr(kday,nth)
         area_zt=fltarr(kday,nth)
         sfile=strarr(kday)
         sfile(icount)=lfile
         dum=transpose(marksf2(*,*,0))
         lon=0.*dum
         lat=0.*dum
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy	; area of each grid point
         endfor
      endif
;
; sum area of gridpoints in Arctic vortex
;
      for thlev=0,nth-1 do begin
          pv=transpose(pv2(*,*,thlev))
          mark=transpose(marksf2(*,*,thlev))
          index=where(lat gt 0. and mark gt 0.0,npts)
          if npts gt 1L then begin
             result=moment(pv(index))
             pv_zt(icount,thlev)=result(0)*1.e6

             a0=total(area(index))/1.e6      ; millions of sqare km
             area_zt(icount,thlev)=a0
          endif
      endfor

skipit:
icount=icount+1L
goto,jump
;
; plot altitude-time series of Arctic vortex area
;
plotit:
yy=strmid(sfile(0),6,2)
if long(yy) lt 90L then y1='20'+yy
if long(yy) gt 90L then y1='19'+yy
yy=strmid(sfile(kday-1),6,2)
if long(yy) lt 90L then y2='20'+yy
if long(yy) gt 90L then y2='19'+yy

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='ukmo_areaxpv_'+y1+'-'+y2+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; vortex strength indicator as product of vortex area and mean vortex pv
;
mpv_zt=0.*pv_zt
th0=500.
for k=0,nth-1 do mpv_zt(*,k)=pv_zt(*,k)*( (th(k)/th0) )^(-9./2.)

;pv_zt=smooth(pv_zt,3,/edge_truncate)
;area_zt=smooth(area_zt,3,/edge_truncate)
index=where(area_zt ne 0. and pv_zt ne 0. and mpv_zt lt 1.e8)
;pv_zt(index)=alog(pv_zt(index))
prod=0.*pv_zt
prod(index)=mpv_zt(index)*area_zt(index)

xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,kday,kday,1,1],[500.,500.,2000.,2000.,500.],min_value=0.,$
      xrange=[1,kday],yrange=[500.,2000.],/nodata,charsize=2,$
      ytitle='Theta (K)',title='MetO Analyses '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,3,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),0,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,450.
    plots,kindex(i)+1,500.,/continue,/data
    xyouts,kindex(i)+1,325.,xlab,/data,alignment=0.5,charsize=3
endfor
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(prod eq 0.)
if index(0) ne -1 then prod(index)=-9999./0.
prod=smooth(prod,3,/edge_truncate,/NaN)
level=2.*findgen(nlvls)
contour,prod,1.+findgen(kday),th,levels=level,/fill,$
        /cell_fill,/overplot,c_color=col1,min_value=0.
contour,prod,1.+findgen(kday),th,levels=level,c_color=0,$
        /follow,/overplot,min_value=0.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
      xtitle='Vortex Area x Mean MPV',charsize=2
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
   spawn,'convert -trim ukmo_areaxpv_'+y1+'-'+y2+'.ps -rotate -90 ukmo_areaxpv_'+y1+'-'+y2+'.jpg'
endif

endfor
end
