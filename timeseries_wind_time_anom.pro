;
; subtract multi-year time mean from data stored in /harvey/UKMKO_means
;
; latitude-altitude sections of UKMO pressure data time anomalies
;
@stddat
@kgmt
@ckday
@kdate
@date2uars
@rd_ukmo
@drawvectors

; set color table
loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
!noeras=1

title2='UKMO Zonal Mean '+['Geopotential Height',$
        'Temperature','Zonal Wind','Meridional Wind','Wind Speed']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
setplot='ps'
nlg=0l
nlat=0l
nlv=0l
lstmn=1
lstdy=1
lstyr=2004
ledmn=5
leddy=1
ledyr=2004
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      UKMO Version '
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
tp_time=fltarr(nday)
sp_time=fltarr(nday)
tp0_time=fltarr(nday)
sp0_time=fltarr(nday)

icount=0L

rlon=999.0
;read,' Enter desired longitude (0-360 by 3.75, 999 for zonal mean)  ',rlon
;ilon=0
ilon=rlon/3.75

; define viewport location 
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.35
cbaryoff=0.08
cbarydel=0.01

if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='timeseries_anom.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif

; Compute initial Julian date
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

;***Read UKMO data
      file='/aura3/data/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
print,file
      rd_ukmo,file,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
              zp,tp,up,vp
      if iflg ne 0 then goto, jump

      spawn,'ls /aura2/harvey/UKMO_means/Datfiles/ukmo_'+mon(imn-1)+'*-*.sav',ifiles
      restore,ifiles(0)

; Declare plotting arrays
      x2d=fltarr(nlat,nlv)
      y2d=fltarr(nlat,nlv)
      for k=0,nlv-1 do x2d(*,k)=alat
      for j=0,nlat-1 do y2d(j,*)=p
      x2dw=fltarr(nlat-1,nlv)
      y2dw=fltarr(nlat-1,nlv)
      for k=0,nlv-1 do x2dw(*,k)=wlat
      for j=0,nlat-2 do y2dw(j,*)=p

      zzm_mean=fltarr(nlat,nlv)
      tzm_mean=fltarr(nlat,nlv)
      uzm_mean=fltarr(nlat-1,nlv)
      vzm_mean=fltarr(nlat-1,nlv)
      z2=fltarr(nlat,nlv)
      t2=fltarr(nlat,nlv)
      u2=fltarr(nlat-1,nlv)
      v2=fltarr(nlat-1,nlv)

; Zonal means
      if rlon eq 999. then begin      
      for il=0,nlg-1 do begin
          z2(0:nlat-1,0:nlv-1)=z2(0:nlat-1,0:nlv-1)+ $
                              zp(il,0:nlat-1,0:nlv-1)/1000.
          zzm_mean(0:nlat-1,0:nlv-1)=zzm_mean(0:nlat-1,0:nlv-1)+ $
                             z_mean_all(il,0:nlat-1,0:nlv-1)

          t2(0:nlat-1,0:nlv-1)=t2(0:nlat-1,0:nlv-1)+ $
                              tp(il,0:nlat-1,0:nlv-1)
          tzm_mean(0:nlat-1,0:nlv-1)=tzm_mean(0:nlat-1,0:nlv-1)+ $
                             t_mean_all(il,0:nlat-1,0:nlv-1)


          u2(0:nlat-2,0:nlv-1)=u2(0:nlat-2,0:nlv-1)+ $
                             up(il,0:nlat-2,0:nlv-1)
          uzm_mean(0:nlat-2,0:nlv-1)=uzm_mean(0:nlat-2,0:nlv-1)+ $
                             u_mean_all(il,0:nlat-2,0:nlv-1)

          v2(0:nlat-2,0:nlv-1)=v2(0:nlat-2,0:nlv-1)+ $
                             vp(il,0:nlat-2,0:nlv-1)
          vzm_mean(0:nlat-2,0:nlv-1)=vzm_mean(0:nlat-2,0:nlv-1)+ $
                             v_mean_all(il,0:nlat-2,0:nlv-1)
      endfor
      s2=sqrt( (u2/nlg)^2.0 + (v2/nlg)^2.0 )
      szm_mean=sqrt ( (uzm_mean/nlg)^2.0 + (vzm_mean/nlg)^2.0 )
      z2=(z2/nlg)	; - (zzm_mean/nlg)
      t2=(t2/nlg)	; - (tzm_mean/nlg)
      u2=(u2/nlg)	; - (uzm_mean/nlg)
      v2=(v2/nlg)	; - (vzm_mean/nlg)

      zzm_mean=(zzm_mean/nlg)
      tzm_mean=(tzm_mean/nlg)
      uzm_mean=(uzm_mean/nlg)
      vzm_mean=(vzm_mean/nlg)

      endif

; Individual longitudes
      if rlon ne 999. then begin
         z2(0:nlat-1,0:nlv-1)=zp(ilon,0:nlat-1,0:nlv-1)/1000. - zzm_mean
         t2(0:nlat-1,0:nlv-1)=tp(ilon,0:nlat-1,0:nlv-1) - tzm_mean
         u2(0:nlat-2,0:nlv-1)=up(ilon,0:nlat-2,0:nlv-1) - uzm_mean
         v2(0:nlat-2,0:nlv-1)=vp(ilon,0:nlat-2,0:nlv-1) - vzm_mean
      endif

tp_time(icount)=min(t2)
sp_time(icount)=max(s2)
tp0_time(icount)=min(tzm_mean)
sp0_time(icount)=max(szm_mean)

index=where(t2 eq min(t2))
print,'min temp ',x2d(index(0)),y2d(index(0))
index=where(s2 eq max(s2))
print,'max wind ',x2dw(index(0)),y2dw(index(0))

icount=icount+1L
goto, jump

plotit:

; Autoscale if scale values for parameter/level are not defined
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,findgen(nday),sp_time,title='Max Zonal Mean Wind Speed',xtitle='2004 Jday',yrange=[0.,120.]
oplot,findgen(nday),sp0_time,psym=1

if setplot eq 'ps' then device, /close
end
