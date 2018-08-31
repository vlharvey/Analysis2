;
; ECMWF version.  Inspired by Gary Thomas to compare SSWs with PMCs
;
; from zonal mean temperature and zonal wind determine
; store 2-D arrays (day vs. altitude) of dT/dy and Ubar 
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
dir='/aura7/harvey/ECMWF_data/Datfiles/ecmwf_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1L & lstdy=1L & lstyr=1979L
ledmn=8L & leddy=31L & ledyr=2002L
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1950 or lstyr gt 2002 then stop,'Year out of range '
if ledyr lt 1950 or ledyr gt 2002 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd=lonarr(nfile)
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
      date=syr+smn+sdy
      yyyymmdd(icount)=long(syr+smn+sdy)
      print,date
;
; read ECMWF pressure data on /aura7
;
; ALAT            FLOAT     = Array[73]
; ALON            FLOAT     = Array[144]
; GP              FLOAT     = Array[144, 73, 23]
; NC              LONG      =          144
; NL              LONG      =           23
; NR              LONG      =           73
; OZ              FLOAT     = Array[144, 73, 23]
; PRESS           FLOAT     = Array[23]
; PV              FLOAT     = Array[144, 73, 23]
; SH              FLOAT     = Array[144, 73, 23]
; TP              FLOAT     = Array[144, 73, 23]
; UU              FLOAT     = Array[144, 73, 23]
; VV              FLOAT     = Array[144, 73, 23]
; WW              FLOAT     = Array[144, 73, 23]
;
      ifile=dir+smn+'_'+sdy+'_'+syr+'_12Z.sav'
      dum1=findfile(ifile)
      if dum1(0) ne '' then restore,ifile
      if dum1(0) eq '' then goto,skip

      if icount eq 0L then begin
         y60n=where(alat eq 60.)
         y60s=where(alat eq -60.)
; 
; dTdy between 60 and 90 for NH and SH.
; Ubar at 60 in NH and SH
;
         nh_dtdy=fltarr(nfile,nl)
         sh_dtdy=fltarr(nfile,nl)
         nh_ubar60=fltarr(nfile,nl)
         sh_ubar60=fltarr(nfile,nl)
      endif
;
; calculate zonal mean temperature and zonal wind
;
      uzm=fltarr(nr,nl)
      tzm=fltarr(nr,nl)
      for k=0,nl-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tp(good,j,k))/float(nc)
              uzm(j,k)=total(uu(good,j,k))/float(nc)
          endfor
          if tzm(nr-1,k) ne 0. and tzm(y60n(0),k) ne 0. then nh_dtdy(icount,k)=tzm(nr-1,k)-tzm(y60n(0),k)
          if tzm(0,k) ne 0. and tzm(y60s(0),k) ne 0. then sh_dtdy(icount,k)=tzm(0,k)-tzm(y60s(0),k)
if abs(nh_dtdy(icount,k)) gt 100. then stop
if abs(sh_dtdy(icount,k)) gt 100. then stop
          if uzm(y60n(0),k) ne 0. then nh_ubar60(icount,k)=uzm(y60n(0),k)
          if uzm(y60s(0),k) ne 0. then sh_ubar60(icount,k)=uzm(y60s(0),k)
      endfor
;
; check
; 
;erase
;level=-40.+4.*findgen(21)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*mcolor/nlvls
;!type=2^2+2^3
;if icount gt 1 then begin
;contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;print,min(nh_dtdy),max(nh_dtdy),min(nh_ubar60),max(nh_ubar60)
;level=-50.+5.*findgen(21)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*mcolor/nlvls
;!type=2^2+2^3
;if icount gt 1 then begin
;if icount lt 300 then begin
;contour,nh_ubar60(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_ubar60(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_ubar60(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;if icount ge 300 then begin
;contour,nh_ubar60(icount-300:icount-1,*),findgen(300),press,/ylog,yrange=[1000.,1.],/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_ubar60(icount-300:icount-1,*),findgen(300),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_ubar60(icount-300:icount-1,*),findgen(300),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;endif
    skip:
    icount=icount+1L
goto,jump
;
; save file
;
saveit:
comment='dT/dy=Tbar90-Tbar60 and Ubar60 is zonal mean zonal wind at 60N/S'
save,file='ECMWF_dTdy_Ubar_SSW_Climo.sav',yyyymmdd,press,nh_dtdy,sh_dtdy,nh_ubar60,sh_ubar60,comment
end
