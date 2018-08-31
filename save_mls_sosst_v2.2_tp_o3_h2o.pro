;****************************************************************************************
; Convert EOS-MLS from he5 to IDL save "SOSST" format.  Interpolate to 0-120 km grid.   *
; Save catalog file each day with date, time, longitude, latitude, fdoy, ltime, lday,   *
; Save "meta" data each day in one file containing all of the original species and      *
; precision data on pressure surfaces as well as status and quality arrays.  Save       *
; interpolated species and temperature, pressure, and atmospheric density in separate   *
; daily files that contain "mix", "error", and "mask" arrays to match SO data format.	*
;											*
; Programed by: V. Lynn Harvey  3/05/07
;
; version 2.2 provisional data	3/05/07
;
; Temperature, Ozone, and Water Vapor only for AIM analysis (need GPH files to 
; interpolate to SOSST altitude grid
;
;               CU/LASP									*
;****************************************************************************************
@stddat
@kgmt
@ckday
@kdate
@readl2gp_std
@aura2date
loadct,38
mcolor=byte(!p.color)
;
; version
;
sver='v2.2'
;
; restore SOSST altitude grid
;
restore,'/aura6/data/MLS_data/Datfiles_SOSST/o3_mls_v2.2_20040921.sav
;restore,'/aura3/data/SAGE_II_data/Datfiles_SOSST/altitude.sav'
nz=n_elements(altitude)
;
; enter dates to convert MLS pressure data
;
lstmn=12L & lstdy=2L & lstyr=8L
ledmn=12L & leddy=2L & ledyr=8L
lstday=0L & ledday=0L 
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 or lstyr gt 2011 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 2011 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
dir='/aura6/data/MLS_data/Datfiles/'
odir='/aura6/data/MLS_data/Datfiles_SOSST/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
;
; loop over dates
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
      z = stddat(imn,idy,iyr,ndays)
      if ndays gt ledday then stop,' Normal termination condition '
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; check for save file
;
      dum=findfile(odir+'cat_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) ne '' then goto,jump
;
; look for EOS-MLS data files for today
;
      spawn,'ls '+dir+'MLS-Aura_L2GP-GPH_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',gpfiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-H2O_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',h2ofiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-O3_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',o3files
      spawn,'ls '+dir+'MLS-Aura_L2GP-Temperature_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',tpfiles
      result4=size(gpfiles)
      result5=size(h2ofiles)
      result9=size(o3files)
      result10=size(tpfiles)
;
; this logic will jump day if any one of the above 10 products are missing
;
      if result4(0) eq 0L or result5(0) eq 0L or result9(0) eq 0L or result10(0) eq 0L then begin
         print,'MLS data missing on '+sdate,iday
         goto,jump
      endif
      print,'MLS data complete on '+sdate
;
; read EOS-MLS data today.  Original data is in the form of structures
;
      gp=readl2gp_std(gpfiles(0),swathName='GPH',variableName=variableName,precisionName=precisionName)
      h2o=readl2gp_std(h2ofiles(0),swathName='H2O',variableName=variableName,precisionName=precisionName)
      o3=readl2gp_std(o3files(0),swathName='O3',variableName=variableName,precisionName=precisionName)
      tp=readl2gp_std(tpfiles(0),swathName='Temperature',variableName=variableName,precisionName=precisionName)
;
; extract some catalog information from MLS ozone structure... dimensions,x,y,t,p
;
      mprof=o3.ntimes
      pmls=o3.pressure		; 37 elements
      pmls2=tp.pressure		; 47 elements
      mlev=n_elements(pmls)
      mlev2=n_elements(pmls2)
;
; common levels where zflag=1
;
      zflag=fltarr(mlev2)
      for k=0L,mlev2-1L do begin
          index=where(pmls2(k) eq pmls)
          if index(0) ne -1L then zflag(k)=1.
      endfor
      kindex=where(zflag eq 1.)
      longitude=o3.longitude
      index=where(longitude lt 0. and longitude ne -999.990)
      if index(0) ne -1L then longitude(index)=longitude(index)+360.
      latitude=o3.latitude
      mtime=float(o3.time)
;
; convert MLS time from elapsed seconds to date (yyyymmddhh) 
;
      aura2date,yyyymmddhh,mtime
;
; convert yyyymmdd to UT time
;
      yyyymmdd=yyyymmddhh/100L
      uttime=mtime
      istime=1993010100L
      ehr=uttime/60./60.       ; convert time from seconds to hours
      hh2=0.*uttime
      for n=0L,mprof-1L do begin
          yy1=istime/1000000L
          if yy1 mod 4 eq 0 then mno(1)=29L
          if yy1 mod 4 ne 0 then mno(1)=28L
          mm1=istime/10000L-yy1*100L
          dd1=istime/100L-yy1*10000L-mm1*100L
          dd2=dd1+long(ehr(n))/24L
          hh1=istime-yy1*1000000L-mm1*10000L-dd1*100L
          yy2=yy1 & mm2=mm1
          while dd2 gt mno(mm2-1) do begin
                dd2=dd2-mno(mm2-1)
                mm2=mm2+1L
                if mm2 gt 12L then begin
                   mm2=mm2-12L
                   yy2=yy2+1L
                   if yy2 mod 4 eq 0 then mno(1)=29
                   if yy2 mod 4 ne 0 then mno(1)=28
                endif
          endwhile
          hh2(n)=ehr(n) mod 24
          if hh2(n) ge 24. then begin
             hh2(n)=hh2(n)-24.
             dd2=dd2+1L
             if dd2 gt mno(mm2-1L) then begin
                dd2=dd2-mno(mm2-1L)
                mm2=mm2+1L
                if mm2 gt 12L then begin
                   mm2=mm2-12L
                   yy2=yy2+1L
                endif
             endif
          endif
      endfor
      uttime=hh2
;
; build fractional day from day of year and UT time
;
      fdoy=0.*uttime
      syyyymmdd=strcompress(yyyymmdd,/remove_all)
      for i=0L,mprof-1L do begin
          kyr=long(strmid(syyyymmdd(i),0,4))
          kmn=long(strmid(syyyymmdd(i),4,2))
          kdy=long(strmid(syyyymmdd(i),6,2))
          z = kgmt(kmn,kdy,kyr,kday)
          fdoy(i)=float(kday)+uttime(i)/24.
      endfor
;
; CALCULATE LOCAL TIME AND LOCAL DAY.  THE TIME DIFFERENCE FROM UTTIME, dH, IS THE EAST 
; LONGITUDE MULTIPLIED BY 24 HOURS/360 DEGREES.  0 LONGITUDE IS UTC, 180 LONGITUDE IS 
; INTERNATIONAL DATE LINE.  THEREFORE, THE LOCAL TIME CAN BE CALCULATED AS: 
;     0 LE LON LE 180:  LTIME=UT+dH.  IF LT>24 THEN LT=LT-24 AND LOCAL DAY=DAY+1.
;   180 LT LON LT 360:  LTIME=UT-dH.  IF LT<0  THEN LT=LT+24 AND LOCAL DAY=DAY-1.
;
      LDAY=FIX(FDOY) & LTIME=UTTIME
      L=LONGITUDE
      X=WHERE(L GT 180.,NX)
      IF NX GT 0 THEN L(X)=360.-L(X)	; L is the Delta Longitude from 0, ranging from 0-180.
      X=WHERE(LONGITUDE LE 180. and LONGITUDE ne -999.990,NX)
      IF NX GT 0 THEN LTIME(X)=UTTIME(X)+L(X)*24./360.
      X=WHERE(LONGITUDE GT 180.,NX)
      IF NX GT 0 THEN LTIME(X)=UTTIME(X)-L(X)*24./360.
      X=WHERE(LTIME GT 24.,NX)
      IF NX GT 0 THEN BEGIN
         LTIME(X)=LTIME(X)-24.
         LDAY(X)=fix(FDOY(X))+1
      ENDIF
      X=WHERE(LTIME LT 0.,NX)
      IF NX GT 0 THEN BEGIN
         LTIME(X)=LTIME(X)+24.
         LDAY(X)=fix(FDOY(X))-1
      ENDIF
;
; extract individual species information from MLS structures
;
      o3mls=transpose(o3.L2GPVALUE)
      o3precision=transpose(o3.L2GPPRECISION)
      o3status=o3.STATUS
      o3quality=o3.QUALITY
      o3convergence=o3.CONVERGENCE
      o3mask=0.*o3mls
      gpmls=transpose(gp.L2GPVALUE)
      gpprecision=transpose(gp.L2GPPRECISION)
      gpstatus=gp.STATUS
      gpquality=gp.QUALITY
      gpconvergence=gp.CONVERGENCE
      gpmask=0.*gpmls
      h2omls=transpose(h2o.L2GPVALUE)
      h2oprecision=transpose(h2o.L2GPPRECISION)
      h2ostatus=h2o.STATUS
      h2oquality=h2o.QUALITY
      h2oconvergence=h2o.CONVERGENCE
      h2omask=0.*h2omls
      tpmls=transpose(tp.L2GPVALUE)
      tpprecision=transpose(tp.L2GPPRECISION)
      tpstatus=tp.STATUS
      tpquality=tp.QUALITY
      tpconvergence=tp.CONVERGENCE
      tpmask=0.*tpmls
;
; use quality, status, and precision flags to "mask" suspect data
;
      o3bad=where(o3precision lt 0.)
      if o3bad(0) ne -1L then o3mask(o3bad)=-99.
      gpbad=where(gpprecision lt 0.)
      if gpbad(0) ne -1L then gpmask(gpbad)=-99.
      h2obad=where(h2oprecision lt 0.)
      if h2obad(0) ne -1L then h2omask(h2obad)=-99.
      tpbad=where(tpprecision lt 0.)
      if tpbad(0) ne -1L then tpmask(tpbad)=-99.
;
; status=0 is good, all odd values are bad
;
      o3bad=where(o3status mod 2 ne 0L)
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(gpstatus mod 2 ne 0L)
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
      h2obad=where(h2ostatus mod 2 ne 0L)
      if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpstatus mod 2 ne 0L)
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose quality limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
      o3bad=where(o3quality lt 0.4)			; do not use if quality < 0.4
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(tpquality lt 0.6)			; do not use if quality < 0.6
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
      h2obad=where(h2oquality lt 0.9)			; do not use if h2oquality < 0.9
      if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpquality lt 0.5)			; do not use if tpquality < 0.5
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose convergence limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
      o3bad=where(o3convergence gt 1.8)                     ; do not use if convergence > 1.8
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
;     h2obad=where(h2oconvergence gt 9999.)                   ; no convergence for water
;     if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      tpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; calculate geometric altitude from geopotential height
; use geometric altitude to interpolate data from pressure levels to SOSST altitude grid
;
      ks=1.931853d-3
      ecc=0.081819
      gamma45=9.80
      rtd=double(180./!pi)
      dtr=1./rtd
      zmls=0.*gpmls
      gpkm=gpmls/1000.
      for j=0L,mprof-1L do begin
          sin2=sin( (latitude(j)*dtr)^2.0 )
          numerator=1.0+ks*sin2
          denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
          gammas=gamma45*(numerator/denominator)
          r=6378.137/(1.006803-(0.006706*sin2))
          zmls(j,*)= (r*gpkm(j,*))/ ( (gammas/gamma45)*r - gpkm(j,*) )
      endfor
      index=where(gpmls eq -999.990)
      if index(0) ne -1L then zmls(index)=-999.990
;
; interpolate from MLS data on pressure levels to to 1 km altitude grid
; (zprof begins at the ground, "altitude" array does too)
;
      approf=alog(pmls)				; interpolate in log pressure
      pgrid=-99.+fltarr(mprof,nz)
      o3grid=-99.+fltarr(mprof,nz)
      o3precgrid=-99.+fltarr(mprof,nz)
      o3maskgrid=-99.+fltarr(mprof,nz)
      h2ogrid=-99.+fltarr(mprof,nz)
      h2oprecgrid=-99.+fltarr(mprof,nz)
      h2omaskgrid=-99.+fltarr(mprof,nz)
      tpgrid=-99.+fltarr(mprof,nz)
      tpprecgrid=-99.+fltarr(mprof,nz)
      tpmaskgrid=-99.+fltarr(mprof,nz)
      for i=0L,mprof-1L do begin
          zprof=reform(zmls(i,kindex))
          zprof2=reform(zmls(i,*))
          o3prof=reform(o3mls(i,*))
          o3precprof=reform(o3precision(i,*))
          o3maskprof=reform(o3mask(i,*))
          h2oprof=reform(h2omls(i,*))
          h2oprecprof=reform(h2oprecision(i,*))
          h2omaskprof=reform(h2omask(i,*))
          tpprof=reform(tpmls(i,*))
          tpprecprof=reform(tpprecision(i,*))
          tpmaskprof=reform(tpmask(i,*))
          for k=0L,nz-1L do begin           
              zp=altitude(k)
              if zp lt min(zprof) then goto,jumplev
              if zp gt max(zprof) then goto,jumplev
              for kk=0L,mlev-2L do begin
                  z0=zprof(kk) & z1=zprof(kk+1L)
                  if z0 le zp and z1 ge zp then begin
                     zscale=(z1-zp)/(z1-z0)
                     pgrid(i,k)=approf(kk+1L)-zscale*(approf(kk+1L)-approf(kk))
                     pgrid(i,k)=exp(pgrid(i,k))

                     o3grid(i,k)=o3prof(kk+1L)-zscale*(o3prof(kk+1L)-o3prof(kk))
                     o3precgrid(i,k)=o3precprof(kk+1L)-zscale*(o3precprof(kk+1L)-o3precprof(kk))
                     if o3maskprof(kk) ne -99. and o3maskprof(kk+1L) ne -99. then $
                        o3maskgrid(i,k)=0.              ; unmask if both are good

                     goto,jumpout
                  endif
              endfor
jumpout:
;
; variables with 47 levels
;
              for kk=0L,mlev2-2L do begin
                  z0=zprof2(kk) & z1=zprof2(kk+1L)
                  if z0 le zp and z1 ge zp then begin
                     zscale=(z1-zp)/(z1-z0)

                     h2ogrid(i,k)=h2oprof(kk+1L)-zscale*(h2oprof(kk+1L)-h2oprof(kk))
                     h2oprecgrid(i,k)=h2oprecprof(kk+1L)-zscale*(h2oprecprof(kk+1L)-h2oprecprof(kk))
                     if h2omaskprof(kk) ne -99. and h2omaskprof(kk+1L) ne -99. then $
                        h2omaskgrid(i,k)=0.              ; unmask if both are good

                     tpgrid(i,k)=tpprof(kk+1L)-zscale*(tpprof(kk+1L)-tpprof(kk))
                     tpprecgrid(i,k)=tpprecprof(kk+1L)-zscale*(tpprecprof(kk+1L)-tpprecprof(kk))
                     if tpmaskprof(kk) ne -99. and tpmaskprof(kk+1L) ne -99. then $
                        tpmaskgrid(i,k)=0.              ; unmask if both are good

                     goto,jumplev
                  endif
              endfor

              jumplev:
          endfor
;
; check
;
;plot,tpprof,zprof2,yrange=[0.,130.],xrange=[180.,400.],thick=2
;oplot,tpgrid(i,*),altitude,psym=2
;stop
;plot,1.e6*o3prof,zprof,yrange=[0.,130.],xrange=[0.,15.],thick=2
;oplot,1.e6*o3grid(i,*),altitude,psym=2
;a=findgen(8)*(2*!pi/8.)
;usersym,cos(a),sin(a),/fill
;x=where(o3maskgrid(i,*) eq -99.)
;oplot,1.e6*o3grid(i,x),altitude(x),psym=8,color=.9*mcolor
;stop

          jumpprof:
      endfor
;
; remove bad profiles
;
      index=where(mtime gt 0. and yyyymmdd eq long(sdate),mprof)
      if index(0) ne -1L then begin
         yyyymmdd=yyyymmdd(index)
         syyyymmdd=syyyymmdd(index)
         uttime=uttime(index)
         longitude=longitude(index)
         latitude=latitude(index)
         fdoy=fdoy(index)
         lday=lday(index)
         ltime=ltime(index)
         pgrid=reform(pgrid(index,*))
         o3grid=reform(o3grid(index,*))
         o3precgrid=reform(o3precgrid(index,*))
         o3maskgrid=reform(o3maskgrid(index,*))
         h2ogrid=reform(h2ogrid(index,*))
         h2oprecgrid=reform(h2oprecgrid(index,*))
         h2omaskgrid=reform(h2omaskgrid(index,*))
         tpgrid=reform(tpgrid(index,*))
         tpprecgrid=reform(tpprecgrid(index,*))
         tpmaskgrid=reform(tpmaskgrid(index,*))
      endif
;
; save daily catalog file containing altitude,date,fdoy,id,latitude,longitude,time
;
      comment=strarr(3)
      comment(0)='date is date (yyyymmdd); time is UT time (fractional hours)'
      comment(1)='fdoy is fractional day; id is yyyymmdd+daily event #'
      comment(2)='altitude = SOSST Altitude grid (km)'
      date=long(sdate)
      time=uttime
      sprof=string(format='(i4.4)',1.+findgen(mprof))
      id=syyyymmdd+'.'+sprof
      catfile=odir+'cat_mls_'+sver+'_'+sdate+'.sav'
      save,file=catfile,date,time,latitude,longitude,id,fdoy,altitude,comment
;
; save daily "SOSST" files with mix,err,mask on 0-120km altitude grid
; convert precision to % error.  force error to be positive even if data and/or precision is negative
; negatives can be subsequently identified by where the mask is -99
;
      comment=strarr(4)
      comment(0)='Errors in % and are always positive.'
      comment(1)='Fill values are set to -99.'
      comment(2)='Mask: set to -99 for all fill values and bad data'
;
; Ozone
;
      mix=o3grid
      err=-99.+0.*o3precgrid
      index=where(mix ne -99. and o3precgrid ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(o3precgrid(index))/abs(mix(index))
      mask=o3maskgrid
      save,file=odir+'o3_mls_'+sver+'_'+sdate+'.sav',altitude,id,mix,err,mask,comment
;
; Water vapor
;
      mix=h2ogrid
      err=-99.+0.*h2oprecgrid
      index=where(mix ne -99. and h2oprecgrid ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(h2oprecgrid(index))/abs(mix(index))
      mask=h2omaskgrid
      save,file=odir+'h2o_mls_'+sver+'_'+sdate+'.sav',altitude,id,mix,err,mask,comment
;
; Temperature, Pressure, and Total Atmospheric density
;
      comment=strarr(7)
      comment(0)='id is yyyymmdd+daily event #'
      comment(1)='altitude = SOSST Altitude grid (km)'
      comment(2)='pressure = Pressure (hPa)'
      comment(3)='temperature = Temperature (K)'
      comment(4)='temperature_error = Temperature Error (%)'
      comment(5)='temperature_mask = Temperature Mask; bad data set to -99'
      comment(6)='dentot = Total atmospheric density (molecules cm^-3)'

      pressure=pgrid
      temperature=tpgrid
      temperature_prec=tpprecgrid
      temperature_mask=tpmaskgrid
      temperature_error=-99.+0.*temperature_prec
      index=where(TEMPERATURE_PREC ne -99. and TEMPERATURE ne -99.)
      if index(0) ne -1L then $
         temperature_error(index)=100.*abs(temperature_prec(index))/abs(temperature(index))
;
; Calculate total density from the ideal gas law: den=(P/T)*(1/k)
; For pressure in mbar and temperature in K, to get density in cm^-3, use k=1.38e-19
;
      k=1.38e-19
      kk=1./k
      dentot=(pressure/temperature)*kk
;     av=6.02e23
;     gasconst=8.3145		; alternate method from B. Gamblin
;     dentot= (100.*pressure*av)/(gasconst*temperature*1.e6)
      x=where(temperature_mask eq -99.)
      if x(0) ne -1L then dentot(x)=-99.
 
      tpdfile=odir+'tpd_mls_'+sver+'_'+sdate+'.sav'
      save,file=tpdfile,altitude,id,pressure,temperature,temperature_error,$
           temperature_mask,dentot,comment
goto,jump
end
