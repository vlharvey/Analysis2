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
; OH and HO2 data for Laura Holt
;  8/7/2009
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
lstmn=12L & lstdy=13L & lstyr=9L
ledmn=6L & leddy=27L & ledyr=10L
lstday=0L & ledday=0L 
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
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
;     dum=findfile(odir+'cat_mls_'+sver+'_'+sdate+'.sav')
;     if dum(0) ne '' then goto,jump
;
; look for EOS-MLS data files for today
;
      spawn,'ls '+dir+'MLS-Aura_L2GP-OH_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',ohfiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-GPH_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',gpfiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-HO2_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',ho2files
      spawn,'ls '+dir+'MLS-Aura_L2GP-Temperature_v02-2*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'.he5',tpfiles
      result1=size(ohfiles)
      result4=size(gpfiles)
      result5=size(ho2files)
      result10=size(tpfiles)
;
; this logic will jump day if any one of the above 4 products are missing
;
      if result1(0) eq 0L or result4(0) eq 0L or result5(0) eq 0L or result10(0) eq 0L then begin
         print,'MLS data missing on '+sdate
         goto,jump
      endif
      print,'MLS data complete on '+sdate
;     print,ohfiles,gpfiles,ho2files,tpfiles
;
; read EOS-MLS data today.  Original data is in the form of structures
;
      oh=readl2gp_std(ohfiles(0),swathName='OH',variableName=variableName,precisionName=precisionName)
      gp=readl2gp_std(gpfiles(0),swathName='GPH',variableName=variableName,precisionName=precisionName)
      ho2=readl2gp_std(ho2files(0),swathName='HO2',variableName=variableName,precisionName=precisionName)
      tp=readl2gp_std(tpfiles(0),swathName='Temperature',variableName=variableName,precisionName=precisionName)
;
; number of levels varies
; tp.pressure <Expression>    FLOAT     = Array[47]
; ho2.pressure <Expression>    FLOAT     = Array[37]
; oh.pressure <Expression>    FLOAT     = Array[49]
;
      mprof=oh.ntimes
      pmls=ho2.pressure		; 37 elements (HO2)
      pmls2=tp.pressure		; 47 elements (Temp and GPH)
      pmls3=oh.pressure		; 49 elements (OH)
      mlev=n_elements(pmls)
      mlev2=n_elements(pmls2)
      mlev3=n_elements(pmls3)
;
; common levels where zflag=1
;
      zflag=fltarr(mlev3)
      for k=0L,mlev3-1L do begin
          index1=where(pmls3(k) eq pmls)
          index2=where(pmls3(k) eq pmls2)
          if index1(0) ne -1L and index2(0) ne 0L then begin
             zflag(k)=1.
;            print,k,pmls3(k),pmls2(index2(0)),pmls(index1(0))
          endif
      endfor
      kindex=where(zflag eq 1.,mlev)
      longitude=oh.longitude
      index=where(longitude lt 0. and longitude ne -999.990)
      if index(0) ne -1L then longitude(index)=longitude(index)+360.
      latitude=oh.latitude
      mtime=float(oh.time)
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
      ohmls=transpose(oh.L2GPVALUE)
      ohprecision=transpose(oh.L2GPPRECISION)
      ohstatus=oh.STATUS
      ohquality=oh.QUALITY
      ohconvergence=oh.CONVERGENCE
      gpmls=transpose(gp.L2GPVALUE)
      gpprecision=transpose(gp.L2GPPRECISION)
      gpstatus=gp.STATUS
      gpquality=gp.QUALITY
      gpconvergence=gp.CONVERGENCE
      ho2mls=transpose(ho2.L2GPVALUE)
      ho2precision=transpose(ho2.L2GPPRECISION)
      ho2status=ho2.STATUS
      ho2quality=ho2.QUALITY
      ho2convergence=ho2.CONVERGENCE
      tpmls=transpose(tp.L2GPVALUE)
      tpprecision=transpose(tp.L2GPPRECISION)
      tpstatus=tp.STATUS
      tpquality=tp.QUALITY
      tpconvergence=tp.CONVERGENCE
;
; remove non-common pressure levels
;
      ohmls=ohmls(*,kindex)
      ohprecision=ohprecision(*,kindex)
      ohstatus=ohstatus(*,kindex)
      ohquality=ohquality(*,kindex)
      ohconvergence=ohconvergence(*,kindex)
      ohmask=0.*ohmls
      gpmls=gpmls(*,kindex)
      gpprecision=gpprecision(*,kindex)
      gpstatus=gpstatus(*,kindex)
      gpquality=gpquality(*,kindex)
      gpconvergence=gpconvergence(*,kindex)
      gpmask=0.*gpmls
      ho2mls=ho2mls(*,kindex)
      ho2precision=ho2precision(*,kindex)
      ho2status=ho2status(*,kindex)
      ho2quality=ho2quality(*,kindex)
      ho2convergence=ho2convergence(*,kindex)
      ho2mask=0.*ho2mls
      tpmls=tpmls(*,kindex)
      tpprecision=tpprecision(*,kindex)
      tpstatus=tpstatus(*,kindex)
      tpquality=tpquality(*,kindex)
      tpconvergence=tpconvergence(*,kindex)
      tpmask=0.*tpmls
;
; use quality, status, and precision flags to "mask" suspect data
;
      ohbad=where(ohprecision lt 0.)
      if ohbad(0) ne -1L then ohmask(ohbad)=-99.
      gpbad=where(gpprecision lt 0.)
      ho2bad=where(ho2precision lt 0.)
      if ho2bad(0) ne -1L then ho2mask(ho2bad)=-99.
      tpbad=where(tpprecision lt 0.)
      if tpbad(0) ne -1L then tpmask(tpbad)=-99.
;
; status=0 is good, all odd values are bad
;
      ohbad=where(ohstatus mod 2 ne 0L)
      if ohbad(0) ne -1L then ohmask(ohbad,*)=-99.
      gpbad=where(gpstatus mod 2 ne 0L)
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
      ho2bad=where(ho2status mod 2 ne 0L)
      if ho2bad(0) ne -1L then ho2mask(ho2bad,*)=-99.
      tpbad=where(tpstatus mod 2 ne 0L)
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose quality limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
;     ohbad=where(ohquality lt 1.1)			; ignore
;     if ohbad(0) ne -1L then ohmask(ohbad,*)=-99.
      gpbad=where(tpquality lt 0.6)			; do not use if quality < 0.6
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
;     ho2bad=where(ho2quality lt 1.1)			; ignore
;     if ho2bad(0) ne -1L then ho2mask(ho2bad,*)=-99.
      tpbad=where(tpquality lt 0.5)			; do not use if tpquality < 0.5
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose convergence limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
      ohbad=where(ohconvergence gt 1.1)                   ; do not use if convergence > 1.1
      if ohbad(0) ne -1L then ohmask(ohbad,*)=-99.
      gpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
      ho2bad=where(ho2convergence gt 9999.)                   ; no not use if convergence > 1.1
      if ho2bad(0) ne -1L then ho2mask(ho2bad,*)=-99.
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
      ohgrid=-99.+fltarr(mprof,nz)
      ohprecgrid=-99.+fltarr(mprof,nz)
      ohmaskgrid=-99.+fltarr(mprof,nz)
      ho2grid=-99.+fltarr(mprof,nz)
      ho2precgrid=-99.+fltarr(mprof,nz)
      ho2maskgrid=-99.+fltarr(mprof,nz)
      for i=0L,mprof-1L do begin
          zprof=reform(zmls(i,*))
          ohprof=reform(ohmls(i,*))
          ohprecprof=reform(ohprecision(i,*))
          ohmaskprof=reform(ohmask(i,*))
          ho2prof=reform(ho2mls(i,*))
          ho2precprof=reform(ho2precision(i,*))
          ho2maskprof=reform(ho2mask(i,*))
          for k=0L,nz-1L do begin           
              zp=altitude(k)
              if zp lt min(zprof) then goto,jumplev
              if zp gt max(zprof) then goto,jumplev
              for kk=0L,mlev-2L do begin
                  z0=zprof(kk) & z1=zprof(kk+1L)
                  if z0 le zp and z1 ge zp then begin
                     zscale=(z1-zp)/(z1-z0)

                     ho2grid(i,k)=ho2prof(kk+1L)-zscale*(ho2prof(kk+1L)-ho2prof(kk))
                     ho2precgrid(i,k)=ho2precprof(kk+1L)-zscale*(ho2precprof(kk+1L)-ho2precprof(kk))
                     if ho2maskprof(kk) ne -99. and ho2maskprof(kk+1L) ne -99. then $
                        ho2maskgrid(i,k)=0.              ; unmask if both are good

                     ohgrid(i,k)=ohprof(kk+1L)-zscale*(ohprof(kk+1L)-ohprof(kk))
                     ohprecgrid(i,k)=ohprecprof(kk+1L)-zscale*(ohprecprof(kk+1L)-ohprecprof(kk))
                     if ohmaskprof(kk) ne -99. and ohmaskprof(kk+1L) ne -99. then $
                        ohmaskgrid(i,k)=0.              ; unmask if both are good

                      goto,jumplev
                   endif
               endfor

              jumplev:
          endfor
;
; check
;
;plot,1.e9*ohprof,zprof,yrange=[0.,130.],xrange=[0.,5.],thick=2
;oplot,1.e9*ohgrid(i,*),altitude,psym=2
;a=findgen(8)*(2*!pi/8.)
;usersym,cos(a),sin(a),/fill
;x=where(ohmaskgrid(i,*) eq -99.)
;oplot,1.e9*ohgrid(i,x),altitude(x),psym=8,color=.9*mcolor
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
         ohgrid=reform(ohgrid(index,*))
         ohprecgrid=reform(ohprecgrid(index,*))
         ohmaskgrid=reform(ohmaskgrid(index,*))
         ho2grid=reform(ho2grid(index,*))
         ho2precgrid=reform(ho2precgrid(index,*))
         ho2maskgrid=reform(ho2maskgrid(index,*))
      endif
;
; save daily catalog file containing altitude,date,fdoy,id,latitude,longitude,time
;
;     comment=strarr(3)
;     comment(0)='date is date (yyyymmdd); time is UT time (fractional hours)'
;     comment(1)='fdoy is fractional day; id is yyyymmdd+daily event #'
;     comment(2)='altitude = SOSST Altitude grid (km)'
;     date=long(sdate)
;     time=uttime
      sprof=string(format='(i4.4)',1.+findgen(mprof))
      id=syyyymmdd+'.'+sprof
;     catfile=odir+'cat_mls_'+sver+'_'+sdate+'.sav'
;     save,file=catfile,date,time,latitude,longitude,id,fdoy,altitude,comment
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
; OH
;
      mix=ohgrid
      err=-99.+0.*ohprecgrid
      index=where(mix ne -99. and ohprecgrid ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(ohprecgrid(index))/abs(mix(index))
      mask=ohmaskgrid
      save,file=odir+'oh_mls_'+sver+'_'+sdate+'.sav',altitude,id,mix,err,mask,comment
;
; HO2
;
      mix=ho2grid
      err=-99.+0.*ho2precgrid
      index=where(mix ne -99. and ho2precgrid ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(ho2precgrid(index))/abs(mix(index))
      mask=ho2maskgrid
      save,file=odir+'ho2_mls_'+sver+'_'+sdate+'.sav',altitude,id,mix,err,mask,comment
goto,jump
end
