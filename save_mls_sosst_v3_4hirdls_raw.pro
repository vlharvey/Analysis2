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
; version 3.3 data 9/21/10
;
; Temperature - only use: Even status, +ve precision, quality > 0.65, convergence < 1.2
; GPH: As for temperature
; O3 - only use: Even status, +ve precision, quality > 0.6, convergence < 1.18
; HNO3 - only use: Even status, +ve precision, quality > 0.5, convergence < 1.4
; O3, reject any profile where any value from 1000-68hPa is less than -0.15 ppmv, or where 316hPa is less than -0.30 ppmv
; HNO3, reject any profile where 316hPa is less than -2 ppbv, or where 316 or 215 are less than -1.6 ppbv
;
; HIRDLS variables: T, O3, HNO3, GPH
; O3, T, GPH have 55 levels. HNO3 has 37 levels.
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
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
;
; version
;
sver='v2.2'
sver='v3.3'
;
; restore SOSST altitude grid
;
restore,'/aura6/data/MLS_data/Datfiles_SOSST/o3_mls_v2.2_20040921.sav
;restore,'/aura3/data/SAGE_II_data/Datfiles_SOSST/altitude.sav'
nz=n_elements(altitude)
;
; enter dates to convert MLS pressure data
;
lstmn=9L & lstdy=1L & lstyr=4L
ledmn=9L & leddy=17L & ledyr=10L
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
;     dum=findfile(odir+'cat_mls_'+sver+'_'+sdate+'.sav')
;     if dum(0) ne '' then goto,jump
;
; look for EOS-MLS data files for today
;
;     spawn,'ls '+dir+'MLS-Aura_L2GP-BrO_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',brofiles
;     spawn,'ls '+dir+'MLS-Aura_L2GP-CO_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',cofiles
;     spawn,'ls '+dir+'MLS-Aura_L2GP-ClO_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',clofiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-GPH_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',gpfiles
;     spawn,'ls '+dir+'MLS-Aura_L2GP-H2O_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',h2ofiles
;     spawn,'ls '+dir+'MLS-Aura_L2GP-HCl_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',hclfiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-HNO3_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',hno3files
;     spawn,'ls '+dir+'MLS-Aura_L2GP-N2O_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',n2ofiles
      spawn,'ls '+dir+'MLS-Aura_L2GP-O3_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',o3files
      spawn,'ls '+dir+'MLS-Aura_L2GP-Temperature_v03-3*'+string(FORMAT='(i4.4,a1,i3.3)',iyr,'d',iday)+'*.he5',tpfiles
;     result1=size(brofiles)
;     result2=size(cofiles)
;     result3=size(clofiles)
      result4=size(gpfiles)
;     result5=size(h2ofiles)
;     result6=size(hclfiles)
      result7=size(hno3files)
;     result8=size(n2ofiles)
      result9=size(o3files)
      result10=size(tpfiles)
;
; this logic will jump day if any one of the above 10 products are missing
;
      if result4(0) eq 0L or $
         result7(0) eq 0L or $
         result9(0) eq 0L or result10(0) eq 0L then begin
         print,'MLS data missing on '+sdate
         goto,jump
      endif
      print,'MLS data complete on '+sdate
;     print,brofiles,cofiles,clofiles,gpfiles,h2ofiles,hclfiles,hno3files,n2ofiles,o3files,tpfiles
;
; read EOS-MLS data today.  Original data is in the form of structures
;
;     bro=readl2gp_std(brofiles(0),swathName='BrO',variableName=variableName,precisionName=precisionName)
;     co=readl2gp_std(cofiles(0),swathName='CO',variableName=variableName,precisionName=precisionName)
;     clo=readl2gp_std(clofiles(0),swathName='ClO',variableName=variableName,precisionName=precisionName)
      gp=readl2gp_std(gpfiles(0),swathName='GPH',variableName=variableName,precisionName=precisionName)
;     h2o=readl2gp_std(h2ofiles(0),swathName='H2O',variableName=variableName,precisionName=precisionName)
;     hcl=readl2gp_std(hclfiles(0),swathName='HCl',variableName=variableName,precisionName=precisionName)
      hno3=readl2gp_std(hno3files(0),swathName='HNO3',variableName=variableName,precisionName=precisionName)
;     n2o=readl2gp_std(n2ofiles(0),swathName='N2O',variableName=variableName,precisionName=precisionName)
      o3=readl2gp_std(o3files(0),swathName='O3',variableName=variableName,precisionName=precisionName)
      tp=readl2gp_std(tpfiles(0),swathName='Temperature',variableName=variableName,precisionName=precisionName)
;
; extract some catalog information from MLS ozone structure... dimensions,x,y,t,p
;
      mprof=o3.ntimes
      pmls=o3.pressure			; 55 elements
      pmls2=hno3.pressure		; 37 elements
      mlev=n_elements(pmls)
      mlev2=n_elements(pmls2)
;
; common levels where zflag=1
;
      zflag=fltarr(mlev)
      for k=0L,mlev-1L do begin		; loop over pressure array with most levels
          index=where(pmls(k) eq pmls2)
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
;     bromls=transpose(bro.L2GPVALUE)
;     broprecision=transpose(bro.L2GPPRECISION)
;     brostatus=bro.STATUS
;     broquality=bro.QUALITY
;     broconvergence=bro.CONVERGENCE
;     bromask=0.*bromls
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
;     comls=transpose(co.L2GPVALUE)
;     coprecision=transpose(co.L2GPPRECISION)
;     costatus=co.STATUS
;     coquality=co.QUALITY
;     coconvergence=co.CONVERGENCE
;     comask=0.*comls
;     hclmls=transpose(hcl.L2GPVALUE)
;     hclprecision=transpose(hcl.L2GPPRECISION)
;     hclstatus=hcl.STATUS
;     hclquality=hcl.QUALITY
;     hclconvergence=hcl.CONVERGENCE
;     hclmask=0.*hclmls
;     clomls=transpose(clo.L2GPVALUE)
;     cloprecision=transpose(clo.L2GPPRECISION)
;     clostatus=clo.STATUS
;     cloquality=clo.QUALITY
;     cloconvergence=clo.CONVERGENCE
;     clomask=0.*clomls
;     h2omls=transpose(h2o.L2GPVALUE)
;     h2oprecision=transpose(h2o.L2GPPRECISION)
;     h2ostatus=h2o.STATUS
;     h2oquality=h2o.QUALITY
;     h2oconvergence=h2o.CONVERGENCE
;     h2omask=0.*h2omls
      hno3mls=transpose(hno3.L2GPVALUE)
      hno3precision=transpose(hno3.L2GPPRECISION)
      hno3status=hno3.STATUS
      hno3quality=hno3.QUALITY
      hno3convergence=hno3.CONVERGENCE
      hno3mask=0.*hno3mls
;     n2omls=transpose(n2o.L2GPVALUE)
;     n2oprecision=transpose(n2o.L2GPPRECISION)
;     n2ostatus=n2o.STATUS
;     n2oquality=n2o.QUALITY
;     n2oconvergence=n2o.CONVERGENCE
;     n2omask=0.*n2omls
      tpmls=transpose(tp.L2GPVALUE)
      tpprecision=transpose(tp.L2GPPRECISION)
      tpstatus=tp.STATUS
      tpquality=tp.QUALITY
      tpconvergence=tp.CONVERGENCE
      tpmask=0.*tpmls
;
; use quality, status, and precision flags to "mask" suspect data
;
;     brobad=where(broprecision lt 0.)
;     if brobad(0) ne -1L then bromask(brobad)=-99.
      o3bad=where(o3precision lt 0.)
      if o3bad(0) ne -1L then o3mask(o3bad)=-99.
      gpbad=where(gpprecision lt 0.)
      if gpbad(0) ne -1L then gpmask(gpbad)=-99.
;     cobad=where(coprecision lt 0.)
;     if cobad(0) ne -1L then comask(cobad)=-99.
;     hclbad=where(hclprecision lt 0.)
;     if hclbad(0) ne -1L then hclmask(hclbad)=-99.
;     clobad=where(cloprecision lt 0.)
;     if clobad(0) ne -1L then clomask(clobad)=-99.
;     h2obad=where(h2oprecision lt 0.)
;     if h2obad(0) ne -1L then h2omask(h2obad)=-99.
      hno3bad=where(hno3precision lt 0.)
      if hno3bad(0) ne -1L then hno3mask(hno3bad)=-99.
;     n2obad=where(n2oprecision lt 0.)
;     if n2obad(0) ne -1L then n2omask(n2obad)=-99.
      tpbad=where(tpprecision lt 0.)
      if tpbad(0) ne -1L then tpmask(tpbad)=-99.
;
; status=0 is good, all odd values are bad
;
;     brobad=where(brostatus mod 2 ne 0L)
;     if brobad(0) ne -1L then bromask(brobad,*)=-99.
      o3bad=where(o3status mod 2 ne 0L)
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(gpstatus mod 2 ne 0L)
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
;     cobad=where(costatus mod 2 ne 0L)
;     if cobad(0) ne -1L then comask(cobad,*)=-99.
;     hclbad=where(hclstatus mod 2 ne 0L)
;     if hclbad(0) ne -1L then hclmask(hclbad,*)=-99.
;     clobad=where(clostatus mod 2 ne 0L)
;     if clobad(0) ne -1L then clomask(clobad,*)=-99.
;     h2obad=where(h2ostatus mod 2 ne 0L)
;     if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      hno3bad=where(hno3status mod 2 ne 0L)
      if hno3bad(0) ne -1L then hno3mask(hno3bad,*)=-99.
;     n2obad=where(n2ostatus mod 2 ne 0L)
;     if n2obad(0) ne -1L then n2omask(n2obad,*)=-99.
      tpbad=where(tpstatus mod 2 ne 0L)
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose quality limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
;     brobad=where(broquality lt 1.2)			; do not use if quality < 1.2
;     if brobad(0) ne -1L then bromask(brobad,*)=-99.
      o3bad=where(o3quality lt 0.6)			; do not use if quality < 0.6	; v3 value
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(tpquality lt 0.65)			; do not use if quality < 0.65	; v3 value
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
;     cobad=where(coquality lt 0.2)			; do not use if quality < 0.2
;     if cobad(0) ne -1L then comask(cobad,*)=-99.
;     hclbad=where(hclquality lt 1.)			; do not use if quality < 1.
;     if hclbad(0) ne -1L then hclmask(hclbad,*)=-99.
;     clobad=where(cloquality lt 0.8)			; do not use if quality < 0.8
;     if clobad(0) ne -1L then clomask(clobad,*)=-99.
;     h2obad=where(h2oquality lt 0.9)			; do not use if h2oquality < 0.9
;     if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      hno3bad=where(hno3quality lt 0.5)			; do not use if quality < 0.5	; v3 value
      if hno3bad(0) ne -1L then hno3mask(hno3bad,*)=-99.
;     n2obad=where(n2oquality lt 0.5)			; do not use if quality < 0.5
;     if n2obad(0) ne -1L then n2omask(n2obad,*)=-99.
      tpbad=where(tpquality lt 0.65)			; do not use if tpquality < 0.65	; v3 value
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; impose convergence limits for v2.2 (see http://avdc.gsfc.nasa.gov/Data/Aura/MLS/MLS_v2.2_provisional.html)
;
;     brobad=where(broconvergence gt 1.5)                   ; do not use if convergence > 1.5
;     if brobad(0) ne -1L then bromask(brobad,*)=-99.
      o3bad=where(o3convergence gt 1.18)                     ; do not use if convergence > 1.18		; v3
      if o3bad(0) ne -1L then o3mask(o3bad,*)=-99.
      gpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2		; v3
      if gpbad(0) ne -1L then gpmask(gpbad,*)=-99.
;     cobad=where(coconvergence gt 1.8)                     ; do not use if convergence > 1.8
;     if cobad(0) ne -1L then comask(cobad,*)=-99.
;     hclbad=where(hclconvergence gt 1.5)                    ; do not use if convergence > 1.5
;     if hclbad(0) ne -1L then hclmask(hclbad,*)=-99.
;     clobad=where(cloconvergence gt 1.5)                   ; do not use if convergence > 1.5
;     if clobad(0) ne -1L then clomask(clobad,*)=-99.
;     h2obad=where(h2oconvergence gt 9999.)                   ; no convergence for water
;     if h2obad(0) ne -1L then h2omask(h2obad,*)=-99.
      hno3bad=where(hno3convergence gt 1.4)                 ; do not use if convergence > 1.4		; v3
      if hno3bad(0) ne -1L then hno3mask(hno3bad,*)=-99.
;     n2obad=where(n2oconvergence gt 1.55)                   ; do not use if convergence > 1.55
;     if n2obad(0) ne -1L then n2omask(n2obad,*)=-99.
      tpbad=where(tpconvergence gt 1.2)                     ; do not use if convergence > 1.2		; v3
      if tpbad(0) ne -1L then tpmask(tpbad,*)=-99.
;
; For O3, reject any profile where any value from 1000-68hPa is less than -0.15 ppmv, or where 316hPa is less than -0.30 ppmv
; For HNO3, reject any profile where 316hPa is less than -2 ppbv, or where 316 or 215 are less than -1.6 ppbv
;
      for iprof=0L,mprof-1L do begin
          o3prof=reform(o3mls(iprof,*))
          bad=where((pmls le 1000. and pmls ge 68. and o3prof lt -0.15) or (long(pmls) eq 316 and o3prof lt -0.3))
;if bad(0) ne -1L then stop,'check bad o3'
          if bad(0) ne -1L then o3mask(iprof,*)=-99.	; reject entire profile

          hno3prof=reform(hno3mls(iprof,*))
          bad=where((long(pmls2) eq 316 and hno3prof lt -2.) or (long(pmls2) eq 316 and hno3prof lt -1.6) or (long(pmls2) eq 215 and hno3prof lt -1.6))
;if bad(0) ne -1L then stop,'check bad hno3'
          if bad(0) ne -1L then hno3mask(iprof,*)=-99.
      endfor
;
; save raw pressure and meta data for long term storage.  save all species into 1 daily file
;
      rawfile=odir+'raw_mls_'+sver+'_'+sdate+'.sav'
;     save,file=rawfile,pmls,pmls2,mtime,bromls,broprecision,bromask,o3mls,o3precision,o3mask,gpmls,gpprecision,gpmask,$
;          comls,coprecision,comask,hclmls,hclprecision,hclmask,clomls,cloprecision,clomask,h2omls,h2oprecision,h2omask,$
;          hno3mls,hno3precision,hno3mask,n2omls,n2oprecision,n2omask,tpmls,tpprecision,tpmask
; 
      save,file=rawfile,pmls,pmls2,mtime,o3mls,o3precision,o3mask,gpmls,gpprecision,gpmask,$
           hno3mls,hno3precision,hno3mask,tpmls,tpprecision,tpmask
;
goto,jump
end
