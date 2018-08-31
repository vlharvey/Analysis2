;
; use subroutines provided by DAAC to read UARS MLS data v5 and v6 for HNO3
; save ozone, temperature, geopotential height, and pressure in SOSST format
; i.e. mix, err, mask arrays
;
@stddat
@kgmt
@ckday
@kdate
@aura2date

dy=2.5
mnr=long((180./dy) +1.0)
latbin=-90.+dy*findgen(mnr)
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
dirm='/aura6/data/MLS_data/Datfiles/'
odir='/aura6/data/MLS_data/Datfiles_SOSST/'
;
; user enters date
;
lstmn=9L & lstdy=1L & lstyr=91L
ledmn=8L & leddy=1L & ledyr=99L
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
; --- Calculate UARS day from (imn,idy,iyr) information.
;
      z = date2uars(imn,idy,iyr,uday)
      if uday lt 8L then goto,jump
      print,iyr,imn,idy,uday
      suday=string(FORMAT='(I4.4)',uday)
;
; look for 6 MLS data files on this day and jump day if any are missing
;
      dum=findfile(dirm+'MLS_L3*D'+suday+'*PROD')
      result=size(dum)
      if result(0) eq 0L or result(1) lt 6L then goto,jump
;print,dum
;stop
;
; read version 5 UARS MLS data from Goddard DAAC
; The UARS pressure array is defined as:
; P = 1000 * 10^(-i/6) where i=0,1,2,... (indices of vertical dimension)
;
      pressure=10.0^(3-findgen(43)/6)	; all pressure levels
      dum=findfile(dirm+'MLS_L3AT_SCLO_D'+suday+'*PROD')
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3at, x, h, s, FILE = infile, SWAP = swap
      s = size(x[0].data) ; s[1] = size of vertical dimension for data
      pclo = pressure(x(0).INDEX_1ST_PT:x(0).INDEX_1ST_PT+s[1]-1L)
      clodata=x.data
      cloprec=x.quality
      nlevclo=s(1)

      dum='/aura6/data/MLS_data/Datfiles/MLS_L3AT_SGPH_D'+suday+'*PROD'
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3at, x, h, s, FILE = infile, SWAP = swap
      s = size(x[0].data)
      pgp = pressure(x(0).INDEX_1ST_PT:x(0).INDEX_1ST_PT+s[1]-1L)
      gpdata=x.data
      gpprec=x.quality
      nlevgp=s(1)

      dum='/aura6/data/MLS_data/Datfiles/MLS_L3AT_SHNO3_D'+suday+'*PROD'
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3at, x, h, s, FILE = infile, SWAP = swap
      s = size(x[0].data)
      phno3 = pressure(x(0).INDEX_1ST_PT:x(0).INDEX_1ST_PT+s[1]-1L)
      hno3data=x.data
      hno3prec=x.quality
      nlevhno3=s(1)

      dum='/aura6/data/MLS_data/Datfiles/MLS_L3AT_STEMP_D'+suday+'*PROD'
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3at, x, h, s, FILE = infile, SWAP = swap
      s = size(x[0].data)
      ptemp = pressure(x(0).INDEX_1ST_PT:x(0).INDEX_1ST_PT+s[1]-1L)
      tempdata=x.data
      tempprec=x.quality
      nlevtemp=s(1)

      dum='/aura6/data/MLS_data/Datfiles/MLS_L3AT_SO3_205_D'+suday+'*PROD'
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3at, x, h, s, FILE = infile, SWAP = swap
      mtime=24.*x.time(1)/86400000.
      mlat=x.lat
      mlon=x.lon
      s = size(x[0].data)
      Po3 = pressure(x(0).INDEX_1ST_PT:x(0).INDEX_1ST_PT+s[1]-1L)
      o3data = x.data
      o3prec = x.quality
      nlevo3=s(1)
;
; mprof=number of daily profiles
;
      s = size(x.data)
      mprof=s(2)

;     dum='/aura6/data/MLS_data/Datfiles/MLS_L3AT_SO3_183_D'+suday+'*PROD'
;     infile=dum(0)
;     read_3at, x, h, s, FILE = infile, SWAP = swap
;     s = size(x[0].data)
;     Po3183 = 10.0^(3-findgen(s[1])/6)
;     o3183data = x.data
;     o3183prec = x.quality
;     nlevo3183=s(1)
;
; read quality and status flags
;
      dum='/aura6/data/MLS_data/Datfiles/MLS_L3TP_SPARAM_L3TP_D'+suday+'*PROD'
print,dum
      infile=dum(0)
if infile eq '' then goto,jump
      read_3tp, x, h, s, FILE = infile, SWAP = swap
      cloqual=x.param.QUAL_CLO
      o3qual=x.param.QUAL_O3_205
      hno3qual=x.param.QUAL_O3_205
      tempqual=x.param.QUAL_TEMP
      gpqual=x.param.QUAL_TEMP
      mmafstat=x.param.MMAF_STAT
;
; create "mask" arrays for each variable
;
      clomask=0.*clodata
      for k=0L,nlevclo-1L do begin
          levdata=reform(clodata(k,*))
          levprec=reform(cloprec(k,*))
          index=where( (mmafstat ne 'G' and mmafstat ne 'T' and mmafstat ne 't') $
                        or cloqual ne 4. or levprec lt 0.,npts)
          if index(0) ne -1L then clomask(k,index)=-99.
      endfor
      hno3mask=0.*hno3data
      for k=0L,nlevhno3-1L do begin
          levdata=reform(hno3data(k,*))
          levprec=reform(hno3prec(k,*))
          index=where( (mmafstat ne 'G' and mmafstat ne 'T' and mmafstat ne 't') $
                        or hno3qual ne 4. or levprec lt 0.,npts) 
          if index(0) ne -1L then hno3mask(k,index)=-99.
      endfor
      tempmask=0.*tempdata
      for k=0L,nlevtemp-1L do begin
          levdata=reform(tempdata(k,*))
          levprec=reform(tempprec(k,*))
          index=where( (mmafstat ne 'G' and mmafstat ne 'T' and mmafstat ne 't') $
                        or tempqual ne 4. or levprec lt 0.,npts) 
          if index(0) ne -1L then tempmask(k,index)=-99.
      endfor
      gpmask=0.*gpdata
      for k=0L,nlevgp-1L do begin
          levdata=reform(gpdata(k,*))
          levprec=reform(gpprec(k,*))
          index=where( (mmafstat ne 'G' and mmafstat ne 'T' and mmafstat ne 't') $
                        or gpqual ne 4. or levprec lt 0.,npts)
          if index(0) ne -1L then gpmask(k,index)=-99.
      endfor
      o3mask=0.*o3data
      for k=0L,nlevo3-1L do begin
          levdata=reform(o3data(k,*))
          levprec=reform(o3prec(k,*))
          index=where( (mmafstat ne 'G' and mmafstat ne 'T' and mmafstat ne 't') $
                        or o3qual ne 4. or levprec lt 0.,npts)
          if index(0) ne -1L then o3mask(k,index)=-99.
      endfor
;
; build fractional day from day of year and UT time
;
      fdoy=0.*mtime
      syyyymmdd=sdate+strarr(mprof)
      for i=0L,mprof-1L do begin
          kyr=long(strmid(syyyymmdd(i),0,4))
          kmn=long(strmid(syyyymmdd(i),4,2))
          kdy=long(strmid(syyyymmdd(i),6,2))
          z = kgmt(kmn,kdy,kyr,kday)
          fdoy(i)=float(kday)+mtime(i)/24.
      endfor
;
; save UARS MLS in pseudo-SOSST format
;
; save daily catalog file containing pressure,date,fdoy,id,latitude,longitude,time
;
      comment=strarr(3)
      comment(0)='date is date (yyyymmdd); time is UT time (fractional hours)'
      comment(1)='fdoy is fractional day; id is yyyymmdd+daily event #'
      comment(2)='pressure = UARS pressure grid (hPa)'
      date=long(sdate)
      time=mtime
      latitude=mlat
      longitude=mlon
      sprof=string(format='(i4.4)',1.+findgen(mprof))
      id=sdate+'.'+sprof
      catfile=odir+'cat_uarsmls_'+sdate+'.sav'
      save,file=catfile,date,time,latitude,longitude,id,fdoy,pressure,comment
;
; save ClO
;
      mix=clodata
      err=-99.+0.*cloprec
      index=where(mix ne -99. and cloprec ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(cloprec(index))/abs(mix(index))
      mask=clomask
      clopress=pclo
      mix=transpose(mix)
      err=transpose(err)
      mask=transpose(mask)
      save,file=odir+'clo_uarsmls_'+sdate+'.sav',clopress,id,mix,err,mask,comment
;help,'ClO ',clopress,id,mix,err,mask
;
; save geopotential height
;
      mix=gpdata
      err=-99.+0.*gpprec
      index=where(mix ne -99. and gpprec ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(gpprec(index))/abs(mix(index))
      mask=gpmask
      mix=transpose(mix)
      err=transpose(err)
      mask=transpose(mask)
      save,file=odir+'gp_uarsmls_'+sdate+'.sav',pressure,id,mix,err,mask,comment
;help,'GP ',pressure,id,mix,err,mask
;
; save nitric acid
;
      mix=hno3data
      err=-99.+0.*hno3prec
      index=where(mix ne -99. and hno3prec ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(hno3prec(index))/abs(mix(index))
      mask=hno3mask
      hno3press=phno3
      mix=transpose(mix)
      err=transpose(err)
      mask=transpose(mask)
      save,file=odir+'hno3_uarsmls_'+sdate+'.sav',hno3press,id,mix,err,mask,comment
;help,'HNO3 ',hno3press,id,mix,err,mask
;
; save temperature
;
      mix=tempdata
      err=-99.+0.*tempprec
      index=where(mix ne -99. and tempprec ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(tempprec(index))/abs(mix(index))
      mask=tempmask
      mix=transpose(mix)
      err=transpose(err)
      mask=transpose(mask)
      save,file=odir+'temp_uarsmls_'+sdate+'.sav',pressure,id,mix,err,mask,comment
;help,'Temp ',pressure,id,mix,err,mask
;
; save ozone
;
      mix=o3data
      err=-99.+0.*o3prec
      index=where(mix ne -99. and o3prec ne -99.)
      if index(0) ne -1L then $
         err(index)=100.*abs(o3prec(index))/abs(mix(index))
      mask=o3mask
      o3press=po3
      mix=transpose(mix)
      err=transpose(err)
      mask=transpose(mask)
      save,file=odir+'o3_uarsmls_'+sdate+'.sav',o3press,id,mix,err,mask,comment
;help,'O3 ',o3press,id,mix,err,mask
goto,jump
end
