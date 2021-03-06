;
; calculate area within NH vortex and save for all levels and years
;
@stddat
@kgmt
@ckday
@kdate

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
lstmn=11
lstdy=1
lstyr=1991
iyr0=1991
ledmn=4
leddy=30
ledyr=2006
lstday=0
ledday=0
nyear=ledyr-lstyr+1L
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

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
kday=366L
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
; build MetO filename and read data
;
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      ifile=mon(imn-1)+sdy+'_'+uyr
      dum=findfile(diru+ifile+'.nc3')
      if dum(0) eq '' then goto,jump
      ncid=ncdf_open(diru+ifile+'.nc3')
      print,ifile
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      mark2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,10,mark2
      ncdf_close,ncid
;
; only on the first day
;
    if icount eq 0L then begin
       area_ave=fltarr(nyear,kday,nth)
       yyyymmdd=lonarr(nyear,kday)
       doy=lonarr(nyear,kday)
       dum=transpose(mark2(*,*,0))
       lon=0.*dum
       lat=0.*dum
       for i=0,nc-1 do lat(i,*)=alat
       for j=0,nr-1 do lon(*,j)=alon
       area=0.*lat
       nrr=91
       yeq=findgen(nrr)
       latcircle=fltarr(nrr)
       latsum=fltarr(nrr)
       hem_frac=fltarr(nrr)
       for j=0,nrr-2 do begin
           hy=re*dtr
           dx=re*cos(yeq(j)*dtr)*360.*dtr
           latcircle(j)=dx*hy	; area in each latitude circle
       endfor
       for j=0L,nrr-1 do latsum(j)=total(latcircle(j:nrr-1))
       for j=0,nrr-1 do begin
           index=where(yeq ge yeq(j))
; fraction of the hemisphere in each latitude circle
           if index(0) ne -1 then $
              hem_frac(j)=100.*total(latcircle(index))/hem_area
           if yeq(j) eq 0. then hem_frac(j)=100.
       endfor
; area in km^2
       deltax=alon(1)-alon(0)
       deltay=alat(1)-alat(0)
       for j=0,nr-1 do begin
           hy=re*deltay*dtr
           dx=re*cos(alat(j)*dtr)*deltax*dtr
           area(*,j)=dx*hy	; area of each grid point
       endfor
       icount=1L
    endif	; if first day
;
; sum area of gridpoints within Tnat isopleth
;
    myr=iyr-iyr0
    yyyymmdd(myr,iday-1L)=long(syr+smn+sdy)
    doy(myr,iday-1L)=iday
;print,long(syr+smn+sdy)
    for thlev=0,nth-1 do begin
        mark=transpose(mark2(*,*,thlev))
        index=where(lat gt 20. and mark eq 1.)
        if index(0) ne -1 then begin
           a0=total(area(index))
           area_ave(myr,iday-1L,thlev)=a0/1.e6		; millions of sqare km
;print,th(thlev),a0/1.e6
        endif
    endfor
;stop
goto,jump
saveit:
save,file='ukmo_nhvortex_area_allyears.sav',area_ave,th,yyyymmdd,doy
end
