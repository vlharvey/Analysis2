;
; save zonal mean water deviation from multi-year mean
;
@stddat
@kgmt
@ckday
@kdate

mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
;start_year=[2007,2008,2009,2010,2011,2012,2013,2014,2015]
;start_date=[-27, -21, -24, -24, -26, -27, -34, -28,-42]
;end_date=[66, 65, 61, 61, 64, 61, 64, 80]
;nyear=n_elements(start_year)

lstmn=8
lstdy=13
lstyr=2004
ledmn=5
leddy=11
ledyr=2015
lstday=0
ledday=0
;
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
;
; longitude grid
;
dx=15.
nc=long(360./dx)+1
longrid=dx*findgen(nc)
nr=91L
latgrid=-90.+2.*findgen(nr)
dy=latgrid(1)-latgrid(0)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      if icount eq 0 then begin
         sdate_all=strarr(kday)
         dfs_all=fltarr(kday)
      endif
      sdate_all(icount)=sdate

if iday ge 60 and iday le 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(6,21,long(syr))
if iday lt 60 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr)-1L)
if iday gt 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr))
      print,sdate,' ',dfs
      dfs_all(icount)=dfs
;
; restore MLS on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v3.3_'+sdate+'.sav'
;
; apply mask
;
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      nlv=n_elements(altitude)
      if icount eq 0 then begin
         tanom=fltarr(kday,nr,nlv)
      endif
;
; compute mean around lat circle
;
    tbar=fltarr(nr,nlv)
    nbar=lonarr(nr,nlv)
    for ii=0L,n_elements(id)-1L do begin
        tmask_prof=reform(MASK(ii,*))
        good=where(tmask_prof ne -99.,ngood)
        if good(0) ne -1L then begin
           ymean=latitude(ii)
           for j=0L,nr-1L do begin
               if ymean ge latgrid(j)-dy/2. and ymean lt latgrid(j)+dy/2. then begin
                  tbar(j,good)=tbar(j,good)+mix(ii,good)
                  nbar(j,good)=nbar(j,good)+1L
               endif
           endfor
        endif
    endfor
    index=where(nbar gt 1.)
    if index(0) ne -1L then tbar(index)=tbar(index)/float(nbar(index))
;
; make multi-year mean and exclude current year
;
    tmean=0.*tbar
    nmean=0*nbar
    smon=strmid(sdate,4,2)
    sday=strmid(sdate,6,2)
    spawn,'ls '+mdir+'cat_mls_v3.3_????'+smon+sday+'.sav',cfiles
    spawn,'ls '+mdir+'h2o_mls_v3.3_????'+smon+sday+'.sav',tfiles
;
; extract years
; 
years=strarr(n_elements(cfiles))
for i=0L,n_elements(cfiles)-1L do begin
result=strsplit(cfiles(i),/extract,'/')
result2=strsplit(result(5),/extract,'.')
result3=strsplit(result2(1),/extract,'_')
years(i)=string(result3(1))
endfor
good=WHERE(STRMATCH(years, sdate) NE 1)
cfiles=cfiles(good)
tfiles=tfiles(good)

    for i=0,n_elements(cfiles)-1L do begin
        restore,cfiles(i)
        restore,tfiles(i)
        print,cfiles(i)
        index=where(mask eq -99.)
        if index(0) ne -1L then mix(index)=-99.
        for ii=0L,n_elements(id)-1L do begin
            tmask_prof=reform(MASK(ii,*))
            good=where(tmask_prof ne -99.,ngood)
            if good(0) ne -1L then begin
               ymean=latitude(ii)
               for j=0L,nr-1L do begin
                   if ymean ge latgrid(j)-dy/2. and ymean lt latgrid(j)+dy/2. then begin
                      tmean(j,good)=tmean(j,good)+mix(ii,good)
                      nmean(j,good)=nmean(j,good)+1L
                   endif
               endfor
            endif
        endfor	; loop over profiles
    endfor	; loop over past years
    index=where(nmean gt 1.)
    if index(0) ne -1L then tmean(index)=tmean(index)/float(nmean(index))
;
; anomaly from multi-year mean
;
tanom(icount,*,*)=tbar-tmean
h2oanom=tanom*1.e6
;

skipmls:
icount=icount+1L
goto,jump

saveit:
save,file='MLS_YZ_H2O_anomaly_2004-2015.sav',sdate_all,dfs_all,h2oanom,latgrid,altitude

end
