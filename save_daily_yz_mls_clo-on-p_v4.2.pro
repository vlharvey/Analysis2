;
; save daily zonal mean ClO and Z on p to go along with existing save file of T, U, V, Z, O3 from gridded data
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
nlvls=18
col1=(findgen(nlvls)/float(nlvls))*mcolor
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_Grid/'

restore,'MLS_YZ_UVTO3_2004-2018.sav'
dfs_old=DFS_ALL         ; FLOAT     = Array[4276]
sdate_old=SDATE_ALL     ; STRING    = Array[4276]
tbar_old=TBAR_ALL       ; FLOAT     = Array[4276, 96, 55]
ubar_old=UBAR_ALL       ; FLOAT     = Array[4276, 96, 55]
vbar_old=VBAR_ALL       ; FLOAT     = Array[4276, 96, 55]
zbar_old=ZBAR_ALL       ; FLOAT     = Array[4276, 96, 55]
kdayold=n_elements(sdate_old)
mindate=min(sdate_old)
lstmn=strmid(mindate,4,2)
lstdy=strmid(mindate,6,2)
lstyr=strmid(mindate,0,4)
maxdate=max(sdate_old)
ledmn=strmid(maxdate,4,2)
leddy=strmid(maxdate,6,2)
ledyr=strmid(maxdate,0,4)
print,'First Day ',lstyr,lstmn,lstdy
print,'Last Day ',ledyr,ledmn,leddy

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
;
      dum=findfile(mdir+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav'
;
; declare arrays on first day
;
      nlv=n_elements(pmls)
      nlv55=n_elements(pmls2)
      nr=n_elements(lat)
      if icount eq 0 then begin
         clobar_all=fltarr(kday,nr,nlv)-9999.
         pressure=pmls
      endif
      clomean=mean(clo,dim=4,/Nan)         ; mean over both nodes
;
; compute zonal means
;
      clobar=mean(clomean,dim=1,/nan)
;
;erase
;!type=2^2+2^3
;contour,clobar,lat,pressure,nlevels=nlvls,/noerase,c_color=col1,/ylog,yrange=[200,1.],xrange=[-90,90],/cell_fill,title=sdate
;contour,clobar,lat,pressure,nlevels=nlvls,/noerase,color=mcolor,/overplot
;
; retain all daily zonal means
;
    clobar_all(icount,*,*)=clobar

skipmls:
icount=icount+1L
goto,jump

saveit:
;
; extract 37 levels from zbar_old
; common levels where zflag=1
;
zflag=fltarr(nlv55)
for k=0L,nlv55-1L do begin
    index=where(pmls2(k) eq pmls)
    if index(0) ne -1L then zflag(k)=1.
endfor
zindex=where(zflag eq 1.)
zbar_all=zbar_old(*,*,zindex)

save,file='MLS_YZ_ClO_2004-2018.sav',sdate_all,dfs_all,clobar_all,zbar_all,lat,pressure

end
