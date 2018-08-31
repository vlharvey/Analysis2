;
; save old marker field
;
@stddat
@kgmt
@ckday
@kdate

loadct,38
mcolor=!p.color
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
diro='/aura3/data/UKMO_data/Datfiles/mark_'
lstmn=10 & lstdy=1 & lstyr=91 & lstday=0
ledmn=12 & leddy=31 & ledyr=5 & ledday=0
;
; Ask interactive questions- get starting/ending date
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
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition'
;
; read old marker field
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      uyr=strmid(syr,2,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      dum1=findfile(diru+ifile+'.nc3')
      if dum1(0) eq '' then goto,jump
      if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
print,dum1
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nl
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nl)
      mark2=fltarr(nr,nc,nl)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,10,mark2
      ncdf_close,ncid
;
; save marker to separate file
;
      ofile=diro+ifile+'.nc'
      ncid = ncdf_create(ofile,/CLOBBER)
      latdimid=ncdf_dimdef(ncid, 'number_of_latitudes' , nr)
      londimid=ncdf_dimdef(ncid, 'number_of_longitudes', nc)
      levdimid=ncdf_dimdef(ncid, 'number_of_levels'    , nl)
      lonsid = ncdf_vardef(ncid, 'longitudes',  londimid)
      latsid = ncdf_vardef(ncid, 'latitudes' ,  latdimid)
      levsid = ncdf_vardef(ncid, 'th_levels' ,  levdimid)
      mksid  = ncdf_vardef(ncid, 'mark'       , [latdimid,londimid,levdimid])
      ncdf_control,ncid,/ENDEF
      ncdf_varput, ncid, lonsid, alon , COUNT=[nc]
      ncdf_varput, ncid, latsid, alat , COUNT=[nr]
      ncdf_varput, ncid, levsid, th   , COUNT=[nl]
      ncdf_varput, ncid, mksid , mark2, COUNT=[nr,nc,nl]
      ncdf_close,ncid

goto,jump
end
