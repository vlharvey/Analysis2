;
; from zonal mean temperature and zonal wind determine
; the dates of major and minor stratospheric warming
; events during the UARS period.
;
; store ssw flag arrays for both hemispheres and entire data record
; in one IDL save file
;
@rd_ukmo_nc3

diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
ifile='                             '
close,1
openr,1,'ssw_finder.fil'
nfile=0L
readf,1,nfile
yyyymmdd=lonarr(nfile)
for n=0,nfile-1 do begin
    readf,1,ifile
    uyr=strmid(ifile,7,2)
    if strmid(uyr,0,1) eq '9' then iyr='19'+uyr
    if strmid(uyr,0,1) ne '9' then iyr='20'+uyr
    tmp=strmid(ifile,0,4)
    index=where(mon eq tmp)
    imn=index(0)+1L
    idy=strmid(ifile,4,2)
    syr=string(FORMAT='(I4)',iyr)
    smn=string(FORMAT='(I2.2)',imn)
    sdy=string(FORMAT='(I2.2)',idy)
    yyyymmdd(n)=long(syr+smn+sdy)
    print,ifile,long(syr+smn+sdy)
    dum1=findfile(diru+ifile+'.nc3')
    if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
    if dum1(0) eq '' then goto,jump
    if n eq 0 then begin
       nr=0L
       nc=0L
       nth=0L
       ncdf_diminq,ncid,0,name,nr
       ncdf_diminq,ncid,1,name,nc
       ncdf_diminq,ncid,2,name,nth
       alon=fltarr(nc)
       alat=fltarr(nr)
       th=fltarr(nth)
       ncdf_varget,ncid,0,alon
       ncdf_varget,ncid,1,alat
       ncdf_varget,ncid,2,th
       u2=fltarr(nr,nc,nth)
       p2=fltarr(nr,nc,nth)
       y60n=where(alat eq 61.25)
       y60s=where(alat eq -61.25)
; 
; set to 1 if minor SSW. set to 2 if major SSW
;
       nh_ssw_flag=fltarr(nfile,nth)
       sh_ssw_flag=fltarr(nfile,nth)
    endif
    ncdf_varget,ncid,4,p2
    ncdf_varget,ncid,6,u2
    ncdf_close,ncid

; Temperature=theta*(p/po)^R/cp and divide by 1000 for km
    t2=0.*p2
    for k=0,nth-1 do $
        t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; calculate zonal mean temperature and zonal wind
;
    uzm=fltarr(nr,nth)
    tzm=fltarr(nr,nth)
    tnp=fltarr(nth)
    unp=fltarr(nth)
    tsp=fltarr(nth)
    usp=fltarr(nth)
    t60n=fltarr(nth)
    u60n=fltarr(nth)
    t60s=fltarr(nth)
    u60s=fltarr(nth)
    for k=0,nth-1 do begin
        for j=0,nr-1 do begin
            tzm(j,k)=total(t2(j,*,k))/float(nc)
            uzm(j,k)=total(u2(j,*,k))/float(nc)
        endfor
        tnp(k)=tzm(nr-1,k)
        tsp(k)=tzm(0,k)
        t60n(k)=total(tzm(y60n,k))
        t60s(k)=total(tzm(y60s,k))
        un=total(uzm(y60n(0):nr-1,k))
        us=total(uzm(0:y60s(0)-1,k))
;
; set to 1 if minor SSW. set to 2 if major SSW
;
        if tnp(k) gt t60n(k) and un gt 0. then nh_ssw_flag(n,k)=1.
        if tnp(k) gt t60n(k) and un lt 0. then nh_ssw_flag(n,k)=2.
        if tsp(k) gt t60s(k) and us gt 0. then sh_ssw_flag(n,k)=1.
        if tsp(k) gt t60s(k) and us lt 0. then sh_ssw_flag(n,k)=2.
    endfor
if max(nh_ssw_flag(n,*)) eq 1. then print,'minor'
if max(nh_ssw_flag(n,*)) eq 2. then print,'major'
    jump:
endfor		; loop over days
;
; save file
;
comment='flag set to 1 if minor SSW. set to 2 if major SSW'
save,file='MetO_SSW_Climo.sav',yyyymmdd,th,nh_ssw_flag,sh_ssw_flag
end
