;
; from zonal mean temperature and zonal wind determine
; store 2-D arrays (day vs. altitude) of dT/dy and Ubar 
; for both hemispheres and entire data record in one IDL save file
;
@rd_ukmo_nc3

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
diru='/Volumes/earth/aura3/data/UKMO_data/Datfiles/ukmo_'
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
; dTdy between 60 and 90 for NH and SH.
; Ubar at 60 in NH and SH
;
       nh_dtdy=fltarr(nfile,nth)
       sh_dtdy=fltarr(nfile,nth)
       nh_ubar60=fltarr(nfile,nth)
       sh_ubar60=fltarr(nfile,nth)
    endif
    ncdf_varget,ncid,4,p2
    ncdf_varget,ncid,6,u2
    ncdf_close,ncid
;
; Temperature=theta*(p/po)^R/cp and divide by 1000 for km
;
    t2=0.*p2
    for k=0,nth-1 do $
        t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; calculate zonal mean temperature and zonal wind
;
    uzm=fltarr(nr,nth)
    tzm=fltarr(nr,nth)
    for k=0,nth-1 do begin
        for j=0,nr-1 do begin
            good=where(t2(j,*,k) ne 0.,nn)
            if good(0) ne -1L then begin
               tzm(j,k)=total(t2(j,good,k))/float(nn)
               uzm(j,k)=total(u2(j,good,k))/float(nn)
            endif
        endfor
        if tzm(nr-1,k) ne 0. and tzm(y60n(0),k) ne 0. then nh_dtdy(n,k)=tzm(nr-1,k)-tzm(y60n(0),k)
        if tzm(0,k) ne 0. and tzm(y60s(0),k) ne 0. then sh_dtdy(n,k)=tzm(0,k)-tzm(y60s(0),k)
if abs(nh_dtdy(n,k)) gt 100. then stop
if abs(sh_dtdy(n,k)) gt 100. then stop
        if uzm(y60n(0),k) ne 0. then nh_ubar60(n,k)=uzm(y60n(0),k)
        if uzm(y60s(0),k) ne 0. then sh_ubar60(n,k)=uzm(y60s(0),k)
    endfor
;
; check
; 
erase
;level=-40.+4.*findgen(21)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*mcolor/nlvls
;!type=2^2+2^3
;if n gt 1 then begin
;contour,nh_dtdy(0:n-1,*),findgen(n),th,/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_dtdy(0:n-1,*),findgen(n),th,/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_dtdy(0:n-1,*),findgen(n),th,/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;print,min(nh_dtdy),max(nh_dtdy),min(nh_ubar60),max(nh_ubar60)
;level=-100.+10.*findgen(21)
;nlvls=n_elements(level)
;col1=1+indgen(nlvls)*mcolor/nlvls
;!type=2^2+2^3
;if n gt 1 then begin
;if n lt 300 then begin
;contour,nh_ubar60(0:n-1,*),findgen(n),th,/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_ubar60(0:n-1,*),findgen(n),th,/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_ubar60(0:n-1,*),findgen(n),th,/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;if n ge 300 then begin
;contour,nh_ubar60(n-300:n-1,*),findgen(300),th,/noeras,levels=level,c_color=col1,/cell_fill,color=0,ytitle='Theta (K)'
;index=where(level lt 0.)
;contour,nh_ubar60(n-300:n-1,*),findgen(300),th,/noeras,levels=level(index),/follow,color=0,c_linestyle=5,/overplot
;index=where(level gt 0.)
;contour,nh_ubar60(n-300:n-1,*),findgen(300),th,/noeras,levels=level(index),/follow,color=mcolor,/overplot
;endif
;endif
    jump:
endfor		; loop over days
;
; save file
;
comment='dT/dy=Tbar90-Tbar60 and Ubar60 is zonal mean zonal wind at 60N/S'
save,file='MetO_dTdy_Ubar_SSW_Climo.sav',yyyymmdd,th,nh_dtdy,sh_dtdy,nh_ubar60,sh_ubar60,comment
end
