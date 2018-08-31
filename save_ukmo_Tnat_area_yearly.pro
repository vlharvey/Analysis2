;
; for Cynthia's climo paper
; calculate area within 2004-2005 NH Tnat and save for all levels and years
; Tnat values at theta levels and on dates in Datfiles/Tnat_SC.dat
; calculate Tnat averages at each level and use that as threshold for all days
;
loadct,38
mcolor=fix(byte(!p.color))
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
;
; loop over all years and restore Tnat values 
;
iyr0=1994L
iyr1=2004L
for iyr=iyr0,iyr1 do begin
;
; restore C. Singleton Tnat save file
; AVGTNAT         FLOAT     = Array[121, 42]
; INT_AVGTNAT     FLOAT     = Array[121, 42]
; JDAY            LONG      = Array[121]
; ZO              INT       = Array[42]
;
if iyr lt 2000L then yy=iyr-1900L
if iyr ge 2000L then yy=iyr-2000L
yy1=string(format='(i2.2)',yy)
yy2=string(format='(i2.2)',yy+1L)
if yy eq 99L then yy2=string(format='(i2.2)',0L)
restore,'/aura7/shaw/MetFiles_ForLynn/Tnat_'+yy1+yy2+'_SC_METO.dat
print,yy1,' ',yy2
caldat,jday,m,d,y       ; convert Julian day to YYYY,MM,DD
sdate=strcompress(y,/remove_all)+string(format='(i2.2)',m)+string(format='(i2.2)',d)
;
; loop over days
;
kday=n_elements(d)
for iday=0L,kday-1L do begin
    idy=d(iday)
    imn=m(iday)
    kyr=y(iday)
    if kyr ge 2000 then uyr=kyr-2000
    if kyr lt 2000 then uyr=kyr-1900
    uyr=string(FORMAT='(I2.2)',uyr)
    syr=string(FORMAT='(I4.4)',kyr)
    smn=string(FORMAT='(I2.2)',imn)
    sdy=string(FORMAT='(I2.2)',idy)
    ifile=mon(imn-1)+sdy+'_'+uyr
    lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
    ncid=ncdf_open(diru+ifile+'.nc3')
    print,iday,ifile
    ncdf_diminq,ncid,0,name,nr
    ncdf_diminq,ncid,1,name,nc
    ncdf_diminq,ncid,2,name,nth
    alon=fltarr(nc)
    alat=fltarr(nr)
    th=fltarr(nth)
    p2=fltarr(nr,nc,nth)
    mark2=fltarr(nr,nc,nth)
    ncdf_varget,ncid,0,alon
    ncdf_varget,ncid,1,alat
    ncdf_varget,ncid,2,th
    ncdf_varget,ncid,4,p2
    ncdf_varget,ncid,10,mark2
    ncdf_close,ncid
;
; MetO temperature
;
    t2=0.*p2
    for k=0,nth-1 do t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
;
; do all of this only on the first day of each year
;
    if iday eq 0L then begin
;
; interpolate INT_AVGTNAT(124, 42) from zo levels to MetO theta levels (th)
;
       avgtnat_meto=fltarr(kday,nth)
       for i=0L,kday-1L do begin
       for kk=0L,nth-1L do begin
           th0=th(kk)
           for k=1L,n_elements(zo)-1L do begin
               thup=zo(k-1) & thlw=zo(k)
               if thup ge th0 and thlw le th0 and $
                  INT_AVGTNAT(i,k-1) ne -99. and INT_AVGTNAT(i,k) ne -99. then begin
                  zscale=(thup-th0)/(thup-thlw)
                  avgtnat_meto(i,kk)=INT_AVGTNAT(i,k-1)+zscale*(INT_AVGTNAT(i,k)-INT_AVGTNAT(i,k-1))
               endif
           endfor
       endfor
       endfor
;
; calculate area 
;
       area_ave=fltarr(kday,nth)
       sfile=strarr(kday)
       dum=transpose(t2(*,*,0))
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
       deltax=alon(1)-alon(0)
       deltay=alat(1)-alat(0)
       for j=0,nr-1 do begin
           hy=re*deltay*dtr
           dx=re*cos(alat(j)*dtr)*deltax*dtr
           area(*,j)=dx*hy	; area of each grid point
       endfor
    endif	; if first day
    sfile(iday)=lfile
;
; extract today's Tnat profile
;
    today=where(sdate eq syr+smn+sdy,nprof)
    if nprof eq 0L then goto,jumpday
    if nprof gt 1L then stop
    tnat_prof=reform(avgtnat_meto(today,*))
;
; check
;
;index=where(tnat_prof ne 0.)
;tmax=max(tnat_prof(index))
;tmin=min(tnat_prof(index))
;plot,tnat_prof(index),th(index),psym=2,xrange=[tmin,tmax],yrange=[300.,1000.],/ynozero
;oplot,tnat_prof(index),th(index),psym=0,thick=4
;tnat_prof_orig=reform(INT_AVGTNAT(today,*))
;index=where(tnat_prof_orig gt 0.)
;oplot,tnat_prof_orig(index),zo(index),thick=3,color=mcolor*.9
;wait,1
;
; sum area of gridpoints within Tnat isopleth
;
    for thlev=0,nth-1 do begin
        temp=transpose(t2(*,*,thlev))
        mark=transpose(mark2(*,*,thlev))
        index=where(lat gt 20. and mark eq 1. and temp le tnat_prof(thlev))	;avgtnat_meto(iday,thlev))
        if index(0) ne -1 then begin
;print,th(thlev),min(temp(index)),tnat_prof(thlev)
           a0=total(area(index))
           area_ave(iday,thlev)=a0/1.e6		; millions of sqare km
        endif
;if index(0) eq -1L then begin
;   kindex=where(lat gt 20. and mark eq 1.)
;   if kindex(0) ne -1L then print,'Not cold enough ',th(thlev),min(temp(kindex)),tnat_prof(thlev)
;endif
    endfor

    jumpday:
endfor		; loop over days
save,file='ukmo_Tnat_area_'+yy1+yy2+'.sav',area_ave,th,sfile,sdate
endfor		; loop over years
end
