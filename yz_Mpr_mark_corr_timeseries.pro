;
; timeseries of MSFpr and ZM anticyclones and YZ correlation
; between 400 and 1000 K and -60 to 60N by 10deg
;
loadct,38
device,decompose=0
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.15,0.55,0.35]
yorig=[0.6,0.6,0.15]
xlen=0.3
ylen=0.3
cbaryoff=0.08
cbarydel=0.01
set_plot,'x'
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
pi2 = 6.2831853071796
dtr=pi2/360.
re=6.37E3
dir='/aura2/harvey/Hovmoller/Datfiles/'
month=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
nmonth=n_elements(month)
ifiles=[$
'ukmo_jan_01_92-dec_31_92_nc3_xt.sav',$
'ukmo_jan_01_93-dec_31_93_nc3_xt.sav',$
'ukmo_jan_01_94-dec_31_94_nc3_xt.sav',$
'ukmo_jan_01_95-dec_31_95_nc3_xt.sav',$
'ukmo_jan_01_96-dec_31_96_nc3_xt.sav',$
'ukmo_jan_01_97-dec_31_97_nc3_xt.sav',$
'ukmo_jan_01_98-dec_31_98_nc3_xt.sav',$
'ukmo_jan_01_99-dec_31_99_nc3_xt.sav',$
'ukmo_jan_01_00-dec_31_00_nc3_xt.sav',$
'ukmo_jan_01_01-dec_31_01_nc3_xt.sav',$
'ukmo_jan_01_02-dec_31_02_nc3_xt.sav',$
'ukmo_jan_01_03-dec_31_03_nc3_xt.sav',$
'ukmo_jan_01_04-dec_31_04_nc3_xt.sav']
nyear=n_elements(ifiles)
nday=365L
for n=0,nyear-1L do begin 
;
; contents of save file***
; ALON            FLOAT     = Array[96]
; MARKXT          FLOAT     = Array[96, 365, 13, 5]
; MSFXT           FLOAT     = Array[96, 365, 13, 5]
; SFILE           STRING    = Array[365]
; TH2             FLOAT     = Array[5]
; TMPXT           FLOAT     = Array[96, 365, 13, 5]
; UXT             FLOAT     = Array[96, 365, 13, 5]
; YMIDS           FLOAT     = Array[13]
;
    restore,dir+ifiles(n)
    print,'restored ',ifiles(n)
    nth=n_elements(th2)
    nr=n_elements(YMIDS)
    nc=n_elements(alon)
    if n eq 0L then begin
markyzmean=fltarr(nday,nr,nth)
msfyzmean=fltarr(nday,nr,nth)
tmpyzmean=fltarr(nday,nr,nth)
markyzall=fltarr(nyear,nday,nr,nth)
msfyzall=fltarr(nyear,nday,nr,nth)
tmpyzall=fltarr(nyear,nday,nr,nth)

       sfileall=strarr(nday*nyear)
       rlat=0.
;      print,ymids
;      read,'Enter desired latitude ',rlat
       index=where(ymids eq rlat)
       ilat=index(0)
       ralt=400.
;      print,th2
;      read,'Enter desired altitude ',ralt
       index=where(th2 eq ralt)
       ialt=index(0)
ilat0=ilat
ialt0=ialt
    endif
    dindex=where(strmid(sfile,0,6) ne 'feb_29',nfile)
    sfileall(n*nday:((n+1L)*nday)-1L)=sfile(dindex)
;
; MARKXT          FLOAT     = Array[96, 365, 13, 5]
;
for j=0L,nr-1L do begin
for k=0L,nth-1L do begin
    for iday=0,nday-1L do begin
        idy=dindex(iday)
        if idy gt 58L then idy=idy-1L	; after 28 Feb subtract 1 day
;
; retain zonal means for all years
;
        msfyzall(n,idy,j,k)=total(msfxt(*,idy,j,k))/float(nc)
        tmpyzall(n,idy,j,k)=total(tmpxt(*,idy,j,k))/float(nc)
;
; zonal mean time (annual) mean
;
        msfyzmean(idy,j,k)=msfyzmean(idy,j,k)+total(msfxt(*,idy,j,k))/float(nc)
        tmpyzmean(idy,j,k)=tmpyzmean(idy,j,k)+total(tmpxt(*,idy,j,k))/float(nc)

        index=where(markxt(*,idy,j,k) lt 0.,nin)
        if index(0) ne -1L then begin
           markyzall(n,idy,j,k)=float(nin)
           markyzmean(idy,j,k)=markyzmean(idy,j,k)+float(nin)
        endif
    endfor
endfor
endfor

endfor	; loop over year
;
; time means in Temperature and MSF
;
markyzmean=markyzmean/float(nyear)
msfyzmean=msfyzmean/float(nyear)
tmpyzmean=tmpyzmean/float(nyear)
;
; deviation from time means in Temperature and MSF
;
msfpryz=fltarr(nyear,nday,nr,nth)
tpryz=fltarr(nyear,nday,nr,nth)

for n=0,nyear-1L do begin
for j=0,nday-1L do begin

for klat=0L,nr-1L do begin
for kth=0L,nth-1L do begin
    msfpryz(n,j,klat,kth)=msfyzall(n,j,klat,kth)-msfyzmean(n,klat,kth)
    tpryz(n,j,klat,kth)=tmpyzall(n,j,klat,kth)-tmpyzmean(n,klat,kth)
endfor
endfor

endfor
endfor
;
; reform indices from (nyear,nday) to (nyear*nday)
;
tpryz2d=fltarr(nday*nyear,nr,nth)
msfpryz2d=fltarr(nday*nyear,nr,nth)
markyz2d=fltarr(nday*nyear,nr,nth)
for j=0L,nr-1L do begin
for k=0L,nth-1L do begin
for n=0L,nyear-1L do begin
    tpryz2d(n*nday:((n+1L)*nday)-1L,j,k)=tpryz(n,*,j,k)
    msfpryz2d(n*nday:((n+1L)*nday)-1L,j,k)=msfpryz(n,*,j,k)
    markyz2d(n*nday:((n+1L)*nday)-1L,j,k)=markyzall(n,*,j,k)
endfor
endfor
endfor

for imon=0L,nmonth-1L do begin
imm0=imon+1
smm0=month(imon)
;
; loop over each latitude and altitude bin
;
yzcorr=fltarr(nr,nth)-99.
for ilat=0L,nr-1L do begin
for ialt=0L,nth-1L do begin
    rlat=ymids(ilat)
    ralt=th2(ialt)
;
; plot Hov at ilat and ialt
;
sth=strcompress(string(fix(th2(ialt))),/remove_all)
if ymids(ilat) eq 0. then slat='Eq'
if ymids(ilat) gt 0. then slat=strcompress(string(ymids(ilat)),/remove_all)+'N'
if ymids(ilat) lt 0. then slat=strcompress(string(ymids(ilat)),/remove_all)+'S'
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_Mprmark_'+smm0+'_'+sth+'K_'+slat+'_uars+mark_corr.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
!type=2^2+2^3
set_viewport,xmn,xmx,ymn,ymx
mpr1d=reform(MSFPRYZ2D(*,ilat,ialt))
index=where(mpr1d eq 0. or abs(mpr1d) gt 10000.)
if index(0) ne -1L then mpr1d(index)=-9999./0.
mpr1d=smooth(mpr1d,30,/NaN,/edge_truncate)		; monthly smoother
imin=-5000.
imax=5000.
;print,sth,'  ',slat,'  ',imin,imax
if abs(imin) eq 0. and abs(imax) eq 0. then goto,jump2
nlev=11
col1=1.0+(findgen(nlev)/nlev)*mcolor
cint=(imax-imin)/(nlev-1.)
level=imin+cint*findgen(nlev)
index=where(strmid(sfileall,0,6) eq 'jan_01',nytick)
plot,findgen(nday*nyear),mpr1d,xrange=[0,nfile*nyear-1L],psym=0,yrange=[-5000.,5000.],$
     xticks=nytick-1,xtickv=index,xtickname=strmid(sfileall(index),7,2),/noeras,title='ZM MSF-MSFavg'
;
; monthly means during each year
;
mark1d=abs(reform(markyz2d(*,ilat,ialt)))
smonth=month(imm0-1)
xyouts,.4,.95,smonth+' MetO '+sth+' K '+slat,/normal,charsize=2.5
smm=strmid(sfileall,0,3)
xindex=where(smm eq smm0 and mark1d ne 0. and finite(mpr1d) eq 1L,nx)
if xindex(0) eq -1L then goto,jump2
mprdata=mpr1d(xindex)
markdata=mark1d(xindex)
date=sfileall(xindex)
syear=strmid(date,7,2)
lyear=0L*lonarr(nx)
index=where(strmid(syear,0,1) eq '9')
if index(0) ne -1L then lyear(index)=1900L+long(syear(index))
index=where(strmid(syear,0,1) eq '0')
if index(0) ne -1L then lyear(index)=2000L+long(syear(index))
nx=max(lyear)-min(lyear)+1
lyear2=min(lyear)+indgen(nx)

mprmon=fltarr(nx)
for iyear=min(lyear2),max(lyear2) do begin
    index=where(lyear eq iyear)
    if index(0) ne -1L then begin
       xx=xindex(index(0))
       o3val=mprdata(index(0))
       if n_elements(index) gt 1L then begin
          o3val=total(mprdata(index))/n_elements(index)
          xx=total(xindex(index))/n_elements(index)
       endif
       mprmon(iyear-min(lyear2))=o3val
       oplot,[xx,xx],[o3val,o3val],psym=8
    endif
endfor
;
; zonal mean timeseries
;
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
!type=2^2+2^3
set_viewport,xmn,xmx,ymn,ymx
imin=0.
imax=1.
;print,sth,'  ',slat,'  ',imin,imax
if abs(imin) eq 0. and abs(imax) eq 0. then goto,jump2
nlev=11
col1=1.0+(findgen(nlev)/nlev)*mcolor
col1=lc+0.*findgen(nlev)
cint=(imax-imin)/(nlev-1.)
level=imin+cint*findgen(nlev)
index0=where(strmid(sfileall,0,6) eq 'jan_01',nytick)
plot,findgen(nday*nyear),mark1d,xrange=[0,nfile*nyear-1L],psym=0,yrange=[0.,100.],$
     xticks=nytick-1,xtickv=index0,xtickname=strmid(sfileall(index0),7,2),/noeras,title='ZM Anticyclones'

markmon=fltarr(nx)
for iyear=min(lyear2),max(lyear2) do begin
    index=where(lyear eq iyear)
    if index(0) ne -1L then begin
       xx=xindex(index(0))
;      o3val=markdata(index(0))
;
; require at least 2 obs
;
       o3val=-99.
       if n_elements(index) gt 1L then begin
          o3val=total(markdata(index))/n_elements(index)
          xx=total(xindex(index))/n_elements(index)
          oplot,xindex(index),markdata(index),psym=8,symsize=0.5,color=mcolor*.8
       endif
       markmon(iyear-min(lyear2))=o3val
       oplot,[xx,xx],[o3val,o3val],psym=8,color=mcolor*.9
    endif
endfor
;
; correlation coefficient for this latitude/altitude bin
;
index=where(finite(mprmon) eq 1 and markmon ne -99.)
if index(0) ne -1L then begin
result=correlate(mprmon(index),markmon(index))
xyouts,50.,90.,'r= '+strcompress(result(0)),/data,charsize=1.25
yzcorr(ilat,ialt)=result(0)
print,rlat,ralt,result(0)
endif

;if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_Mprmark_'+smm0+'_'+sth+'K_'+slat+'_uars+mark_corr.ps '+$
         '-rotate -90 yz_Mprmark_'+smm0+'_'+sth+'K_'+slat+'_uars+mark_corr.jpg'
   spawn,'/usr/bin/rm yz_Mprmark_'+smm0+'_'+sth+'K_'+slat+'_uars+mark_corr.ps'
endif
jump2:
endfor	; loop over theta
endfor	; loop over latitude
endfor	; loop over month
;
; latitude/altitude plot of correlation coefficient
;
;erase
;!type=2^2+2^3
;set_viewport,0.15,0.85,0.15,0.85
;nlvls=21
;level=-1.0+0.1*findgen(nlvls)
;col1=1.0+(findgen(nlvls)/nlvls)*mcolor
;contour,yzcorr,ymids,th2,/noeras,levels=level,/cell_fill,c_color=col1,$
;        xrange=[min(ymids),max(ymids)],yrange=[min(th2),max(th2)],$
;        title=smonth+' MetO Correlation',/normal,charsize=2.5,$
;        xtitle='Latitude',ytitle='Theta'
;index=where(level lt 0.)
;contour,yzcorr,ymids,th2,/noeras,levels=level(index),/follow,color=mcolor,c_labels=0*level(index)
;index=where(level gt 0.)
;contour,yzcorr,ymids,th2,/noeras,levels=level(index),/follow,color=0,c_labels=0*level(index)

end
