;
; plot entire 14 years in 1 longitude-time section
; 4 panel, includes Tpr hov, mark hov, and timeseries of both
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
xorig=[0.15,0.55,0.15,0.55]
yorig=[0.6,0.6,0.15,0.15]
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
    restore,dir+ifiles(n)
    print,'restored ',ifiles(n)
    nth=n_elements(th2)
    nr=n_elements(YMIDS)
    nc=n_elements(alon)
    if n eq 0L then begin
       uxtall=fltarr(nc,nday,nyear)
       tmpxtall=fltarr(nc,nday,nyear)
       msfxtall=fltarr(nc,nday,nyear)
       markxtall=fltarr(nc,nday,nyear)
       uxtmean=fltarr(nc,nday)
       tmpxtmean=fltarr(nc,nday)
       msfxtmean=fltarr(nc,nday)
       markxtmean=fltarr(nc,nday)
       sfileall=strarr(nday*nyear)
       rlat=0.
       print,ymids
       read,'Enter desired latitude ',rlat
       index=where(ymids eq rlat)
       ilat=index(0)
       ralt=400.
       print,th2
       read,'Enter desired altitude ',ralt
       index=where(th2 eq ralt)
       ialt=index(0)
    endif
    index=where(strmid(sfile,0,6) ne 'feb_29',nfile)
    uxtall(*,*,n)=uxt(*,index,ilat,ialt)
    tmpxtall(*,*,n)=tmpxt(*,index,ilat,ialt)
    msfxtall(*,*,n)=msfxt(*,index,ilat,ialt)
    markxtall(*,*,n)=markxt(*,index,ilat,ialt)
    uxtmean=uxtmean+uxt(*,index,ilat,ialt)
    tmpxtmean=tmpxtmean+tmpxt(*,index,ilat,ialt)
    msfxtmean=msfxtmean+msfxt(*,index,ilat,ialt)
    markxtmean=markxtmean+markxt(*,index,ilat,ialt)
    sfileall(n*nday:((n+1L)*nday)-1L)=sfile(index)
endfor
;
; deviation from time mean in Temperature and MSF
;
tmpxtmean=tmpxtmean/float(nyear)
msfxtmean=msfxtmean/float(nyear)
tprxt=-9999.+0.*tmpxtall
msfprxt=-9999.+0.*msfxtall
for n=0,nyear-1L do begin
for j=0,nday-1L do begin
for i=0,nc-1L do begin
    if tmpxtall(i,j,n) ne 0. then begin
       tprxt(i,j,n)=tmpxtall(i,j,n)-tmpxtmean(i,j)
       msfprxt(i,j,n)=msfxtall(i,j,n)-msfxtmean(i,j)
    endif
endfor
endfor
endfor
;
; reform tprxt to be a 2-D array over entire time period
;
tprxt2d=fltarr(nc,nday*nyear)
msfprxt2d=fltarr(nc,nday*nyear)
markxt2d=fltarr(nc,nday*nyear)
for n=0L,nyear-1L do begin
    tprxt2d(*,n*nday:((n+1L)*nday)-1L)=tprxt(*,*,n)
    msfprxt2d(*,n*nday:((n+1L)*nday)-1L)=msfprxt(*,*,n)
    markxt2d(*,n*nday:((n+1L)*nday)-1L)=markxtall(*,*,n)
endfor
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
   device,/landscape,bits=8,filename='hov_Tprmark_'+sth+'K_'+slat+'_uars+mark_timeseries.ps'
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
xyouts,.4,.95,'MetO '+sth+' K '+slat,/normal,charsize=2.5
plt=reform(tprxt2d,nc,nfile*nyear)
;index=where(plt eq 0.)
;if index(0) ne -1L then plt(index)=-9999./0.
;plt=smooth(plt,30,/NaN,/edge_truncate)
imin=-10.
imax=10.
print,sth,'  ',slat,'  ',imin,imax
if abs(imin) eq 0. and abs(imax) eq 0. then goto,jump2
nlev=11
col1=1.0+(findgen(nlev)/nlev)*mcolor
cint=(imax-imin)/(nlev-1.)
level=imin+cint*findgen(nlev)
index=where(strmid(sfileall,0,6) eq 'jan_01',nytick)
contour,plt,alon,findgen(nfile*nyear),xrange=[0.,360.],/fill,$
        /cell_fill,yrange=[0,nfile*nyear-1L],xstyle=1,ystyle=1,xticks=6,$
        yticks=nytick-1,ytickv=index,ytickname=strmid(sfileall(index),7,2),$
        xtitle='Longitude',ytitle='Time',c_color=col1,$
        title='Tpr',/noeras,levels=level,charsize=2,min_value=-9999.
ymnb=yorig(0)-cbaryoff
ymxb=ymnb +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlev)
for jj=0,nlev-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(jj)
    x1=x1+dx
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
plt1d=fltarr(nday*nyear)
for i=0L,nday*nyear-1L do plt1d(i)=total(plt(*,i))/float(nc)
plot,findgen(nday*nyear),plt1d,xrange=[0,nfile*nyear-1L],psym=3,yrange=[-10.,10.],$
     xticks=nytick-1,xtickv=index,xtickname=strmid(sfileall(index),7,2),/noeras,title='ZM Tpr'

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
!type=2^2+2^3
set_viewport,xmn,xmx,ymn,ymx
plt=abs(reform(markxt2d,nc,nfile*nyear))
index=where(plt eq 0.)
if index(0) ne -1L then plt(index)=-9999.
;plt=smooth(plt,10,/NaN,/edge_truncate)
imin=0.
imax=1.
print,sth,'  ',slat,'  ',imin,imax
if abs(imin) eq 0. and abs(imax) eq 0. then goto,jump2
nlev=11
col1=1.0+(findgen(nlev)/nlev)*mcolor
col1=lc+0.*findgen(nlev)
cint=(imax-imin)/(nlev-1.)
level=imin+cint*findgen(nlev)
index0=where(strmid(sfileall,0,6) eq 'jan_01',nytick)
contour,plt,alon,findgen(nfile*nyear),xrange=[0.,360.],/fill,$
        /cell_fill,yrange=[0,nfile*nyear-1L],xstyle=1,ystyle=1,xticks=6,$
        yticks=nytick-1,ytickv=index0,ytickname=strmid(sfileall(index0),7,2),$
        xtitle='Longitude',ytitle='Time',c_color=col1,$
        title='Anticyclones',/noeras,levels=level,charsize=2,min_value=-9999.
;ymnb=yorig(2)-cbaryoff
;ymxb=ymnb +cbarydel
;set_viewport,xmn,xmx,ymnb,ymxb
;!type=2^2+2^3+2^6
;plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,charsize=1.5
;ybox=[0,10,10,0,0]
;x1=imin
;dx=(imax-imin)/float(nlev)
;for jj=0,nlev-1 do begin
;    xbox=[x1,x1,x1+dx,x1+dx,x1]
;    polyfill,xbox,ybox,color=col1(jj)
;    x1=x1+dx
;endfor
;
; zonal mean timeseries
;
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
!type=2^2+2^3
set_viewport,xmn,xmx,ymn,ymx
plt1d=fltarr(nday*nyear)
for i=0L,nday*nyear-1L do begin
    index=where(plt(*,i) ne -9999. and plt(*,i) ne 0.,nc2)
    if index(0) ne -1L then plt1d(i)=total(plt(index,i))	;/float(nc2)
endfor
plot,findgen(nday*nyear),plt1d,xrange=[0,nfile*nyear-1L],psym=0,yrange=[0.,50.],$
     xticks=nytick-1,xtickv=index0,xtickname=strmid(sfileall(index0),7,2),/noeras,title='ZM Anticyclones'

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim hov_Tprmark_'+sth+'K_'+slat+'_uars+mark_timeseries.ps '+$
         '-rotate -90 hov_Tprmark_'+sth+'K_'+slat+'_uars+mark_timeseries.jpg'
endif
jump2:
end
