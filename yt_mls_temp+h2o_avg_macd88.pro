;
; annual cycle average from 2005-2012
; altitude-time series of MLS temperature + H2O
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.8
ylen=0.275
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
dir='/Volumes/earth/aura6/data/MLS_data/Datfiles/'
lstmn=1
lstdy=1
lstyr=2005
ledmn=12
leddy=31
ledyr=2012
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
minyear=lstyr
maxyear=ledyr
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nr=91L
latbin=-90.+2.*findgen(nr)
;goto,quick

z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=365L
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

for iyear=lstyr,ledyr do begin
    syear=strcompress(iyear,/remove_all)
    if kcount eq 0L then begin
       mlsh2o_yt=fltarr(kday,nr)
       mlstp_yt=fltarr(kday,nr)
       nmlsh2o_yt=fltarr(kday,nr)
       nmlstp_yt=fltarr(kday,nr)
       kcount=1
    endif
    restore,'yt_mls_temp_h2o_pmc_'+syear+'.sav'	;,h2o,temp,latbin,sdates
;
; remove leap day
;
    index=where(sdates ne '20080229')
    temp=temp(index,*)
    h2o=h2o(index,*)
    if iyear eq lstyr+1 then sdates_orig=sdates
    if iyear eq lstyr then sdate_all=sdates(index)
    if iyear gt lstyr then sdate_all=[sdate_all,sdates(index)]
;
; sum
;
   good=where(temp ne 0.)
   if good(0) ne -1L then mlstp_yt(good)=mlstp_yt(good)+temp(good)
   if good(0) ne -1L then nmlstp_yt(good)=nmlstp_yt(good)+1.
   good=where(h2o ne 0.)
   if good(0) ne -1L then mlsh2o_yt(good)=mlsh2o_yt(good)+h2o(good)
   if good(0) ne -1L then nmlsh2o_yt(good)=nmlsh2o_yt(good)+1.

endfor
;
; average
;
index=where(nmlstp_yt ne 0.)
if index(0) ne -1L then mlstp_yt(index)=mlstp_yt(index)/nmlstp_yt(index)
index=where(nmlsh2o_yt ne 0.)
if index(0) ne -1L then mlsh2o_yt(index)=mlsh2o_yt(index)/nmlsh2o_yt(index)
;
; save temp, h2o, etc
;
h2o=mlsh2o_yt
temp=mlstp_yt
save,file='yt_mls_temp_h2o_pmc_'+yearlab+'.sav',h2o,temp,latbin,sdates
quick:
restore,'yt_mls_temp_h2o_pmc_'+yearlab+'.sav'
mlsh2o_yt=h2o
mlstp_yt=temp
sdate_all=sdates_orig
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)	;+'/'+sday(xindex)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yt_mls_temp+h2o_'+yearlab+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot latitude-time temperature and water
;
erase
xyouts,0.3,0.95,yearlab+' MLS 0.005 hPa',/normal,color=0,charsize=2
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(mlsh2o_yt eq 0.)
if index(0) ne -1L then mlsh2o_yt(index)=0./0.
mlsh2o_yt=smooth(mlsh2o_yt,3,/NaN)
index=where(mlstp_yt eq 0.)
if index(0) ne -1L then mlstp_yt(index)=0./0.
mlstp_yt=smooth(mlstp_yt,3,/NaN)
tlevel=130.+4.*findgen(23)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlstp_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='Temperature',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlstp_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
contour,mlstp_yt,1.+findgen(kday),latbin,levels=[160.],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=0.25*findgen(24)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlsh2o_yt,1.+findgen(kday),latbin,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/cell_fill,c_color=col1,title='Water Vapor',$
      levels=tlevel,yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
contour,mlsh2o_yt,1.+findgen(kday),latbin,levels=[4.],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(ppmv)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor


    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim yt_mls_temp+h2o_'+yearlab+'.ps -rotate -90 yt_mls_temp+h2o_'+yearlab+'.jpg'
    endif
end
