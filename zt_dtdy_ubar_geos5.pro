;
; plot average annual cycle as a function of altitude of the number of times
; major and minor SSW conditions are satisfied
;
; VLH 2/3/12
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=700
nydim=700
xorig=[0.1,0.5,0.1,0.5]
yorig=[0.5,0.5,0.2,0.2]
xlen=0.3
ylen=0.2
cbaryoff=0.08
cbarydel=0.02
setplot='ps'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
mon=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
;
; restore climo of SSW proxy arrays
;
;COMMENT         STRING    = 'dT/dy=Tbar90-Tbar60 and Ubar60 is zonal mean zonal wind at 60N/S'
;NH_DTDY         FLOAT     = Array[6720, 22]
;NH_UBAR         FLOAT     = Array[6720, 22]
;SH_DTDY         FLOAT     = Array[6720, 22]
;SH_UBAR         FLOAT     = Array[6720, 22]
;TH              FLOAT     = Array[22]
;YYYYMMDD        LONG      = Array[6720]
;
restore,file='GEOS5_dTdy_Ubar_SSW_Climo.sav
nl=n_elements(th)
yyyymmdd_all=YYYYMMDD
syyyymmdd_all=strcompress(yyyymmdd_all,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
sday=strmid(syyyymmdd_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
xindex=where(smon eq '07' and sday eq '15' and long(syear) mod 2 eq 0,nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
xlabs=syear(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)

nmonth=12L
majorsswnh=0./0.*fltarr(nmonth,nl)
minorsswnh=0./0.*fltarr(nmonth,nl)
majorsswsh=0./0.*fltarr(nmonth,nl)
minorsswsh=0./0.*fltarr(nmonth,nl)
imonth=strmid(strcompress(yyyymmdd_all,/remove_all),4,2)
for i=0L,nmonth-1L do begin		; loop over month
    thismonth=string(format='(i2.2)',i+1)
    index=where(imonth eq thismonth,idays)
;
; strip out this month
; 
    nh_minor=NH_DTDY(index, *)
    nh_major=NH_UBAR(index,*)
    sh_minor=SH_DTDY(index, *)
    sh_major=SH_UBAR(index,*)
;
; count number of warmings at each altitude
;
    for k=0L,nl-1L do begin
        if i+1 ge 10 or i+1 le 3 then begin
           index=where(nh_major(*,k) lt 0.,ndays)
           if index(0) ne -1L then majorsswnh(i,k)=100.*ndays/float(idays)
           index=where(nh_minor(*,k) gt 0.,ndays)
           if index(0) ne -1L then minorsswnh(i,k)=100.*ndays/float(idays)
        endif
        if i+1 ge 4 and i+1 le 10 then begin
           index=where(sh_major(*,k) lt 0.,ndays)
           if index(0) ne -1L then majorsswsh(i,k)=100.*ndays/float(idays)
           index=where(sh_minor(*,k) gt 0.,ndays)
           if index(0) ne -1L then minorsswsh(i,k)=100.*ndays/float(idays)
        endif

    endfor
endfor
majorsswnh_save=majorsswnh
minorsswnh_save=minorsswnh
for k=0L,nl-1L do begin
    majorsswnh(0:5,k)=majorsswnh_save(6:11,k)
    minorsswnh(0:5,k)=minorsswnh_save(6:11,k)
    majorsswnh(6:11,k)=majorsswnh_save(0:5,k)
    minorsswnh(6:11,k)=minorsswnh_save(0:5,k)
endfor

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_dtdy_ubar_geos5_'+yearlab+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
erase
level=5.*findgen(21)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*mcolor/nlvls
!type=2^2+2^3+2^7
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nhlabs=['7','8','9','10','11','12','1','2','3','4','5','6']
xyouts,.2,.775,'GEOS5 '+yearlab,/normal,charsize=2.5,charthick=2,color=0
contour,minorsswnh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,ytitle='Theta (K)',title='NH Minor SSWs %',$
        xrange=[1.,nmonth],xticks=nmonth-1,xtickname=nhlabs,xtickv=1+findgen(nmonth),$
        charsize=1.25,charthick=2,xticklen=-0.05
contour,minorsswnh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        /follow,color=mcolor,/overplot,c_labels=0*index
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,minorsswsh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,title='SH Minor SSWs %',$
        xrange=[1.,nmonth],xticks=nmonth-1,xtickname=strcompress(long(1+findgen(nmonth)),/remove_all),xtickv=1+findgen(nmonth),$
        charsize=1.25,charthick=2,xticklen=-0.05
contour,minorsswsh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        /follow,color=mcolor,/overplot,c_labels=0*index
imin=min(level)
imax=max(level)
xmnb=xorig(1)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(1)+cbarydel,yorig(1)+ylen-cbarydel
!type=2^2+2^3+2^5+2^7
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,/noeras,charsize=1.25,charthick=2
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

!type=2^2+2^3+2^7
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=2.*findgen(21)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*mcolor/nlvls
bad=where(majorsswnh eq 0.)
if bad(0) ne -1L then majorsswnh(bad)=0./0.
contour,majorsswnh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,ytitle='Theta (K)',title='NH Major SSWs %',$
        xrange=[1.,nmonth],xticks=nmonth-1,xtickname=nhlabs,xtickv=1+findgen(nmonth),$
        charsize=1.25,charthick=2,xticklen=-0.05
contour,majorsswnh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=2.*level,$
        /follow,color=mcolor,/overplot,c_labels=0*index
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,majorsswsh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
        c_color=col1,/cell_fill,color=0,title='SH Major SSWs %',$
        xrange=[1.,nmonth],xticks=nmonth-1,xtickname=strcompress(long(1+findgen(nmonth)),/remove_all),xtickv=1+findgen(nmonth),$
        charsize=1.25,charthick=2,xticklen=-0.05
;contour,majorsswsh,1.+findgen(nmonth),th,yrange=[400.,2000.],/noeras,levels=level,$
;        /follow,color=mcolor,/overplot,c_labels=0*index
imin=min(level)
imax=max(level)
xmnb=xorig(3)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(3)+cbarydel,yorig(3)+ylen-cbarydel
!type=2^2+2^3+2^5+2^7
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],$
     color=0,/noeras,charsize=1.25,charthick=2
xbox=[0,10,10,0,0]
y1=imin
dy=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
;
; save jpg
;
if setplot eq 'ps' then begin
device,/close
spawn,'convert -trim zt_dtdy_ubar_geos5_'+yearlab+'.ps -rotate -90 zt_dtdy_ubar_geos5_'+yearlab+'.jpg'
endif
end
