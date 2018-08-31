;
; save areas in equivalent latitude
; save GEOS-5 vortex area as a function of altitude and time
; save GEOS-5 vortex area based on Nash definition
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto
@calcelat2d

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
latcircle=fltarr(nrr)
hem_frac=fltarr(nrr)
for j=0,nrr-2 do begin
    hy=re*dtr
    dx=re*cos(yeq(j)*dtr)*360.*dtr
    latcircle(j)=dx*hy
endfor
for j=0,nrr-1 do begin
    if yeq(j) ge 0. then index=where(yeq ge yeq(j))
    if index(0) ne -1 then hem_frac(j)=100.*total(latcircle(index))/hem_area
    if yeq(j) eq 0. then hem_frac(j)=100.
endfor

lstmn=11L & lstdy=1L & lstyr=6L 
ledmn=3L & leddy=1L & ledyr=7L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15,0.15,0.15]
yorig=[0.75,0.50,0.15]
xlen=0.8
ylen=0.2
cbaryoff=0.04
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
goto,plotit
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto, saveit

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump

      ncid=ncdf_open(dir+sdate+stimes(0)+'nc5')
      marknew2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,marknew2
      ncdf_close,ncid

      if kcount eq 0L then begin
         area_zt_nc4=fltarr(kday,nth)
         area_zt_nc5=fltarr(kday,nth)
         area_zt_nash=fltarr(kday,nth)
         harea_zt=fltarr(kday,nth)
         sfile=strarr(kday)
         dum=transpose(mark2(*,*,0))
         lon=0.*dum
         lat=0.*dum
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy    ; area of each grid point
         endfor
         kcount=1L
      endif
      sfile(icount)=lfile
;
; loop over theta
;
      for thlev=0,nth-1 do begin
          mark1=transpose(mark2(*,*,thlev))
          index=where(lat gt 0. and mark1 gt 0.0,nn)
          if index(0) ne -1 then begin
             phem=100.*total(area(index))/hem_area	; % of hemisphere
             for jj=1L,nrr-1L do begin
                 if hem_frac(jj) le phem and hem_frac(jj-1) gt phem then begin
                    scale=(hem_frac(jj)-phem)/(hem_frac(jj)-hem_frac(jj-1))
                    area_zt_nc4(icount,thlev)=yeq(jj)-scale*(yeq(jj)-yeq(jj-1))	; Equivalent latitude
;print,hem_frac(jj-1),phem,hem_frac(jj),scale
;print,yeq(jj-1),area_zt_nc4(icount,thlev),yeq(jj)
;stop
                 endif
             endfor
          endif
              
          index=where(lat gt 0. and mark1 lt 0.0,nn)
          if index(0) ne -1 then begin
             phem=100.*total(area(index))/hem_area      ; % of hemisphere
             for jj=1L,nrr-1L do begin
                 if hem_frac(jj) le phem and hem_frac(jj-1) gt phem then begin
                    scale=(hem_frac(jj)-phem)/(hem_frac(jj)-hem_frac(jj-1))
                    harea_zt(icount,thlev)=yeq(jj)-scale*(yeq(jj)-yeq(jj-1))        ; Equivalent latitude
                 endif
             endfor
          endif

          mark1=transpose(marknew2(*,*,thlev))
          index=where(lat gt 0. and mark1 gt 0.0,nn)
          if index(0) ne -1 then begin
             phem=100.*total(area(index))/hem_area      ; % of hemisphere
             for jj=1L,nrr-1L do begin
                 if hem_frac(jj) le phem and hem_frac(jj-1) gt phem then begin
                    scale=(hem_frac(jj)-phem)/(hem_frac(jj)-hem_frac(jj-1))
                    area_zt_nc5(icount,thlev)=yeq(jj)-scale*(yeq(jj)-yeq(jj-1)) ; Equivalent latitude
                 endif
             endfor
          endif
;
; area of Nash vortex
;
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          speed1=sqrt(u1^2+v1^2)
          pv1=transpose(pv2(*,*,thlev))
          elat1=calcelat2d(pv1,alon,alat)
;
; integrate wind speed and PV in Elat bins
;
          nbins=37
          dy=2.5
          latmin=0.
          latmax=90.
          elatbin=latmin+dy*findgen(nbins)
          speedbin=-999.+0.*fltarr(nbins)                               ; average windspeed per elat bin
          pvbin=0.*elatbin
          for n=0,nbins-2 do begin
              t=where(pv1 lt 1000. and lat ge latmin and elat1 ge elatbin(n) and elat1 lt elatbin(n+1),it)
              if it gt 2 then begin
                  result=moment(pv1(t))
                  pvbin(n)=result(0)
                  if min(lat(t))-latmin le dy then begin ; make sure bins are resolved (do not intersect latmin)
                     speedbin(n)=-999.
                     goto,jumpnhbin
                  endif
                  speedbin(n)=total(speed1(t))/float(it)
              endif
              jumpnhbin:
          endfor                                                        ; loop over Elat bins
          s=where(lat ge latmin and elat1 ge elatbin(nbins-1),is)
          if is gt 2 then begin
             result=moment(pv1(s))
             pvbin(n)=result(0)
             if min(lat(s))-latmin gt dy then speedbin(nbins-1)=total(speed1(s))/float(is)
          endif
;
; compute PV gradient wrt Equivalent latitude
;
          dpvbin=0.*pvbin
          for i=0,nbins-2L do dpvbin(i)=pvbin(i+1)-pvbin(i)
          dpvbin(nbins-1)=pvbin(nbins-1)-pvbin(nbins-2)
;
; impose Nash filter poleward of 80deg (and add new one Equatorward of lat0)
;
          lat0=70.
          index=where(elatbin ge lat0)                                  ; filter down poleward of 80deg
          speedbin(index)=speedbin(index)*(90.-elatbin(index))/30.
          dpvbin(index)=dpvbin(index)*(90.-elatbin(index))/30.
          lat0=25.
          if th(thlev) lt 600. then lat0=45.
          index=where(elatbin le lat0)                                  ; filter down equatorward of lat0
          speedbin(index)=speedbin(index)*(elatbin(index))/(2.*lat0)
          dpvbin(index)=dpvbin(index)*(elatbin(index))/(2.*lat0)
          dpvbin=dpvbin/max(dpvbin)                                     ; normalise
;
; vortex edge is where dPV/dElat multiplied by the wind speed integrated in Elat bins is maximum
; and integrated wind speed must be greater than 15.2 m/s
;
          prod=dpvbin*speedbin
          index=where(prod eq max(prod))
          if index(0) ne -1L then edgepv=pvbin(index)
;print,th(thlev),edgepv(0),speedbin(index(0))
          if speedbin(index(0)) lt 15.2 then goto,skiplev

          index=where(lat gt 0. and pv1 ge edgepv(0),nn)
          if index(0) ne -1 then begin
             phem=100.*total(area(index))/hem_area      ; % of hemisphere
             for jj=1L,nrr-1L do begin
                 if hem_frac(jj) le phem and hem_frac(jj-1) gt phem then begin
                    scale=(hem_frac(jj)-phem)/(hem_frac(jj)-hem_frac(jj-1))
                    area_zt_nash(icount,thlev)=yeq(jj)-scale*(yeq(jj)-yeq(jj-1)) ; Equivalent latitude
;print,hem_frac(jj-1),phem,hem_frac(jj),scale
;print,yeq(jj-1),area_zt_nash(icount,thlev),yeq(jj)
;stop
                 endif
             endfor
          endif
;print,th(thlev),area_zt_nc4(icount,thlev),area_zt_nash(icount,thlev)
skiplev:
      endfor
;stop

skipit:
icount=icount+1L
goto,jump

saveit:
;
; plot altitude-time series of Arctic vortex area
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
y1='20'+string(format='(i2.2)',long(min(yy(index))))
y2='20'+string(format='(i2.2)',long(max(yy(index))))
;
; save file
;
save,file='vortex_area_'+y1+'_'+y2+'_geos5_vs_nash.sav',area_zt_nc4,area_zt_nc5,area_zt_nash,$
     harea_zt,th,sfile,y1,y2

plotit:
restore,'vortex_area_2006_2007_geos5_vs_nash.sav'
kday=n_elements(SFILE)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
   device,/landscape,/bold,/color,bits_per_pixel=8,/helvetica,filename='zt_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,kday,kday,1,1],[500.,500.,4600.,4600.,500.],min_value=0.,$
      xrange=[1,kday],yrange=[500.,4600.],/nodata,charsize=1.5,color=0,$
      ytitle='Theta (K)',title=y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,3,2) eq '15',nxtick)
if kindex(0) ne -1L then begin
xmon=long(strmid(sfile(kindex),0,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,300.
    plots,kindex(i)+1,500.,/continue,/data,color=0
    xyouts,kindex(i)+1,0.,xlab,/data,alignment=0.5,charsize=1.5,color=0
endfor
endif
xyouts,xmn+xlen+0.02,0.06+(ymx+ymn)/2.,'Kinematic',orientation=-90.,/normal,color=0,charsize=1.5
nlvls=21
col1=reverse(1+indgen(nlvls)*icolmax/nlvls)
level=30.+2.5*findgen(nlvls)
level2=[30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.]
area_zt=smooth(area_zt_nc4,5,/edge_truncate)
harea_zt=smooth(harea_zt,5,/edge_truncate)
contour,area_zt,1.+findgen(kday),th,levels=level,/fill,/cell_fill,/overplot,c_color=col1
contour,area_zt,1.+findgen(kday),th,levels=level2,c_color=0,/follow,/overplot
;contour,harea_zt,1.+findgen(kday),th,levels=level2,c_color=mcolor,/follow,/overplot,thick=4
area_zt_nc4=area_zt

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,kday,kday,1,1],[500.,500.,4600.,4600.,500.],min_value=0.,$
      xrange=[1,kday],yrange=[500.,4600.],/nodata,charsize=1.5,color=0,$
      ytitle='Theta (K)',xtickname=[' ',' '],xticks=1
if kindex(0) ne -1L then begin
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,300.
    plots,kindex(i)+1,500.,/continue,/data,color=0
    xyouts,kindex(i)+1,0.,xlab,/data,alignment=0.5,charsize=1.5,color=0
endfor
endif
xyouts,xmn+xlen+0.02,0.02+(ymx+ymn)/2.,'Nash',orientation=-90.,/normal,color=0,charsize=1.5
area_zt=smooth(area_zt_nash,5,/edge_truncate)
harea_zt=smooth(harea_zt,5,/edge_truncate)
contour,area_zt,1.+findgen(kday),th,levels=level,/fill,/cell_fill,/overplot,c_color=col1
contour,area_zt,1.+findgen(kday),th,levels=level2,c_color=0,/follow,/overplot
;contour,harea_zt,1.+findgen(kday),th,levels=level2,c_color=mcolor,/follow,/overplot,thick=4
area_zt_nash=area_zt
imin=min(level)
imax=max(level)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Vortex Equivalent Latitude',charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,kday,kday,1,1],[500.,500.,4600.,4600.,500.],min_value=0.,$
      xrange=[1,kday],yrange=[500.,4600.],/nodata,charsize=1.5,color=0,$
      ytitle='Theta (K)',xtickname=[' ',' '],xticks=1
if kindex(0) ne -1L then begin
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,500.
    plots,kindex(i)+1,500.,/continue,/data,color=0
    xyouts,kindex(i)+1,0.,xlab,/data,alignment=0.5,charsize=1.5,color=0
endfor
endif
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
nlvls=n_elements(col2)
area_zt=0.*area_zt_nc4
index=where(abs(area_zt_nc4) gt 1. and abs(area_zt_nash) gt 1.)
area_zt(index)=area_zt_nc4(index)-area_zt_nash(index)
level=-25.+5.*findgen(nlvls)
contour,area_zt,1.+findgen(kday),th,levels=level,/fill,/cell_fill,/overplot,c_color=col2
index=where(level gt 0. or level lt 0.)
contour,area_zt,1.+findgen(kday),th,levels=level(index),c_labels=1.+0.*level(index),c_color=0,/follow,/overplot

imin=min(level)
imax=max(level)
ymnb=yorig(2) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Equivalent Latitude Difference',charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.ps -rotate -90 '+$
         'zt_geos5_vortex_area_'+y1+'-'+y2+'_vs_nash.jpg'
endif
end
