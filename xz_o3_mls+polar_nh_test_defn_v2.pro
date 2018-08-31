;
; test that save file exists with MLS and marker interpolation done
; ***MLS only***
; test LOP definition.  v2 defn: look at anticyclone mean/sigma and ambient mean/sigma
; within X km from anticyclone center of mass.
; VLH 9/7/2006
;
@aura2date
@stddat
@kgmt
@ckday
@kdate
@range_ring

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.10]
yorig=[0.10]
cbaryoff=0.015
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=0
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
ifile='                             '
lstmn=1 & lstdy=10 & lstyr=2006 & lstday=0
ledmn=1 & leddy=20 & ledyr=2006 & ledday=0
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
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
      if ndays gt ledday then stop,' Normal termination condition '
;
; read UKMO data
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      idate=long(smn+sdy)
      uyr=strmid(syr,2,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
dum=findfile('/Users/harvey/Analysis/Datfiles_MLS/mls+mark_'+syr+smn+sdy+'.sav')
if dum(0) eq '' then goto,jump
      print,iyr,imn,idy
;
; postscript file
;
      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,$
                 filename='MLS_LOPs/xz_o3_mls+polar_nh_'+lfile+'_v2.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; restore save file with MLS and MetO marker
;
restore,'/Users/harvey/Analysis/Datfiles_MLS/mls+mark_'+syr+smn+sdy+'.sav'
nc=n_elements(alon)
nr=n_elements(alat)
nth=n_elements(th)
mlon=reform(mlon2(0,*))
mlat=reform(mlat2(0,*))
muttime=reform(mtime2(0,*))
;
; polar projection of MetO data at 1000 K
;
      rtheta=1000.
      index=where(rtheta eq th)
      itheta=index(0)
      stheta=strcompress(string(fix(th(itheta))),/remove_all)
      sf1=transpose(sf2(*,*,itheta))
      mark1=transpose(mark2(*,*,itheta))
      sf=fltarr(nc+1,nr)
      sf(0:nc-1,0:nr-1)=sf1
      sf(nc,*)=sf(0,*)
      mark=fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
;
; polar plot of ozone
;
      erase
      !type=2^2+2^3
      xyouts,.35,.975,syr+smn+sdy+'  '+stheta+' K',/normal,color=0,charsize=2
      set_viewport,.15,.45,.65,.95
      MAP_SET,90,0,-180,/ortho,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5,title='Ozone'
      oplot,findgen(361),0.1+0.*findgen(361),psym=0,color=0
      contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark gt 0.05 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark lt -0.05 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      a=findgen(8)*(2*!pi/8.)
      usersym,cos(a),sin(a),/fill
      index=where(mtheta2 gt rtheta-80. and mtheta2 le rtheta+80. and mlat2 gt 0.,npt)
      if index(0) eq -1L then goto,jump
      xdata=mlon2(index) & ydata=mlat2(index)
      o3data=mo3(index)
      o3mask=mo3mask(index)
      o3mark=mmark2(index)
      index=where(xdata lt 0.)
      if index(0) ne -1 then xdata(index)=xdata(index)+360.
      omin=2.0
      omax=10.
      for i=0L,n_elements(o3data)-1L do begin
          if o3mask(i) ne -99. then $
             oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                   color=((o3data(i)-omin)/(omax-omin))*mcolor,symsize=1.5
      endfor
      loadct,0
      contour,mark,x,alat,levels=[.05],c_color=mcolor*.7,/overplot,/follow,c_labels=0,/noeras,thick=5
      loadct,38
      contour,mark,x,alat,levels=[-.05],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=5
;
; determine number of NH swaths and print time of swaths on polar plot
;
      index=where(mlat ge 81.7,npt)
      tswath=muttime(index)
      flag=0.*tswath
      nswath=50L
      tsave=fltarr(nswath)
      mcount=0L
      for i=0L,npt-1L do begin
          index=where(abs(tswath(i)-tswath) lt 1. and flag eq 0.)
          if index(0) ne -1L then begin
             flag(index)=1.0
             kindex=where(abs(muttime-tswath(index(0))) le 0.5 and mlat gt 0.)
;            oplot,mlon(kindex),mlat(kindex),psym=8,symsize=1.25,color=(float(i+1)/float(npt))*mcolor
             stime=string(FORMAT='(I2.2)',long(tswath(index(0))))
             tsave(mcount)=tswath(index(0))
             mcount=mcount+1L
             xtmp=mlon(kindex)
             ytmp=mlat(kindex)
             index=where(ytmp eq min(ytmp))
             xyouts,xtmp(index(0)),30.,stime,/data,charsize=1.5,alignment=0.5,color=0,charthick=3
          endif
      endfor
      tsave=tsave(0:mcount-1L)
      MAP_SET,90,0,-180,/ortho,/contin,/grid,/noeras,color=0,/noborder,charsize=1.5
;
; polar plot of marker
;
      set_viewport,.55,.85,.65,.95
      MAP_SET,90,0,-180,/ortho,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5,title='Marker'
      oplot,findgen(361),0.1+0.*findgen(361),psym=0,color=0
      contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark gt 0.05 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark lt -0.05 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      a=findgen(8)*(2*!pi/8.)
      usersym,cos(a),sin(a),/fill
      index=where(mtheta2 gt rtheta-80. and mtheta2 le rtheta+80. and mlat2 gt 0.,npt)
      if index(0) eq -1L then goto,jump
      xdata=mlon2(index) & ydata=mlat2(index)
      o3data=mo3(index)
      o3mask=mo3mask(index)
      o3mark=mmark2(index)
      o3theta=mtheta2(index)
      index=where(xdata lt 0.)
      if index(0) ne -1 then xdata(index)=xdata(index)+360.
      mmin=min(o3mark)
      nlvls=abs(mmin)
      col1=1+0.2*indgen(nlvls)*mcolor
      for i=0L,n_elements(o3data)-1L do begin
          if o3mask(i) ne -99. then begin
             if o3mark(i) eq 1.0 then oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                   color=0.9*mcolor,symsize=1.5
          endif
          if o3mask(i) ne -99. then begin
             if abs(o3mark(i)) eq 0.0 then oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                   color=0.65*mcolor,symsize=1.5
          endif
          if o3mask(i) ne -99. then begin
             for ihigh=1L,abs(mmin) do begin
             if o3mark(i) eq -1L*ihigh then $
                oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,color=col1(ihigh-1L),symsize=1.5
             endfor
          endif
      endfor
      MAP_SET,90,0,-180,/ortho,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5
      loadct,0
      contour,mark,x,alat,levels=[.05],c_color=mcolor*.7,/overplot,/follow,c_labels=0,/noeras,thick=5
      loadct,38
      contour,mark,x,alat,levels=[-.05],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=5
;
; loop over theta levels ranging from 2000-500K
;
ho3_save=-99.
ao3_save=-99.
hth_save=-99.
ath_save=-99.
hx_save=-99.
ax_save=-99.
hy_save=-99.
ay_save=-99.

ho3mean_save=fltarr(nth)
ao3mean_save=fltarr(nth)
hthmean_save=fltarr(nth)
athmean_save=fltarr(nth)
ho3sigma_save=fltarr(nth)
ao3sigma_save=fltarr(nth)
hthsigma_save=fltarr(nth)
athsigma_save=fltarr(nth)
nlvls=12L
col1=1+indgen(nlvls)*mcolor/(float(nlvls))
for k=0L,11L do begin
    rtheta=th(k)
    index=where(rtheta eq th)
    itheta=index(0)
    if itheta eq 0 then dtheta=(th(0)-th(1))/2.0
    if itheta eq nth-1L then dtheta=(th(nth-2)-th(nth-1))/2.0
    if itheta gt 0 and itheta le nth then dtheta=(th(itheta-1)-th(itheta+1))/2.0
    index=where(mtheta2 le th(0) and mtheta2 ge th(11) and $
                mtheta2 gt rtheta-dtheta and mtheta2 le rtheta+dtheta and mlat2 gt 0.,npt)
    if index(0) eq -1L then goto,jumplev
    xdata=mlon2(index) & ydata=mlat2(index)
;
; normalise by area of each gridpoint.  without this logic, the center of mass is shifted poleward
;
    adata=0.*xdata
    deltax=alon(1)-alon(0)
    deltay=alat(1)-alat(0)
    for jj=0L,npt-1L do begin
        hy=re*deltay*dtr
        dx=re*cos(ydata(jj)*dtr)*deltax*dtr
        adata(jj)=1.	;dx*hy    ; area of each grid point
    endfor
    adata=adata/max(adata)

    o3data=mo3(index)
    o3mask=mo3mask(index)
    o3mark=mmark2(index)
    o3theta=mtheta2(index)
    index=where(xdata lt 0.)
    if index(0) ne -1 then xdata(index)=xdata(index)+360.
;
; look for anticyclone "center of mass" in Pacific sector
;
index=where(xdata ge 90. and xdata le 270.)
mmin=min(o3mark(index))
hindex=where(o3mark eq mmin)
resultx=moment(xdata(hindex))
x0=resultx(0)
resulty=moment(ydata(hindex)*adata(hindex))	; normalise by area
y0=resulty(0)
if y0 le 30. then goto,jumplev		; ignore tropical anticyclones
;
; draw anticyclone edges
;
mark1=transpose(mark2(*,*,itheta))
mark=fltarr(nc+1,nr)
mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
mark(nc,*)=mark(0,*)
contour,mark,x,alat,levels=[-.05],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=5
;
; draw range ring 1000 km beyond anticyclone edge
;
dxf=re*abs(x0-xdata(hindex))*dtr*cos(y0*dtr)
dyf=re*abs(y0-ydata(hindex))*dtr
dist=sqrt(dxf*dxf+dyf*dyf)
dist0=max(dist)+1000.
print,'theta= ',rtheta,'  dtheta= ',dtheta,' center of mass= ',x0,y0,' dist outward= ',dist0
range_ring,y0,x0,dist0,360,bear,lats,lons
oplot,[x0,x0],[y0,y0],psym=8,color=col1(k),symsize=3
oplot,lons,lats,color=col1(k),thick=3
;
; extract ambient ozone within dist0 of anticyclone and plot
;
dxf=re*abs(x0-xdata)*dtr*cos(y0*dtr)
dyf=re*abs(y0-ydata)*dtr
dist=sqrt(dxf*dxf+dyf*dyf)
agood=where(o3mark eq 0. and dist le dist0)
hgood=where(o3mark lt 0. and dist le dist0)
;
; save ambient and anticyclone ozone and theta
;
if hgood(0) ne -1L then begin
;  oplot,xdata(hgood),ydata(hgood),psym=8,color=mcolor,symsize=2
   ho3_save=[ho3_save,o3data(hgood)]
   hth_save=[hth_save,o3theta(hgood)]
   hx_save=[hx_save,xdata(hgood)]
   hy_save=[hy_save,ydata(hgood)]
   if n_elements(hgood) ge 3L then begin
      result=moment(o3data(hgood))
      ho3mean_save(k)=result(0)
      ho3sigma_save(k)=sqrt(result(1))
      result=moment(o3theta(hgood))
      hthmean_save(k)=result(0)
      hthsigma_save(k)=sqrt(result(1))
   endif
endif
if agood(0) ne -1L then begin
;  oplot,xdata(agood),ydata(agood),psym=1,color=0,symsize=2
   ao3_save=[ao3_save,o3data(agood)]
   ath_save=[ath_save,o3theta(agood)]
   ax_save=[ax_save,adata(agood)]
   ay_save=[ay_save,adata(agood)]
   if n_elements(agood) ge 3L then begin
      result=moment(o3data(agood))
      ao3mean_save(k)=result(0)
      ao3sigma_save(k)=sqrt(result(1))
      result=moment(o3theta(agood))
      athmean_save(k)=result(0)
      athsigma_save(k)=sqrt(result(1))
   endif
endif
jumplev:
endfor	; end loop over theta
;
; plot anticyclone and nearby ambient ozone values
;
set_viewport,.25,.75,.1,.6
plot,1.+findgen(9),[th(11),th(0)],/nodata,color=0,charsize=1.5,xtitle='Ozone (ppmv)',ytitle='Theta (K)'
xyouts,6.,1900.,'Anticyclone',color=mcolor*.9,/data,charsize=1.5,charthick=2
xyouts,6.,1800.,'Nearby Ambient',color=0,/data,charsize=1.5,charthick=2

agood=where(ao3_save ne -99.)
if agood(0) ne -1L then begin
   ao3_save=ao3_save(agood)
   ath_save=ath_save(agood)
   ax_save=ax_save(agood)
   ay_save=ay_save(agood)
   oplot,ao3_save(agood),ath_save(agood),psym=8,color=0,symsize=0.5
endif
hgood=where(ho3_save ne -99.)
if hgood(0) ne -1L then begin
   ho3_save=ho3_save(hgood)
   hth_save=hth_save(hgood)
   hx_save=hx_save(hgood)
   hy_save=hy_save(hgood)
   oplot,ho3_save(hgood),hth_save(hgood),psym=8,color=.9*mcolor,symsize=0.5
endif
;
; plot means and sigmas
;
loadct,0
agood=where(ao3mean_save ne 0.)
if agood(0) ne -1L then begin
   ao3mean_save=ao3mean_save(agood)
   athmean_save=athmean_save(agood)
   ao3sigma_save=ao3sigma_save(agood)
   athsigma_save=athsigma_save(agood)
   for ii=0L,n_elements(agood)-1L do begin
       oplot,[ao3mean_save(ii),ao3mean_save(ii)],[athmean_save(ii),athmean_save(ii)],psym=8,color=0,symsize=2
       plots,ao3mean_save(ii)-ao3sigma_save(ii),athmean_save(ii)
       plots,ao3mean_save(ii)+ao3sigma_save(ii),athmean_save(ii),thick=6,/continue,color=mcolor*.7
   endfor
endif
loadct,38

hgood=where(ho3mean_save ne 0.)
if hgood(0) ne -1L then begin
   ho3mean_save=ho3mean_save(hgood)
   hthmean_save=hthmean_save(hgood)
   ho3sigma_save=ho3sigma_save(hgood)
   hthsigma_save=hthsigma_save(hgood)
   for ii=0L,n_elements(hgood)-1L do begin
       oplot,[ho3mean_save(ii),ho3mean_save(ii)],[hthmean_save(ii),hthmean_save(ii)],psym=8,color=mcolor*.8,symsize=2
       plots,ho3mean_save(ii)-ho3sigma_save(ii),hthmean_save(ii)
       plots,ho3mean_save(ii)+ho3sigma_save(ii),hthmean_save(ii),thick=4,/continue,color=mcolor*.7
   endfor
endif
;
; mark anticyclone points with ozone less than ambient-1sigma
;
for k=0L,11L do begin
    rtheta=th(k)
    index=where(rtheta eq th)
    itheta=index(0)
    if itheta eq 0 then dtheta=(th(0)-th(1))/2.0
    if itheta eq nth-1L then dtheta=(th(nth-2)-th(nth-1))/2.0
    if itheta gt 0 and itheta le nth then dtheta=(th(itheta-1)-th(itheta+1))/2.0
index=where(athmean_save le th(0) and athmean_save ge th(11) and $
            athmean_save gt rtheta-dtheta and athmean_save le rtheta+dtheta)
if index(0) eq -1L then goto,jumplev2
athmean=athmean_save(index(0))
ao3thresh=ao3mean_save(index(0))-1.5*ao3sigma_save(index(0))
    index=where(hth_save le th(0) and hth_save ge th(11) and $
                hth_save gt rtheta-dtheta and hth_save le rtheta+dtheta and $
                ho3_save lt ao3thresh,npt)
if index(0) eq -1L then goto,jumplev2
oplot,ho3_save(index),hth_save(index),psym=8,color=mcolor*.1,symsize=1.5

jumplev2:
endfor

      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim MLS_LOPs/xz_o3_mls+polar_nh_'+lfile+'_v2.ps -rotate -90 '+$
               ' MLS_LOPs/xz_o3_mls+polar_nh_'+lfile+'_v2.jpg'
      endif
      icount=icount+1L
goto,jump
end
