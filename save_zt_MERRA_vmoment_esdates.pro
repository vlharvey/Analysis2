;-----------------------------------------------------------------------------------------------------------------------------
; Reads in MERRA data and plots number of cyclonic lobes in the NH on each level each day
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3
@range_ring
@vortexshape

re=40000000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0

px1a = .22
px1b = .73
px2a = .52
px2b = .95
py1a = .50
py1b = .95
py2a = .45
py2b = .66
py3a = .15
py3b = .35

SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.6,0.1,0.6,0.1,0.6]
yorig=[0.7,0.7,0.4,0.4,0.1,0.1]
xlen=0.25
ylen=0.25
cbaryoff=0.1
cbarydel=0.01

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=!p.color
icolmax=255
mcolor=icolmax
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nzdim,ysize=nydim,retain=2,colors=162
endif
days = 0
months = ['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
MONTH = ['04','05','06','07','08','09','10','11']
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
;
; Read ES day zeros
;
restore, '/Users/harvey/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
kcount=0L
result=size(MAXHEIGHTTHETA)
nevents=result(1)
for iES = 0L, nevents - 1L do begin
    sevent=string(format='(i2.2)',ies+1)
    sevent=strtrim(strcompress(string(format='(I3.2)',ies+1)),2)
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_zt_nvort_ES_event_'+sevent+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

    icount=0L
    ndays=61
    for iday =0L, ndays-1L do begin
        if iday lt 30L then sday=string(format='(I3.2)',iday-30)
        if iday ge 30L then sday=string(format='(I2.2)',iday-30)
        sdate=esdate(iday+60*ies)
        if sdate eq '' then goto,jumpday        ; missing day
print,iday-30,' ',sdate
;
; read daily file
;
        dum=findfile(dir+sdate+'.nc3')
        if dum ne '' then ncfile0=dir+sdate+'.nc3'
        rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        if iflag ne 0L then goto,jumpday
        tmp2=0.*p2
        for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286

marker_USLM = make_array(nc,nr,nth)
for k=0,nth-1 do begin
  marker_USLM(*,*,k) = transpose(mark2(*,*,k))
endfor
shape = vortexshape(marker_USLM, alat, alon)
centroid=shape.nhcentroid
centroidx=reform(centroid(0,*))
centroidy=reform(centroid(1,*))
axis=shape.axis
majoraxis=reform(axis(0,*))
minoraxis=reform(axis(1,*))

;
; wrap around point
;
        lons = fltarr(nc+1L)
        lons[0:nc-1L] = alon
        lons[nc] = alon[0L]
;
; coordinate transformation
;
        nr2=nr/2
        x2d=fltarr(nc+1,nr/2)
        y2d=fltarr(nc+1,nr/2)
        for i=0,nc do y2d(i,*)=alat(nr/2:nr-1)
        for j=0,nr/2-1 do x2d(*,j)=lons
        xcn=fltarr(nc+1,nr2)
        ycn=fltarr(nc+1,nr2)
        for j=nr2,nr-1 do begin
            ANG = (90. - alat(j)) * RADG * 0.5
            FACTOR = TAN(ANG) * FAC20
            for i=0,nc do begin
                THETA = (lons(i) - 90.) * RADG
                xcn(i,j-nr2) = FACTOR * COS(THETA)
                ycn(i,j-nr2) = FACTOR * SIN(THETA)
            endfor
        endfor
;
; theta below stratopause max height
;
        if maxheighttheta[ies,iday] le 0. or maxheighttheta[ies,iday] gt 10000. then maxheighttheta[ies,iday]= maxheighttheta[ies,iday-1]
        xtheta = th - maxheighttheta[ies,iday]
        x = where(th ge 350. and th le 10000., nz)
        if x(0) eq -1L then goto,jumpday
        th2=th(x)
        th2r=reverse(th(x))
        mark2=mark2(*,*,x)
        sf2=sf2(*,*,x)
        z2=z2(*,*,x)
        xphase=0*th2
        yphase=0*th2
;
; declare nvortex array on first day
; 
        if iday eq 0L then begin
           nvort_all=fltarr(ndays,nz)
           xcentroid_all=fltarr(ndays,nz)
           ycentroid_all=fltarr(ndays,nz)
           ellip_all=fltarr(ndays,nz)
           sdate_all=strarr(ndays)
        endif
        sdate_all(iday)=sdate
;
; max theta below stratopause max
;
        thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(5000.) - min(alog(th2),/nan))  *  254.
        maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(5000.)))*254.
        if maxheighttheta[ies,iday] gt 5000. then begin
           thlevel = reverse(alog(th2) - min(alog(th2),/nan))/(alog(maxheighttheta[ies,iday]) - min(alog(th2),/nan))*254.
           maxthetalevel = (alog(maxheighttheta[ies,iday])/(alog(maxheighttheta[ies,iday])))*254.
        endif
;
; plot
;
erase
        altarray=fltarr(nz)
	for ii = 0L, nz - 1L do begin
        map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,position = [0.25,0.45,0.75,0.95],color=0,$
                title='MERRA ES event ' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)+', Day '+ sday+' ('+sdate+')'

            marker=fltarr(nc+1L,nr2)
            marker[0L:nc-1L,*] = transpose(reform(mark2[nr2:nr-1,*,nz-1L-ii]))
            marker[nc,*] = marker(0,*)
            sf=fltarr(nc+1L,nr2)
            sf[0L:nc-1L,*] = transpose(reform(sf2[nr2:nr-1,*,nz-1L-ii]))
            sf[nc,*] = sf(0,*)
            z=fltarr(nc+1L,nr2)
            z[0L:nc-1L,*] = transpose(reform(z2[nr2:nr-1,*,nz-1L-ii]))
            z[nc,*] = z(0,*)
;
; contour vortex edge
;
             contour,marker,lons,alat(nr2:nr-1),levels=[.1],/overplot,color=thlevel[ii],thick=8

;
            smin=min(sf)
            smax=max(sf)
            slevel=smin+((smax-smin)/39.)*findgen(40)
;           contour,sf,lons,alat(nr2:nr-1),levels=slevel,/overplot,color=0
;
; superimpose center of vortex
;
            index=where(marker gt 0.)
            if index(0) ne -1L then begin
               skm=strcompress(long(mean(z(index))),/remove_all)+' km'
               if ii mod 2 eq 0 then xyouts,px1a-0.05,py1a+ii*((py1b-py1a)/nz),skm,color=thlevel(ii),/normal,charsize=1.25,charthick=3
               altarray(nz-1-ii)=mean(z(index))
;
xm0=centroidx(nz-1L-ii)
ym0=centroidy(nz-1L-ii)
oplot,[xm0,xm0],[ym0,ym0],psym=8,color=0,symsize=3
xcentroid_all(iday,ii)=xm0
ycentroid_all(iday,ii)=ym0
;
; histograms of vortex latitude and longitude
;
!type=2^2+2^3
py1=fltarr(nc)
px1=fltarr(nr)
for i=0L,nc-1L do begin
    index=where(marker(i,*) gt 0.,nn)
    py1(i)=float(nn)
endfor
for j=0L,nr-1L do begin
    index=where(mark2(j,*,nz-1L-ii) gt 0.,nn)
    px1(j)=float(nn)
endfor
;
; are there two cyclonic vortices?
;
n0=findgen(nc)
n1=1.+findgen(nc)
vortlon=0.*alon
index=where(py1 ne 0.)
if index(0) ne -1L then vortlon(index)=1.
index=where(abs(vortlon(n0)-vortlon(n1)) gt 0.,nv)
nextra=1
if nv eq 0L then nv=2			; circumpolar has no zeros. set to 2 to get 1 vortex
index=where(vortlon eq 1.)
if min(alon(index)) eq min(alon) and max(alon(index)) ne max(alon) then nextra=0        ; GM edge
if min(alon(index)) ne min(alon) and max(alon(index)) eq max(alon) then nextra=0        ; GM edge
if nv gt 2L then begin
   nextra=0.5*nv		; each vortex results in 2 edge points - unless it lies exactly along the GM
;  if min(alon(index)) eq min(alon) and max(alon(index)) ne max(alon) then nextra=0        ; GM edge
;  if min(alon(index)) ne min(alon) and max(alon(index)) eq max(alon) then nextra=0        ; GM edge
endif
nv=round(nv-nextra)
print,'number of cyclonic lobes ',nv

;set_viewport,0.1,0.45,0.1,0.4
;plot,alat,100.*px1/float(nc),color=0,psym=8,xrange=[0.,90.],yrange=[0.,100],ytitle='% Lons f(y)'
;oplot,alat,100.*px1/float(nc),color=thlevel(ii),psym=8
;set_viewport,0.55,0.9,0.1,0.4
;plot,alon,100.*py1/float(nr),color=0,psym=8,yrange=[0.,60],ytitle='% Lats f(x)'
;oplot,alon,100.*py1/float(nr),color=thlevel(ii),psym=8
;xyouts,300.,50.,strcompress(long(nv)),charsize=2,charthick=2,color=0,/data
nvort_all(iday,ii)=nv
;if nv gt 1 then stop

            endif

	endfor		; loop over altitude
;
set_viewport,0.3,0.7,0.1,0.4
ellipticity=minoraxis/majoraxis
;
; set moments to NaN where there are multiple vortices
;
index=where(nvort_all(iday,*) gt 1.)
if index(0) ne -1L then begin
ellipticity(index)=0./0.
xcentroid_all(iday,index)=0./0.
ycentroid_all(iday,index)=0./0.
endif
ellip_all(iday,*)=reverse(ellipticity)

plot,ellipticity,th,color=0,thick=4,xrange=[0.,1.0],yrange=[400.,5000.]
oplot,ellipticity,th,color=0,psym=8
index=where(ellipticity eq min(ellipticity))
;print,'min ellipticity= ',ellipticity(index),th(index)

        jumpday:
        icount=icount+1L

endfor		; loop over days 0 to +30

plotit:

altarray=reverse(altarray)
        x2d=fltarr(ndays,nz)
        y2d=fltarr(ndays,nz)
        for i=0,ndays-1 do y2d(i,*)=altarray
        for j=0,nz-1 do x2d(*,j)=-30+findgen(ndays)

;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_zt_nvort_ES_event_'+sevent+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

erase
set_viewport,0.2,0.9,0.3,0.7
contour,nvort_all,-30+findgen(ndays),altarray,/nodata,xrange=[-30,30],yrange=[10.,70.],$
        xtitle='Days since ES onset ('+strmid(sdate_all(30),0,4)+')',ytitle='Approximate Altitude (km)',/noeras,color=0,charsize=2,charthick=2
;index=where(nvort_all eq 1)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=0,symsize=0.5
;index=where(nvort_all eq 2)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=250,symsize=1
;index=where(nvort_all eq 3)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=200,symsize=1
;index=where(nvort_all eq 4)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=170,symsize=1
;index=where(nvort_all eq 5)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=50,symsize=1
contour,smooth(ellip_all,5,/edge_truncate),-30+findgen(ndays),altarray,/overplot,level=0.2*findgen(6),/follow,thick=5,color=0,/nodata
contour,smooth(ellip_all,5,/edge_truncate,/Nan),-30+findgen(ndays),altarray,/overplot,level=0.2,/follow,thick=5,color=100           
contour,smooth(ellip_all,5,/edge_truncate,/Nan),-30+findgen(ndays),altarray,/overplot,level=0.4,/follow,thick=5,color=150
contour,smooth(ellip_all,5,/edge_truncate,/Nan),-30+findgen(ndays),altarray,/overplot,level=0.6,/follow,thick=5,color=200
contour,smooth(ellip_all,5,/edge_truncate,/Nan),-30+findgen(ndays),altarray,/overplot,level=0.8,/follow,thick=5,color=250

index=where(nvort_all gt 1)
if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=0,symsize=1


;imin=1.
;imax=5.
;ymnb=0.3 -cbaryoff
;ymxb=ymnb  +cbarydel
;set_viewport,0.2,0.9,ymnb,ymxb
;!type=2^2+2^3+2^6
;plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
;      xtitle='MERRA # Cyclonic Vortices',charthick=2,charsize=1.5
;ybox=[0,10,10,0,0]
;x1=imin
;col1=[0,250,200,170,50]
;nlvls=n_elements(col1)
;dx=(imax-imin)/float(nlvls)
;for j=0,nlvls-1 do begin
;xbox=[x1,x1,x1+dx,x1+dx,x1]
;polyfill,xbox,ybox,color=col1(j)
;x1=x1+dx
;endfor

        if setplot ne 'ps' then stop
        if setplot eq 'ps' then begin
           device, /close
           spawn,'convert -trim ../Figures/merra_zt_nvort_ES_event_'+sevent+'.ps -rotate -90 ../Figures/merra_zt_nvort_ES_event_'+sevent+'.png'
           spawn,'rm -f ../Figures/merra_zt_nvort_ES_event_'+sevent+'.ps'
        endif
;
; save file
;
save,filename='merra_zt_vmoment_ES_event_'+sevent+'.sav',nvort_all,ellip_all,xcentroid_all,ycentroid_all,th2r,sdate_all,altarray

endfor		; loop over ES events
end
