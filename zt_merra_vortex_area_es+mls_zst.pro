;
; save MERRA vortex area as a function of altitude and time about ES days zeros
; overplot MLS stratopause height and height of max temp
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

sver='v3.3'
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=9L & lstdy=1L & lstyr=2009L 
ledmn=5L & leddy=1L & ledyr=2010L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)

lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,39
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
xorig=[0.15]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirm='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
;restore, '/Users/harvey/Desktop/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
;result=size(MAXHEIGHTTHETA)
;nevents=result(1)

mlsesdates=['20060130','20090205','20120130','20130123']
nevents=n_elements(mlsesdates)

for iES = 0L, nevents - 1L do begin
    sevent=strtrim(strcompress(string(format='(I3.2)',ies+1)),2)
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_vortex_area_ES_event_'+sevent+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

    esdate0=mlsesdates(ies)
    iyr=long(strmid(esdate0,0,4))
    imn=long(strmid(esdate0,4,2))
    idy=long(strmid(esdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
koff=60
    jday0=jday-koff
    jday1=jday+koff
    CALDAT, jday0, lstmn ,lstdy , lstyr
    CALDAT, jday1, ledmn ,leddy , ledyr

lstday=0L & ledday=0L
if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
;goto,quick
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=long(ledday-lstday+1L)
;
; Compute initial Julian date
;
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
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
print,sdate
;
; read daily file
;
        dum=findfile(dir+sdate+'.nc3')
        if dum ne '' then ncfile0=dir+sdate+'.nc3'
        rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        if iflag ne 0L then goto,jump
        tmp2=0.*p2
        for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286

      if icount eq 0 then begin
         sfile=strarr(kday)
         area_zt_nc4=fltarr(kday,nth)
         ztavg=fltarr(kday)	; height of polar cap average stratopause
         ztmax=fltarr(kday)	; height of maximum temperature anywhere in the polar cap
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
      sfile(icount)=sdate
;
; loop over theta
;
      for thlev=0,nth-1 do begin
          mark1=transpose(mark2(*,*,thlev))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=100.*total(area(index))/hem_area
      endfor
;
; MLS
;
      dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) eq '' then goto,jumpday
      restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
      restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
      nz=n_elements(altitude)
      mprof=n_elements(longitude)
      mlev=n_elements(altitude)
      mlat=latitude
      mlon=longitude
      mtemp=temperature
      mpress=pressure
;
; retain NH polar cap
;

      index=where(mlat ge 70.)
      if index(0) eq -1L then goto,jumpday
      mlat=reform(mlat(index))
      mlon=reform(mlon(index))
      mtemp=reform(mtemp(index,*))
      mpress=reform(mpress(index,*))
      mtheta=mtemp*(1000./mpress)^0.286
      index=where(mtemp lt 0.)
      if index(0) ne -1L then mtheta(index)=-99.
;
; polar cap mean and max T profile
;
      thavg=fltarr(mlev)
      tavg=fltarr(mlev)
      tmax=fltarr(mlev)
      for k=0L,mlev-1L do begin
          good=where(mtemp(*,k) gt 0.)
          if good(0) ne -1L then thavg(k)=mean(mtheta(good,k))
          if good(0) ne -1L then tavg(k)=mean(mtemp(good,k))
          if good(0) ne -1L then tmax(k)=max(mtemp(good,k))
      endfor
;erase
;plot,tavg,altitude,xrange=[150,350],color=0
;oplot,tmax,altitude,color=0

zindex=where(altitude gt 15. and altitude lt 100.)
zprof=altitude(zindex)
thprof=thavg(zindex)
tavgprof=tavg(zindex)
tmaxprof=tmax(zindex)
index=where(tavgprof eq max(tavgprof))
ztavg(icount)=thprof(index)

;oplot,[max(tavgprof),max(tavgprof)],[zprof(index),zprof(index)],psym=8,color=0
index=where(tmaxprof eq max(tmaxprof))
ztmax(icount)=thprof(index)

;oplot,[max(tmaxprof),max(tmaxprof)],[zprof(index),zprof(index)],psym=4,symsize=2,color=0

      jumpday:
      icount=icount+1L
goto,jump

saveit:
;
; plot altitude-time series of Arctic vortex area
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
if long(min(yy(index))) lt 90L then y1='20'+string(FORMAT='(I2.2)',min(yy(index)))
if long(min(yy(index))) gt 90L then y1='19'+string(FORMAT='(I2.2)',min(yy(index)))
if long(max(yy(index))) lt 90L then y2='20'+string(FORMAT='(I2.2)',max(yy(index)))
if long(max(yy(index))) gt 90L then y2='19'+string(FORMAT='(I2.2)',max(yy(index)))
;
; save file
;
save,file='vortex_area_ES_'+esdate0+'_merra.sav',area_zt_nc4,th,sfile,y1,y2,koff,ztavg,ztmax
index=where(area_zt_nc4 le 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.

plotit:
restore,'vortex_area_ES_'+esdate0+'_merra.sav'
kday=n_elements(SFILE)
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/merra_vortex_area_ES_event_'+sevent+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot zt vortex area
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=18
col1=1+indgen(nlvls)*icolmax/nlvls
level=2.5*findgen(nlvls)
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
contour,area_zt_nc4,-1.*koff+findgen(kday),th,color=0,xtitle='Days From ES Onset',thick=6,yrange=[500.,8000.],/noeras,ytitle='Theta (K)',/cell_fill,$
     title='MLS ES event '+sevent+' ('+esdate0+')',levels=level,c_color=col1
contour,area_zt_nc4,-1.*koff+findgen(kday),th,levels=level,c_color=mcolor,/follow,/overplot,c_labels=1+0*level,min_value=0
plots,0,500
plots,0,8000.,/continue,color=0,thick=3
oplot,-1.*koff+findgen(kday),ztavg,psym=8,color=mcolor*.9
oplot,-1.*koff+findgen(kday),ztmax,psym=8,color=0
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='MERRA Vortex Area (% NH)',charthick=2,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ../Figures/merra_vortex_area_ES_event_'+sevent+'.ps -rotate -90 ../Figures/merra_vortex_area_ES_event_'+sevent+'.jpg'
endif

endfor	; loop over ES events
end
