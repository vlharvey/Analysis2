;
; contour CO and the vortex edge as a function of latitude and time
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

sver='v3.3'

loadct,38
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=1000
nydim=700
xorig=[.1,0.1,0.55,0.55]
yorig=[.25,.1,.55,.1]
xlen=0.8
ylen=0.6
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
!noeras=1
dirm='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
lstmn=1L & lstdy=1L & lstyr=2006L
ledmn=3L & leddy=15L & ledyr=2006L
lstday=0L & ledday=0L
if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
goto,quick
;
; get date range
;
print, ' '
print, '      GEOS-5 Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=(ledday-lstday+1L)
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
      if ndays gt ledday then goto,plotit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
;
; read MLS data
;
      dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipit
      restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
      restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
      restore,dirm+'co_mls_'+sver+'_'+sdate+'.sav'              ; mix
      nz=n_elements(altitude)
      nthlev=n_elements(thlev)
      mprof=n_elements(longitude)
      mlev=n_elements(altitude)
      muttime=time
      mlat=latitude
      mlon=longitude
      bad=where(mask eq -99.)
      if bad(0) ne -1L then mix(bad)=-99.
      good=where(mix ne -99.)
      if good(0) eq -1L then goto,skipit
      mco=mix
      restore,dirm+'h2o_mls_'+sver+'_'+sdate+'.sav'              ; water vapor mix
      bad=where(mask eq -99.)
      if bad(0) ne -1L then mix(bad)=-99.
      good=where(mix ne -99.)
      if good(0) eq -1L then goto,skipit
      mh2o=mix
      mtemp=temperature
      mpress=pressure
;
; eliminate bad uttimes and SH
;
      index=where(muttime gt 0.)
      if index(0) eq -1L then goto,skipit
      muttime=reform(muttime(index))
      mlat=reform(mlat(index))
      mlon=reform(mlon(index))
      mtemp=reform(mtemp(index,*))
      mpress=reform(mpress(index,*))
      mco=reform(mco(index,*))
      mh2o=reform(mh2o(index,*))
      mtheta=mtemp*(1000./mpress)^0.286
      index=where(mtemp lt 0.)
      if index(0) ne -1L then mtheta(index)=-99.
;
; construct 2d MLS arrays to match CO
;
      mpress2=mpress
      mtime2=0.*mco
      mlat2=0.*mco
      mlon2=0.*mco
      for i=0L,mlev-1L do begin
          mtime2(*,i)=muttime
          mlat2(*,i)=mlat
          mlon2(*,i)=mlon
      endfor
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,skipit
;
; read new vortex
;
      ncid=ncdf_open(dir+sdate+stimes(0)+'nc5')
      marknew2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,marknew2
      ncdf_close,ncid
      index=where(marknew2 lt 0.)
      if index(0) ne -1L then marknew2(index)=-1.0
      speed2=sqrt(u2^2+v2^2)
      if icount eq 0 then begin
         sdate_all=strarr(kday)
         x2d=fltarr(nc,nr)
         y2d=fltarr(nc,nr)
         for i=0L,nc-1L do y2d(i,*)=alat
         for j=0L,nr-1L do x2d(*,j)=alon
         ytco=-9999.+fltarr(kday,nr)
         yth2o=-9999.+fltarr(kday,nr)
         ytspeed=-9999.+fltarr(kday,nr)
         ytmark=-9999.+fltarr(kday,nr)
      endif
      sdate_all(icount)=sdate
;
; loop over theta
;
index=where(th eq 3000.)
ilev=index(0)
      for kk=ilev,ilev do begin
          rlev=th(kk)
          slev=strcompress(long(rlev),/remove_all)+'K'
          kindex=where(mlat2 ge 0. and abs(mtheta-rlev) le 50. and mco ne -99. and mh2o ne -99.,mprof)
          if kindex(0) eq -1L then goto,skipit
          codata=mco(kindex)*1.e6
          h2odata=mh2o(kindex)*1.e6
          xdata=mlon2(kindex)
          ydata=mlat2(kindex)
;
; zonal mean MLS
;
          dlat=alat(1)-alat(0)
          for iprof=0L,mprof-1L do begin
          for j=0L,nr-1L do begin
              index=where(ydata ge alat(j)-dlat/2. and ydata lt alat(j)+dlat/2.)
              if index(0) ne -1L then begin
                 ytco(icount,j)=mean(codata(index))
                 yth2o(icount,j)=mean(h2odata(index))
              endif
          endfor
          endfor
          pv1=transpose(pv2(*,*,kk))
          sf1=transpose(sf2(*,*,kk))
          speed1=transpose(speed2(*,*,kk))
          mark1=transpose(marknew2(*,*,kk))
;
; bin max marker and CO in lat
;
          for j=0L,nr-1L do begin
              ytmark(icount,j)=mean(mark1(*,j))
              ytspeed(icount,j)=mean(speed1(*,j))
          endfor	; loop over latitude
     endfor		; loop over altitude
skipit:
            icount=icount+1L
;     endfor    ; loop over 4 daily times

goto,jump
plotit:
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
save,file='yt_geos5_co+mark_'+yearlab+'.sav',alat,kday,sdate_all,ytco,yth2o,ytmark,ytspeed,slev

quick:
restore,'yt_geos5_co+mark_'+yearlab+'.sav'
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
; save postscript version
;
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename='yt_geos5_co+mark_'+sdate+'_'+slev+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.thick=2.0                   ;Plotted lines twice as thick
         !p.charsize=2.0
      endif

erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
index=where(ytco ne -9999.)
imin=min(ytco(index))
imax=max(ytco(index))
iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
;ytco=smooth(ytco,3)
;ytmarkco=smooth(ytmark,3)
contour,ytco,1.+findgen(kday),alat,levels=level,color=0,c_color=col1,/noeras,charsize=2,charthick=2,$
        xrange=[1,kday],yrange=[25,85],/fill,title=yearlab+' CO, Arctic Vortex, and Wind Speed',min_value=-9999.,yticks=6,$
        ytitle='Latitude',xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
contour,ytco,1.+findgen(kday),alat,levels=level,color=0,/noeras,/follow,/overplot,min_value=-9999.
contour,ytmark,1.+findgen(kday),alat,/overplot,levels=[0.1,0.3,0.5,0.7,0.9],color=0,thick=10,/noeras,/follow,min_value=-9999.
loadct,0
contour,ytspeed,1.+findgen(kday),alat,/overplot,levels=30.+20*findgen(20),color=200,thick=6,/noeras,/follow,min_value=-9999.
loadct,39
omin=min(level)
omax=max(level)
set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
!type=2^2+2^3+2^6
plot,[omin,omax],[0,0],yrange=[0,10],charsize=2,charthick=2,$
      xrange=[omin,omax],xtitle=slev+' MLS CO (ppmv)',/noeras,xstyle=1,color=0
ybox=[0,10,10,0,0]
x1=omin
dx=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim yt_geos5_co+mark_'+sdate+'_'+slev+'.ps -rotate -90 '+$
         'yt_geos5_co+mark_'+sdate+'_'+slev+'.jpg'
endif
end
