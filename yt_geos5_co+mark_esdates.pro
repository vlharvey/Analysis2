;
; contour CO and the vortex edge as a function of latitude and time
; +/- 30 days around all ES events
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

sver='v3.3'

loadct,39
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
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
mlsesdates=['20060130','20090205','20120130','20130123']
;sabesdates=['20060128','20090205','20120128','20130123']

nevent=n_elements(mlsesdates)
for ievent=0,nevent-1L do begin
    esdate0=mlsesdates(ievent)
    iyr=long(strmid(esdate0,0,4))
    imn=long(strmid(esdate0,4,2))
    idy=long(strmid(esdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
    jday0=jday-30
    jday1=jday+30
    CALDAT, jday0, lstmn ,lstdy , lstyr
    CALDAT, jday1, ledmn ,leddy , ledyr

lstday=0L & ledday=0L
if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
goto,quick
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
; construct 2d MLS latitude
;
      mlat2=0.*mco
      for i=0L,mlev-1L do mlat2(*,i)=mlat
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
;
; declare arrays on first day of event
;
      if icount eq 0 then begin
         sdate_all=strarr(kday)
         x2d=fltarr(nc,nr)
         y2d=fltarr(nc,nr)
         for i=0L,nc-1L do y2d(i,*)=alat
         for j=0L,nr-1L do x2d(*,j)=alon
         zindex=where(th eq 500. or th eq 1000. or th eq 2000. or th eq 3000. or th eq 4000.,nth2)
         ytco=-9999.+fltarr(kday,nr,nth2)
         yth2o=-9999.+fltarr(kday,nr,nth2)
         ytspeed=-9999.+fltarr(kday,nr,nth2)
         ytmark=-9999.+fltarr(kday,nr,nth2)
      endif
      sdate_all(icount)=sdate
;
; loop over theta
;
      for kk=0,nth2-1L do begin
          ilev=zindex(kk)
          rlev=th(ilev)
          slev=strcompress(long(rlev),/remove_all)+'K'
          print,slev
          kindex=where(mlat2 ge 0. and abs(mtheta-rlev) le 50. and mco ne -99. and mh2o ne -99.,mprof)
          if kindex(0) eq -1L then goto,skipit
          codata=mco(kindex)*1.e6
          h2odata=mh2o(kindex)*1.e6
          ydata=mlat2(kindex)
;
; zonal mean MLS
;
          dlat=alat(1)-alat(0)
          for iprof=0L,mprof-1L do begin
          for j=0L,nr-1L do begin
              index=where(ydata ge alat(j)-dlat/2. and ydata lt alat(j)+dlat/2.)
              if index(0) ne -1L then begin
                 ytco(icount,j,kk)=mean(codata(index))
                 yth2o(icount,j,kk)=mean(h2odata(index))
              endif
          endfor
          endfor
          pv1=transpose(pv2(*,*,ilev))
          sf1=transpose(sf2(*,*,ilev))
          speed1=transpose(speed2(*,*,ilev))
          mark1=transpose(marknew2(*,*,ilev))
;
; bin max marker and CO in lat
;
          for j=0L,nr-1L do begin
              ytmark(icount,j,kk)=mean(mark1(*,j))
              ytspeed(icount,j,kk)=mean(speed1(*,j))
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
save,file='yt_geos5_co+mark_'+yearlab+'.sav',alat,kday,sdate_all,ytco,yth2o,ytmark,ytspeed,slev,th,zindex,nth2

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
; loop over theta
;
      for kk=0,nth2-1L do begin
          ilev=zindex(kk)
          rlev=th(ilev)
          slev=strcompress(long(rlev),/remove_all)+'K'

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
         device,/landscape,bits=8,filename='yt_geos5_co+mark_'+yearlab+'_'+slev+'.ps'
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
plotarray=ytco
stitle=yearlab+' CO, Arctic Vortex, and Wind Speed'
cbartitle=slev+' MLS CO (ppmv)'
if rlev le 2000. then begin
   plotarray=yth2o
   stitle=yearlab+' H2O, Arctic Vortex, and Wind Speed'
   cbartitle=slev+' MLS H2O (ppmv)'
endif
index=where(plotarray ne -9999.)
imin=min(plotarray(index))
if rlev le 2000. then imin=5.
imax=max(plotarray(index))
iint=(imax-imin)/float(nlvls)
level=imin+iint*findgen(nlvls)
contour,smooth(plotarray(*,*,kk),3),1.+findgen(kday),alat,levels=level,color=0,c_color=col1,/noeras,charsize=2,charthick=2,$
        xrange=[1,kday],yrange=[30,80],/fill,title=stitle,min_value=-9999.,$
        ytitle='Latitude',xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
;contour,smooth(ytco(*,*,kk),3),1.+findgen(kday),alat,levels=level,color=0,/noeras,/follow,/overplot,min_value=-9999.
contour,smooth(ytmark(*,*,kk),3),1.+findgen(kday),alat,/overplot,levels=[0.1,0.2,0.4,0.6,0.8,1.0],color=0,thick=10,/noeras,/follow,min_value=-9999.
loadct,0
contour,smooth(ytspeed(*,*,kk),3),1.+findgen(kday),alat,/overplot,levels=30.+10*findgen(20),color=200,thick=6,/noeras,/follow,min_value=-9999.
plots,31,30
plots,31,80,/continue,thick=5,color=100
loadct,39
omin=min(level)
omax=max(level)
set_viewport,xmn,xmx,ymn-cbaryoff,ymn-cbaryoff+cbarydel
!type=2^2+2^3+2^6
plot,[omin,omax],[0,0],yrange=[0,10],charsize=2,charthick=2,$
      xrange=[omin,omax],xtitle=cbartitle,/noeras,xstyle=1,color=0
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
   spawn,'convert -trim yt_geos5_co+mark_'+yearlab+'_'+slev+'.ps -rotate -90 '+$
         'yt_geos5_co+mark_'+yearlab+'_'+slev+'.jpg'
endif

endfor	; loop over altitude
endfor	; loop over ES event
end
