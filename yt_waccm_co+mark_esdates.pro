;
; WACCM version
; contour CO and the vortex edge as a function of latitude and time
; +/- 30 days around all ES events
;
@stddat
@kgmt
@ckday
@kdate

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
dir='/Volumes/Data/WACCM/WACCM4/mee00fpl_FW2/mee00fpl_FW2.cam2.h3.Year'
dir2='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
mlsesdates=[20041218,20060130,20151221,20231223,20261210,20320204,20331220,20390226,20420104]
mlsesdates=strcompress(mlsesdates,/remove_all)

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
      if ndays gt ledday then goto,plotit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
;
; read WACCM
;
if smn eq '02' and sdy eq '29' then goto,jump
print,sdate
;
; read daily files
;
        print,'reading ',strmid(sdate,2,2),smn+sdy
        ncfile0=dir+strmid(sdate,2,2)+'_'+smn+sdy+'_CO.sav'
        ncfile1=dir+strmid(sdate,2,2)+'_'+smn+sdy+'_TP.sav'
        restore,ncfile0
        restore,ncfile1
        ncfile0=dir2+sdate+'_3D_dyn.nc3'
        print,ncfile0
        ncid=ncdf_open(ncfile0)
        result0=ncdf_inquire(ncid)
        for idim=0,result0.ndims-1 do begin
            ncdf_diminq,ncid,idim,name,dim
            if name eq 'number_of_latitudes' then nr=dim
            if name eq 'number_of_longitudes' then nc=dim
            if name eq 'number_of_levels' then nth=dim
        endfor
        for ivar=0,result0.nvars-1 do begin
            result=ncdf_varinq(ncid,ivar)
            if result.name eq 'latitude' or result.name eq 'longitude' or result.name eq 'theta' or result.name eq 'MARK' or $
               result.name eq 'U' or result.name eq 'V' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
            if result.name eq 'latitude' then alat=data
            if result.name eq 'longitude' then alon=data
            if result.name eq 'theta' then th=data
            if result.name eq 'MARK' then mark2=data
            if result.name eq 'U' then u2=data
            if result.name eq 'V' then v2=data
;           print,ivar,result.name,min(data),max(data)
        endfor
        ncdf_close,ncid
        speed2=sqrt(u2^2+v2^2)
;
; declare arrays on first day of event
;
      if icount eq 0 then begin
         nl=n_elements(lev)
         nr=n_elements(lat)
         sdate_all=strarr(kday)
         zindex=where(th eq 500. or th eq 1000. or th eq 2000. or th eq 3000. or th eq 4000.,nth2)

         ytmark=fltarr(kday,nr,nth2)
         ytspeed=fltarr(kday,nr,nth2)
         ytco=fltarr(kday,nr,nth2)
         yth2o=fltarr(kday,nr,nth2)
         kcount=1L
      endif
      sdate_all(icount)=sdate
;
; zonal mean WACCM theta
;
      wtheta=0.*tp
      for j=0L,nr-1L do wtheta(j,*)=tp(j,*)*(1000./lev)^0.286
;
; 2d latitude array
;
      mlat2=0.*tp
      for k=0L,nl-1L do mlat2(*,k)=lat
;
; use wtheta to interpolate WACCM to th levels
;
      for kk=0,nth2-1L do begin
          ilev=zindex(kk)
          rlev=th(ilev)
          slev=strcompress(long(rlev),/remove_all)+'K'
          print,slev

          for k=1L,nl-1L do begin
          km1=k-1L				; waccm theta arrays are top down
          for j=0L,nr-1L do begin
          if wtheta(j,km1) gt rlev and wtheta(j,k) le rlev then begin
             zscale=(wtheta(j,km1)-rlev)/(wtheta(j,km1)-wtheta(j,k))
             ytco(icount,j,kk)=co(j,km1)-zscale*(co(j,km1)-co(j,k))
;print,wtheta(j,k),rlev,wtheta(j,km1),zscale
;print,co(j,k), ytco(icount,j,kk),co(j,km1)
;stop
          endif
          endfor
          endfor
          speed1=transpose(speed2(*,*,ilev))
          mark1=transpose(mark2(*,*,ilev))
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
save,file='yt_waccm_co+mark_'+yearlab+'.sav',alat,kday,sdate_all,ytco,yth2o,ytmark,ytspeed,slev,th,zindex,nth2

quick:
restore,'yt_waccm_co+mark_'+yearlab+'.sav'
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
         device,/landscape,bits=8,filename='yt_waccm_co+mark_'+yearlab+'_'+slev+'.ps'
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
plotarray=ytco*1.e6
stitle=yearlab+' CO, Arctic Vortex, and Wind Speed'
cbartitle=slev+' WACCM CO (ppmv)'
;if rlev le 2000. then begin
;   plotarray=yth2o
;   stitle=yearlab+' H2O, Arctic Vortex, and Wind Speed'
;   cbartitle=slev+' WACCM H2O (ppmv)'
;endif
index=where(plotarray ne -9999.)
imin=min(plotarray(index))
;if rlev le 2000. then imin=5.
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
   spawn,'convert -trim yt_waccm_co+mark_'+yearlab+'_'+slev+'.ps -rotate -90 '+$
         'yt_waccm_co+mark_'+yearlab+'_'+slev+'.jpg'
endif

endfor	; loop over altitude
endfor	; loop over ES event
end
