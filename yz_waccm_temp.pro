;
; zonal mean temperature WACCM
;
@stddat
@kgmt
@ckday
@kdate

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
diro='/Users/harvey/Harvey_etal_2014/Post_process/'
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
esdates=[$	; SSW dates
20030208,$
20061220,$
20081209,$
20090104,$
20150213,$
20260214,$
20271219,$
20300110,$
20331231,$
20401219,$
20411227]
;20041218,20060130,20151221,20231223,20261210,20320204,20331220,20390226,20420104]	; ES dates
esdates=strcompress(esdates,/remove_all)
;result=size(MAXHEIGHTTHETA)
;nevents=result(1)
nevents=n_elements(esdates)

for iES = 0L, nevents - 1L do begin
    sevent=strtrim(strcompress(string(format='(I3.2)',ies+1)),2)
    esdate0=esdates(ies)
    iyr=long(strmid(esdate0,0,4))
    imn=long(strmid(esdate0,4,2))
    idy=long(strmid(esdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
koff=30
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
      if ndays gt ledday then stop,' Normal termination condition'
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
if smn eq '02' and sdy eq '29' then goto,jump
print,sdate
;
; read daily file
;
        ncfile0=dir+sdate+'_3D_dyn.nc3'
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
            if result.name eq 'latitude' or result.name eq 'longitude' or result.name eq 'theta' or result.name eq 'P' or result.name eq 'MARK' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
            if result.name eq 'latitude' then alat=data
            if result.name eq 'longitude' then alon=data
            if result.name eq 'theta' then th=data
            if result.name eq 'P' then p2=data
            if result.name eq 'MARK' then mark2=data
;           print,ivar,result.name,min(data),max(data)
        endfor
        ncdf_close,ncid
        temp2=0.*p2
        for k=0L,nth-1L do temp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^.286

      if icount eq 0 then begin
         ztavg=fltarr(kday)     ; height of polar cap average stratopause
         ztmax=fltarr(kday)     ; height of maximum temperature anywhere in the polar cap
         yindex=where(alat ge 70.)
         sfile=strarr(kday)
         area_zt_nc4=fltarr(kday,nth)
         ztavg=fltarr(kday)     ; height of polar cap average stratopause
         ztmax=fltarr(kday)     ; height of maximum temperature anywhere in the polar cap
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
; polar cap mean and max T profile
;
      tzm=fltarr(nr,nth)
      tavg=fltarr(nth)
      tmax=fltarr(nth)
      for k=0L,nth-1L do begin
          tavg(k)=mean(temp2(yindex,*,k))
          tmax(k)=max(temp2(yindex,*,k))
          for j=0L,nr-1L do tzm(j,k)=mean(temp2(j,*,k))
      endfor

good=where(th le 10000. and th ge 500.)
tavgprof=tavg(good)
tmaxprof=tmax(good)
thprof=th(good)
index=where(tavgprof eq max(tavgprof))
ztavg(icount)=thprof(index)
index=where(tmaxprof eq max(tmaxprof))
ztmax(icount)=thprof(index)

;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/waccm_vortex_area_ssw_event_'+sevent+'.ps'
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
nlvls=26
col1=1+indgen(nlvls)*icolmax/nlvls
level=150.+5.*findgen(nlvls)
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
contour,tzm,alat,th,color=0,xtitle='Days From SSW Onset',thick=6,yrange=[300,10000],/noeras,ytitle='Theta (K)',/cell_fill,$
     title='WACCM SSW event '+sevent+' ('+esdate0+')',levels=level,c_color=col1
contour,tzm,alat,th,levels=level,c_color=mcolor,/follow,/overplot,c_labels=1+0*level,min_value=0
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='WACCM Temperature (K)',charthick=2,charsize=1.5
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
   spawn,'convert -trim ../Figures/waccm_vortex_area_ssw_event_'+sevent+'.ps -rotate -90 ../Figures/waccm_vortex_area_ssw_event_'+sevent+'.jpg'
endif

      jumpday:
      icount=icount+1L
goto,jump

endfor	; loop over SSW events
end
