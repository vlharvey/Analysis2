;
; zonal mean temperature deviation from multi-year mean
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
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
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.12
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
;start_year=[2007,2008,2009,2010,2011,2012,2013,2014,2015]
;start_date=[-27, -21, -24, -24, -26, -27, -34, -28,-42]
;end_date=[66, 65, 61, 61, 64, 61, 64, 80]
;nyear=n_elements(start_year)

lstmn=8         ; NH
lstdy=1
lstyr=2004
ledmn=5         ; NH
leddy=11
ledyr=2015
lstday=0
ledday=0
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
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
; longitude grid
;
dx=15.
nc=long(360./dx)+1
longrid=dx*findgen(nc)
nr=91L
latgrid=-90.+2.*findgen(nr)
dy=latgrid(1)-latgrid(0)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition'
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
if iday ge 60 and iday le 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(6,21,long(syr))
if iday lt 60 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr)-1L)
if iday gt 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr))
      print,sdate,' ',dfs
;
; restore MLS on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
;
; apply mask
;
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      nlv=n_elements(altitude)
      tbar=fltarr(nr,nlv)
;
; compute mean temp around lat circle
;
    tbar=fltarr(nr,nlv)
    nbar=lonarr(nr,nlv)
    for ii=0L,n_elements(id)-1L do begin
        tmask_prof=reform(TEMPERATURE_MASK(ii,*))
        good=where(tmask_prof ne -99.,ngood)
        if good(0) ne -1L then begin
           ymean=latitude(ii)
           for j=0L,nr-1L do begin
               if ymean ge latgrid(j)-dy/2. and ymean lt latgrid(j)+dy/2. then begin
                  tbar(j,good)=tbar(j,good)+temperature(ii,good)
                  nbar(j,good)=nbar(j,good)+1L
               endif
           endfor
        endif
    endfor
    index=where(nbar gt 1.)
    if index(0) ne -1L then tbar(index)=tbar(index)/float(nbar(index))
;
; make multi-year mean and exclude current year
;
    tmean=0.*tbar
    nmean=0*nbar
    smon=strmid(sdate,4,2)
    sday=strmid(sdate,6,2)
    spawn,'ls '+mdir+'cat_mls_v3.3_????'+smon+sday+'.sav',cfiles
    spawn,'ls '+mdir+'tpd_mls_v3.3_????'+smon+sday+'.sav',tfiles
;
; extract years
; 
years=strarr(n_elements(cfiles))
for i=0L,n_elements(cfiles)-1L do begin
result=strsplit(cfiles(i),/extract,'/')
result2=strsplit(result(5),/extract,'.')
result3=strsplit(result2(1),/extract,'_')
years(i)=string(result3(1))
endfor
good=WHERE(STRMATCH(years, sdate) NE 1)
cfiles=cfiles(good)
tfiles=tfiles(good)

    for i=0,n_elements(cfiles)-1L do begin
        restore,cfiles(i)
        restore,tfiles(i)
        print,cfiles(i)
        index=where(temperature_mask eq -99.)
        if index(0) ne -1L then temperature(index)=-99.
        for ii=0L,n_elements(id)-1L do begin
            tmask_prof=reform(TEMPERATURE_MASK(ii,*))
            good=where(tmask_prof ne -99.,ngood)
            if good(0) ne -1L then begin
               ymean=latitude(ii)
               for j=0L,nr-1L do begin
                   if ymean ge latgrid(j)-dy/2. and ymean lt latgrid(j)+dy/2. then begin
                      tmean(j,good)=tmean(j,good)+temperature(ii,good)
                      nmean(j,good)=nmean(j,good)+1L
                   endif
               endfor
            endif
        endfor	; loop over profiles
    endfor	; loop over past years
    index=where(nmean gt 1.)
    if index(0) ne -1L then tmean(index)=tmean(index)/float(nmean(index))
;
; anomaly from multi-year mean
;
tanom=tbar-tmean
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='MLS_YZ_Anomalies/yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps'
   !p.charsize=1.55
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean temperature and CO
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=-20+4.*findgen(11)
nlvls=n_elements(col2)
index=where(tanom eq 0.0)
if index(0) ne -1L then tanom(index)=0./0.
contour,tanom,latgrid,altitude,/noeras,xrange=[-90,90],yrange=[10.,90.],ytitle='Altitude (km)',charsize=1.5,color=0,xtitle='Latitude',xticks=6,charthick=2,title=sdate+' DFS='+strcompress(dfs,/r),$
        levels=level,c_color=col2,/cell_fill
index=where(level gt 0.)
contour,tanom,latgrid,altitude,/noeras,levels=level(index),color=0,/follow,/overplot
index=where(level lt 0.)
contour,tanom,latgrid,altitude,/noeras,levels=level(index),color=mcolor,/follow,/overplot
contour,tanom,latgrid,altitude,/noeras,levels=[0.],color=0,thick=3,/follow,/overplot

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.2,charthick=2,xtitle='MLS Temperature Anomaly (K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col2(jj)
x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert MLS_YZ_Anomalies/yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps -rotate -90 MLS_YZ_Anomalies/yz_mls_temp_minus_multiyear_mean_'+sdate+'.png'
   spawn,'rm -f MLS_YZ_Anomalies/yz_mls_temp_minus_multiyear_mean_'+sdate+'.ps'
endif

goto,jump

end
