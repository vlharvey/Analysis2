;
; read CIPS orbits on one day
;
;spawn,'ls /aura7/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_*_2007-148*_v03.20.nc',fnames20_all	; if running on aura
spawn,'ls /Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_*_2007-148*_v03.20.nc',fnames20_all
norbit=n_elements(fnames20_all)
restore,'read_cips_file.sav

ALBLIM=1.	;LOWER LIMIT FOR ALBEDO -- ANYTHING SMALLER THAN THIS IS ASSUMED NOT TO BE A CLOUD.^M
ALBLIM=0.	;LOWER LIMIT FOR ALBEDO -- ANYTHING SMALLER THAN THIS IS ASSUMED NOT TO BE A CLOUD.^M
ERRLIM=1.0	;MAXIMUM ALLOWED RATIO OF ALBEDO_ERR/ALBEDO^M
SZALIM_HI=91.	;DATA WITH SZA > SZALIM_HI ARE SUSPECT^M
SZALIM_LO=50.	;DATA WITH SZA < SZALIM_LO ARE SUSPECT^M

korbit=0L
for iorbit=0L,norbit-1L do begin
fname20=fnames20_all(iorbit)
print,'fname20= ',fname20
;
; strip out orbit number and YYYY-DOY to find corresponding file for v3.22
;
result=strsplit(fname20,'/',/extract)
result2=strsplit(result(n_elements(result)-1),'_',/extract)
sorbit=result2(4)
sdate=result2(5)
dum=findfile('/Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_'+sorbit+'_'+sdate+'_v03.22.nc')
dum2=findfile('/Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_ozo_orbit_'+sorbit+'_'+sdate+'_v03.20.nc')
dum3=findfile('/Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_ozo_orbit_'+sorbit+'_'+sdate+'_v03.22.nc')
if dum(0) eq '' or dum2(0) eq '' or dum3(0) eq '' then goto,skiporbit
fname22='/Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_'+sorbit+'_'+sdate+'_v03.22.nc'
print,sorbit+' '+sdate

data20=read_cips_file(fname20,/full_path,attributes=attributes)
data22=read_cips_file(fname22,/full_path,attributes=attributes)
;
; extract variables from data structures
;
AIM_ORBIT_NUMBER20=data20.AIM_ORBIT_NUMBER                ;INT Cumulative mission orbit number
UTDATE20=DATA20.UT_DATE                                  ;LONG       2009001 UTC date of this orbit
UTTIME20=(*DATA20[0].UT_TIME)                            ;POINTER  Number of seconds elapsed since orbit_start_time_ut
nlayers20=(*DATA20[0].nlayers)
PERCENT_CLOUDS20=DATA20.PERCENT_CLOUDS                    ;FLOAT    Percentage of tiles that have an identified cloud
CLOUD_INDEX20=(*DATA20[0].CLOUD_PRESENCE_MAP)             ;1 FOR CLOUD, 0 FOR NO CLOUD
LATITUDE20=(*DATA20[0].LATITUDE)                          ;latitude for every pixel
LONGITUDE20=(*DATA20[0].LONGITUDE)                        ;longitude for every pixel
X=WHERE(LATITUDE20 GT 90,NX)
IF NX GT 0 THEN LATITUDE20(X)=180-LATITUDE20(X)           ;correct latitude for crossing over the NP
X=WHERE(LATITUDE20 lt -90.,nx)
if nx gt 0L then latitude20(x)=-90.-(latitude20(x)+90.)   ;correct latitude for crossing over the SP
X=WHERE(LONGITUDE20 LT 0,NX)
IF NX GT 0 THEN LONGITUDE20(X)=LONGITUDE20(X)+360
SZA20 = (*DATA20[0].ZENITH_ANGLE_RAY_PEAK)                ;Zenith angle at the peak of the Rayleigh contribution
ALB20 = (*DATA20[0].cld_albedo)                           ;Cloud albedo in Garys (10^-6 sr^-1) i.e., alb x 1.e6
ALB_ERR20 = (*DATA20[0].CLD_ALBEDO_UNC)                   ;1 sigma formal uncertainty of cld_albedo
IWC20 = (*DATA20[0].ICE_WATER_CONTENT)
IWC_ERR20 = (*DATA20[0].ICE_WATER_CONTENT_UNC)
RAD20=(*DATA20[0].PARTICLE_RADIUS)                        ;Particle radius for each tile
RAD_ERR20=(*DATA20[0].PARTICLE_RADIUS_UNC)                ;Particle radius uncertainty for each tile
ozo=DATA20.INDICATORS
RATALL20=(*ozo[0].ratall)

AIM_ORBIT_NUMBER22=data22.AIM_ORBIT_NUMBER                ;INT Cumulative mission orbit number
if AIM_ORBIT_NUMBER22 ne AIM_ORBIT_NUMBER20 then stop,'different orbit numbers,AIM_ORBIT_NUMBER20,AIM_ORBIT_NUMBER22
UTDATE22=DATA22.UT_DATE                                  ;LONG       2009001 UTC date of this orbit
UTTIME22=(*DATA22[0].UT_TIME)                            ;POINTER  Number of seconds elapsed since orbit_start_time_ut
nlayers22=(*DATA22[0].nlayers)
PERCENT_CLOUDS22=DATA22.PERCENT_CLOUDS                    ;FLOAT    Percentage of tiles that have an identified cloud
CLOUD_INDEX22=(*DATA22[0].CLOUD_PRESENCE_MAP)             ;1 FOR CLOUD, 0 FOR NO CLOUD
LATITUDE22=(*DATA22[0].LATITUDE)                          ;latitude for every pixel
LONGITUDE22=(*DATA22[0].LONGITUDE)                        ;longitude for every pixel
X=WHERE(LATITUDE22 GT 90,NX)
IF NX GT 0 THEN LATITUDE22(X)=180-LATITUDE22(X)           ;correct latitude for crossing over the NP
X=WHERE(LATITUDE22 lt -90.,nx)
if nx gt 0L then latitude22(x)=-90.-(latitude22(x)+90.)   ;correct latitude for crossing over the SP
X=WHERE(LONGITUDE22 LT 0,NX)
IF NX GT 0 THEN LONGITUDE22(X)=LONGITUDE22(X)+360
SZA22 = (*DATA22[0].ZENITH_ANGLE_RAY_PEAK)                ;Zenith angle at the peak of the Rayleigh contribution
ALB22 = (*DATA22[0].cld_albedo)                           ;Cloud albedo in Garys (10^-6 sr^-1) i.e., alb x 1.e6
ALB_ERR22 = (*DATA22[0].CLD_ALBEDO_UNC)                   ;1 sigma formal uncertainty of cld_albedo
IWC22 = (*DATA22[0].ICE_WATER_CONTENT)
IWC_ERR22 = (*DATA22[0].ICE_WATER_CONTENT_UNC)
RAD22=(*DATA22[0].PARTICLE_RADIUS)                        ;Particle radius for each tile
RAD_ERR22=(*DATA22[0].PARTICLE_RADIUS_UNC)                ;Particle radius uncertainty for each tile
ozo=DATA22.INDICATORS
RATALL22=(*ozo[0].ratall)
;
; remove bad values
;
GOOD20 = WHERE(FINITE(ALB20) EQ 1 AND FINITE(ALB_ERR20) EQ 1 AND $
             SZA20 LT SZALIM_HI AND SZA20 GT SZALIM_LO AND $
            (CLOUD_INDEX20 EQ 1 OR CLOUD_INDEX20 EQ 0) AND $
             ALB20 GE ALBLIM and nlayers20 gt 6.,NGOOD20)

GOOD22 = WHERE(FINITE(ALB22) EQ 1 AND FINITE(ALB_ERR22) EQ 1 AND $
             SZA22 LT SZALIM_HI AND SZA22 GT SZALIM_LO AND $
            (CLOUD_INDEX22 EQ 1 OR CLOUD_INDEX22 EQ 0) AND $
             ALB22 GE ALBLIM and nlayers22 gt 6.,NGOOD22)

if good20(0) eq -1L or good22(0) eq -1L then goto,skiporbit

UTDATE20=UTDATE20(good20)
UTTIME20=UTTIME20(good20)
PERCENT_CLOUDS20=PERCENT_CLOUDS20(good20)
CLOUD_INDEX20=CLOUD_INDEX20(good20)
LATITUDE20=LATITUDE20(good20)
LONGITUDE20=LONGITUDE20(good20)
SZA20=SZA20(good20)
ALB20=ALB20(good20)
ALB_ERR20=ALB_ERR20(good20)
IWC20=IWC20(good20)
IWC_ERR20=IWC_ERR20(good20)
RAD20=RAD20(good20)
RAD_ERR20=RAD_ERR20(good20)
RATALL20=RATALL20(good20)

UTDATE22=UTDATE22(good22)
UTTIME22=UTTIME22(good22)
PERCENT_CLOUDS22=PERCENT_CLOUDS22(good22)
CLOUD_INDEX22=CLOUD_INDEX22(good22)
LATITUDE22=LATITUDE22(good22)
LONGITUDE22=LONGITUDE22(good22)
SZA22=SZA22(good22)
ALB22=ALB22(good22)
ALB_ERR22=ALB_ERR22(good22)
IWC22=IWC22(good22)
IWC_ERR22=IWC_ERR22(good22)
RAD22=RAD22(good22)
RAD_ERR22=RAD_ERR22(good22)
RATALL22=RATALL22(good22)

sdate=strcompress(min(UTDATE20),/remove_all)
sorbit=strcompress(AIM_ORBIT_NUMBER20,/remove_all)
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,0.2*cos(a),0.2*sin(a),/fill
if korbit eq 0L then begin
   setplot='ps'
   read,'setplot=',setplot
   korbit=1L
endif
nxdim=750
nydim=750
xorig=[0.15,0.55]
yorig=[0.35,0.35]
xlen=0.3
ylen=0.5
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='cips_l4_ratall_v3.20_v3.22_'+sdate+'_'+sorbit+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xyouts,.325,.9,sdate+' Orbit '+sorbit,/normal,color=0,charsize=1.75
plot,sza20,ratall20,psym=8,color=0,title='v3.20',$
     ytitle='RATALL',xtitle='SZA',xrange=[30.,100.],yrange=[0.96,1.01]
index=where(cloud_index20 eq 1 and alb20 ge 1.,ncloud20)	; for comparison to 1G summary file
if ncloud20 gt 0L then oplot,sza20(index),ratall20(index),psym=8,color=.9*mcolor
plots,30,0.992250
plots,100,0.992250,/continue,color=mcolor*.9,thick=3
xyouts,32,0.962,strcompress(ncloud20,/remove_all)+' Clouds',/data,color=0
xyouts,32,0.966,strcompress(ngood20,/remove_all)+' Obs',/data,color=0
xyouts,32,1.008,'Threshold = 0.992250',/data,color=0
if ncloud20 gt 0L then xyouts,32,1.004,'Max RATALL = '+strcompress(max(ratall20(index)),/remove_all),/data,color=0

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,sza22,ratall22,psym=8,color=0,title='v3.22',$
     xtitle='SZA',xrange=[30.,100.],yrange=[0.96,1.01]
index=where(cloud_index22 eq 1 and alb22 ge 1.,ncloud22)
if ncloud22 gt 0L then oplot,sza22(index),ratall22(index),psym=8,color=.9*mcolor
plots,30,0.994369
plots,100,0.994369,/continue,color=mcolor*.9,thick=3
xyouts,32,0.962,strcompress(ncloud22,/remove_all)+' Clouds',/data,color=0
xyouts,32,0.966,strcompress(ngood22,/remove_all)+' Obs',/data,color=0
xyouts,32,1.008,'Threshold = 0.994369',/data,color=0
if ncloud22 gt 0L then xyouts,32,1.004,'Max RATALL = '+strcompress(max(ratall22(index)),/remove_all),/data,color=0

if ncloud20 gt ncloud22 then stop,'v3.20 has more clouds!!!'

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim cips_l4_ratall_v3.20_v3.22_'+sdate+'_'+sorbit+'.ps'+$
         ' -rotate -90 compare_cips_l4_v3.20_v3.22_'+sdate+'_'+sorbit+'.jpg'
   spawn,'rm -f cips_l4_ratall_v3.20_v3.22_'+sdate+'_'+sorbit+'.ps'
endif

HEAP_GC
skiporbit:
endfor	; loop over orbits
end
