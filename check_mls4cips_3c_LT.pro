;
; check MLS and CIPS LT sampling for both data versions 4 and 5
; VLH 11/4/2017
;
@stddat
@kgmt
@ckday
@kdate
@mkltime

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,0.5*cos(a),0.5*sin(a),/fill
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15,0.55,0.55]
yorig=[0.55,0.15,0.55,0.15]
xlen=0.225
ylen=0.3
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
cdir='/atmos/harvey/CIPS_data/Datfiles/Level_3c_Summary/cips_3c_'
odir='/Volumes/Data/CIPS_data/Datfiles_MLS_DMP/cips_3c_mls_'
version='v04.20_r05'
;version='v05.10_r01'

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
           /bold,/color,bits_per_pixel=8,/helvetica,filename='check_mls4cips_3c_LT_SZA_'+version+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif

ialb=2
salb=string(format='(i2.2)',ialb)+'G'
spawn,'ls '+cdir+'*'+version+'_'+salb+'_all.sav',cfiles
nfile=n_elements(cfiles)
;
; restore CIPS season
;
; ALB             FLOAT     = Array[589, 70]
; ALB_STD         FLOAT     = Array[589, 70]
; DFS             LONG      = Array[589]
; DOY             INT       = Array[589]
; IWC             FLOAT     = Array[589, 70]
; IWC_STD         FLOAT     = Array[589, 70]
; LATHI           INT       = Array[70]
; LATLO           INT       = Array[70]
; LON             FLOAT     = Array[589, 70]
; LTIME           FLOAT     = Array[589, 70]
; NBIN            INT       =       70
; NREV            LONG      =          589
; NUM_CLD         INT       = Array[589, 70]
; NUM_OBS         INT       = Array[589, 70]
; RAD             FLOAT     = Array[589, 70]
; RAD_STD         FLOAT     = Array[589, 70]
; REV             INT       = Array[589]
; SZA             FLOAT     = Array[589, 70]
; UT              FLOAT     = Array[589, 70]
; YEAR            INT       = Array[589]
;
allyears=7+indgen(10)	; 7-16
nyear=n_elements(allyears)
col1=1+indgen(nyear)*mcolor/float(nyear)

for ifile=0L,nfile-1L do begin
    restore,cfiles(ifile)
    print,'restored '+cfiles(ifile)
    result=strsplit(cfiles(ifile),'_',/extract)
    year=long(result(-5))
    shem=result(-6)
    if shem eq 'south' then latlo=-1.*latlo
    if shem eq 'south' then lathi=-1.*lathi
    if shem eq 'south' then year=long(strmid(result(-5),0,2))
;
; since MLS LT sampling is the same each day, read in any day
;
    restore,mdir+'cat_mls_v4.2_20100101.sav'
    mprof=n_elements(longitude)
    mlev=n_elements(altitude)
    muttime=time
    mlat=latitude
    mlon=longitude
    mltime=0.*time
    mldoy=0.*time
    mkltime,muttime,mlon,mltime,mldoy
    if ifile eq 0 then begin
       erase
       ymn=0.3
       ymx=0.7
       set_viewport,0.2,0.8,ymn,ymx
       !type=2^2+2^3
       plot,mltime,mlat,psym=2,color=0,xrange=[0,24],yrange=[-90,90],ytitle='Geographic Latitude',xtitle='Local Time',charsize=2,charthick=2,title=version
       yinc=(ymx-ymn)/float(nyear)
    endif

print,min(latlo),max(latlo)
help,latlo,ltime
latlo2d=0.*ltime
for i=0,nrev-1L do latlo2d(i,*)=latlo
nbin2=nbin/2
LATLO2D1=reform(LATLO2D(*,0:nbin2-1))
LATLO2D2=reform(LATLO2D(*,nbin2:nbin-1))
ltime1=reform(ltime(*,0:nbin2-1))
ltime2=reform(ltime(*,nbin2:nbin-1))

index=where(year eq allyears)
thiscolor=col1(index)
oplot,ltime2,latlo2d2,psym=1,color=thiscolor,symsize=0.5
oplot,ltime1,latlo2d1,psym=4,color=thiscolor
xyouts,18.5,-75,'MLS',/data,color=0,charsize=2,charthick=2
xyouts,19,45,'CIPS ASC',/data,color=0,charsize=1.5,charthick=2
xyouts,12,20,'CIPS DES',/data,color=0,charsize=1.5,charthick=2
xyouts,1.5,-50,'CIPS ASC',/data,color=0,charsize=1.5,charthick=2
xyouts,12,-30,'CIPS DES',/data,color=0,charsize=1.5,charthick=2
if version eq 'v04.20_r05' and ifile lt nyear-1 then xyouts,0.8+0.02,ymx-0.03-ifile*yinc,strcompress(2000L+year,/r),color=thiscolor,/normal,charsize=2,charthick=2
if version eq 'v05.10_r01' and ifile lt nyear then xyouts,0.8+0.02,ymx-0.03-ifile*yinc,strcompress(2000L+year,/r),color=thiscolor,/normal,charsize=2,charthick=2

endfor	; loop over NH and SH yearly files
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim check_mls4cips_3c_LT_SZA_'+version+'.ps -rotate -90 check_mls4cips_3c_LT_SZA_'+version+'.jpg'
;  spawn,'rm -f check_mls4cips_3c_LT_SZA_'+version+'.ps'
endif

end
