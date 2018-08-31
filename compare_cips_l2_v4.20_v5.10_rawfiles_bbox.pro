;
; compare CIPS level 2 v4.2 and v5.10
; plot BBOX values for recent seasons
; VLH 1/10/2017
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
xorig=[0.15,0.15,0.15]
yorig=[0.7,0.4,0.1]
xlen=0.7
ylen=0.225
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
pth='/atmos/harvey/CIPS_data/Datfiles/Level_2/cips_sci_2_orbit_'
;
; loop over years
;
for iyear=2015,2015 do begin		; 2015/2016 is the last SH v4.2 season
syear=strcompress(long(iyear),/r)

;goto,quick

lstmn=11
lstdy=11	; dfs -40
lstyr=iyear
ledmn=3
leddy=10		; dfs +80
ledyr=iyear+1
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date ',lstmn,lstdy,lstyr
;read,' Enter ending date ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
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
kcount=0L
;
; USE THE CLOUD PRESENCE MAP ARRAY TO CALCULATE FREQUENCIES. AND USES THE CPM=1 VALUE TO GET THE ALBEDOS.
; LOWER LIMIT FOR ALBEDO -- ANYTHING SMALLER THAN THIS IS ASSUMED NOT TO BE A CLOUD.
; USING LIM=-99 ESSENTIALLY INCLUDES ALL POINTS THAT ARE FOUND WITH CLOUD_PRESENCE_MAP,
; EVEN IF THE ALBEDO IS NEGATIVE (WHICH DOES HAPPEN) -- BUT THEN THE ALB/ALB_ERR TEST MIGHT CATCH IT.
ALBLIM=2.
;ERRLIM=1.0      ;MAXIMUM ALLOWED RATIO OF ALBEDO_ERR/ALBEDO - not used here
SZALIM_HI=92.      ;DATA WITH SZA > SZALIM ARE SUSPECT
SZALIM_LO=42.

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotyear
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sday=string(FORMAT='(I3.3)',iday)
      sdate=syr+smn+sdy
;
; convert YYYYMMDD to DFS
;
;     jday=julday(smn,sdy,syr)
;     yyyy=strmid( strcompress(sdate_all(0),/remove_all), 0,4)		; year on DFS -40
;     dfs_all(icount)=jday-julday(12,21,yyyy)
;     print,sdate,' ',dfs_all(icount)
;
; get nc filenames on this day (I get an error if I try to read the gzipped file)
;
      spawn,'ls '+pth+'*'+syr+'-'+sday+'_v04.20_r05_cld.nc',fnames20
      spawn,'ls '+pth+'*'+syr+'-'+sday+'_v04.20_r05_cat.nc',fnamescat20
      spawn,'ls '+pth+'*'+syr+'-'+sday+'_v05.10_r01_cld.nc',fnames51
      spawn,'ls '+pth+'*'+syr+'-'+sday+'_v05.10_r01_cat.nc',fnamescat51
      if fnames20(0) eq '' or fnames51(0) eq '' then begin
         print,'no orbit files for one of the versions'
         goto,skipcips
      endif
      norbit20=n_elements(fnames20)
      norbit51=n_elements(fnames51)
      if norbit20 ne norbit51 then begin
         print,'different number of orbit files'
         goto,skipcips
      endif
      norbit=norbit20
;
; loop over orbits
;
      norbit=2L       ; testing purposes
      FOR iorbit = 0,norbit-1 DO BEGIN
          FNAME=FNAMESCAT20(iorbit)
          print,fname
;
; read 4.2 catalog file
;
          ncid=ncdf_open(fname)
          result=ncdf_inquire(ncid)
          for idim=0,result.ndims-1 do begin
              ncdf_diminq,ncid,idim,name,dim
              if name eq 'dim1_NLAYERS' then dim1_NLAYERS=dim
              if name eq 'dim2_NLAYERS' then dim2_NLAYERS=dim
              if name eq 'dim1_RATALL' then dim1_RATALL=dim
              if name eq 'dim2_RATALL' then dim2_RATALL=dim
              if name eq 'dim1_BBOX' then dim1_BBOX=dim
              if name eq 'dim1_LATITUDE' then dim1_LATITUDE=dim
              if name eq 'dim2_LATITUDE' then dim2_LATITUDE=dim
              if name eq 'dim1_LONGITUDE' then dim1_LONGITUDE=dim
              if name eq 'dim2_LONGITUDE' then dim2_LONGITUDE=dim
              if name eq 'dim1_ZENITH_ANGLE_RAY_PEAK' then dim1_ZENITH_ANGLE_RAY_PEAK=dim
              if name eq 'dim2_ZENITH_ANGLE_RAY_PEAK' then dim2_ZENITH_ANGLE_RAY_PEAK=dim
              if name eq 'dim1_COMMON_VOLUME_MAP' then dim1_COMMON_VOLUME_MAP=dim
              if name eq 'dim2_COMMON_VOLUME_MAP' then dim2_COMMON_VOLUME_MAP=dim
;             print,'read ',name,' dimension ',dim
          endfor
          for ivar=0,result.nvars-1 do begin
              result=ncdf_varinq(ncid,ivar)
              ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
              if result.name eq 'AIM_ORBIT_NUMBER' then AIM_ORBIT_NUMBER=data           ; Integer orbit number
              if result.name eq 'VERSION' then VERSION=data                             ; Data version number
              if result.name eq 'REVISION' then REVISION=data                           ; Data revision number
              if result.name eq 'PRODUCT_CREATION_TIME' then PRODUCT_CREATION_TIME=data ; String containing UT time at which data file was produced
              if result.name eq 'DEPENDENT_1B_VERSION' then DEPENDENT1BVERSION=data	; Version of lower level 1B data used to produce this data set
              if result.name eq 'UT_DATE' then UTDATE=data				; UT time for each element (fractional hour)
              if result.name eq 'UT_TIME' then UTTIME=data				; UT time for each element (fractional hour)
              if result.name eq 'HEMISPHERE' then HEMISPHERE=data                       ; N (north) or S (south)
              if result.name eq 'ORBIT_START_TIME' then ORBIT_START_TIME=data           ; GPS start time of orbit (microseconds from 0000 UT on 6 Jan 1980)
              if result.name eq 'ORBIT_START_TIME_UT' then ORBIT_START_TIME_UT=data     ; Start time of orbit in yyyymmdd-hr:min:sec format
              if result.name eq 'ORBIT_END_TIME' then ORBIT_END_TIME=data               ; GPS end time of orbit (microseconds from 0000 UT on 6 Jan 1980)
              if result.name eq 'STACK_ID' then STACK_ID=data                           ; Obsolete
              if result.name eq 'XDIM' then XDIM=data                                   ; number of along track elements
              if result.name eq 'YDIM' then YDIM=data                                   ; number of cross track elements
              if result.name eq 'NLAYERS' then NLAYERS=data                             ; Number of observations at the location of each element; each observation corresponds to a different observing geometry and thus SCA in the phase function 
              if result.name eq 'RATALL' then RATALL=data                               ; Indicator of forward vs. backward scattering ratio [see Bailey et al., 2009] / [1933,412], range: 0 to 1.28827
              if result.name eq 'QUALITY_FLAGS' then QUALITY_FLAGS=data                 ; Indicators of data quality for each element.  In v4.20 for NLayers > 5, QF=0.  NLayers = 4 or 5, QF=1. NLayers < 4, QF=2.
              if result.name eq 'KM_PER_PIXEL' then KM_PER_PIXEL=data                   ; Linear dimension of square pixel occupying area of CIPS resolution element.
              if result.name eq 'BBOX' then BBOX=data                                   ; Bounding Box: Bottom-Left and Top-Right indices of the smallest rectangle which both circumscribes a set of cells on a grid and is parallel to the grid axes
              if result.name eq 'CENTER_LON' then CENTER_LON=data                       ; Center longitude of the orbit 
              if result.name eq 'LATITUDE' then LATITUDE=data                           ; Latitude of each element; Latitudes greater (less) than 90 (-90) indicate ascending node data
              if result.name eq 'LONGITUDE' then LONGITUDE=data                         ; Longitude of each element. Ranges from -180 to 180
              if result.name eq 'ZENITH_ANGLE_RAY_PEAK' then SZA=data			; Solar zenith angle (SZA) of each element. 
              if result.name eq 'COMMON_VOLUME_MAP' then COMMON_VOLUME_MAP=data         ; Indicator for whether this location is within the single "Common Volume" where both CIPS and SOFIE observe each orbit.  1 = in ; 0 = not in
              if result.name eq 'NOTES' then NOTES=data                                 ; NOTES
;             print,'read variable ',result.name
          endfor
          ncdf_close,ncid
          sorbit=strcompress(AIM_ORBIT_NUMBER,/remove_all)
;         latitude_orig=latitude                                ;save original lats to determine asc/desc
;         longitude_orig=longitude                              ;save original lons to determine asc/desc
;         X=WHERE(LATITUDE GT 90,NX)
;         IF NX GT 0 THEN LATITUDE(X)=180-LATITUDE(X)           ;correct latitude for crossing over the NP
;         X=WHERE(LATITUDE lt -90.,nx)
;         if nx gt 0L then latitude(x)=-90.-(latitude(x)+90.)   ;correct latitude for crossing over the SP
;         X=WHERE(LONGITUDE LT 0,NX)
;         IF NX GT 0 THEN LONGITUDE(X)=LONGITUDE(X)+360
;
;---DO NOT--- GET RID OF ANY INFINITE DATA YET
;
;         good=WHERE(finite(uttime) eq 1,ngood)	; all data
;         IF NGOOD GT 0 THEN BEGIN
;    
; retain all 4.2 orbit data
;
             if icount eq 0L then begin
;               LAT_ALL20=latitude(good)
;               LON_ALL20=longitude(good)
;               LATORIG_ALL20=latitude_orig(good)
;               LONORIG_ALL20=longitude_orig(good)

                orbit_all20=AIM_ORBIT_NUMBER
                date_all20=long(string(utdate))
                bbox_all20=lonarr(kday*30,4)			; build a 2D array (so first element matches date to be plotted) - truncate used x-elements at the end
                clon_all20=CENTER_LON
             endif
             if icount gt 0L then begin
;               LAT_ALL20=[LAT_ALL20,latitude(good)]
;               LON_ALL20=[LON_ALL20,longitude(good)]
;               LATORIG_ALL20=[LATORIG_ALL20,latitude_orig(good)]
;               LONORIG_ALL20=[LONORIG_ALL20,longitude_orig(good)]

                orbit_all20=[orbit_all20,AIM_ORBIT_NUMBER]
                date_all20=[date_all20,long(string(utdate))]
                clon_all20=[clon_all20,CENTER_LON]			; center_lon is 1 number/orbit. need array of same for plotting seasonal results
             endif

             bbox_all20(icount,*)=bbox                       

;         ENDIF
;
; read v5.10 catalogue file
;
          FNAME=FNAMESCAT51(iorbit)
          print,fname
          ncid=ncdf_open(fname)
          result=ncdf_inquire(ncid)
          for idim=0,result.ndims-1 do begin
              ncdf_diminq,ncid,idim,name,dim
              if name eq 'dim1_UT_TIME' then dim1_UT_TIME=dim
              if name eq 'dim2_UT_TIME' then dim2_UT_TIME=dim
              if name eq 'dim1_NLAYERS' then dim1_NLAYERS=dim
              if name eq 'dim2_NLAYERS' then dim2_NLAYERS=dim
              if name eq 'dim1_LATITUDE' then dim1_LATITUDE=dim
              if name eq 'dim2_LATITUDE' then dim2_LATITUDE=dim
              if name eq 'dim1_LONGITUDE' then dim1_LONGITUDE=dim
              if name eq 'dim2_LONGITUDE' then dim2_LONGITUDE=dim
              if name eq 'dim1_ZENITH_ANGLE_RAY_PEAK' then dim1_ZENITH_ANGLE_RAY_PEAK=dim
              if name eq 'dim2_ZENITH_ANGLE_RAY_PEAK' then dim2_ZENITH_ANGLE_RAY_PEAK=dim
              if name eq 'dim1_COMMON_VOLUME_MAP' then dim1_COMMON_VOLUME_MAP=dim
              if name eq 'dim2_COMMON_VOLUME_MAP' then dim2_COMMON_VOLUME_MAP=dim
;             print,'read ',name,' dimension ',dim
          endfor
          for ivar=0,result.nvars-1 do begin
              result=ncdf_varinq(ncid,ivar)
              ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
              if result.name eq 'AIM_ORBIT_NUMBER' then AIM_ORBIT_NUMBER=data           ; short AIM_ORBIT_NUMBER(structure_elements) ;
              if result.name eq 'VERSION' then VERSION=data                             ; char VERSION(structure_elements, string) ;
              if result.name eq 'REVISION' then REVISION=data                           ; char REVISION(structure_elements, string) ;
              if result.name eq 'PRODUCT_CREATION_TIME' then PRODUCT_CREATION_TIME=data ; char PRODUCT_CREATION_TIME(structure_elements, string) ;
              if result.name eq 'DEPENDENT_1B_VERSION' then DEPENDENT1BVERSION=data	; char DEPENDENT_1B_VERSION(structure_elements, string) ;
              if result.name eq 'UT_DATE' then UTDATE=data				; int UT_DATE(structure_elements) ;
              if result.name eq 'UT_TIME' then UTTIME=data				; double UT_TIME(structure_elements) ;
              if result.name eq 'HEMISPHERE' then HEMISPHERE=data                       ; char HEMISPHERE(structure_elements, string) ;
              if result.name eq 'ORBIT_START_TIME' then ORBIT_START_TIME=data           ; double ORBIT_START_TIME(structure_elements) ;
              if result.name eq 'ORBIT_START_TIME_UT' then ORBIT_START_TIME_UT=data     ; double ORBIT_START_TIME_UT(structure_elements) ;
              if result.name eq 'ORBIT_END_TIME' then ORBIT_END_TIME=data               ; double ORBIT_END_TIME(structure_elements) ;
;             if result.name eq 'STACK_ID' then STACK_ID=data                           ; short STACK_ID(structure_elements) ;
              if result.name eq 'XDIM' then XDIM=data                                   ; int XDIM(structure_elements) ;
              if result.name eq 'YDIM' then YDIM=data                                   ; int YDIM(structure_elements) ;
              if result.name eq 'NLAYERS' then NLAYERS=data                             ; short NLAYERS(structure_elements, dim2_NLAYERS, dim1_NLAYERS) ;
              if result.name eq 'RATALL' then RATALL=data                               ; float RATALL(structure_elements, dim2_RATALL, dim1_RATALL) ;
              if result.name eq 'QUALITY_FLAGS' then QUALITY_FLAGS=data                 ; int QUALITY_FLAGS(structure_elements) ;
              if result.name eq 'KM_PER_PIXEL' then KM_PER_PIXEL=data                   ; float KM_PER_PIXEL(structure_elements) ;
              if result.name eq 'BBOX' then BBOX=data                                   ; int BBOX(structure_elements, dim1_BBOX) ;
              if result.name eq 'CENTER_LON' then CENTER_LON=data                       ; double CENTER_LON(structure_elements) ;
              if result.name eq 'LATITUDE' then LATITUDE=data                           ; float LATITUDE(structure_elements, dim2_LATITUDE, dim1_LATITUDE) ;
              if result.name eq 'LONGITUDE' then LONGITUDE=data                         ; float LONGITUDE(structure_elements, dim2_LONGITUDE, dim1_LONGITUDE) ;
              if result.name eq 'ZENITH_ANGLE_RAY_PEAK' then SZA=data			; float ZENITH_ANGLE_RAY_PEAK(structure_elements, dim2_ZENITH_ANGLE_RAY_PEAK, dim1_ZENITH_ANGLE_RAY_PEAK) ;
              if result.name eq 'COMMON_VOLUME_MAP' then COMMON_VOLUME_MAP=data         ; byte COMMON_VOLUME_MAP(structure_elements, dim2_COMMON_VOLUME_MAP, dim1_COMMON_VOLUME_MAP) ;
              if result.name eq 'NOTES' then NOTES=data                                 ; char NOTES(structure_elements, string) ;
;             print,'read variable ',result.name
          endfor
          ncdf_close,ncid
;         sorbit=strcompress(AIM_ORBIT_NUMBER,/remove_all)
;         latitude_orig=latitude                                ;save original lats to determine asc/desc
;         longitude_orig=longitude                              ;save original lons to determine asc/desc
;         X=WHERE(LATITUDE GT 90,NX)
;         IF NX GT 0 THEN LATITUDE(X)=180-LATITUDE(X)           ;correct latitude for crossing over the NP
;         X=WHERE(LATITUDE lt -90.,nx)
;         if nx gt 0L then latitude(x)=-90.-(latitude(x)+90.)   ;correct latitude for crossing over the SP
;         X=WHERE(LONGITUDE LT 0,NX)
;         IF NX GT 0 THEN LONGITUDE(X)=LONGITUDE(X)+360
;
;---DO NOT--- GET RID OF ANY INFINITE DATA YET
;
;         good=WHERE(finite(uttime) eq 1,ngood)
;         IF NGOOD GT 0 THEN BEGIN
;
; retain 5.1 orbit data
;
             if icount eq 0L then begin
;               LAT_ALL51=latitude(good)
;               LON_ALL51=longitude(good)
;               LATORIG_ALL51=latitude_orig(good)
;               LONORIG_ALL51=longitude_orig(good)

                orbit_all51=AIM_ORBIT_NUMBER
                date_all51=long(string(utdate))
                bbox_all51=lonarr(kday*30,4)                    ; build a 2D array (so first element matches date to be plotted) - truncate used x-elements at the end
                clon_all51=CENTER_LON
             endif
             if icount gt 0L then begin
;               LAT_ALL51=[LAT_ALL51,latitude(good)]
;               LON_ALL51=[LON_ALL51,longitude(good)]
;               LATORIG_ALL51=[LATORIG_ALL51,latitude_orig(good)]
;               LONORIG_ALL51=[LONORIG_ALL51,longitude_orig(good)]

                orbit_all51=[orbit_all51,AIM_ORBIT_NUMBER]
                date_all51=[date_all51,long(string(utdate))]
                clon_all51=[clon_all51,CENTER_LON]                      ; center_lon is 1 number/orbit. need array of same for plotting seasonal results
             endif
;         ENDIF

             bbox_all51(icount,*)=bbox

          icount=icount+1L
      endfor  ; loop over orbits

      skipcips:
;     icount=icount+1L
goto,jump

plotyear:
;
; truncate first element of bbox arrays to match length of clon,date,orbit arrays
;
npts=n_elements(clon_all51)
bbox_all20=reform(bbox_all20(0:npts-1,*))	; assumes that there are same number of orbits in 4.2 as 5.1
bbox_all51=reform(bbox_all51(0:npts-1,*))
stop

quick:

;
; postscript file
;
;     if setplot eq 'ps' then begin
;        lc=0
;        set_plot,'ps'
;        xsize=nxdim/100.
;        ysize=nydim/100.
;        !p.font=0
;        device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
;               /bold,/color,bits_per_pixel=8,/helvetica,filename='compare_cips_l4_v4.20_v5.10_rawfiles_bbox_'+sdate+'.ps'
;        !p.charsize=1.25
;        !p.thick=2
;        !p.charthick=5
;        !p.charthick=5
;        !y.thick=2
;        !x.thick=2
;     endif

;     erase
;     xmn=xorig(0)
;     xmx=xorig(0)+xlen
;     ymn=yorig(0)
;     ymx=yorig(0)+ylen
;     set_viewport,xmn,xmx,ymn,ymx
;     !type=2^2+2^3
;     plot,sza20,alb20,psym=1,color=0,ytitle='ALB',xtitle='SZA',xrange=[SZALIM_LO,SZALIM_HI],yrange=[0.,100.],title=sdate+' (v4.2 black, v5.1 color)',charsize=1.5,charthick=2

;   if setplot ne 'ps' then stop
;     if setplot eq 'ps' then begin
;        device, /close
;        spawn,'convert -trim compare_cips_l4_v4.20_v5.10_rawfiles_bbox_'+sdate+'.ps -rotate -90 compare_cips_l4_v4.20_v5.10_rawfiles_bbox_'+sdate+'.jpg'
;        spawn,'rm -f compare_cips_l4_v4.20_v5.10_rawfiles_bbox_'+sdate+'.ps'
;     endif


endfor	; loop over years
end
