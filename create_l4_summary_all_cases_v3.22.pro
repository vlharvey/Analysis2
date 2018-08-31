; =======================================================================================
; Generates the Level 4 summary data binned in latitude.  Binning is done orbit-by-orbit
; using a defined albedo threshold.  
; This version of the code generates all three versions of the summary data, corresponding 
; to binning for:
;                 all pixels
;                 cloud pixels only
;                 non-cloud pixels only
;
; Results are saved to separate IDL save files for all three cases.
; =======================================================================================

; ----------------------------------------------------
; Define the season and hemisphere to analyze
; ----------------------------------------------------

yr='2008'
yr='2007'

hem='south'
hem='north'

season=hem+'_'+yr

; -----------------------------------------
; Set albedo threshold to use for binning.
; -----------------------------------------

threshold='10G'
alb_thresh=10.
threshold='01G'
alb_thresh=1.
threshold='05G'
alb_thresh=5.

; ------------------------------------------------------------------
; Read in list of Level 4 file names for this season (pre-compiled)
; ------------------------------------------------------------------

datpath='/aim/data/cips/v3.20/'+season+'/level_4/'
filein='./file_lists/'+season+'.files'

files=read_ascii_file(filein)
rev=fix(strmid(files,17,5))
year=fix(strmid(files,23,4))
doy=fix(strmid(files,28,3))

; ------------------------------------------------------------------------
; Set min/max revs for each season and cut input files to this range.
; ------------------------------------------------------------------------

if(season eq 'north_2007') then revminmax=[435,2149]
if(season eq 'south_2007') then revminmax=[2537,4280]
if(season eq 'north_2008') then revminmax=[5579,7590]
if(season eq 'south_2008') then revminmax=[8188,10300]

x=where(rev ge revminmax[0] and rev le revminmax[1])
files=files(x)
rev=rev(x)
year=year(x)
doy=doy(x)
nrev=n_elements(rev)

; -----------------------------
; Calculate days from solstice 
; -----------------------------

if(hem eq 'north') then begin
  sol=julday(6,21,year(0))-julday(1,1,year(0))
  dfs=doy-sol
endif else begin
  sol=julday(12,21,year(0))-julday(1,1,year(0))
  dfs=doy-sol
  newyear=where(doy lt 180)
  dfs(newyear)=dfs(newyear)+365
endelse

; ----------------------
; Define latitude bins..
; ----------------------

latlo=[50,55,60,65,70,75,80,90,80,75,70,65,60,55]
lathi=[55,60,65,70,75,80,90,80,75,70,65,60,55,50]
nbin=n_elements(latlo)

; ---------------------------------------------------------------
; Define which bins correspond to ascending and descending nodes
; (this is strictly by convention)
; ---------------------------------------------------------------

astart=0      ; Ascending node bins
astop=6       ;          "
dstart=7      ; Descending node bins
dstop=nbin-1  ;          "

; -------------------------------
; Define and initialize arrays 
; -------------------------------

def=-999.    ; default bad value
nmin=1       ; Minimum # of good points needed to calculate and average in a lat bin.
sza_min=50.  ; Minimum solar zenith angle allowed in CIPS data.
sza_max=91.  ; Maximum solar zenith angle allowed in CIPS data.


; ---------------------------------------------
; Define arrays for all three binning scenarios
; ---------------------------------------------

  freq=fltarr(nrev,nbin)+def

  lon_all=freq
  sza_all=freq
ltime_all=freq
   ut_all=freq
  alb_all=freq
  rad_all=freq
  iwc_all=freq
  alb_std_all=freq
  rad_std_all=freq
  iwc_std_all=freq
  num_cld_all=intarr(nrev,nbin)
  num_obs_all=num_cld_all

  lon_cld=freq
  sza_cld=freq
ltime_cld=freq
   ut_cld=freq
  alb_cld=freq
  rad_cld=freq
  iwc_cld=freq
  alb_std_cld=freq
  rad_std_cld=freq
  iwc_std_cld=freq
  num_cld_cld=num_cld_all
  num_obs_cld=num_cld_all

  lon_nocld=freq
  sza_nocld=freq
ltime_nocld=freq
   ut_nocld=freq
  alb_nocld=fltarr(nrev,nbin)
  rad_nocld=freq
  iwc_nocld=fltarr(nrev,nbin)
  alb_std_nocld=freq
  rad_std_nocld=freq
  iwc_std_nocld=freq
  num_cld_nocld=num_cld_all
  num_obs_nocld=num_cld_all

; ---------------------------------------------------------------
; Loop over orbits, read  Level 4 data and calculate averages.
; ---------------------------------------------------------------

for irev=0,nrev-1 do begin      ; Loop over orbits

  rev0=rev(irev)
  file=files(irev)

  dat=read_cips_file(datpath+file,/full_path)

  nx=dat.xdim
  ny=dat.ydim
  ut_in=*dat.ut_time
  lat_in=*dat.latitude
  lon_in=*dat.longitude
  sza_in=*dat.zenith_angle_ray_peak
  alb_in=*dat.cld_albedo
  rad_in=*dat.particle_radius
  iwc_in=*dat.ice_water_content
       c=*dat.ozone_col_density
   sigma=*dat.scale_height_ratio
  cld_presence=*dat.cloud_presence_map

  asc=where(abs(lat_in) gt 90.,complement=dsc)

; ---------------------------
; Reorder the lat/lon arrays
; ---------------------------

  g=where(lon_in lt 0,ng)
  if(ng gt 0) then lon_in(g)+=360.

  if(dat.hemisphere eq 'N' and asc[0] ne -1) then lat_in(asc)=180.-lat_in(asc)
  if(dat.hemisphere eq 'S' and asc[0] ne -1) then lat_in(asc)=-180.-lat_in(asc)

; ------------------------
; Calculate local time
; ------------------------

  ltime_in=fltarr(nx,ny)*!values.f_nan
;  nadir_rayleigh=fltarr(nx,ny)*!values.f_nan
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
       if(finite(ut_in[ix,iy]) eq 1 and finite(lon_in[ix,iy]) eq 1) then begin
         ltime_in(ix,iy)=ut_in(ix,iy)-(360.-lon_in(ix,iy))*24./360.
         if(ltime_in(ix,iy) lt 0.) then ltime_in(ix,iy)=ltime_in(ix,iy)+24.
       endif
; testing......
;      if(finite(c[ix,iy]) eq 1) then begin
;        no3=c[ix,iy]*1.e15
;        sig=sigma[ix,iy]
;        zenith=sza_in[ix,iy]
;        rp=rayleigh_phase([180.-zenith])
;        nadir_rayleigh[ix,iy]=calc_albedo_smb(no3,sig,rp,[0.],[zenith])/1.e-6
;      endif
;.................
  endfor
  endfor


; ------------------------
; Loop over latitude bins
; ------------------------

  for k=0,nbin-1 do begin

    if(k le astop) then begin   ; we're in ascending node
      latmin=latlo(k)
      latmax=lathi(k)
      acdc=asc
    endif else begin            ; we're in descending node (note - order is switched)
      latmin=lathi(k)       
      latmax=latlo(k)
      acdc=dsc
    endelse
    if(acdc[0] eq -1) then goto,skip

;-----------------------------------------------------------------------
;    Find points in this bin. Discriminate cloud and non-cloud pixels,
;    and impose albedo and solar zenith angle screens. 
;-----------------------------------------------------------------------

    all=where(finite(sza_in(acdc)) eq 1 and  $
              abs(lat_in(acdc)) ge latmin and abs(lat_in(acdc)) lt latmax ,nall)

    cld=where(finite(sza_in(acdc)) eq 1 and  $
              sza_in(acdc) ge sza_min and sza_in(acdc) le sza_max and $
              abs(lat_in(acdc)) ge latmin and abs(lat_in(acdc)) lt latmax and $
              cld_presence(acdc) eq 1 and alb_in(acdc) ge alb_thresh,ncld)

    nocld=where(finite(sza_in(acdc)) eq 1 and $
              sza_in(acdc) ge sza_min and sza_in(acdc) le sza_max and $
              abs(lat_in(acdc)) ge latmin and abs(lat_in(acdc)) lt latmax and $
              cld_presence(acdc) eq 0 ,nnocld)

    low=where(finite(sza_in(acdc)) eq 1 and $
              sza_in(acdc) ge sza_min and sza_in(acdc) le sza_max and $
              abs(lat_in(acdc)) ge latmin and abs(lat_in(acdc)) lt latmax and $
              cld_presence(acdc) eq 1 and alb_in(acdc) lt alb_thresh,nlow)

    total=where(finite(sza_in(acdc)) eq 1 and $
              sza_in(acdc) ge sza_min and sza_in(acdc) le sza_max and $
              abs(lat_in(acdc)) ge latmin and abs(lat_in(acdc)) lt latmax and $
              cld_presence(acdc) ge 0,ntot)

;    not calculating frequency anymore
;    if(ntot ne 0) then freq[irev,k]=100.*float(ncld)/(float(ntot))

;   -------------------
;   Cloud points option
;   -------------------

    if(ncld ge nmin) then begin

      num_cld_cld[irev,k]=ncld
      num_obs_cld[irev,k]=ncld

      sza_cld[irev,k]=median(sza_in(acdc(cld)))
      ltime_cld[irev,k]=median(ltime_in(acdc(cld)))
      ut_cld[irev,k]=median(ut_in(acdc(cld)))

      dm=lon_in(acdc(cld))
      df=dm-mean(dm)
      if(max(dm)-min(dm) gt 300) then begin  ;crossed 0 lon line
        p=where(dm gt 180)
        dm(p)-=360.
        lonavg=median(dm)
        if(lonavg lt 0) then lonavg+=360.
      endif else begin
        lonavg=median(dm)
      endelse
      lon_cld[irev,k]=lonavg

      alb_cld[irev,k]=median(alb_in(acdc(cld)))
      rad_cld[irev,k]=median(rad_in(acdc(cld)))
      iwc_cld[irev,k]=median(iwc_in(acdc(cld)))

      if(ncld ge 5) then begin
        alb_std_cld[irev,k]=stddev(alb_in(acdc(cld)))
        rad_std_cld[irev,k]=stddev(rad_in(acdc(cld)))
        iwc_std_cld[irev,k]=stddev(iwc_in(acdc(cld)))
      endif


    endif

;   -----------------------
;   No cloud points option
;   -----------------------

    if(nnocld ge nmin) then begin

      num_obs_nocld[irev,k]=nnocld

      sza_nocld[irev,k]=median(sza_in(acdc(nocld)))
      ltime_nocld[irev,k]=median(ltime_in(acdc(nocld)))
      ut_nocld[irev,k]=median(ut_in(acdc(nocld)))

      dm=lon_in(acdc(nocld))
      df=dm-mean(dm)
      if(max(dm)-min(dm) gt 300) then begin  ;crossed 0 lon line
        p=where(dm gt 180)
        dm(p)-=360.
        lonavg=median(dm)
        if(lonavg lt 0) then lonavg+=360.
      endif else begin
        lonavg=median(dm)
      endelse
      lon_nocld[irev,k]=lonavg

;     albedo, radius and IWC fields have no dynamic values for this case

    endif

;   -------------------
;   All points option
;   -------------------

    if(ntot ge nmin) then begin

      num_cld_all[irev,k]=ncld
      num_obs_all[irev,k]=ntot

      sza_all[irev,k]=median(sza_in(acdc(total)))
      ltime_all[irev,k]=median(ltime_in(acdc(total)))
      ut_all[irev,k]=median(ut_in(acdc(total)))

      dm=lon_in(acdc(total))
      df=dm-mean(dm)
      if(max(dm)-min(dm) gt 300) then begin  ;crossed 0 lon line
        p=where(dm gt 180)
        dm(p)-=360.
        lonavg=median(dm)
        if(lonavg lt 0) then lonavg+=360.
      endif else begin
        lonavg=median(dm)
      endelse
      lon_all[irev,k]=lonavg

      if(ncld eq 0) then begin
        alb_all[irev,k]=0.0
        iwc_all[irev,k]=0.0
      endif else begin
        dum_alb=alb_in[acdc[cld]]
        dum_iwc=iwc_in[acdc[cld]]
        if(nnocld+nlow gt 0) then begin
           dum_alb=[dum_alb,fltarr(nnocld+nlow)]
           dum_iwc=[dum_iwc,fltarr(nnocld+nlow)]
        endif
        alb_all[irev,k]=mean(dum_alb)  
        iwc_all[irev,k]=mean(dum_iwc)  
      endelse

    endif

  endfor

  heap_gc
  skip:
endfor   ; end loop over orbits...

; ----------------------
; Save data to files
; ----------------------

savfile='cips_l4_summary_'+threshold+'_'+season+'_v03.20.sav.all'
lon=lon_all
sza=sza_all
alb=alb_all
rad=rad_all
iwc=iwc_all
alb_std=alb_std_all
rad_std=rad_std_all
iwc_std=iwc_std_all
ltime=ltime_all
ut=ut_all
num_cld=num_cld_all
num_obs=num_obs_all
save,file=savfile, nrev,nbin,rev,year,doy,dfs,latlo,lathi, $
                   alb,rad,iwc,alb_std,rad_std,iwc_std,ltime,ut,lon,sza,num_cld,num_obs

savfile='cips_l4_summary_'+threshold+'_'+season+'_v03.20.sav.cld'
lon=lon_cld
sza=sza_cld
alb=alb_cld
rad=rad_cld
iwc=iwc_cld
alb_std=alb_std_cld
rad_std=rad_std_cld
iwc_std=iwc_std_cld
ltime=ltime_cld
ut=ut_cld
num_cld=num_cld_cld
num_obs=num_obs_cld
save,file=savfile, nrev,nbin,rev,year,doy,dfs,latlo,lathi, $
                   alb,rad,iwc,alb_std,rad_std,iwc_std,ltime,ut,lon,sza,num_cld,num_obs

savfile='cips_l4_summary_'+season+'_v03.20.sav.nocld'
lon=lon_nocld
sza=sza_nocld
alb=alb_nocld
rad=rad_nocld
iwc=iwc_nocld
alb_std=alb_std_nocld
rad_std=rad_std_nocld
iwc_std=iwc_std_nocld
ltime=ltime_nocld
ut=ut_nocld
num_cld=num_cld_nocld
num_obs=num_obs_nocld
save,file=savfile, nrev,nbin,rev,year,doy,dfs,latlo,lathi, $
                   alb,rad,iwc,alb_std,rad_std,iwc_std,ltime,ut,lon,sza,num_cld,num_obs

end
