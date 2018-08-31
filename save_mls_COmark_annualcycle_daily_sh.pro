;
; save daily zonal mean CO gradient marker for each day (to hopefully match the mls_daily_zonal_means_from_gridded TUV file).
; SH version. Also save daily altitude profile of marker area (% of the hemisphere).
;
loadct,39
mcolor=byte(!p.color)
icmm1=mcolor-1B
icmm2=mcolor-2B
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!NOERAS=-1
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.1
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/Users/harvey/Harvey_etal_2018/Code/Save_files/daily_mls_coelatedge+merra2_sfelatedge_'
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
nmon=n_elements(smonth)

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
index=where(yeq mod 2 eq 0,nrr2)
yeq2=yeq(index)                         ; 2 degree bins

latcircle=fltarr(nrr)
hem_frac=fltarr(nrr)
for j=0,nrr-2 do begin
    hy=re*dtr
    dx=re*cos(yeq(j)*dtr)*360.*dtr
    latcircle(j)=dx*hy
endfor
for j=0,nrr-1 do begin
    if yeq(j) ge 0. then index=where(yeq ge yeq(j))
    if index(0) ne -1 then hem_frac(j)=100.*total(latcircle(index))/hem_area
    if yeq(j) eq 0. then hem_frac(j)=100.
endfor

;
; get lon,lat information
;
restore,'smidemax_300-year_TUmark_djf_jja.sav
;
; MLS daily zonal means
; KDAY            FLOAT     =       4664.00
; NLV             LONG      =           55
; NLV2            LONG      =           37
; NR              LONG      =           96
; LAT             DOUBLE    = Array[96]
; PRESS37         FLOAT     = Array[37]
; PRESS55         FLOAT     = Array[55]
; SDATE_ALL       STRING    = Array[4664]
; TBAR            FLOAT     = Array[96, 55, 4664]
; UBAR            FLOAT     = Array[96, 55, 4664]
; VBAR            FLOAT     = Array[96, 55, 4664]
; ZBAR            FLOAT     = Array[96, 55, 4664]
; H2OBAR          FLOAT     = Array[96, 55, 4664]
; O3BAR           FLOAT     = Array[96, 55, 4664]
; N2OBAR          FLOAT     = Array[96, 37, 4664]
; COBAR           FLOAT     = Array[96, 37, 4664]
;
restore,'/atmos/aura6/data/MLS_data/Pre_process/mls_daily_zonal_means_from_gridded_20040808-20170515.sav

;goto,quick
;
; file listing of MLS CO marker files
;
spawn,'ls '+dir+'*_2d_sh.sav',ifiles
nfile=n_elements(ifiles)
for ifile=0L,nfile-1L do begin
    restore,ifiles(ifile)
    print,ifiles(ifile)
    if ifile eq 0L then sdate_all=SDATE_TIME
    if ifile gt 0L then sdate_all=[sdate_all,SDATE_TIME]
endfor
nday=n_elements(sdate_all)
;
; loop over files
;
icount=0L
kcount=0L
for ifile=0L,nfile-1L do begin
    result=strsplit(ifiles(ifile),'_',/extract)
    yyyymm=result(-3)
    smon=strmid(yyyymm,4,2)
    imon=long(smon)
;
;   if imon ne 1 and imon ne 2 and imon ne 12 and imon ne 6 and imon ne 7 and imon ne 8 then goto,skipmonth		; DJF and JJA only
;   if imon ne 1 and imon ne 7 then goto,skipmonth		; Jan/Jul only
;
; restore monthly data
;
; DELATLEVS3D     FLOAT     = Array[30, 91, 37]
; HLATPDF_TIME_3D FLOAT     = Array[30, 46, 37]
; LLATPDF_TIME_3D FLOAT     = Array[30, 46, 37]
; LOWLAT_ELATEDGE_2D FLOAT     = Array[30, 37]
; LOWLAT_ELATINNER_2D FLOAT     = Array[30, 37]
; LOWLAT_ELATOUTER_2D FLOAT     = Array[30, 37]
; MARKMLS4D       FLOAT     = Array[144, 96, 37, 30]
; MARKSFELATEDGE_2D FLOAT     = Array[30, 37]
; NASHELATEDGE_2D FLOAT     = Array[30, 37]
; NASHINNER_2D    FLOAT     = Array[30, 37]
; NASHOUTER_2D    FLOAT     = Array[30, 37]
; NOVORTEX_FLAG_2D FLOAT     = Array[30, 37]
; PMLS            FLOAT     = Array[37]
; SDATE_TIME      STRING    = Array[30]
; SFELATEDGE_2D   FLOAT     = Array[30, 37]
; SFMARKEDGE_2D   FLOAT     = Array[30, 37]
; SPBIN3D         FLOAT     = Array[30, 91, 37]
; YEQ             FLOAT     = Array[91]
;
        print,ifiles(ifile)
        restore,ifiles(ifile)
; 
; declare arrays
;
        if icount eq 0L then begin
           np=n_elements(pmls)
           zm_markco_mls_daily_sh=fltarr(nr,np,nday)
           area_markco_mls_prof_sh=fltarr(np,nday)

           nr=n_elements(lat)
           nc=144L
           lon=2.5*findgen(nc)
           x2d=0.*fltarr(nc,nr)
           y2d=0.*fltarr(nc,nr)
           for i=0,nc-1 do y2d(i,*)=lat
           for j=0,nr-1 do x2d(*,j)=lon
           area=0.*y2d
           deltax=lon(1)-lon(0)
           deltay=lat(1)-lat(0)
           for j=0,nr-1 do begin
               hy=re*deltay*dtr
               dx=re*cos(lat(j)*dtr)*deltax*dtr
               area(*,j)=dx*hy    ; area of each grid point
           endfor

           icount=1L
        endif
;
; daily zonal means
;
        zmmark=mean(MARKMLS4D,dim=1,/Nan)		; daily zonal mean 
        zm_markco_mls_daily_sh(*,*,kcount:kcount+n_elements(SDATE_TIME)-1L)=zm_markco_mls_daily_sh(*,*,kcount:kcount+n_elements(SDATE_TIME)-1L)+zmmark
;
; save profile of Arctic vortex area - ADD THIS LATER
;
        for ii=0L,n_elements(SDATE_TIME)-1L do begin
            mark3d=reform(MARKMLS4D(*,*,*,ii))
            for kk=0L,np-1L do begin
                mark2d=reform(mark3d(*,*,kk))
                index=where(mark2d gt 0. and y2d lt 0.)
                if index(0) ne -1L then area_markco_mls_prof_sh(kk,kcount+ii)=100.*total(area(index))/hem_area
            endfor	; loop over pressure levels
;print,SDATE_TIME(ii),' ',area_markco_mls_prof_sh(*,ii)
        endfor	; loop over days

        kcount=kcount+n_elements(SDATE_TIME)

skipmonth:
endfor  ; loop over files
;
; need height information
;
restore,'mls_djf_jja.sav
zindex=0.*press55
for k = 0,n_elements(PRESS55)-1 do begin
    index = where(press37 eq press55(k))
    if index(0) ne -1 then zindex(k) = 1.0
endfor
good=where(zindex eq 1.0)
zbar=ZBAR_DJF(*,good)/1000.
zprof=mean(zbar,dim=1)	; altitude profile
;
; save daily zonal means and daily profiles of vortex area (SH)
;
save,file='mls_COmark_daily_zm_sh.sav',nr,np,nday,lat,pmls,zprof,zm_markco_mls_daily_sh,area_markco_mls_prof_sh,sdate_all

end
