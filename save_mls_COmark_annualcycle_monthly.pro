;
; multi-year monthly mean annual cycle of T, U, CO gradient marker
; to save daily annual cycle need to accommodate missing days
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
;
; get lon,lat information
;
restore,'smidemax_300-year_TUmark_djf_jja.sav

;goto,quick
;
; file listing of MLS CO marker files
;
spawn,'ls '+dir+'*.sav',ifiles
nfile=n_elements(ifiles)
;
; loop over files
;
icount=0L
for ifile=0L,nfile-1L do begin
    result=strsplit(ifiles(ifile),'_',/extract)
    yyyymm=result(-3)
    smon=strmid(yyyymm,4,2)
    imon=long(smon)

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
        mavg2=mean(MARKMLS4D,dim=4)	; average over days in the month
; 
; declare DJF and JJA arrays
;
        if icount eq 0L then begin
           np=n_elements(pmls)
           zm_markco_mls=fltarr(nmon,nr,np)
           nzm_markco_mls=fltarr(nmon)
           icount=1L
        endif
;
; monthly mean zonal mean
;
        dum=mean(MARKMLS4D,dim=1,/Nan)		; daily zonal mean 
        zm_markco_mls(imon-1,*,*)=zm_markco_mls(imon-1,*,*)+mean(dum,dim=3,/Nan)		; monthly mean 
        nzm_markco_mls(imon-1)=nzm_markco_mls(imon-1)+1L
;
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
; save 3d DJF and JJA CO-Mark
;
save,file='mls_COmark_annualcycle_monthly.sav',nr,np,alat,pmls,zprof,zm_markco_mls,nzm_markco_mls
quick:
restore,'mls_COmark_annualcycle_monthly.sav'
;
; average multi-year
;
for i=0L,nmon-1L do begin
if nzm_markco_mls(i) gt 0. then zm_markco_mls(i,*,*)=zm_markco_mls(i,*,*)/float(nzm_markco_mls(i))
endfor
;
; select level
;
ralt=50.
print,zprof
read,'Enter desired altitude ',ralt
index=where(abs(ralt-zprof) eq min(abs(ralt-zprof)))
ialt=index(0)
salt=strcompress(long(zprof(ialt)),/r)+'km'
zm_markco_mls_lev=reform(zm_markco_mls(*,*,ialt))
;
; plot
;
if setplot eq 'ps' then begin
   lc=0
   !p.font=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
           /bold,/color,bits_per_pixel=8,/helvetica,filename='annualcycle_COMark_mls.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; DJF
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=10
tlevel=0.1+0.1*findgen(nlvls)
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/float(nlvls-1))*mcolor
col1(-1)=col1(-1)-1
contour,zm_markco_mls_lev,1+findgen(nmon),alat,/noera,/cell_fill,color=0,c_color=col1,levels=tlevel,xrange=[1,12],yrange=[-90,90],$
        ytitle='Latitude',charsize=1.5,charthick=2,xticks=nmon-1,xtickname=smonth
contour,zm_markco_mls_lev,1+findgen(nmon),alat,/noera,/foll,color=0,levels=tlevel,/overplot

xmnb=xmx +cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(0)+0.01,yorig(0)+ylen-0.01
!type=2^2+2^3+2^5
imin=min(tlevel)
imax=max(tlevel)
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],/noeras,color=0,charsize=1,ytitle='MLS CO Gradient Vortex'
xbox=[0,10,10,0,0]
y2=imin
dy=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    ybox=[y2,y2,y2+dy,y2+dy,y2]
    polyfill,xbox,ybox,color=col1(j)
    y2=y2+dy
endfor
;
; Close PostScript file and return control to X-windows
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim annualcycle_COMark_mls.ps -rotate -90 annualcycle_COMark_mls.jpg'
;  spawn,'rm -f annualcycle_COMark_mls.ps'
endif

end
