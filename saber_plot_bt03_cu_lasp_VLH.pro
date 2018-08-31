;
; plot SABER wave data
;
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
xorig=[0.1,0.55,0.1,0.55,0.1,0.55]
yorig=[.7,.7,.4,.4,.1,.1]
xlen=0.3
ylen=0.2
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif

year=2009
READ,'ENTER YEAR (yyyy): ', YEAR
yearstring    = STRTRIM( STRING(year,'(I4)') )
dirheader='/atmos/harvey/SABER_data/Datfiles/'
;
; list all SABER files for the given year
; /atmos/harvey/SABER_data/Datfiles/SABER_FLUXES_2012_10_71.SAV		: EPDIV_ALL,EPY_ALL,EPZ_ALL,latitude,longitude,zgrid,yyyyddd
; /atmos/harvey/SABER_data/Datfiles/SABER_GEOP_2012_10_71.SAV		: AV_AGEO,ZM_AGEO,ZONAL_GSPECTRA,latitude,longitude,zgrid,yyyyddd
; /atmos/harvey/SABER_data/Datfiles/SABER_PROFS_2012_10_71.SAV		: UPROFILES,VPROFILES,GPROFILES,TPROFILES,latitude,longitude,zgrid,yyyyddd
; /atmos/harvey/SABER_data/Datfiles/SABER_STAWAVES_2012_10_71.SAV	: STAWAVES_G,STAWAVES_T,latitude,longitude,zgrid,yyyyddd
; /atmos/harvey/SABER_data/Datfiles/SABER_TEMP_2012_10_71.SAV		: AV_ATEMP,ZM_ATEMP,ZONAL_TSPECTRA,latitude,longitude,zgrid,yyyyddd
; /atmos/harvey/SABER_data/Datfiles/SABER_WINDS_2012_10_71.SAV		: UWIND,UPERTWIND,VPERTWIND,latitude,longitude,zgrid,yyyyddd
;
; ARRAYS:
; LATITUDE        FLOAT     = Array[35]
; LONGITUDE       FLOAT     = Array[12]
; YYYYDDD         STRING    = Array[62]
; ZGRID           FLOAT     = Array[120]

; AV_AGEO         FLOAT     = Array[62, 12, 35, 120]
; AV_ATEMP        FLOAT     = Array[62, 12, 35, 120]
; EPDIV_ALL       FLOAT     = Array[35, 62, 120]
; EPY_ALL         FLOAT     = Array[35, 62, 120]
; EPZ_ALL         FLOAT     = Array[35, 62, 120]
; GPROFILES       FLOAT     = Array[62, 12, 35, 120]
; STAWAVES_G      FLOAT     = Array[35, 120, 7, 2, 62]
; STAWAVES_T      FLOAT     = Array[35, 120, 7, 2, 62]
; TPROFILES       FLOAT     = Array[62, 12, 35, 120]
; UPERTWIND       COMPLEX   = Array[35, 62, 6, 120]
; UPROFILES       FLOAT     = Array[62, 12, 35, 120]
; UWIND           FLOAT     = Array[35, 62, 120]
; VPERTWIND       COMPLEX   = Array[35, 62, 6, 120]
; VPROFILES       FLOAT     = Array[62, 12, 35, 120]
; ZM_AGEO         FLOAT     = Array[62, 35, 120]
; ZM_ATEMP        FLOAT     = Array[62, 35, 120]
; ZONAL_GSPECTRA  COMPLEX   = Array[35, 62, 6, 120]
; ZONAL_TSPECTRA  COMPLEX   = Array[35, 62, 6, 120]
;
spawn,'ls '+dirheader+'*_'+yearstring+'_*SAV',ifiles
for i=0L,n_elements(ifiles)-1L do begin
    restore,ifiles(i)
    print,ifiles(i)
endfor
days=float(strmid(yyyyddd,4,3))
;
; PLOT
;
yindex=where(latitude eq 75.)
slat=strcompress(long(latitude(yindex(0))),/remove_all)
conarr=reform(EPDIV_ALL(yindex(0),*,*))
conarr=reform(UWIND(yindex(0),*,*))
nz=n_elements(zgrid)
nfiles=n_elements(yyyyddd)

H = 6.5
altgrid=H*zgrid
NUM_CON  =   11
CONVALS  =   FLTARR(num_con)
baseval  = -100.
coninc   = 20.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls

MISSING = WHERE(CONARR EQ 0.)
IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='SABER_6pan_'+yearstring+'_'+slat+'.ps'
   !p.charsize=1.
   !p.thick=1.5
   !p.charthick=3
   !p.charthick=3
   !y.thick=1.5
   !x.thick=1.5
endif
erase
conarr=reform(ZM_ATEMP(*,yindex(0),*))
NUM_CON  =   21
CONVALS  =   FLTARR(num_con)
baseval  = 100.
coninc   = 10.
level=[160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.]      ;,baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls

MISSING = WHERE(CONARR EQ 0.)
IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.

xyouts,.3,.95,'SABER '+yearstring+' at '+slat+' Latitude',charsize=1.5,charthick=1.5,color=0,/normal
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,3,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='Tbar',XTITLE = 'DOY', YTITLE = 'Altitude (km)'
index=where(level le 260.)
CONTOUR, CONARR, days, altgrid, LEVELS = level, /FOLLOW, /overplot,color=0
omin=min(level)
omax=max(level)
xmnb=xorig(0)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(0),yorig(0)+ylen
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(0)+ylen/3.0,'(K)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0

conarr=reform(UWIND(yindex(0),*,*))
NUM_CON  =   9
baseval  = -80.
coninc   = 20.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls

MISSING = WHERE(CONARR EQ 0.)
IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,3,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='Ubar', XTITLE = 'DOY'
index=where(level gt 0.)
CONTOUR, CONARR, days, altgrid, LEVELS = level(index), /FOLLOW, /overplot,color=0
index=where(level le 0.)
contour,CONARR, days, altgrid, LEVELS = level(index), c_linestyle=3,/overplot,color=mcolor
omin=min(level)
omax=max(level)
xmnb=xorig(1)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(1),yorig(1)+ylen
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(1)+ylen/3.0,'(m/s)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0

conarr=2.*reform(abs(ZONAL_GSPECTRA(yindex(0),*,1,*)))*1000.
NUM_CON  =   11
baseval  = 200.
coninc   = 200.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls

MISSING = WHERE(CONARR EQ 0.)
;IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.
;
; interpolate small gaps in time
;
for k=0,nz-1 do begin
    dlev=reform(conarr(*,k))
    for i=1,nfiles-1 do begin
        if dlev(i) eq 0. or dlev(i-1) eq 0. then begin
           for ii=i+1,nfiles-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    conarr(*,k)=dlev
endfor
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,3,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='Wave 1 GPH',XTITLE = 'DOY',YTITLE = 'Altitude (km)'
index=where(level gt 0.)
CONTOUR, CONARR, days, altgrid, LEVELS = level(index), /FOLLOW, /overplot,color=0
omin=min(level)
omax=max(level)
xmnb=xorig(2)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(2)+ylen/3.0,'(km)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0

conarr=2.*reform(abs(ZONAL_GSPECTRA(yindex(0),*,2,*)))*1000.
NUM_CON  =  11
baseval  = 200.
coninc   = 200.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls

MISSING = WHERE(CONARR EQ 0.)
;IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.
;
; interpolate small gaps in time
;
for k=0,nz-1 do begin
    dlev=reform(conarr(*,k))
    for i=1,nfiles-1 do begin
        if dlev(i) eq 0. or dlev(i-1) eq 0. then begin
           for ii=i+1,nfiles-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump2
               endif
           endfor
jump2:
        endif
    endfor
    conarr(*,k)=dlev
endfor
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,3,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='Wave 2 GPH', XTITLE = 'DOY'
index=where(level gt 0.)
CONTOUR, CONARR, days, altgrid, LEVELS = level(index), /FOLLOW, /overplot,color=0
omin=min(level)
omax=max(level)
xmnb=xorig(3)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(3)+ylen/3.0,'(km)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0


conarr=reform(EPY_ALL(yindex(0),*,*))
NUM_CON  =   11
index=where(finite(conarr) eq 1)
baseval  = -20.
coninc   = 4.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
MISSING = WHERE(CONARR EQ 0.)
;IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.
;
; interpolate small gaps in time
;
for k=0,nz-1 do begin
    dlev=reform(conarr(*,k))
    for i=1,nfiles-1 do begin
        if dlev(i) eq 0. or dlev(i-1) eq 0. then begin
           for ii=i+1,nfiles-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump3
               endif
           endfor
jump3:
        endif
    endfor
    conarr(*,k)=dlev
endfor
xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,5,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='Meridional EPF', XTITLE = 'DOY', YTITLE = 'Altitude (km)'
index=where(level gt 0.)
CONTOUR, CONARR, days, altgrid, LEVELS = level(index), /FOLLOW, /overplot,color=0
index=where(level lt 0.)
contour,CONARR, days, altgrid, LEVELS = level(index), c_linestyle=3,/overplot,color=mcolor
omin=min(level)
omax=max(level)
xmnb=xorig(4)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(4)+ylen/3.0,'(kg/ms!u2!n)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0

conarr=reform(EPDIV_ALL(yindex(0),*,*))
NUM_CON  =   11
index=where(finite(conarr) eq 1)
baseval  = -200.
coninc   = 40.
level=baseval+coninc*findgen(num_con)
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
MISSING = WHERE(CONARR EQ 0.)
;IF MISSING(0) NE -1 THEN CONARR(MISSING) = 0./0.
;
; interpolate small gaps in time
;
for k=0,nz-1 do begin
    dlev=reform(conarr(*,k))
    for i=1,nfiles-1 do begin
        if dlev(i) eq 0. or dlev(i-1) eq 0. then begin
           for ii=i+1,nfiles-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump4
               endif
           endfor
jump4:
        endif
    endfor
    conarr(*,k)=dlev
endfor
xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
conarr=smooth(conarr,5,/nan,/edge_truncate)
CONTOUR, CONARR, days, altgrid, LEVELS = level,/cell_fill,c_color=col1,/noeras,yrange = [10,130], Xrange = [5,75], color=0, $
         title='EPF DIV', XTITLE = 'DOY'
index=where(level lt 0.)
CONTOUR, CONARR, days, altgrid, LEVELS = level(index), /FOLLOW, /overplot,color=0,thick=3
index=where(level eq 0.)
contour,CONARR, days, altgrid, LEVELS = level(index), c_linestyle=3,/overplot,color=mcolor
omin=min(level)
omax=max(level)
xmnb=xorig(5)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=mcolor
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
xyouts,xmxb+0.035,yorig(5)+ylen/3.0,'(m/s/day)',/normal,charsize=1.25,orientation=90.,color=0
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim SABER_6pan_'+yearstring+'_'+slat+'.ps -rotate -90 SABER_6pan_'+yearstring+'_'+slat+'.jpg'
endif
end
