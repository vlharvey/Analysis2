;
; multi-year DJF and JJA YZ of T, U, CO gradient marker
; SmidEmax 300 years
; monthly mean of daily averages data 
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
xorig=[0.15,0.55,0.15,0.55]
yorig=[0.55,0.55,0.15,0.15]
xlen=0.325
ylen=0.325
cbaryoff=0.1
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO_Vortex_data/daily_waccm-smidemax_coelatedge+sfelatedge_'
smonth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
nmonth=n_elements(smonth)
;
; get lon,lat information
;
restore,'smidemax_300-year_TUmark_djf_jja.sav

;goto,quick
;
; loop over years and months
;
icount=0L
for iyear=1,300 do begin
    syear=string(FORMAT='(I3.3)',iyear)
    for imonth=0L,nmonth-1L do begin
        imon=imonth+1
        smon=string(FORMAT='(I2.2)',imon)
        if imon ne 1 and imon ne 2 and imon ne 12 and imon ne 6 and imon ne 7 and imon ne 8 then goto,skipmonth		; DJF and JJA only
;       if imon ne 1 and imon ne 7 then goto,skipmonth		; Jan/Jul only
;
; restore monthly data
;
        ofile=dir+syear+smon+'_2d_?h.sav'
        dum=findfile(ofile)
        print,'reading '+dum
;
; DELATLEVS3D     FLOAT     = Array[31, 91, 22]
; HLATPDF_TIME_3D FLOAT     = Array[31, 46, 22]
; LLATPDF_TIME_3D FLOAT     = Array[31, 46, 22]
; LOWLAT_ELATEDGE_2D FLOAT     = Array[31, 22]
; LOWLAT_ELATINNER_2D FLOAT     = Array[31, 22]
; LOWLAT_ELATOUTER_2D FLOAT     = Array[31, 22]
; MARKCO4D        FLOAT     = Array[144, 96, 22, 31]
; NASHELATEDGE_2D FLOAT     = Array[31, 22]
; NASHINNER_2D    FLOAT     = Array[31, 22]
; NASHOUTER_2D    FLOAT     = Array[31, 22]
; NOVORTEX_FLAG_2D FLOAT     = Array[31, 22]
; SDATE_TIME      STRING    = Array[31]
; SFELATEDGE_2D   FLOAT     = Array[31, 22]
; SPBIN3D         FLOAT     = Array[31, 91, 22]
; TH              FLOAT     = Array[22]
; YEQ             FLOAT     = Array[91]
;
        restore,dum
        mavg2=mean(MARKCO4D,dim=4)	; average over days in the month
; 
; declare DJF and JJA arrays
;
        if icount eq 0L then begin
           djf_markco=0.*mavg2
           jja_markco=0.*mavg2
           ndjf_markco=0.*mavg2
           njja_markco=0.*mavg2
           icount=1L
        endif
;
; DJF
;
    if smon eq '12' or smon eq '01' or smon eq '02' then begin
       djf_markco=djf_markco+mavg2
       ndjf_markco=ndjf_markco+1.
    endif
;
; JJA
;
    if smon eq '06' or smon eq '07' or smon eq '08' then begin
       jja_markco=jja_markco+mavg2
       njja_markco=njja_markco+1.
    endif

skipmonth:
endfor  ; loop over months
endfor  ; loop over years

djf_markco=djf_markco/ndjf_markco
jja_markco=jja_markco/njja_markco
;
; save 3d DJF and JJA CO-Mark
;
save,file='smidemax_300-year_COmark_djf_jja.sav',nc,nr,nth,alon,alat,th,djf_markco,jja_markco
quick:
restore,'smidemax_300-year_COmark_djf_jja.sav'
;
; calculate zonal means
;
djf_markcoyz=mean(djf_markco,dim=1)
jja_markcoyz=mean(jja_markco,dim=1)
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
           /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_multi-year_djf_jja_COMark_smidemax.ps'
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
col1=(findgen(nlvls)/float(nlvls))*mcolor
myz=jja_markcoyz
zyz=mean(JJA_Z,dim=2)
contour,myz,alat,zyz,/noera,/fill,color=0,c_color=col1,levels=tlevel,xrange=[-90,90],yrange=[30,100],ytitle='Altitude (km)',charsize=1.5,charthick=2,title='JJA',xticks=6
contour,myz,alat,zyz,/noera,/foll,color=0,levels=tlevel,/overplot
myz=mean(JJA_MARK,dim=2)
contour,myz,alat,zyz,/noera,/foll,color=mcolor*.9,levels=[0.1,0.5,0.9],/overplot,thick=3

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
myz=djf_markcoyz
zyz=mean(DJF_Z,dim=2)
contour,myz,alat,zyz,/noera,/fill,color=0,c_color=col1,levels=tlevel,xrange=[-90,90],yrange=[30,100],charsize=1.5,charthick=2,title='DJF',xticks=6
contour,myz,alat,zyz,/noera,/foll,color=0,levels=tlevel,/overplot
myz=mean(DJF_MARK,dim=2)
contour,myz,alat,zyz,/noera,/foll,color=mcolor*.9,levels=[0.1,0.5,0.9],/overplot,thick=3

xmnb=xmx +cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,yorig(1)+0.01,yorig(1)+ylen-0.01
!type=2^2+2^3+2^5
imin=min(tlevel)
imax=max(tlevel)
plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],/noeras,color=0,charsize=1,ytitle='Vortex Marker'
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
   spawn,'convert -trim yz_multi-year_djf_jja_COMark_smidemax.ps -rotate -90 yz_multi-year_djf_jja_COMark_smidemax.jpg'
;  spawn,'rm -f yz_multi-year_djf_jja_COMark_smidemax.ps'
endif

end
