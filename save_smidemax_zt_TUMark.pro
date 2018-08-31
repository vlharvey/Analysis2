;
; save average annual cycle of T, U, TEM, GW, etc...and CO Mark at all latitudes and altitudes
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
xorig=[0.15,0.15]
yorig=[0.6,0.15]
xlen=0.7
ylen=0.35
cbaryoff=0.05
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO2x1SmidEmax_yBWCN/3d_CO2x1SmidEmax_yBWCN_'
dir1='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO_Vortex_data/daily_waccm-smidemax_coelatedge+sfelatedge_'
dir2='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO2x1SmidEmax_yBWCN/ZM_CO2x1SmidEmax_yBWCN_'

smonth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
lstdoy=0*mno
leddoy=0*mno
lstdoy(0)=0L
leddoy(0)=mno(0)-1
for i=1L,n_elements(mno)-1 do lstdoy(i)=total(mno(0:i-1))		; DOY of the first day of each month (IDL index)
for i=1L,n_elements(mno)-1 do leddoy(i)=total(mno(0:i))-1		; DOY of the last day of each month (IDL index)
nmonth=n_elements(smonth)

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
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
;goto,quick
;goto,quick1
;
; build MMDD dates and read in multi-year daily means to build average annual cycle in isobaric and isentropic data
;
year1files=file_search(dir+'001????.nc3')
ndays=n_elements(year1files)
mmdd=strarr(ndays)
for ii=0L,ndays-1L do begin
    dum=strsplit(year1files(ii),'.',/extract)
    dum2=strsplit(dum(0),'_',/extract)
    mmdd(ii)=strmid(dum2(-1),3,4)
;
; read multi-year daily mean zonal means: omega, qjoule, qrs_aur, qrs_euv, ttgw, t, utgw, v*, w*, z
;
    restore,dir2+mmdd(ii)+'.sav'      ; zm and sig: omega, qjoule, qrs_aur, qrs_euv, ttgw, t, utgw, v*, w*, z
    restore,dir+mmdd(ii)+'.sav'      ; zm and sig: co, ipv, mark, p, qdf, sf, U, V, z 
    print,'reading '+dir2+mmdd(ii)+'.sav'
stop

    if ii eq 0L then begin
       nr=n_elements(lat)
       nz=n_elements(lev)
       EKGWSPEC_DAILY_MEAN_3d=fltarr(nr,nz+1,ndays)	; save mean annual cycle at all latitudes and altitudes
       OMEGA_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       QJOULE_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       QRS_AUR_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       QRS_EUV_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       TTGW_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       T_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       UTGWSPEC_DAILY_MEAN_3d=fltarr(nr,nz,ndays)
       VSTAR_DAILY_MEAN_3d=fltarr(nr,nz+1,ndays)
       WSTAR_DAILY_MEAN_3d=fltarr(nr,nz+1,ndays)
       Z3_DAILY_MEAN_3d=fltarr(nr,nz,ndays)

       nth=n_elements(th)
       co_on_th_3d=fltarr(nr,nth,ndays)
       pv_on_th_3d=fltarr(nr,nth,ndays)
       mark_on_th_3d=fltarr(nr,nth,ndays)
markarea_on_th_3d_nh=fltarr(nth,ndays)
markarea_on_th_3d_sh=fltarr(nth,ndays)
       p_on_th_3d=fltarr(nr,nth,ndays)
       u_on_th_3d=fltarr(nr,nth,ndays)
       v_on_th_3d=fltarr(nr,nth,ndays)
       z_on_th_3d=fltarr(nr,nth,ndays)
    endif
    EKGWSPEC_DAILY_MEAN_3d(*,*,ii)=EKGWSPEC_DAILY_MEAN
    OMEGA_DAILY_MEAN_3d(*,*,ii)=OMEGA_DAILY_MEAN
    QJOULE_DAILY_MEAN_3d(*,*,ii)=QJOULE_DAILY_MEAN
    QRS_AUR_DAILY_MEAN_3d(*,*,ii)=QRS_AUR_DAILY_MEAN
    QRS_EUV_DAILY_MEAN_3d(*,*,ii)=QRS_EUV_DAILY_MEAN
    TTGW_DAILY_MEAN_3d(*,*,ii)=TTGW_DAILY_MEAN
    T_DAILY_MEAN_3d(*,*,ii)=T_DAILY_MEAN
    UTGWSPEC_DAILY_MEAN_3d(*,*,ii)=UTGWSPEC_DAILY_MEAN
    VSTAR_DAILY_MEAN_3d(*,*,ii)=VSTAR_DAILY_MEAN
    WSTAR_DAILY_MEAN_3d(*,*,ii)=WSTAR_DAILY_MEAN
    Z3_DAILY_MEAN_3d(*,*,ii)=Z3_DAILY_MEAN

    co_on_th_3d(*,*,ii)=mean(COAVG,dim=2)
    pv_on_th_3d(*,*,ii)=mean(IPVAVG,dim=2)
    mark_on_th_3d(*,*,ii)=mean(MAVG,dim=2)
    p_on_th_3d(*,*,ii)=mean(PAVG,dim=2)
    u_on_th_3d(*,*,ii)=mean(UAVG,dim=2)
    v_on_th_3d(*,*,ii)=mean(VAVG,dim=2)
    z_on_th_3d(*,*,ii)=mean(ZAVG,dim=2)

;markarea_on_th_3d_nh=fltarr(*,ii)=
;markarea_on_th_3d_sh=fltarr(*,ii)=

endfor
save,file='3d_zm_annual_cycle_smidemax.sav',nc,nr,nz,nth,alon,alat,lev,ilev,th,ndays,mmdd,EKGWSPEC_DAILY_MEAN_3d,OMEGA_DAILY_MEAN_3d,$
     QJOULE_DAILY_MEAN_3d,QRS_AUR_DAILY_MEAN_3d,QRS_EUV_DAILY_MEAN_3d,TTGW_DAILY_MEAN_3d,T_DAILY_MEAN_3d,UTGWSPEC_DAILY_MEAN_3d,$
     VSTAR_DAILY_MEAN_3d,WSTAR_DAILY_MEAN_3d,Z3_DAILY_MEAN_3d,co_on_th_3d,pv_on_th_3d,mark_on_th_3d,p_on_th_3d,u_on_th_3d,v_on_th_3d,z_on_th_3d
quick1:
restore,'3d_zm_annual_cycle_smidemax.sav'
;
; loop over years and months and retain all zonal mean CO gradient vortex positions
;
comark_on_th_4d=fltarr(nr,nth,ndays,300)
comarkarea_on_th_3d_nh=fltarr(nth,ndays,300)
comarkarea_on_th_3d_sh=fltarr(nth,ndays,300)
;
; area
;
lon=0.*fltarr(nc,nr)
lat=0.*fltarr(nc,nr)
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
x2d=fltarr(nc,nr)
y2d=fltarr(nc,nr)
for i=0,nc-1 do y2d(i,*)=alat
for j=0,nr-1 do x2d(*,j)=alon
;
; loop over years
;
for iyear=1,2 do begin				; change to 300 for final analysis
    syear=string(FORMAT='(I3.3)',iyear)
    for imonth=0L,nmonth-1L do begin
        imon=imonth+1
        smon=string(FORMAT='(I2.2)',imon)
;
; restore monthly data
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
        spawn,'ls '+dir1+syear+smon+'_2d_?h.sav',ifiles
        if ifiles(0) eq '' then goto,skipmonth

        if n_elements(ifiles) eq 1L then begin					; one file per month
           print,ifiles(0)
           restore,ifiles(0)
           cozm=mean(MARKCO4D,dim=1)						; Array[nr,nth,nday]
           comark_on_th_4d(*,*,lstdoy(imonth):leddoy(imonth),iyear-1L)= cozm

result=strsplit(ifiles(0),'_',/extract)
hemlab=strmid(result(-1),0,2)
for iday=0L,n_elements(SDATE_TIME)-1L do begin
    area_prof=fltarr(nth)
    for kk=0L,nth-1L do begin
        marklev=reform(MARKCO4D(*,*,kk,iday))
        if hemlab eq 'nh' then index=where(marklev gt 0. and y2d gt 0.)
        if hemlab eq 'sh' then index=where(marklev gt 0. and y2d lt 0.)
        if index(0) ne -1L then area_prof(kk)=total(area(index))/hem_area
    endfor
    if hemlab eq 'nh' then comarkarea_on_th_3d_nh(*,lstdoy(imonth)+iday,iyear-1)=area_prof
    if hemlab eq 'sh' then comarkarea_on_th_3d_sh(*,lstdoy(imonth)+iday,iyear-1)=area_prof
endfor
        endif
        if n_elements(ifiles) eq 2L then begin                                  ; two files per month in Oct and Nov
           for ii=0L,n_elements(ifiles)-1L do begin
               print,ifiles(ii)
               restore,ifiles(ii)
               cozm=mean(MARKCO4D,dim=1)                                            ; Array[nr,nth,nday]
               comark_on_th_4d(*,*,lstdoy(imonth):leddoy(imonth),iyear-1L)=comark_on_th_4d(*,*,lstdoy(imonth):leddoy(imonth),iyear-1L) + cozm

result=strsplit(ifiles(ii),'_',/extract)
hemlab=strmid(result(-1),0,2)
for iday=0L,n_elements(SDATE_TIME)-1L do begin
    area_prof=fltarr(nth)
    for kk=0L,nth-1L do begin
        marklev=reform(MARKCO4D(*,*,kk,iday))
        if hemlab eq 'nh' then index=where(marklev gt 0. and y2d gt 0.)
        if hemlab eq 'sh' then index=where(marklev gt 0. and y2d lt 0.)
        if index(0) ne -1L then area_prof(kk)=total(area(index))/hem_area
    endfor
    if hemlab eq 'nh' then comarkarea_on_th_3d_nh(*,lstdoy(imonth)+iday,iyear-1)=area_prof
    if hemlab eq 'sh' then comarkarea_on_th_3d_sh(*,lstdoy(imonth)+iday,iyear-1)=area_prof
endfor
           endfor
        endif
    endfor		; loop over months
skipmonth:
endfor		; loop over years
;
; average annual cycle
; when all 300 years are read change MAX to MEAN
;
comark_on_th_3d=max(comark_on_th_4d,dim=4)			; comark_on_th_4d=fltarr(nr,nth,ndays,300)
comarkarea_on_th_2d_nh=transpose(max(comarkarea_on_th_3d_nh,dim=3))	; comarkarea_on_th_3d_nh=fltarr(nth,ndays,300)
comarkarea_on_th_2d_sh=transpose(max(comarkarea_on_th_3d_sh,dim=3))	; comarkarea_on_th_3d_sh=fltarr(nth,ndays,300)
;
; save file
;
save,file='3d_comark_smidemax.sav',nc,nr,nth,alon,alat,th,mmdd,ndays,comark_on_th_3d,comarkarea_on_th_2d_nh,comarkarea_on_th_2d_sh

quick:
restore,'3d_comark_smidemax.sav
;
; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_NH_SH_vortex_smidemax_2pan.ps'
   !p.charsize=1
   !p.thick=2
   !p.charthick=2
   !y.thick=2
   !x.thick=2
endif
;
; choose latitude - from there decide whether you want NH or SH files
;
;print,alat
rlat=60.
nhindex=where(alat ge rlat)
shindex=where(alat le -1.*rlat)
slat=strcompress(long(rlat),/r)
;
zt_mark_nh=transpose(mean(mark_on_th_3d(nhindex,*,*),dim=1))		; mark_on_th_3d=fltarr(nr,nth,ndays)
zt_z_nh=transpose(mean(z_on_th_3d(nhindex,*,*),dim=1))
zt_mark_sh=transpose(mean(mark_on_th_3d(shindex,*,*),dim=1))               ; mark_on_th_3d=fltarr(nr,nth,ndays)
zt_z_sh=transpose(mean(z_on_th_3d(shindex,*,*),dim=1))
;
;comarkarea_on_th_2d_nh
;comarkarea_on_th_2d_sh

; shift NH to put winter in the middle
;
;znhshift=0.*znh
;mnhshift=0.*mnh
;znhshift(0:183,*)=reform(znh(181:364,*))        ; July-Dec
;znhshift(184:364,*)=reform(znh(0:180,*))        ; Jan-July
;mnhshift(0:183,*)=reform(mnh(181:364,*))        ; July-Dec
;mnhshift(184:364,*)=reform(mnh(0:180,*))        ; Jan-July
;
; plot
;
    smon=strmid(mmdd,0,2)
    sday=strmid(mmdd,2,2)
    xindex=where(sday eq '15',nxticks)
    xlabs=smon(xindex)
    xlabshift=xlabs
;   xlabshift=[xlabs(6:-1),xlabs(0:5)]
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
tlevel=[0.01,0.05,0.1+0.1*findgen(10)]
nlvls=n_elements(tlevel)
col1=(findgen(nlvls)/float(nlvls-1))*mcolor
col1(-1)=col1(-1)-1
    contour,comarkarea_on_th_2d_nh,findgen(ndays),zt_z_nh,/noera,/fill,color=0,c_color=col1,levels=tlevel,xrange=[0,ndays-1],yrange=[30,125],ytitle='Altitude (km)',charsize=1,charthick=2,title='NH',$
            xticks=nxticks-1,xtickname=xlabshift,xtickv=xindex
    contour,comarkarea_on_th_2d_nh,findgen(ndays),zt_z_nh,/noera,/follow,color=0,levels=tlevel(0:-1:2),/overplot
    contour,zt_mark_nh,findgen(ndays),zt_z_nh,/noera,/foll,color=250,thick=5,levels=[0.1,0.5,0.9],/overplot

    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,comarkarea_on_th_2d_sh,findgen(ndays),zt_z_sh,/noera,/fill,color=0,c_color=col1,levels=tlevel,xrange=[0,ndays-1],yrange=[30,125],ytitle='Altitude (km)',charsize=1,charthick=2,title='SH',$
            xticks=nxticks-1,xtickname=xlabshift,xtickv=xindex
    contour,comarkarea_on_th_2d_sh,findgen(ndays),zt_z_sh,/noera,/follow,color=0,levels=tlevel(0:-1:2),/overplot
    contour,zt_mark_sh,findgen(ndays),zt_z_sh,/noera,/foll,color=250,thick=5,levels=[0.1,0.5,0.9],/overplot
;
imin=min(tlevel)
imax=max(tlevel)
ymnb=min(yorig) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Polar Vortex Area (Frac Hem)',/noeras,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;
; Close PostScript file and return control to X-windows
;
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim zt_NH_SH_vortex_smidemax_2pan.ps -rotate -90 zt_NH_SH_vortex_smidemax_2pan.jpg
;      spawn,'rm -f zt_NH_SH_vortex_smidemax_2pan.ps'
    endif

end
