;
; SAD timeseries
; read multi-year monthly means of WACCM Surface Area Density of NAT
;
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=fix(byte(!p.color))
if mcolor ne 255 then mcolor=255
icmm1=mcolor-1B
icmm2=mcolor-2B
device,decompose=0
nlvls=19
col1=1+indgen(nlvls)*mcolor/nlvls
!NOERAS=-1
!P.FONT=1
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.4,0.7]
xorig=[0.2,0.4,0.7]
yorig=[0.2,0.25,0.25]
xlen=0.6
ylen=0.6
cbaryoff=0.05
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif

dir='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_WACCM_1_T3I_sad_nat_'
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
for imon=0L,nmon-1L do begin

    restore,file=dir+smon(imon)+'.sav'	;,alon,alat,lev,wsadgrd,nwsadgrd
    print,'read CCMVal2_REF-B1_WACCM_1_T3I_sad_nat_'+smon(imon)+'.sav'
;
; compute zonal mean
;
    nr=n_elements(alat)
    nc=n_elements(alon)
    nl=n_elements(lev)
    wsadzm=fltarr(nr,nl)
    for k=0L,nl-1L do begin
    for j=0L,nr-1L do begin
        index=where(wsadgrd(*,j,k) ne 0.,npts)
        if npts ge 3L then begin
           wsadzm(j,k)=total(wsadgrd(index,j,k))/float(npts)
        endif
    endfor
    endfor
;
; timeseries
;
    if imon eq 0L then wsadtime=fltarr(nmon,nr,nl)
    wsadtime(imon,*,*)=wsadzm

endfor
rlev=29.7280 & rlat=71.0526
print,lev
read,'Enter pressure level ',rlev
index=where(abs(rlev-lev) eq min(abs(rlev-lev)))
ilev=index(0)
slev=strcompress(rlev,/remove_all)
;print,alat
;read,'Enter latitude ',rlat
index=where(abs(rlat-alat) eq min(abs(rlat-alat)))
ilat=index(0)
slat=strcompress(rlat,/remove_all)
wsad=reform(wsadtime(*,ilat,ilev))
;
; postscript
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='step3_plot_timeseries_monthly_waccm_ccmval_sad_nat_'+slat+'_'+slev+'hPa.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot
;
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;   index=where(wsadzm ge 0.001)
;   imin=min(wsadzm(index))
;   imax=max(wsadzm(index))
;   iint=(imax-imin)/float(nlvls)
;   level=imin+iint*findgen(nlvls)

level=[0.0001,0.0005,0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5]/10.
nlvls=n_elements(level)
col1=1+indgen(nlvls)*mcolor/nlvls
wsad=wsad*1.e8
    plot,1+findgen(nmon),wsad,title='WACCM SAD NAT '+slat+' '+slev+' hPa',xrange=[1.,12.],yrange=[0.0001,0.1],/ylog,$
         charsize=1.5,ytitle='10^8 * (m-1)',xtitle='Month',color=0,psym=8,xticks=nmon-1,xtickname=smon
for j=0L,nr-1L do begin
wsad=1.e8*reform(wsadtime(*,j,ilev))
oplot,1+findgen(nmon),wsad,psym=8,color=(float(j)/float(nr))*mcolor
endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim step3_plot_timeseries_monthly_waccm_ccmval_sad_nat_'+slat+'_'+slev+'hPa.ps -rotate -90 '+$
                            'step3_plot_timeseries_monthly_waccm_ccmval_sad_nat_'+slat+'_'+slev+'hPa.jpg'
        spawn,'/usr/bin/rm step3_plot_timeseries_monthly_waccm_ccmval_sad_nat_'+slat+'_'+slev+'hPa.ps'
     endif
end
