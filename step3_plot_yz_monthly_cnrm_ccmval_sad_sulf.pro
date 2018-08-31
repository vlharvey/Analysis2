;
; read multi-year monthly means of Surface Area Density of sulf and plot zonal mean
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

dir='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_sulf_'
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
for imon=0L,nmon-1L do begin

    restore,file=dir+smon(imon)+'.sav'	;,alon,alat,lev,wsadgrd,nwsadgrd
    print,'read CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_sulf_'+smon(imon)+'.sav'
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
; postscript
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='step3_plot_yz_monthly_cnrm_ccmval_sad_sulf_'+smon(imon)+'.ps'
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
print,smon(imon),' ',max(wsadzm)

wsadzm=wsadzm*1.e8
    contour,wsadzm,alat,lev,/ylog,c_color=col1,/cell_fill,title='CNRM SAD Sulfur '+smon(imon),xrange=[-90.,90.],yrange=[500.,1.],$
         charsize=1.5,ytitle='Pressure (hPa)',xtitle='Latitude',color=0,levels=level
    contour,wsadzm,alat,lev,/ylog,color=0,levels=level,/overplot,/follow,c_charsize=1.5,c_labels=0*level

    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    slab=' '+strarr(n_elements(level))
    plot,[min(level),max(level)],[0,0],yrange=[0,10],color=0,xrange=[min(level),max(level)],$
          xticks=n_elements(level)-1L,xtickname=slab,/noeras,charsize=1.5,xtitle='10^8 * (m!u-1!n)'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

slab=strcompress(string(format='(f6.4)',level),/remove_all)
x2=xmn
dx=(xmx-xmn)/(float(nlvls)-1)
for i=0L,n_elements(slab)-1L do begin
    slab0=slab(i)
    flab0=float(slab(i))
    if flab0 ge 1. then slab0=strcompress(long(slab0),/remove_all)
    if flab0 lt 1. then slab0=strcompress(string(format='(f3.1)',flab0),/remove_all)
    if flab0 lt 0.1 then slab0=strcompress(string(format='(f4.2)',flab0),/remove_all)
    if flab0 lt 0.01 then slab0=strcompress(string(format='(f5.3)',flab0),/remove_all)
    if flab0 lt 0.001 then slab0=strcompress(string(format='(f6.4)',flab0),/remove_all)
    if flab0 lt 0.0001 then slab0=strcompress(string(format='(f7.5)',flab0),/remove_all)
    xyouts,x2,ymnb-0.05,slab0,charsize=1.1,/normal,orientation=90.,color=0,charthick=2
    x2=x2+dx
endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim step3_plot_yz_monthly_cnrm_ccmval_sad_sulf_'+smon(imon)+'.ps -rotate -90 '+$
                            'step3_plot_yz_monthly_cnrm_ccmval_sad_sulf_'+smon(imon)+'.jpg'
        spawn,'/usr/bin/rm step3_plot_yz_monthly_cnrm_ccmval_sad_sulf_'+smon(imon)+'.ps'
     endif

endfor
end
