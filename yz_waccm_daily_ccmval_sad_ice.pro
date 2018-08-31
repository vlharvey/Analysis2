;
; zonal mean daily SAD
; restore save files from Doug Kinnison of WACCM Surface Area Density of ICE
; /aura7/harvey/CCMval_data/Datfiles
;
; CCMVal2_REF-B1_WACCM_1_T3I_sad_ice_1960-1989.nc
; CCMVal2_REF-B1_WACCM_1_T3I_sad_ice_1990-1999.nc
; CCMVal2_REF-B1_WACCM_1_T3I_sad_ice_2000-2005.nc
;
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

loadct,39
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=1
!p.charsize=1
!p.charthick=2
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.25]
yorig=[0.15]
xlen=0.5
ylen=0.7
cbaryoff=0.05
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif

dir='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_WACCM_1_T3I_sad_ice_'
spawn,'ls '+dir+'*sav',ncfiles
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
    ncfile=ncfiles(ifile)
    print,'opening '+ncfile
    restore,ncfile	; alon,alat,lev,sadgrd
;
; set date
;
result=strsplit(ncfile,'_',/extract)
result2=strsplit(result(8),'.',/extract)
sdate=result2(0)
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='yz_waccm_daily_ccmval_sad_ice_'+sdate+'.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; compute zonal mean SAD
;
nc=n_elements(alon)
nr=n_elements(alat)
nl=n_elements(lev)
    sadzm=-99.+0.*fltarr(nr,nl)
    for j=0L,nr-1L do begin
    for k=0L,nl-1L do begin
        index=where(sadgrd(*,j,k) gt 1.e-10,nx)
        if nx ge 3L then sadzm(j,k)=mean(sadgrd(index,j,k))*1.e6
    endfor
    endfor
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
    index=where(sadzm ne -99.)
    imin=min(sadzm(index))
    imax=max(sadzm(index))
    int=(imax-imin)/float(nlvls)
    level=imin+int*findgen(nlvls)
    contour,sadzm,alat,lev,levels=level,/ylog,/cell_fill,c_color=col1,color=0,title='WACCM SAD ICE '+sdate,$
            xtitle='Latitude',xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[200.,10.],$
            ytitle='Pressure (hPa)',min_value=-99.
    contour,sadzm,alat,lev,levels=level,/follow,color=0,/overplot

    imin=min(level)
    imax=max(level)
    ymnb=yorig(0)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
         xtitle='m(-1)',charsize=2,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
    endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim yz_waccm_daily_ccmval_sad_ice_'+sdate+'.ps -rotate -90 '+$
                            'yz_waccm_daily_ccmval_sad_ice_'+sdate+'.jpg'
        spawn,'/usr/bin/rm yz_waccm_daily_ccmval_sad_ice_'+sdate+'.ps'
     endif
endfor	; loop over files
end
