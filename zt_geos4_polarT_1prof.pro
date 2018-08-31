;
; compute daily Arctic GEOS-4 Z, T, U, V, Q, PV for time period at one 
; random gridpoint poleward of 70 N
; plot altitude-time section 
; save entire 3-D (y,z,t) arrays for time period
;
@kdate
@rd_ukmo_nc3
@rd_geos5_dat
@rd_geos5_nc3_meto

loadct,38
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
xorig=[0.15]
yorig=[0.225]
xlen=0.7
ylen=0.5
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smon=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
dir='/aura7/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
;dir='/Users/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
gdir='/aura7/harvey/GEOS4_data/Datfiles/'
;gdir='/Users/harvey/GEOS4_data/Datfiles/'
goto,plotit
spawn,'ls '+gdir+'DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.origp.200601*_1200.V01.sav',ncfiles0
spawn,'ls '+gdir+'DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.origp.200602*_1200.V01.sav',ncfiles1
ncfiles=[ncfiles0,ncfiles1]
icount=0L
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
;
; restore GEOS data
;
    restore,file=ncfiles(ifile)        ; alonold,alatold,press,znew,unew,vnew,tnew,pvnew,qnew
    glon=alonold
    glat=alatold
    p=press
    nlg=n_elements(glon)
    nlat=n_elements(glat)
    nlv=n_elements(p)
;
; single profile location
;
ilat=5  ; 77.5 N
slat=strcompress(string(format='(f5.1)',glat(ilat))+'N',/remove_all)
;ilon=5  ; 18.75 E
ilon=25  ; 93.75 E
;ilon=55  ; 206.25 E
;ilon=75  ; 281.25 E
slon=strcompress(string(format='(f6.2)',glon(ilon))+'E',/remove_all)
;
; daily profiles of Z, T, U, V, Q, PV
;
    zbar=reform(znew(ilon,ilat,*))
    tbar=reform(tnew(ilon,ilat,*))
    ubar=reform(unew(ilon,ilat,*))
    vbar=reform(vnew(ilon,ilat,*))
    qbar=reform(qnew(ilon,ilat,*))
    pvbar=reform(pvnew(ilon,ilat,*))

    result=strsplit(ncfiles(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(7),'_',/extract)
    sdate=result3(0)
    print,sdate
;
; declare time period arrays on first day
;
    if icount eq 0L then begin
       zbar_zt=fltarr(nlv,nfile)
       tbar_zt=fltarr(nlv,nfile)
       ubar_zt=fltarr(nlv,nfile)
       vbar_zt=fltarr(nlv,nfile)
       qbar_zt=fltarr(nlv,nfile)
       pvbar_zt=fltarr(nlv,nfile)
       sdate_all=strarr(nfile)
       icount=1
    endif
    sdate_all(ifile)=sdate
;
; save Arctic means on each day
;
    zbar_zt(*,ifile)=zbar
    tbar_zt(*,ifile)=tbar
    ubar_zt(*,ifile)=ubar
    vbar_zt(*,ifile)=vbar
    qbar_zt(*,ifile)=qbar
    pvbar_zt(*,ifile)=pvbar

    jump:
endfor          ; loop over time steps

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='zt_geos4_polarT_1prof_'+sdate0+'-'+sdate1+'_'+slon+'_'+slat+'.sav',zbar_zt,tbar_zt,ubar_zt,$
     vbar_zt,qbar_zt,pvbar_zt,sdate_all,p,nfile,slon,slat,glon,glat,ilon,ilat

plotit:
;restore,'zt_geos4_polarT_1prof_20060101-20060228_18.75E_77.5N.sav
restore,'zt_geos4_polarT_1prof_20060101-20060228_93.75E_77.5N.sav
;restore,'zt_geos4_polarT_1prof_20060101-20060228_206.25E_77.5N.sav
;restore,'zt_geos4_polarT_1prof_20060101-20060228_281.25E_77.5N.sav
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)

    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='zt_geos4_polarT_1prof_'+sdate0+'-'+sdate1+'_'+slon+'_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot time-altitude section of T at a single Arctic location
;
    erase
    level=175.+5.*findgen(24)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*mcolor/nlvls
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    syr=strmid(sdate_all,0,4)
    smn=strmid(sdate_all,4,2)
    sdy=strmid(sdate_all,6,2)
    xindex=where(sdy eq '15',nxticks)
    xlabs=smn(xindex)
    contour,transpose(tbar_zt),findgen(nfile),p,charsize=2,/noeras,yrange=[110,min(p)],/ylog,charthick=2,$
            ytitle='Pressure (hPa)',levels=level,c_color=col1,/cell_fill,color=0,title='GEOS-4 Temperature '+slon+' '+slat,$
            xrange=[0.,nfile-1],zticks=6,xticks=nxticks-1,xtickname=' '+strarr(nxticks+1)
    contour,transpose(tbar_zt),findgen(nfile),p,levels=level(0:*:2),color=0,/follow,/noeras,/overplot,min_value=0.,thick=2
;   contour,transpose(ubar_zt),findgen(nfile),p,charsize=1.5,/noeras,/overplot,levels=[10,20,40,60,80,100,120],color=0
;   contour,transpose(ubar_zt),findgen(nfile),p,charsize=1.5,/noeras,/overplot,levels=-1.*[120,100,80,60,40,20,10],color=mcolor,c_linestyle=5
;   x2d=0.*transpose(pbar_zt)
    for ii=0L,nxticks-1L do xyouts,xindex(ii),200.,smon(long(xlabs(ii))-1),/data,color=0,charsize=2,charthick=2,alignment=0.5
    for ii=0L,nxticks-1L do begin
        plots,xindex(ii),110.
        plots,xindex(ii),140.,/continue,color=0,thick=2,/data
    endfor
    xyouts,5.,200.,syr(0),/data,color=0,charsize=2,charthick=2,alignment=0.5

    imin=min(level)
    imax=max(level)
    ymnb=yorig(0)-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=2,charthick=2
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim zt_geos4_polarT_1prof_'+sdate0+'-'+sdate1+'_'+slon+'_'+slat+'.ps -rotate -90 '+$
                           'zt_geos4_polarT_1prof_'+sdate0+'-'+sdate1+'_'+slon+'_'+slat+'.jpg'
    endif
end
