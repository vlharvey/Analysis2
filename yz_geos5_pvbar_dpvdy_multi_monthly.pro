;
; plot monthly mean zonal mean PV and dPV/dy in GEOS-5
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.15]
yorig=[0.15]
xlen=0.7
ylen=0.7
cbaryoff=0.065
cbarydel=0.02
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
!noeras=1
gdir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
stimes=[$
'_AVG.V01.']
month=['January','February','March','April','May','June','July','August','September','October','November','December']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']

month=['January','July']
nmon=['01','07']
for imon=0L,n_elements(nmon)-1L do begin
    spawn,'ls yz_geos5_pvbar_dpvdy_monthly_*'+nmon(imon)+'.sav',ifiles
    nn=n_elements(ifiles)
    for ii=0L,nn-1L do begin
;
; pvzm,dpvdyzm,sfzm,tzm,uzm,markzm,alat,th,icount
;
        restore,ifiles(ii)
        print,'restored '+ifiles(ii)
        if ii eq 0L then begin
           pvzm_tot=0.*pvzm
           dpvdyzm_tot=0.*pvzm
           sfzm_tot=0.*pvzm
           tzm_tot=0.*pvzm
           uzm_tot=0.*pvzm
           markzm_tot=0.*pvzm
        endif
        pvzm_tot=pvzm_tot+pvzm
        dpvdyzm_tot=dpvdyzm_tot+dpvdyzm
        sfzm_tot=sfzm_tot+sfzm
        tzm_tot=tzm_tot+tzm
        uzm_tot=uzm_tot+uzm
        markzm_tot=markzm_tot+markzm
    endfor
    pvzm=pvzm_tot/float(nn)
    dpvdyzm=dpvdyzm_tot/float(nn)
    sfzm=sfzm_tot/float(nn)
    tzm=tzm_tot/float(nn)
    uzm=uzm_tot/float(nn)
    markzm=markzm_tot/float(nn)
;
; normalize by PV
;
;    dpvdyzm=dpvdyzm/pvzm
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_geos5_pvbar_dpvdy_multi_monthly_'+nmon(imon)+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
       !p.thick=2.0
       !p.charsize=2.0
    endif
erase
loadct,39
xyouts,.25,.9,'GEOS-5 '+month(imon)+' 2004 - 2008',/normal,charsize=2,charthick=2,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
tlevel=[-1.,-0.5,-0.1,-0.05,-0.01,-0.005,-0.001,-0.0005,-0.0001,-0.00005,-0.00001,0.,$
        0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.]
;tlevel=-.5+0.05*findgen(21)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
xmin=30 & xmax=90.
if imon eq 1L then begin
   xmin=-90. & xmax=-30.
endif
contour,dpvdyzm,alat,th,/noeras,xrange=[xmin,xmax],yrange=[500.,4000.],charsize=1.5,color=0,$
        ytitle='Theta (K)',title='Zonal Mean dPV/dy',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel lt 0.)
contour,dpvdyzm,alat,th,levels=tlevel(index),color=mcolor,/follow,/overplot
index=where(tlevel gt 0.)
contour,dpvdyzm,alat,th,levels=tlevel(index),color=0,/follow,/overplot
ulevel=-140.+20.*findgen(25)
nlvls=n_elements(ulevel)
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(ulevel lt 0.)
contour,uzm,alat,th,levels=ulevel(index),color=0,/follow,/overplot,c_linestyle=5,thick=2
index=where(ulevel gt 0.)
contour,uzm,alat,th,levels=ulevel(index),color=mcolor,/follow,/overplot,thick=2
;contour,tzm,alat,th,levels=[240.,250.,260.,270.],color=0,/follow,/overplot,thick=3
;contour,markzm,alat,th,levels=[0.5],color=mcolor,/follow,/overplot,thick=5
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Km!nkg!u-1!ns!u-1!n'
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
       device,/close
       spawn,'convert -trim yz_geos5_pvbar_dpvdy_multi_monthly_'+nmon(imon)+'.ps '+$
             ' -rotate -90 yz_geos5_pvbar_dpvdy_multi_monthly_'+nmon(imon)+'.jpg'
    endif
jumpmonth:
endfor	; loop over months
end
