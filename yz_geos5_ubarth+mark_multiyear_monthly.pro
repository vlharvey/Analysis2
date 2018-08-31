;
; multi-year monthly mean zonal mean
; GEOS-5 version
;
; plot daily Ubar from theta data
; superimpose zonal mean marker
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.125
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS5?0.MetO.'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS5?0.MetO.'
month=['January','February','March','April','May','June','July','August','September','October','November','December']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
years=[2003,2004,2005,2006,2007,2008,2009,2010,2011]

for ii=0L,11L do begin
;
;***Read GEOS data
;
    smn=string(FORMAT='(I2.2)',ii+1)
    print,smn
goto,quick

    spawn,'ls '+dir+'????'+smn+'??_AVG.V01.nc3',ifiles
    if ifiles(0) eq '' then goto,jumpmonth
    nday=n_elements(ifiles)
    for icount=0,nday-1L do begin
        rd_geos5_nc3_meto,ifiles(icount),nc,nr,nth,alon,alat,th,$
                 pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
        index=where(mark2 lt 0.)
        if index(0) ne -1L then mark2(index)=0.	; zero out anticyclones
;
; compute temperature
;
        t2=0.*u2
        for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
;
; zonal mean zonal wind
;
        if icount eq 0L then begin
           tzm=fltarr(nr,nth)
           tzm_all=fltarr(nday,nr,nth)
           uzm=fltarr(nr,nth)
           uzm_all=fltarr(nday,nr,nth)
           markzm=fltarr(nr,nth)
           markzm_all=fltarr(nday,nr,nth)
        endif

        for k=0L,nth-1L do begin
        for j=0L,nr-1L do begin
            index=where(t2(j,*,k) ne 0.,np)
            if index(0) ne -1L then begin
            tzm(j,k)=tzm(j,k)+total(t2(j,index,k))/float(np)
            tzm_all(icount,j,k)=total(t2(j,index,k))/float(np)
            uzm(j,k)=uzm(j,k)+total(u2(j,index,k))/float(np)
            uzm_all(icount,j,k)=total(u2(j,index,k))/float(np)
;           markzm(j,k)=markzm(j,k)+max(mark2(j,*,k))
;           markzm_all(icount,j,k)=max(mark2(j,*,k))
            markzm(j,k)=markzm(j,k)+total(mark2(j,index,k))/float(np)
            markzm_all(icount,j,k)=total(mark2(j,index,k))/float(np)
            endif
        endfor
        endfor
     endfor	; loop over daily files
;
; average
;
tzm=tzm/float(nday)
uzm=uzm/float(nday)
markzm=markzm/float(nday)
;
quick:
restore,'yz_geos5_ubarth+mark_multiyear_'+smn+'.sav'	;,alat,th,tzm,uzm,markzm

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_geos5_ubarth+mark_multiyear_'+smn+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot zonal mean zonal wind
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
level=-100.+10.*findgen(nlvls)
index=where(uzm eq 0.)
if index(0) ne -1L then uzm(index)=-9999.
contour,uzm,alat,th,/noeras,xrange=[-90.,90.],yrange=[400.,3000.],charsize=2,color=0,$
      ytitle='Potential Temperature (K)',xticks=6,/cell_fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title='GEOS-5 '+month(ii),charthick=2
index=where(level gt 0.)
contour,uzm,alat,th,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,thick=2
index=where(level lt 0.)
contour,uzm,alat,th,levels=level(index),color=mcolor,/follow,/overplot,min_value=-9999.,c_linestyle=5,thick=2
contour,markzm,alat,th,levels=[0.1,.5,.9],color=0,/follow,/overplot,thick=15,c_labels=1+0*findgen(3)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Ubar (m/s)',charsize=1.5,charthick=2
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
   spawn,'convert -trim yz_geos5_ubarth+mark_multiyear_'+smn+'.ps -rotate -90 yz_geos5_ubarth+mark_multiyear_'+smn+'.jpg'
endif
;save,file='yz_geos5_ubarth+mark_multiyear_'+smn+'.sav',alat,th,tzm,uzm,markzm

jumpmonth:
endfor	; loop over months
end
