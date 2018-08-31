;
; save 3-event average GEOS-5 vortex area as a function of altitude and time about ES days zeros
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=9L & lstdy=1L & lstyr=2009L 
ledmn=5L & leddy=1L & ledyr=2010L
y1=strcompress(lstyr,/remove_all)
y2=strcompress(ledyr,/remove_all)

lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.65,0.175]
xlen=0.8
ylen=0.3
cbaryoff=0.1
cbarydel=0.01
set_plot,'ps'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=mcolor
   !p.background=mcolor
endif
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
restore, '/Users/harvey/Desktop/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
dir1='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
result=size(MAXHEIGHTTHETA)
nevents=result(1)
for iES = 0L, nevents - 1L do begin
    sevent=strtrim(strcompress(string(format='(I3.2)',ies+1)),2)
    icount=0L
    ndays=61
    sdays=strarr(ndays)
    for iday =0L, ndays-1L do begin
        sdays(iday)=string(format='(I3.2)',iday-30)
        sdate=esdate(iday+60*ies)
        if sdate eq '' then goto,jumpday        ; missing day
;
; read daily file
;
        dum=findfile(dir1+sdate+'_AVG.V01.nc3')
        if dum ne '' then ncfile0=dir1+sdate+'_AVG.V01.nc3'
        if dum eq '' then ncfile0=dir+sdate+'_AVG.V01.nc3'
print,ncfile0
        ncid=ncdf_open(ncfile0)
        ncdf_diminq,ncid,0,name,nr
        ncdf_diminq,ncid,1,name,nc
        ncdf_diminq,ncid,2,name,nth
        alon=fltarr(nc)
        alat=fltarr(nr)
        th=fltarr(nth)
        mark2=fltarr(nr,nc,nth)
        ncdf_varget,ncid,0,alon
        ncdf_varget,ncid,1,alat
        ncdf_varget,ncid,2,th
        ncdf_varget,ncid,10,mark2
        ncdf_close,ncid

      if ies eq 0L and iday eq 0L then begin
         print,th
         ilev=0L
         rlev=2000.
;        read,'Enter theta surface ',rlev
         index=where(rlev eq th)
         ilev=index(0)
         slev=strcompress(long(rlev),/remove_all)
;
; save postscript version
;
         if setplot eq 'ps' then begin
            set_plot,'ps'
            xsize=nxdim/100.
            ysize=nydim/100.
            device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                   /bold,/color,bits_per_pixel=8,/times,filename='../Figures/geos5_vortex_area_ES_event_avg_'+slev+'K.ps'
            !p.charsize=1.25
            !p.thick=2
            !p.charthick=5
            !p.charthick=5
            !y.thick=2
            !x.thick=2
         endif
 
         area_zt_nc4=fltarr(ndays,nth)
         narea_zt_nc4=fltarr(ndays,nth)
         sfile=strarr(ndays)
         dum=transpose(mark2(*,*,0))
         lon=0.*dum
         lat=0.*dum
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
         kcount=1L
      endif
      sfile(icount)=sdate
;
; loop over theta
;
      for thlev=0,nth-1 do begin
          mark1=transpose(mark2(*,*,thlev))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=area_zt_nc4(icount,thlev)+100.*total(area(index))/hem_area
          if index(0) ne -1 then narea_zt_nc4(icount,thlev)=narea_zt_nc4(icount,thlev)+1.
      endfor
      jumpday:
      icount=icount+1L
endfor          ; loop over days
endfor  ; loop over ES events
;
; plot altitude-time series of mean Arctic vortex area
;
area_zt_nc4=area_zt_nc4/narea_zt_nc4
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=17
col1=1+indgen(nlvls)*icolmax/nlvls
level=2.*findgen(nlvls)
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
contour,area_zt_nc4,1.+findgen(ndays),th,color=0,xtitle='Days From ES Onset',thick=6,yrange=[500.,6000.],/noeras,ytitle='Theta (K)',/cell_fill,$
     xticks=ndays/10,xtickname=sdays(0:ndays-1:10),title='Composite MLS ES Event',levels=level,c_color=col1
contour,area_zt_nc4,1.+findgen(ndays),th,levels=level,c_color=mcolor,/follow,/overplot,c_labels=1+0*level,min_value=0
plots,31,500
plots,31,6000.,/continue,color=0,thick=3
compositestrat=fltarr(ndays)
ncompositestrat=fltarr(ndays)
for ies=0,2 do begin
    meantheta=reform(dailymeanTHETA(ies,*))
    index=where(meantheta ne 0.)
    compositestrat(index)=compositestrat(index)+meantheta(index)
    ncompositestrat(index)=ncompositestrat(index)+1.
endfor
compositestrat=compositestrat/ncompositestrat
smoothcompositestrat=smooth(compositestrat,5,/nan)
compositestrat(31:60)=smoothcompositestrat(31:60)
oplot,1.+findgen(ndays),compositestrat,psym=8,color=mcolor*.9
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='GEOS5 Vortex Area (% NH)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;
; scatterplot
;
xmn=xorig(1)
xmx=xorig(1)+0.34
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=30
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
area_lev=reform(area_zt_nc4(*,ilev))
esday=-30.+findgen(ndays)
plot,findgen(31),th,color=0,xtitle='GEOS5 Vortex Area (% NH) at '+slev+' K',xrange=[0.,30.],yrange=[500.,6000.],/noeras,ytitle='Stratopause Theta (K)',/nodata
plots,0,500.
plots,30,6000.,/continue,color=0,thick=2
for k=0L,ndays-1 do begin
    kk=k-30L
    if abs(kk) gt 10L then oplot,[area_lev(k),area_lev(k)],[compositestrat(k),compositestrat(k)],color=(k/float(ndays))*mcolor,psym=8
endfor
index=where(abs(esday) le 10)
area_lev(index)=0.
maxtheta=compositestrat
index=where(maxtheta gt 0. and area_lev gt 0.)
result=correlate(area_lev(index),maxtheta(index))
xyouts,20.,1000.,'r= '+string(format='(f4.2)',result),/data,color=0
imin=-30
imax=30
ymnb=yorig(1)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='Days From ES Onset'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

set_viewport,0.61,0.95,ymn,ymx
!type=2^2+2^3
rprofile=fltarr(nth)
for ilev=0L,nth-1L do begin
    area_lev=reform(area_zt_nc4(*,ilev))
    index=where(abs(esday) le 10)
    area_lev(index)=0.
    index=where(maxtheta gt 0. and area_lev gt 0.)
    rprofile(ilev)=correlate(area_lev(index),maxtheta(index))
endfor
plot,rprofile,th,color=0,xtitle='Correlation Coefficient',xrange=[-1.,1.],yrange=[500.,6000.],/noeras,ytitle='Theta (K)',thick=5
plots,0,500
plots,0,6000.,/continue,color=0,thick=3,linestyle=5

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ../Figures/geos5_vortex_area_ES_event_avg_'+slev+'K.ps -rotate -90 ../Figures/geos5_vortex_area_ES_event_avg_'+slev+'K.jpg'
endif

end
