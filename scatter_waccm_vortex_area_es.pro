;
; save GEOS-5 vortex area as a function of altitude and time about ES days zeros
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
xorig=[0.15]
yorig=[0.25]
xlen=0.8
ylen=0.5
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
restore, '/Users/harvey/Desktop/Harvey_etal_2014/Post_process/WACCM_ES_daily_max_T_Z.sav'
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
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
        ncfile0=dir+sdate+'_3D_dyn.nc3'
        print,ncfile0
        ncid=ncdf_open(ncfile0)
        result0=ncdf_inquire(ncid)
        for idim=0,result0.ndims-1 do begin
            ncdf_diminq,ncid,idim,name,dim
            if name eq 'number_of_latitudes' then nr=dim
            if name eq 'number_of_longitudes' then nc=dim
            if name eq 'number_of_levels' then nth=dim
        endfor
        for ivar=0,result0.nvars-1 do begin
            result=ncdf_varinq(ncid,ivar)
            if result.name eq 'latitude' or result.name eq 'longitude' or result.name eq 'theta' or result.name eq 'MARK' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
            if result.name eq 'latitude' then alat=data
            if result.name eq 'longitude' then alon=data
            if result.name eq 'theta' then th=data
            if result.name eq 'MARK' then mark2=data
;           print,ivar,result.name,min(data),max(data)
        endfor
        ncdf_close,ncid

      if iday eq 0L then begin
         area_zt_nc4=fltarr(ndays,nth)
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
          if index(0) ne -1 then area_zt_nc4(icount,thlev)=100.*total(area(index))/hem_area
      endfor
      jumpday:
      icount=icount+1L
endfor          ; loop over days

saveit:
;
; plot altitude-time series of Arctic vortex area
;
yy=strmid(sfile,6,2)
index=where(yy ne '')
y1='20'+string(FORMAT='(I2.2)',min(yy(index)))
y2='20'+string(FORMAT='(I2.2)',max(yy(index)))
;
; save file
;
save,file='vortex_area_ES_'+esdate(30+ndays*(ies))+'_waccm.sav',area_zt_nc4,th,sfile,y1,y2
index=where(area_zt_nc4 le 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.

plotit:
restore,'vortex_area_ES_'+esdate(30+ndays*(ies))+'_waccm.sav'
kday=n_elements(SFILE)
if ies eq 0L then begin
   print,th
   ilev=0L
   rlev=0.
   read,'Enter theta surface ',rlev
   index=where(rlev eq th)
   ilev=index(0)
   slev=strcompress(long(rlev),/remove_all)
endif
;
; save postscript version
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/times,filename='../Figures/scatter_waccm_vortex_area_ES_event_'+sevent+'_'+slev+'K.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot vortex area scatter plot
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=30
col1=1+indgen(nlvls)*icolmax/nlvls
index=where(area_zt_nc4 eq 0.)
if index(0) ne -1L then area_zt_nc4(index)=0./0.
area_zt_nc4=smooth(area_zt_nc4,5,/nan)
area_lev=reform(area_zt_nc4(*,ilev))
plot,1.+findgen(40),th,color=0,xtitle='WACCM Vortex Area (% NH) '+slev+' K',xrange=[0.,40.],yrange=[500.,10000.],/noeras,ytitle='Polar Cap Mean Stratopause Theta (K)',$
     title='WACCM ES event '+sevent+' ('+esdate(30+ndays*(ies))+')',/nodata
for k=0L,ndays-1 do oplot,[area_lev(k),area_lev(k)],[dailymeanTHETA(ies,k),dailymeantHETA(ies,k)],color=((ndays-k)/float(ndays))*mcolor,psym=8
maxtheta=reform(dailymeanTHETA(ies,*))
index=where(maxtheta gt 0. and area_lev gt 0.)
result=correlate(area_lev(index),maxtheta(index))
xyouts,20.,1000.,'r= '+strcompress(result),/data,charsize=2,charthick=2,color=0
imin=-30
imax=30
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
      xtitle='Days From ES Onset',charthick=2,charsize=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim ../Figures/scatter_waccm_vortex_area_ES_event_'+sevent+'_'+slev+'K.ps -rotate -90 ../Figures/scatter_waccm_vortex_area_ES_event_'+sevent+'_'+slev+'K.jpg'
endif

endfor	; loop over ES events
end
