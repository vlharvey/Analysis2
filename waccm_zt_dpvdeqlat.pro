;
; save the product of the average wind speed at the vortex edge
; multiplied by the normalised gradient of PV wrt Equivalent latitude
; as a function of altitude and time
;
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
lstmn=11L & lstdy=1L & lstyr=0L 
ledmn=5L & leddy=1L & ledyr=0L
lstday=0L & ledday=0L
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.35]
xlen=0.8
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
dir='/aura3/data/WACCM_data/Datfiles/waccm_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
ifiles=[$
'waccm_nh_files_03-04.fil']
nyear=n_elements(ifiles)
for iyear=0L,nyear-1L do begin
ifile=''
close,1
openr,1,ifiles(iyear)
ndays=0L
readf,1,ndays
sfile=strarr(ndays)
for iday=0L,ndays-1L do begin
    readf,1,ifile
    if iday gt 0L then sfile(iday)=ifile
    dum=findfile(dir+ifile+'.nc3')
    if dum(0) eq '' then goto,skipit
    ncid=ncdf_open(dir+ifile+'.nc3')
    print,'opening ',dir+ifile+'.nc3'
    if iday eq 0L then begin
       nr=0L & nc=0L & nth=0L
       ncdf_diminq,ncid,0,name,nr
       ncdf_diminq,ncid,1,name,nc
       ncdf_diminq,ncid,2,name,nth
       alon=fltarr(nc)
       alat=fltarr(nr)
       th=fltarr(nth)
       pv2=fltarr(nr,nc,nth)
       u2=fltarr(nr,nc,nth)
       v2=fltarr(nr,nc,nth)
       ncdf_varget,ncid,0,alon
       ncdf_varget,ncid,1,alat
       ncdf_varget,ncid,2,th

       pv_zt=fltarr(ndays,nth)
       pvgavg_zt=fltarr(ndays,nth)
       prod_zt=fltarr(ndays,nth)
       sfile=strarr(ndays)
       sfile(iday)=ifile
       dum=transpose(u2(*,*,0))
       lon=0.*dum
       lat=0.*dum
       for i=0,nc-1 do lat(i,*)=alat
       for j=0,nr-1 do lon(*,j)=alon
    endif
    ncdf_varget,ncid,3,pv2
    ncdf_varget,ncid,5,u2
    ncdf_varget,ncid,6,v2
    ncdf_close,ncid
;
; loop over theta
;
    for thlev=0,nth-1 do begin
        u1=transpose(u2(*,*,thlev))
        v1=transpose(v2(*,*,thlev))
        speed1=sqrt(u1^2+v1^2)
        pv1=transpose(pv2(*,*,thlev))
        index=where(abs(pv1) gt 1000.)
        if index(0) ne -1 then pv1(index)=pv1(index)/0.
        pv1=smooth(pv1,3,/NaN)
        elat1=calcelat2d(pv1,alon,alat)
if max(pv1) gt 1000. then goto,skiplev
;
; integrate wind speed around PV isopleths
;
        nbins=37
        dy=2.5
        latmin=0.
        latmax=90.
        elatbin=latmin+dy*findgen(nbins)
        speedbin=-999.+0.*fltarr(nbins)               ; average windspeed per bin
        speedsigbin=-999.+0.*fltarr(nbins)               ; average windspeed per bin
        pvbin=-999.+0.*elatbin
        pvsigbin=0.*elatbin
        for n=0,nbins-2 do begin
            t=where(abs(pv1) lt 1000. and $
                    lat ge latmin and elat1 ge elatbin(n) and elat1 lt elatbin(n+1),it)
; check latmin.  make sure bins are resolved (do not intersect latmin)
            if (it gt 2) then begin
                result=moment(pv1(t))
                pvbin(n)=result(0)
                pvsigbin(n)=sqrt(result(1))

                result=moment(speed1(t))
                speedbin(n)=total(speed1(t))/float(it)
                speedsigbin(n)=sqrt(result(1))

                if min(lat(t))-latmin le dy then begin
                   pvbin(n)=-999.
                   pvsigbin(n)=-999.
                   speedbin(n)=-999.
                   speedsigbin(n)=-999.
                   goto,jumpnhbin
                endif
            endif
            jumpnhbin:
        endfor        ; loop over bins
        s=where(lat ge latmin and elat1 ge elatbin(nbins-1),is)
        if is gt 2 then begin
           result=moment(pv1(s))
           pvbin(nbins-1)=result(0)
           pvsigbin(nbins-1)=sqrt(result(1))

           result=moment(speed1(s))
           speedbin(nbins-1)=total(speed1(s))/float(is)
           speedsigbin(nbins-1)=sqrt(result(1))

           if min(lat(s))-latmin le dy then begin
              pvbin(nbins-1L)=-999. 
              pvsigbin(nbins-1L)=-999. 
              speedbin(nbins-1)=-999.
              speedsigbin(nbins-1)=-999.
           endif
        endif
        dpvbin=-999.+0.*pvbin
        for i=0,nbins-2L do $
            if pvbin(i) ne -999. and pvbin(i+1) ne -999. then dpvbin(i)=pvbin(i+1)-pvbin(i)
        if pvbin(nbins-1) ne -999. and pvbin(nbins-2) ne -999. then $
           dpvbin(nbins-1)=pvbin(nbins-1)-pvbin(nbins-2)
;
; impose Nash filter poleward of 80deg (and add new one Equatorward of lat0)
;
        lat0=70.
        index=where(elatbin ge lat0 and speedbin gt 0. and dpvbin ne -999.) 
        if index(0) ne -1L then begin
           speedbin(index)=speedbin(index)*(90.-elatbin(index))/30.
           dpvbin(index)=dpvbin(index)*(90.-elatbin(index))/30.
           pvsigbin(index)=pvsigbin(index)*(90.-elatbin(index))/30.
        endif
        lat0=25.
        if th(thlev) lt 600. then lat0=45.
        index=where(elatbin le lat0 and speedbin gt 0. and dpvbin ne -999.)
        if index(0) ne -1L then begin
           speedbin(index)=speedbin(index)*(elatbin(index))/(2.*lat0)
           dpvbin(index)=dpvbin(index)*(elatbin(index))/(2.*lat0)
           pvsigbin(index)=pvsigbin(index)*(elatbin(index))/(2.*lat0)
        endif
;
; normalise
;
        index=where(dpvbin ne -999.)
        if index(0) eq -1L then goto,skiplev
        dpvbin(index)=dpvbin(index)/max(dpvbin(index))
        pvsigbin(index)=pvsigbin(index)/max(pvsigbin(index))
;
; current index is the maximum in the product of the average wind speed
; and the average PV gradient, both are a function of Equivalent latitude
;
        prod=0.*dpvbin
        index=where(speedbin ne -999. and dpvbin ne -999.)
        prod(index)=dpvbin(index)*speedbin(index)
        index=where(prod eq max(prod))
        if index(0) ne -1L then begin
           pv_zt(iday,thlev)=pvbin(index(0))
           pvgavg_zt(iday,thlev)=dpvbin(index(0))
           prod_zt(iday,thlev)=prod(index(0))
        endif

if thlev lt 0 then begin
erase
!type=2^2+2^3
set_viewport,.1,.4,.6,.9
plot,elatbin,pvbin,xrange=[0.,90.],yrange=[0.,max(pvbin)],psym=1,/noeras,$
     title='PV '+string(th(thlev)),charsize=2
set_viewport,.5,.8,.6,.9
smax=max(speedbin)
index=where(speedbin gt 0.)
plot,elatbin(index),speedbin(index),xrange=[0.,90.],title='Avg Speed',$
     yrange=[0.,smax],color=icolmax*.3,thick=2
oplot,elatbin(index),speedsigbin(index),psym=1
set_viewport,.1,.4,.2,.5
smax=max(dpvbin)
index=where(dpvbin gt 0.)
plot,elatbin(index),dpvbin(index),xrange=[0.,90.],title='dPV/dEqlat',$
     yrange=[0.,1.0],color=icolmax*.9,thick=2
oplot,elatbin,pvsigbin
set_viewport,.5,.8,.2,.5
smax=max(prod)
index=where(prod gt 0.)
plot,elatbin(index),prod(index),xrange=[0.,90.],title='dPV*Speed',$
     yrange=[0.,smax],thick=2
endif
skiplev:
      endfor
skipit:
endfor          ; loop over days
y1=strmid(sfile(0),0,4)
y2=strmid(sfile(ndays-1),0,4)
;
; save file
;
save,file='waccm_vortex_strength_'+y1+'_'+y2+'.sav',prod_zt,th,sfile

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='waccm_pvg_'+y1+'-'+y2+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot altitude-time section of vortex strength indicator, look at prod=dpv*speed
;
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
plot,[1,ndays,ndays,1,1],[400.,400.,2000.,2000.,400.],min_value=0.,$
      xrange=[1,ndays],yrange=[400.,2000.],/nodata,charsize=2,$
      ytitle='Theta (K)',title='WACCM '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
kindex=where(strmid(sfile,6,2) eq '15',nxtick)
xmon=long(strmid(sfile(kindex),4,2))
for i=0,nxtick-1 do begin
    xlab=smon(xmon(i)-1)
    plots,kindex(i)+1,350.
    plots,kindex(i)+1,400.,/continue,/data
    xyouts,kindex(i)+1,225.,xlab,/data,alignment=0.5,charsize=3
endfor
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=5.*findgen(nlvls)
index=where(prod_zt eq 0.)
if index(0) ne -1 then prod_zt(index)=-999.
contour,prod_zt,1.+findgen(ndays),th,levels=level,/fill,$
        /cell_fill,/overplot,c_color=col1,min_value=-999.
contour,prod_zt,1.+findgen(ndays),th,levels=level,c_color=0,$
        /follow,/overplot,min_value=-999.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
      xtitle='Edge Avg Wind Speed x Normalised PV Gradient',charsize=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim waccm_pvg_'+y1+'-'+y2+'.ps -rotate -90 waccm_pvg_'+y1+'-'+y2+'.jpg'
endif

endfor
end
