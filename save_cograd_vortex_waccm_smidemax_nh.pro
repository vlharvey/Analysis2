;
; this version is based on free-running WACCM - Ethan's 300 year "SmidEmax" run
; modified from:
; upgrades: MLS v4 and MERRA2
; DJF
; choose equatorward-most CO Elat local maximum in the first derivative
; CO elat
;
@rd_waccm_smidemax_nc3
@calcelat2d

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
index=where(yeq mod 2 eq 0,nrr2)
yeq2=yeq(index)				; 2 degree bins

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

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.1,0.55,0.1,0.55]
yorig=[0.55,0.55,0.15,0.15]
cbaryoff=0.02
cbarydel=0.01
xlen=0.3
ylen=0.3
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
;syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017']
iyears=1+indgen(300)
syear=string(format='(i3.3)',iyears)
nyear=n_elements(syear)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; get file listing
;
dir='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO2x1SmidEmax_yBWCN/3d_CO2x1SmidEmax_yBWCN_'
;
; loop over years and months
;
for iyear=0L,nyear-1L do begin
for imon=0L,nmon-1L do begin
;   if imon gt 2 and imon lt 9 then goto,skipmon		; ONDJFM 
    if imon ne 3 and imon ne 8 then goto,skipmon		; April and September
    ifiles=file_search(dir+syear(iyear)+smon(imon)+'??.nc3',count=nfile)
    if nfile eq 0L then goto,skipmon
;
; declare (nfile) arrays
;
    highlat_time=fltarr(nfile)	; PDF peak of high Elat geographical lats (use these to determine if CO can be used to mark the vortex today)
    lowlat_time=fltarr(nfile)	; PDF peak of low Elat geographical lats
    hlatpdf_time=fltarr(nfile,nrr2)
    llatpdf_time=fltarr(nfile,nrr2)
    elatedge_time=fltarr(nfile)
    coedge_time=fltarr(nfile)
    lowlat_elatedge_time=fltarr(nfile)
    lowlat_elatinner_time=fltarr(nfile)
    lowlat_elatouter_time=fltarr(nfile)
    lowlat_coedge_time=fltarr(nfile)
    sfelatedge_time=fltarr(nfile)
    markcoelatedge_time=fltarr(nfile)
    marksfelatedge_time=fltarr(nfile)
    pvelatedge_time=fltarr(nfile)
    nashedge_time=fltarr(nfile)
    nashinner_time=fltarr(nfile)
    nashouter_time=fltarr(nfile)
    novortex_flag=-9999.+0.*fltarr(nfile)
    sdate_time=strarr(nfile)
    firstday_flag=0L
;
; loop over files
;
    icount=0L
    FOR n=0l,nfile-1l DO BEGIN
        result=strsplit(ifiles(n),'.',/extract)
        result2=result(-2)
        result3=strsplit(result2,'_',/extract)
        sdate=result3(-1)
        sdate_time(n)=sdate
;
; read WACCM netcdf theta data to get both CO and PV and don't need MERRA. everything is more consistent now with both on theta

        ifile=ifiles(n)
print,ifile
        dum=findfile(ifile)
        if dum(0) eq '' then goto,jumpstep
        rd_waccm_smidemax_nc3,ifile,nc,nr,nth,alon,alat,th,pv2,p2,u2,v2,qdf2,co2,z2,sf2,mark2
;
; wind speed
;
        sp2=sqrt(u2^2.+v2^2.)
;
; chores
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
; on first day of each month
; declare CO marker fields, both 3D and 4D
; declare vortex edge arrays
;
        delatlevs2d=fltarr(nrr,nth)
        spbin2d=fltarr(nrr,nth)
        markco=0.*co2
        if n eq 0L then markco4d=fltarr(nc,nr,nth,nfile)
        if n eq 0L then begin
           lowlat_elatedge_2d=fltarr(nfile,nth)
           lowlat_elatinner_2d=fltarr(nfile,nth)
           lowlat_elatouter_2d=fltarr(nfile,nth)
           sfelatedge_2d=fltarr(nfile,nth)
           novortex_flag_2d=fltarr(nfile,nth)
           delatlevs3d=fltarr(nfile,nrr,nth)
           spbin3d=fltarr(nfile,nrr,nth)
           hlatpdf_time_3d=fltarr(nfile,nrr2,nth)
           llatpdf_time_3d=fltarr(nfile,nrr2,nth)
           nashelatedge_2d=fltarr(nfile,nth)
           nashinner_2d=fltarr(nfile,nth)
           nashouter_2d=fltarr(nfile,nth)
        endif
;
; loop over WACCM theta levels
;
kcount=0L
for kk=nth-1L,0L,-1L do begin   ; all theta levels from the bottom-up
;
; strip out theta level
;
    pv1=transpose(reform(pv2(0:nr-1,*,kk)))
    sp1=transpose(reform(sp2(0:nr-1,*,kk)))
    mark1=transpose(reform(mark2(0:nr-1,*,kk)))
    qdf1=transpose(reform(qdf2(0:nr-1,*,kk)))
    sf1=transpose(reform(sf2(0:nr-1,*,kk)))
    z1=transpose(reform(z2(0:nr-1,*,kk)))
    speed1=transpose(reform(sp2(0:nr-1,*,kk)))
;
; horizontal CO
;
    co1=transpose(reform(co2(0:nr-1,*,kk)))*1.e6
    mark1co=0.*co1
;
; CO Elat - set SH values to -1
;
    index=where(lat lt 0.)
    cosave=co1
    cosave(index)=-1.*cosave(index)
    elat=calcelat2d(cosave,alon,alat)
;
; SF Elat - multiply NH by -1 to get max in vortex and set SH values to -1.e12
;
    index=where(lat lt 0.)
    sfsave=sf1+abs(min(sf1))
    index=where(lat gt -100.)
    sfsave(index)=-1.*sfsave(index)	; want max in the NH vortex like PV
    index=where(lat lt 0.)		; eliminate SH 
    sfsave(index)=-1.e12
    sfelat=calcelat2d(sfsave,alon,alat)
;
; PV elat
;
    pvelat=calcelat2d(pv1,alon,alat)
;
; retain rlev for polar plot
;
;   rlev=3000.
    rlev=5000.
;   rlev=6000.
    slev=strcompress(long(rlev),/r)+'K'
    if abs(th(kk)-rlev) eq min(abs(th-rlev)) then begin
       co1save=co1
       elatsave=elat
       pvelatsave=pvelat
       speed1save=speed1
       qdf1save=qdf1
       sf1save=sf1
       sfelatsave=sfelat
       sp1save=sp1
       mark1save=mark1
    endif
;
; save elat min associated with SF vortex previously marked
;
    markcoelatedge_time(n)=-99.
    marksfelatedge_time(n)=-99
    index=where(y2d gt 0. and mark1 gt 0.)
    if index(0) ne -1L then begin
       markcoelatedge_time(n)=min(elat(index))
       marksfelatedge_time(n)=min(sfelat(index))
    endif
;
; Geographic Latitude PDFs for Elat >80 and <60
;
    index=where(y2d gt 0. and elat ge 80.,nn)
    y2h=histogram(y2d(index),min=0,max=90,binsize=2)/float(nn)
    hlatpdf_time(n,*)=y2h
    highindex=where(y2h eq max(y2h))
    highco=yeq2(highindex(0))
    highlat_time(n)=highco
    index=where(y2d gt 0. and elat gt 0. and elat le 60.,nn)
    y2=histogram(y2d(index),min=0,max=90,binsize=2)/float(nn)
    llatpdf_time(n,*)=y2
    lowindex=where(y2 eq max(y2))
    lowco=yeq2(lowindex(0))
    lowlat_time(n)=lowco
    d=lowco/highco
hlatpdf_time_3d(n,*,kk)=y2h
llatpdf_time_3d(n,*,kk)=y2
;
; if the range of geographic latitudes for the vortex is the entire hemisphere - or if the highest CO values are confined to low lats - then the vortex does not exist
;erase
;!type=2^2+2^3
;plot,yeq2,y2h,color=0,thick=5,position=[.1,.1,.9,.9],title=string(th(kk))
;oplot,yeq2,y2,color=250,thick=3
index=where(y2h ne 0.)			; PDF of geographic vortex latitudes
;if index(0) ne -1L then print,th(kk),mean(z1)/1000.,' Geographic Vortex Latitudes ',min(yeq2(index)),max(yeq2(index))
;
; where low elat latitudes eq 0. and high elat latitudes ne 0.
;
    index=where(y2 eq 0. and y2h ne 0.)		; will give latitudes where two distributions do not overlap
;   if index(0) eq -1L then print,'NO VORTEX TODAY?'
;   if index(0) ne -1L then print,'vortex home? ',yeq2(index)
    if index(0) ne -1L then novortex_flag(n)=yeq2(index(0))
    novortex_flag_2d(n,kk)=novortex_flag(n)
;
; mean CO within 1 degree spaced elat bins
;
    cobin=fltarr(nrr)
    for i=0L,nrr-1L do begin
        ip1=i+1
        if i eq nrr-1 then ip1=nrr-1
        im1=i-1
        if i eq 0L then im1=0
        index=where(y2d gt 0. and elat ge yeq(im1) and elat lt yeq(ip1))
        if index(0) ne -1L then cobin(i)=mean(co1(index))
    endfor
    if abs(th(kk)-rlev) eq min(abs(th-rlev)) then cobinsave=cobin
;
; this needs to be AFTER cobinsave is initialized!
; new logic - this needs to be added to MLS/MERRA codes too
; if the range of geographic latitudes for the vortex is the entire hemisphere - or if the highest CO values are confined to low lats - then the vortex does not exist
;
index=where(y2h ne 0.)
if min(yeq2(index)) lt 20. then goto,jumplev
if max(yeq2(index)) lt 45. then goto,jumplev
;
; derivative of Elat wrt CO
;
    delatlevs=smooth(deriv(yeq,cobin),7,/edge_truncate)
    d2elatlevs=smooth(deriv(yeq,delatlevs),7,/edge_truncate)
;
; mean PV within 1 degree spaced elat
;
    pvbin=fltarr(nrr)
    spbin=fltarr(nrr)
    for i=0L,nrr-1L do begin
        ip1=i+1
        if i eq nrr-1 then ip1=nrr-1
        im1=i-1
        if i eq 0L then im1=0
        index=where(y2d gt 0. and pvelat ge yeq(im1) and pvelat lt yeq(ip1))
        if index(0) ne -1L then pvbin(i)=mean(pv1(index))
        if index(0) ne -1L then spbin(i)=mean(sp1(index))
    endfor
;
; mean QDF and speed within SF bins
;
       sfbin=fltarr(nrr)
       qdfbin=fltarr(nrr)
       spbinsf=fltarr(nrr)
       for i=0L,nrr-1L do begin
           ip1=i+1
           if i eq nrr-1 then ip1=nrr-1
           im1=i-1
           if i eq 0L then im1=0
           index=where(y2d gt 0. and sfelat ge yeq(im1) and sfelat lt yeq(ip1))
           if index(0) ne -1L then sfbin(i)=mean(sf1(index))
           if index(0) ne -1L then qdfbin(i)=mean(qdf1(index))
           if index(0) ne -1L then spbinsf(i)=mean(sp1(index))
       endfor
       qdfbin=smooth(qdfbin,3,/edge_truncate)
       spbinsf=smooth(spbinsf,3,/edge_truncate)
       spbin2d(*,kk)=spbinsf	;/max(spbinsf)
;
; flag nodes in QDFbin
;
       flag=0*qdfbin
       for i=0L,nrr-2L do if qdfbin(i)*qdfbin(i+1) lt 0. then flag(i)=1.
;
; smooth to allow for overlaps within -1 to +1
;
       flag=smooth(flag,5) 
;
; flag speed maxima
;
       spflag=0*qdfbin
       spbinsf=smooth(spbinsf,7,/edge_truncate)
       for i=1L,nrr-2L do if spbinsf(i) gt spbinsf(i-1) and spbinsf(i) gt spbinsf(i+1) then spflag(i)=1.
;
; where both QDF node and local max in speed
; use only poleward most jet - do not require node in QDF
;
       index=where(spflag eq 1. and spbinsf eq max(spbinsf))	; and flag gt 0.)
;      if th(kk) gt min(th) then index=where(spflag eq 1. and abs(yeq-sfelatedge_2d(n,kk-1)) le 10.)    ; within 10 degrees of edge below
       sfedge=-99.
       if index(0) ne -1L then sfedge=max(yeq(index))
;      if index(0) ne -1L then print,'SF edge ',sfedge
       if index(0) ne -1L then sfelatedge_time(n)=max(yeq(index))
       if index(0) ne -1L then sfelatedge_2d(n,kk)=sfedge
;
; yikes, why is PV vs Elat so noisy?
;
;      pvbin=smooth(pvbin,3,/edge_truncate)
;
; derivative of PV Elat wrt PV
;
       dpvelatlevs=smooth(deriv(yeq,pvbin),7,/edge_truncate)
       d2pvelatlevs=smooth(deriv(yeq,dpvelatlevs),7,/edge_truncate)
;
; apply "sloping filter" equal to 1 at 80 and 0 at 90 (1.125-yeq/80.)/0.125
;
       index=where(yeq ge 80.)
       delatlevs(index)=delatlevs(index)*((1.125-yeq(index)/80.)/0.125)
       dpvelatlevs(index)=dpvelatlevs(index)*((1.125-yeq(index)/80.)/0.125)
       spbin(index)=spbin(index)*((1.125-yeq(index)/80.)/0.125)
;
; add a sloping filter in the tropics
;
       index=where(yeq le 30.)
       delatlevs(index)=delatlevs(index)*(yeq(index)/30.)
       dpvelatlevs(index)=dpvelatlevs(index)*(yeq(index)/30.)
       spbin(index)=spbin(index)*(yeq(index)/30.)
;
; retain all levels
;
       delatlevs2d(*,kk)=delatlevs/max(delatlevs)
;
; Nash product of wind and dpv/elat
;
       prod=dpvelatlevs*spbin
;
; absolute maximum in the first derivative - as Nash
;
       index=where(delatlevs eq max(delatlevs))
       elatedge=yeq(index(0))		; Elat value of the edge
       coedge=cobin(index(0))		; CO value of the edge
       index=where(dpvelatlevs eq max(dpvelatlevs))
       pvelatedge=yeq(index(0))		; PV Elat value of the edge
       pvedge=pvbin(index(0))		; PV value of the edge
       index=where(prod eq max(prod))
       nashedge=yeq(index(0))           ; Elat value of the "Nash" edge
;
; inner and outer: look within 10 degrees of max in the gradient
;
       index=where(abs(yeq-nashedge) gt 10.)
       d2pvelatlevs(index)=0.
       index=where(d2pvelatlevs ne 0.)
       if index(0) eq -1L then stop,'pv second derivative all zero'
       dum=d2pvelatlevs(index)
       ydum=yeq(index)
       index=where(dum eq max(dum))
       nashouter_time(n)=ydum(index(0))
       index=where(dum eq min(dum))
       nashinner_time(n)=ydum(index(0))
;
; save 2d
;
       nashelatedge_2d(n,kk)=nashedge
       nashinner_2d(n,kk)=nashinner_time(n)
       nashouter_2d(n,kk)=nashouter_time(n)
;
; most equatorward local maximum in the CO-Elat first derivative
; require positive slope over 2 points prior to max and negative slope following max
;
       dflag=0.*delatlevs
       ilim=2
       for i=ilim,n_elements(delatlevs)-ilim-1L do begin
           if delatlevs(i) gt delatlevs(i-1) and delatlevs(i) gt delatlevs(i+1) then begin
              if delatlevs(i) gt delatlevs(i-ilim) and delatlevs(i) gt delatlevs(i+ilim) then dflag(i)=1.
           endif
       endfor
;
; and of the CO contours co-located with maximum CO gradients, require those gradients to be some fraction of the strength of the maximum gradient
;
;      index=where(dflag eq 1 and delatlevs ge max(delatlevs)*0.5)						; Harvey et al (2015) criteria
       index=where(dflag eq 1 and delatlevs ge max(delatlevs)*0.4 and yeq gt 30.)				; relax 50% of max criteria and reinstate min-elat threshhold

;      index=where(dflag eq 1 and delatlevs ge max(delatlevs)*0.95)						; max co vs elat
;      if th(kk) gt min(th) then index=where(dflag eq 1. and abs(yeq-lowlat_elatedge_2d(n,kk-1)) le 10.)	; within 10 degrees of edge below

       edgecandidates=[-99.]
       coedgecandidates=[-99.]
       if index(0) ne -1L then edgecandidates=yeq(index)
       if index(0) ne -1L then coedgecandidates=cobin(index)
;print,'Candidates ',edgecandidates
;
; broadcast edges
;
;print,'CO Elat Edge ',elatedge
;print,'CO Elat Edge lowlat ',min(edgecandidates)
;print,'PV Elat Edge ',pvelatedge
;print,'Nash Edge ',nashedge
;if abs(th(kk)-rlev) eq min(abs(th-rlev)) then begin
;erase
;!type=2^2+2^3
;plot,yeq,cobin,color=0,thick=5,position=[.1,.1,.9,.9],title=string(th(kk))
;axis,/yax,yrange=[min(delatlevs),max(delatlevs)],/save,color=250
;oplot,yeq,delatlevs,color=250,thick=3
;oplot,yeq,max(delatlevs)*0.4+0.*delatlevs,color=0,psym=2
;stop
;endif
;
; save edges
;
       elatedge_time(n)=elatedge
       coedge_time(n)=coedge
       index=where(edgecandidates eq min(edgecandidates))
       lowlat_elatedge_time(n)=edgecandidates(index)
       lowlat_coedge_time(n)=coedgecandidates(index)
;print,'CO Value for Elat edge: ',coedge
;print,'CO Value for low lat Elat edge: ',coedgecandidates(index)
       pvelatedge_time(n)=pvelatedge
       nashedge_time(n)=nashedge
;
; inner and outer boundaries
;
       if lowlat_elatedge_time(n) ne -99. then begin
          index=where(abs(yeq-lowlat_elatedge_time(n)) gt 10.)
          d2elatlevs(index)=0.
          index=where(d2elatlevs ne 0.)
          if index(0) eq -1L then stop,'second derivative all zero'
          dum=d2elatlevs(index)
          ydum=yeq(index)
          index=where(dum eq max(dum))
          lowlat_elatouter_time(n)=ydum(index)
          index=where(dum eq min(dum))
          lowlat_elatinner_time(n)=ydum(index)
       endif
;
; retain edges at all levels
;
       lowlat_elatedge_2d(n,kk)=lowlat_elatedge_time(n)
       lowlat_elatinner_2d(n,kk)=lowlat_elatinner_time(n)
       lowlat_elatouter_2d(n,kk)=lowlat_elatouter_time(n)
;
; fill marker array
;
       index=where(elat ge lowlat_elatedge_time(n))
       if index(0) ne -1L then mark1co(index)=1.
       markco(*,*,kk)=mark1co
;
; end loop over levels
;
jumplev:
       kcount=kcount+1L
endfor
; 
; save postscript version
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='Polar_plots/polar_daily_coelat+cograd_waccm-smidemax_'+sdate+'_2d_nh.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
erase
xyouts,.4,.475,sdate,color=0,charsize=2,charthick=2,/normal
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=min(cobinsave)
imax=max(cobinsave)
nlvls=26L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
map_set,90,0,-90,/stereo,/contin,/grid,/noerase,color=0,charsize=1.5
dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=co1save
dum(nc,*)=dum(0,*)
x=fltarr(nc+1)
x(0:nc-1)=alon
x(nc)=x(0)+360.

contour,dum,x,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot
index=where(abs(th-rlev) eq min(abs(th-rlev)))
dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=speed1save
dum(nc,*)=dum(0,*)
;contour,dum,x,alat,levels=60+20*findgen(4),/noeras,charsize=2,c_color=[190,210,250,0],/foll,/overplot,thick=3
dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=elatsave
dum(nc,*)=dum(0,*)
contour,dum,x,alat,levels=lowlat_elatedge_2d(n,index(0)),/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=8
;
; add co marker
;
;index=where(abs(th-rlev) eq min(abs(th-rlev)))
;mark=reform(markco(*,*,index(0)))
;contour,mark,alon,alat,levels=[0.1],/noeras,charsize=2,color=mcolor*.9,/foll,/overplot,c_labels=[1],thick=8
;stop

markco4d(*,*,*,n)=markco

ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,charthick=2,xtitle='WACCM CO '+slev
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor

dum=mean(z2,dim=1)
zmean=mean(dum,dim=1)/1000.     ; global averaged altitude (km)

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)+0.05
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(delatlevs2d gt 0. and finite(delatlevs2d) eq 1)
imin=min(delatlevs2d(index))
imax=max(delatlevs2d)
level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
if max(level) ne 0. then begin
contour,delatlevs2d,yeq,zmean,levels=level,c_color=col1,/fill,/noeras,charsize=1,charthick=2,ytitle='Approx. Alt (km)',xtitle='Equivalent Latitude',xrange=[20,90],$
        yrange=[min(zmean),100.],thick=5,color=0
dum=reform(lowlat_elatedge_2d(n,*))
index=where(dum gt 10.)
if index(0) ne -1L then oplot,dum(index),zmean(index),psym=8,color=0
if index(0) ne -1L then oplot,dum(index),zmean(index),color=0,thick=3
ymnb=ymn -cbaryoff-0.05
ymxb=ymnb+cbarydel
set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,charthick=2,xtitle='2D CO PDF'
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor
endif
;
; QDF/Mark
;
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=-400.
imax=400.       ;max(max(qdf1))
level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
map_set,90,0,-90,/stereo,/contin,/grid,/noerase,color=0,charsize=1.5
dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=qdf1save
dum(nc,*)=dum(0,*)
contour,dum,x,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot

dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=sf1save
dum(nc,*)=dum(0,*)
contour,dum,x,alat,nlevels=20,/noeras,charsize=2,color=0,/foll,/overplot,thick=2,c_labels=0*level
;map_set,90,0,-90,/stereo,/contin,/noerase,color=mcolor
;contour,sfelatsave,alon,alat,levels=10+10*findgen(8),/noeras,charsize=2,color=.9*mcolor,/foll,/overplot,c_labels=[1,1,1,1,1,1,1,1]                        ; all contours

dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=sp1save
dum(nc,*)=dum(0,*)
contour,dum,x,alat,levels=60+20*findgen(4),/noeras,charsize=2,c_color=[190,210,250,0],/foll,/overplot,thick=3
;contour,mark1save,alon,alat,levels=[0.1],/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=5,c_linestyle=5
index=where(abs(th-rlev) eq min(abs(th-rlev)))
dum=fltarr(nc+1,nr)
dum(0:nc-1,*)=sfelatsave
dum(nc,*)=dum(0,*)
contour,dum,x,alat,levels=sfelatedge_2d(n,index(0)),/noeras,charsize=2,color=mcolor,/foll,/overplot,thick=8
;
; Nash edge
;
loadct,0
dum(0:nc-1,*)=pvelatsave
dum(nc,*)=dum(0,*)
contour,dum,x,alat,levels=nashelatedge_2d(n,index(0)),/noeras,charsize=2,color=150,/foll,/overplot,thick=8
loadct,39

ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,charthick=2,xtitle='WACCM QDF '+slev
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor
;
; normalized speed f(elat)
;
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)+0.05
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(spbin2d gt 0. and finite(spbin2d) eq 1)
imin=0.         ; min(spbin2d(index))
imax=100.       ; max(spbin2d)
level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
contour,spbin2d,yeq,zmean,levels=level,c_color=col1,/fill,/noeras,charsize=1,charthick=2,ytitle='Approx. Alt (km)',xtitle='Equivalent Latitude',xrange=[20,90],$
        yrange=[min(zmean),100.],thick=5,color=0
dum=reform(sfelatedge_2d(n,*))
index=where(dum gt 10.)
if index(0) ne -1L then oplot,dum(index),zmean(index),psym=8,color=0
if index(0) ne -1L then oplot,dum(index),zmean(index),color=0,thick=3

loadct,0
dum=reform(nashelatedge_2d(n,*))
index=where(dum gt 10.)
if index(0) ne -1L then oplot,dum(index),zmean(index),psym=8,color=150
if index(0) ne -1L then oplot,dum(index),zmean(index),color=150,thick=3
loadct,39

ymnb=ymn -cbaryoff-0.05
ymxb=ymnb+cbarydel
set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,charthick=2,xtitle='Wind Speed'
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor

index=where(th ge 500. and th le 100000.)
;print,'Vortex Flag ',reform(novortex_flag_2d(n,index))
;
; retain elat binned dCO and speed
;
delatlevs3d(n,*,*)=delatlevs2d
spbin3d(n,*,*)=spbin2d

; Close PostScript file and return control to X-windows
       if setplot ne 'ps' then stop
       if setplot eq 'ps' then begin
          device,/close
          spawn,'convert -trim Polar_plots/polar_daily_coelat+cograd_waccm-smidemax_'+sdate+'_2d_nh.ps -rotate -90 '+$
                              'Polar_plots/polar_daily_coelat+cograd_waccm-smidemax_'+sdate+'_2d_nh.jpg'
          spawn,'rm -f Polar_plots/polar_daily_coelat+cograd_waccm-smidemax_'+sdate+'_2d_nh.ps'
       endif
jumpday:
;endfor
icount=icount+1L
jumpstep:
endfor	; loop over files
;
; save monthly file
;
save,filename='/atmos/harvey/WACCM_data/Datfiles/Datfiles_Ethan_600yr/CO_Vortex_data/daily_waccm-smidemax_coelatedge+sfelatedge_'+syear(iyear)+smon(imon)+'_2d_nh.sav',$
		lowlat_elatedge_2d,$		; CO Elat Edge (equatorward-most maximum in 1st derivative that is half the magnitude of absolute maximum) *** Use this as CO Edge ***
 		lowlat_elatinner_2d,$		; CO Elat inner boundary
 		lowlat_elatouter_2d,$		; CO Elat outer
 		sfelatedge_2d,$			; SF Elat Edge (stongest jet in SF Elat space) *** Use this as SF Edge ***
 		delatlevs3d,$			; dCO/Elat on which the lowlat elatedge is based
 		spbin3d,$			; speed/Elat on which the SF elatedge is based
 		sdate_time,th,yeq,$
 		novortex_flag_2d,$		; -9999 when low and high lat PDFs of Elat distribution overlap. Set to Elat when they don't
                markco4d,$			; 4d MLS marker array based on CO
                hlatpdf_time_3d,$		; geographic lat PDF of CO > 80
                llatpdf_time_3d,$		; geographic lat PDF of CO < 60
                nashelatedge_2d,nashinner_2d,nashouter_2d	; Nash edge, inner, outer as a function of day and level
skipmon:
endfor	; loop over months
endfor	; loop over years
end
