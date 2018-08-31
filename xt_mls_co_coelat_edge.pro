;
; Hovmoller of CO + CO Elat Edge
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3
@rd_sdwaccm4_nc3
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
nxdim=750
nydim=750
xorig=[0.25]
yorig=[0.25]
xlen=0.5
ylen=0.7
cbaryoff=0.15
cbarydel=0.01
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
nyear=n_elements(syear)
start_year=long(syear)
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
dirm='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirw='/Volumes/cloud/data/WACCM_data/Datfiles_SD/'
dir='/Volumes/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_'

for nn=start_year(0),start_year(nyear-2L) do begin

lstmn=11         ; NH
lstdy=1
lstyr=nn
ledmn=3         ; NH
leddy=1
ledyr=nn
lstday=0
ledday=0
;
; loop over years
;
for iyear=lstyr,ledyr do begin
for ilat=60,60 do begin ;50,80,5 do begin

kcount=0
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

z = stddat(lstmn,lstdy,iyear,lstday)
z = stddat(ledmn,leddy,iyear+1,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
;
; from old code
;
;
; declare (nfile) arrays
;
nrr=91L
yeq=findgen(nrr)
index=where(yeq mod 2 eq 0,nrr2)
yeq2=yeq(index)                         ; 2 degree bins
nfile=kday
highlat_time=fltarr(nfile)      ; PDF peak of high Elat geographical lats (use these to determine if CO can be used to mark the vortex today)
lowlat_time=fltarr(nfile)       ; PDF peak of low Elat geographical lats
hlatpdf_time=fltarr(nfile,nrr2)
llatpdf_time=fltarr(nfile,nrr2)
elatedge_time=fltarr(nfile)
coedge_time=fltarr(nfile)
lowlat_elatedge_time=fltarr(nfile)
lowlat_coedge_time=fltarr(nfile)
sfelatedge_time=fltarr(nfile)
markcoelatedge_time=fltarr(nfile)
marksfelatedge_time=fltarr(nfile)
pvelatedge_time=fltarr(nfile)
nashedge_time=fltarr(nfile)
sdate_time=strarr(nfile)
;
; Compute initial Julian date
;
iyr = iyear
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
n=icount

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      sdate_all(icount)=sdate
      print,sdate

      ifiles=file_search(dir+sdate+'.sav',count=nfile)
      if nfile eq 0L then goto,skipmls
;
; read gridded MLS data on pressure
; CO_GRID         FLOAT     = Array[144, 96, 37]
; GP_GRID         FLOAT     = Array[144, 96, 55]
; H2O_GRID        FLOAT     = Array[144, 96, 55]
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; N2O_GRID        FLOAT     = Array[144, 96, 37]
; O3_GRID         FLOAT     = Array[144, 96, 55]
; PMLS            FLOAT     = Array[37]
; PMLS2           FLOAT     = Array[55]
; TP_GRID         FLOAT     = Array[144, 96, 55]
;
        ifile=dir+sdate+'.sav'
        dum=findfile(ifile)
        if dum(0) eq '' then goto,skipmls
        restore,ifile
        print,'restored '+ifile
;
; MERRA
;
        ifile=dirm+sdate+'.nc3'
        rd_merra_nc3,ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
            u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
;       print,iflag,ifile
        sp2=sqrt(u2^2.+v2^2.)
;
; SD-WACCM
;
;       ifile=dirw+'f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h5.'+sdate+'.nc3'
;       rd_sdwaccm4_nc3,ifile,nc,nr,nth,alon,alat,th,$
;          pv2,p2,z2,u2,v2,q2,qdf2,mark2,sf2,h2o2,n2o2,o32,iflg
;
; declare time period arrays on first day
;
        if kcount eq 0L then begin
          rth=1800.
;         print,th
;         read,'Enter desired theta ',rth
          index=where(abs(rth-th) eq min(abs(rth-th)))
          ith=index(0)
          mlsco_xt=fltarr(nc,kday)
          mlselat_xt=fltarr(nc,kday)
          kcount=1
;
; postscript file
;
        if setplot eq 'ps' then begin
           lc=0
           xsize=nxdim/100.
           ysize=nydim/100.
           set_plot,'ps'
           device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                  /bold,/color,bits_per_pixel=8,/helvetica,filename='xt_mls_co_coelat_edge_'+syr+'_'+sth+'.ps'
           !p.charsize=1.25
           !p.thick=2
           !p.charthick=5
           !y.thick=2
           !x.thick=2
        endif

        endif
        sth=strcompress(long(th(ith)),/remove_all)+'K'
        pv1=transpose(reform(pv2(0:nr-1,*,ith)))
        sp1=transpose(reform(sp2(0:nr-1,*,ith)))
        mark1=transpose(reform(mark2(0:nr-1,*,ith)))
        qdf1=transpose(reform(qdf2(0:nr-1,*,ith)))
        sf1=transpose(reform(sf2(0:nr-1,*,ith)))
        q1=transpose(reform(q2(0:nr-1,*,ith)))
;
; choose pressure
;
        rpress=2.0
;       print,pmls
;       read,'Enter desired pressure ',rpress
        index=where(abs(rpress-pmls) eq min(abs(rpress-pmls)))
        ilev=index(0)
        spress=strcompress(pmls(ilev),/remove_all)+'hPa'
        nc=n_elements(lon)
        nr=n_elements(lat)
        nl=n_elements(pmls)
;       alat=lat
;       alon=lon
        dum=co_grid(*,*,0)
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
        x2d=fltarr(nc,nr)
        y2d=fltarr(nc,nr)
        for i=0,nc-1 do y2d(i,*)=alat
        for j=0,nr-1 do x2d(*,j)=alon
;
; horizontal CO gradient
;
        co1=reform(co_grid(*,*,ilev))*1.e6
;       co1=smooth(co1,3,/edge_truncate)
;
; CO at rlat
;
        yindex=where(abs(rlat-alat) eq min(abs(rlat-alat)))
        mlsco_xt(*,icount)=reform(co1(*,yindex(0)))



;
; old code
;
        comax=max(abs(co1(*,nr/2:-1)))
        cograd1=co1*0.0
        for j = 0, nr-1 do begin
            jm1=j-1
            jp1=j+1
            if j eq 0 then jm1=0
            if j eq 0 then dy2=(alat(1)-alat(0))*!pi/180.
            if j eq nr-1 then jp1=nr-1
            if j eq nr-1 then dy2=(alat(nr-1)-alat(nr-2))*!pi/180.
            if (j gt 0 and j lt nr-1) then dy2=(alat(jp1)-alat(jm1))*!pi/180.
            csy=cos(alat(j)*!pi/180.)
            for i = 0, nc-1 do begin
                ip1 = i+1
                im1 = i-1
                if i eq 0 then im1 = nc-1
                if i eq 0 then dx2 = (alon(1)-alon(0))*!pi/180.
                if i eq nc-1 then ip1 = 0
                if i eq nc-1 then dx2 = (alon(0)-alon(nc-1))*!pi/180.
                if (i gt 0 and i lt nc-1) then dx2=(alon(ip1)-alon(im1))*!pi/180.
        
                dqdx = (co1(ip1,j)-co1(im1,j))/(dx2*csy)
                dqdy = (co1(i,jp1)-co1(i,jm1))/dy2
                cograd1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)	;/abs(co1(i,j))	;/comax
;
; allowing negative on vortex interior is useful to isolate edge but if vortex is offset from the pole dCO/dy<0 should be >0
;
;               if (dqdy le 0.0) then cograd1(i,j) = -sqrt(dqdx*dqdx+dqdy*dqdy)$
;                                                    *abs(co1(i,j))/comax
;
; without normalization
;
;               cograd1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)       ;/abs(co1(i,j))
;               if (dqdy le 0.0) then cograd1(i,j) = -1.0*cograd1(i,j)

            endfor
        endfor
;
; poles are bad and neighbooring lats
;
        cograd1(*,0)=0./0.
        cograd1(*,1)=0./0.
        cograd1(*,nr-1)=0./0.
        cograd1(*,nr-2)=0./0.
;
; horizontal PV gradient
;
        pvgrad1=pv1*0.0
        for j = 0, nr-1 do begin
            jm1=j-1
            jp1=j+1
            if j eq 0 then jm1=0
            if j eq 0 then dy2=(alat(1)-alat(0))*!pi/180.
            if j eq nr-1 then jp1=nr-1
            if j eq nr-1 then dy2=(alat(nr-1)-alat(nr-2))*!pi/180.
            if (j gt 0 and j lt nr-1) then dy2=(alat(jp1)-alat(jm1))*!pi/180.
            csy=cos(alat(j)*!pi/180.)
            for i = 0, nc-1 do begin
                ip1 = i+1
                im1 = i-1
                if i eq 0 then im1 = nc-1
                if i eq 0 then dx2 = (alon(1)-alon(0))*!pi/180.
                if i eq nc-1 then ip1 = 0
                if i eq nc-1 then dx2 = (alon(0)-alon(nc-1))*!pi/180.
                if (i gt 0 and i lt nc-1) then dx2=(alon(ip1)-alon(im1))*!pi/180.
        
                dqdx = (pv1(ip1,j)-pv1(im1,j))/(dx2*csy)
                dqdy = (pv1(i,jp1)-pv1(i,jm1))/dy2
                pvgrad1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)	;/abs(pv1(i,j))  ;/pvmax
;
; allowing negative on vortex interior is useful to isolate edge but if vortex is offset from the pole dCO/dy<0 should be >0
;
;               if (dqdy le 0.0) then pvgrad1(i,j) = -sqrt(dqdx*dqdx+dqdy*dqdy)$
;                                                    *abs(pv1(i,j))/pvmax
;
; without normalization
;
;               pvgrad1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)       ;/abs(pv1(i,j))
;               if (dqdy le 0.0) then pvgrad1(i,j) = -1.0*pvgrad1(i,j)

            endfor
        endfor
;
; poles are bad and neighbooring lats
;
        pvgrad1(*,0)=0./0.
        pvgrad1(*,1)=0./0.
        pvgrad1(*,nr-1)=0./0.
        pvgrad1(*,nr-2)=0./0.
;
; CO Elat - set SH values to -1
;
        index=where(lat lt 0.)
        cosave=co1
        cosave(index)=-1.*cosave(index)
        elat=calcelat2d(cosave,alon,alat)
        mlselat_xt(*,icount)=reform(elat(*,yindex(0)))

;
; SF Elat - multiply NH by -1 to get max in vortex and set SH values to -1.e12
;
        index=where(lat lt 0.)
        sfsave=sf1+abs(min(sf1))
;set_viewport,.1,.9,.1,.9
;erase
;map_set,0,0,0,/contin,/grid,/noeras,color=0
;imin=min(sfsave)
;imax=max(sfsave)
;nlvls=26L
;col1=1+(indgen(nlvls)/float(nlvls))*mcolor
;contour,sfsave,alon,alat,levels=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls),c_color=col1,/follow,/noeras
;stop
;       sfsave(index)=-1.*sfsave(index)	; make SH negative
;       sfsave=1./sfsave

index=where(lat gt -100.)
sfsave(index)=-1.*sfsave(index)	; want max in the NH vortex like PV
index=where(lat lt 0.)		; eliminate SH 
sfsave(index)=-1.e12
set_viewport,.1,.9,.5,.9
erase
;map_set,0,0,0,/contin,/grid,/noeras,color=0
;imin=min(sfsave)
;imax=max(sfsave)
;nlvls=26L
;col1=1+(indgen(nlvls)/float(nlvls))*mcolor
;contour,sfsave,alon,alat,levels=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls),c_color=col1,/follow,/noeras
;stop

        sfelat=calcelat2d(sfsave,alon,alat)
;erase
;map_set,0,0,0,/contin,/grid,/noeras,color=0
;imin=min(elat)
;imax=max(elat)
;nlvls=26L
;col1=1+(indgen(nlvls)/float(nlvls))*mcolor
;contour,elat,alon,alat,levels=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls),c_color=col1,/follow,/noeras
;stop
;
; PV elat
;
        pvelat=calcelat2d(pv1,alon,alat)
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
index=where(y2h ne 0.)
;plot,yeq2(index),y2h(index),color=mcolor*.9,thick=5
        highindex=where(y2h eq max(y2h))
        highco=yeq2(highindex(0))
        highlat_time(n)=highco
print,'High Elat Lat ',highco
        index=where(y2d gt 0. and elat gt 0. and elat le 60.,nn)
oplot,lon(index),lat(index),psym=8,color=0
        y2=histogram(y2d(index),min=0,max=90,binsize=2)/float(nn)
        llatpdf_time(n,*)=y2
;index=where(y2 ne 0.)
;set_viewport,.1,.9,.1,.4
;plot,yeq2,y2h,color=mcolor*.9,thick=5
;oplot,yeq2,y2,color=mcolor*.3,thick=5


        lowindex=where(y2 eq max(y2))
        lowco=yeq2(lowindex(0))
        lowlat_time(n)=lowco
print,'Low Elat Lat ',lowco
        d=lowco/highco
;stop
erase
;
; where low elat latitudes eq 0. and high elat latitudes ne 0.
;
index=where(y2 eq 0. and y2h ne 0.)
if index(0) eq -1L then print,'NO VORTEX TODAY?'
if index(0) ne -1L then print,'vortex home? ',yeq2(index)
;d_time(n)=d
print,'Low Elat lats/High Elat lats = ',d
;
; mean CO and CO gradient within 1 degree spaced elat bins
; mean CO gradient
;
       cobin=fltarr(nrr)
       cogradbin=fltarr(nrr)
       for i=0L,nrr-1L do begin
           ip1=i+1
           if i eq nrr-1 then ip1=nrr-1
           im1=i-1
           if i eq 0L then im1=0
           index=where(y2d gt 0. and elat ge yeq(im1) and elat lt yeq(ip1))
           if index(0) ne -1L then cobin(i)=mean(co1(index))
           if index(0) ne -1L then cogradbin(i)=mean(cograd1(index))
       endfor
;
; derivative of Elat wrt CO
;
       delatlevs=smooth(deriv(yeq,cobin),7,/edge_truncate)
;
; mean PV within 1 degree spaced elat
; mean PV gradient
;
       pvbin=fltarr(nrr)
       pvgradbin=fltarr(nrr)
       spbin=fltarr(nrr)
       for i=0L,nrr-1L do begin
           ip1=i+1
           if i eq nrr-1 then ip1=nrr-1
           im1=i-1
           if i eq 0L then im1=0
           index=where(y2d gt 0. and pvelat ge yeq(im1) and pvelat lt yeq(ip1))
           if index(0) ne -1L then pvbin(i)=mean(pv1(index))
           if index(0) ne -1L then pvgradbin(i)=mean(pvgrad1(index))
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
;
; flag nodes in QDFbin
;
       flag=0*qdfbin
       for i=0L,nrr-2L do if qdfbin(i)*qdfbin(i+1) lt 0. then flag(i)=1.
;
; flag speed maxima
;
       spflag=0*qdfbin
       for i=1L,nrr-2L do if spbinsf(i) gt spbinsf(i-1) and spbinsf(i) gt spbinsf(i+1) then spflag(i)=1.
;
; where both QDF node and local max in speed
;
       index=where(flag eq 1. and spflag eq 1.)
       sfedge=-99.
       if index(0) ne -1L then sfedge=max(yeq(index))
       if index(0) ne -1L then sfelatedge_time(n)=max(yeq(index))
;
; yikes, why is PV vs Elat so noisy?
;
;      pvbin=smooth(pvbin,3,/edge_truncate)
;
; derivative of PV Elat wrt PV
;
       dpvelatlevs=smooth(deriv(yeq,pvbin),7,/edge_truncate)
;
; apply "sloping filter" equal to 1 at 80 and 0 at 90 (1.125-yeq/80.)/0.125
;
       index=where(yeq ge 80.)
       delatlevs(index)=delatlevs(index)*((1.125-yeq(index)/80.)/0.125)
       dpvelatlevs(index)=dpvelatlevs(index)*((1.125-yeq(index)/80.)/0.125)
       spbin(index)=spbin(index)*((1.125-yeq(index)/80.)/0.125)
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
; most equatorward local maximum in the CO-Elat first derivative
; require positive slope over 2 points prior to max and negative slope following max
;
flag=0.*delatlevs
ilim=2
for i=ilim,n_elements(delatlevs)-ilim-1L do begin
    if delatlevs(i) gt delatlevs(i-1) and delatlevs(i) gt delatlevs(i+1) then begin
       if delatlevs(i) gt delatlevs(i-ilim) and delatlevs(i) gt delatlevs(i+ilim) then flag(i)=1.
    endif
endfor
;
; and of the CO contours co-located with maximum CO gradients, require those gradients to be some fraction of the strength of the maximum gradient
;
index=where(flag eq 1 and delatlevs ge max(delatlevs)*0.5)	;yeq gt 30.)
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
; end old code
;

skipmls:
      icount=icount+1L
      n=n+1
goto,jump

plotit:

sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
;xindex=where(sday eq '01' or sday eq '05' or sday eq '10' or sday eq '15' or sday eq '20' or sday eq '25',nxticks)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7       ; ticks outward

       imin=min(mlsco_xt)
       imax=max(mlsco_xt)
       nlvls=26L
       col1=1+(indgen(nlvls)/float(nlvls))*mcolor
       level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
mlsco_xtshift=0.*mlsco_xt
mlselat_xtshift=0.*mlsco_xt
index=where(alon gt 180.)
mlsco_xtshift(index-min(index),*)=mlsco_xt(index,*)
mlselat_xtshift(index-min(index),*)=mlselat_xt(index,*)
index=where(alon lt 180.)
mlsco_xtshift(index+max(index),*)=mlsco_xt(index,*)
mlselat_xtshift(index+max(index),*)=mlselat_xt(index,*)

contour,mlsco_xtshift,alon-180.,findgen(kday),/noeras,yrange=[kday-1,0.],xrange=[-180.,180.],$
      charsize=2,color=0,xtitle='Longitude',/fill,c_color=col1,charthick=2,$
      levels=level,yticks=nxticks-1,ytickname=xlabs,ytickv=xindex,min_value=-99.
;contour,mlselat_xtshift,alon-180,findgen(kday),levels=30+10*findgen(6),color=mcolor,/follow,/overplot,thick=1
;
; superimpose COElat edge
;
for ii=0,kday-1L do begin
    today=reform(mlselat_xt(*,ii))
    index=where(abs(today-lowlat_elatedge_time(ii)) le 1.)
    if index(0) ne -1L then oplot,alon(index)-180,ii+0*index,psym=8,color=mcolor
endfor

       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle=sth+' CO (ppmv)',charthick=2
       ybox=[0,10,10,0,0]
       x2=imin
       dx=(imax-imin)/(float(nlvls)-1)
       for j=1,nlvls-1 do begin
           xbox=[x2,x2,x2+dx,x2+dx,x2]
           polyfill,xbox,ybox,color=col1(j)
           x2=x2+dx
       endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim xt_mls_co_coelat_edge_'+syr+'_'+sth+'.ps -rotate -90 xt_mls_co_coelat_edge_'+syr+'_'+sth+'.png'
    endif
endfor  ; loop over latitudes
endfor  ; loop over years

endfor
end
