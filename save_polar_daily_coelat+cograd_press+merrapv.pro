;
; DJF
; choose equatorward-most CO Elat local maximum in the first derivative
; CO elat
; CO horizontal gradient - dCO/dx + dCO/dy
;
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
;index=where(yeq mod 2 eq 0,nrr)
;yeq=yeq(index)				; 2 degree bins

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
yorig=[0.55,0.55,0.05,0.05]
cbaryoff=0.02
cbarydel=0.01
xlen=0.4
ylen=0.4
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
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
dirm='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
dirw='/Volumes/cloud/data/WACCM_data/Datfiles_SD/'
dir='/Volumes/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_'

for iyear=0L,nyear-1L do begin
for imon=0L,nmon-1L do begin
     if imon gt 3 and imon lt 10 then goto,skipmon		; NDJFMA
;    if imon ne 2 then goto,skipmon		; Mar
;    if imon ne 0 then goto,skipmon		; Jan only
;    if imon ne 11 then goto,skipmon		; Dec only
    ifiles=file_search(dir+syear(iyear)+smon(imon)+'??.sav',count=nfile)
    if nfile eq 0L then goto,skipmon
;
; declare (nfile) arrays
;
    hlatpdf_time=fltarr(nfile,nrr)
    llatpdf_time=fltarr(nfile,nrr)
    elatedge_time=fltarr(nfile)
    coedge_time=fltarr(nfile)
    lowlat_elatedge_time=fltarr(nfile)
    lowlat_coedge_time=fltarr(nfile)
    pvelatedge_time=fltarr(nfile)
    nashedge_time=fltarr(nfile)
    sdate_time=strarr(nfile)
;
; loop over files
;
    icount=0L
    kcount=0L
    FOR n=0l,nfile-1l DO BEGIN
        result=strsplit(ifiles(n),'.',/extract)
        result2=result(1)
        result3=strsplit(result2,'_',/extract)
        sdate=result3(-1)
        print,sdate,icount
        sdate_time(n)=sdate
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
        restore,ifile
        print,'restored '+ifile
        dum=findfile(ifile)
        if dum(0) eq '' then goto,jumpstep
;
; MERRA
;
        ifile=dirm+sdate+'.nc3'
        rd_merra_nc3,ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
            u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        print,iflag,ifile
        sp2=sqrt(u2^2.+v2^2.)
;
; SD-WACCM
;
;       ifile=dirw+'f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h5.'+sdate+'.nc3'
;       rd_sdwaccm4_nc3,ifile,nc,nr,nth,alon,alat,th,$
;          pv2,p2,z2,u2,v2,q2,qdf2,mark2,sf2,h2o2,n2o2,o32,iflg
;
; choose theta surface
;
        rth=800.
;       print,th
;       read,'Enter desired theta ',rth
        index=where(abs(rth-th) eq min(abs(rth-th)))
        ith=index(0)
        sth=strcompress(long(th(ith)),/remove_all)+'K'
        pv1=transpose(reform(pv2(0:nr-1,*,ith)))
        sp1=transpose(reform(sp2(0:nr-1,*,ith)))
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
; postscript file
;
        if setplot eq 'ps' then begin
           lc=0
           xsize=nxdim/100.
           ysize=nydim/100.
           set_plot,'ps'
           device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                  /bold,/color,bits_per_pixel=8,/helvetica,filename='polar_daily_coelat+cograd_press_merrapv_'+sdate+'_'+spress+'.ps'
           !p.charsize=1.25
           !p.thick=2
           !p.charthick=5
           !y.thick=2
           !x.thick=2
        endif
;
; horizontal CO gradient
;
        co1=reform(co_grid(*,*,ilev))*1.e6
;       co1=smooth(co1,3,/edge_truncate)
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
;
; PV elat
;
        pvelat=calcelat2d(pv1,alon,alat)
;
; Geographic Latitude PDFs for Elat >60 and <80
;
        index=where(y2d gt 0. and elat ge 80.,nn)
        y2h=histogram(y2d(index),min=0,max=90,binsize=1)/float(nn)
hlatpdf_time(n,*)=y2h
highindex=where(y2h eq max(y2h))
highco=yeq(highindex(0))
print,'High Elat Lat ',highco
erase
        xmn=xorig(1)+0.075
        xmx=xorig(1)+xlen-0.075
        ymn=yorig(1)+ylen/2.
        ymx=yorig(1)+ylen
        set_viewport,xmn,xmx,ymn,ymx
        !type=2^2+2^3
        plot,yeq,y2h,color=0,/noeras,charsize=1.5,ytitle='PDF',xtitle='Latitude',xrange=[0,90],yrange=[0,max(y2h)],/nodata
index=where(y2h ne 0.)
        oplot,yeq(index),y2h(index),color=mcolor*.9,thick=5
        index=where(y2d gt 0. and elat le 60.,nn)
        y2=histogram(y2d(index),min=0,max=90,binsize=1)/float(nn)
llatpdf_time(n,*)=y2
index=where(y2 ne 0.)
        oplot,yeq(index),y2(index),color=mcolor*.3,thick=5
lowindex=where(y2 eq max(y2))
lowco=yeq(lowindex(0))
print,'Low ELat Lat ',lowco
d=lowco/highco
;
; where low elat latitudes eq 0. and high elat latitudes ne 0.
;
index=where(y2 eq 0. and y2h ne 0.)
if index(0) eq -1L then print,'NO VORTEX TODAY?'
if index(0) ne -1L then print,'vortex home? ',yeq(index)
;d_time(n)=d
;print,'d= ',d
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
print,'Candidates ',edgecandidates
;
; broadcast edges
;
print,'CO Elat Edge ',elatedge
print,'CO Elat Edge lowlat ',min(edgecandidates)
print,'PV Elat Edge ',pvelatedge
print,'Nash Edge ',nashedge
;
; save edges
;
       elatedge_time(n)=elatedge
       coedge_time(n)=coedge
       index=where(edgecandidates eq min(edgecandidates))
       lowlat_elatedge_time(n)=edgecandidates(index)
       lowlat_coedge_time(n)=coedgecandidates(index)
print,'CO Value for Elat edge: ',coedge
print,'CO Value for low lat Elat edge: ',coedgecandidates(index)
       pvelatedge_time(n)=pvelatedge
       nashedge_time(n)=nashedge
; 
; save postscript version
;
       xyouts,.7,.95,sdate,color=0,charsize=2,charthick=2,/normal
       xmn=xorig(0)
       xmx=xorig(0)+xlen
       ymn=yorig(0)
       ymx=yorig(0)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       imin=min(cobin)
       imax=max(cobin)
       nlvls=26L
       col1=1+(indgen(nlvls)/float(nlvls))*mcolor
       level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
       map_set,90,0,-90,/ortho,/contin,/grid,/noerase,color=0,charsize=1.5,title='CO '+spress	;,limit=[40.,0.,90.,360.]
       contour,co1,alon,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot
       map_set,90,0,-90,/ortho,/contin,/noerase,color=mcolor
;      contour,elat,alon,alat,levels=10+10*findgen(8),/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1,1,1,1,1,1,1,1]			; all contours
;      contour,elat,alon,alat,levels=50,/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=5
       contour,elat,alon,alat,levels=elatedge,/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=10
       contour,elat,alon,alat,levels=min(edgecandidates),/noeras,charsize=2,color=mcolor*.9,/foll,/overplot,c_labels=[1],thick=8
       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
       ybox=[0,10,10,0,0]
       x2=imin
       dx=(imax-imin)/(float(nlvls)-1)
       for j=1,nlvls-1 do begin
           xbox=[x2,x2,x2+dx,x2+dx,x2]
           polyfill,xbox,ybox,color=col1(j)
           x2=x2+dx
       endfor

       xmn=xorig(1)+0.075
       xmx=xorig(1)+xlen-0.075
       ymn=yorig(1)-0.05
       ymx=yorig(1)+ylen/2.-0.05
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       plot,yeq,cobin,color=0,/noeras,charsize=1.5,ytitle='Mean CO',xtitle='Elat',xrange=[20,90],thick=5,title='MLS'
       axis,/yax,yrange=[min(delatlevs),max(delatlevs)],/save,ytitle='dCO/elat',charsize=1.5,color=mcolor*.9
       oplot,yeq,delatlevs,color=mcolor*.9,thick=4
;      imin=min(cogradbin)
;      imax=max(cogradbin)
;      level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
;      map_set,90,0,-90,/ortho,/contin,/grid,/noerase,color=0,charsize=1.5,title='CO Gradient'      ;,limit=[40.,0.,90.,360.]
;      contour,cograd1,alon,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot
;      map_set,90,0,-90,/ortho,/contin,/noerase,color=mcolor
;      ymnb=ymn -cbaryoff
;      ymxb=ymnb+cbarydel
;      set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
;      !type=2^2+2^3+2^6
;      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1
;      ybox=[0,10,10,0,0]
;      x2=imin
;      dx=(imax-imin)/(float(nlvls)-1)
;      for j=1,nlvls-1 do begin
;          xbox=[x2,x2,x2+dx,x2+dx,x2]
;          polyfill,xbox,ybox,color=col1(j)
;          x2=x2+dx
;      endfor
;
; PV
;
       xmn=xorig(2)
       xmx=xorig(2)+xlen
       ymn=yorig(2)
       ymx=yorig(2)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       imin=min(pvbin)
       imax=max(pvbin)
       level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
       map_set,90,0,-90,/ortho,/contin,/grid,/noerase,color=0,charsize=1.5,title='PV '+sth  ;,limit=[40.,0.,90.,360.]
       contour,pv1,alon,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot
       map_set,90,0,-90,/ortho,/contin,/noerase,color=mcolor
;      contour,pvelat,alon,alat,levels=10+10*findgen(8),/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1,1,1,1,1,1,1,1]                       ; all contours
       contour,pvelat,alon,alat,levels=pvelatedge,/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=10
       contour,pvelat,alon,alat,levels=nashedge,/noeras,charsize=2,color=mcolor*.9,/foll,/overplot,c_labels=[1],thick=8
;      contour,elat,alon,alat,levels=elatedge,/noeras,charsize=2,color=mcolor,/foll,/overplot,c_labels=[1],thick=5,c_linestyle=5
       contour,sp1,alon,alat,levels=40+20*findgen(4),/noeras,charsize=2,c_color=[150,250,0],/foll,/overplot,thick=3
       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='PV'
       ybox=[0,10,10,0,0]
       x2=imin
       dx=(imax-imin)/(float(nlvls)-1)
       for j=1,nlvls-1 do begin
           xbox=[x2,x2,x2+dx,x2+dx,x2]
           polyfill,xbox,ybox,color=col1(j)
           x2=x2+dx
       endfor
;
; PV Gradient
;
       xmn=xorig(3)+0.075
       xmx=xorig(3)+xlen-0.075
       ymn=yorig(3)+0.05
       ymx=yorig(3)+ylen-0.05
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
       plot,yeq,pvbin,color=0,/noeras,charsize=1.5,ytitle='Mean PV',xtitle='Elat',xrange=[20,90],thick=5,title='MERRA'
       axis,/yax,yrange=[min(spbin),max(spbin)],/save,ytitle='Speed/elat',charsize=1.5,color=mcolor*.3
       oplot,yeq,spbin,color=mcolor*.3,thick=4
       axis,/yax,yrange=[min(dpvelatlevs),max(dpvelatlevs)],/save,ytitle='dPV/elat',charsize=1.5,color=mcolor*.9
       oplot,yeq,dpvelatlevs,color=mcolor*.9,thick=4

;      imin=min(pvgradbin)
;      imax=max(pvgradbin)
;      level=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
;      map_set,90,0,-90,/ortho,/contin,/grid,/noerase,color=0,charsize=1.5,title='PV Gradient'  ;,limit=[40.,0.,90.,360.]
;      contour,pvgrad1,alon,alat,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,/overplot
;      map_set,90,0,-90,/ortho,/contin,/noerase,color=mcolor
;      ymnb=ymn -cbaryoff
;      ymxb=ymnb+cbarydel
;      set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
;      !type=2^2+2^3+2^6
;      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='PV Gradient'
;      ybox=[0,10,10,0,0]
;      x2=imin
;      dx=(imax-imin)/(float(nlvls)-1)
;      for j=1,nlvls-1 do begin
;          xbox=[x2,x2,x2+dx,x2+dx,x2]
;          polyfill,xbox,ybox,color=col1(j)
;          x2=x2+dx
;      endfor

; Close PostScript file and return control to X-windows
       if setplot ne 'ps' then stop	;wait,1
       if setplot eq 'ps' then begin
          device, /close
          spawn,'convert -trim polar_daily_coelat+cograd_press_merrapv_'+sdate+'_'+spress+'.ps -rotate -90 '+$
                              'polar_daily_coelat+cograd_press_merrapv_'+sdate+'_'+spress+'.jpg'
          spawn,'rm -f polar_daily_coelat+cograd_press_merrapv_'+sdate+'_'+spress+'.ps'
       endif
jumplev:
;endfor
icount=icount+1L
jumpstep:
endfor	; loop over files
;
; save monthly file
;
save,filename='daily_mls_coelatedge+merra_pvelatedge_'+syear(iyear)+smon(imon)+'_'+spress+'.sav',elatedge_time,coedge_time,lowlat_elatedge_time,lowlat_coedge_time,pvelatedge_time,nashedge_time,sdate_time,spress,sth,hlatpdf_time,llatpdf_time

skipmon:
endfor	; loop over months
endfor	; loop over years
end
