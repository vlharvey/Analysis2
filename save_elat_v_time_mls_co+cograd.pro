;
; save daily mean CO and CO gradient vs Elat
; CO horizontal gradient - dCO/dx + dCO/dy
;
@calcelat2d

re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
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
xorig=[0.55,0.55,0.1,0.1]
yorig=[0.05,0.05,0.55,0.55]
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
dir='/Volumes/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_'

for iyear=0L,nyear-1L do begin
ifiles=file_search(dir+syear(iyear)+'????.sav',count=nfile)
if nfile eq 0L then goto,skipmon
;
; declare time based arrays of CO, CO gradient, retain edges
;
ytco=fltarr(nfile,nrr)
ytdco=fltarr(nfile,nrr)
ytprod=fltarr(nfile,nrr)
ytcograd=fltarr(nfile,nrr)
ytcoedge=fltarr(nfile)
ytelatedge=fltarr(nfile)
sdate_all=strarr(nfile)
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
    sdate_all(n)=sdate
    print,sdate,icount
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
;
; first day
;
;if kcount eq 0L then begin
   rpress=1.
;  print,pmls
;  read,'Enter desired pressure ',rpress
   index=where(abs(rpress-pmls) eq min(abs(rpress-pmls)))
   ilev=index(0)
   spress=strcompress(pmls(ilev),/remove_all)+'hPa'
   nc=144			;n_elements(lon)
   nr=96			;n_elements(lat)
   nl=n_elements(pmls)
   alat=-90+1.89474*findgen(nr)	;lat	(files in 2004 and others? have latitude(nprof) arrays)
   alat(nr-1)=90.
   alon=2.5*findgen(nc)		;lon
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
   kcount=1
;endif
;
; normalized CO gradient
;
cograd2=0.*co_grid
for k=0,nl-1 do begin
co1=co_grid(*,*,k)
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
        cograd1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)/abs(co1(i,j))	;/comax
;
; allowing negative on vortex interior is useful to isolate edge but if vortex is offset from the pole dCO/dy<0 should be >0
;
;       if (dqdy le 0.0) then cograd1(i,j) = -sqrt(dqdx*dqdx+dqdy*dqdy)$
;                                            *abs(co1(i,j))/comax
;
; without normalization
;
;       cograd1(i,j) = sqrt(dqdx*dqdx+dqdy*dqdy)       ;/abs(co1(i,j))
;       if (dqdy le 0.0) then cograd1(i,j) = -1.0*cograd1(i,j)

    endfor
endfor
cograd2(*,*,k)=smooth(cograd1*1.e6,3,/edge_truncate)
endfor
;
; loop over levels
;
;for ith=20L,20L do begin	;0,nl-1L do begin
;
; extract level
;
    n2o=n2o_grid(*,*,ilev)*1.e9
    co=co_grid(*,*,ilev)*1.e6
;
; area of CO "Nash" vortex
;
    index=where(lat lt 0.)
    cosave=co
    cosave(index)=-1.*cosave(index)
    elat=calcelat2d(cosave,alon,alat)
;
    cograd=cograd2(*,*,ilev)
index=where(elat le 10.)	; zero tropical gradients that got enhanced by divide by co
cograd(index)=0.
cograd=smooth(cograd,5,/edge_truncate)
    index=where(abs(pmls2-pmls(ilev)) eq min(abs(pmls2-pmls(ilev))))
    gp=gp_grid(*,*,index(0))/1000.
;
; set CO levels for hemisphere
;
nlvls=26L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
nhindex=where(y2d gt 0.,nn)
imin=min(co(nhindex))
imax=max(co(nhindex))
colevel=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
if max(colevel) eq 0. then goto,jumplev			; bad data
if finite(max(colevel)) eq 0 then goto,jumplev
;
; elat as a function of CO
;
elatlevs=0.*colevel
for i=1L,nlvls-1L do begin
    index=where(y2d gt 0. and co gt colevel(i-1) and co le colevel(i))
    if index(0) ne -1L then elatlevs(i)=mean(elat(index))
endfor
;
; mean CO within 1 degree spaced elat
; mean CO gradient
;
cobins=fltarr(nrr)
cogradbins=fltarr(nrr)
for i=0L,nrr-1L do begin
    ip1=i+1
    if i eq nrr-1 then ip1=nrr-1
    im1=i-1
    if i eq 0L then im1=0
    index=where(y2d gt 0. and elat ge yeq(im1) and elat lt yeq(ip1))
    if index(0) ne -1L then cobins(i)=mean(co(index))
    if index(0) ne -1L then cogradbins(i)=mean(cograd(index))
endfor
;
; retain today co values and mean gradient values per elat bin at this lev
;
ytco(n,*)=cobins
ytcograd(n,*)=cogradbins

delatlevs=smooth(deriv(yeq,cobins),7,/edge_truncate)
ytdco(n,*)=delatlevs		; save rate of change of CO with Elat

d2elatlevs=deriv(delatlevs)
;
; product of 1) mean COgrad as a function of elat and 2) dCO/delat
;
prod=delatlevs*cogradbins
ytprod(n,*)=prod
;
; extract node from d2Elat/dCO2
;
elatsave=yeq
d2elatsave=d2elatlevs
;
; absolute maximum in the first derivative - as Nash
;
;colevsave=colevel(1:-1)
colevsave=cobins
index=where(delatlevs eq max(delatlevs))
elatedge=elatsave(index(0))
coedge=colevsave(index(0))
; 
; absolute maximum in the product
;
index=where(prod eq max(prod))
elatedge2=elatsave(index(0))
coedge2=colevsave(index(0))
ytcoedge(n)=coedge2
ytelatedge(n)=elatedge2
;
; most equatorward local maximum in the first derivative
; IS THIS NECESSARY?
;
flag=0.*delatlevs
for i=1L,n_elements(delatlevs)-2L do begin
    if delatlevs(i) gt delatlevs(i-1) and delatlevs(i) gt delatlevs(i+1) then flag(i)=1.
endfor
nhindex=where(y2d gt 0.,nn)
cohemmean=mean(co(nhindex))
cohemsig=stdev(co(nhindex))
index=where(flag eq 1 and cobins gt cohemmean-0.5*cohemsig)
elatedge=yeq(index(0))
coedge=colevsave(index(0))
;
; mean gradient within CO bins
;
cogradbin=0.*colevel
for i=1L,nlvls-1L do begin
    index=where(co gt colevel(i-1) and co le colevel(i),nn)
    if index(0) ne -1L then cogradbin(i)=mean(cograd(index))
endfor
cogradbin=smooth(cogradbin,3,/edge_truncate)	; want edge influence to reduce highest-most CO bins
delatgradlevs=smooth(deriv(yeq,cogradbins),7,/edge_truncate)
d2elatgradlevs=deriv(delatgradlevs)     ;,cobins)
;
elatsave=yeq
d2elatgradsave=d2elatgradlevs
;
; superimpose all local maxima - bimodal distributions usually want lower CO value
;
cogradnodes=0./0.
for i=1L,nlvls-2L do begin
    if cogradbin(i-1) lt cogradbin(i) and cogradbin(i+1) lt cogradbin(i) then begin
       cogradnodes=[cogradnodes,colevel(i)]
    endif
endfor
index=where(finite(cogradnodes) eq 1.)
cogradnodes=cogradnodes(index)

jumplev:
;endfor
icount=icount+1L
jumpstep:
endfor	; loop over files
;
; save yearly file at 1 hPa
;
save,filename='elat_v_time_mls_co+cograd_'+syear(iyear)+'_'+spress+'.sav',ytco,ytdco,ytprod,ytcograd,ytcoedge,ytelatedge,yeq,sdate_all,nfile,spress
;
; postscript file
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='elat_v_time_mls_co+cograd_'+syear(iyear)+'_'+spress+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
; 
; plot
;
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=26L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
imin=min(ytco(nhindex))
imax=max(ytco(nhindex))
colevel=imin+((imax-imin)/float(nlvls-1))*findgen(nlvls)
contour,ytco,findgen(nfile),yeq,levels=colevel,xrange=[0,nfile],yrange=[0,90],/noeras,color=0,c_color=col1




;
; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim elat_v_time_mls_co+cograd_'+syear(iyear)+'_'+spress+'.ps -rotate -90 '+$
                            'elat_v_time_mls_co+cograd_'+syear(iyear)+'_'+spress+'.jpg'
;       spawn,'rm -f elat_v_time_mls_co+cograd_'+syear(iyear)+'_'+spress+'.ps'
     endif
;
skipmon:
endfor	; loop over years
end
