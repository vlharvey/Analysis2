;
; read MERRA IDL save files.
; Input data:  IDL> restore,'/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_YYYYMMDD.sav
; Input data:  IDL> restore,'/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_rad_YYYYMMDD.sav
; 
; LATITUDE_WACCM  FLOAT     = Array[96]
; LONGITUDE_WACCM FLOAT     = Array[144]
; PRESSURE        FLOAT     = Array[41]
; PSGRD           FLOAT     = Array[144, 96]
; QVGRD           FLOAT     = Array[144, 96, 41]
; TGRD            FLOAT     = Array[144, 96, 41]
; UGRD            FLOAT     = Array[144, 96, 41]
; VGRD            FLOAT     = Array[144, 96, 41]
; ZGRD            FLOAT     = Array[144, 96, 41]
; add QGRD            FLOAT     = Array[144, 96, 41]

; ALAT            FLOAT     = Array[96]
; ALON            FLOAT     = Array[144]
; PGRD            FLOAT     = Array[144, 96, 72]
; TGRD            FLOAT     = Array[144, 96, 72]
; UGRD            FLOAT     = Array[144, 96, 72]
; VGRD            FLOAT     = Array[144, 96, 72]
;		
; calculate PV, QDF and 
; interpolate to isentropic surfaces.  Output daily theta .nc files.
; next step: calculation of stream function
;
@compvort

loadct,38
device,decompose=0
mcolor=byte(!p.color)
nlvls=30L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
thlev=[5000.,4800.,4600.,4400.,4200.,4000.,3800.,3600.,$
       3400.,3200.,3000.,2800.,2600.,2400.,2200.,2000.,$
       1800.,1600.,1400.,1200.,1000., 900., 800., 700.,$
        600., 550., 500., 450., 400., 350.]
nth=n_elements(thlev)
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
nru=72L & ncu=144L
ulat=-88.75+2.5*findgen(nru)
ulon=2.5*findgen(ncu)
dirw='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_'
ifiles=file_search(dirw+'press_2015????.sav')
nfile=n_elements(ifiles)
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN

result=strsplit(ifiles(n),'_',/extract)
result2=strsplit(result(3),'.',/extract)
sdate=result2(0)
print,sdate
dum=findfile(dirw+'theta_'+sdate+'.nc3')
if dum(0) ne '' then goto,jumpfile
restore,dirw+'press_'+sdate+'.sav'
print,'reading '+dirw+'press_'+sdate+'.sav'
restore,dirw+'rad_'+sdate+'.sav'
alon=LONGITUDE_WACCM
alat=LATITUDE_WACCM
nc=n_elements(alon)
nr=n_elements(alat)
nl=n_elements(pressure)
tl=tgrd
vl=vgrd
ul=ugrd
zl=zgrd
qvl=qvgrd
;QNEW
;
; remove tropospheric levels for diabatic heating rate for which there is no Nv data
;
index=where(levels_rad le max(pressure))
levels_rad=levels_rad(index)
qnew=reform(qnew(*,*,index))
;
; interpoloate qnew to ql with same pressure levels as Nv data (fill above 0.1 with NaN)
;print,levels_rad
;      0.10000000      0.30000001      0.40000001      0.50000000      0.69999999       1.0000000       2.0000000       3.0000000       4.0000000       5.0000000       7.0000000       10.000000       20.000000
;       30.000000       40.000000       50.000000       70.000000       100.00000       150.00000       200.00000       250.00000       300.00000       350.00000       400.00000       450.00000       500.00000
;IDL> print,pressure
;    0.0150000    0.0200000    0.0300000    0.0400000    0.0500000    0.0600000    0.0700000    0.0800000     
;     0.100000     0.150000     0.200000     0.300000     0.400000     0.500000     0.600000     0.700000
;     0.800000      1.00000      1.50000      2.00000      3.00000      4.00000      5.00000      6.00000      7.00000      8.00000      10.0000      15.0000      20.0000      30.0000      40.0000      50.0000
;      60.0000      70.0000      80.0000      100.000      150.000      200.000      300.000      400.000      500.000
;
ql=0.*tl
for k=0L,nl-1L do begin
    zp=pressure(k)
    if zp lt min(levels_rad) then ql(*,*,k)=0./0.
    if zp lt min(levels_rad) then goto,jumplev
    for kk=0L,n_elements(levels_rad)-2L do begin
        if levels_rad(kk) le zp and levels_rad(kk+1L) ge zp then begin
           zscale=(alog(zp)-alog(levels_rad(kk)))/(alog(levels_rad(kk+1L))-alog(levels_rad(kk)))
           ql(*,*,k)=qnew(*,*,kk)-zscale*(qnew(*,*,kk)-qnew(*,*,kk+1))
;print,levels_rad(kk+1),zp,levels_rad(kk),zscale
;print,qnew(10,10,kk+1),ql(10,10,k),qnew(10,10,kk)
;stop
        endif
    endfor
    jumplev:
endfor
;
; initialize 3d pressure array
;
prl=0.*tgrd
for i=0L,nc-1L do $
for j=0L,nr-1L do $
    prl(i,j,*)=pressure
;
; save memory
;
tgrd=0 & ugrd=0 & vgrd=0 & pgrd=0 & zgrd=0 & qvgrd=0
;
; calculate theta, absolute vorticity, potential vorticity
;
thl=0.*tl
for L=0L,NL-1L do $
    THL(*,*,L)=TL(*,*,L)*(1000./PRL(*,*,L))^.286
eta=0.*tl
compvort,ul,vl,eta,alon,alat,nc,nr
;
; SKIP AS LOWEST LEVEL IS 500 HPA
; initialize surface arrays
;
;tsfc=reform(tl(*,*,nl-1L))
;thsfc=reform(thl(*,*,nl-1L))
;usfc=reform(ul(*,*,nl-1L))
;vsfc=reform(vl(*,*,nl-1L))
;qdfsfc=0.*tsfc
;pvsfc=0.*tsfc
;
; initialize isentropic arrays
;
UGRD=fltarr(nr,nc,nth)
VGRD=fltarr(nr,nc,nth)
IPVGRD=fltarr(nr,nc,nth)
PGRD=fltarr(nr,nc,nth)
ZGRD=fltarr(nr,nc,nth)
QDFGRD=fltarr(nr,nc,nth)
QVGRD=fltarr(nr,nc,nth)
QGRD=fltarr(nr,nc,nth)
;
; LOOP OVER LATITUDES
;
pv=0.*tl
qdf=0.*tl
for LAT=0L,NR-1L do begin
    JP1=LAT-1
    JM1=LAT+1
    JP2=LAT-2
    JM2=LAT+2
    IF LAT EQ 0 THEN begin
       JP1=0
       JP2=0
    ENDIF
    IF LAT EQ NR-1L THEN begin
       JM1=NR-1L
       JM2=NR-1L
    ENDIF
    IF LAT EQ NR-2 then JM2=NR-1
    IF LAT EQ 1 then JP2=0
    DY1=RADEA*(ALAT(JP1)-ALAT(JM1))*DTR
    DY2=RADEA*(ALAT(JP2)-ALAT(JM2))*DTR
    DX1=RADEA*COS(ALAT(LAT)*DTR)*PI2/(.5*NC)
    DX2=RADEA*COS(ALAT(LAT)*DTR)*PI2/(.25*NC)
;
; LOOP OVER LONGITUDES
;
    for I=0,NC-1L do begin
        IP1=I+1
        IM1=I-1
        IP2=I+2
        IM2=I-2
        IF I EQ 0 THEN begin
           IM1=NC-1
           IM2=NC-2
        ENDIF
        IF I EQ NC-1 THEN begin
           IP1=0
           IP2=1
        ENDIF
        IF I EQ 1 then IM2=NC-1
        IF I EQ NC-2 then IP2=0
;
; COMPUTE ISENTROPIC POTENTIAL VORTICITY ON PRESSURE SURFACE
;
        for K=0,NL-1L do begin
            LM1=K+1
            LP1=K-1
            IF K EQ NL-1L then LM1=NL-1L
            IF K EQ 0L then LP1=0L
            DTHDP=(THL(I,LAT,LP1)-THL(I,LAT,LM1))/(PRL(I,LAT,LP1)-PRL(I,LAT,LM1))
            DUDP=(ul(I,LAT,LP1)-ul(I,LAT,LM1))/(PRL(I,LAT,LP1)-PRL(I,LAT,LM1))
            DVDP=(vl(I,LAT,LP1)-vl(I,LAT,LM1))/(PRL(I,LAT,LP1)-PRL(I,LAT,LM1))
            DTHDX=(4./3.)*(THL(IP1,LAT,K)-THL(IM1,LAT,K))/DX1 - $
                  (1./3.)*(THL(IP2,LAT,K)-THL(IM2,LAT,K))/DX2
            IF LAT LE 1 OR LAT GE NR-2 THEN begin
               DTHDY=(THL(I,JP1,K)-THL(I,JM1,K))/DY1
            endif
            IF LAT gt 1 and LAT lt NR-2 THEN begin
               DTHDY=(4./3.)*(THL(I,JP1,K)-THL(I,JM1,K))/DY1 - $
                     (1./3.)*(THL(I,JP2,K)-THL(I,JM2,K))/DY2
            ENDIF
            IF DTHDP GE 0. THEN begin
               PV(I,LAT,K)=1.E12
            endif
            IF DTHDP lt 0. THEN begin
               PV(I,LAT,K)=eta(I,LAT,K)-(DUDP*DTHDY-DVDP*DTHDX)/DTHDP
               PV(I,LAT,K)=-9.8*DTHDP*PV(I,LAT,K)/100.
            ENDIF

; normalized by RADEA. The signed sqrt of Q is taken.
            arg1 = (ul(IP1,LAT,K)-ul(IM1,LAT,K))/DX1 $
                  - vl(I,LAT,K)*TAN(ALAT(LAT)*DTR)/RADEA
            arg2 = (vl(IP1,LAT,K)-vl(IM1,LAT,K))/DX1 $
                  + ul(I,LAT,K)*TAN(ALAT(LAT)*DTR)/RADEA
            DVDY = (vl(I,JP1,K)-vl(I,JM1,K))/DY1
            DUDY = (ul(I,JP1,K)-ul(I,JM1,K))/DY1
            if abs(arg1) gt 1.e12 or abs(arg2) gt 1.e12 then QDF(I,LAT,K) = 1.e12
            if abs(arg1) lt 1.e12 and abs(arg2) lt 1.e12 then begin
               qtemp=(0.5*(arg1*arg1+DVDY*DVDY)+arg2*DUDY)*RADEA*RADEA
               if qtemp ge 0.0 then QDF(I,LAT,K) = sqrt(qtemp)
               if qtemp lt 0.0 then QDF(I,LAT,K) = -sqrt(-qtemp)
            endif
; end calculation of QDF deformation diagnostic
        endfor
;       qdfsfc(i,lat)=qdf(i,lat,nl-1)
;       pvsfc(i,lat)=pv(i,lat,nl-1)

; LINEARLY INTERPOLATE IN THETA TO ISENTROPIC LEVEL
      for m=0,nth-1L do begin
      UGRD(LAT,I,M)=0.
      VGRD(LAT,I,M)=0.
      IPVGRD(LAT,I,M)=1.E12
      PGRD(LAT,I,M)=0.
      ZGRD(LAT,I,M)=0.
      QDFGRD(LAT,I,M)=0.
      QVGRD(LAT,I,M)=0.
      QGRD(LAT,I,M)=0.

      for K=1,NL-1L do begin
      LM1=K-1
      LP1=K+1
      IF K EQ NL-1L then LP1=NL-1L

      IF THLEV(M) LE THL(I,LAT,LM1) AND THLEV(M) GT THL(I,LAT,K) THEN begin
      PLM1=PRL(I,LAT,LM1)
      PL=PRL(I,LAT,K)
      SCALE=(THLEV(M)-THL(I,LAT,K))/(THL(I,LAT,LM1)-THL(I,LAT,K))
      UGRD(LAT,I,m)=UL(I,LAT,K)+SCALE*(UL(I,LAT,LM1)-UL(I,LAT,K))
      VGRD(LAT,I,m)=VL(I,LAT,K)+SCALE*(VL(I,LAT,LM1)-VL(I,LAT,K))
      IPVGRD(LAT,I,m)=PV(I,LAT,K)+SCALE*(PV(I,LAT,LM1)-PV(I,LAT,K))
      PGRD(LAT,I,m)=PL^.286 + SCALE*(PLM1^.286-PL^.286)
      PGRD(LAT,I,m)=PGRD(LAT,I,m)^(1./.286)
      ZGRD(LAT,I,M)=ZL(I,LAT,K)+SCALE*(ZL(I,LAT,LM1)-ZL(I,LAT,K))
      QDFGRD(LAT,I,M)=QDF(I,LAT,K)+SCALE*(QDF(I,LAT,LM1)-QDF(I,LAT,K))
      QVGRD(LAT,I,M)=QVL(I,LAT,K)+SCALE*(QVL(I,LAT,LM1)-QVL(I,LAT,K))
      QGRD(LAT,I,M)=QL(I,LAT,K)+SCALE*(QL(I,LAT,LM1)-QL(I,LAT,K))
      ENDIF
      endfor

; if desired theta surface intersects the ground then
; isentropic quantities are equal to the sfc values (at psfc)
;     if thlev(m) le thsfc(i,lat) then begin
;     PGRD(LAT,I,M)=psfc(i,lat)
;     UGRD(LAT,I,M)=usfc(i,lat)
;     VGRD(LAT,I,M)=vsfc(i,lat)
;     IPVGRD(LAT,I,M)=pvsfc(i,lat)
;     QDFGRD(LAT,I,M)=qdfsfc(i,lat)
;     endif

      endfor

    endfor
endfor

; average QDF for polar rows
for M = 0, NTH-1 do begin
    xlst=0.0
    frst=0.0
    for I = 0, nc-1 do begin
       xlst = xlst + QDFGRD(nr-2,i,m)/float(nc)
       frst = frst + QDFGRD(1,i,m)/float(nc)
    endfor
    for I = 0, nc-1 do begin
       QDFGRD(nr-1,i,m) = xlst
       QDFGRD(0,i,m) = frst
    endfor
endfor
;
; check
; 
rlev=2000.
;print,thlev
;read,'Enter theta surface ',rlev
index=where(thlev eq rlev)
ilev=index(0)
slev=string(rlev)
pp=transpose(qdfgrd(*,*,ilev))
pp=transpose(zgrd(*,*,ilev))
;pp=transpose(ipvgrd(*,*,ilev))
level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
!type=2^2+2^3
erase
contour,pp,alon,alat,levels=level,/cell_fill,c_color=col1,/noeras,xrange=[0.,360.],$
        yrange=[-90.,90.],title=sdate+'  '+slev+' K'
contour,pp,alon,alat,levels=level,/follow,c_color=0,/noeras,/overplot
contour,pp,alon,alat,levels=[0.],/follow,c_color=0,thick=3,/noeras,/overplot
pp=transpose(qgrd(*,*,ilev))
level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
index=where(level lt 0.)
if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=mcolor,c_linestyle=5,/noeras,/overplot
index=where(level gt 0.)
;if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=0,/noeras,/overplot
;
; write theta file
;
ofile=dirw+'theta_'+sdate+'.nc.sav'
print,'writing ',ofile
save,file=ofile,nr,nc,nth,alon,alat,thlev,ipvgrd,pgrd,ugrd,vgrd,zgrd,qdfgrd,qvgrd,qgrd

;nocid = ncdf_create(ofile,/CLOBBER)
;latdimid=ncdf_dimdef(nocid, 'number_of_latitudes' , nr)
;londimid=ncdf_dimdef(nocid, 'number_of_longitudes', nc)
;levdimid=ncdf_dimdef(nocid, 'number_of_levels'    , nth)
;lonsid = ncdf_vardef(nocid, 'longitudes',  londimid)
;latsid = ncdf_vardef(nocid, 'latitudes' ,  latdimid)
;levsid = ncdf_vardef(nocid, 'th_levels' ,  levdimid)
;ipvid  = ncdf_vardef(nocid, 'ipv'       , [latdimid,londimid,levdimid])
;prsid  = ncdf_vardef(nocid, 'press'     , [latdimid,londimid,levdimid])
;uuuid  = ncdf_vardef(nocid, 'u_wind'    , [latdimid,londimid,levdimid])
;vvvid  = ncdf_vardef(nocid, 'v_wind'    , [latdimid,londimid,levdimid])
;qdfid  = ncdf_vardef(nocid, 'qdf'       , [latdimid,londimid,levdimid])
;qdfid  = ncdf_vardef(nocid, 'qv'        , [latdimid,londimid,levdimid])
;zid  = ncdf_vardef(nocid, 'gph'       , [latdimid,londimid,levdimid])
;ncdf_control,nocid,/ENDEF
;ncdf_varput, nocid, lonsid, alon  , COUNT=[nc]
;ncdf_varput, nocid, latsid, alat  , COUNT=[nr]
;ncdf_varput, nocid, levsid, thlev , COUNT=[nth]
;ncdf_varput, nocid, ipvid , ipvgrd, COUNT=[nr,nc,nth]
;ncdf_varput, nocid, prsid , pgrd  , COUNT=[nr,nc,nth]
;ncdf_varput, nocid, uuuid , ugrd  , COUNT=[nr,nc,nth]
;ncdf_varput, nocid, vvvid , vgrd  , COUNT=[nr,nc,nth]
;ncdf_varput, nocid, qdfid , qdfgrd, COUNT=[nr,nc,nth]
;ncdf_varput, nocid, qdfid , qvgrd , COUNT=[nr,nc,nth]
;ncdf_varput, nocid, zid , zgrd , COUNT=[nr,nc,nth]
;ncdf_close,nocid
jumpfile:
ENDFOR		; LOOP OVER files
end
