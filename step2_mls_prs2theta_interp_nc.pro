;
; read MLS gridded IDL save files.
; Input data:  IDL> restore,'/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_YYYYMMDD.sav
; Input data:  IDL> restore,'/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_U_V_v3.3_YYYYMMDD.sav
; 
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
; U               FLOAT     = Array[144, 96, 55]
; V               FLOAT     = Array[144, 96, 55]
; 
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
sver = 'v3.3'
loadct,39
device,decompose=0
mcolor=byte(!p.color)
nlvls=30L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
thlev=[10000.,9000.,8000.,7000.,6000.,5000.,4800.,4600.,4400.,4200.,4000.,3800.,3600.,$
       3400.,3200.,3000.,2800.,2600.,2400.,2200.,2000.,$
       1800.,1600.,1400.,1200.,1000., 900., 800., 700.,$
        600., 550., 500., 450., 400., 350.]
nth=n_elements(thlev)
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
nru=72L & ncu=144L
dirw='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_'
ifiles=file_search(dirw+'ALL_'+sver+'_20??????.sav')
ufiles=file_search(dirw+'U_V_'+sver+'_20??????.sav')
nfile=n_elements(ifiles)
nfile2=n_elements(ufiles)
if nfile ne nfile2 then stop,'Check number of all and uv files'
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN

result=strsplit(ifiles(n),'_',/extract)
result2=strsplit(result(-1),'.',/extract)
sdate=result2(0)
print,sdate
dum=findfile(dirw+'theta_'+sdate+'.nc.sav')
if dum(0) ne '' then goto,jumpfile
;if long(strmid(sdate,4,2)) ne 12L then goto,jumpfile	; December
;dum=findfile(dirw+'theta_'+sdate+'.nc3')
;if dum(0) ne '' then goto,jumpfile
restore,dirw+'ALL_'+sver+'_'+sdate+'.sav'
print,'reading '+dirw+'ALL_'+sver+'_'+sdate+'.sav'
restore,dirw+'U_V_'+sver+'_'+sdate+'.sav'
print,'reading '+dirw+'U_V_'+sver+'_'+sdate+'.sav'
alon=LON
alat=LAT
pressure=pmls2
press37=pmls
nc=n_elements(alon)
nr=n_elements(alat)
nl=n_elements(pressure)
nl37=n_elements(press37)
tl=tp_grid
vl=smooth(v,3,/edge_truncate,/Nan)
ul=smooth(u,3,/edge_truncate,/Nan)
zl=gp_grid
h2ol=h2o_grid
;
; interpoloate co_grid to col with same pressure as temperature
;IDL> print,pmls
;      1000.00      681.292      464.159      316.228      215.443      146.780      100.000      68.1292      46.4159      31.6228      21.5443      14.6780
;      10.0000      6.81292      4.64159      3.16228      2.15443      1.46780      1.00000     0.681292     0.464159     0.316228     0.215443     0.146780
;     0.100000    0.0464159    0.0215443    0.0100000   0.00464159   0.00215443   0.00100000  0.000464159  0.000215443  0.000100000  4.64159e-05  2.15443e-05
;  1.00000e-05
;IDL> print,pmls2
;      1000.00      825.404      681.292      562.341      464.159      383.119      316.228      261.016      215.443      177.828      146.780      121.153
;      100.000      82.5404      68.1292      56.2341      46.4159      38.3119      31.6228      26.1016      21.5443      17.7828      14.6780      12.1153
;      10.0000      8.25404      6.81292      5.62341      4.64159      3.83119      3.16228      2.61016      2.15443      1.77828      1.46780      1.21153
;      1.00000     0.681292     0.464159     0.316228     0.215443     0.146780     0.100000    0.0464159    0.0215443    0.0100000   0.00464159   0.00215443
;   0.00100000  0.000464159  0.000215443  0.000100000  4.64159e-05  2.15443e-05  1.00000e-05
;
col=0.*tl
for k=0L,nl-1L do begin
    zp=pressure(k)
    for kk=0L,nl37-2L do begin
        if pmls(kk) ge zp and pmls(kk+1L) le zp then begin
           zscale=(alog(zp)-alog(pmls(kk)))/(alog(pmls(kk+1L))-alog(pmls(kk)))
           col(*,*,k)=co_grid(*,*,kk)-zscale*(co_grid(*,*,kk)-co_grid(*,*,kk+1))
;print,pmls(kk+1),zp,pmls(kk),zscale
;print,co_grid(10,10,kk+1),col(10,10,k),co_grid(10,10,kk)
;if finite(co_grid(10,10,kk)) eq 1 then stop
        endif
    endfor
    jumplev:
endfor
;
; flip pressure data to run from top to bottom
;
pressure=reverse(pressure)
for i=0L,nc-1L do begin
for j=0L,nr-1L do begin
    tmp=reform(tl(i,j,*))
    tl(i,j,*)=reverse(tmp)
    tmp=reform(vl(i,j,*))
    vl(i,j,*)=reverse(tmp)
    tmp=reform(ul(i,j,*))
    ul(i,j,*)=reverse(tmp)
    tmp=reform(zl(i,j,*))
    zl(i,j,*)=reverse(tmp)
    tmp=reform(h2ol(i,j,*))
    h2ol(i,j,*)=reverse(tmp)
    tmp=reform(col(i,j,*))
    col(i,j,*)=reverse(tmp)
endfor
endfor
;
; initialize 3d pressure array
;
prl=0.*tl
for i=0L,nc-1L do $
for j=0L,nr-1L do $
    prl(i,j,*)=pressure
;
; save memory
;
tp_grid=0 & u=0 & v=0 & gp_grid=0 & h2o_grid=0 & n2o_grid=0 & o3_grid=0 & co_grid=0
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
H2OGRD=fltarr(nr,nc,nth)
COGRD=fltarr(nr,nc,nth)
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
      H2OGRD(LAT,I,M)=0.
      COGRD(LAT,I,M)=0.

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
      H2OGRD(LAT,I,M)=H2OL(I,LAT,K)+SCALE*(H2OL(I,LAT,LM1)-H2OL(I,LAT,K))
      COGRD(LAT,I,M)=COL(I,LAT,K)+SCALE*(COL(I,LAT,LM1)-COL(I,LAT,K))
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
rlev=3000.
;print,thlev
;read,'Enter theta surface ',rlev
index=where(thlev eq rlev)
ilev=index(0)
slev=string(rlev)
pp=transpose(qdfgrd(*,*,ilev))
;pp=transpose(ipvgrd(*,*,ilev))
imin=-1000.
imax=1000.
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
!type=2^2+2^3
erase
MAP_SET,90,0,-90,/stereo,/noeras,/grid,/contin,title=sdate+' '+slev+' K',charsize=2.0
contour,pp,alon,alat,levels=level,/cell_fill,c_color=col1,/noeras,/overplot
;contour,pp,alon,alat,levels=level,/follow,c_color=0,/noeras,/overplot
;contour,pp,alon,alat,levels=[0.],/follow,c_color=0,thick=3,/noeras,/overplot
;pp=transpose(cogrd(*,*,ilev))
;level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
index=where(level lt 0.)
if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=mcolor,c_linestyle=5,/noeras,/overplot
index=where(level gt 0.)
if index(0) ne -1L then contour,pp,alon,alat,levels=level(index),/follow,c_color=0,/noeras,/overplot
pp=transpose(zgrd(*,*,ilev))
level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
contour,pp,alon,alat,levels=level,/follow,c_color=0,/noeras,/overplot,thick=3
;
; write theta file
;
ofile=dirw+'theta_'+sdate+'.nc.sav'
print,'writing ',ofile
save,file=ofile,nr,nc,nth,alon,alat,thlev,ipvgrd,pgrd,ugrd,vgrd,zgrd,qdfgrd,h2ogrd,cogrd

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
