;
; read WACCM U, V, T monthly netCDF file. calculate PV, QDF and 
; interpolate to isentropic surfaces.  Output daily theta .nc files.
;
; interpolate to 72 latitudes as MetO -88.75 to 88.75 by 2.5
; and 144 longitudes 0 to 357.5 by 2.5
; for calculation of stream function
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
ncfile='/aura3/randall/waccm/U_V_T_S1.nc'
ncid=ncdf_open(ncfile)
result0=ncdf_inquire(ncid)
for idim=0,result0.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'lon' then nc=dim
    if name eq 'lat' then nr=dim
    if name eq 'lev' then nl=dim
    if name eq 'time' then nt=dim
;   print,'read ',name,' dimension ',dim
endfor
ncfile2='/aura3/randall/waccm/O3_CH4_NOY_S1.nc'
ncid2=ncdf_open(ncfile2)
result2=ncdf_inquire(ncid2)
;
; loop over timesteps
;
FOR ITIME=0l,NT-1l DO BEGIN
;
; loop over variables
;
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    if result.name ne 'T' and result.name ne 'U' and result.name ne 'V' then $
       ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'P0' then p0=data
    if result.name eq 'lat' then alat=data
    if result.name eq 'lon' then alon=data
    if result.name eq 'lev' then lev=data
    if result.name eq 'ilev' then ilev=data
    if result.name eq 'time' then time=data
    if result.name eq 'hyai' then hyai=data
    if result.name eq 'hybi' then hybi=data
    if result.name eq 'hyam' then hyam=data
    if result.name eq 'hybm' then hybm=data
    if result.name eq 'date' then date=data
    if result.name eq 'PS' then psfc=data/100.
    if result.name eq 'T' or result.name eq 'U' or result.name eq 'V' then $
       ncdf_varget,ncid,ncdf_varid(ncid,result.name),data, OFFSET=[0,0,0,itime], COUNT=[nc,nr,nl,1]
    if result.name eq 'T' then tl=data
    if result.name eq 'U' then ul=data
    if result.name eq 'V' then vl=data
;   print,ivar,result.name,min(data),max(data)
endfor
for ivar=0,result2.nvars-1 do begin
    result=ncdf_varinq(ncid2,ivar)
    if result.name ne 'CH4' and result.name ne 'NOY' and result.name ne 'O3' then $
       ncdf_varget,ncid2,ncdf_varid(ncid2,result.name),data
    if result.name eq 'CH4' or result.name eq 'NOY' or result.name eq 'O3' then $
       ncdf_varget,ncid2,ncdf_varid(ncid2,result.name),data, OFFSET=[0,0,0,itime], COUNT=[nc,nr,nl,1]
    if result.name eq 'CH4' then ch4=data
    if result.name eq 'NOY' then noy=data
    if result.name eq 'O3' then  o3=data
;   print,ivar,result.name,min(data),max(data)
endfor
;
;============================================================
; Calculate Pressure
;============================================================
prl        = fltarr(nc,nr,nl)
Pzero       = P0/100.
; prl(i,j,k) = A(k)*PO + B(k)*PS(i,j)
FOR ilon = 0, nc-1 DO $
    FOR ilat = 0, nr-1 DO $
        FOR ialt = 0, nl-1 DO $
            prl(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat,itime)
;
; calculate theta, absolute vorticity, potential vorticity
;
thl=0.*tl
for L=0L,NL-1L do $
    THL(*,*,L)=TL(*,*,L)*(1000./PRL(*,*,L))^.286

eta=0.*tl
compvort,ul,vl,eta,alon,alat,nc,nr
;
; initialize surface arrays
;
tsfc=reform(tl(*,*,nl-1L))
thsfc=reform(thl(*,*,nl-1L))
usfc=reform(ul(*,*,nl-1L))
vsfc=reform(vl(*,*,nl-1L))
o3sfc=reform(o3(*,*,nl-1L))
ch4sfc=reform(ch4(*,*,nl-1L))
noysfc=reform(noy(*,*,nl-1L))
qdfsfc=0.*noysfc
pvsfc=0.*noysfc
;
; initialize isentropic arrays
;
UGRD=fltarr(nr,nc,nth)
VGRD=fltarr(nr,nc,nth)
IPVGRD=fltarr(nr,nc,nth)
PGRD=fltarr(nr,nc,nth)
QDFGRD=fltarr(nr,nc,nth)
O3GRD=fltarr(nr,nc,nth)
CH4GRD=fltarr(nr,nc,nth)
NOYGRD=fltarr(nr,nc,nth)
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
            LM1=K-1
            LP1=K+1
            IF K EQ 0 then LM1=0
            IF K EQ NL-1L then LP1=NL-1L
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
        qdfsfc(i,lat)=qdf(i,lat,nl-1)
        pvsfc(i,lat)=pv(i,lat,nl-1)

; LINEARLY INTERPOLATE IN THETA TO ISENTROPIC LEVEL
      for m=0,nth-1L do begin
      UGRD(LAT,I,M)=0.
      VGRD(LAT,I,M)=0.
      IPVGRD(LAT,I,M)=1.E12
      PGRD(LAT,I,M)=0.
      QDFGRD(LAT,I,M)=0.
      O3GRD(LAT,I,M)=0.
      CH4GRD(LAT,I,M)=0.
      NOYGRD(LAT,I,M)=0.

      for K=1,NL-1L do begin
      LM2 = K-2
      LM1 = K-1
      LP1 = K+1
      IF K EQ 1 then LM2 = 0
      IF K EQ NL-1L then LP1 = NL-1L

      IF THLEV(M) LE THL(I,LAT,LM1) AND THLEV(M) GT THL(I,LAT,K) THEN begin
      PLM1=PRL(I,LAT,LM1)
      PL=PRL(I,LAT,K)
      SCALE=(THLEV(M)-THL(I,LAT,K))/(THL(I,LAT,LM1)-THL(I,LAT,K))
      UGRD(LAT,I,m)=UL(I,LAT,K)+SCALE*(UL(I,LAT,LM1)-UL(I,LAT,K))
      VGRD(LAT,I,m)=VL(I,LAT,K)+SCALE*(VL(I,LAT,LM1)-VL(I,LAT,K))
      IPVGRD(LAT,I,m)=PV(I,LAT,K)+SCALE*(PV(I,LAT,LM1)-PV(I,LAT,K))
      PGRD(LAT,I,m)=PL^.286 + SCALE*(PLM1^.286-PL^.286)
      PGRD(LAT,I,m)=PGRD(LAT,I,m)^(1./.286)
      QDFGRD(LAT,I,M)=QDF(I,LAT,K)+SCALE*(QDF(I,LAT,LM1)-QDF(I,LAT,K))
      O3GRD(LAT,I,M)=O3(I,LAT,K)+SCALE*(O3(I,LAT,LM1)-O3(I,LAT,K))
      CH4GRD(LAT,I,M)=CH4(I,LAT,K)+SCALE*(CH4(I,LAT,LM1)-CH4(I,LAT,K))
      NOYGRD(LAT,I,M)=NOY(I,LAT,K)+SCALE*(NOY(I,LAT,LM1)-NOY(I,LAT,K))
      ENDIF

      endfor

; if desired theta surface intersects the ground then
; isentropic quantities are equal to the sfc values (at psfc)
      if thlev(m) le thsfc(i,lat) then begin
      PGRD(LAT,I,M)=psfc(i,lat)
      UGRD(LAT,I,M)=usfc(i,lat)
      VGRD(LAT,I,M)=vsfc(i,lat)
      IPVGRD(LAT,I,M)=pvsfc(i,lat)
      QDFGRD(LAT,I,M)=qdfsfc(i,lat)
      O3GRD(LAT,I,M)=o3sfc(i,lat)
      CH4GRD(LAT,I,M)=ch4sfc(i,lat)
      NOYGRD(LAT,I,M)=noysfc(i,lat)
      endif

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
; interpolate to MetO latitudes
;
PGRD2=fltarr(nru,ncu,nth)
UGRD2=fltarr(nru,ncu,nth)
VGRD2=fltarr(nru,ncu,nth)
IPVGRD2=fltarr(nru,ncu,nth)
QDFGRD2=fltarr(nru,ncu,nth)
O3GRD2=fltarr(nru,ncu,nth)
CH4GRD2=fltarr(nru,ncu,nth)
NOYGRD2=fltarr(nru,ncu,nth)
for m=0,nth-1 do begin
for i=0,ncu-1 do begin
    xp=ulon(i)
    for ii=0,nc-1 do begin
        ip1=ii+1
        if ii eq nc-1L then ip1=0L
        xlon=alon(ii)
        xlonp1=alon(ip1)
        if ii eq nc-1L then xlonp1=alon(ip1)+360.
        IF xp ge xlon AND xp le xlonp1 THEN begin
           xscale=(xp-xlon)/(xlonp1-xlon)
           for j=0,nru-1 do begin
               yp=ulat(j)
               for jj = 0, nr-2 do begin
                   jp1 = jj+1
                   xlat=alat(jj)
                   xlatp1=alat(jp1)
                   IF yp ge xlat AND yp le xlatp1 THEN begin
                   YSCALE=(yp-xlat)/(xlatp1-xlat)
;
; interpolate in longitude at 2 bounding latitudes
;
                   uj=UGRD(jj,ii,m)+XSCALE*(UGRD(jj,ip1,m)-UGRD(jj,ii,m))
                   vj=VGRD(jj,ii,m)+XSCALE*(VGRD(jj,ip1,m)-VGRD(jj,ii,m))
                   pvj=IPVGRD(jj,ii,m)+XSCALE*(IPVGRD(jj,ip1,m)-IPVGRD(jj,ii,m))
                   pj=PGRD(jj,ii,m)+XSCALE*(PGRD(jj,ip1,m)-PGRD(jj,ii,m))
                   qj=QDFGRD(jj,ii,m)+XSCALE*(QDFGRD(jj,ip1,m)-QDFGRD(jj,ii,m))
                   o3j=O3GRD(jj,ii,m)+XSCALE*(O3GRD(jj,ip1,m)-O3GRD(jj,ii,m))
                   ch4j=CH4GRD(jj,ii,m)+XSCALE*(CH4GRD(jj,ip1,m)-CH4GRD(jj,ii,m))
                   noyj=NOYGRD(jj,ii,m)+XSCALE*(NOYGRD(jj,ip1,m)-NOYGRD(jj,ii,m))

                   ujp1=UGRD(jp1,ii,m)+XSCALE*(UGRD(jp1,ip1,m)-UGRD(jp1,ii,m))
                   vjp1=VGRD(jp1,ii,m)+XSCALE*(VGRD(jp1,ip1,m)-VGRD(jp1,ii,m))
                   pvjp1=IPVGRD(jp1,ii,m)+XSCALE*(IPVGRD(jp1,ip1,m)-IPVGRD(jp1,ii,m))
                   pjp1=PGRD(jp1,ii,m)+XSCALE*(PGRD(jp1,ip1,m)-PGRD(jp1,ii,m))
                   qjp1=QDFGRD(jp1,ii,m)+XSCALE*(QDFGRD(jp1,ip1,m)-QDFGRD(jp1,ii,m))
                   o3jp1=O3GRD(jp1,ii,m)+XSCALE*(O3GRD(jp1,ip1,m)-O3GRD(jp1,ii,m))
                   ch4jp1=CH4GRD(jp1,ii,m)+XSCALE*(CH4GRD(jp1,ip1,m)-CH4GRD(jp1,ii,m))
                   noyjp1=NOYGRD(jp1,ii,m)+XSCALE*(NOYGRD(jp1,ip1,m)-NOYGRD(jp1,ii,m))
;
; interpolate in latitude from interpolated longitudes
;
                   UGRD2(j,i,m)=uj+yscale*(ujp1-uj)
                   VGRD2(j,i,m)=vj+yscale*(vjp1-vj)
                   IPVGRD2(j,i,m)=pvj+yscale*(pvjp1-pvj)
                   PGRD2(j,i,m)=pj+yscale*(pjp1-pj)
                   QDFGRD2(j,i,m)=qj+yscale*(qjp1-qj)
                   O3GRD2(j,i,m)=o3j+yscale*(o3jp1-o3j)
                   CH4GRD2(j,i,m)=ch4j+yscale*(ch4jp1-ch4j)
                   NOYGRD2(j,i,m)=noyj+yscale*(noyjp1-noyj)

                   ENDIF
               endfor
           endfor
        ENDIF
    endfor
endfor
endfor
;
; check
; 
sdate=strcompress(string(date(itime)),/remove_all)
rlev=2000.
;print,thlev
;read,'Enter theta surface ',rlev
index=where(thlev eq rlev)
ilev=index(0)
slev=string(rlev)
pp=transpose(qdfgrd2(*,*,ilev))
level=min(pp)+((max(pp)-min(pp))/float(nlvls))*findgen(nlvls)
!type=2^2+2^3
erase
contour,pp,ulon,ulat,levels=level,/cell_fill,c_color=col1,/noeras,xrange=[0.,360.],$
        yrange=[-90.,90.],title=sdate+'  '+slev+' K'
contour,pp,ulon,ulat,levels=level,/follow,c_color=0,/noeras,/overplot
contour,pp,ulon,ulat,levels=[0.],/follow,c_color=0,thick=3,/noeras,/overplot
;
; write daily theta file
;
ofile='/aura3/data/WACCM_data/Datfiles/waccm_'+sdate+'.nc'
print,'writing ',ofile
nocid = ncdf_create(ofile,/CLOBBER)
latdimid=ncdf_dimdef(nocid, 'number_of_latitudes' , nru)
londimid=ncdf_dimdef(nocid, 'number_of_longitudes', ncu)
levdimid=ncdf_dimdef(nocid, 'number_of_levels'    , nth)
lonsid = ncdf_vardef(nocid, 'longitudes',  londimid)
latsid = ncdf_vardef(nocid, 'latitudes' ,  latdimid)
levsid = ncdf_vardef(nocid, 'th_levels' ,  levdimid)
ipvid  = ncdf_vardef(nocid, 'ipv'       , [latdimid,londimid,levdimid])
prsid  = ncdf_vardef(nocid, 'press'     , [latdimid,londimid,levdimid])
uuuid  = ncdf_vardef(nocid, 'u_wind'    , [latdimid,londimid,levdimid])
vvvid  = ncdf_vardef(nocid, 'v_wind'    , [latdimid,londimid,levdimid])
qdfid  = ncdf_vardef(nocid, 'qdf'       , [latdimid,londimid,levdimid])
o3id   = ncdf_vardef(nocid, 'ozone'     , [latdimid,londimid,levdimid])
ch4id  = ncdf_vardef(nocid, 'methane'   , [latdimid,londimid,levdimid])
noyid  = ncdf_vardef(nocid, 'NOY'       , [latdimid,londimid,levdimid])
ncdf_control,nocid,/ENDEF
ncdf_varput, nocid, lonsid, ulon  , COUNT=[ncu]
ncdf_varput, nocid, latsid, ulat  , COUNT=[nru]
ncdf_varput, nocid, levsid, thlev , COUNT=[nth]
ncdf_varput, nocid, ipvid , ipvgrd2, COUNT=[nru,ncu,nth]
ncdf_varput, nocid, prsid , pgrd2  , COUNT=[nru,ncu,nth]
ncdf_varput, nocid, uuuid , ugrd2  , COUNT=[nru,ncu,nth]
ncdf_varput, nocid, vvvid , vgrd2  , COUNT=[nru,ncu,nth]
ncdf_varput, nocid, qdfid , qdfgrd2, COUNT=[nru,ncu,nth]
ncdf_varput, nocid, o3id  , o3grd2 , COUNT=[nru,ncu,nth]
ncdf_varput, nocid, ch4id , ch4grd2, COUNT=[nru,ncu,nth]
ncdf_varput, nocid, noyid , noygrd2, COUNT=[nru,ncu,nth]
ncdf_close,nocid
ENDFOR		; LOOP OVER TIMESTEPS
ncdf_close,ncid
ncdf_close,ncid2
end
