;
; /aura7/harvey/WACCM_data/Datfiles/Datfiles_Liu
; WA3548T08CO_2x.cam2.h2.2002-01-28-00000.nc files from Hanli
; read WACCM pressure files and store in binary format at
; the horizontal and vertical resolution of the UARS MetO data
; with extension to ~140 km
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
!p.background=mcolor
nlvls=30
col1=1+indgen(nlvls)*mcolor/nlvls
!noeras=1
;
; Ask interactive questions- get starting/ending dates
;
lstmn=1 & lstdy=1 & lstyr=2002
ledmn=1 & leddy=31 & ledyr=2002
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_Liu/WAX3548T08CO_2x.cam2.h2.'
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
;
; UARS Met Office grid (+ 24 pressure levels above 0.1 hPa)
;
ncold=96L
alonold=3.75*findgen(ncold)
wlonold=1.875+3.75*findgen(ncold)
nrold=72L
nr1old=73L
alatold=90.-2.5*findgen(nr1old)
wlatold=88.75-2.5*findgen(nrold)
p=[1000.,681.3,464.2,316.2,215.4,146.8,100.,68.13,46.42,31.62,21.54,14.68,10.,6.813,$
   4.642,3.162,2.154,1.468,1.,0.681,0.464,0.316,0.215,0.147,0.1,0.0681,0.0464,0.0316,0.0215,$
   0.0147,0.01,0.00681,0.00464,0.00316,0.00215,0.00147,0.001,0.000681,0.000464,0.000316,0.000215,$
   0.000147,0.0001,0.0000681,0.0000464,0.0000316,0.0000215,0.0000147,0.00001]
nlvold=n_elements(p)
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(i4)',iyr)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      date=syr+smn+sdy
      sdate=syr+'-'+smn+'-'+sdy
;
; read WACCM data
;
      ifile=dir+sdate+'-00000.nc'
      ofile=dir+sdate+'-00000.dat'
      ncid=ncdf_open(ifile)
      result0=ncdf_inquire(ncid)
      for idim=0,result0.ndims-1 do begin
          ncdf_diminq,ncid,idim,name,dim
          if name eq 'lon' then nc=dim
          if name eq 'lat' then nr=dim
          if name eq 'lev' then nl=dim
          if name eq 'time' then nt=dim
;         print,'read ',name,' dimension ',dim
      endfor
;
; loop over variables
;
      for ivar=0,result0.nvars-1 do begin
          result=ncdf_varinq(ncid,ivar)
          ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
          if result.name eq 'P0' then p0=data
          if result.name eq 'lat' then alat=data
          if result.name eq 'lon' then alon=data
          if result.name eq 'lev' then lev=data
          if result.name eq 'ilev' then ilev=data
          if result.name eq 'time' then time=data
          if result.name eq 'datesec' then time=data/86400.
          if result.name eq 'hyai' then hyai=data
          if result.name eq 'hybi' then hybi=data
          if result.name eq 'hyam' then hyam=data
          if result.name eq 'hybm' then hybm=data
          if result.name eq 'date' then date=data
          if result.name eq 'PS' then psfc=data     ;/100.
          if result.name eq 'T' then t3d=data
          if result.name eq 'U' then u3d=data
          if result.name eq 'V' then v3d=data
          if result.name eq 'e' then e3d=data
          if result.name eq 'CO' then co3d=data
          if result.name eq 'NO' then no3d=data
          if result.name eq 'QRL_TOT' then qrl3d=data
          if result.name eq 'QRS_TOT' then qrs3d=data
          if result.name eq 'O3' then  o3d=data
          if result.name eq 'Z3' then  z3d=data

          print,ivar,result.name,min(data),max(data)
      endfor
      ncdf_close,ncid
;
; make alon,alat,lev float arrays
;
      alon=float(alon)
      alat=float(alat)
      lev=float(lev)
;
;============================================================
; Calculate Pressure : p3d(i,j,k) = A(k)*PO + B(k)*PS(i,j)
;============================================================
      p3d        = fltarr(nc,nr,nl)
      Pzero      = P0
      FOR ilon = 0, nc-1 DO $
          FOR ilat = 0, nr-1 DO $
              FOR ialt = 0, nl-1 DO $
                  p3d(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat)
;
; convert dT/dt to dTh/dt
;
      q3d=qrl3d+qrs3d
      q3d=q3d*(100000./p3d)^0.286
;
; compute atmospheric density
; p=rho R T -> rho=P/RT where R=287 J/K kg. Pressure in Pascals.
;
      rho=p3d/(t3d*287.)
;
; convert species from (NO molecules/air molecules) to NO molecules/cm3
; assume the molecular weight of one molecule of air is 29 grams 
; (weight of O is 16, weight of N is 14, atm is 80% N2 and 20% O2)
; (molecules NO/m^3 air) = (mol NO/mol air) * (1 molecule air/29 grams) * 
;                          (1000 g air/1 kg air) * (AIR DENSITY/m^3 air) * (Avagadros #/1 mole air)
; Avagadros # = 6.022e23
;
      no3d=no3d * (1./29.) * (1000./1.) * rho * 6.022e23
      no3d=no3d/1.e6                                 ; divide by 1.e6 for m-3 to cm-3
      co3d=co3d * (1./29.) * (1000./1.) * rho * 6.022e23
      co3d=co3d/1.e6 
      o3d=o3d * (1./29.) * (1000./1.) * rho * 6.022e23
      o3d=o3d/1.e6
;
; make sure that 3-d arrays go from top down
;
      wlat=alat & wlon=alon
      nc=n_elements(alon) & nc1=n_elements(alon)
      nr=n_elements(alat) & nr1=n_elements(alat)
      nlv=n_elements(lev)
;
; check zonal mean temperature and U
;
tbarg=fltarr(nr1,nlv)
ubarg=fltarr(nr1,nlv)
for k=0L,nlv-1L do begin
for j=0L,nr1-1L do begin
    tbarg(j,k)=total(t3d(*,j,k))/float(nc)
    ubarg(j,k)=total(u3d(*,j,k))/float(nc)
endfor
endfor
erase
!type=2^2+2^3
set_viewport,.15,.85,.15,.85
contour,tbarg,alat,lev,levels=130.+10.*findgen(30),/noerase,title=sdate+' T and U',/fill,$
        c_color=col1,/ylog,yrange=[max(lev),min(lev)],color=0,charsize=2
contour,tbarg,alat,lev,levels=130.+10.*findgen(30),/noerase,/follow,color=0,/overplot
contour,ubarg,alat,lev,levels=-100.+10.*findgen(10),/noerase,/overplot,color=mcolor,thick=2,c_linestyle=5
contour,ubarg,alat,lev,levels=10.+10.*findgen(10),/noerase,/overplot,color=0,thick=2
;
; convert pressure to hPa, create log pressure array, convert dTh/dt to (K/day).
;
      p3d=p3d/100.
      ap3d=alog(p3d)
      q3d=q3d*86400.
;
; 2d latitude longitude arrays for Q,Z,T,NO,CO,e
;
      x2d=fltarr(ncold,nr1old)
      y2d=fltarr(ncold,nr1old)
      for i=0L,ncold-1 do y2d(i,*)=alatold
      for j=0L,nr1old-1 do x2d(*,j)=alonold
;
; interpolate new grid to old resolution (alat->alatold, wlat->wlatold, lev->p)
; new pressure goes from bottom up
;
      qnew=fltarr(ncold,nr1old,nlvold)
      znew=fltarr(ncold,nr1old,nlvold)
      tnew=fltarr(ncold,nr1old,nlvold)
      nonew=fltarr(ncold,nr1old,nlvold)
      conew=fltarr(ncold,nr1old,nlvold)
      o3new=fltarr(ncold,nr1old,nlvold)

      unew=fltarr(ncold,nrold,nlvold)
      vnew=fltarr(ncold,nrold,nlvold)
      for kk=0L,nlvold-1L do begin
          pp=alog(p(kk))
          for ii=0L,ncold-1L do begin
          xp=alonold(ii)
          for jj=0L,nr1old-1L do begin
          yp=alatold(jj)
          for i=0L,nc-1L do begin
              ip1=i+1
              if i eq nc-1L then ip1=0
              xlon=alon(i)
              xlonp1=alon(ip1)
              if i eq nc-1L then xlonp1=360.+alon(ip1)
              if xp ge xlon and xp le xlonp1 then begin
                 xscale=(xp-xlon)/(xlonp1-xlon)
                 for j=0L,nr1-2L do begin
                     jp1=j+1
                     xlat=alat(j)
                     xlatp1=alat(jp1)
                     if yp ge xlat and yp le xlatp1 then begin
                        yscale=(yp-xlat)/(xlatp1-xlat)
                        for k=0L,nlv-2L do begin
                            kp1=k+1
;
; impose a more rigorous vertical interpolation scale factor based on
; ALL 8 surrounding gridpoints, not just 2: (j,i,k) and (j,i,kp1)
;
                            pj1=ap3d(i,j,k)+xscale*(ap3d(ip1,j,k)-ap3d(i,j,k))
                            pjp1=ap3d(i,jp1,k)+xscale*(ap3d(ip1,jp1,k)-ap3d(i,jp1,k))
                            pj2=ap3d(i,j,kp1)+xscale*(ap3d(ip1,j,kp1)-ap3d(i,j,kp1))
                            pjp2=ap3d(i,jp1,kp1)+xscale*(ap3d(ip1,jp1,kp1)-ap3d(i,jp1,kp1))
                            p1=pj1+yscale*(pjp1-pj1)	; top
                            p2=pj2+yscale*(pjp2-pj2)	; bottom
                            if pp ge p1 and pp le p2 then begin
                               pscale=(pp-p1)/(p2-p1)
                               qj1=q3d(i,j,k)+xscale*(q3d(ip1,j,k)-q3d(i,j,k))
                               qjp1=q3d(i,jp1,k)+xscale*(q3d(ip1,jp1,k)-q3d(i,jp1,k))
                               qj2=q3d(i,j,kp1)+xscale*(q3d(ip1,j,kp1)-q3d(i,j,kp1))
                               qjp2=q3d(i,jp1,kp1)+xscale*(q3d(ip1,jp1,kp1)-q3d(i,jp1,kp1))
                               q1=qj1+yscale*(qjp1-qj1)
                               q2=qj2+yscale*(qjp2-qj2)
                               qnew(ii,jj,kk)=q1+pscale*(q2-q1)
                               tj1=t3d(i,j,k)+xscale*(t3d(ip1,j,k)-t3d(i,j,k))
                               tjp1=t3d(i,jp1,k)+xscale*(t3d(ip1,jp1,k)-t3d(i,jp1,k))
                               tj2=t3d(i,j,kp1)+xscale*(t3d(ip1,j,kp1)-t3d(i,j,kp1))
                               tjp2=t3d(i,jp1,kp1)+xscale*(t3d(ip1,jp1,kp1)-t3d(i,jp1,kp1))
                               t1=tj1+yscale*(tjp1-tj1)
                               t2=tj2+yscale*(tjp2-tj2)
                               tnew(ii,jj,kk)=t1+pscale*(t2-t1)
;if yp eq 90 then stop,tnew(ii,jj,kk),t1,t2,pscale,t1+pscale*(t2-t1)
                               zj1=z3d(i,j,k)+xscale*(z3d(ip1,j,k)-z3d(i,j,k))
                               zjp1=z3d(i,jp1,k)+xscale*(z3d(ip1,jp1,k)-z3d(i,jp1,k))
                               zj2=z3d(i,j,kp1)+xscale*(z3d(ip1,j,kp1)-z3d(i,j,kp1))
                               zjp2=z3d(i,jp1,kp1)+xscale*(z3d(ip1,jp1,kp1)-z3d(i,jp1,kp1))
                               z1=zj1+yscale*(zjp1-zj1)
                               z2=zj2+yscale*(zjp2-zj2)
                               znew(ii,jj,kk)=z1+pscale*(z2-z1)
;
; add NO, CO, and e
;
                               zj1=no3d(i,j,k)+xscale*(no3d(ip1,j,k)-no3d(i,j,k))
                               zjp1=no3d(i,jp1,k)+xscale*(no3d(ip1,jp1,k)-no3d(i,jp1,k))
                               zj2=no3d(i,j,kp1)+xscale*(no3d(ip1,j,kp1)-no3d(i,j,kp1))
                               zjp2=no3d(i,jp1,kp1)+xscale*(no3d(ip1,jp1,kp1)-no3d(i,jp1,kp1))
                               z1=zj1+yscale*(zjp1-zj1)
                               z2=zj2+yscale*(zjp2-zj2)
                               nonew(ii,jj,kk)=z1+pscale*(z2-z1)
                               zj1=co3d(i,j,k)+xscale*(co3d(ip1,j,k)-co3d(i,j,k))
                               zjp1=co3d(i,jp1,k)+xscale*(co3d(ip1,jp1,k)-co3d(i,jp1,k))
                               zj2=co3d(i,j,kp1)+xscale*(co3d(ip1,j,kp1)-co3d(i,j,kp1))
                               zjp2=co3d(i,jp1,kp1)+xscale*(co3d(ip1,jp1,kp1)-co3d(i,jp1,kp1))
                               z1=zj1+yscale*(zjp1-zj1)
                               z2=zj2+yscale*(zjp2-zj2)
                               conew(ii,jj,kk)=z1+pscale*(z2-z1)
                               zj1=o3d(i,j,k)+xscale*(o3d(ip1,j,k)-o3d(i,j,k))
                               zjp1=o3d(i,jp1,k)+xscale*(o3d(ip1,jp1,k)-o3d(i,jp1,k))
                               zj2=o3d(i,j,kp1)+xscale*(o3d(ip1,j,kp1)-o3d(i,j,kp1))
                               zjp2=o3d(i,jp1,kp1)+xscale*(o3d(ip1,jp1,kp1)-o3d(i,jp1,kp1))
                               z1=zj1+yscale*(zjp1-zj1)
                               z2=zj2+yscale*(zjp2-zj2)
                               o3new(ii,jj,kk)=z1+pscale*(z2-z1)

                               goto,jumplev1
                            endif
                        endfor
jumplev1:
;
; is MetO pressure level below lowest GEOS-5 level? (is p gt p3d(*,*,nlv-1)?)
;
                        k=nlv-1L
                        pj1=ap3d(i,j,k)+xscale*(ap3d(ip1,j,k)-ap3d(i,j,k))
                        pjp1=ap3d(i,jp1,k)+xscale*(ap3d(ip1,jp1,k)-ap3d(i,jp1,k))
                        up=pj1+yscale*(pjp1-pj1)
                        if pp gt up then begin
                           qj1=q3d(i,j,k)+xscale*(q3d(ip1,j,k)-q3d(i,j,k))
                           qjp1=q3d(i,jp1,k)+xscale*(q3d(ip1,jp1,k)-q3d(i,jp1,k))
                           qnew(ii,jj,kk)=qj1+yscale*(qjp1-qj1)
                           tj1=t3d(i,j,k)+xscale*(t3d(ip1,j,k)-t3d(i,j,k))
                           tjp1=t3d(i,jp1,k)+xscale*(t3d(ip1,jp1,k)-t3d(i,jp1,k))
                           tnew(ii,jj,kk)=tj1+yscale*(tjp1-tj1)
                           zj1=z3d(i,j,k)+xscale*(z3d(ip1,j,k)-z3d(i,j,k))
                           zjp1=z3d(i,jp1,k)+xscale*(z3d(ip1,jp1,k)-z3d(i,jp1,k))
                           znew(ii,jj,kk)=zj1+yscale*(zjp1-zj1)
;
; add NO, CO, and e
;
                           zj1=no3d(i,j,k)+xscale*(no3d(ip1,j,k)-no3d(i,j,k))
                           zjp1=no3d(i,jp1,k)+xscale*(no3d(ip1,jp1,k)-no3d(i,jp1,k))
                           nonew(ii,jj,kk)=zj1+yscale*(zjp1-zj1)
                           zj1=co3d(i,j,k)+xscale*(co3d(ip1,j,k)-co3d(i,j,k))
                           zjp1=co3d(i,jp1,k)+xscale*(co3d(ip1,jp1,k)-co3d(i,jp1,k))
                           conew(ii,jj,kk)=zj1+yscale*(zjp1-zj1)
                           zj1=o3d(i,j,k)+xscale*(o3d(ip1,j,k)-o3d(i,j,k))
                           zjp1=o3d(i,jp1,k)+xscale*(o3d(ip1,jp1,k)-o3d(i,jp1,k))
                           o3new(ii,jj,kk)=zj1+yscale*(zjp1-zj1)

                        endif
                     endif
                 endfor
              endif
          endfor
          endfor	; z and T latitudes
;
; interpolate winds to staggered grid
;
          xp=wlonold(ii)
          for jj=0L,nrold-1L do begin
          yp=wlatold(jj)
          for i=0L,nc-1L do begin
              ip1=i+1
              if i eq nc-1L then ip1=0
              xlon=wlon(i)
              xlonp1=wlon(ip1)
              if i eq nc-1L then xlonp1=360.+wlon(ip1)
              if xp ge xlon and xp le xlonp1 then begin
                 xscale=(xp-xlon)/(xlonp1-xlon)
                 for j=0L,nr-2L do begin
                     jp1=j+1
                     xlat=wlat(j)
                     xlatp1=wlat(jp1)
                     if yp ge xlat and yp le xlatp1 then begin
                        yscale=(yp-xlat)/(xlatp1-xlat)
                        for k=0L,nlv-2L do begin
                            kp1=k+1
                            pj1=ap3d(i,j,k)+xscale*(ap3d(ip1,j,k)-ap3d(i,j,k))
                            pjp1=ap3d(i,jp1,k)+xscale*(ap3d(ip1,jp1,k)-ap3d(i,jp1,k))
                            pj2=ap3d(i,j,kp1)+xscale*(ap3d(ip1,j,kp1)-ap3d(i,j,kp1))
                            pjp2=ap3d(i,jp1,kp1)+xscale*(ap3d(ip1,jp1,kp1)-ap3d(i,jp1,kp1))
                            p1=pj1+yscale*(pjp1-pj1)    ; top
                            p2=pj2+yscale*(pjp2-pj2)    ; bottom
                            if pp ge p1 and pp le p2 then begin
                               pscale=(pp-p1)/(p2-p1)
                               uj1=u3d(i,j,k)+xscale*(u3d(ip1,j,k)-u3d(i,j,k))
                               ujp1=u3d(i,jp1,k)+xscale*(u3d(ip1,jp1,k)-u3d(i,jp1,k))
                               uj2=u3d(i,j,kp1)+xscale*(u3d(ip1,j,kp1)-u3d(i,j,kp1))
                               ujp2=u3d(i,jp1,kp1)+xscale*(u3d(ip1,jp1,kp1)-u3d(i,jp1,kp1))
                               u1=uj1+yscale*(ujp1-uj1)
                               u2=uj2+yscale*(ujp2-uj2)
                               unew(ii,jj,kk)=u1+pscale*(u2-u1)
                               vj1=v3d(i,j,k)+xscale*(v3d(ip1,j,k)-v3d(i,j,k))
                               vjp1=v3d(i,jp1,k)+xscale*(v3d(ip1,jp1,k)-v3d(i,jp1,k))
                               vj2=v3d(i,j,kp1)+xscale*(v3d(ip1,j,kp1)-v3d(i,j,kp1))
                               vjp2=v3d(i,jp1,kp1)+xscale*(v3d(ip1,jp1,kp1)-v3d(i,jp1,kp1))
                               v1=vj1+yscale*(vjp1-vj1)
                               v2=vj2+yscale*(vjp2-vj2)
                               vnew(ii,jj,kk)=v1+pscale*(v2-v1)
                               goto,jumplev2
                            endif
                        endfor
jumplev2:
                        k=nlv-1L
                        pj1=ap3d(i,j,k)+xscale*(ap3d(ip1,j,k)-ap3d(i,j,k))
                        pjp1=ap3d(i,jp1,k)+xscale*(ap3d(ip1,jp1,k)-ap3d(i,jp1,k))
                        up=pj1+yscale*(pjp1-pj1)
                        if pp gt up then begin
                           uj1=u3d(i,j,k)+xscale*(u3d(ip1,j,k)-u3d(i,j,k))
                           ujp1=u3d(i,jp1,k)+xscale*(u3d(ip1,jp1,k)-u3d(i,jp1,k))
                           unew(ii,jj,kk)=uj1+yscale*(ujp1-uj1)
                           vj1=v3d(i,j,k)+xscale*(v3d(ip1,j,k)-v3d(i,j,k))
                           vjp1=v3d(i,jp1,k)+xscale*(v3d(ip1,jp1,k)-v3d(i,jp1,k))
                           vnew(ii,jj,kk)=vj1+yscale*(vjp1-vj1)
                        endif
                     endif
                 endfor
              endif
          endfor
          endfor
          endfor	; loop over longitude
          print,'interpolated to ',p(kk),znew(10,10,kk),tnew(10,10,kk),unew(10,10,kk),vnew(10,10,kk)
;
; check
;
;erase
;nlvls=21
;col1=1+indgen(nlvls)*mcolor/nlvls
;slev=strtrim(string(p(kk)),2)+' hPa'
;
;set_viewport,0.1,0.45,0.5,0.75
;grd1=reform(znew(*,*,kk))
;print,'Gp ',min(grd1),max(grd1)
;index=where(grd1 ne 0.)
;imin=min(grd1(index)) & imax=max(grd1)
;level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
;map_set,0,180,0,/contin,/grid,title='Gp '+slev,color=0,/noeras
;contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
;contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
;map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
;set_viewport,0.55,0.9,0.5,0.75
;grd1=reform(tnew(*,*,kk))
;print,'Tp ',min(grd1),max(grd1)
;index=where(grd1 ne 0.)
;imin=min(grd1(index)) & imax=max(grd1)
;level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
;map_set,0,180,0,/contin,/grid,title='Tp '+slev,color=0,/noeras
;contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
;contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
;map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
;set_viewport,0.1,0.45,0.15,0.4
;grd1=reform(nonew(*,*,kk))
;print,'NO ',min(grd1),max(grd1)
;imin=min(grd1) & imax=max(grd1)
;level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
;map_set,0,180,0,/contin,/grid,title='NO '+slev,color=0,/noeras
;contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
;contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
;map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
;set_viewport,0.55,0.9,0.15,0.4
;grd1=reform(enew(*,*,kk))
;print,'e ',min(grd1),max(grd1)
;imin=min(grd1) & imax=max(grd1)
;level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
;map_set,0,180,0,/contin,/grid,title='e '+slev,color=0,/noeras
;contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
;contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
;map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
      endfor		; loop over altitude
;
; loop over pressure levels and output /f77 binary to match MetO format
;
      close,1
      openw,1,ofile,/f77
      for kk=0L,nlvold-1L do begin
          plevel=p(kk)
          writeu,1,plevel
          z=reform(znew(*,*,kk))
          t=reform(tnew(*,*,kk))
          u=reform(unew(*,*,kk))
          v=reform(vnew(*,*,kk))
          q=reform(qnew(*,*,kk))
          no=reform(nonew(*,*,kk))
          co=reform(conew(*,*,kk))
          o3=reform(o3new(*,*,kk))
          writeu,1,u,v,t,z,q,no,co,o3
          print,plevel,u(0,0),v(0,0),t(0,0),z(0,0),q(0,0),no(0,0),co(0,0),o3(0,0)
      endfor
      close,1
goto,jump
end
