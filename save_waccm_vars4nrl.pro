; 
; save NO, NO2, O3, T at every timestep in 1 file per variable
; save species in molecules/cm3
; send to NRL in support of CHARM
; VLH 11/15/10
;
@stddat
@kgmt
@ckday
@kdate

lstmn=1
lstdy=1
lstyr=2002
ledmn=1
leddy=10
ledyr=2002
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1
ntimes=long(kday*48.)
dir='/Volumes/data/WACCM/TEMAOA.cam2.h1.'

rtd=double(180./!pi)
dtr=1./rtd
ks=1.931853d-3
ecc=0.081819
gamma45=9.80

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
;     syr=string(FORMAT='(i4)',iyr)
      syr='0001'
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
;
; read WACCM data
;
      spawn,'ls '+dir+syr+'-'+smn+'-'+sdy+'-*.nc',ncfiles
      for ifile=0L,n_elements(ncfiles)-1L do begin
          result=strsplit(ncfiles(ifile),'-',/extract)
          result2=strsplit(result(3),'.',/extract)
          sec=result2(0)
          shrs=string(FORMAT='(f5.1)',float(sec)/3600.+((idy-1)*24.))
          sdate=syr+'-'+smn+'-'+sdy+'-'+sec
          print,sdate
          ncfile=ncfiles(ifile)
          ncid=ncdf_open(ncfile)
          result0=ncdf_inquire(ncid)
          for idim=0,result0.ndims-1 do begin
              ncdf_diminq,ncid,idim,name,dim
              if name eq 'lon' then nc=dim
              if name eq 'lat' then nr=dim
              if name eq 'lev' then nl=dim
              if name eq 'time' then nt=dim
;             print,'read ',name,' dimension ',dim
          endfor
;
; loop over variables
;
          for ivar=0,result0.nvars-1 do begin
              result=ncdf_varinq(ncid,ivar)

              if result.name eq 'P0' or result.name eq 'lat' or result.name eq 'lon' or result.name eq 'lev' or result.name eq 'hyai' or $
                 result.name eq 'hybi' or result.name eq 'hyam' or result.name eq 'hybm' or result.name eq 'date' or $
                 result.name eq 'PS' or result.name eq 'T' or result.name eq 'U' or result.name eq 'V' or result.name eq 'O' or $
                 result.name eq 'NO2' or result.name eq 'CO' or result.name eq 'NO' or result.name eq 'O3' or $
                 result.name eq 'O2' or result.name eq 'Z3' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data

              if result.name eq 'P0' then p0=data
              if result.name eq 'lat' then lat=data
              if result.name eq 'lon' then lon=data
              if result.name eq 'lev' then lev=data
              if result.name eq 'hyai' then hyai=data
              if result.name eq 'hybi' then hybi=data
              if result.name eq 'hyam' then hyam=data
              if result.name eq 'hybm' then hybm=data
              if result.name eq 'date' then date=data
              if result.name eq 'PS' then psfc=data	;/100.
              if result.name eq 'T' then tgrd=data
              if result.name eq 'U' then ugrd=data
              if result.name eq 'V' then vgrd=data
              if result.name eq 'NO2' then no2grd=data
              if result.name eq 'CO' then cogrd=data
              if result.name eq 'NO' then nogrd=data
              if result.name eq 'O3' then  o3grd=data
              if result.name eq 'O2' then  o2grd=data
              if result.name eq 'O' then  ogrd=data
              if result.name eq 'Z3' then  ggrd=data/1000.

;             print,ivar,result.name,min(data),max(data)
jumpvar:
          endfor
          ncdf_close,ncid
;
; convert geopotential to geometric height
;  
          zgrd=0.*ggrd
          for k=0L,nl-1L do begin
              for j=0L,nr-1L do begin
                  sin2=sin( (lat(j)*dtr)^2.0 )
                  numerator=1.0+ks*sin2
                  denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
                  gammas=gamma45*(numerator/denominator)
                  r=6378.137/(1.006803-(0.006706*sin2))
                  zgrd(*,j,k)=(r*ggrd(*,j,k))/ ( (gammas/gamma45)*r - ggrd(*,j,k) )
              endfor
          endfor
;
; Calculate Pressure : pgrd(i,j,k) = A(k)*PO + B(k)*PS(i,j)
;
          pgrd        = fltarr(nc,nr,nl)
          Pzero       = P0      ;/100.
          FOR ilon = 0, nc-1 DO $
              FOR ilat = 0, nr-1 DO $
                  FOR ialt = 0, nl-1 DO $
                      pgrd(ilon,ilat,ialt) = hyam(ialt)*Pzero + hybm(ialt)*PSFC(ilon,ilat)
;
; compute atmospheric density
; p=rho R T -> rho=P/RT where R=287 J/K kg. Pressure in Pascals.
;
          rho=pgrd/(tgrd*287.)
;
; to convert species from (NO molecules/air molecules) to NO molecules/cm3
; assume the molecular weight of one molecule of air is 29 grams (weight of O is 16, weight of N is 14, atm is 80% N2 and 20% O2)
; (mol NO/mol air) * (1 molecule air/29 grams) * (1000 g air/1 kg air) * (AIR DENSITY/m^3 air) * (Avagadros #/1 mole NO) = (molecules NO/m^3 air)
; Avagadros # = 6.022e23
;
          no_conc=nogrd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6	; divide by 1.e6 for m-3 to cm-3
          no2_conc=no2grd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6
          o3_conc=o3grd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6
          co_conc=cogrd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6
          o2_conc=o2grd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6
          o_conc=ogrd * (1./29.) * (1000./1.) * rho * 6.022e23 / 1.e6
;
; retain all data
;
          if icount eq 0L then begin
             sdate_all=strarr(ntimes)
             no_data=fltarr(nc,nr,nl,ntimes)
             no2_data=fltarr(nc,nr,nl,ntimes)
             o3_data=fltarr(nc,nr,nl,ntimes)
             o2_data=fltarr(nc,nr,nl,ntimes)
             gp_data=fltarr(nc,nr,nl,ntimes)
             co_data=fltarr(nc,nr,nl,ntimes)
             tp_data=fltarr(nc,nr,nl,ntimes)
             u_data=fltarr(nc,nr,nl,ntimes)
             v_data=fltarr(nc,nr,nl,ntimes)
             z_data=fltarr(nc,nr,nl,ntimes)
             o_data=fltarr(nc,nr,nl,ntimes)
          endif
          sdate_all(icount)=sdate
          no_data(*,*,*,icount)=no_conc
          no2_data(*,*,*,icount)=no2_conc
          o3_data(*,*,*,icount)=o3_conc
          o2_data(*,*,*,icount)=o2_conc
          o_data(*,*,*,icount)=o_conc
          gp_data(*,*,*,icount)=ggrd
          co_data(*,*,*,icount)=co_conc
          tp_data(*,*,*,icount)=tgrd
          u_data(*,*,*,icount)=ugrd
          v_data(*,*,*,icount)=vgrd
          z_data(*,*,*,icount)=zgrd
          icount=icount+1L
      endfor	; loop over time steps
goto, jump
;
; save files
;
saveit:
save,file='waccm_no_data.sav',no_data,lon,lat,lev,sdate_all
save,file='waccm_no2_data.sav',no2_data,lon,lat,lev,sdate_all
save,file='waccm_o3_data.sav',o3_data,lon,lat,lev,sdate_all
save,file='waccm_o2_data.sav',o2_data,lon,lat,lev,sdate_all
save,file='waccm_o_data.sav',o_data,lon,lat,lev,sdate_all
save,file='waccm_gp_data.sav',gp_data,lon,lat,lev,sdate_all
save,file='waccm_co_data.sav',co_data,lon,lat,lev,sdate_all
save,file='waccm_tp_data.sav',tp_data,lon,lat,lev,sdate_all
save,file='waccm_u_data.sav',u_data,lon,lat,lev,sdate_all
save,file='waccm_v_data.sav',v_data,lon,lat,lev,sdate_all
save,file='waccm_z_data.sav',z_data,lon,lat,lev,sdate_all
end
