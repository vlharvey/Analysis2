;
; read nc files from Doug Kinnison of Surface Area Density of Ice
; /aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_CNRM-ACM_2_T3I_ice_*nc
;
; NOTE: time:units = "hours since 1990-01-01 00:00" and units are m-1
;
; 10 day instantaneous fields
;
; restore pressure conversion variables and compute pressure levels
;
ncfile='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_CNRM-ACM_T1_hyblev60.nc'
ncid=ncdf_open(ncfile)
result0=ncdf_inquire(ncid)
for idim=0,result0.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'lev' then nl=dim
;   print,'read ',name,' dimension ',dim
endfor
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    if result.name eq 'hyam' then hyam=data
    if result.name eq 'hybm' then hybm=data
    if result.name eq 'p0' then p0=data
endfor
plev=hyam*p0 + hybm*air_pressure_at_surface

dir='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_ice'
spawn,'ls '+dir+'*nc',ncfiles
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin
    ncfile=ncfiles(ifile)
    print,'opening '+ncfile
    ncid=ncdf_open(ncfile)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'lon' then nc=dim
        if name eq 'lat' then nr=dim
        if name eq 'plev' then nl=dim
        if name eq 'time' then nt=dim
;       print,'read ',name,' dimension ',dim
    endfor
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)
;       print,result.name
        if result.name ne 'sad_ice' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if result.name eq 'lat' then alat=data
        if result.name eq 'lon' then alon=data
        if result.name eq 'plev' then lev=data
        if result.name eq 'time' then begin
           time=float(data)/24.
;
; set time = 0 to 1/1/1990 and build date array
;
           lstdy=1 & lstmn=1 & lstyr=1990
           leddy=1 & ledmn=1 & ledyr=2010
           lstday=0
           ledday=0
           z = stddat_noleap_cnrm(lstmn,lstdy,lstyr,lstday)
           z = stddat_noleap_cnrm(ledmn,leddy,ledyr,ledday)
           if ledday lt lstday then stop,' Wrong dates! '
           kday=ledday-lstday+1L
           sdate=strarr(nt)
           rdate=fltarr(nt)
;
; Compute initial Julian date
;
           iyr = lstyr
           idy = lstdy
           imn = lstmn
           z = kgmt(imn,idy,iyr,iday)
           iday = iday - 1

; --- Loop over all days --------
           jump: iday = iday + 1
                 kdate,float(iday),iyr,imn,idy
                 ckday,iday,iyr
                 z = stddat_noleap_cnrm(imn,idy,iyr,ndays)
                 if imn eq 2 and idy eq 29 then begin
                    print,'skip leap day ',iday,iyr,imn,idy,ndays
                    goto,skipleap
                 endif
                 if ndays gt max(time) then goto,jumpout

                 good=where(time eq ndays)
                 if n_elements(good) gt 1L then print,'mult times',time(good)
                 if good(0) ne -1L then begin
print,iday,iyr,imn,idy,' CNRM day ',ndays
                    syr=string(FORMAT='(I4)',iyr)
                    smn=string(FORMAT='(I2.2)',imn)
                    sdy=string(FORMAT='(I2.2)',idy)
                    sdate(good)=syr+smn+sdy
                    rdate(good)=time(good)
                 endif
skipleap:
           goto,jump
           jumpout:
           index=where(sdate ne '',nn)
           sdate=sdate(index)
           rdate=rdate(index)
;
; check. in the 2000-2005 file there are non-unique times. these are filered out above
;
           if nn ne nt then begin
              for ii=0L,nt-1L do begin
                  index=where(rdate eq time(ii))
;print,ii,time(ii),rdate(index(0)),index(0)
                  if index(0) eq -1L then print,time(ii)
              endfor
              stop
           endif
        endif   ; if time variable
;
; read SAD data one time step at a time and save daily file. units = "m-1"
;
        if result.name eq 'sad_ice' then begin
           count = [nc,nr,nl,1]
           for n=0L,nt-1L do begin
               sdate0=sdate(n)
               ofile=dir+'_'+sdate0+'.sav'
               print,ofile,time(n)
               dum=findfile(ofile)
               if dum(0) ne '' then begin
                  print,'duplicate ',ofile
                  goto,jumpday
               endif
               offset = [0,0,0,n]
               ncdf_varget,ncid,ncdf_varid(ncid,result.name),sadgrd,count=count,offset=offset
               save,file=ofile,alon,alat,lev,sadgrd
jumpday:
           endfor	; loop over time steps
        endif
    endfor	; loop over variables
    ncdf_close,ncid
endfor	; loop over files
end
