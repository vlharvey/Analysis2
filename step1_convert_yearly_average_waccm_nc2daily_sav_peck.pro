;
; yearly average files
; read WACCM netcdf data from Ethan Peck
; convert yearly to daily files
;
; /aura3/data/WACCM_data/Datfiles/noaurfcoaverageYear.nc
;
;dir='/aura3/data/WACCM_data/Datfiles/noaurfcoaverageYear'
dir='/Volumes/earth/aura3/data/WACCM_data/Datfiles/noaurfplaverageYear'
;
; read WACCM data
;
    ncfile0='/Volumes/earth/aura3/data/WACCM_data/Datfiles/noaurfpl_FWaverageYear.nc'
    print,ncfile0
    ncid=ncdf_open(ncfile0)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'lon' then nc=dim
        if name eq 'lat' then nr=dim
        if name eq 'lev' then nl=dim
        if name eq 'time' then nt=dim
        print,'read ',name,' dimension ',dim
    endfor
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)
        if result.name eq 'NOX' or result.name eq 'O3' or result.name eq 'NOY' or result.name eq 'VTH3d' or $
           result.name eq 'UV3d' or result.name eq 'UW3d' or result.name eq 'TH' or result.name eq 'CO' or $
           result.name eq 'CH4' or result.name eq 'NO2' or result.name eq 'CLONO2' or result.name eq 'OMEGA' then goto,jumpvar1
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if result.name eq 'P0' then p0=data
        if result.name eq 'hyai' then hyai=data
        if result.name eq 'hybi' then hybi=data
        if result.name eq 'hyam' then hyam=data
        if result.name eq 'hybm' then hybm=data
        if result.name eq 'PS' then ps=data     ;/100.
        if result.name eq 'lat' then alat=data
        if result.name eq 'lon' then alon=data
        if result.name eq 'lev' then lev=data
        if result.name eq 'time' then time=data
        if result.name eq 'date' then date=data
        if result.name eq 'T' then t4d=data
        if result.name eq 'U' then u4d=data
        if result.name eq 'V' then v4d=data
        if result.name eq 'Z3' then g4d=data/1000.
        print,ivar,result.name,min(data),max(data)
jumpvar1:
    endfor
    ncdf_close,ncid
    sdate=strcompress(date,/remove_all)
    index=where(date lt 1001)
    sdate(index)='0'+sdate(index)
;
; Calculate 3d Pressure: p(i,j,k,n) = A(k)*PO + B(k)*PS(i,j,n) in Pascals
;
    P0=100000.
    p4d=fltarr(nc,nr,nl,nt)
    Pzero=P0
    FOR ilon=0,nc-1 DO $
        FOR ilat=0,nr-1 DO $
            FOR ialt=0,nl-1 DO $
                p4d(ilon,ilat,ialt,*)=(hyam(ialt)*Pzero + hybm(ialt)*PS(ilon,ilat,*)) / 100.
;
; IDL save file for each day
;
    for iday=0L,nt-1L do begin 
        ofile=dir+'_'+sdate(iday)+'.sav'
        print,ofile
        pgrd=reform(p4d(*,*,*,iday))
        tgrd=reform(t4d(*,*,*,iday))
        ugrd=reform(u4d(*,*,*,iday))
        vgrd=reform(v4d(*,*,*,iday))
        save,file=ofile,alon,alat,pgrd,tgrd,ugrd,vgrd
    endfor	; loop over days
end
