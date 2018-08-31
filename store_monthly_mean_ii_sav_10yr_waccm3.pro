;
; store multi-year monthly mean inertial instability frequency
;
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_TNV3/wa3_tnv3_'
spawn,'rm -f '+dir+'*ii_*-*'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
spawn,'ls '+dir+'1*nc3 '+dir+'2*nc3',allfiles
;
; extract month numbers from filenames
;
wmon=strarr(n_elements(allfiles))
for n=0L,n_elements(allfiles)-1L do begin
    result=strsplit(allfiles(n),'_',/extract)
    wmon(n)=strmid(result(4),4,2)
endfor
for m=0,nmonths-1 do begin
;
; pull out days in this month from all years
;
    smon=string(format='(i2.2)',m+1L)
    index=where(wmon eq smon,nfile)
    ifiles=allfiles(index)
    print,nfile,' ',mon(m),' days'
    for n=0L,nfile-1 do begin
        ifile=ifiles(n)
        print,ifile
        ncid=ncdf_open(ifile)
        nr=0L & nc=0L & nth=0L
        ncdf_diminq,ncid,0,name,nr
        ncdf_diminq,ncid,1,name,nc
        ncdf_diminq,ncid,2,name,nth
        alon=fltarr(nr) & alat=fltarr(nc) & th=fltarr(nth)
        ncdf_varget,ncid,0,alon
        ncdf_varget,ncid,1,alat
        ncdf_varget,ncid,2,th
        pv2=fltarr(nr,nc,nth)
        ncdf_varget,ncid,3,pv2
        ncdf_close,ncid
;
; on first day reset monthly mean arrays
;
        if n eq 0L then begin
           y3d=0.*pv2
           for k=0,nth-1 do $
               for i=0,nc-1 do y3d(*,i,k)=alat
           ii_mean=0.0*pv2
        endif
;
; summate gridpoints with anomalous PV
;
        index=where(pv2*y3d lt 0. and pv2 ne 1.00000e+12)
        if index(0) ne -1 then ii_mean(index)=ii_mean(index)+1L
    endfor          ; loop over files
    ii_mean=ii_mean/float(nfile)
    save,file=dir+'ii_'+mon(m)+'avg.sav',nr,nc,nth,alon,alat,th,ii_mean
endfor
end
