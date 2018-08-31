;
; save multi-year monthly means of CNRM Surface Area Density of sulf
; /aura7/harvey/CCMval_data/Datfiles
;
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
for imon=0L,nmon-1L do begin

    dir='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_sulf_????'+smon(imon)+'??.sav'
    spawn,'ls '+dir,ncfiles
    nfile=n_elements(ncfiles)
    for ifile=0L,nfile-1L do begin
        ncfile=ncfiles(ifile)
        restore,ncfile	;,alon,alat,lev,sadgrd
        result=strsplit(ncfile,'_',/extract)
        result2=strsplit(result(8),'.',/extract)
        sdate=result2(0)
        print,'opening '+ncfile
        if ifile eq 0L then begin
           wsadgrd=0.*sadgrd
           nwsadgrd=long(0.*sadgrd)
        endif
        index=where(sadgrd ne 0.)
        if index(0) ne -1L then begin
           wsadgrd(index)=wsadgrd(index)+sadgrd(index)
           nwsadgrd(index)=nwsadgrd(index)+1L
        endif
    endfor	; loop over files

    index=where(nwsadgrd gt 0L)
    if index(0) ne -1L then begin
       wsadgrd(index)=wsadgrd(index)/float(nwsadgrd(index))
       erase
       plot,wsadgrd(index),psym=3,title=sdate
    endif
;
; save monthly mean
;
    save,file='/aura7/harvey/CCMval_data/Datfiles/CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_sulf_'+smon(imon)+'.sav',$
          alon,alat,lev,wsadgrd,nwsadgrd
    print,'saved CCMVal2_REF-B1_CNRM-ACM_2_T3I_sad_sulf_'+smon(imon)+'.sav'
endfor
end
