;
; store multi-year monthly mean GEOS-5 in all fields 
;
@rd_geos5_nc3_meto

loadct,38
mcolor=!p.color
mcolor=byte(!p.color)
device,decompose=0
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
spawn,'ls '+dir+'*'+stimes(0)+'nc3',allfiles
;
; extract month numbers from filenames
;
wmon=strarr(n_elements(allfiles))
for n=0L,n_elements(allfiles)-1L do begin
    result=strsplit(allfiles(n),'.',/extract)
    wmon(n)=strmid(result(6),4,2)
endfor
spawn,'/usr/bin/rm '+dir+'*avg.sav'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
for m=0,nmonths-1 do begin
;
; pull out days in this month from all years
;
    smon=string(format='(i2.2)',m+1L)
    index=where(wmon eq smon,nfile)
    ifiles=allfiles(index)
    kfile=0L
    for n=0L,nfile-1 do begin
        ifile=ifiles(n)
;
; read data
;
        rd_geos5_nc3_meto,ifile,nc,nr,nth,alon,alat,th,$
                 pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
;
; make anticyclones all -1
;
        index=where(mark2 lt 0.)
        if index(0) ne -1 then mark2(index)=-1.0*mark2(index)/mark2(index)
           
        if n eq 0 then begin
           pv_mean=0.0*pv2
           npv_mean=0.0*pv2
           p_mean=0.0*pv2
           u_mean=0.0*pv2
           v_mean=0.0*pv2
           mark_mean=0.0*pv2
           sf_mean=0.0*pv2
        endif
;
; summate
;
        index=where(abs(pv2) le 10.)
        if index(0) ne -1L then begin
           pv_mean(index)=pv_mean(index)+pv2(index)
           npv_mean(index)=npv_mean(index)+1.0
        endif
        p_mean=p_mean+p2
        u_mean=u_mean+u2
        v_mean=v_mean+v2
        mark_mean=mark_mean+mark2
        sf_mean=sf_mean+sf2
        kfile=kfile+1L
    endfor          ; loop over files
;
; divide by total number of days summed (except for PV)
;
    index=where(npv_mean gt 0.)
    pv_mean(index)=pv_mean(index)/npv_mean(index)
    p_mean=p_mean/float(kfile)
    u_mean=u_mean/float(kfile)
    v_mean=v_mean/float(kfile)
    mark_mean=mark_mean/float(kfile)
    sf_mean=sf_mean/float(kfile)

; write monthly means 
    save,file=dir+mon(m)+'avg.sav',nc,nr,nth,alon,alat,th,$
         pv_mean,p_mean,u_mean,v_mean,mark_mean,sf_mean

    thlev=20
    plt1=transpose(reform(pv_mean(*,*,thlev),nr,nc))
    plt=fltarr(nc+1,nr)
    plt(0:nc-1,0:nr-1)=plt1
    plt(nc,*)=plt(0,*)
    mplt1=transpose(reform(mark_mean(*,*,thlev),nr,nc))
    mplt=fltarr(nc+1,nr)
    mplt(0:nc-1,0:nr-1)=mplt1
    mplt(nc,*)=mplt(0,*)
    x=fltarr(nc+1)
    x(0:nc-1)=alon
    x(nc)=x(0)
    !type=2^2+2^3
    nlvls=30
    erase
    col1=3+indgen(nlvls)*mcolor/nlvls
    MAP_SET,0,0,0,/contin,/grid,/noeras,title=mon(m)
    contour,plt,x,alat,nlevels=30,/fill,/cell_fill,/overplot,c_color=col1,/noeras
    contour,plt,x,alat,nlevels=30,/follow,/overplot,c_color=0,/noeras,c_labels=0*indgen(30)
    contour,mplt,x,alat,levels=[0.1],/follow,color=0,thick=5,/noeras,/overplot
    contour,mplt,x,alat,levels=[-0.1],/follow,color=mcolor,thick=5,/noeras,/overplot
    MAP_SET,0,0,0,/contin,/noeras
endfor  ; loop over months
end
