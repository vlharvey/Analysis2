;
; store multi-year monthly mean SD-WACCM in T, U, Mark
; /Volumes/cloud/data/WACCM_data/Datfiles_SD_New/f.e11.FWTREFC1SD.f19.f19.ccmi30.001.cam.h5.????????.nc3
;
@rd_sdwaccm4_nc3_new

loadct,39
mcolor=!p.color
mcolor=byte(!p.color)
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.1
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif

dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD_New/f.e11.FWTREFC1SD.f19.f19.ccmi30.001.cam.h5.'
allfiles1=file_search(dir+'199*.nc3')
allfiles2=file_search(dir+'2*.nc3')
allfiles=[allfiles1,allfiles2]
;
; extract month numbers from filenames
;
nday=n_elements(allfiles)
wmon=strarr(nday)
for n=0L,nday-1L do begin
    result=strsplit(allfiles(n),'.',/extract)
    wmon(n)=strmid(result(-2),4,2)
endfor

;spawn,'/usr/bin/rm '+dir+'*avg.sav'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
for m=0,nmonths-1 do begin
;
; pull out days in this month from all years
;
    smon=string(format='(i2.2)',m+1L)
    index=where(wmon eq smon,nfile)	; all january files
    ifiles=allfiles(index)		; ifiles=all january days
    kfile=0L
    print,smon

    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,font_size=9
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='sdwaccm_Ubar_markbar_'+smon+'.ps'
       !p.charsize=2
       !p.thick=2
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

    for n=0L,nfile-1 do begin
        ifile=ifiles(n)
        print,ifile
;
; read data
;
;       rd_sdwaccm4_nc3_new,ifile,nc,nr,nth,alon,alat,th,$
;          pv2,p2,z2,u2,v2,q2,qdf2,mark2,sf2,h2o2,n2o2,o32,co22,no2,no22,iflag
        ncid=ncdf_open(ifile)
        result0=ncdf_inquire(ncid)
        for idim=0,result0.ndims-1 do begin
            ncdf_diminq,ncid,idim,name,dim
            if name eq 'number_of_latitudes' then nr=dim
            if name eq 'number_of_longitudes' then nc=dim
            if name eq 'number_of_levels' then nth=dim
;           print,'read ',name,' dimension ',dim
        endfor
        for ivar=0,result0.nvars-1 do begin
            result=ncdf_varinq(ncid,ivar)
            if result.name eq 'latitude' or result.name eq 'longitude' or result.name eq 'theta' or result.name eq 'P' or $
               result.name eq 'GPH' or result.name eq 'U' or result.name eq 'MARK' then ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
            if result.name eq 'latitude' then alat=data             ; latitude  (deg)
            if result.name eq 'longitude' then alon=data    ; longitude (deg E)
            if result.name eq 'theta' then th=data           ; potential temperature (K)
            if result.name eq 'P' then p2=data              ; pressure  (hPa)
            if result.name eq 'U' then u2=data              ; zonal wind (m/s)
            if result.name eq 'GPH' then z2=data            ; geopotential height (m)
            if result.name eq 'MARK' then mark2=data              ; vortex marker
;           print,ivar,result.name,min(data),max(data)
        endfor
        ncdf_close,ncid
;
; make anticyclones all -1
;
        index=where(mark2 lt 0.)
        if index(0) ne -1 then mark2(index)=-1.0*mark2(index)/mark2(index)
           
        if n eq 0 then begin
           p_mean=0.0*p2
           z_mean=0.0*p2
           u_mean=0.0*p2
           v_mean=0.0*p2
           mark_mean=0.0*p2
        endif
;
; summate
;
        p_mean=p_mean+p2
        z_mean=z_mean+z2
        u_mean=u_mean+u2
        mark_mean=mark_mean+mark2
        kfile=kfile+1L
    endfor          ; loop over files
;
; divide by total number of days summed
;
    p_mean=p_mean/float(kfile)
    z_mean=z_mean/float(kfile)
    u_mean=u_mean/float(kfile)
    mark_mean=mark_mean/float(kfile)

; write monthly means 
    save,file=dir+mon(m)+'avg.sav',nc,nr,nth,alon,alat,th,$
         p_mean,z_mean,u_mean,mark_mean
;
; plot multi-year monthly mean zonal mean
;
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    plt=mean(u_mean,dim=2)
    mplt=mean(mark_mean,dim=2)
    zplt=mean(z_mean,dim=2)
    nlvls=21
    ulevel=-100+10.*findgen(nlvls)
    col1=3+indgen(nlvls)*mcolor/nlvls
    contour,plt,alat,zplt/1000.,levels=ulevel,/fill,/cell_fill,c_color=col1,/noeras,title=strupcase(strmid(mon(m),0,3)),xticks=6,xtitle='Latitude',ytitle='Altitude (km)',color=0,yrange=[10,110]
    index=where(ulevel gt 0.)
    contour,plt,alat,zplt/1000.,levels=ulevel(index),/follow,/overplot,c_color=0,/noeras
    index=where(ulevel lt 0.)
    contour,plt,alat,zplt/1000.,levels=ulevel(index),/follow,/overplot,c_color=mcolor,c_linestyle=5,/noeras
    contour,mplt,alat,zplt/1000.,levels=0.1+0.1*findgen(9),/follow,color=0,thick=10,/noeras,/overplot
    contour,mplt,alat,zplt/1000.,levels=-0.9+0.1*findgen(9),/follow,color=mcolor,thick=10,/noeras,/overplot
    imin=min(ulevel)
    imax=max(ulevel)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle='SDWACCM Zonal Mean Wind Speed (m/s)'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim sdwaccm_Ubar_markbar_'+smon+'.ps -rotate -90 sdwaccm_Ubar_markbar_'+smon+'.png'
    endif
endfor  ; loop over months
end
