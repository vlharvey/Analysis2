;
; read MERRA2 netcdf data and interpolate to WACCM grid
; include Specific Humidity and delP
;
longitude_waccm=2.5*findgen(144)
latitude_waccm=-90.+1.89474*findgen(96)
levels_waccm=[0.015, 0.0263, 0.0401, 0.0567, 0.0776, 0.1045, 0.1395, 0.1854,$
    0.2449, 0.3217, 0.4204, 0.5462, 0.7059, 0.9072, 1.1599, 1.4756, 1.8678, 2.3525,$
    2.9483, 3.6765, 4.5616, 5.6318, 6.9183, 8.4563, 10.2849, 12.4601,$
    15.0502, 18.1243, 21.7610, 26.0491, 31.0889, 36.9927, 43.9096, 52.0159,$
    61.4956, 72.5578, 85.4390, 100.514, 118.25, 139.115, 163.661, 192.541,$
    226.513, 266.479, 312.791, 356.250, 393.750, 431.250, 468.750, 506.250,$
    543.750, 581.250, 618.750, 656.250, 687.500, 712.500, 737.500, 762.500,$
    787.500, 810.000, 827.500, 842.500, 857.500, 872.500, 887.500, 902.500,$
    917.500, 932.500, 947.500, 962.500, 977.500, 992.500]
nc_waccm=n_elements(longitude_waccm)
nr_waccm=n_elements(latitude_waccm)
latitude_waccm(nr_waccm-1)=90.
nl_waccm=n_elements(levels_waccm)
;
; /atmos/harvey/MERRA2_data/Datfiles/MERRA2_200.inst6_3d_ana_Nv.19920801.nc4
;
print,'starting'
spawn,'date'
dir='/atmos/harvey/MERRA2_data/Datfiles/'
pre1='MERRA2*.inst6_3d_ana_Nv.????????'

dum1=file_search(dir+pre1+'.nc4',count=count1)
for ifile=0L,count1-1L do begin
;
; extract date
;
    result=strsplit(dum1(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    sdate=result2(-2)
    print,sdate
;
; check for processed data
;
    dum=file_search(dir+'MERRA2-on-WACCM_'+sdate+'00.sav')  ; if 00Z is there then skip day
    if dum(0) ne '' then goto,jumpday
;
; read MERRA U, V, T, PS, DELP, QV
;
    ifilename=dum1(ifile)
    print,ifilename
    ncid=ncdf_open(ifilename)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'lon' then nc_assim=dim
        if name eq 'lat' then nr_assim=dim
        if name eq 'lev' then nl_assim=dim
        if name eq 'time' then nt_assim=dim
        print,'read ',name,' dimension ',dim
    endfor
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data

        if result.name eq 'lon' then longitude_assim=data
        if result.name eq 'lat' then latitude_assim=data
        if result.name eq 'lev' then levels_assim=data		; model levels range from 1 to 72
        if result.name eq 'time' then time=data
        if result.name eq 'PS' then ps_all=data			; surface pressure
        if result.name eq 'DELP' then delp_all=data		; layer pressure thickness
        if result.name eq 'T' then t_all=data			; air temperature
        if result.name eq 'U' then u_all=data			; zonal wind
        if result.name eq 'V' then v_all=data			; meridional wind
        if result.name eq 'O3' then o3_all=data			; ozone
        if result.name eq 'QV' then qv_all=data			; specific humidity
        print,ivar,result.name,min(data),max(data)
    endfor
    ncdf_close,ncid
;
; rearrange longitudes from -180 to 180 to 0 to 360
;
    index=where(longitude_assim lt 0.)
    if index(0) ne -1L then longitude_assim(index)=longitude_assim(index)+360.
    index=sort(longitude_assim)
    longitude_assim=longitude_assim(index)
    ps_all=ps_all(index,*,*)/100.						; hPa
    delp_all=delp_all(index,*,*,*)
    t_all=t_all(index,*,*,*)
    u_all=u_all(index,*,*,*)
    v_all=v_all(index,*,*,*)
    o3_all=o3_all(index,*,*,*)
    qv_all=qv_all(index,*,*,*)
;
; loop over 4x per day
;
    for itime=0L,nt_assim-1L do begin

    stime=string(format='(i2.2)',time(itime)/60L)
    ps=reform(ps_all(*,*,itime))
    delp=reform(delp_all(*,*,*,itime))
    t=reform(t_all(*,*,*,itime))
    u=reform(u_all(*,*,*,itime))
    v=reform(v_all(*,*,*,itime))
    o3=reform(o3_all(*,*,*,itime))
    qv=reform(qv_all(*,*,*,itime))
;
; interpolate to WACCM horizontal grid
; retain original model levels
;
    psnew=fltarr(nc_waccm,nr_waccm)
    dpnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    tnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    unew=fltarr(nc_waccm,nr_waccm,nl_assim)
    vnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    o3new=fltarr(nc_waccm,nr_waccm,nl_assim)
    qvnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    qnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    for ii=0L,nc_waccm-1L do begin
        xp=longitude_waccm(ii)
        if xp gt max(longitude_assim) then xp=xp-360.
        if xp lt min(longitude_assim) then xp=xp+360.
        for jj=0L,nr_waccm-1L do begin
            yp=latitude_waccm(jj)

            for i=0L,nc_assim-1L do begin
                ip1=i+1
                if i eq nc_assim-1L then ip1=0
                xlon=longitude_assim(i)
                xlonp1=longitude_assim(ip1)
                if i eq nc_assim-1L then xlonp1=360.+longitude_assim(ip1)
                if xp ge xlon and xp le xlonp1 then begin
                   xscale=(xp-xlon)/(xlonp1-xlon)
                   for j=0L,nr_assim-1L do begin
                       jm1=j
                       jp1=j+1
                       if j eq nr_assim-1L then begin
                          jm1=j-1
                          jp1=j
                       endif
                       xlat=latitude_assim(jm1)
                       xlatp1=latitude_assim(jp1)

                       if yp ge xlat and yp le xlatp1 then begin
                          yscale=(yp-xlat)/(xlatp1-xlat)

                          pj1=ps(i,jm1)+xscale*(ps(ip1,jm1)-ps(i,jm1))
                          pjp1=ps(i,jp1)+xscale*(ps(ip1,jp1)-ps(i,jp1))
                          psnew(ii,jj)=pj1+yscale*(pjp1-pj1)

                          dj1=delp(i,jm1,*)+xscale*(delp(ip1,jm1,*)-delp(i,jm1,*))
                          djp1=delp(i,jp1,*)+xscale*(delp(ip1,jp1,*)-delp(i,jp1,*))
                          dpnew(ii,jj,*)=dj1+yscale*(djp1-dj1)

                          tj1=t(i,jm1,*)+xscale*(t(ip1,jm1,*)-t(i,jm1,*))
                          tjp1=t(i,jp1,*)+xscale*(t(ip1,jp1,*)-t(i,jp1,*))
                          tnew(ii,jj,*)=tj1+yscale*(tjp1-tj1)

                          uj1=u(i,jm1,*)+xscale*(u(ip1,jm1,*)-u(i,jm1,*))
                          ujp1=u(i,jp1,*)+xscale*(u(ip1,jp1,*)-u(i,jp1,*))
                          unew(ii,jj,*)=uj1+yscale*(ujp1-uj1)

                          vj1=v(i,jm1,*)+xscale*(v(ip1,jm1,*)-v(i,jm1,*))
                          vjp1=v(i,jp1,*)+xscale*(v(ip1,jp1,*)-v(i,jp1,*))
                          vnew(ii,jj,*)=vj1+yscale*(vjp1-vj1)

                          o3j1=o3(i,jm1,*)+xscale*(o3(ip1,jm1,*)-o3(i,jm1,*))
                          o3jp1=o3(i,jp1,*)+xscale*(o3(ip1,jp1,*)-o3(i,jp1,*))
                          o3new(ii,jj,*)=o3j1+yscale*(o3jp1-o3j1)

                          qvj1=qv(i,jm1,*)+xscale*(qv(ip1,jm1,*)-qv(i,jm1,*))
                          qvjp1=qv(i,jp1,*)+xscale*(qv(ip1,jp1,*)-qv(i,jp1,*))
                          qvnew(ii,jj,*)=qvj1+yscale*(qvjp1-qvj1)

                          goto,jumplev1
                       endif
                   endfor
                endif
            endfor
jumplev1:
        endfor        ; loop over WACCM latitudes
    endfor            ; loop over WACCM longitudes
;
; plot 
;
;   erase
;   loadct,39
;   device,decompose=0
;;   tbar=mean(t,dim=1)
;;   contour,tbar,latitude_assim,levels_waccm,/ylog,yrange=[1000.,0.01],levels=150+10.*findgen(30),thick=3
;;   tbar_new=mean(tnew,dim=1)
;;   contour,tbar_new,latitude_waccm,levels_waccm,/overplot,levels=150+10.*findgen(30),thick=3,color=250

;ilev=0
;xyouts,.4,.8,sdate+stime+'  '+strcompress(levels_waccm(ilev))+' hPa',/normal,charsize=2
;set_viewport,.1,.45,.2,.6
;map_set,0,180,0,/contin,/grid,color=255,/noeras,title='Original Grid'
;;contour,reform(ps),longitude_assim,latitude_assim,/noeras,levels=520.+20.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;;contour,reform(ps),longitude_assim,latitude_assim,/overplot,levels=520.+20.*findgen(20),/foll,color=0
;;contour,reform(delp(*,*,ilev)),longitude_assim,latitude_assim,/noeras,levels=500.+20.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;;contour,reform(delp(*,*,ilev)),longitude_assim,latitude_assim,/overplot,levels=500.+20.*findgen(20),/foll,color=0
;contour,reform(t(*,*,ilev)),longitude_assim,latitude_assim,/noeras,levels=160.+5.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;contour,reform(t(*,*,ilev)),longitude_assim,latitude_assim,/overplot,levels=160.+5.*findgen(20),/foll,color=0
;map_set,0,180,0,/contin,/grid,color=255,/noeras
;
;set_viewport,.55,.9,.2,.6
;map_set,0,180,0,/contin,/grid,color=255,/noeras,title='WACCM Grid'
;;contour,reform(psnew),longitude_waccm,latitude_waccm,/noeras,levels=520.+20.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;;contour,reform(psnew),longitude_waccm,latitude_waccm,/overplot,levels=520.+20.*findgen(20),/foll,color=0
;;contour,reform(dpnew(*,*,ilev)),longitude_waccm,latitude_waccm,/noeras,levels=500.+20.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;;contour,reform(dpnew(*,*,ilev)),longitude_waccm,latitude_waccm,/overplot,levels=500.+20.*findgen(20),/foll,color=0
;contour,reform(tnew(*,*,ilev)),longitude_waccm,latitude_waccm,/noeras,levels=160.+5.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;contour,reform(tnew(*,*,ilev)),longitude_waccm,latitude_waccm,/overplot,levels=160.+5.*findgen(20),/foll,color=0
;map_set,0,180,0,/contin,/grid,color=255,/noeras
;
; IDL save file for each day
; levels_assim are integer 1-72 model levels whereas in MERRA they were approximate pressure values
;
    ofile=dir+'MERRA2-on-WACCM_'+sdate+stime+'.sav'
    print,ofile
    save,file=ofile,longitude_waccm,latitude_waccm,levels_waccm,levels_assim,tnew,unew,vnew,psnew,dpnew,qvnew,o3new
    endfor      ; loop over times/day
;
; remove original file
;
;   spawn,'rm -f '+ifilename
jumpday:
endfor	; loop over daily files
print,'ending'
spawn,'date'
end
