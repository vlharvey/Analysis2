;
; read MERRA netcdf data and interpolate to WACCM grid
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
; /Volumes/Data/MERRA_data/Datfiles/MERRA301.prod.assim.inst6_3d_ana_Nv.20111231.SUB.nc
; /Volumes/Data/MERRA_data/Datfiles/MERRA300.prod.assim.tavg3_3d_rad_Cp.20111231.SUB.nc
;
dir='/Volumes/Data/MERRA_data/Datfiles/'
pre1='MERRA*.prod.assim.inst6_3d_ana_Nv.'

;pre2='MERRA*.prod.assim.tavg3_3d_rad_Cp.'
dum1=file_search(dir,pre1+'*SUB.nc.gz',count=count1)
;dum2=file_search(dir,pre2+'*.nc',count=count2)
;if count1 ne count2 then stop,'check number of files'
for ifile=0L,count1-1L do begin
;
; extract date
;
    result=strsplit(dum1(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    sdate=result2(4)
    print,sdate
;
; read MERRA radiation products
;
;    print,dum2(ifile)
;    ncid=ncdf_open(dum2(ifile))
;    result0=ncdf_inquire(ncid)
;    for idim=0,result0.ndims-1 do begin
;        ncdf_diminq,ncid,idim,name,dim
;        if name eq 'longitude' then nc_rad=dim
;        if name eq 'latitude' then nr_rad=dim
;        if name eq 'levels' then nl_rad=dim
;        if name eq 'time' then nt_rad=dim
;        print,'read ',name,' dimension ',dim
;    endfor
;    for ivar=0,result0.nvars-1 do begin
;        result=ncdf_varinq(ncid,ivar)
;        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
;
;        if result.name eq 'longitude' then longitude_rad=data
;        if result.name eq 'latitude' then latitude_rad=data
;        if result.name eq 'levels' then levels_rad=data
;        if result.name eq 'dtdtlwr' then dtdtlwr=data		; T tendency from terrestrial radiation
;        if result.name eq 'dtdtswr' then dtdtswr=data		; T tendency from solar radiation
;        if result.name eq 'dtdtswrclr' then dtdtswrclr=data	; T tendency from solar radiation (clear sky)
;        print,ivar,result.name,min(data),max(data)
;    endfor
;    ncdf_close,ncid
;
; sum to get net radiation
;
;    dtdtnet=dtdtlwr+dtdtswr+dtdtswrclr
;    bad=where(dtdtlwr eq 1.00000e+15 or dtdtswr eq 1.00000e+15 or dtdtswrclr eq 1.00000e+15)
;    if bad(0) ne -1L then dtdtnet(bad)=0./0.
;    good=where(finite(dtdtnet) eq 1)
;    if good(0) ne -1L then dtdtnet(good)=dtdtnet(good)*86400.	; K/day
;
; rearrange longitudes from -180 to 180 to 0 to 360
;
;   index=where(longitude_rad lt 0.)
;    if index(0) ne -1L then longitude_rad(index)=longitude_rad(index)+360.
;    index=sort(longitude_rad)
;    longitude_rad=longitude_rad(index)
;    dtdtnet=dtdtnet(index,*,*)
;
; reverse levels to be top down
;
;    index=sort(levels_rad)
;    levels_rad=levels_rad(index)
;    dtdtnet=dtdtnet(*,*,index)
;
; read MERRA U, V, T, PS, DELP, QV
;
    spawn,'gunzip '+dum1(ifile)
    result=strsplit(dum1(ifile),'.',/extract)
    ifilename=result(0)+'.'+result(1)+'.'+result(2)+'.'+result(3)+'.'+result(4)+'.'+result(5)+'.'+result(6)
    print,ifilename
    ncid=ncdf_open(ifilename)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'longitude' then nc_assim=dim
        if name eq 'latitude' then nr_assim=dim
        if name eq 'levels' then nl_assim=dim
        if name eq 'time' then nt_assim=dim
;       print,'read ',name,' dimension ',dim
    endfor
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data

        if result.name eq 'longitude' then longitude_assim=data
        if result.name eq 'latitude' then latitude_assim=data
        if result.name eq 'levels' then levels_assim=data
        if result.name eq 'ps' then ps=data			; surface pressure
        if result.name eq 'delp' then delp=data			; layer pressure thickness
        if result.name eq 't' then t=data			; air temperature
        if result.name eq 'u' then u=data			; zonal wind
        if result.name eq 'v' then v=data			; meridional wind
;       if result.name eq 'o3' then o3=data			; ozone
        if result.name eq 'qv' then qv=data			; specific humidity
;       print,ivar,result.name,min(data),max(data)
    endfor
    ncdf_close,ncid
;
; rearrange longitudes from -180 to 180 to 0 to 360
;
    index=where(longitude_assim lt 0.)
    if index(0) ne -1L then longitude_assim(index)=longitude_assim(index)+360.
    index=sort(longitude_assim)
    longitude_assim=longitude_assim(index)
    ps=ps(index,*)/100.						; hPa
    delp=delp(index,*,*)
    t=t(index,*,*)
    u=u(index,*,*)
    v=v(index,*,*)
;   o3=o3(index,*,*)
    qv=qv(index,*,*)
;
; interpolate to WACCM horizontal grid
; retain original pressure grid
;
    psnew=fltarr(nc_waccm,nr_waccm)
    dpnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    tnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    unew=fltarr(nc_waccm,nr_waccm,nl_assim)
    vnew=fltarr(nc_waccm,nr_waccm,nl_assim)
;   o3new=fltarr(nc_waccm,nr_waccm,nl_assim)
    qvnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    qnew=fltarr(nc_waccm,nr_waccm,nl_assim)
    for ii=0L,nc_waccm-1L do begin
        xp=longitude_waccm(ii)
        if xp gt max(longitude_waccm) then xp=xp-360.
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

;                         o3j1=o3(i,jm1,*)+xscale*(o3(ip1,jm1,*)-o3(i,jm1,*))
;                         o3jp1=o3(i,jp1,*)+xscale*(o3(ip1,jp1,*)-o3(i,jp1,*))
;                         o3new(ii,jj,*)=o3j1+yscale*(o3jp1-o3j1)

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
;;
;; interpolate q
;;
;    for ii=0L,nc_waccm-1L do begin
;        xp=longitude_waccm(ii)
;        if xp gt max(longitude_waccm) then xp=xp-360.
;        for jj=0L,nr_waccm-1L do begin
;            yp=latitude_waccm(jj)
;
;            for ir=0L,nc_rad-1L do begin
;                im1=ir
;                ip1=ir+1
;                if ir eq nc_rad-1L then ip1=0
;                xlon=longitude_rad(im1)
;                xlonp1=longitude_rad(ip1)
;                if ir eq nc_rad-1L then xlonp1=360.+longitude_rad(ip1)
;                if xp ge xlon and xp le xlonp1 then begin
;                   xscale=(xp-xlon)/(xlonp1-xlon)
;                   for jr=0L,nr_rad-1L do begin
;                       jm1=jr
;                       jp1=jr+1
;                       if jr eq nr_rad-1L then begin
;                          jm1=jr-1
;                          jp1=jr
;                       endif
;                       xlat=latitude_rad(jm1)
;                       xlatp1=latitude_rad(jp1)
;                       if yp ge xlat and yp le xlatp1 then begin
;                          yscale=(yp-xlat)/(xlatp1-xlat)
;
;                          for kk=1L,nl_waccm-1L do begin
;                              zp=levels_waccm(kk)
;                              for k=0L,nl_rad-15L do begin			; NaN below 600 hPa
;                                  kp1=k+1
;                                  z0=levels_rad(k)
;                                  z1=levels_rad(k+1)
;                                  if zp ge z0 and zp le z1 then begin
;                                     zscale=(zp-z0)/(z1-z0)
;                                     qj1=dtdtnet(im1,jm1,k)+xscale*(dtdtnet(ip1,jm1,k)-dtdtnet(im1,jm1,k))
;                                     qjp1=dtdtnet(im1,jp1,k)+xscale*(dtdtnet(ip1,jp1,k)-dtdtnet(im1,jp1,k))
;                                     qj2=dtdtnet(im1,jm1,kp1)+xscale*(dtdtnet(ip1,jm1,kp1)-dtdtnet(im1,jm1,kp1))
;                                     qjp2=dtdtnet(im1,jp1,kp1)+xscale*(dtdtnet(ip1,jp1,kp1)-dtdtnet(im1,jp1,kp1))
;                                     q1=qj1+yscale*(qjp1-qj1)
;                                     q2=qj2+yscale*(qjp2-qj2)
;                                     qnew(ii,jj,kk)=q1+zscale*(q2-q1)
;;print,z0,zp,z1,zscale,latitude_waccm(jj)
;;print,q1,qnew(ii,jj,kk),q2
;                                     goto,jumplev2
;                                  endif
;jumplev2:
;                              endfor
;                          endfor
;                       endif
;                   endfor
;                endif
;            endfor
;;
;; radiation lats do not extend to poles
;;
;            qnew(*,0,*)=qnew(*,1,*)
;            qnew(*,nr_waccm-1,*)=qnew(*,nr_waccm-2,*)
;
;        endfor        ; loop over WACCM latitudes
;    endfor            ; loop over WACCM longitudes
;
; test
;
   loadct,39
   device,decompose=0
   tbar=fltarr(nr_assim,nl_assim)
;  qbar=fltarr(nr_rad,nl_rad)
   for j=0,nr_assim-1L do $
       for k=0,nl_assim-1L do tbar(j,k)=mean(t(*,j,k))
;  for j=0,nr_rad-1L do $
;      for k=0,nl_rad-1L do qbar(j,k)=mean(dtdtnet(*,j,k))
   contour,tbar,latitude_assim,levels_assim,/ylog,yrange=[1000.,0.01],levels=150+10.*findgen(30),thick=3
;  contour,ps,longitude_assim,latitude_assim,levels=500+10.*findgen(50)
;  contour,qbar,latitude_rad,levels_rad,/ylog,yrange=[1000.,0.01],levels=-20.+findgen(20),thick=3,c_linestyle=5
;  contour,qbar,latitude_rad,levels_rad,/overplot,levels=1.+findgen(19),thick=3
   tbar_new=fltarr(nr_waccm,nl_assim)
;  qbar_new=fltarr(nr_waccm,nl_assim)
   for j=0,nr_waccm-1L do $
       for k=0,nl_assim-1L do tbar_new(j,k)=mean(tnew(*,j,k))
;  for j=0,nr_waccm-1L do $
;      for k=0,nl_assim-1L do qbar_new(j,k)=mean(qnew(*,j,k))
   contour,tbar_new,latitude_waccm,levels_assim,/overplot,levels=150+10.*findgen(30),thick=3,color=250
;  contour,psnew,longitude_waccm,latitude_waccm,/overplot,levels=500+10.*findgen(50),color=250
;  contour,qbar_new,latitude_waccm,levels_assim,/overplot,levels=-20.+findgen(20),thick=3,c_linestyle=5,color=250
;  contour,qbar_new,latitude_waccm,levels_assim,/overplot,levels=1.+findgen(19),thick=3,color=250
;
; IDL save file for each day
;
    ofile=dir+'MERRA-on-WACCM_'+sdate+'.sav'
    print,ofile
    save,file=ofile,longitude_waccm,latitude_waccm,levels_assim,tnew,unew,vnew,psnew,dpnew,qvnew		;,qnew
;
; remove original file
;
    spawn,'rm -f '+ifilename
endfor	; loop over daily files
end
