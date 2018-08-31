;
; read MERRA2 netcdf data and interpolate to WACCM grid
; -- radiation diagnostics
;
longitude_waccm=2.5*findgen(144)
latitude_waccm=-90.+1.89474*findgen(96)
;levels_merra=[0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,10.,$
;              15.,20.,30.,40.,50.,60.,70.,80.,100.,150.,200.,300.,400.,500.]
;levels_rad=[0.1, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10., 20., 30., 40., 50., 70., 100., 150., 200., 250., 300.,$
;            350., 400., 450., 500., 550., 600., 650., 700., 725., 750., 775., 800., 825., 850., 875., 900., 950., 975., 1000.]
nc_waccm=n_elements(longitude_waccm)
nr_waccm=n_elements(latitude_waccm)
latitude_waccm(nr_waccm-1)=90.
print,'starting '
spawn,'date'
;
; /atmos/harvey/MERRA2_data/Datfiles/MERRA2_100.tavg3_3d_rad_Np.19910701.nc4	(100 may also be 200, 300, 400)
;
dir='/atmos/harvey/MERRA2_data/Datfiles/'
pre2='MERRA2*.tavg3_3d_rad_Np.????????'
dum2=file_search(dir+pre2+'.nc4',count=count2)
for ifile=0L,count2-1L do begin
;
; extract date
;
    result=strsplit(dum2(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    sdate=result2(-2)
    print,sdate
    dum=file_search(dir+'MERRA2-on-WACCM_rad_'+sdate+'00.sav')	; if 00Z is there then skip day
    if dum(0) ne '' then goto,jumpday
;
; read MERRA radiation products
;
     ifilename=dum2(ifile)	;strmid(dum2(ifile),0,85)
     print,ifilename
     ncid=ncdf_open(ifilename)
     result0=ncdf_inquire(ncid)
     for idim=0,result0.ndims-1 do begin
         ncdf_diminq,ncid,idim,name,dim
         if name eq 'lon' then nc_rad=dim
         if name eq 'lat' then nr_rad=dim
         if name eq 'lev' then nl_rad=dim
         if name eq 'time' then nt_rad=dim
         print,'read ',name,' dimension ',dim
     endfor
     for ivar=0,result0.nvars-1 do begin
         result=ncdf_varinq(ncid,ivar)
         ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
 
         if result.name eq 'lon' then longitude_rad=data
         if result.name eq 'lat' then latitude_rad=data
         if result.name eq 'lev' then levels_rad=data
         if result.name eq 'time' then time=data
         if result.name eq 'DTDTLWR' then dtdtlwr_all=data		; T tendency from terrestrial radiation
         if result.name eq 'DTDTLWRCLR' then dtdtlwrclr_all=data	; T tendency from terrestrial radiation (clear sky)
         if result.name eq 'DTDTSWR' then dtdtswr_all=data		; T tendency from solar radiation
         if result.name eq 'DTDTSWRCLR' then dtdtswrclr_all=data	; T tendency from solar radiation (clear sky)
         print,ivar,result.name,min(data),max(data)
     endfor
     ncdf_close,ncid
;
; only retain 0, 360, 720, 1080 times to match Nv files
;
     index=where(time eq 0L or time eq 360L or time eq 720L or time eq 1080L,nt_rad)
     time=time(index)
     dtdtlwr_all=dtdtlwr_all(*,*,*,index)
     dtdtlwrclr_all=dtdtlwrclr_all(*,*,*,index)
     dtdtswr_all=dtdtswr_all(*,*,*,index)
     dtdtswrclr_all=dtdtswrclr_all(*,*,*,index)
;
;  rearrange longitudes from -180 to 180 to 0 to 360
;
     index=where(longitude_rad lt 0.)
     if index(0) ne -1L then longitude_rad(index)=longitude_rad(index)+360.
     index=sort(longitude_rad)
     longitude_rad=longitude_rad(index)
     dtdtlwr_all=dtdtlwr_all(index,*,*,*)
     dtdtlwrclr_all=dtdtlwrclr_all(index,*,*,*)
     dtdtswr_all=dtdtswr_all(index,*,*,*)
     dtdtswrclr_all=dtdtswrclr_all(index,*,*,*)
;
; reverse levels to be top down
;
     index=sort(levels_rad)
     levels_rad=levels_rad(index)
     dtdtlwr_all=dtdtlwr_all(*,*,index,*)
     dtdtlwrclr_all=dtdtlwrclr_all(*,*,index,*)
     dtdtswr_all=dtdtswr_all(*,*,index,*)
     dtdtswrclr_all=dtdtswrclr_all(*,*,index,*)
;
; loop over times/day
;
     for itime=0L,nt_rad-1L do begin

     stime=string(format='(i2.2)',time(itime)/60L)
     dtdtlwr=reform(dtdtlwr_all(*,*,*,itime))
     dtdtlwrclr=reform(dtdtlwrclr_all(*,*,*,itime))
     dtdtswr=reform(dtdtswr_all(*,*,*,itime))
     dtdtswrclr=reform(dtdtswrclr_all(*,*,*,itime))
;
; sum to get net radiation
;
     dtdtnet=dtdtlwr+dtdtswr+dtdtswrclr+dtdtlwrclr
     bad=where(dtdtlwr eq 1.00000e+15 or dtdtswr eq 1.00000e+15 or dtdtlwrclr eq 1.00000e+15 or dtdtswrclr eq 1.00000e+15)
     if bad(0) ne -1L then dtdtnet(bad)=0./0.
     good=where(finite(dtdtnet) eq 1)
     if good(0) ne -1L then dtdtnet(good)=dtdtnet(good)*86400.	; K/day
;
; interpolate to WACCM horizontal grid
; retain original pressure grid
;
    qnew=fltarr(nc_waccm,nr_waccm,nl_rad)
    for ii=0L,nc_waccm-1L do begin
        xp=longitude_waccm(ii)
        if xp gt max(longitude_rad) then xp=xp-360.
        if xp lt min(longitude_rad) then xp=xp+360.
        for jj=0L,nr_waccm-1L do begin
            yp=latitude_waccm(jj)

            for i=0L,nc_rad-1L do begin
                ip1=i+1
                if i eq nc_rad-1L then ip1=0
                xlon=longitude_rad(i)
                xlonp1=longitude_rad(ip1)
                if i eq nc_rad-1L then xlonp1=360.+longitude_rad(ip1)
                if xp ge xlon and xp le xlonp1 then begin
                   xscale=(xp-xlon)/(xlonp1-xlon)
                   for j=0L,nr_rad-1L do begin
                       jm1=j
                       jp1=j+1
                       if j eq nr_rad-1L then begin
                          jm1=j-1
                          jp1=j
                       endif
                       xlat=latitude_rad(jm1)
                       xlatp1=latitude_rad(jp1)

                       if yp ge xlat and yp le xlatp1 then begin
                          yscale=(yp-xlat)/(xlatp1-xlat)

                          for k=0L,nl_rad-1L do begin
                          if finite(dtdtnet(i,jm1,k)) eq 1 and finite(dtdtnet(ip1,jm1,k)) eq 1 and finite(dtdtnet(i,jp1,k)) eq 1 and finite(dtdtnet(ip1,jp1,k)) eq 1 then begin
                          qj1=dtdtnet(i,jm1,k)+xscale*(dtdtnet(ip1,jm1,k)-dtdtnet(i,jm1,k))
                          qjp1=dtdtnet(i,jp1,k)+xscale*(dtdtnet(ip1,jp1,k)-dtdtnet(i,jp1,k))
                          qnew(ii,jj,k)=qj1+yscale*(qjp1-qj1)
                          endif
                          endfor

                          goto,jumplev1
                       endif
                   endfor
                endif
            endfor
jumplev1:
 
        endfor        ; loop over WACCM latitudes
    endfor            ; loop over WACCM longitudes
;
; test
;
;erase
;loadct,39
;!type=2^2+2^3
;device,decompose=0
;set_viewport,.1,.9,.2,.8
;qbar=fltarr(nr_rad,nl_rad)
;for j=0,nr_rad-1L do $
;    for k=0,nl_rad-1L do qbar(j,k)=mean(dtdtnet(*,j,k))
;contour,qbar,latitude_rad,levels_rad,/ylog,xrange=[-90,90],yrange=[1000.,0.01],levels=-60.+3.*findgen(20),thick=3,c_linestyle=5,title=sdate+stime
;contour,qbar,latitude_rad,levels_rad,/overplot,levels=1.+findgen(19),thick=3
;qbar_new=fltarr(nr_waccm,nl_rad)
;for j=0,nr_waccm-1L do $
;    for k=0,nl_rad-1L do qbar_new(j,k)=mean(qnew(*,j,k))
;contour,qbar_new,latitude_waccm,levels_rad,/overplot,levels=-60.+3.*findgen(20),thick=3,c_linestyle=5,color=250
;contour,qbar_new,latitude_waccm,levels_rad,/overplot,levels=1.+findgen(19),thick=3,color=250
;stop
;erase
;xyouts,.4,.8,sdate+stime,/normal,charsize=2
;set_viewport,.1,.45,.2,.6
;contour,reform(dtdtnet(*,*,0)),longitude_rad,latitude_rad,/noeras,levels=-60.+4.*findgen(20),/cell,c_color=1+findgen(20)*255/20.
;contour,reform(dtdtnet(*,*,0)),longitude_rad,latitude_rad,/overplot,levels=-60.+4.*findgen(20),/foll,color=0
;contour,reform(dtdtnet(*,*,0)),longitude_rad,latitude_rad,/overplot,levels=0,/foll,color=0,thick=10
;
;set_viewport,.55,.9,.2,.6
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/noeras,levels=-60.+4.*findgen(20),/cell,c_color=1+findgen(20)*255/20.
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/overplot,levels=-60.+4.*findgen(20),/foll,color=0
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/overplot,levels=0,/foll,color=0,thick=10
;
;set_viewport,.1,.9,.2,.8
;map_set,0,180,0,/contin,/grid,color=255,/noeras,title=sdate+stime+'  '+strcompress(levels_rad(0))+' hPa'
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/noeras,levels=-60.+4.*findgen(20),/cell,c_color=1+findgen(20)*255/20.,/overplot
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/overplot,levels=-60.+4.*findgen(20),/foll,color=0
;contour,reform(qnew(*,*,0)),longitude_waccm,latitude_waccm,/overplot,levels=0,/foll,color=0,thick=10
;map_set,0,180,0,/contin,/grid,color=255,/noeras
;xnoon=(360./24.)*(12.-float(stime))		; longitude of local noon is 360/24*(LT-UT) where LT is noon
;if xnoon gt 360. then xnoon=xnoon-360
;if xnoon lt 0. then xnoon=xnoon+360
;oplot,xnoon+0.*findgen(181),-90.+findgen(181),thick=5,color=255
;
; IDL save file for each day
;
    ofile=dir+'MERRA2-on-WACCM_rad_'+sdate+stime+'.sav'
    print,ofile
    save,file=ofile,longitude_waccm,latitude_waccm,levels_rad,qnew

    endfor	; loop over times/day
;
; remove original file
;
;   spawn,'rm -f '+ifilename
jumpday:
endfor	; loop over daily files
print,'ending '
spawn,'date'
end
