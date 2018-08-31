;
; ERA40 zonal mean temperature and zonal wind
; store 2-D arrays (day vs. altitude) of dT/dy and Ubar for SSWs
; for both hemispheres and Ubar at the Equator for QBO.
; save yearly sav files
;
; VLH 9/9/09
;
@stddat
@kgmt
@ckday
@kdate
@rd_era40_nc

;a=findgen(8)*(2*!pi/8.)
;usersym,cos(a),sin(a),/fill
;loadct,39
;icolmax=byte(!p.color)
;icolmax=fix(icolmax)
;if icolmax eq 0 then icolmax=255
;mcolor=icolmax
;nxdim=700
;nydim=700
;xorig=[0.2,0.2,0.2]
;yorig=[0.65,0.375,0.075]
;xlen=0.6
;ylen=0.2
;cbaryoff=0.08
;cbarydel=0.02
;device,decompose=0
;!NOERAS=-1
;setplot='ps'
;read,'setplot=',setplot
;if setplot ne 'ps' then begin
;   set_plot,'x'
;   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
;   !p.background=icolmax
;endif
dir='/aura7/harvey/ERA40_data/Datfiles/era40_ua_12Z_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
comment='dT/dy=Tbar85-Tbar60 and Ubar65 is zonal mean zonal wind at 65N/S'
lstmn=9L & lstdy=1L & lstyr=1957L
ledmn=8L & leddy=31L & ledyr=2002L
;
; Get start and end dates
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1950 or lstyr gt 2002 then stop,'Year out of range '
if ledyr lt 1950 or ledyr gt 2002 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd_all=lonarr(nfile)
syyyymmdd_all=strarr(nfile)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' Starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
      yyyymmdd_all(icount)=long(syr+smn+sdy)
      syyyymmdd_all(icount)=date
;
; read ERA40 pressure data, i.e. /aura7/harvey/ERA40_data/Datfiles/era40_ua_12Z_19600131.nc
;
      rd_era40_nc,dir+date+'.nc',nc,nr,nl,alon,alat,press,tp,uu,vv,gp,iflg
      if iflg eq 1 then goto,skip

      if kcount eq 0L then begin
         y60n=where(alat eq 60.)
         y60s=where(alat eq -60.)
         y65n=where(alat eq 65.)
         y65s=where(alat eq -65.)
         y85n=where(alat eq 85.)
         y85s=where(alat eq -85.)
         yeq=where(alat eq 0.)
; 
; WMO minor warming is T85 - T60 > 0. major is U65 < 0. Both defined at 10 hPa.
;
         nh_dtdy_all=-9999.+0.*fltarr(nfile,nl)
         sh_dtdy_all=-9999.+0.*fltarr(nfile,nl)
         nh_ubar65_all=-9999.+0.*fltarr(nfile,nl)
         sh_ubar65_all=-9999.+0.*fltarr(nfile,nl)
         ubar_eq_all=-9999.+0.*fltarr(nfile,nl)
         kcount=1
      endif
;
; calculate zonal mean temperature and zonal wind
;
      uzm=-9999.+0.*fltarr(nr,nl)
      tzm=-9999.+0.*fltarr(nr,nl)
      for k=0,nl-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tp(*,j,k))/float(nc)
              uzm(j,k)=total(uu(*,j,k))/float(nc)
          endfor
          if tzm(y85n(0),k) lt 100. or tzm(y60n(0),k) lt 100. then stop,'check temp NH'
          if tzm(y85s(0),k) lt 100. or tzm(y60s(0),k) lt 100. then stop,'check temp SH'
          nh_dtdy_all(icount,k)=tzm(y85n(0),k)-tzm(y60n(0),k)
          sh_dtdy_all(icount,k)=tzm(y85s(0),k)-tzm(y60s(0),k)
          if abs(nh_dtdy_all(icount,k)) gt 100. then stop,'check dtdy NH'
          if abs(sh_dtdy_all(icount,k)) gt 100. then stop,'check dtdy SH'

          nh_ubar65_all(icount,k)=uzm(y65n(0),k)
          sh_ubar65_all(icount,k)=uzm(y65s(0),k)
          ubar_eq_all(icount,k)=uzm(yeq(0),k)
      endfor
;
; to plot along the way, change logic to if icount gt 0
; 
;      syear=strmid(syyyymmdd,0,4)
;      smon=strmid(syyyymmdd,4,2)
;      sday=strmid(syyyymmdd,6,2)
;      if icount lt 0L then begin	; do not plot each day	;gt 1 then begin
;         erase
;         level=-40.+5.*findgen(13) ; -40 to +20
;         nlvls=n_elements(level)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         good=where(long(syear) ne 0L)
;         minyear=long(min(long(syear(good))))
;         maxyear=long(max(long(syear)))
;         yearlab=strcompress(maxyear,/remove_all)
;         if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;         if icount lt 10L then xindex=findgen(icount)
;         if icount lt 10L then nxticks=icount
;         if icount gt 10L then xindex=where(sday eq '01' or sday eq '15',nxticks)
;         if icount gt 100L then xindex=where(sday eq '15',nxticks)
;         xlabs=smon(xindex)+'/'+sday(xindex)
;         if icount gt 100L then xlabs=smon(xindex)
;         xyouts,.38,.9,yearlab,/normal,color=0,charsize=3,charthick=3
;         xmn=xorig(0)
;         xmx=xorig(0)+xlen
;         ymn=yorig(0)
;         ymx=yorig(0)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Tbar 85N - 60N',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(level lt 0.)
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         index=where(level gt 0.)
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;;        contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0],$
;;                /follow,color=0,/overplot,thick=3,min_value=-9999.
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(level)
;         imax=max(level)
;         xmnb=xorig(0)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(0)+cbarydel,yorig(0)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K',color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;         ulevel=-50.+10.*findgen(13)
;         nlvls=n_elements(ulevel)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         xmn=xorig(1)
;         xmx=xorig(1)+xlen
;         ymn=yorig(1)
;         ymx=yorig(1)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Ubar at 65N',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(ulevel gt 0.)
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         index=where(ulevel lt 0.)
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;;        contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0.],$
;;                /follow,color=0,/overplot,thick=3,min_value=-9999.
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(ulevel)
;         imax=max(ulevel)
;         xmnb=xorig(1)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(1)+cbarydel,yorig(1)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='m/s',$
;              color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;         ulevel=-50.+10.*findgen(9)
;         nlvls=n_elements(ulevel)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         xmn=xorig(2)
;         xmx=xorig(2)+xlen
;         ymn=yorig(2)
;         ymx=yorig(2)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Ubar at the Equator',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(ulevel gt 0.)
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         index=where(ulevel lt 0.)
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[5,15],$
;                 c_labels=[0,0],/follow,color=0,/overplot,min_value=-9999.
;;        contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0.],$
;;                /follow,color=0,/overplot,thick=3,min_value=-9999.
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(ulevel)
;         imax=max(ulevel)
;         xmnb=xorig(2)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(2)+cbarydel,yorig(2)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='m/s',$
;              color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;      endif
;;index=where(ubar_eq ne -9999.)
;;print,min(ubar_eq(index)),max(ubar_eq)
;;
;; save yearly plot and sav file on July 15th
;;
;      if icount gt 1L and smon(icount) eq '07' and sday(icount) eq '15' then begin
;         if icount lt 10L then xindex=findgen(icount)
;         if icount lt 10L then nxticks=icount
;         if icount gt 10L then xindex=where(sday eq '01' or sday eq '15',nxticks)
;         if icount gt 100L then xindex=where(sday eq '15',nxticks)
;         xlabs=smon(xindex)+'/'+sday(xindex)
;         if icount gt 100L then xlabs=smon(xindex)
;         good=where(long(syear) ne 0L)
;         minyear=long(min(long(syear(good))))
;         maxyear=long(max(long(syear)))
;         yearlab=strcompress(maxyear,/remove_all)
;         if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
;         set_plot,'ps'
;         xsize=nxdim/100.
;         ysize=nydim/100.
;         !p.font=0
;         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
;                /bold,/color,bits_per_pixel=8,/helvetica,filename='ERA40_dTdy_Ubar_QBO_'+yearlab+'.ps'
;         !p.charsize=1.25
;         !p.thick=2
;         !p.charthick=5
;         !p.charthick=5
;         !y.thick=2
;         !x.thick=2
;         erase
;         level=-40.+5.*findgen(13) ; -40 to +20
;         nlvls=n_elements(level)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         xyouts,.38,.9,yearlab,/normal,color=0,charsize=2,charthick=3
;         xmn=xorig(0)
;         xmx=xorig(0)+xlen
;         ymn=yorig(0)
;         ymx=yorig(0)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Tbar 85N - 60N',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(level lt 0.)
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         index=where(level gt 0.)
;         contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=level(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;;        contour,nh_dtdy(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0],$
;;                /follow,color=0,/overplot,thick=3,min_value=-9999.
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(level)
;         imax=max(level)
;         xmnb=xorig(0)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(0)+cbarydel,yorig(0)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='K',color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;         ulevel=-50.+10.*findgen(13)
;         nlvls=n_elements(ulevel)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         xmn=xorig(1)
;         xmx=xorig(1)+xlen
;         ymn=yorig(1)
;         ymx=yorig(1)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Ubar at 65N',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(ulevel gt 0.)
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         index=where(ulevel lt 0.)
;         contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;;        contour,nh_ubar65(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0.],$
;;                /follow,color=0,/overplot,thick=3,min_value=-9999.
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(ulevel)
;         imax=max(ulevel)
;         xmnb=xorig(1)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(1)+cbarydel,yorig(1)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='m/s',$
;              color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;         ulevel=-50.+10.*findgen(9)
;         nlvls=n_elements(ulevel)
;         col1=1+indgen(nlvls)*mcolor/nlvls
;         !type=2^2+2^3+2^7
;         xmn=xorig(2)
;         xmx=xorig(2)+xlen
;         ymn=yorig(2)
;         ymx=yorig(2)+ylen
;         set_viewport,xmn,xmx,ymn,ymx
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel,$
;                 c_color=col1,/cell_fill,color=0,ytitle='Pressure (hPa)',title='Ubar at the Equator',xrange=[0.,icount-1],$
;                 xticks=nxticks-1,xtickname=' '+strarr(nxticks+1),charsize=1.25,charthick=2,xticklen=-0.05,min_value=-9999.
;         index=where(ulevel gt 0.)
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 c_labels=0*index,/follow,color=0,/overplot,min_value=-9999.
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[5,15],$
;                 c_labels=[0,0],/follow,color=0,/overplot,min_value=-9999.
;         index=where(ulevel lt 0.)
;         contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=ulevel(index),$
;                 /follow,color=mcolor,/overplot,c_linestyle=5,c_labels=0*index,min_value=-9999.
;;        contour,ubar_eq(0:icount-1,*),findgen(icount),press,/ylog,yrange=[1000.,1.],/noeras,levels=[0.],$
;;                /follow,color=mcolor,/overplot,thick=3,min_value=-9999.,c_labels=[0]
;         for ii=0L,nxticks-1L do xyouts,xindex(ii),3000.,xlabs(ii),/data,color=0,charsize=1.25,charthick=2,alignment=0.5
;         imin=min(ulevel)
;         imax=max(ulevel)
;         xmnb=xorig(2)+xlen+cbaryoff
;         xmxb=xmnb+cbarydel
;         set_viewport,xmnb,xmxb,yorig(2)+cbarydel,yorig(2)+ylen-cbarydel
;         !type=2^2+2^3+2^5+2^7
;         plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='m/s',$
;              color=0,/noeras,charsize=1.25,charthick=2
;         xbox=[0,10,10,0,0]
;         y1=imin
;         dy=(imax-imin)/float(nlvls)
;         for j=0,nlvls-1 do begin
;             ybox=[y1,y1,y1+dy,y1+dy,y1]
;             polyfill,xbox,ybox,color=col1(j)
;             y1=y1+dy
;         endfor
;;
;; save yearly file
;
;         nh_dtdy=reform(nh_dtdy(0:icount-1L,*))
;         sh_dtdy=reform(sh_dtdy(0:icount-1L,*))
;         nh_ubar65=reform(nh_ubar65(0:icount-1L,*))
;         sh_ubar65=reform(sh_ubar65(0:icount-1L,*))
;         ubar_eq=reform(ubar_eq(0:icount-1L,*))
;         yyyymmdd=yyyymmdd(0:icount-1L)
;         nday=icount
;         save,file='ERA40_dTdy_Ubar_QBO_'+yearlab+'.sav',nl,nday,yyyymmdd,press,nh_dtdy,sh_dtdy,$
;              nh_ubar65,sh_ubar65,ubar_eq,comment
;         icount=0L
;;
;; save jpg
;;
;         device,/close
;         spawn,'convert -trim ERA40_dTdy_Ubar_QBO_'+yearlab+'.ps -rotate -90 ERA40_dTdy_Ubar_QBO_'+yearlab+'.jpg'
;         goto,jump
;      endif
      skip:
      icount=icount+1L
goto,jump
saveit:

syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
sday=strmid(syyyymmdd_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nday=icount
save,file='ERA40_dTdy_Ubar_QBO_'+yearlab+'.sav',nl,nday,yyyymmdd_all,press,nh_dtdy_all,sh_dtdy_all,$
            nh_ubar65_all,sh_ubar65_all,ubar_eq_all,comment
end
