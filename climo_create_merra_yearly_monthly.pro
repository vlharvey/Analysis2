;
; create MERRA yearly files. monthly means
;
dirh='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
spawn,'ls '+dirh+'*daily_zm.sav',ifiles
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
model_years=1979+indgen(36)
model_years=string(FORMAT='(i4)',long(model_years))
nyears=n_elements(model_years)
;
; loop over model years
;
for iyear=0L,nyears-1L do begin

      dum=findfile(dirh+model_years(iyear)+'_daily_zm.sav')
      restore,dum(0)
;
; monthly averages
;
      nr=n_elements(latbin)
      nl=n_elements(pressure)
      kday=12L
      UGRD_MONAVG=fltarr(nr,nl,kday)
      VGRD_MONAVG=fltarr(nr,nl,kday)
      GPGRD_MONAVG=fltarr(nr,nl,kday)
      H2OGRD_MONAVG=fltarr(nr,nl,kday)
      TGRD_MONAVG=fltarr(nr,nl,kday)
      for n=0L,kday-1L do begin
          monindex=where(strmid(SDATE_ALL,0,2) eq string(FORMAT='(I2.2)',n+1))
          if monindex(0) ne -1L then begin
             UGRD_MON=reform(UGRD_AVG(*,*,monindex))
             VGRD_MON=reform(VGRD_AVG(*,*,monindex))
             GPGRD_MON=reform(GPGRD_AVG(*,*,monindex))
             H2OGRD_MON=reform(H2OGRD_AVG(*,*,monindex))
             TGRD_MON=reform(TGRD_AVG(*,*,monindex))
             for j=0L,nr-1L do begin
             for k=0L,nl-1L do begin
                 index=where(UGRD_MON(j,k,*) ne 0.)
                 if index(0) ne -1L then UGRD_MONAVG(j,k,n)=mean(UGRD_MON(j,k,index))
                 index=where(VGRD_MON(j,k,*) ne 0.)
                 if index(0) ne -1L then VGRD_MONAVG(j,k,n)=mean(VGRD_MON(j,k,index))
                 index=where(GPGRD_MON(j,k,*) ne 0.)
                 if index(0) ne -1L then GPGRD_MONAVG(j,k,n)=mean(GPGRD_MON(j,k,index))
                 index=where(H2OGRD_MON(j,k,*) ne 0.)
                 if index(0) ne -1L then H2OGRD_MONAVG(j,k,n)=mean(H2OGRD_MON(j,k,index))
                 index=where(TGRD_MON(j,k,*) ne 0.)
                 if index(0) ne -1L then TGRD_MONAVG(j,k,n)=mean(TGRD_MON(j,k,index))
             endfor
             endfor
         endif
      endfor
      print,max(TGRD_MONAVG),max(H2OGRD_MONAVG),max(GPGRD_MONAVG),max(VGRD_MONAVG),max(UGRD_MONAVG)
;
; check
;
erase
!type=2^2+2^3
loadct,39
mcolor=byte(!p.color)
mcolor=fix(mcolor)
device,decompose=0
if mcolor eq 0 then mcolor=255
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
ilat=where(min(abs(latbin-0.94736842)) eq abs(latbin-0.94736842))
plotarray=transpose(reform(TGRD_MONAVG(ilat,*,*)))
omin=min(plotarray)
omax=max(plotarray)
level=omin+((omax-omin)/nlvls)*findgen(nlvls+1)
contour,plotarray,1+findgen(kday),pressure,levels=level,c_color=col1,/cell_fill,/noeras,yrange=[max(pressure),min(pressure)],/ylog,$
        xrange=[1.,kday],xticks=11,ytitle='Pressure',xtitle='Month',charsize=1.5,min_value=-99.,title='MERRA Temp '+model_years(iyear)
contour,plotarray,1+findgen(kday),pressure,levels=level,c_color=0,/follow,/noeras,/overplot
;
; save yearly file
;
ofile=dirh+model_years(iyear)+'_monthly_zm.sav'
save,file=ofile,latbin,pressure,UGRD_MONAVG,VGRD_MONAVG,GPGRD_MONAVG,H2OGRD_MONAVG,TGRD_MONAVG
endfor	; loop over years
end
