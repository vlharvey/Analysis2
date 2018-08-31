;---------------------------------------------------------------------------------------------------
; pressure data
; plot zt of polar cap temperature for MERRA, MERRA2, and the difference
;-----------------------------------------------------

nxdim=750
nydim=750
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=!p.color
icolmax=255
mcolor=icolmax
device,decompose=0
!NOERAS=-1
nlvls=21
col1=(findgen(nlvls)/float(nlvls))*mcolor
xorig=[0.15,0.15,0.15]
yorig=[0.7,0.4,0.1]
xlen=0.7
ylen=0.225
cbaryoff=0.08
cbarydel=0.01
!NOERAS=-1
!p.font=1
a=findgen(8)*(2*!pi/8.)
usersym,0.3*cos(a),0.3*sin(a),/fill
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,font_size=9
   device,/landscape,bits=8,filename='zt_merra_merra2_diff_zt.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; restore
;
yearlab='19800101-20160131'
restore,'/Volumes/Data/MERRA_data/Pre_process/merra_ZM-Z_'+yearlab+'.sav'   ;,sdates,pressure,lat,zbar2d
restore,'/Volumes/Data/MERRA_data/Pre_process/merra_ZM-T_'+yearlab+'.sav'   ;,sdates,pressure,lat,tbar2d
;index=where(sdates ge 20020101 and sdates le 20080101)
;index=where(sdates ge 20060101 and sdates le 20061201)
;index=where(sdates ge 20040101 and sdates le 20041201)
;sdates=sdates(index)
;TBAR2D=TBAR2D(index,*,*)
;ZBAR2D=ZBAR2D(index,*,*)

sdates_merra=sdates
yindex=where(lat gt 60.)
tpc_merra=mean(TBAR2D(*,yindex,*),dim=2)
zpc_merra=mean(ZBAR2D(*,yindex,*),dim=2)
syear=strmid(sdates_merra,0,4)
xindex=where(strmid(sdates_merra,4,4) eq '0101' and long(syear) mod 5 eq 0,nxticks)
yearlabel=' '
if nxticks gt 1 then xlab=strmid(sdates_merra(xindex),0,4)
if nxticks le 1 then begin
   xindex=where(strmid(sdates_merra,6,2) eq '01',nxticks)
   xlab=strmid(sdates_merra(xindex),4,2)
   yearlabel=strmid(sdates_merra(xindex(0)),0,4)
endif

!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
tlevel=170+5.*findgen(nlvls)
contour,tpc_merra,findgen(n_elements(SDATES_merra)),zpc_merra,levels=tlevel,c_color=col1,ytitle='Altitude (km)',/cell,xtickname=xlab,xticks=nxticks-1,xtickv=xindex,yrange=[30,80],title='MERRA',charsize=1.5,color=0
;contour,tpc_merra,findgen(n_elements(SDATES_merra)),zpc_merra,levels=[270],color=0,/overplot,/follow,thick=3
xyouts,5,70,yearlabel,charsize=3,color=0,charthick=2

restore,'merra2_ZM-Z_'+yearlab+'.sav'	;,sdates,pressure,lat,zbar2d
restore,'merra2_ZM-T_'+yearlab+'.sav'	;,sdates,pressure,lat,tbar2d
;index=where(sdates ge 20020101 and sdates le 20080101)
;index=where(sdates ge 20060101 and sdates le 20061201)
;index=where(sdates ge 20040101 and sdates le 20041201)
;sdates=sdates(index)
;TBAR2D=TBAR2D(index,*,*)
;ZBAR2D=ZBAR2D(index,*,*)

sdates_merra2=sdates
yindex=where(lat gt 60.)
tpc_merra2=mean(TBAR2D(*,yindex,*),dim=2)
zpc_merra2=mean(ZBAR2D(*,yindex,*),dim=2)

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,tpc_merra2,findgen(n_elements(SDATES_merra2)),zpc_merra2,levels=tlevel,c_color=col1,ytitle='Altitude (km)',/cell,xtickname=xlab,xticks=nxticks-1,xtickv=xindex,yrange=[30,80],title='MERRA-2',charsize=1.5,color=0
;contour,tpc_merra2,findgen(n_elements(SDATES_merra2)),zpc_merra2,levels=[270],color=0,/overplot,/follow,thick=3

level  = tlevel
nlvls  = n_elements(level)
slab=' '+strarr(n_elements(level))
!type=2^2+2^3+2^5+2^6
plot,[0,0],[0,0],xrange=[0,10],yrange=[0,1],/noeras,yticks=n_elements(level)-1L,$
      position = [.9,.4,.94,.95],ytickname=slab,/nodata
xyouts,.97,.5,'60-90N Temperature (K)',/normal,orientation=90,color=0,charsize=1.5,charthick=2
xbox=[0,10,10,0,0]
y2=0
dy= 1./(n_elements(level)-1.)
for j=0,n_elements(col1)-2 do begin
    ybox=[y2,y2,y2+dy,y2+dy,y2]
    polyfill,xbox,ybox,color=col1[j]
    y2=y2+dy
endfor
loadct,0
slab=strcompress(string(format='(i3)',level),/remove_all)
slabcolor = fltarr(n_elements(level))*0.
slabcolor[0:8] = 255        ; set first few labels to white so they are visible
y1=dy/2 ; center of first color level
for i=0L,n_elements(slab)-2L do begin
    slab0=slab[i]
    if i mod 2 eq 0 then xyouts,5,y1-dy/2.,slab0,charsize=1.3,/data,color=slabcolor[i],align = .5 ; This should place the label on the left side of each color level
    y1=y1+dy
endfor

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
dlevel=[-50,-20,-10,-5,-2,0,2,5,10,20,50]
;dlevel=-10+2.*findgen(11)
;dlevel=-25+5.*findgen(11)
;dlevel=-50+10.*findgen(11)
diff=tpc_merra2-tpc_merra
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,diff,findgen(n_elements(SDATES_merra2)),zpc_merra2,levels=dlevel,c_color=col2,ytitle='Altitude (km)',/cell,xtickname=xlab,xticks=nxticks-1,xtickv=xindex,yrange=[30,80],title='MERRA-2 minus MERRA',charsize=1.5,color=0
;contour,tpc_merra2,findgen(n_elements(SDATES_merra2)),zpc_merra2,levels=[270],color=0,/overplot,/follow,thick=3

level  = dlevel
nlvls  = n_elements(level)
slab=' '+strarr(n_elements(level))
!type=2^2+2^3+2^5+2^6
plot,[0,0],[0,0],xrange=[0,10],yrange=[0,1],/noeras,yticks=n_elements(level)-1L,$
      position = [.9,.1,.94,.35],ytickname=slab,/nodata
xyouts,.97,.25,'dT (K)',/normal,orientation=90,color=0,charsize=1.5,charthick=2
xbox=[0,10,10,0,0]
y2=0
dy= 1./(n_elements(level)-1.)
for j=0,n_elements(col2)-2 do begin
    ybox=[y2,y2,y2+dy,y2+dy,y2]
    polyfill,xbox,ybox,color=col2[j]
    y2=y2+dy
endfor
loadct,0
slab=strcompress(string(format='(i3)',level),/remove_all)
slabcolor = fltarr(n_elements(level))*0.
slabcolor[0:2] = mcolor        ; set first few labels to white so they are visible
y1=dy/2 ; center of first color level
for i=0L,n_elements(slab)-2L do begin
    slab0=slab[i]
    xyouts,5,y1-dy/2.,slab0,charsize=1.3,/data,color=slabcolor[i],align = .5 ; This should place the label on the left side of each color level
    y1=y1+dy
endfor


if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_merra_merra2_diff_zt.ps -rotate -90 zt_merra_merra2_diff_zt.png'
endif

end
