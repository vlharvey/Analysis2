;
; contour maximum temperature timeseries in each hemisphere
;
loadct,38
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,.5*cos(a),.5*sin(a),/fill
!NOERAS=-1
!P.FONT=1
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.20,0.20]
yorig=[0.55,0.15]
xlen=0.6
ylen=0.3
cbaryoff=0.08
cbarydel=0.01
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='tz_maxt.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; restore max temperature from MetO
;
restore,'MetO_MaxT_Climo.sav'	;,yyyymmdd,th,nh_maxt_flag,sh_maxt_flag,comment
nday=n_elements(yyyymmdd)
syyyymmdd=strcompress(yyyymmdd,/remove_all)
syr=strmid(syyyymmdd,2,2)
smn=strmid(syyyymmdd,4,2)
sdy=strmid(syyyymmdd,6,2)
xindex=where(smn eq '01' and sdy eq '01',nxticks)
xlabs=syr(xindex)
n0=findgen(nxticks)
n1=1.+findgen(nxticks)
diff=abs(xindex(n0)-xindex(n1))
index=where(diff eq 0. or (diff ge 360. and diff le 367.),nxticks)
xindex=xindex(index)
;
; plot
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
xyouts,.2,.9,'MetO Maximum Daily Polar Temperature',/normal,charsize=2
!type=2^2+2^3
nlvls=21
level=200.+5.*findgen(nlvls)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,nh_maxt_flag,findgen(nday),th,levels=level,/fill,/cell_fill,c_color=col1,ytitle='Theta',$
        xtickname=syr(xindex),xtickv=xindex,xticks=nxticks-1,title='Northern Hemisphere',charsize=1.5,$
        yrange=[400.,2000.]
contour,nh_maxt_flag,findgen(nday),th,levels=[280.],/overplot,/follow,c_color=[0],c_labels=[0]
contour,nh_maxt_flag,findgen(nday),th,levels=[300.],/overplot,/follow,c_color=mcolor*.9,c_labels=[0],thick=2
contour,nh_maxt_flag,findgen(nday),th,levels=[320.],/overplot,/follow,c_color=mcolor*.1,c_labels=[0],thick=2

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,sh_maxt_flag,findgen(nday),th,levels=level,/fill,/cell_fill,c_color=col1,ytitle='Theta',$
        xtickname=syr(xindex),xtickv=xindex,xticks=nxticks-1,title='Southern Hemisphere',charsize=1.5,$
        yrange=[400.,2000.]
contour,sh_maxt_flag,findgen(nday),th,levels=[280.],/overplot,/follow,c_color=[0],c_labels=[0]
contour,sh_maxt_flag,findgen(nday),th,levels=[300.],/overplot,/follow,c_color=mcolor*.9,c_labels=[0],thick=2
contour,sh_maxt_flag,findgen(nday),th,levels=[320.],/overplot,/follow,c_color=mcolor*.1,c_labels=[0],thick=2
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,charsize=1.5
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
   spawn,'convert -trim tz_maxt.ps -rotate -90 tz_maxt.jpg'
endif
end
