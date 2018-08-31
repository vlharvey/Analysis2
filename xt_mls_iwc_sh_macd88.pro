;
; Hovmoller of MLS IWC derived from T and H2O using Hervig's equation
; SH version
;
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.25]
yorig=[0.20]
xlen=0.5
ylen=0.7
cbaryoff=0.12
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
for lstyr=2007,2012 do begin
slstyr=strcompress(lstyr,/remove_all)
ledyr=lstyr+1L
sledyr=strcompress(ledyr,/remove_all)
lday=28
if ledyr mod 4 eq 0 then lday=29
kday=92L+lday
restore,'avg_mls_iwc_nov_'+slstyr+'.sav'
icenov=ice
restore,'avg_mls_iwc_dec_'+slstyr+'.sav'
icedec=ice
restore,'avg_mls_iwc_jan_'+sledyr+'.sav'
icejan=ice
restore,'avg_mls_iwc_feb_'+sledyr+'.sav'
icefeb=ice
;
; create NDJF ICE
;
nr=n_elements(latbin)
nc=n_elements(lonbin)
icendjf=fltarr(nc,nr,kday)
icendjf(*,*,0:29)=icenov
icendjf(*,*,30:30+31-1)=icedec
icendjf(*,*,30+31:30+31+31-1)=icejan
icendjf(*,*,30+31+31:30+31+31+lday-1)=icefeb
;
; year date label
;
yearlab=slstyr+'-'+sledyr
;
rlat=-55.
print,latbin
;read,'Enter latitude ',rlat
slat=strcompress(long(rlat),/remove_all)
index=where(latbin eq rlat)
ilat=index(0)
xtice=reform(icendjf(*,ilat,*))
;
; interpolate small gaps in time
;
for k=0,kday-1 do begin
    dlev=reform(xtice(*,k))
    for i=1,nc-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,nc-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    xtice(*,k)=dlev
endfor

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='xt_mls_iwc_'+yearlab+'_'+slat+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; add wrap-around point in longitude
;
lonbin2=30.*findgen(nc+1)
xtice2=fltarr(nc+1,kday)
for k=0L,kday-1L do begin
    for i=0L,nc-2L do begin
        if xtice(i,k) gt 0. and xtice(i+1,k) gt 0. then xtice2(i+1,k)=(xtice(i,k)+xtice(i+1,k))/2.
    endfor
    if xtice(0,k) gt 0. and xtice(nc-1,k) gt 0. then xtice2(0,k)=(xtice(0,k)+xtice(nc-1,k))/2.
    xtice2(nc,k)=xtice(0,k)
endfor
;
; plot xt IWC
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
index=where(xtice2 le 0.)
if index(0) ne -1L then xtice2(index)=0./0.
xtice2=smooth(xtice2,3,/NaN,/edge_truncate)
if index(0) ne -1L then xtice2(index)=0./0.

tlevel=12.5*findgen(20)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,xtice2,lonbin2,1.+findgen(kday),/noeras,yrange=[kday,1.],xrange=[0.,360.],xticks=6,$
      charsize=1.5,color=0,xtitle='Longitude',/cell_fill,c_color=col1,title=yearlab+' Lat= '+slat+' 83km',$
      levels=tlevel,yticks=3,ytickname='15 '+['Nov','Dec','Jan','Feb'],ytickv=[15.,15.+30.,15.+30.+31.,15.+30.+31.+31.],min_value=-99.
contour,xtice2,lonbin2,1.+findgen(kday),levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Derived IWC (ug/m2)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim xt_mls_iwc_'+yearlab+'_'+slat+'.ps -rotate -90 xt_mls_iwc_'+yearlab+'_'+slat+'.png'
    endif
endfor
end
