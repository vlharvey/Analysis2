;
; ACE and SLIMCAT for Cynthia's JGR 2006 paper
; plot time-theta section of vortex methane
;
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,0.3*cos(a),0.3*sin(a),/fill
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=600 & nydim=600
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
month=['                      Dec',$
       '                      Jan',$
       '                      Feb',$
       '                      Mar',$
       '                      Apr',$
       '                         ']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/Theta/'
;
; restore SLIMCAT methane at ACE locations during 2004-2005
; SC_ACE          FLOAT     = Array[923, 42]
; SC_DATE_ACE     LONG      = Array[923]
; SC_ID_ACE       LONG      = Array[923]
; SC_LAT_ACE      FLOAT     = Array[923]
; SC_LON_ACE      FLOAT     = Array[923]
; ZO              INT       = Array[42]
restore,'CH4SC_interptoACE_2004_2005.sav
min_sc_id=min(SC_ID_ACE)
max_sc_id=max(SC_ID_ACE)
;
; flip in the vertical
;
SC_CH4=0.*SC_ACE
for i=0L,n_elements(SC_DATE_ACE)-1L do begin
    prof=reverse(reform(SC_ACE(i,*)))
    sc_ch4(i,*)=prof*1.e6
endfor
sc_ch4_save=sc_ch4
for i=0,39 do begin
    smoothit,sc_ch4,sc_ch4smooth
    sc_ch4=sc_ch4smooth
endfor
index=where(sc_ch4_save lt 0.)
;if index(0) ne -1L then sc_ch4(index)=-99.
;
; restore ACE data for 2004 and 2005
;
restore,'/aura3/data/ACE_data/Datfiles_SOSST/v2.2/cat_ace_v2.2.2004
DATE4=date
FDOY4=fdoy
ID4=id
LAT4=latitude
LON4=longitude
restore,'/aura3/data/ACE_data/Datfiles_SOSST/v2.2/Theta/ch4_ace_v2.2_theta.2004
mix4=MIX

restore,'/aura3/data/ACE_data/Datfiles_SOSST/v2.2/cat_ace_v2.2.2005
DATE5=date
FDOY5=fdoy
ID5=id
LAT5=latitude
LON5=longitude
restore,'/aura3/data/ACE_data/Datfiles_SOSST/v2.2/Theta/ch4_ace_v2.2_theta.2005
mix5=MIX
;
; concatenate ACE years
;
date=[date4,date5]
fdoy=[fdoy4,fdoy5]
id=[id4,id5]
lat=[lat4,lat5]
lon=[lon4,lon5]
mix=[mix4,mix5]
;
; remove ACE data outside SLIMCAT date range
;
index=where(id ge min_sc_id and id le max_sc_id,nday)
date=date(index)
fdoy=fdoy(index)
id=id(index)
lat=lat(index)
lon=lon(index)
mix=mix(index,*)
;
; remove ACE data that does not have corresponding SLIMCAT data
;
flag=fltarr(nday)
for i=0L,nday-1L do begin
    x=where(id(i) eq sc_id_ace)
    if x(0) eq -1L then flag(i)=1.
;   if x(0) ne -1L then print,i,id(i),sc_id_ace(x)
endfor
index=where(flag eq 0.,nday)
date=date(index)
fdoy=fdoy(index)
id=id(index)
lat=lat(index)
lon=lon(index)
mix=mix(index,*)*1.e6
;
; hardwired fix to fdoy for plotting
;
index=where(fdoy lt 330.)
fdoy(index)=fdoy(index)+max(fdoy)

mix_save=mix
for i=0,39 do begin
    smoothit,mix,mixsmooth
    mix=mixsmooth
endfor
index=where(mix_save lt 0.)
;if index(0) ne -1L then mix(index)=-99.


if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,font_size=9
   device,/landscape,bits=8,filename='zt_ace_ch4+slimcat.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
nlvls=n_elements(level)
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
contour,sc_ch4,fdoy,theta,/noeras,c_color=col1,/cell_fill,levels=level,$
        yrange=[400.,700.],xrange=[min(fdoy),max(fdoy)],xticks=n_elements(month)-1,xtickname=month,$
        charsize=2,ytitle='Alitude (km)',color=0,xticklen=-0.02,$
        title='Methane',min_value=-99.
contour,sc_ch4,fdoy,theta,/noeras,color=0,/follow,levels=level,/overplot,$
        c_labels=1*level,min_value=-99.,thick=2,c_charsize=2
contour,mix,fdoy,theta,/noeras,color=mcolor,/follow,levels=level,/overplot,$
        c_labels=1*level,min_value=-99.,thick=3,c_charsize=2

sdate=strcompress(date,/remove_all)
smmdd=strmid(sdate,4,4)
xx=where(strmid(sdate,7,1) eq '0' or strmid(sdate,7,1) eq '5',nn)
for ii=0L,nn-1L do begin
    if ii eq 0L then begin
       plots,fdoy(xx(ii)),320.
       plots,fdoy(xx(ii)),700.,color=0,/continue
       xyouts,fdoy(xx(ii)),320.,smmdd(xx(ii)),orientation=90,color=0.,charsize=1.75,alignment=0
    endif
    if ii gt 0L then begin
       if long(fdoy(xx(ii))) ne long(fdoy(xx(ii-1))) then begin
          plots,fdoy(xx(ii)),320.
          plots,fdoy(xx(ii)),700.,color=0,/continue
          xyouts,fdoy(xx(ii)),320.,smmdd(xx(ii)),orientation=90,color=0.,charsize=1.75,alignment=0
       endif
    endif
endfor
!type=2^2+2^3+2^5
xmnb=xorig(0)+xlen+cbaryoff
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
slab=strcompress(string(FORMAT='(F5.2)',level))
plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
     yticks=n_elements(level)-1L,ytickname=slab,$
     yrange=[min(level),max(level)],charsize=1.5,title='ppmv'
xbox=[0,10,10,0,0]
y1=min(level)
dy=(max(level)-min(level))/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_ace_ch4+slimcat.ps -rotate -90 zt_ace_ch4+slimcat.jpg'
endif
end
