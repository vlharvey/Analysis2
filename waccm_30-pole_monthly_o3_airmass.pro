;
; average individual monthly averages to get multi-year monthly averages
; plot annual cycle from 30-pole to compare to occultation plot
;
@fillit
@smoothit

setplot='x'
read,'setplot?',setplot
loadct,38
device,decompose=0
mcolor=byte(!p.color)
mcolor=fix(mcolor)
if mcolor eq 0 then mcolor=255
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
icmm1=mcolor-1
icmm2=mcolor-2
!noeras=1
a=findgen(6)*(2*!pi/6.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
yorig=[0.5,0.5,0.5,0.125,0.125,0.125]
xorig=[0.1,0.4,0.7,0.1,0.4,0.7]
ylen=0.3
xlen=0.25
cbaryoff=0.06
cbarydel=0.02
nhxlabels=['J','A','S','O','N','D','J','F','M','A','M','J']
shxlabels=['J','F','M','A','M','J','J','A','S','O','N','D']
if setplot ne 'ps' then begin
   lc=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/portrait,bits=8,filename='waccm_30-pole_o3_airmass.ps'
   device,/color
   device,/inch,xoff=0.05,yoff=.1,xsize=xsize,ysize=ysize
endif
dir='/aura3/data/WACCM_data/Datfiles/waccm_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
nlat=18
ybin=-85.+10.*findgen(nlat)
nmonth=12L
for imonth=0L,nmonth-1L do begin
    spawn,'ls -1 '+dir+'*airmass*'+mon(imonth)+'*.sav',ifiles
    for i=0L,n_elements(ifiles)-1L do begin
        restore,ifiles(i)
print,'restored '+ifiles(i)
        if i eq 0L then begin
           HIGH_O3_AVG_ALL=0.*high_o3_avg
           HIGH_O3_NUM_ALL=0L*high_o3_num
           HIGH_O3_SIG_ALL=0.*high_o3_avg
           OUT_O3_AVG_ALL=0.*out_o3_avg
           OUT_O3_NUM_ALL=0L*out_o3_num
           OUT_O3_SIG_ALL=0.*out_o3_avg
           VORTEX_O3_AVG_ALL=0.*vortex_o3_avg
           VORTEX_O3_NUM_ALL=0L*vortex_o3_num
           VORTEX_O3_SIG_ALL=0.*vortex_o3_avg
nth=n_elements(th)
x2d=0.*vortex_o3_avg
y2d=0.*vortex_o3_avg
for ii=0,nlat-1 do y2d(ii,*)=th
for j=0,nth-1 do x2d(*,j)=ybin
;
; annual cycle arrays
;
           if imonth eq 0L then begin
           out_annual_30np=fltarr(nmonth,nth)
           out_annual_num_30np=fltarr(nmonth,nth)
           high_annual_30np=fltarr(nmonth,nth)
           high_annual_num_30np=fltarr(nmonth,nth)
           vortex_annual_30np=fltarr(nmonth,nth)
           vortex_annual_num_30np=fltarr(nmonth,nth)
           out_annual_30sp=fltarr(nmonth,nth)
           out_annual_num_30sp=fltarr(nmonth,nth)
           high_annual_30sp=fltarr(nmonth,nth)
           high_annual_num_30sp=fltarr(nmonth,nth)
           vortex_annual_30sp=fltarr(nmonth,nth)
           vortex_annual_num_30sp=fltarr(nmonth,nth)
           endif
        endif
        index=where(high_o3_avg gt 0.)
        if index(0) ne -1L then begin
        HIGH_O3_AVG_ALL(index)=HIGH_O3_AVG_ALL(index)+high_o3_avg(index)*float(high_o3_num(index))
        HIGH_O3_NUM_ALL(index)=HIGH_O3_NUM_ALL(index)+high_o3_num(index)
        HIGH_O3_SIG_ALL(index)=HIGH_O3_SIG_ALL(index)+high_o3_sig(index)*float(high_o3_num(index))
        endif
        index=where(out_O3_avg gt 0.)
        if index(0) ne -1L then begin
        OUT_O3_AVG_ALL(index)=OUT_O3_AVG_ALL(index)+out_o3_avg(index)*float(out_o3_num(index))
        OUT_O3_NUM_ALL(index)=OUT_O3_NUM_ALL(index)+out_o3_num(index)
        OUT_O3_SIG_ALL(index)=OUT_O3_SIG_ALL(index)+out_o3_sig*float(out_o3_num(index))
        endif
        index=where(vortex_O3_avg gt 0.)
        if index(0) ne -1L then begin
        VORTEX_O3_AVG_ALL(index)=VORTEX_O3_AVG_ALL(index)+vortex_o3_avg(index)*float(vortex_o3_num(index))
        VORTEX_O3_NUM_ALL(index)=VORTEX_O3_NUM_ALL(index)+vortex_o3_num(index)
        VORTEX_O3_SIG_ALL(index)=VORTEX_O3_SIG_ALL(index)+vortex_o3_sig(index)*float(vortex_o3_num(index))
        endif
    endfor
;
; multi-year monthly averages
;
    index=where(HIGH_O3_AVG_ALL gt 0.)
    HIGH_O3_AVG_ALL(index)=HIGH_O3_AVG_ALL(index)/HIGH_O3_NUM_ALL(index)
    HIGH_O3_SIG_ALL(index)=HIGH_O3_SIG_ALL(index)/HIGH_O3_NUM_ALL(index)
    index=where(OUT_O3_AVG_ALL gt 0.)
    OUT_O3_AVG_ALL(index)=OUT_O3_AVG_ALL(index)/OUT_O3_NUM_ALL(index)
    OUT_O3_SIG_ALL(index)=OUT_O3_SIG_ALL(index)/OUT_O3_NUM_ALL(index)
    index=where(VORTEX_O3_AVG_ALL gt 0.)
    VORTEX_O3_AVG_ALL(index)=VORTEX_O3_AVG_ALL(index)/VORTEX_O3_NUM_ALL(index)
    VORTEX_O3_SIG_ALL(index)=VORTEX_O3_SIG_ALL(index)/VORTEX_O3_NUM_ALL(index)
;
; loop over altitude
;
    for k=0L,nth-1L do begin
;       dum=reform(out_O3_AVG_ALL(*,k))
;       ndum=reform(out_O3_num_ALL(*,k))
;       index=where(ybin gt 30. and dum gt 0.,klat)
;       if index(0) ne -1 then out_annual_30np(imonth,k)=total(dum(index))/float(klat)
;       if index(0) ne -1 then out_annual_num_30np(imonth,k)=total(ndum(index))
;       index=where(ybin lt -30. and dum gt 0.,klat)
;       if index(0) ne -1 then out_annual_30sp(imonth,k)=total(dum(index))/float(klat)
;       if index(0) ne -1 then out_annual_num_30sp(imonth,k)=total(ndum(index))

;       dum=reform(high_O3_AVG_ALL(*,k))
;       ndum=reform(high_O3_num_ALL(*,k))
;       index=where(ybin gt 30. and dum gt 0.,klat)
;       if index(0) ne -1 then high_annual_30np(imonth,k)=total(dum(index))/float(klat)
;       if index(0) ne -1 then high_annual_num_30np(imonth,k)=total(ndum(index))
;       index=where(ybin lt -30. and dum gt 0.,klat)
;       if index(0) ne -1 then high_annual_30sp(imonth,k)=total(dum(index))/float(klat)
;       if index(0) ne -1 then high_annual_num_30sp(imonth,k)=total(ndum(index))


        odum=reform(out_O3_AVG_ALL(*,k))
        nodum=reform(out_O3_num_ALL(*,k))
        hdum=reform(high_O3_AVG_ALL(*,k))
        nhdum=reform(high_O3_num_ALL(*,k))
        index=where(ybin gt 30. and odum gt 0. and hdum gt 0.,klat)
        if index(0) ne -1 then out_annual_30np(imonth,k)=total(odum(index))/float(klat)
        if index(0) ne -1 then out_annual_num_30np(imonth,k)=total(nodum(index))
        if index(0) ne -1 then high_annual_30np(imonth,k)=total(hdum(index))/float(klat)
        if index(0) ne -1 then high_annual_num_30np(imonth,k)=total(nhdum(index))
        index=where(ybin lt -30. and odum gt 0. and hdum gt 0.,klat)
        if index(0) ne -1 then out_annual_30sp(imonth,k)=total(odum(index))/float(klat)
        if index(0) ne -1 then out_annual_num_30sp(imonth,k)=total(nodum(index))
        if index(0) ne -1 then high_annual_30sp(imonth,k)=total(hdum(index))/float(klat)
        if index(0) ne -1 then high_annual_num_30sp(imonth,k)=total(nhdum(index))

        dum=reform(vortex_O3_AVG_ALL(*,k))
        ndum=reform(vortex_O3_num_ALL(*,k))
        index=where(ybin gt 30. and dum gt 0.,klat)
        if index(0) ne -1 then vortex_annual_30np(imonth,k)=total(dum(index))/float(klat)
        if index(0) ne -1 then vortex_annual_num_30np(imonth,k)=total(ndum(index))
        index=where(ybin lt -30. and dum gt 0.,klat)
        if index(0) ne -1 then vortex_annual_30sp(imonth,k)=total(dum(index))/float(klat)
        if index(0) ne -1 then vortex_annual_num_30sp(imonth,k)=total(ndum(index))
    endfor
endfor		; loop over months

nhroto=0.*out_annual_30np
nhroth=0.*high_annual_30np
nhroth2=0.*high_annual_30np
tmpo=0.*out_annual_30np
tmph=0.*high_annual_30np
tmph2=0.*high_annual_30np
index=where(out_annual_30np gt 0.)
tmpo(index)=out_annual_30np(index)
index=where(out_annual_30np le 0.)
if index(0) ne -1 then tmpo(index)=-9999.
index=where(high_annual_30np gt 0.)
tmph(index)=high_annual_30np(index)
index=where(high_annual_30np le 0.)
if index(0) ne -1 then tmph(index)=-9999.
index=where(high_annual_30np gt 0. and out_annual_30np gt 0.)
tmph2(index)=high_annual_30np(index)-out_annual_30np(index)
index=where(high_annual_30np le 0.)
if index(0) ne -1 then tmph2(index)=-9999.
for k=0,nth-1 do begin
    nhroto(0:5,k)=out_annual_30np(6:11,k)
    nhroto(6:11,k)=out_annual_30np(0:5,k)
    nhroth(0:5,k)=tmph(6:11,k)
    nhroth(6:11,k)=tmph(0:5,k)
    nhroth2(0:5,k)=tmph2(6:11,k)
    nhroth2(6:11,k)=tmph2(0:5,k)
endfor
erase
!type=2^2+2^3
xyouts,.3,.9,'WACCM',/normal,charsize=3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=findgen(13)
nlev=n_elements(level)
col1=1+indgen(nlev)*mcolor/nlev
yindex=where((th le 1600. and th ge 1000.) or th eq 800. or $
              th eq 600. or th eq 400.,nyticks)
sth=strcompress(string(fix(th(yindex))),/remove_all)
contour,nhroto,findgen(nmonth),th,/noeras,title='30-NP Ambient',$
     ytitle='Potential Temperature (K)',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=nhxlabels,$
     max_value=9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,yticks=nyticks-1,ytickname=sth,ytickv=th(yindex),charsize=1.25
contour,nhroto,findgen(nmonth),th,/noeras,/follow,levels=level,$
        c_labels=1+0*level,max_value=9999.,/overplot,color=0
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,nhroth,findgen(nmonth),th,/noeras,$
     ytitle=' ',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=nhxlabels,$
     min_value=-9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,ytickname=' '+strarr(nyticks),title='30N-NP Anticyclones',charsize=1.25
contour,nhroth,findgen(nmonth),th,/noeras,/follow,levels=level,$
        c_labels=1+0*level,max_value=9999.,/overplot,color=0
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6

!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=-2.0+.2*findgen(21)
nlev=n_elements(level)
col1=1+indgen(nlev)*mcolor/nlev
contour,nhroth2,findgen(nmonth),th,/noeras,$
     ytitle=' ',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=nhxlabels,$
     min_value=-9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,ytickname=' '+strarr(nyticks),title='Anticyclones-Ambient',charsize=1.25
index=where(level lt 0.)
contour,nhroth2,findgen(nmonth),th,/noeras,/follow,levels=level(index),$
        c_labels=0*level(index),min_value=-9999.,/overplot,color=mcolor,charsize=2
index=where(level gt 0.)
contour,nhroth2,findgen(nmonth),th,/noeras,/follow,levels=level(index),$
        c_labels=0*level(index),min_value=-9999.,/overplot,color=0,charsize=2
contour,nhroth2,findgen(nmonth),th,/noeras,/follow,levels=[0],c_labels=[0],$
        min_value=-9999.,/overplot,color=0,thick=3
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6
xyouts,12.5,1000.,'Altitude (km)',/data,alignment=.5,orientation=-90.,charsize=1.25
xyouts,11.3,1580.,'45',/data,charsize=1.25
xyouts,11.3,1400.,'40',/data,charsize=1.25
xyouts,11.3,1200.,'36',/data,charsize=1.25
xyouts,11.3,1000.,'32',/data,charsize=1.25
xyouts,11.3,800.,'30',/data,charsize=1.25
xyouts,11.3,600.,'23',/data,charsize=1.25
xyouts,11.3,400.,'15',/data,charsize=1.25

tmpo=0.*out_annual_30sp
tmph=0.*high_annual_30sp
tmph2=0.*high_annual_30sp
index=where(out_annual_30sp gt 0.)
tmpo(index)=out_annual_30sp(index)
index=where(out_annual_30sp le 0.)
if index(0) ne -1 then tmpo(index)=-9999.
index=where(high_annual_30sp gt 0.)
tmph(index)=high_annual_30sp(index)
index=where(high_annual_30sp le 0.)
if index(0) ne -1 then tmph(index)=-9999.
index=where(high_annual_30sp gt 0.)
tmph2(index)=high_annual_30sp(index)-out_annual_30sp(index)
index=where(high_annual_30sp le 0.)
if index(0) ne -1 then tmph2(index)=-9999.

!type=2^2+2^3
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=findgen(13)
nlev=n_elements(level)
col1=1+indgen(nlev)*mcolor/nlev
contour,tmpo,findgen(nmonth),th,/noeras,title='30-SP Ambient',$
     ytitle='Potential Temperature (K)',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=shxlabels,$
     max_value=9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,yticks=nyticks-1,ytickname=sth,ytickv=th(yindex),charsize=1.25
contour,tmpo,findgen(nmonth),th,/noeras,/follow,levels=level,$
        c_labels=1+0*level,max_value=9999.,/overplot,color=0
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6
!psym=0
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xorig(3),xorig(3)+xlen,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xtitle='Ozone (ppmv)',$
    xrange=[imin,imax],/noeras,charsize=1.25
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/float(nlev)
for j=0,nlev-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor

!type=2^2+2^3
xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,tmph,findgen(nmonth),th,/noeras,$
     ytitle=' ',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=shxlabels,$
     min_value=-9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,ytickname=' '+strarr(nyticks),title='30-SP Anticyclones',charsize=1.25
contour,tmph,findgen(nmonth),th,/noeras,/follow,levels=level,$
        c_labels=1+0*level,max_value=9999.,/overplot,color=0
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6
!psym=0
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xorig(4),xorig(4)+xlen,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xtitle='Ozone (ppmv)',$
    xrange=[imin,imax],/noeras,charsize=1.25
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/float(nlev)
for j=0,nlev-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor

!type=2^2+2^3
xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
level=-2.0+.2*findgen(21)
;level=-50.+5.*findgen(21)
nlev=n_elements(level)
col1=1+indgen(nlev)*mcolor/nlev
col1(nlev-1)=col1(nlev-2)
contour,tmph2,findgen(nmonth),th,/noeras,$
     ytitle=' ',xrange=[0.,nmonth-1],yrange=[400.,1600.],xtickname=shxlabels,$
     min_value=-9999.,/cell_fill,c_color=col1,xticks=nmonth-1,levels=level,$
     xticklen=-0.02,ytickname=' '+strarr(nyticks),title='Anticyclones-Ambient',charsize=1.25
index=where(level lt 0.)
contour,tmph2,findgen(nmonth),th,/noeras,/follow,levels=level(index),$
        c_labels=1+0*level(index),min_value=-9999.,/overplot,color=mcolor,charsize=2
index=where(level gt 0.)
contour,tmph2,findgen(nmonth),th,/noeras,/follow,levels=level(index),$
        c_labels=1+0*level(index),min_value=-9999.,/overplot,color=0,charsize=2
contour,tmph2,findgen(nmonth),th,/noeras,/follow,levels=[0],c_labels=[0],$
        min_value=-9999.,/overplot,color=0,thick=3
plots,6,400.
plots,6,320.,/continue,color=lc,thick=6
xyouts,12.5,1000.,'Altitude (km)',/data,alignment=.5,orientation=-90.,charsize=1.25
xyouts,11.3,1580.,'45',/data,charsize=1.25
xyouts,11.3,1400.,'40',/data,charsize=1.25
xyouts,11.3,1200.,'36',/data,charsize=1.25
xyouts,11.3,1000.,'32',/data,charsize=1.25
xyouts,11.3,800.,'30',/data,charsize=1.25
xyouts,11.3,600.,'23',/data,charsize=1.25
xyouts,11.3,400.,'15',/data,charsize=1.25
!psym=0
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xorig(5),xorig(5)+xlen,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xtitle='Ozone Difference (ppmv)',$
    xrange=[imin,imax],/noeras,charsize=1.25
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/float(nlev)
for j=0,nlev-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor
if setplot eq 'ps' then device,/close
if setplot ne 'ps' then stop
end
