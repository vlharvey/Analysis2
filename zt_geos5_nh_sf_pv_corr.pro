;
; plot GEOS-5 altitude-time section of the correlation coefficient
; between SF and PV within each of the vortex definitions
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,38
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,0.15]
yorig=[0.55,0.15]
xlen=0.75
ylen=0.35
cbaryoff=0.04
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
stimes=[$
'_AVG.V01.']
slabs=['AVG']
ntimes=n_elements(stimes)
!noeras=1
dirm='/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
lstmn=1L & lstdy=17L & lstyr=2007L
ledmn=4L & leddy=1L & ledyr=2007L
lstday=0L & ledday=0L
;
; get date range
;
print, ' '
print, '      GEOS-5 Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdates=strarr(kday)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
      sdates(kcount)=sdate
      kcount=kcount+1L
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+sdate+stimes(0)+'nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump

;
; read new vortex
;
      ncid=ncdf_open(dir+sdate+stimes(0)+'nc5')
      marknew2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,3,marknew2
      ncdf_close,ncid
;
; zt arrays of correlation coefficients
;
      if icount eq 0L then begin
         rcorr3=-999.+0.*fltarr(kday,nth)
         rcorr5=-999.+0.*fltarr(kday,nth)
         y2d=fltarr(nc,nr)
         for i=0,nc-1 do y2d(i,*)=alat
         icount=1L
      endif
;
; loop over theta surfaces
;
      for ilev=0L,nth-1L do begin
          mark=transpose(mark2(*,*,ilev))
          marknew=transpose(marknew2(*,*,ilev))
          sf=transpose(sf2(*,*,ilev))
          pv=transpose(pv2(*,*,ilev))
;
; PV vs SF correlations in vortex
; 
          index=where(y2d gt 40.)	; and mark gt 0.)
sf1=sf(index) & pv1=pv(index)
index=sort(sf1)
          if index(0) ne -1L then rcorr3(kcount-1,ilev)=correlate(sf1(index),pv1(index))
          index=where(y2d gt 0. and marknew gt 0.)
;          if index(0) ne -1L then rcorr5(kcount-1,ilev)=correlate(sf(index),pv(index))
print,th(ilev),rcorr3(kcount-1,ilev)
      endfor
erase
plot,rcorr3(kcount-1,*),th,color=0,title=sdate,xrange=[-1.,1.]
stop
goto,jump
;
; plot
;
plotit:
;
; save postscript version
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_geos5_nh_sf_pv_corr_'+syr+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2
endif
;
; Isotachs + vortex edge
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
xindex=where(strmid(sdates,6,2) eq '15')
xindex=[0,xindex,kday-1]
xlabs=sdates(xindex)
xindex=1.+xindex
nxticks=n_elements(xindex)
nlvls=21
col1=reverse(1+indgen(nlvls)*icolmax/nlvls)
level=-1.+0.1*findgen(nlvls)
index=where(rcorr3 eq -999.)
if index(0) ne -1L then rcorr3(index)=0./0.
rcorr3=smooth(rcorr3,3,/edge_truncate,/nan)
contour,rcorr3,1.+findgen(kday),th,levels=level,/fill,/cell_fill,c_color=col1,$
        title='Arctic Vortex PV SF Correlation (Strongest Jet)',min_value=-999.,xrange=[1.,kday],$
        yrange=[min(th),max(th)],ytitle='Theta (K)',/noeras,color=0,xticks=nxticks-1,xtickv=xindex,$
        xtickname=' '+strarr(nxticks)
index=where(level lt 0.)
contour,rcorr3,1.+findgen(kday),th,levels=level(index),c_color=mcolor,/follow,/overplot,min_value=-999.
index=where(level gt 0.)
contour,rcorr3,1.+findgen(kday),th,levels=level(index),c_color=0,/follow,/overplot,min_value=-999.
contour,rcorr3,1.+findgen(kday),th,levels=[0.],c_color=0,/follow,/overplot,thick=3,min_value=-999.
set_viewport,xmx+cbaryoff,xmx+cbaryoff+cbarydel,ymn,ymx
!type=2^2+2^3+2^5
omin=min(level)
omax=max(level)
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
index=where(rcorr5 eq -999.)
if index(0) ne -1L then rcorr5(index)=0./0.
rcorr5=smooth(rcorr5,3,/edge_truncate,/nan)
contour,rcorr5,1.+findgen(kday),th,levels=level,/fill,/cell_fill,c_color=col1,$
        title='Arctic Vortex PV SF Correlation (Poleward-Most Jet)',min_value=-999.,xrange=[1.,kday],$
        yrange=[min(th),max(th)],ytitle='Theta (K)',/noeras,color=0,xticks=nxticks-1,xtickv=xindex,$
        xtickname=' '+strarr(nxticks)
index=where(level lt 0.)
contour,rcorr5,1.+findgen(kday),th,levels=level(index),c_color=mcolor,/follow,/overplot,min_value=-999.
index=where(level gt 0.)
contour,rcorr5,1.+findgen(kday),th,levels=level(index),c_color=0,/follow,/overplot,min_value=-999.
contour,rcorr5,1.+findgen(kday),th,levels=[0.],c_color=0,/follow,/overplot,thick=3,min_value=-999.
for i=0,nxticks-1 do xyouts,xindex(i),-200.,xlabs(i),alignment=0.5,orientation=90.,color=0,/data
set_viewport,xmx+cbaryoff,xmx+cbaryoff+cbarydel,ymn,ymx
!type=2^2+2^3+2^5
omin=min(level)
omax=max(level)
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim zt_geos5_nh_sf_pv_corr_'+syr+'.ps -rotate -90 '+$
         'zt_geos5_nh_sf_pv_corr_'+syr+'.jpg'
endif
end
