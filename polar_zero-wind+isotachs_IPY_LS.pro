;
; Lower Stratosphere (LS)
; reads in .nc3 and plots polar projection zero wind line at a given theta level
; plus isotachs at 350 K in the upper troposphere
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra2_nc3

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'

mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=9L & lstdy=20L & lstyr=2009L
ledmn=10L & leddy=20L & ledyr=2009L
lstday=0L & ledday=0L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

stime=['00','06','12','18']
ntime=n_elements(stime)

kday=ledday-lstday+1L
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition '

      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      for itime=0L,ntime-1L do begin

      rd_merra2_nc3,dir+sdate+stime(itime)+'.nc3',nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag

      if iflag eq 1 then goto,jump
      if icount eq 0 then begin
         theta=1800.
;        print,th
;        read,'Enter theta ',theta
         index=where(theta eq th)
         if index(0) eq -1 then stop,'Invalid theta level '
         thlev=index(0)
         icount=1
      endif
      theta=th(thlev)
      stheta=strcompress(string(fix(theta)),/remove_all)
      qdf1=transpose(qdf2(*,*,thlev))
      sf1=transpose(sf2(*,*,thlev))
      pv1=transpose(pv2(*,*,thlev))
      z1=transpose(z2(*,*,thlev))
print,'mean height 45-90S of ',th(thlev),' is ',mean(z1(0:nr/4,*)),' km'
ustrat=transpose(u2(*,*,thlev))
vstrat=transpose(v2(*,*,thlev))
spstrat1=ustrat		;sqrt(ustrat^2.+vstrat^2.)

      u1=transpose(u2(*,*,-1))
      v1=transpose(v2(*,*,-1))	; wind speed at 350 K
      sp1=sqrt(u1^2.+v1^2.)
      mark1=transpose(mark2(*,*,thlev))
      qdf=0.*fltarr(nc+1,nr)
      qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
      qdf(nc,*)=qdf(0,*)
      sf=0.*fltarr(nc+1,nr)
      sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
      sf(nc,*)=sf(0,*)
      spstrat=0.*fltarr(nc+1,nr)
      spstrat(0:nc-1,0:nr-1)=spstrat1(0:nc-1,0:nr-1)
      spstrat(nc,*)=spstrat(0,*)

      sp=0.*fltarr(nc+1,nr)
      sp(0:nc-1,0:nr-1)=sp1(0:nc-1,0:nr-1)
      sp(nc,*)=sp(0,*)
      pv=0.*fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
      pv(nc,*)=pv(0,*)
      mark=0.*fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
      x=fltarr(nc+1)
      x(0:nc-1)=alon
      x(nc)=alon(0)+360.
      lon=0.*sf
      lat=0.*sf
      for i=0,nc   do lat(i,*)=alat
      for j=0,nr-1 do lon(*,j)=x

      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,font_size=9
         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                /bold,/color,bits_per_pixel=8,/times,filename='Figures/polar_sh_zerowind_'+sdate+stime(itime)+'_LS.ps'
         !p.charsize=1.25
         !p.thick=2
         !p.charthick=5
         !y.thick=2
         !x.thick=2
      endif
      erase
      !psym=0
      !type=2^2+2^3
      !p.charthick=2
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      MAP_SET,-90,0,90,/ortho,/noeras,/grid,/contin,/noborder,title=sdate+stime(itime),charsize=2.0,latdel=10,color=0,limit=[-90.,0,-40.,360.]
      map_continents,color=mcolor*.3,/fill_continents,/coasts
;
; indicate wedge of Justin's CIPS/AIRS data -90 to -45 longitude and -45 to -75 S
;
loadct,39
      oplot,-90+0*findgen(41),-80.+findgen(41),thick=10,color=250 
      oplot,-45+0*findgen(41),-80.+findgen(41),thick=10,color=250 
      oplot,-90+findgen(45),-40.+0*findgen(45),thick=10,color=250 
      oplot,-90+findgen(45),-80.+0*findgen(45),thick=10,color=250 

      index=where(lat lt 0.)
;     imin=min(pv(index))
imin=-0.0008
      imax=0
      nlvls=20
      pvint=(imax-imin)/float(nlvls)
      level=imin+pvint*findgen(nlvls)
;print,level
loadct,0
      col1=1+indgen(nlvls)*icolmax/float(nlvls)
;     contour,pv,x,alat,/overplot,levels=level,c_color=col1,/fill,/noeras
      contour,spstrat,x,alat,/overplot,levels=[0],/follow,c_labels=0*level,/noeras,color=0,thick=15
;
; vortex edge at all levels
;
;     contour,mark,x,alat,/overplot,levels=[-0.1],thick=10,color=20
;     loadct,0
;     contour,mark,x,alat,/overplot,levels=0.1*findgen(10),thick=10,color=0
;     for k=nth-1,0 do begin
;           theta=th(k)
;           stheta=strcompress(string(fix(theta)),/remove_all)
;           mark1=transpose(mark2(*,*,k))
;           mark=0.*fltarr(nc+1,nr)
;           mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
;           mark(nc,*)=mark(0,*)
;           contour,mark,x,alat,/overplot,levels=[0.1],thick=10,color=(float(k)/float(nth))*mcolor,c_label=stheta,c_charsize=1.5,c_charthick=2
;     endfor
      loadct,39
      contour,sp,x,alat,/overplot,levels=[35.],thick=9,color=.35*mcolor
      contour,sp,x,alat,/overplot,levels=[40.],thick=9,color=.6*mcolor
      contour,sp,x,alat,/overplot,levels=[45.],thick=9,color=.7*mcolor
      contour,sp,x,alat,/overplot,levels=[50.],thick=9,color=.75*mcolor
      contour,sp,x,alat,/overplot,levels=[55.],thick=9,color=.8*mcolor
      contour,sp,x,alat,/overplot,levels=[60.],thick=9,color=.85*mcolor
      contour,sp,x,alat,/overplot,levels=[66.],thick=9,color=.9*mcolor
      contour,sp,x,alat,/overplot,levels=[70.],thick=9,color=.95*mcolor
      MAP_SET,-90,0,90,/ortho,/noeras,/grid,/contin,/noborder,charsize=2.0,latdel=10,color=0,limit=[0.,-30,360.,-90]
;
loadct,0
;     ymnb=ymn -cbaryoff
;     ymxb=ymnb+cbarydel
;     set_viewport,xmn,xmx,ymnb,ymxb
;     !type=2^2+2^3+2^6
;     plot,[imin,imax],[0,0],yrange=[0,10],$
;           xrange=[imin,imax],xtitle=stheta+' K (30 km, 10 hPa) Potential Vorticity',/noeras,$
;           charsize=1.5,color=0,charthick=2
;     ybox=[0,10,10,0,0]
;     x1=imin
;     dx=(imax-imin)/float(nlvls)
;     for j=0,nlvls-1 do begin
;       xbox=[x1,x1,x1+dx,x1+dx,x1]
;       polyfill,xbox,ybox,color=col1(j)
;       x1=x1+dx
;     endfor

      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim Figures/polar_sh_zerowind_'+sdate+stime(itime)+'_LS.ps -rotate -90 Figures/polar_sh_zerowind_'+sdate+stime(itime)+'_LS.png'
         spawn,'rm -f Figures/polar_sh_zerowind_'+sdate+stime(itime)+'_LS.ps'
      endif

      endfor	; loop over times/day
goto,jump
end
