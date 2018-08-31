;
; base on .dat files
; plot daily zonal mean temperature and wind speeds
;
@stddat
@kgmt
@ckday
@kdate
@rd_nogaps_dat

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.2]
yorig=[0.3]
xlen=0.7
ylen=0.55
cbaryoff=0.065
cbarydel=0.02
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
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
dir='/aura7/harvey/NOGAPS_Alpha/Datfiles/NOGAPSA_'
dir='/aura7/harvey/NOGAPS_Alpha/Datfiles/'
slabs=['00','06','12','18']
slabs=['12']
ntimes=n_elements(slabs)

lstmn=12L & lstdy=2L & lstyr=2007L
ledmn=2L & leddy=28L & ledyr=2009L
lstday=0L & ledday=0L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
if lstyr lt 2006 then stop,'Year out of range '
if ledyr lt 2006 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
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
      if ndays gt ledday then stop,' Normal termination condition '
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
print,sdate
;
; loop over output times
;
      for itime=0L,ntimes-1L do begin

      sdate=syr+smn+sdy+slabs(itime)
      ifile='NOGAPSA_'+sdate+'_MetO_sabmls_aim9c.dat'
      dum=findfile(dir+ifile)
      if dum(0) eq '' then goto,jump
      rd_nogaps_dat,dir+ifile,iflg,nc,nr,nlv,alon,alat,p,g3d,t3d,u3d,v3d,$
         pv3d,o33d,h2o3d
print,'Read NOGAPS'
;
; zonal mean zonal wind
;
      tzm=fltarr(nr,nlv)
      uzm=fltarr(nr,nlv)
      for j=0L,nr-1L do begin
          for k=0L,nlv-1L do begin
              tzm(j,k)=total(t3d(*,j,k))/float(nc)
              uzm(j,k)=total(u3d(*,j,k))/float(nc)
          endfor
      endfor

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='yz_ubar_tbar_nogaps_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
       !p.thick=2.0
       !p.charsize=2.0
    endif
erase
loadct,39
xyouts,.45,.9,sdate,/normal,charsize=2,charthick=2,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
tlevel=120.+10.*findgen(nlvls)
index=where(tzm eq 0.)
if index(0) ne -1L then tzm(index)=0./0.
contour,tzm,alat,p,/noeras,xrange=[-90.,90.],yrange=[100,min(p)],/ylog,charsize=1.5,color=0,$
      ytitle='Pressure (hPa)',title='NOGAPS-ALPHA Tbar + Ubar',xticks=6,/cell_fill,c_color=col1,levels=tlevel
index=where(tlevel mod 10. eq 0)
contour,tzm,alat,p,levels=tlevel(index),color=0,/follow,/overplot,c_labels=1+0*index
slevel=-200.+10.*findgen(41)
index=where(uzm eq 0.)
if index(0) ne -1L then uzm(index)=0./0.
index=where(slevel gt 0.)
contour,uzm,alat,p,levels=slevel(index),color=0,/follow,/overplot,thick=5
index=where(slevel lt 0.)
contour,uzm,alat,p,levels=slevel(index),color=mcolor,/follow,/overplot,c_linestyle=5,thick=5


imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim yz_ubar_tbar_nogaps_'+sdate+'.ps -rotate -90 yz_ubar_tbar_nogaps_'+sdate+'.jpg'
       spawn,'/usr/bin/rm yz_ubar_tbar_nogaps_'+sdate+'.ps'
    endif
 
    endfor	; loop over output times
goto, jump
end
