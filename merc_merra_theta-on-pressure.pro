;
; plot mercator projection of MERRA theta on pressure
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=750
nydim=750
xorig=[0.20]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.07
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mdir='/Volumes/Data/MERRA_data/Datfiles/'
lstmn=6
lstdy=1
lstyr=2011
ledmn=8
leddy=1
ledyr=2011
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 70 then lstyr=lstyr+2000
if ledyr lt 70 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
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
icount=0

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; restore MERRA data
;
; LATITUDE_WACCM  FLOAT     = Array[96]
; LEVELS_ASSIM    DOUBLE    = Array[72]
; LONGITUDE_WACCM FLOAT     = Array[144]
; PSNEW           FLOAT     = Array[144, 96]
; TNEW            FLOAT     = Array[144, 96, 72]
; UNEW            FLOAT     = Array[144, 96, 72]
; VNEW            FLOAT     = Array[144, 96, 72]
;
      dum=file_search(mdir+'MERRA-on-WACCM_'+sdate+'.sav')
      if dum(0) ne '' then restore,mdir+'MERRA-on-WACCM_'+sdate+'.sav'
      if dum(0) eq '' then begin
         print,'missing data on '+sdate
         goto,jump
      endif
;
      if icount eq 0L then begin
         print,LEVELS_ASSIM
         rpress=163.66100
         read,' Enter pressure level ',rpress
         index=where(abs(rpress-LEVELS_ASSIM) eq min(abs(rpress-LEVELS_ASSIM)))
         ilev=index(0)
         slev=strcompress(rpress,/remove_all)+'hPa'
         nr=n_elements(latitude_waccm)
         nc=n_elements(longitude_waccm)
         nl=n_elements(levels_assim)
         alat=latitude_waccm
         alon=longitude_waccm
         pressure=levels_assim
      endif
      thnew=0.*tnew
      for k=0L,nl-1L do thnew(*,*,k)=tnew(*,*,k)*(1000./pressure(k))^.286
      th2d=reform(thnew(*,*,ilev))
      print,'min/max theta ',min(th2d),max(th2d)

      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                /bold,/color,bits_per_pixel=8,/helvetica,filename='merc_merra_theta-on-'+sdate+'_'+slev+'.ps'
         !p.charsize=1.25
         !p.thick=2
         !p.charthick=5
         !p.charthick=5
         !y.thick=2
         !x.thick=2
      endif
;
; plot 
;
      erase
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
      imax=500.
      imin=300.
      nlvls=20
      int=(imax-imin)/nlvls
      level=imin+int*findgen(nlvls+1)
      col1=1+indgen(nlvls+1)*mcolor/nlvls
      !psym=0
      contour,th2d,alon,alat,xrange=[0.,360.],yrange=[-90.,90.],levels=level,c_color=col1,/cell_fill,/noeras,title=sdate+'  ('+slev+')',color=0,$
              xtitle='Longitude',ytitle='Latitude',yticks=6,xticks=6
      contour,th2d,alon,alat,/overplot,levels=level,color=0,c_labels=1+intarr(nlvls)
      map_set,0,180,0,/contin,/grid,/noeras
      ymnb=min(yorig) -cbaryoff-0.01
      ymxb=ymnb+cbarydel
      set_viewport,min(xorig)+0.01,max(xorig)+xlen-0.01,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.25,charthick=2,$
           xtitle='MERRA Potential Temperature'
      ybox=[0,10,10,0,0]
      x2=imin
      dx=(imax-imin)/(float(nlvls)-1)
      for j=1,nlvls-1 do begin
          xbox=[x2,x2,x2+dx,x2+dx,x2]
          polyfill,xbox,ybox,color=col1(j)
          x2=x2+dx
      endfor
;
      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim merc_merra_theta-on-'+sdate+'_'+slev+'.ps -rotate -90 merc_merra_theta-on-'+sdate+'_'+slev+'.jpg'
      endif
      icount=icount+1L
goto,jump
end
