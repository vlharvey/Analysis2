;
; contour theta on pressure

@stddat
@kgmt
@ckday
@kdate
@rd_ukmo
@date2uars
@drawvectors

device,decompose=0

runtitle='UKMO '
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
setplot='x'
read,'setplot ',setplot
nlg=0l
nlat=0l
nlv=0l
lstmn=0
lstdy=0
lstyr=0
ledmn=0
leddy=0
ledyr=0
lsfc=0
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      UKMO Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

; define viewport location 
nxdim=750
nydim=750
xorig=0.15
yorig=0.15
xlen=0.7
ylen=0.7
cbaryoff=0.1
cbarydel=0.02

!p.thick=1
!p.charsize=1.0

if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

; set color table
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
!noeras=1

; Compute initial Julian date
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
      if ndays gt ledday then stop,' Normal termination condition '

;***Read UKMO data
      file='/aura3/data/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr+'_m'+imn+'_d'+idy+'_h12.pp.dat')
      print,file
      rd_ukmo,file,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
              zp,tp,up,vp
      if iflg ne 0 then goto, jump

      if icount eq 0L then begin 
         rsfc=0.0
         print,p
         read,' Enter desired pressure surface ',rsfc
         index=where(fix(100.*rsfc) eq fix(100.*p))
         lsfc=index(0)
         ssfc=strcompress(string(rsfc),/remove_all)
         ssfc=strmid(ssfc,0,7)
      endif

      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename=$
             (string(FORMAT='(a5,i2.2,a1,i2.2,a1,i4,a1,a7,a5)',$
             'ukmo_',imn,'_',idy,'_',iyr,'_',ssfc,'mb.ps'))
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.charsize=1.0
      endif

; Declare plotting arrays
      z2=fltarr(nlg+1,nlat)
      t2=fltarr(nlg+1,nlat)
      x=fltarr(nlg+1)
      x(0:nlg-1)=alon(0:nlg-1)
      x(nlg)=x(0)+360.
;
; shift longitudes for plotting
;
      z2(0:nlg-1,0:nlat-1)=zp(0:nlg-1,0:nlat-1,lsfc)/1000.
      z2(nlg,*)=z2(0,*)
      t2(0:nlg-1,0:nlat-1)=tp(0:nlg-1,0:nlat-1,lsfc)
      t2(nlg,*)=t2(0,*)
      th2=t2*(1000./rsfc)^.286

; Find data range for autoscaling contours
      minval=min(z2)
      maxval=max(z2)

; Autoscale if scale values for parameter/level are not defined
      nlvls=20
      level=minval+((maxval-minval)/nlvls)*findgen(nlvls)
      col1=1+indgen(nlvls)*icolmax/nlvls
      mb=strcompress(string(FORMAT='(F7.2,A4)',p(lsfc),' mb '))
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
      mtitle=runtitle+'Theta at'+mb+'on '+date

; Mercator projection
      xmn=xorig
      xmx=xorig+xlen
      ymn=yorig
      ymx=yorig+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !psym=0
      erase
      MAP_SET,90,180,0,/stereo,/contin,/grid,/noeras,title='!6UKMO '+mb+$
             'Geopotential Height on '+date
      contour,z2,x,alat,levels=level,/cell_fill,/noeras,$
              title='!6'+mtitle,c_color=col1,charsize=1.5,$
              xticks=6,yticks=6,xtitle='!6Longitude',$
              ytitle='!6Latitude',/overplot
      contour,z2,x,alat,levels=level,/overplot,c_color=[0],$
              c_labels=1+0*level,c_linestyle=level lt 0,/noeras,$
              charsize=3.,c_annotation=string(fix(level))
      MAP_SET,90,180,0,/stereo,/contin,/grid,/noeras,color=0

      u=fltarr(nlg,nlat-1)
      u(0:nlg-1,0:nlat-2)=up(0:nlg-1,0:nlat-2,lsfc)
      u(nlg-1,*)=u(0,*)
      v=fltarr(nlg,nlat-1)
      v(0:nlg-1,0:nlat-2)=vp(0:nlg-1,0:nlat-2,lsfc)
      v(nlg-1,*)=v(0,*)
      drawvectors,nlg,nlat-1,wlon,wlat,u,v,5,0

; Draw color bar
      imin=min(level)
      imax=max(level)
      ymnb=yorig -cbaryoff
      ymxb=ymnb  +cbarydel
      set_viewport,xmn,xmx,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],$
            xrange=[imin,imax],xtitle='!6Theta (K)'
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(icmm1)
      for j=1,icmm1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=j
          x1=x1+dx
      endfor
;stop
    if setplot eq 'x' then begin
       save=assoc(3,bytarr(nxdim,nydim))
       img=bytarr(nxdim,nydim)
       img(0,0)=TVRD(0,0,nxdim,nydim)
       write_gif,(string(FORMAT='(a5,i2.2,a1,i2.2,a1,i4,a1,a7,a6)',$
                 'ukmo_',imn,'_',idy,'_',iyr,'_',ssfc,'mb.gif')),img
    endif

; Close PostScript file and return control to X-windows
      if setplot eq 'ps' then begin
         device, /close
         set_plot, 'x'
         !p.font=0
         !p.thick=1.0
      endif

      icount=icount+1L
goto, jump

end
