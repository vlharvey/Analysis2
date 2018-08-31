;
; GEOS-5 version
; plot polar projections and yz cross polar sections

@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
loadct,39
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=1
!p.charsize=1.5
!p.charthick=2
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.07,0.255,0.44,0.625,0.81,0.07,0.255,0.44,0.625,0.81,0.07,0.255,0.44,0.625,0.81]
yorig=[0.7,0.7,0.7,0.7,0.7,0.425,0.425,0.425,0.425,0.425,0.15,0.15,0.15,0.15,0.15]
npan=n_elements(xorig)
xlen=0.175
ylen=0.175
cbaryoff=0.02
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir2='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'

lstmn=1
lstdy=1
lstyr=2012
ledmn=2
leddy=25
ledyr=2012
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      GEOS Version '
;print, ' '
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
;
; Compute initial Julian date
;
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
      kcount=0L
;
; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; Read GEOS data
;
      dum=findfile(dir+sdate+'_AVG.V01.nc3')
      if dum(0) eq '' then dum=findfile(dir2+sdate+'_AVG.V01.nc3')
      if dum(0) ne '' then begin
         rd_geos5_nc3_meto,dum(0),nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
;        ncid=ncdf_open(dir+sdate+'_AVG.V01.nc4')
;        mark2=fltarr(nr,nc,nth)
;        ncdf_varget,ncid,3,mark2
;        ncdf_close,ncid
      endif
      print,sdate
;
; wind speed
;
      sp2=sqrt(u2^2.+v2^2.)
      index=where(u2 lt 0.)
      if index(0) ne -1L then sp2(index)=-1.*sp2(index)
;
; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
;
      t2=0.*p2
      z2=0.*p2
      for k=0,nth-1 do begin
          t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
          z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
      endfor
;
; declare arrays 
;
      if icount eq 0 then begin
         xyz=fltarr(nr,nth)
         yyz=fltarr(nr,nth)
         for i=0,nr-1 do yyz(i,*)=th
         for j=0,nth-1 do xyz(*,j)=alat 
      endif
;
; plot
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='yz_cross_pole_geos_nh_'+sdate+'_15pan_temp.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; loop over panels
;
    erase
    for ipan=0,nc-1,nc/npan do begin
         rlon1=alon(ipan)
         index1=where(alon eq rlon1)
         ilon1=index1(0)
         rlon2=rlon1+180.
         if rlon2 gt max(alon) then rlon2=rlon2-360.
         index2=where(alon eq rlon2)
         ilon2=index2(0)
         slon1=string(format='(f7.3)',rlon1)+'E'
         slon2=string(format='(f7.3)',rlon2)+'E'

         tyz=fltarr(nr,nth)
         markyz=fltarr(nr,nth)
         speedyz=fltarr(nr,nth)
         for k=0,nth-1 do begin
             tyz(0:nr/2-1,k)=t2(nr/2:nr-1,ilon1,k)
             tyz(nr/2:nr-1,k)=reverse(t2(nr/2:nr-1,ilon2,k))
             markyz(0:nr/2-1,k)=mark2(nr/2:nr-1,ilon1,k)
             markyz(nr/2:nr-1,k)=reverse(mark2(nr/2:nr-1,ilon2,k))
             speedyz(0:nr/2-1,k)=sp2(nr/2:nr-1,ilon1,k)
             speedyz(nr/2:nr-1,k)=reverse(sp2(nr/2:nr-1,ilon2,k))
         endfor

         xyouts,.35,.925,'GEOS-5  '+sdate,/normal,color=0,charsize=3
         !type=2^2+2^3
         xmn=xorig(kcount)
         xmx=xorig(kcount)+xlen
         ymn=yorig(kcount)
         ymx=yorig(kcount)+ylen
         set_viewport,xmn,xmx,ymn,ymx
         level=180.+5.*findgen(nlvls)
         ylab=' '
         yticklab=[' ',' ']
         ynticks=1
         if kcount eq 0 or kcount eq 5 or kcount eq 10 then ylab='Theta (K)'
         if kcount eq 0 or kcount eq 5 or kcount eq 10 then ynticks=4
         if kcount eq 0 or kcount eq 5 or kcount eq 10 then yticklab=['1000','1900','2800','3700','4600']
         contour,tyz,alat,th,levels=level,/cell_fill,c_color=col1,color=0,$
                 xtitle=slon1+'       '+slon2,ytitle=ylab,yticks=ynticks,ytickname=yticklab,$
                 xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6,yrange=[1000.,max(th)]
         contour,tyz,alat,th,levels=level,/overplot,/follow,color=0
         contour,speedyz,alat,th,levels=[30,50,70,90],/overplot,/follow,color=mcolor,thick=2,c_labels=['']
         contour,speedyz,alat,th,levels=[-90,-70,-50,-30],/overplot,/follow,color=mcolor,thick=2,c_labels=[''],c_linestyle=3
         contour,smooth(markyz,3),alat,th,levels=[0.1],color=0,thick=5,/overplot
         kcount=kcount+1
if kcount gt npan-1 then goto,jumpout
    endfor      ; loop over panels
    icount=icount+1
jumpout:
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff-0.05
    ymxb=ymnb+cbarydel
    set_viewport,min(xorig),max(xorig)+xlen,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle='Temperature (K)'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim yz_cross_pole_geos_nh_'+sdate+'_15pan_temp.ps -rotate -90 '+$
                            'yz_cross_pole_geos_nh_'+sdate+'_15pan_temp.jpg'
        spawn,'/usr/bin/rm yz_cross_pole_geos_nh_'+sdate+'_15pan_temp.ps'
     endif
goto,jump
end
