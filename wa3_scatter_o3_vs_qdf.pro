@stddat
@kgmt
@ckday
@kdate
@rd_waccm3_nc3

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
!NOERAS=-1
device,decompose=0
nxdim=700
nydim=700
xorig=[0.1]
yorig=[0.35]
xlen=0.8
ylen=0.8
cbaryoff=0.04
cbarydel=0.02
!NOERAS=-1
lstmn=1
lstdy=31
lstyr=1990
ledmn=1
leddy=4
ledyr=2004
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '     WACCM Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 57 then lstyr=lstyr+2000
if ledyr lt 57 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1957 then stop,'Year out of range '
if ledyr lt 1957 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

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
dir='/aura7/harvey/WACCM_data/Datfiles/TNV3_files/wa3_tnv3_'

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
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
      ifile=string(FORMAT='(i4.4,i2.2,i2.2,a4)',iyr,imn,idy,'.nc3')
      rd_waccm3_nc3,dir+ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,sf2,o32,ch42,no22,h2o2,iflag
      if iflag eq 1 then goto,jump

; select theta level
    if icount eq 0L then begin
       rlev=1000.
;      print,th
;      read,'Enter theta surface ',rlev
       zindex=where(th eq rlev)
       ilev=zindex(0)
       slev=strcompress(string(fix(th(ilev))),/remove_all)+'K'
       icount=1L
    endif

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='Figures/'+ifile+'_wa3_scatter_o3_vs_qdf_'+slev+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

    mark1=transpose(mark2(*,*,ilev))
    sf1=transpose(sf2(*,*,ilev))
    qdf1=transpose(qdf2(*,*,ilev))
    o31=transpose(o32(*,*,ilev))
    sf=0.*fltarr(nc+1,nr)
    sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
    sf(nc,*)=sf(0,*)
    mark=0.*fltarr(nc+1,nr)
    mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
    mark(nc,*)=mark(0,*)
    o3=0.*fltarr(nc+1,nr)
    o3(0:nc-1,0:nr-1)=o31(0:nc-1,0:nr-1)
    o3(nc,*)=o3(0,*)
    qdf=0.*fltarr(nc+1,nr)
    qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
    qdf(nc,*)=qdf(0,*)
    x=fltarr(nc+1)
    x(0:nc-1)=alon
    x(nc)=alon(0)+360.
x2d=0.*qdf
y2d=0.*qdf
for i=0L,nc do y2d(i,*)=alat
for j=0L,nr-1 do x2d(*,j)=x
    erase
    xyouts,.35,.95,'WACCM 3: '+strmid(ifile,0,8),charsize=2,/normal,color=0
    xmn=0.1
    xmx=0.45
    ymn=0.55
    ymx=0.9
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    MAP_SET,0,180,0,/noeras,/grid,/contin,title='Ozone + QDF '+slev,charsize=1.5,color=0,$
            limit=[0,0,90,360],label=1
    nlvls=19
    level=1.0+0.5*findgen(nlvls)
    qlevel=-300.+20.*findgen(31)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,o3,x,alat,/overplot,levels=level,c_color=col1,$
           /cell_fill,/noeras
    contour,o3,x,alat,/overplot,levels=[7],/follow,c_labels=0*level,/noeras,color=0,thick=4
;   contour,sf,x,alat,/overplot,nlevels=30,color=0
loadct,0
    contour,mark,x,alat,/overplot,levels=[0.1],thick=10,color=150
    contour,mark,x,alat,/overplot,levels=[-0.1],thick=10,color=0
loadct,38
    index=where(qlevel lt 0.)
    contour,qdf,x,alat,levels=qlevel(index),/follow,/overplot,color=mcolor
    index=where(qlevel gt 0.)
    contour,qdf,x,alat,levels=qlevel(index),/follow,/overplot,color=0
    MAP_SET,0,180,0,/noeras,/grid,/contin,charsize=1.5,color=0,limit=[0,0,90,360],label=1
    rlat=51.25
    oplot,findgen(361),rlat+0.*findgen(361),color=0,psym=8,symsize=0.2
    imin=min(level)
    imax=max(level)
    xmnb=xmx+cbaryoff
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='(ppmv)',charsize=1.5,color=0
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=0.6
    xmx=0.95
    ymn=0.55
    ymx=0.9
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    plot,findgen(10),findgen(10),/nodata,color=0,xrange=[3.,8.],yrange=[-200.,200.]
    index=where(o3 le 7. and mark ne 1. and y2d gt 0.,ngood)
    ymin=min(y2d(index))-1 & ymax=max(y2d(index))+1
    for i=0L,ngood-1L do begin
    oplot,[o3(index(i)),o3(index(i))],[qdf(index(i)),qdf(index(i))],psym=8,$
           color=((y2d(index(i))-ymin)/(ymax-ymin))*mcolor,symsize=0.75
    endfor
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a)
    index=where(o3 le 7. and mark lt 0. and y2d gt 0.)
    if index(0) ne -1L then oplot,o3(index),qdf(index),psym=8,color=0,symsize=0.75
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
symin=strcompress(long(ymin),/remove_all)
symax=strcompress(long(ymax),/remove_all)
    imin=ymin
    imax=ymax
    xmnb=xmx
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='Lat    ',charsize=1.5,color=0
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
if j eq 0L then xyouts,xbox(0),ybox(0),symin,charsize=1.25,color=0,/data,charthick=2
if j eq nlvls-1L then xyouts,xbox(0),ybox(0),symax,charsize=1.25,color=0,/data,charthick=2
        y1=y1+dy
    endfor

    zindex=where(rlat eq alat)
    ilat=zindex(0)
    slat=strcompress(string(fix(alat(ilat))),/remove_all)
    markxz=0.*fltarr(nc+1,nth)
    markxz(0:nc-1,0:nth-1)=reform(mark2(ilat,0:nc-1,0:nth-1))
    markxz(nc,*)=markxz(0,*)
    o3xz=0.*fltarr(nc+1,nth)
    o3xz(0:nc-1,0:nth-1)=reform(o32(ilat,0:nc-1,0:nth-1))
    o3xz(nc,*)=o3xz(0,*)
    qdfxz=0.*fltarr(nc+1,nth)
    qdfxz(0:nc-1,0:nth-1)=reform(qdf2(ilat,0:nc-1,0:nth-1))
    qdfxz(nc,*)=qdfxz(0,*)
    xmn=0.1
    xmx=0.45
    ymn=0.1
    ymx=0.45
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    contour,o3xz,x,th,levels=level,/cell_fill,title='Ozone + QDF '+slat+' N',c_color=col1,$
            min_value=0.,xticks=4,xrange=[0.,360.],yrange=[500.,1800.],$
            ytitle='Theta',xtitle='Longitude',color=0,charsize=1.5
    contour,o3xz,x,th,levels=findgen(nlvls),/follow,/overplot,color=0,min_value=0.
    contour,o3xz,x,th,levels=[3],/follow,/overplot,color=0,min_value=0.,thick=3
loadct,0
    contour,smooth(markxz,3),x,th,levels=[0.1],/follow,/overplot,color=150,thick=10
    contour,smooth(markxz,3),x,th,levels=[-0.1],/follow,/overplot,color=0,thick=10
loadct,38
    index=where(qlevel lt 0.)
    contour,qdfxz,x,th,levels=qlevel(index),/follow,/overplot,color=mcolor
    index=where(qlevel gt 0.)
    contour,qdfxz,x,th,levels=qlevel(index),/follow,/overplot,color=0
;   contour,qdfxz,x,th,levels=[0.],/follow,/overplot,color=0,thick=3
    imin=min(level)
    imax=max(level)
    xmnb=xmx+cbaryoff
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='(ppmv)',charsize=1.5,color=0
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    xmn=0.6
    xmx=0.95
    ymn=0.1
    ymx=0.45
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    plot,findgen(10),findgen(10),/nodata,color=0,xrange=[3.,8.],yrange=[-200.,200.]
    th2d=0.*qdfxz
    for i=0L,nc do th2d(i,*)=th
    index=where(o3xz ge 3. and markxz ne 1. and th2d ge 600. and th2d le 1600.,ngood)
    thmin=550. & thmax=1650.
    for i=0L,ngood-1L do begin
    oplot,[o3xz(index(i)),o3xz(index(i))],[qdfxz(index(i)),qdfxz(index(i))],psym=8,$
           color=((thmax-th2d(index(i)))/(thmax-thmin))*mcolor,symsize=0.75
    endfor
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a)
    index=where(o3xz ge 3. and markxz lt 0. and th2d ge 600. and th2d le 1600.,ngood)
    if index(0) ne -1L then oplot,o3xz(index),qdfxz(index),psym=8,color=0,symsize=0.75
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
    imin=thmin
    imax=thmax
    xmnb=xmx
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='Theta    ',charsize=1.5,color=0
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim Figures/'+ifile+'_wa3_scatter_o3_vs_qdf_'+slev+'.ps -rotate -90'+$
             ' Figures/'+ifile+'_wa3_scatter_o3_vs_qdf_'+slev+'.jpg'
    endif
goto, jump
end
