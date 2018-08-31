
; yz zonal mean t and q

@rd_ukmo_nc3

a=findgen(8)*(2*!pi/8.)
usersym,.5*cos(a),.5*sin(a),/fill
loadct,38
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=31
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
; define viewport location
nxdim=750
nydim=750
xorig=[0.20,0.20]
yorig=[0.60,0.15]
xlen=0.6
ylen=0.3
cbaryoff=0.08
cbarydel=0.01
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'yz_cross_pole.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    print,ifile
    if n eq 0 then ifile0=ifile
    iflag=0
    rd_ukmo_nc3,dir+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump

; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
    t2=0.*p2
    z2=0.*p2
    for k=0,nth-1 do begin
        t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
        z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
    endfor

    xyz=fltarr(nr,nth)
    yyz=fltarr(nr,nth)
    for i=0,nr-1 do yyz(i,*)=th
    for j=0,nth-1 do xyz(*,j)=alat 
;
; zonal mean
    uyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    zyz=fltarr(nr,nth)
    qyz=fltarr(nr,nth)
    markyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
    for j=0,nr-1 do begin
        uyz(j,k)=total(u2(j,*,k))/float(nc)
        tyz(j,k)=total(t2(j,*,k))/float(nc)
        zyz(j,k)=total(z2(j,*,k))/float(nc)
        qyz(j,k)=total(q2(j,*,k))/float(nc)
        markyz(j,k)=total(marksf2(j,*,k))/float(nc)
    endfor
    endfor

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='Figures/zonal_mean_qt_'+ifile+'.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

; plot
    !noeras=1
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    xyouts,.4,.95,ifile,/normal,charsize=3
    mtitle='Zonal Mean Temperature'
    !type=2^2+2^3
    level=170.+4.*findgen(nlvls)
    contour,tyz,alat,th,levels=level,/fill,/cell_fill,c_color=col1,ytitle='Theta',$
            xtitle='Latitude',title=mtitle,charsize=1.5
    contour,tyz,alat,th,levels=level,/follow,c_color=mcolor,/overplot
    contour,tyz,alat,th,levels=[170.,180.,190.,200.],/overplot,/follow,c_color=[0],thick=2
    contour,markyz,alat,th,levels=[0.01],/overplot,/follow,c_color=[mcolor],thick=3,c_label=[0]
    contour,markyz,alat,th,levels=[-0.01],/overplot,/follow,c_color=[0],thick=3,c_label=[0]
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    level=-100.+5.*findgen(nlvls)
    print,min(qyz),max(qyz)
    mtitle='Zonal Theta Dot'
    contour,qyz,alat,th,levels=level,/fill,/cell_fill,c_color=col1,ytitle='Theta',$
            xtitle='Latitude',title=mtitle,charsize=1.5
    contour,qyz,alat,th,levels=-100.+5.*findgen(nlvls),/overplot,/follow,c_color=[0]
    contour,qyz,alat,th,levels=5.+5.*findgen(nlvls),/overplot,/follow,c_color=[mcolor]
    contour,markyz,alat,th,levels=[0.01],/overplot,/follow,c_color=[mcolor],thick=3,c_label=[0]
    contour,markyz,alat,th,levels=[-0.01],/overplot,/follow,c_color=[0],thick=3,c_label=[0]
    imin=min(level)
    imax=max(level)
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

; Close PostScript file and return control to X-windows
    if setplot eq 'ps' then device, /close

    stop
    jump:
endfor		; loop over days
end
