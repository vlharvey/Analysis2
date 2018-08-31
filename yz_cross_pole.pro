
; plot polar projection and yz cross polar section

@rd_ukmo_nc3

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

loadct,38
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
; define viewport location
nxdim=750
nydim=750
xorig=[0.30,0.30]
yorig=[0.55,0.13]
xlen=0.4
ylen=0.4
cbaryoff=0.06
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

    if n eq 0 then begin
;   print,th
    th0=1200.
;   read,' Enter desired theta surface ',th0
    thlev=where(th eq th0)
    thlev=thlev(0)
    sth=strcompress(string(fix(th0)),/remove_all)
;   print,alon
    rlon1=181.875
;   read,' Enter longitude ',rlon1
    index1=where(alon eq rlon1)
    if index1(0) eq -1 then stop,'Bad Longitude'
    ilon1=index1(0)
    rlon2=rlon1+180.
    if rlon2 gt max(alon) then rlon2=rlon2-360.
    index2=where(alon eq rlon2)
    ilon2=index2(0)
    slon1=' '+strcompress(string(rlon1),/remove_all)+' E'
    slon2=' '+strcompress(string(rlon2),/remove_all)+' E'
    print,'longitudes ',rlon1,rlon2
    x=fltarr(nc+1)
    x(0:nc-1)=alon(0:nc-1)
    x(nc)=alon(0)+360.
    x2d=fltarr(nc+1,nr)
    y2d=fltarr(nc+1,nr)
    for i=0,nc do y2d(i,*)=alat
    for j=0,nr-1 do x2d(*,j)=x
    xyz=fltarr(nr,nth)
    yyz=fltarr(nr,nth)
    for i=0,nr-1 do yyz(i,*)=th
    for j=0,nth-1 do xyz(*,j)=alat 
    endif
    u1=transpose(u2(*,*,thlev))
    t1=transpose(t2(*,*,thlev))
    z1=transpose(z2(*,*,thlev))
    l1=transpose(marksf2(*,*,thlev))
    h1=transpose(marksf2(*,*,thlev))
    u=fltarr(nc+1,nr)
    u(0:nc-1,0:nr-1)=u1(0:nc-1,0:nr-1)
    u(nc,*)=u(0,*)
    t=fltarr(nc+1,nr)
    t(0:nc-1,0:nr-1)=t1(0:nc-1,0:nr-1)
    t(nc,*)=t(0,*)
    z=fltarr(nc+1,nr)
    z(0:nc-1,0:nr-1)=z1(0:nc-1,0:nr-1)
    z(nc,*)=z(0,*)
    l=fltarr(nc+1,nr)
    l(0:nc-1,0:nr-1)=l1(0:nc-1,0:nr-1)
    l(nc,*)=l(0,*)
    h=fltarr(nc+1,nr)
    h(0:nc-1,0:nr-1)=h1(0:nc-1,0:nr-1)
    h(nc,*)=h(0,*)
    uyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    zyz=fltarr(nr,nth)
    lyz=fltarr(nr,nth)
    hyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
        uyz(0:nr/2-1,k)=u2(nr/2:nr-1,ilon1,k)
        uyz(nr/2:nr-1,k)=reverse(u2(nr/2:nr-1,ilon2,k))
        tyz(0:nr/2-1,k)=t2(nr/2:nr-1,ilon1,k)
        tyz(nr/2:nr-1,k)=reverse(t2(nr/2:nr-1,ilon2,k))
        zyz(0:nr/2-1,k)=z2(nr/2:nr-1,ilon1,k)
        zyz(nr/2:nr-1,k)=reverse(z2(nr/2:nr-1,ilon2,k))
        lyz(0:nr/2-1,k)=marksf2(nr/2:nr-1,ilon1,k)
        lyz(nr/2:nr-1,k)=reverse(marksf2(nr/2:nr-1,ilon2,k))
        hyz(0:nr/2-1,k)=marksf2(nr/2:nr-1,ilon1,k)
        hyz(nr/2:nr-1,k)=reverse(marksf2(nr/2:nr-1,ilon2,k))
    endfor

    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       device,/color,/landscape,bits=8,filename='Figures/ukmo_'+ifile+'_'+sth+'_'+slon1+'.ps'
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

; plot
    !noeras=1
    !type=2^2+2^3
    !p.thick=1
    erase
    ipan=0
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    mtitle='!6'+sth+' K UKMO Zonal Wind Temperature and Vortices on '+ifile
    !type=2^2+2^3
    level=-100.+10.*findgen(nlvls)
    MAP_SET,90,0,-90,/stereo,/noeras,title=mtitle,color=lc
    contour,u,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
    index=where(level gt 0.)
    contour,u,x,alat,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,u,x,alat,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
    tlevel=190.+5.*findgen(nlvls)
    for ii=0,nlvls-1L do $
    contour,t,x,alat,levels=tlevel(ii),/overplot,/follow,c_charsize=1.5,$
            c_labels=[1],c_color=col1(ii),thick=3
    loadct,0
    contour,l,x,alat,levels=[0.1],color=icolmax*.5,thick=8,/overplot
    contour,h,x,alat,levels=[-0.1],color=icolmax*.9,thick=8,/overplot
    loadct,38
    MAP_SET,90,0,-90,/stereo,/grid,/contin,/noeras,/noborder,color=0
    oplot,rlon1+0.*alat,alat,thick=2,color=icolmax
    oplot,rlon2+0.*alat,alat,thick=2,color=icolmax
    ipan=1
    xmn=xorig(ipan)
    xmx=xorig(ipan)+xlen
    ymn=yorig(ipan)
    ymx=yorig(ipan)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    contour,uyz,alat,th,levels=level,/fill,/cell_fill,c_color=col1,$
            xtitle='!6'+slon1+'       Latitude      '+slon2,ytitle='!6Theta',$
            xtickname=['Eq','30','60','NP','60','30','Eq'],xticks=6
    index=where(level gt 0.)
    contour,uyz,alat,th,levels=level(index),/overplot,/follow,color=0
    index=where(level lt 0.)
    contour,uyz,alat,th,levels=level(index),/overplot,/follow,color=icolmax,c_linestyle=5
    for ii=0,nlvls-1L do $
    contour,tyz,alat,th,levels=tlevel(ii),/follow,c_color=col1(ii),/overplot,$
            c_labels=[1],thick=3
    loadct,0
    contour,lyz,alat,th,levels=[0.1],color=icolmax*.5,thick=8,/overplot
    contour,hyz,alat,th,levels=[-0.1],color=icolmax*.9,thick=8,/overplot
    loadct,38
;   oplot,alat,th0+0.*alat,thick=3,color=icolmax
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
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then device, /close

;    stop
    jump:
endfor		; loop over days
end
