
; Plot profiles of PV=0 at each longitude.
; VLH 8/31/03

@rd_ukmo_nc3
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
mcolor=icolmax
a=findgen(4)*(2*!pi/4.)
usersym,cos(a),sin(a),/fill
!noeras=1
runtitle='!6UKMO PV=0'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.6
cbaryoff=0.08
cbarydel=0.02
setplot='x'
read,'SETPLOT ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_th_pv.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
   xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'yz_th.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    print,ifile
    iflag=0
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump
    x2d=fltarr(nr,nth)
    y2d=fltarr(nr,nth)
    for k=0,nth-1 do x2d(*,k)=alat
    for j=0,nr-1 do y2d(j,*)=th
    t2=0.*pv2
    z2=0.*pv2
    for k=0,nth-1 do begin
        t2(0:nr-1,0:nc-1,k)=th(k)*((p2(0:nr-1,0:nc-1,k)/1000.)^(.286) )
        z2(0:nr-1,0:nc-1,k)=(msf2(0:nr-1,0:nc-1,k)-1004.* $
                               t2(0:nr-1,0:nc-1,k))/(9.86*1000.)
    endfor
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    date=ifile
    mtitle=runtitle+' on '+date
    !type=2^2+2^3
thmax=2000.
    plot,[-15,65,65,-15,-15],[330.,330.,thmax,thmax,330.],$
          xrange=[-15,65] ,yrange=[330.,thmax],xtitle='!6Latitude',$
          ytitle='!6Pressure (mb)',xticks=6,title=mtitle
    plots,0.,330.
    plots,0.,thmax,/continue,thick=3
    pyz=fltarr(nr,nth)
    pvyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    for i=0,nc-1 do begin
        pyz=reform(p2(*,i,*),nr,nth)
        pvyz=reform(pv2(*,i,*),nr,nth)
        tyz=reform(t2(*,i,*),nr,nth)
        contour,pvyz,alat,th,levels=[0.],c_labels=[0],/follow,/overplot,$
                color=((i+1.)/(nc+1.))*mcolor,thick=2
;       index=where(pvyz*x2d lt 0.)
;       if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=8,color=((i+1.)/(nc+1.))*mcolor
    endfor
; Draw color bar
    imin=0.
    imax=360.
    ymnb=yorig -cbaryoff
    ymxb=ymnb  +cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],$
         xrange=[imin,imax],xtitle='!6Longitude'
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nc-1)
    for i=0,nc-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=((i+1.)/(nc+1.))*mcolor
        x1=x1+dx
    endfor
;   stop
    if setplot eq 'ps' then device, /close
    if setplot eq 'x' then begin
       save=assoc(3,bytarr(nxdim,nydim))
       img=bytarr(nxdim,nydim)
       img(0,0)=TVRD(0,0,nxdim,nydim)
       write_gif,ifile+'_yz_pv.gif',img
    endif

    jump:
endfor
end
