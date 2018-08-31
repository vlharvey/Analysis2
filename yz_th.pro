
; Read UKMO isentropic netCDF data and plot zonal mean
; temperature temperature as a function of log pressure.
; VLH 3/31/03

@rd_ukmo_nc3
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
!noeras=1
runtitle='!6UKMO Zonal Mean Temperature Tendency'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
nxdim=750
nydim=750
xorig=[0.1]
yorig=[0.15]
xlen=0.8
ylen=0.8
cbaryoff=0.08
cbarydel=0.01
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
   device,/landscape,bits=8,filename='ukmo_th.ps'
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
    y2d=fltarr(nr,nth)
    for k=0,nth-1 do y2d(*,k)=alat
;
; convert diabatic heating rate to temperature tendency
; T = theta* (p/po)^R/cp and divide by 1000 for km
;
    dtyz=fltarr(nr,nth)
    pyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
        for j=0,nr-1 do begin
            pave=moment(p2(j,*,k))
            pyz(j,k)=pave(0)
            qave=moment(q2(j,*,k))
            dtyz(j,k)=qave(0)*((pave(0)/1000.)^(.286))
        endfor
    endfor
;
; plot
;
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                     month(2),' ',15,', ',2001))
    mtitle=runtitle+' on '+date
    !type=2^2+2^3
    plot,[-90,90,90,-90,-90],[1000.,1000.,1.,1.,1000.],/ylog,$
          xrange=[-90,90] ,yrange=[1000.,1.],xtitle='!6Latitude',$
          ytitle='!6Pressure (mb)',xticks=6,title=mtitle
    level=[-7.,-6.,-5.,-4.,-3.,-2.,-1.,-0.75,-0.5,-0.25,-0.1,0.,$
           0.1,0.25,0.5,0.75,1.]
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls

    contour,dtyz,y2d,pyz,levels=level,/overplot,/fill,/cell_fill,c_color=col1
    contour,dtyz,y2d,pyz,levels=level,c_charsize=1.0,$
           /follow,/overplot,c_linestyle=(level lt 0.0),c_color=0
    contour,dtyz,y2d,pyz,levels=[0.],c_charsize=1.0,$
           /follow,/overplot,thick=3,c_color=0

      imin=min(level)
      imax=max(level)
      ymnb=yorig -cbaryoff
      ymxb=ymnb  +cbarydel
      set_viewport,xmn,xmx,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax]
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(icmm1)
      for j=1,icmm1 do begin
            xbox=[x1,x1,x1+dx,x1+dx,x1]
            polyfill,xbox,ybox,color=j
            x1=x1+dx
      endfor
stop
      if setplot eq 'ps' then device, /close
      jump:
endfor
end
