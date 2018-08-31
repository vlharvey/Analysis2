
; Read UKMO isentropic netCDF data and plot zonal mean
; temperature as a function of log pressure.
; VLH 3/31/03

@rd_ukmo_nc3
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
!noeras=1
runtitle='!6UKMO Zonal Mean TTL'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
nxdim=750
nydim=750
xorig=[0.1]
yorig=[0.25]
xlen=0.8
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
   device,/landscape,bits=8,filename='yz_th_ttl.ps'
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
;
; T = theta* (p/po)^R/cp and divide by 1000 for km
;
    t2=0.*pv2
    z2=0.*pv2
    for k=0,nth-1 do begin
        t2(0:nr-1,0:nc-1,k)=th(k)*((p2(0:nr-1,0:nc-1,k)/1000.)^(.286) )
        z2(0:nr-1,0:nc-1,k)=(msf2(0:nr-1,0:nc-1,k)-1004.* $
                               t2(0:nr-1,0:nc-1,k))/(9.86*1000.)
    endfor

    pyz=fltarr(nr,nth)
    tyz=fltarr(nr,nth)
    for k=0,nth-1 do begin
        for j=0,nr-1 do begin
            pave=moment(p2(j,*,k))
            pyz(j,k)=pave(0)
            tave=moment(t2(j,*,k))
            tyz(j,k)=th(k)*((pave(0)/1000.)^(.286))
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
    plot,[-45,45,45,-45,-45],[330.,330.,400.,400.,330.],$
          xrange=[-45,45] ,yrange=[330.,400.],xtitle='!6Latitude',$
          ytitle='!6Pressure (mb)',xticks=6,title=mtitle
    level=180.+4.*findgen(21)
    nlvls=n_elements(level)
    col1=1+indgen(nlvls)*icolmax/nlvls

    contour,tyz,alat,th,levels=level,/overplot,/fill,/cell_fill,c_color=col1
    contour,tyz,alat,th,levels=level,c_charsize=1.0,$
           /follow,/overplot,c_linestyle=(level lt 0.0),c_color=0
    contour,pyz,alat,th,levels=[150.],c_charsize=1.0,$
           /follow,/overplot,thick=3
ttop=fltarr(nr)
ttl=fltarr(nr,nth)
for j=0,nr-1 do begin
    tmp=reform(tyz(j,*))
    index=where(tmp eq min(tmp))
    ttop(j)=th(index(0))
index=where(th le th(index(0)))
    ttl(j,index)=th(index)
endfor
oplot,alat,ttop,psym=0,thick=3
index=where(pyz le 150. and ttl gt 0.)
oplot,x2d(index),y2d(index),psym=2

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
