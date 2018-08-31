;
; reads in .nc3 and .nc4 marker field and plots mercator plots
;
@rd_ukmo_nc3

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.15,0.575]
xlen=0.7
ylen=0.3
cbaryoff=0.03
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'polar_pvgrad+mark.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    print,ifile
    iflag=0
    rd_ukmo_nc3,dir+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump
pv2=pv2*1.e6	; convert to PVU
;
; meridional PV gradient calculation
;
pvgrad2=0*pv2
pvgrad4=0*pv2
dy2=2.*(alat(1)-alat(0))*!pi/180.
dy4=12.*(alat(1)-alat(0))*!pi/180.
for k=0,nth-1 do begin
for j=2,nr-3 do begin
    jm1=j-1
    jp1=j+1
    jm2=j-2
    jp2=j+2
    for i=0,nc-1 do begin
        pvgrad2(j,i,k) = (pv2(jp1,i,k)-pv2(jm1,i,k))/dy2	; 2nd order
        pvgrad4(j,i,k) = (-1.*pv2(jp2,i,k)+8.*pv2(jp1,i,k) $
                          -8.*pv2(jm1,i,k)+1.*pv2(jm2,i,k))/dy4	; 4th order
    endfor
endfor
pvgrad2(0,*,k)=pvgrad2(2,*,k)
pvgrad2(1,*,k)=pvgrad2(2,*,k)
pvgrad4(nr-2,*,k)=pvgrad2(nr-3,*,k)
pvgrad4(nr-1,*,k)=pvgrad2(nr-3,*,k)
endfor

if n eq 0 then begin
theta=0.
;print,th
;read,'Enter theta ',theta
theta=1600.
index=where(theta eq th)
if index(0) eq -1 then stop,'Invalid theta level '
thlev=index(0)
endif
        stheta=strcompress(string(fix(theta)),/remove_all)
        qdf1=transpose(qdf2(*,*,thlev))
        sf1=transpose(sf2(*,*,thlev))
        pv1=transpose(pv2(*,*,thlev))
        pvgrad1=transpose(pvgrad4(*,*,thlev))
        marksfl1=transpose(marksf2(*,*,thlev))
        marksfh1=transpose(marksf2(*,*,thlev))
        qdf=0.*fltarr(nc+1,nr)
        qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
        qdf(nc,*)=qdf(0,*)
        sf=0.*fltarr(nc+1,nr)
        sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
        sf(nc,*)=sf(0,*)
        pv=0.*fltarr(nc+1,nr)
        pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
        pv(nc,*)=pv(0,*)
        pvgrad=0.*fltarr(nc+1,nr)
        pvgrad(0:nc-1,0:nr-1)=pvgrad1(0:nc-1,0:nr-1)
        pvgrad(nc,*)=pvgrad(0,*)
        marksfl=0.*fltarr(nc+1,nr)
        marksfl(0:nc-1,0:nr-1)=marksfl1(0:nc-1,0:nr-1)
        marksfl(nc,*)=marksfl(0,*)
        marksfh=0.*fltarr(nc+1,nr)
        marksfh(0:nc-1,0:nr-1)=marksfh1(0:nc-1,0:nr-1)
        marksfh(nc,*)=marksfh(0,*)
        x=fltarr(nc+1)
        x(0:nc-1)=alon
        x(nc)=alon(0)+360.
        lon=0.*sf
        lat=0.*sf
        for i=0,nc   do lat(i,*)=alat
        for j=0,nr-1 do lon(*,j)=x
index=where(lat gt 85.)
pv(index)=0./0.
pvgrad(index)=0./0.

        if setplot eq 'ps' then begin
           set_plot,'ps'
           xsize=nxdim/100.
           ysize=nydim/100.
         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                /bold,/color,bits_per_pixel=8,filename='merc_'+ifile+'_'+stheta+'K_pvgrad+mark.ps'
!p.font=1
!p.charsize=1.75
!p.charthick=1.5
        endif

; Set plot boundaries
        erase
        !type=2^2+2^3
        xmn=xorig(0)
        xmx=xorig(0)+xlen
        ymn=yorig(0)
        ymx=yorig(0)+ylen
        set_viewport,xmn,xmx,ymn,ymx
;       xyouts,.35,.95,'13 February 2002 '+stheta+' K',/normal,color=0
        MAP_SET,0,0,0,/noeras,/contin,/noborder,color=0,limit=[40.,0.,90.,360.]
        index=where(lat gt 0.)
        sfmin=min(pvgrad(index))
        sfmax=max(pvgrad(index))
sfmax=0.05*1.e6
sfmin=-0.05*1.e6
        nlvls=20
        sfint=(sfmax-sfmin)/nlvls
        sflevel=sfmin+sfint*findgen(nlvls)
        col1=1+indgen(nlvls)*mcolor/nlvls
        !psym=0
        contour,pvgrad,x,alat,/overplot,levels=sflevel,c_color=col1,/cell_fill,/noeras
        index=where(sflevel gt 0.)
        contour,pvgrad,x,alat,/overplot,levels=sflevel(index),color=0
        index=where(sflevel lt 0.)
        contour,pvgrad,x,alat,/overplot,levels=sflevel(index),color=mcolor,c_linestyle=5
        contour,marksfl,x,alat,/overplot,levels=[0.1],thick=10,color=0
        contour,marksfl,x,alat,/overplot,levels=[-0.1],thick=10,color=mcolor
        MAP_SET,0,0,0,/noeras,/contin,/grid,title='PV Gradient',$
              color=0,limit=[40.,0.,90.,360.]
xyouts,-180.,40.,'40',/data,color=0,alignment=1
xyouts,-180.,50.,'50',/data,color=0,alignment=1
xyouts,-180.,60.,'60',/data,color=0,alignment=1
xyouts,-180.,70.,'70',/data,color=0,alignment=1
xyouts,-180.,80.,'80',/data,color=0,alignment=1
xyouts,-180.,90.,'90',/data,color=0,alignment=1
xyouts,xmn-0.06,ymn+0.1,'Latitude',/normal,color=0,orientation=90
xyouts,0.,35.,'0',/data,color=0,alignment=0.5
xyouts,90.,35.,'90',/data,color=0,alignment=0.5
xyouts,270.,35.,'-90',/data,color=0,alignment=0.5
xyouts,xmn+0.275,ymn-0.05,'Longitude',/normal,color=0

; vertical color bar
      xmnb=xorig(0)+xlen+0.1
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(sflevel)
      omax=max(sflevel)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor

      !type=2^2+2^3
      xmn=xorig(1)
      xmx=xorig(1)+xlen
      ymn=yorig(1)
      ymx=yorig(1)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      MAP_SET,0,0,0,/noeras,/contin,title='13 February 2002 '+stheta+' K   PV',$
              color=0,limit=[40.,0.,90.,360.]
      index=where(lat gt 0.)
      sfmin=min(pv(index))
      sfmax=max(pv(index))
      nlvls=20
      sfint=(sfmax-sfmin)/nlvls
      sflevel=sfmin+sfint*findgen(nlvls)
      col1=1+indgen(nlvls)*mcolor/nlvls
      contour,pv,x,alat,/overplot,levels=sflevel,c_color=col1,/cell_fill,/noeras
      contour,pv,x,alat,/overplot,levels=sflevel,color=0,c_linestyle=sflevel lt 0
      contour,pv,x,alat,/overplot,levels=[0],color=0,thick=5
      contour,marksfl,x,alat,/overplot,levels=[0.1],thick=10,color=0
      contour,marksfl,x,alat,/overplot,levels=[-0.1],thick=10,color=mcolor
      MAP_SET,0,0,0,/noeras,/contin,/grid,color=0,limit=[40.,0.,90.,360.]
xyouts,-180.,40.,'40',/data,color=0,alignment=1
xyouts,-180.,50.,'50',/data,color=0,alignment=1
xyouts,-180.,60.,'60',/data,color=0,alignment=1
xyouts,-180.,70.,'70',/data,color=0,alignment=1
xyouts,-180.,80.,'80',/data,color=0,alignment=1
xyouts,-180.,90.,'90',/data,color=0,alignment=1
xyouts,xmn-0.06,ymn+0.1,'Latitude',/normal,color=0,orientation=90
xyouts,0.,35.,'0',/data,color=0,alignment=0.5
xyouts,90.,35.,'90',/data,color=0,alignment=0.5
xyouts,270.,35.,'-90',/data,color=0,alignment=0.5
xyouts,xmn+0.275,ymn-0.05,'Longitude',/normal,color=0

; vertical color bar
      xmnb=xorig(0)+xlen+0.1
      xmxb=xmnb+cbarydel
      set_viewport,xmnb,xmxb,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(sflevel)
      omax=max(sflevel)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],color=0
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor
      !type=2^2+2^3+2^5

      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim merc_'+ifile+'_'+stheta+'K_pvgrad+mark.ps -rotate -90 merc_'+ifile+'_'+stheta+'K_pvgrad+mark.jpg'
      endif

    jump:
endfor		; loop over files
end
