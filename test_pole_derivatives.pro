
; Plot markers.  cyclones = 1;  anticyclones = -1

@rd_ukmo_nc3
@compvort_v2
@compvort

nlvls=20
loadct,38
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.05,0.55]
yorig=[0.25,0.25]
xlen=0.4
ylen=0.4
cbaryoff=0.02
cbarydel=0.01
device,decompose=0
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='pole_deriv.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
   xsize=xsize,ysize=ysize
endif

diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'ukmo_files.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin

; Read UKMO isentropic data
    iflag=0
    readf,1,ifile
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump
    x=fltarr(nc+1)
    x(0:nc-1)=alon(0:nc-1)
    x(nc)=alon(0)+360.

; loop over theta from top down
    FOR thlev=0,nth-1 DO BEGIN
        theta=th(thlev)
        print,'theta level=',theta
        u1=transpose(u2(*,*,thlev))
        v1=transpose(v2(*,*,thlev))
        
; introduce relative vorticity
        zeta1old=u1*0.0
        zeta1new=u1*0.0
        compvort,u1,v1,zeta1old,alon,alat,nc,nr
        compvort_v2,u1,v1,zeta1new,alon,alat,nc,nr
        zetaold=fltarr(nc+1,nr)
        zetaold(0:nc-1,0:nr-1)=zeta1old(0:nc-1,0:nr-1)*1.e5
        zetaold(nc,*)=zetaold(0,*)
        zetanew=fltarr(nc+1,nr)
        zetanew(0:nc-1,0:nr-1)=zeta1new(0:nc-1,0:nr-1)*1.e5
        zetanew(nc,*)=zetanew(0,*)
; PLOT
        erase
        !type=2^2+2^3
        xyouts,.1,.8,string(fix(theta))+'K Relative Vorticity on '+ifile,/normal,charsize=2
        xmn=xorig(0)
        xmx=xorig(0)+xlen
        ymn=yorig(0)
        ymx=yorig(0)+ylen
        set_viewport,xmn,xmx,ymn,ymx
        MAP_SET,90,0,0,/ortho,/CONTIN,/noeras,limit=[85.,0.,90.,360.],title='OLD',charsize=2
        level=-25.+2.5*findgen(19)
        contour,zetaold,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
        contour,zetaold,x,alat,levels=level,/overplot,/follow,c_color=0,c_linestyle=level lt 0
        MAP_SET,90,0,0,/ortho,/CONTIN,/noeras,limit=[85.,0.,90.,360.],color=0
        oplot,findgen(361),alat(nr-1)+0.*findgen(361),psym=3
        oplot,findgen(361),alat(nr-2)+0.*findgen(361),psym=3
        velovect,u1,v1,alon,alat,/overplot
        imin=min(level)
        imax=max(level)
        ymnb=ymn -cbaryoff
        ymxb=ymnb+cbarydel
        set_viewport,xmn,xmx,ymnb,ymxb
        !type=2^2+2^3+2^6
        plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
        ybox=[0,10,10,0,0]
        x2=imin
        dx=(imax-imin)/float(icmm1)
        for j=1,icmm1 do begin
            xbox=[x2,x2,x2+dx,x2+dx,x2]
            polyfill,xbox,ybox,color=j
            x2=x2+dx
        endfor

        xmn=xorig(1)
        xmx=xorig(1)+xlen
        ymn=yorig(1)
        ymx=yorig(1)+ylen
        set_viewport,xmn,xmx,ymn,ymx
        MAP_SET,90,0,0,/ortho,/CONTIN,/noeras,limit=[85.,0.,90.,360.],title='NEW',charsize=2
        contour,zetanew,x,alat,levels=level,/overplot,/cell_fill,c_color=col1
        contour,zetanew,x,alat,levels=level,/overplot,/follow,c_color=0,c_linestyle=level lt 0
        MAP_SET,90,0,0,/ortho,/CONTIN,/noeras,limit=[85.,0.,90.,360.],color=0
        oplot,findgen(361),alat(nr-1)+0.*findgen(361),psym=3
        oplot,findgen(361),alat(nr-2)+0.*findgen(361),psym=3
        velovect,u1,v1,alon,alat,/overplot
        imin=min(level)
        imax=max(level)
        ymnb=ymn -cbaryoff
        ymxb=ymnb+cbarydel
        set_viewport,xmn,xmx,ymnb,ymxb
        !type=2^2+2^3+2^6
        plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras
        ybox=[0,10,10,0,0]
        x2=imin
        dx=(imax-imin)/float(icmm1)
        for j=1,icmm1 do begin
            xbox=[x2,x2,x2+dx,x2+dx,x2]
            polyfill,xbox,ybox,color=j
            x2=x2+dx
        endfor
       stop
    ENDFOR  ; loop over theta
jump:
endfor; loop over days
end
