;
; check smoothed edge
;
@rd_ukmo_nc3

ipan=0
npp=1
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.1]
yorig=[0.15]
xlen=0.7
ylen=0.7
cbaryoff=0.03
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'ukmo_sh_files_00.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
print,ifile
    iflag=0
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump
    index=where(marksf2 ne 0.)
    if index(0) ne -1 then $
       marksf2(index)=marksf2(index)/abs(marksf2(index))
    savemark=marksf2
    marksmooth2=0.*marksf2
;
; smooth marker field
;
    for i=0,nc-1 do begin
        ip1=i+1
        im1=i-1
        if i eq 0 then im1=nc-1
        if i eq nc-1 then ip1=0
        for j=0,nr-1 do begin
            jp1=j+1
            jm1=j-1
            if j eq 0 then jm1=j
            if j eq nr-1 then jp1=j
            for k=0,nth-1 do begin
                kp1=k+1
                km1=k-1
                if k eq 0 then km1=k
                if k eq nth-1 then kp1=k
                marksmooth2(j,i,k)=( (1./4.)*(savemark(j,im1,k) + $
                               2.0*savemark(j,i,k)+savemark(j,ip1,k)) + $
                                  (1./4.)*(savemark(jm1,i,k) + $
                               2.0*savemark(j,i,k)+savemark(jp1,i,k)) + $
                                  (1./4.)*(savemark(j,i,km1) + $
                               2.0*savemark(j,i,k)+savemark(j,i,kp1)))/3.0
             endfor
        endfor
    endfor


if n eq 0 then begin
;theta=0.
;print,th
;read,'Enter theta ',theta
;index=where(theta eq th)
;if index(0) eq -1 then stop,'Invalid theta level '
;thlev=index(0)
endif
for thlev=0,nth-1 do begin
    theta=th(thlev)

; Height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
    speed2=sqrt(u2^2+v2^2)
    temp2=0.*pv2
    zth2=0.*pv2
    for k=0,nth-1 do begin
        temp2(0:nr-1,0:nc-1,k)=th(k)*( (p2(0:nr-1,0:nc-1,k)/1000.)^(.286) )
        zth2(0:nr-1,0:nc-1,k)=(msf2(0:nr-1,0:nc-1,k)-1004.* $
                              temp2(0:nr-1,0:nc-1,k))/(9.86*1000.)
    endfor

    stheta=strcompress(string(fix(theta)),/remove_all)
    qdf1=transpose(qdf2(*,*,thlev))
    sf1=transpose(sf2(*,*,thlev))
    pv1=transpose(pv2(*,*,thlev))
    temp1=transpose(temp2(*,*,thlev))
    zth1=transpose(zth2(*,*,thlev))
    marksfl1=transpose(marksf2(*,*,thlev))
    marksfh1=transpose(marksf2(*,*,thlev))
    marksmooth1=transpose(marksmooth2(*,*,thlev))
print,stheta+' K'
index=where(marksfl1 eq 1.)
if index(0) ne -1 then print,'Min in Vortex ',min(marksmooth1(index))
index=where(marksfh1 eq -1.) 
if index(0) ne -1 then print,'Max in Highs  ',max(marksmooth1(index))
index=where(marksfl1 eq 0.)
print,'Outside Range ',min(marksmooth1(index)),max(marksmooth1(index))

    qdf=0.*fltarr(nc+1,nr)
    qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
    qdf(nc,*)=qdf(0,*)
    sf=0.*fltarr(nc+1,nr)
    sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
    sf(nc,*)=sf(0,*)
    pv=0.*fltarr(nc+1,nr)
    pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
    pv(nc,*)=pv(0,*)
    marksfl=0.*fltarr(nc+1,nr)
    marksfl(0:nc-1,0:nr-1)=marksfl1(0:nc-1,0:nr-1)
    marksfl(nc,*)=marksfl(0,*)
    marksfh=0.*fltarr(nc+1,nr)
    marksfh(0:nc-1,0:nr-1)=marksfh1(0:nc-1,0:nr-1)
    marksfh(nc,*)=marksfh(0,*)
    marksmooth=0.*fltarr(nc+1,nr)
    marksmooth(0:nc-1,0:nr-1)=marksmooth1(0:nc-1,0:nr-1)
    marksmooth(nc,*)=marksmooth(0,*)
    temp=0.*fltarr(nc+1,nr)
    temp(0:nc-1,0:nr-1)=temp1(0:nc-1,0:nr-1)
    temp(nc,*)=temp(0,*)
    zth=0.*fltarr(nc+1,nr)
    zth(0:nc-1,0:nr-1)=zth1(0:nc-1,0:nr-1)
    zth(nc,*)=zth(0,*)
    x=fltarr(nc+1)
    x(0:nc-1)=alon
    x(nc)=alon(0)+360.
    lon=0.*sf
    lat=0.*sf
    for i=0,nc   do lat(i,*)=alat
    for j=0,nr-1 do lon(*,j)=x

    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename=ifile+'_'+stheta+'K_smoothed_edge.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
               xsize=xsize,ysize=ysize
    endif

; Set plot boundaries
    erase
    !psym=0
    MAP_SET,90,0,-90,/stereo,/noeras,/grid,/contin,/noborder,$
            title='!6'+ifile+'  '+stheta+' K',charsize=2.0,latdel=10
    !psym=3
    oplot,findgen(361),0.1+0.*findgen(361)
    nlvls=21
    level=-1.+0.1*findgen(nlvls)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,marksmooth,x,alat,/overplot,levels=level,c_color=col1,/cell_fill,/noeras
    contour,marksmooth,x,alat,/overplot,levels=level,c_color=0,/follow,/noeras
index=where(marksfh ne -1.)
marksfh(index)=-9999.
    contour,marksfh,x,alat,/overplot,levels=[-1.0],c_color=icolmax,thick=2,/follow,/noeras
    contour,marksfl,x,alat,/overplot,levels=[1.0],c_color=icolmax,thick=2,/follow,/noeras
    MAP_SET,90,0,-90,/stereo,/noeras,/grid,/contin,/noborder,charsize=2.0,latdel=10
;   if setplot ne 'ps' then stop
    if setplot eq 'ps' then device, /close
    jump:
endfor		; loop over levels
endfor		; loop over files
end
