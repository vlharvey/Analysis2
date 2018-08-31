;
; calculate the area enclosed by PV contours
; expressed in equivalent latitude and plotted
; in latitude and time for each level
;
@rd_ukmo_nc3

ipan=0
npp=4
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
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'yt_pv_area.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    print,ifile
    iflag=0
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump

    if n eq 0 then begin
       theta=0.
       print,th
       read,'Enter theta ',theta
       index=where(theta eq th)
       if index(0) eq -1 then stop,'Invalid theta level '
       thlev=index(0)
       stheta=strcompress(string(fix(theta)),/remove_all)
       yt_eqlat=fltarr(nfile,nr)

       pv=transpose(pv2(*,*,thlev))
       lon=0.*pv
       lat=0.*pv
       for i=0,nc-1 do lat(i,*)=alat
       for j=0,nr-1 do lon(*,j)=alon
       area=0.*lat
       yeq=findgen(nr)
       latcircle=fltarr(nr)
       hem_frac=fltarr(nr)
       for j=0,nr-2 do begin
           hy=re*dtr
           dx=re*cos(yeq(j)*dtr)*360.*dtr
           latcircle(j)=dx*hy	; area in each latitude circle
       endfor
       for j=0,nr-1 do begin
           if yeq(j) ge 0. then index=where(yeq ge yeq(j))

; fraction of the hemisphere of each latitude circle
           if index(0) ne -1 then $
              hem_frac(j)=100.*total(latcircle(index))/hem_area
           if yeq(j) eq 0. then hem_frac(j)=100.
       endfor
       deltax=alon(1)-alon(0)
       deltay=alat(1)-alat(0)
       for j=0,nr-1 do begin
           hy=re*deltay*dtr
           dx=re*cos(alat(j)*dtr)*deltax*dtr
           area(j,*)=dx*hy	; area of each grid point
       endfor
    endif
    pv=transpose(pv2(*,*,thlev))

    if ipan eq 0 and setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/portrait,bits=8,$
              filename='yt_pv_area_'+stheta+'.ps'
       device,/color
       device,/inch,xoff=3.,yoff=3.,$
              xsize=xsize,ysize=ysize
    endif


    if setplot eq 'ps' then device, /close
    jump:
endfor		; loop over files
end
