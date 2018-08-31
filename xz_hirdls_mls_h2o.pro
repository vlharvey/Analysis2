@aura2date
@rd_ukmo_nc3

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
nlvls=20
col1=1+indgen(nlvls)*mcolor/nlvls
device,decompose=0
icmm1=icolmax-1
icmm2=icolmax-2
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.15,0.15]
yorig=[0.15,0.60]
xlen=0.7
ylen=0.25
cbaryoff=0.08
cbarydel=0.02
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='xz_hirdls_mls_h2o.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
SpeciesNames = ['Temperature',$
                'H2O', $
                'O3',  $
                'N2O', $
                'HNO3']
GeoLoc = ['Pressure',$
          'Time',$
          'Latitude',$
          'Longitude',$
          'SolarZenithAngle',$
          'LocalSolarTime']
hdir='/aura3/data/HIRDLS_data/Datfiles/'
Hfile=hdir+'HIRDLS2_2000d276_MZ3_c1.he5'

;; load HIRDLS data all at once
hirdls=LoadAuraData(Hfile, [GeoLoc, SpeciesNames])

;; file header and tail for MLS
Mfileh=hdir+'MLS-Aura_L2GP-'
Mfilet='_sAura2c--t_2000d276.he5'

;; loop over species
FOR is = 0,N_ELEMENTS(SpeciesNames)-1 DO BEGIN
    SpeciesName = SpeciesNames(is)

;; construct MLS filename
    Mfile=Mfileh + SpeciesName + Mfilet

;; load geolocation data for MLS from 1st species file
    IF is EQ 0 THEN mls = LoadAuraData(Mfile, GeoLoc)

;; extend the mls structure with the species data
    mls=LoadAuraData(Mfile, SpeciesName, mls)
ENDFOR
;
; extract mls and hirdls variables
; time is elapsed seconds since midnight 1 Jan 1993
;
mpress=mls.p		; P               FLOAT     Array[37]
mlev=n_elements(mpress)
mtime=mls.time		; TIME            DOUBLE    Array[3495]
mlat=mls.lat		; LAT             FLOAT     Array[3495]
mlon=mls.lon		; LON             FLOAT     Array[3495]
msza=mls.sza		; SZA             FLOAT     Array[3495]
mlst=mls.lst		; LST             FLOAT     Array[3495]
mprof=n_elements(mlst)
mtemp=mls.t		; T               FLOAT     Array[37, 3495]
mh2o=mls.h2o		; H2O             FLOAT     Array[37, 3495]
mo3=mls.o3		; O3              FLOAT     Array[37, 3495]
mn2o=mls.n2o		; N2O             FLOAT     Array[37, 3495]
mhno3=mls.hno3		; HNO3            FLOAT     Array[37, 3495]

hpress=hirdls.p		;   P               FLOAT     Array[145]
hlev=n_elements(hpress)
htime=hirdls.time	;   TIME            DOUBLE    Array[7848]
hlat=hirdls.lat		;   LAT             FLOAT     Array[7848]
hlon=hirdls.lon		;   LON             FLOAT     Array[7848]
hsza=hirdls.sza		;   SZA             FLOAT     Array[7848]
hlst=hirdls.lst		;   LST             FLOAT     Array[7848]
hprof=n_elements(hlst)
htemp=hirdls.t		;   T               FLOAT     Array[145, 7848]
hh2o=hirdls.h2o		;   H2O             FLOAT     Array[145, 7848]
ho3=hirdls.o3		;   O3              FLOAT     Array[145, 7848]
hn2o=hirdls.n2o		;   N2O             FLOAT     Array[145, 7848]
hno3=hirdls.hno3	;   HNO3            FLOAT     Array[145, 7848]
;
; convert elapsed seconds to dates (yyyymmddhh)
;
aura2date,mdate,mtime
aura2date,hdate,htime
;
; make press,lat,lon 2d
;
mpress2=0.*mo3
mlat2=0.*mo3
mlon2=0.*mo3
for i=0L,mprof-1L do mpress2(*,i)=mpress
for i=0L,mlev-1L do begin
    mlat2(i,*)=mlat
    mlon2(i,*)=mlon
endfor
hpress2=0.*ho3
hlat2=0.*ho3
hlon2=0.*ho3
for i=0L,hprof-1L do hpress2(*,i)=hpress
for i=0L,hlev-1L do begin
    hlat2(i,*)=hlat
    hlon2(i,*)=hlon
endfor
sdate=strcompress(string(mdate(mprof-1)),/remove_all)
sdate=strmid(sdate,0,8)
;
; read UKMO data
;
syr=strmid(sdate,0,4)
uyr=strmid(syr,2,2)
smn=strmid(sdate,4,2)
imn=fix(smn)
sdy=strmid(sdate,6,2)
ifile=mon(imn-1)+sdy+'_'+uyr
print,ifile
rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
            pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
x2d=fltarr(nc,nth)
y2d=fltarr(nc,nth)
for k=0,nth-1 do x2d(*,k)=alon
for j=0,nc-1 do y2d(j,*)=th
rlat=1.25000
slat=strcompress(string(rlat),/remove_all)
index=where(alat eq rlat)
ilat=index(0)
v=reform(v2(ilat,*,*))
u=reform(u2(ilat,*,*))
p=reform(p2(ilat,*,*))
pmn=fltarr(nth)
for k=0,nth-1 do pmn(k)=total(p(*,k))/float(nc)
;
; plot
;
!type=2^2+2^3
erase
xyouts,.1,.93,sdate+' '+' Equatorial Water Vapor and Zonal Wind',charsize=1.5,/normal
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
omin=2.
omax=7.
level=-50.+5.*findgen(21)
contour,u,alon,pmn,levels=level,xrange=[0.,360.],yrange=[1000.,.1],/ylog,$
        /noeras,title='MLS',charsize=2,/follow,/nodata
index=where(abs(mlat2) le 2.,npt)
mx=mlon2(index) & my=mpress2(index) & mwater=mh2o(index)*1.e6
index=where(mx lt 0.)
if index(0) ne -1 then mx(index)=mx(index)+360.
print,'mls ',min(mwater),max(mwater)
for i=0L,npt-1L do $
    oplot,[mx(i),mx(i)],[my(i),my(i)],psym=8,$
          color=((mwater(i)-omin)/(omax-omin))*icolmax
contour,u,alon,pmn,levels=level,/overplot,/follow,color=lc,c_linestyle=level lt 0,thick=2
ymnb=ymn-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[omin,omax],[0,0],yrange=[0,10],xrange=[omin,omax],xtitle='(ppmv)'
ybox=[0,10,10,0,0]
x1=omin
dx=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
!type=2^2+2^3
set_viewport,xmn,xmx,ymn,ymx
contour,u,alon,pmn,levels=level,xrange=[0.,360.],yrange=[1000.,.1],/ylog,$
        /noeras,title='HIRDLS',charsize=2,/follow,/nodata
index=where(abs(hlat2) le 2.,npt)
hx=hlon2(index) & hy=hpress2(index) & hwater=hh2o(index)*1.e6
index=where(hx lt 0.)
if index(0) ne -1 then hx(index)=hx(index)+360.
print,'hirdls ',min(hwater),max(hwater)
for i=0L,npt-1L do $
    oplot,[hx(i),hx(i)],[hy(i),hy(i)],psym=8,$
          color=((hwater(i)-omin)/(omax-omin))*icolmax
contour,u,alon,pmn,levels=level,/overplot,/follow,color=lc,c_linestyle=level lt 0,thick=2
ymnb=ymn-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[omin,omax],[0,0],yrange=[0,10],xrange=[omin,omax],xtitle='(ppmv)'
ybox=[0,10,10,0,0]
x1=omin
dx=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert xz_hirdls_mls_h2o.ps -rotate -90 xz_hirdls_mls_h2o.jpg'
endif

END
