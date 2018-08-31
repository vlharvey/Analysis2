;*********************************************
pro calcelengthhem, q, lon, lat, mlm, xi, elat
;*********************************************

;*************************************************************************
; This code calculates the equivalent length over one hemisphere. This is
; a modified version of calcelength.pro, which is a global calculation.
; Care must be taken in only using one hemisphere, as there may be edge
; effects at the equator.
;
; Last modified by Doug Allen: 7 July 2005.
; Modified: Cynthia Singleton 7/14/05
;*************************************************************************

nelat=n_elements(elat)       ; elat grid
nlon=n_elements(lon)         ; regularly-gridded latitude (degrees)
nlat=n_elements(lat)         ; regularly-gridded longitude (degrees)
dlat=lat(1)-lat(0)           ; latitude grid spacing
dlon=lon(1)-lon(0)
elat2d=FltArr(nlon,nlat)     ; equivalent latitude on 2D grid
a=6.37e3                     ; radius of earth in km
hemarea=2.*!pi*a^2           ; surface area of hemisphere in km^2 (Surface area of Globe = 4.!pi*r^2
latarea=FltArr(nlat)         ; area of grid point at each latitude

dy=2*!pi*a*dlat/360.                 ; great circle distance between two latitudes
dx=2*!pi*a*dlon/360.*cos(!dtor*lat)  ; great circle distance between two longitudes
dbdx=FltArr(nlon,nlat)               ; zonal derivative
dbdy=FltArr(nlon,nlat)               ; meridional derivative
grad2=FltArr(nlon,nlat)              ; gradient squared
xi=FltArr(nelat)                     ; normalized equivalent length
mlm=FltArr(nelat)                    ; mlm
q2int=FltArr(nelat)                  ; integrated gradient squared
qtotarea=FltArr(nelat)               ; integrated area

;------ Calculate area for each grid point ------

if(min(lat) eq -90.) then begin
  ;------ End grid points on the poles ------

  phi1=!pi/180.*(-90.)
  phi2=!pi/180.*(-90.+dlat/2.)
  latarea(0) = hemarea*abs(sin(phi1)-sin(phi2))/float(nlon)

  for i=1,nlat-2 do begin
    phi1=!pi/180.*(lat(i)-dlat/2.)
    phi2=!pi/180.*(lat(i)+dlat/2.)
    latarea(i) = hemarea*abs(sin(phi1)-sin(phi2))/float(nlon)

  endfor

  phi1=!pi/180.*(90.-dlat/2.)
  phi2=!pi/180.*(90.)

  latarea(nlat-1) = hemarea*abs(sin(phi1)-sin(phi2))/float(nlon)

endif else begin

  ;------ End grid points 1/2 grid off the poles ------

  for i=0,nlat-1 do begin
    phi1=!pi/180.*(lat(i)-dlat/2.)
    phi2=!pi/180.*(lat(i)+dlat/2.)
    latarea(i) = hemarea*abs(sin(phi1)-sin(phi2))/float(nlon)
stop
  endfor
endelse

totalarea=total(latarea)*nlon ; Note: this is for a hemisphere!!!

;------ calculate dbdx ------------------------------
; here we extend the data 2 grid points around the GM
;----------------------------------------------------

for j=0,nlat-1 do begin
  dd=FltArr(nlon+4)
  dd(0)=q(nlon-2,j)
  dd(1)=q(nlon-1,j)
  dd(2:nlon+1)=q(*,j)
  dd(nlon+2)=q(1,j)
  dd(nlon+3)=q(2,j)
  dderiv=deriv(dd)/dx(j)

  dbdx(*,j) = dderiv(2:nlon+1)
endfor

;------ calculate dbdy ---------------------------------
; here we extend the data 2 grid points beyond each pole
;-------------------------------------------------------

for i=0,nlon-1 do begin
  ;------ Get longitude 180 degrees from given longitude ------
  merid=nlon/2
  if(i lt merid) then ii=i+merid
  if(i ge merid) then ii=i-merid

  dd=FltArr(nlat+2)

  dd(0)=q(ii,1)
  dd(1)=q(ii,0)
  dd(2:nlat+1)=q(i,*)

  dderiv=deriv(dd)/dy
  dbdy(i,*) = dderiv(2:nlat+1)
endfor

grad2=dbdx^2+dbdy^2

;------ order the data points by increasing value ------

npoints=long(nlon)*long(nlat)

q1d=FltArr(npoints)                ; Data placed in 1-D array
a1d=FltArr(npoints)                ; Area placed in 1-D array
el1d=FltArr(npoints)               ; Elat placed in 1-D array
grad21d=FltArr(npoints)

index=long(0)
for ilon=0,nlon-1 do begin
  for ilat=0,nlat-1 do begin
    q1d(index)=q(ilon,ilat)
    a1d(index)=latarea(ilat)
    grad21d(index)=grad2(ilon,ilat)

    index=index+1
  endfor
endfor
qsort=q1d(sort(q1d))               ; Sorted data
asort=a1d(sort(q1d))               ; Sorted area
grad2sort=grad21d(sort(q1d))       ; Sorted gradient

tasort=FltArr(npoints)             ; Total area for sorted mixing ratio
elsort=FltArr(npoints)             ; Equivalent latitude for sorted mr
grad2intsort=Fltarr(npoints)       ; Integrated gradient

ta=0.                              ; Total area
g2a=0.

for index=long(0),npoints-1 do begin
  ta=ta+asort(index)
  g2a=g2a+asort(index)*grad2sort(index)

  if(ta gt hemarea) then ta=hemarea ; keep area le hemarea
  tasort(index)=ta
  elsort(index)=asin(1.-ta/hemarea)/!dtor
  grad2intsort(index)=g2a
endfor

;------ Interpolate to constant equivalent latitude grid ------

mlm1d=interpol(qsort,elsort,elat)
area1d=interpol(tasort,elsort,elat)
grad2int1d=interpol(grad2intsort,elsort,elat)
dqda=deriv(area1d,mlm1d)
dq2da=deriv(area1d,grad2int1d)

Le2=dqdA^(-2)*dq2dA
Le2=abs(Le2)

;------ Normalize the equivalent length ------

Lo2=(2.*!pi*a*cos(!dtor*elat))^2
xi1d=alog(Le2/Lo2)

xi=xi1d
mlm=mlm1d

end

;*****************************************************************
;                         Main Program
;*****************************************************************
;
; This program calculates the equivalent length on a grid from
; one hemisphere of data. This was written to accomodate the
; output used by Cynthia Singleton. The hemispheric calculation
; is necessary to account for tracers with non-monotonic gradients.
;
; Last modified: 7 July 2005
;
;*****************************************************************

loadct,38
device,decompose=0
icolmax=byte(!p.color)
programstr='~/Programs/Elength/write_elength_grid_hem.pro'

;****
menu:
;****
nchoice=0
print, '1. Read data and calculate equivalent length'
print, '2. Plot data'
print, '3. Quit'
read, nchoice

case(nchoice) of
  1: goto, readdata
  2: goto, plotdata
  3: goto, finish
endcase

;********
readdata:
;********

;------ Equivalent latitude grid for output ------

nelat=45.
delat=90./nelat
elat=90.-delat/2.-Findgen(nelat)*delat

;------ Restore data ------

date=long(20021202)
datestr=string(date,'(I8)')
path='C:\Documents and Settings\shaw.LASP_ENGR\EffectiveDiffusivity\'
restore,path+'SC_NH_CH4'+datestr+'.DAT'
; LAT             FLOAT     = Array[32]
; LON             FLOAT     = Array[128]
; SC_NH           FLOAT     = Array[128, 32, 24]
; THETA           FLOAT     = Array[24]
; YYYYMMDD        LONG      =     20021202

nlat=n_elements(lat)
nlon=n_elements(lon)
ntheta=n_elements(theta)

dlat=90./nlat
latgrid=90.-dlat/2.-Findgen(nlat)*dlat
longrid=lon

data3d=FltArr(nlon,nlat,ntheta) ; data on regular latitude grid

xi2d=FltArr(nelat,ntheta)

mlm2d=FltArr(nelat,ntheta)

;------ Interpolate to regular latitude grid ------

for itheta=0,ntheta-1 do begin
  for ilon=0,nlon-1 do begin
    d=FltArr(nlat)
    d(0:nlat-1)=sc_nh(ilon,*,itheta)
    dd=interpol(d,lat,latgrid)

    data3d(ilon,*,itheta)=dd
  endfor
endfor

for itheta=0,ntheta-1 do begin
  data2d=data3d(*,*,itheta)
  calcelengthhem, data2d, longrid, latgrid, mlm, xi, elat

  xi2d(*,itheta)=xi
  mlm2d(*,itheta)=mlm
endfor

;------ Calculate zonal mean for later use ------

datazm=FltArr(nlat,ntheta)
for ilat=0,nlat-1 do begin
  for itheta=0,ntheta-1 do begin
    d=data3d(*,ilat,itheta)
    dd=mean(d)
    datazm(ilat,itheta)=dd
  endfor
endfor

goto, menu

;********
plotdata:
;********

!P.multi=[0,2,2,0,0]
!P.font=0
set_plot, 'ps'
device, filename=path+'test2.ps', /landscape, /color

nlev=20
dmin=min(datazm*1.e6,/Nan)
dmax=max(datazm*1.e6,/Nan)
dlev=float(dmax-dmin)/(nlev-1)
levels=dmin+Findgen(nlev)*dlev

contour, datazm*1.e6, latgrid, theta, xrange=[0,90], xstyle=1, $
  xtitle='Latitude', xticks=3, xminor=3, yrange=[min(theta),2000.], $
  ystyle=1, ytitle='Potential Temperature (K)', levels= levels, /foll, $
  title='Zonal Mean CH4 for '+datestr

contour, mlm2d*1.e6, elat, theta, xrange=[0,90], xstyle=1, $
  xtitle='Equiv. Latitude', xticks=3, xminor=3, yrange=[min(theta),2000.], $
  ystyle=1, ytitle='Potential Temperature (K)', levels= levels, /foll, $
  title='Modified Lagrangian Mean CH4 for '+datestr

xilevels=Findgen(18)*0.2
nlvls=n_elements(xilevels)
col1=1+indgen(nlvls)*icolmax/nlvls
contour, xi2d, elat, theta, xrange=[0,90], xstyle=1, $
  xtitle='Equiv. Latitude', xticks=3, xminor=3, yrange=[min(theta),2000.], $
  ystyle=1, ytitle='Potential Temperature (K)', levels= xilevels,  $
  title='Normalized Eq. Length for '+datestr, /cell, c_color=col1

contour, xi2d, elat, theta, levels= xilevels, /foll, /overplot

oplot, [0,90], [theta(10),theta(10)]

plot, elat, xi2d(*,10), xrange=[0,90], xstyle=1, $
  xtitle='Equiv. Latitude', xticks=3, xminor=3, $
  yrange=[0,3], ystyle=1, ytitle='Norm. Eq. Length', $
  title='Normalized Eq. Length at '+string(theta(10),'(I3)')+' K'

device, /close
set_plot, 'win'

goto, menu

;******
finish:
;******

end
