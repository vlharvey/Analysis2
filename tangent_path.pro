pro tangent_path,n,nclus,nth,thsfc,x,y,xsat,ysat,xn,yn,zn

; This subroutine passes back the 3d location of the tangent path
; given the occultation lon,lat and the satellite lon,lat.
; HALOE tangent path is 300 km

PI2 = 6.2831853071796
DTR = ASIN(1.0)/90.
RAD = 6.37E6
s   = 300.e3
ds  = s/nclus

for k=0,nth-1 do begin

; co-latitudes
clathal=dtr*(90.-y(n,k))
clatsat=dtr*(90.-ysat(n,k))
dlamba=(xsat(n,k)-x(n,k))
if dlamba gt  180. then dlamba=dlamba-360.
if dlamba lt -180. then dlamba=dlamba+360.
dlamba=dtr*dlamba

; law of cosines for sphere to get angle between sat and hal
term=cos(clatsat)*cos(clathal)+$
     sin(clatsat)*sin(clathal)*cos(dlamba)
angle1=acos(term)
 
; great circle distance (meters) between sat and hal
zterm=angle1*rad
 
; longitudinal projection (sx) of the tangent path (s)
xterm=dlamba*rad*cos(y(n,k)*dtr)
sx=ds*xterm/zterm
dx=sx*360./(pi2*rad*cos(y(n,k)*dtr))

; latitude is determined by equation for great circle in spherical
; coordinates.  first, x,y,z coordinates of hal and sat
xm=sin(clathal)*cos(dtr*x(n,k))
ym=sin(clathal)*sin(dtr*x(n,k))
zm=cos(clathal)
xs=sin(clatsat)*cos(dtr*xsat(n,k))
ys=sin(clatsat)*sin(dtr*xsat(n,k))
zs=cos(clatsat)

; second, coeff for z=f(x,y) plane defined by sat and hal
; and center of the earth in 3d
term=(xs*ym-ys*xm)
if term eq 0. then term=.0000001
cm=(zs*ym-ys*zm)/term
cn=(xs*zm-zs*xm)/term

xn(0,k)=x(n,k)
yn(0,k)=y(n,k)
zn(0,k)=thsfc(k)
xmid=x(n,k)
ymid=y(n,k)
zmid=thsfc(k)

for i=0,nclus/2-2 do begin
    xn(i+1,k)=xn(i,k)+dx
    zn(i+1,k)=thsfc(k)
    term1=(1.+cm^2.0)*cos(xn(i+1,k)*dtr)^2.0
    term2=(1.+cn^2.0)*sin(xn(i+1,k)*dtr)^2.0
    term3=2.*cn*cm*cos(xn(i+1,k)*dtr)*sin(xn(i+1,k)*dtr)
    term=sqrt(1./(term1+term2+term3))
    yn(i+1,k)=90. - asin(term)/dtr
    if y(n,k) lt 0. then yn(i+1,k)=-yn(i+1,k)
    IF xn(i+1,k) LT 0. then xn(i+1,k)=360.+xn(i+1,k)
    IF yn(i+1,k) GT 90.then begin
        yn(i+1,k)=180.-yn(i+1,k)
        xn(i+1,k)=180.+xn(i+1,k)
    ENDIF
    IF yn(i+1,k) LT -90. THEN begin
       yn(i+1,k)=-180.-yn(i+1,k)
       xn(i+1,k)= 180.+xn(i+1,k)
    ENDIF
    IF xn(i+1,k) GT 360. THEN begin
       MULT=xn(i+1,k)/360.
       xn(i+1,k)=xn(i+1,k)-MULT*360.
    ENDIF
endfor
for i=nclus/2-1,nclus-2 do begin
    xn(i+1,k)=xn(i,k)-dx
    zn(i+1,k)=thsfc(k)
    if i eq nclus/2-1 then xn(i+1,k)=xmid-dx
    if i eq nclus/2-1 then zn(i+1,k)=zmid
    term1=(1.+cm^2.0)*cos(xn(i+1,k)*dtr)^2.0
    term2=(1.+cn^2.0)*sin(xn(i+1,k)*dtr)^2.0
    term3=2.*cn*cm*cos(xn(i+1,k)*dtr)*sin(xn(i+1,k)*dtr)
    term=sqrt(1./(term1+term2+term3))
    yn(i+1,k)=90. - asin(term)/dtr
    if y(n,k) lt 0. then yn(i+1,k)=-yn(i+1,k)
    IF xn(i+1,k) LT 0. then xn(i+1,k)=360.+xn(i+1,k)
    IF yn(i+1,k) GT 90.then begin
       yn(i+1,k)=180.-yn(i+1,k)
       xn(i+1,k)=180.+xn(i+1,k)
    ENDIF
    IF yn(i+1,k) LT -90. THEN begin
       yn(i+1,k)=-180.-yn(i+1,k)
       xn(i+1,k)= 180.+xn(i+1,k)
    ENDIF
    IF xn(i+1,k) GT 360. THEN begin
       MULT=xn(i+1,k)/360.
       xn(i+1,k)=xn(i+1,k)-MULT*360.
    ENDIF
endfor
endfor
end
