;
; Read in 2.5x2.5 topograpy data and use geopotential height with it
; to calculate the pressure at the station height.  Mean sea level
; pressure is computed hydrostatically.
;
@rd_ecmwf
@store_ecmwf_pp
;
; read global topography (geopotential meters)
; topg2.5 matches daily ECMWF data (-90 to 90 and 0-360)
;
close,10
openr,10,'topg2.5'
nc=0L & nr=0L
readu,10,nc,nr
topg=fltarr(nc,nr)
readu,10,topg
close,10
;
; convert geopotential to geopotential height
;
topg=topg/9.806
dir='/aura3/data/ECMWF_data/Datfiles/'
ifile='                             '
close,3
openr,3,'ecmwf_press_pp.fil'
nfile=0L
readf,3,nfile
for n=0,nfile-1 do begin
    readf,3,ifile
    print,ifile
    rd_ecmwf,dir+ifile,iflg,nlg,nlat,nlv,alon,alat,p,pv,g3d,t3d,u3d,v3d,ww,sh,oz
;
; first calculate array of mslp (in mb) from Z*g at 1000mb (100000Pa)
; Is there any improvement using virtual temp?  Tv=T(1+.61w) or
; Tv=T/( 1-e/p(1-.622) ) where w is in kg/kg and e is vapor pressure
;
    g1000=reform(g3d(*,*,0))
    t1000=reform(t3d(*,*,0))
    mslp=(100000.+((100000.*g1000*9.806)/(287.*t1000)))/100.
    psfc=0.*mslp
;
; for each horizontal grid point, interpolate in altitude ECMWF p
; using bounding ECMWF Z about topg Z.
;
    for i=0,nlg-1 do begin
        for j=0,nlat-1 do begin
            gtp=topg(i,j)
            if gtp le 0. then begin	; if topg below sea level
               psfc(i,j)=mslp(i,j)
               goto,jumpout
            endif
            for k=0,nlv-2 do begin
                gup=g3d(i,j,k+1)
                glw=g3d(i,j,k)
                if glw gt gtp then begin			; 1000 mb Z above topg
                   psfc(i,j)=1000.-((glw-gtp)/glw)*(1000.-mslp(i,j))
                   goto,jumpout
                endif
                if gup ge gtp AND glw le gtp then begin	; ECMWF Z's bounding topg
                   psfc(i,j)=p(k+1)-((gup-gtp)/(gup-glw))*(p(k+1)-p(k))
                   goto,jumpout
                endif
            endfor
            jumpout:
;           if psfc(i,j) eq 0. then print,gtp,g3d(i,j,0),mslp(i,j),i,j
        endfor
    endfor

;erase
;map_set,0,0,0,/contin,/grid
;contour,mslp,alon,alat,nlevels=20,/overplot
;
; store daily files
;
    store_ecmwf_pp,dir+ifile+'.psfc',nlg,nlat,alon,alat,mslp,psfc
endfor
end
