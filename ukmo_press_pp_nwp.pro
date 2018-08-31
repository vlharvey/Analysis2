; Read in 2.5x3.75 topograpy data and use geopotential height with it
; to calculate the pressure at the station height.  Mean sea level
; pressure is computed hydrostatically.
;
; output pressure in Pascals

@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_pp_nwp
@store_ukmo_pp

lstmn=3 & lstdy=16 & lstyr=8 & lstday=0
ledmn=3 & leddy=24 & ledyr=8 & ledday=0

nlg =96L
nlat=73L
nlv =27L
g3d=fltarr(nlg,nlat,nlv)
t3d=fltarr(nlg,nlat,nlv)
u3d=fltarr(nlg,nlat,nlv)
v3d=fltarr(nlg,nlat,nlv)
psfc=fltarr(nlg,nlat)
mslp=fltarr(nlg,nlat)
g1000=fltarr(nlg,nlat)
t1000=fltarr(nlg,nlat)
gup=0.
glw=0.
gtp=0.
y=fltarr(nlat)
y=90.-2.5*findgen(nlat)
x=fltarr(nlg+1)
x=0.+3.75*findgen(nlg+1)
;
; read global topography (geopotential meters)
; topg3.75 matches daily UKMO z,t format 90- -90 and 0-360
;
topg=fltarr(nlg,nlat)
close,10
openr,10,'topg3.75'
readu,10,topg
close,10
;
; convert geopotential to geopotential height
;
topg=topg/9.806
diri='/aura7/harvey/UKMO_data/Datfiles/'
diro='/aura3/data/UKMO_data/Datfiles/'
;
; Ask interactive questions- get starting/ending date
;
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays gt ledday then stop,' normal termination condition'
;
; read UKMO data
;
      if iyr ge 2000 then iyr1=iyr-2000
      if iyr lt 2000 then iyr1=iyr-1900
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      uyr=strmid(syr,2,2)
      ifile=string(FORMAT='(a8,i2.2,a2,i2.2,a2,i2.2,a11)',$
            'ppassm_y',iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo_pp_nwp,diri+ifile,iflg,nlg,nlat,nlv,alon,alat,wlon,wlat,p,$
                g3d,t3d,u3d,v3d
      if iflg ne 0 then goto, jump
      pa=p*100.

; first calculate array of mslp (in mb) from Z*g at 1000mb (100000Pa)
; Is there any improvement using virtual temp?  Tv=T(1+.61w) or
; Tv=T/( 1-e/p(1-.622) ) where w is in kg/kg and e is vapor pressure
      g1000(0:nlg-1,0:nlat-1)=g3d(0:nlg-1,0:nlat-1,0)
      t1000(0:nlg-1,0:nlat-1)=t3d(0:nlg-1,0:nlat-1,0)
      mslp=100000.+((100000.*g1000*9.806)/(287.*t1000))

; for each horizontal grid point, interpolate in altitude UKMO p
; using bounding UKMO Z about topg Z.
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
                     psfc(i,j)=100000.-((glw-gtp)/glw)*(100000.-mslp(i,j))
                     goto,jumpout
                  endif
                  if gup ge gtp AND glw le gtp then begin	; UKMO Z's bounding topg
                     psfc(i,j)=pa(k+1)-((gup-gtp)/(gup-glw))*(pa(k+1)-pa(k))
                     goto,jumpout
                  endif
              endfor
              jumpout:
;             if psfc(i,j) eq 0. then print,gtp,g3d(i,j,0),mslp(i,j),i,j
           endfor
       endfor
;
; store daily files
;
      store_ukmo_pp,diro+ifile+'.psfc',nlg,nlat,alon,alat,mslp,psfc

goto,jump
end
