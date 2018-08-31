;pro plt_tem, calc = calc, print = print
;
; adapted from code to read in files like f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4.2011-12-23-00000.nc
; to work on update files that look like sdwaccm2012-2014_1_2_2.cam.h4.2012-01-01-00000.nc
;
; need variables th2d, v2d, w2d, vth2d along with date, time, lev, lat, P0
;
dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD/'
spawn,'ls '+dir+'sdwaccm2012-2014_1_2_2.cam.h4.*-00000.nc',ifiles
nfile=n_elements(ifiles)
for ifile=0L,nfile-1L do begin

;h4file  = 'f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4.2002-2004.nc'
;savfile = 'f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4tem.2002-2004.nc
h4file  = ifiles(ifile)
result=strsplit(ifiles(ifile),'.',/extract)
savfile = result(0)+'.'+result(1)+'.'+result(2)+'tem.'+result(3)+'.sav'
print,savfile

;if keyword_set(print) then begin
  set_plot, 'PS'
  device, file='temh4vlh_'+result(3)+'.eps', xsize=7.0, ysize=7.0, $
     encap = 1, /inches, /color, bits=8
;endif

!p.multi = [0,1,2]
!p.charsize=1.2
!x.style=1
!y.style=1
!y.title='Pressure level (hPa)'
!x.margin=[14,2]
!y.margin=[3,5]
loadct, 0

;if keyword_set(calc) then begin
  calc_tem, h4file, lat, ilev, vstar, wstar, date
  save, lat, ilev, vstar, wstar, date, file = savfile
;endif else begin
;  restore, savfile
;endelse
;
; subset dates
;
;ind0 = where(date eq 20031201)
;print, ind0
;ind1 = where(date eq 20040315)
;print, ind1
;
;date = date(ind0:ind1)
;
; take full date range
;
ind0=0
ind1=n_elements(date)-1L
nt = n_elements(date)

vstar = vstar(*,*,ind0:ind1)
wstar = wstar(*,*,ind0:ind1)

nz = n_elements(ilev)
days = indgen(nt)
xran = [1,nt-2]
;xtickv = [0,31,62,91] 
xtickv = [0,365/4.,365/2.,3.*365/4.,365]
xtickn = date(xtickv)
xticks = 3

wts = cos(!pi*lat/180.)

lat0 = 56.
lat1 = 69.
min0 =  min(abs(lat-lat0), ilat0)
min0 =  min(abs(lat-lat1), ilat1)
ilat0 = ilat0(0)
ilat1 = ilat1(0)

print, ilat0, lat(ilat0), 'to', ilat1, lat(ilat1)
ny = ilat1-ilat0+1

vs = reform(vstar(ilat0:ilat1,*,*),ny,nz*nt)
vs = vs##wts(ilat0:ilat1)
vs = reform(vs,nz,nt)/total(wts(ilat0:ilat1))
vs = transpose(vs)
vs = smooth(vs, [5,1])

vlev = [-100,-2,0,100] 
ccol = [180,230,255,255]

contour, vs, days, ilev, lev=vlev, c_col = ccol, /follow, /ylog, $
  yrange=[100,1e-5], ystyle=1, /cell, title = 'V* 55-70N (m/s)', $
  xticks = xticks, xtickn = xtickn, xtickv = xtickv, xran = xran

contour, vs, days, ilev, lev=[-5,-3,-2,-1,0,1,2,3,4,5], /follow, /over, $
  color=0, c_thick = [2,2,2,2,3,2,2,2,2,2], c_linestyle= [1,1,1,1,0,0,0,0,0,0]

lat0 = 70.
lat1 = 90.
min0 =  min(abs(lat-lat0), ilat0)
min0 =  min(abs(lat-lat1), ilat1)
ilat0 = ilat0(0)
ilat1 = ilat1(0)

print, ilat0, lat(ilat0), 'to', ilat1, lat(ilat1)
ny = ilat1-ilat0+1

ws = reform(wstar(ilat0:ilat1,*,*),ny,nz*nt)
ws = ws##wts(ilat0:ilat1)
ws = reform(ws,nz,nt)/total(wts(ilat0:ilat1))
ws = transpose(ws)
ws = smooth(ws, [5,1])

wlev = (indgen(25)-12)*0.5

wlev = [-100,-1,0,100] 
ccol = [180,230,255,255]

contour, 100*ws, days, ilev, lev=wlev, c_col = ccol, /follow, /ylog, $
  yrange=[100,1e-5], ystyle=1, /cell, title = 'W* 75-90N (cm/s)', $
  xticks = xticks, xtickn = xtickn, xtickv = xtickv, xran = xran

contour, 100*ws, days, ilev, lev=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3], /follow, /over, $
  color=0, c_thick = [2,2,2,2,3,2,2,2,2,2], c_linestyle= [1,1,1,1,0,0,0,0,0,0]

xyouts, /norm, .01, .96, 'a)'
xyouts, /norm, .01, .46, 'b)'

!p.multi = 0

;if keyword_set(print) then begin
  device, /close
  set_plot, 'X'
;endif

endfor	; loop over yearly files
;return
end
