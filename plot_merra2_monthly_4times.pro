;
; plot MERRA2 3D monthly means 4 UT
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
mcolor=fix(mcolor)
device,decompose=0
if mcolor eq 0 then mcolor=255
!type=2^2+2^3
!p.charthick=2
!p.charsize=2
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
model_years=1990+indgen(26)	; 1990-2015
nyear=n_elements(model_years)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmonth=n_elements(smon)
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
stime=['00','06','12','18']
ntime=n_elements(stime)
;
; loop over model years
;
for itime=0L,ntime-1L do begin
for imon=0L,0L do begin	;nmonth-1L do begin
for iyear=model_years(0),model_years(nyear-1L) do begin
    if iyear mod 4 eq 0L then mday(1)=29
    if iyear mod 4 ne 0L then mday(1)=28
;
; read monthly file
;
    ofile=dir+'AVG_'+strcompress(string(FORMAT='(I4,I2.2,A1,A2)',iyear,imon+1,'_',stime(itime)))+'.sav'
    restore,ofile	;,alat,alon,th,pv2avg,p2avg,u2avg,v2avg,qdf2avg,mark2avg,qv2avg,z2avg,sf2avg,q2avg,o32avg

if iyear eq model_years(0) then begin
   p2all=0.*p2avg
   u2all=0.*u2avg
   v2all=0.*v2avg
   mark2all=0.*mark2avg
   z2all=0.*z2avg
endif
p2all=p2all+p2avg
u2all=u2all+u2avg
v2all=v2all+v2avg
mark2all=mark2all+mark2avg
z2all=z2all+z2avg

endfor  ; loop over years
;
; average
;
p2all=p2all/float(nyear)
u2all=u2all/float(nyear)
v2all=v2all/float(nyear)
mark2all=mark2all/float(nyear)
z2all=z2all/float(nyear)

if itime eq 0L then erase
!type=2^2+2^3
zindex=where(th eq 4000.)
ubar=mean(u2all,dim=2)
ubar4000=reform(ubar(*,zindex(0)))
zplotarray=mean(z2avg,dim=2)
if itime eq 0L then plot,ubar4000,alat,/noeras,thick=8,color=100,xrange=[30,50],ytitle='Latitude',yrange=[20,50],xtitle='Ubar at 70 km',title=string(imon+1)
if itime eq 0L then oplot,ubar4000,alat,thick=8,color=100
if itime eq 0L then xyouts,.85,.9,stime(itime)+'Z',/normal,color=100

if itime eq 1L then oplot,ubar4000,alat,thick=8,color=150
if itime eq 1L then xyouts,.85,.85,stime(itime)+'Z',/normal,color=150

if itime eq 2L then oplot,ubar4000,alat,thick=8,color=200
if itime eq 2L then xyouts,.85,.8,stime(itime)+'Z',/normal,color=200

if itime eq 3L then oplot,ubar4000,alat,thick=8,color=250
if itime eq 3L then xyouts,.85,.75,stime(itime)+'Z',/normal,color=250


;
endfor	; loop over months
endfor	; loop over UT times
end
