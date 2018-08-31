;
; check .dat file
; /aura7/harvey/WACCM_data/Datfiles/Datfiles_Liu
; WA3548T08CO_2x.cam2.h2.2002-01-28-00000.nc files from Hanli
; read WACCM pressure files and store in binary format at
; the horizontal and vertical resolution of the UARS MetO data
; with extension to ~140 km
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
!p.background=mcolor
nlvls=30
col1=1+indgen(nlvls)*mcolor/nlvls
!noeras=1
;
; Ask interactive questions- get starting/ending dates
;
lstmn=1 & lstdy=1 & lstyr=2002
ledmn=1 & leddy=1 & ledyr=2002
lstday=0 & ledday=0
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_Liu/WAX3548T08CO_2x.cam2.h2.'
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
;
; UARS Met Office grid (+ 24 pressure levels above 0.1 hPa)
;
ncold=96L
alonold=3.75*findgen(ncold)
wlonold=1.875+3.75*findgen(ncold)
nrold=72L
nr1old=73L
alatold=90.-2.5*findgen(nr1old)
wlatold=88.75-2.5*findgen(nrold)
p=[1000.,681.3,464.2,316.2,215.4,146.8,100.,68.13,46.42,31.62,21.54,14.68,10.,6.813,$
   4.642,3.162,2.154,1.468,1.,0.681,0.464,0.316,0.215,0.147,0.1,0.0681,0.0464,0.0316,0.0215,$
   0.0147,0.01,0.00681,0.00464,0.00316,0.00215,0.00147,0.001,0.000681,0.000464,0.000316,0.000215,$
   0.000147,0.0001,0.0000681,0.0000464,0.0000316,0.0000215,0.0000147,0.00001]
nlvold=n_elements(p)
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(i4)',iyr)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      date=syr+smn+sdy
      sdate=syr+'-'+smn+'-'+sdy
;
; read WACCM data
;
      ofile=dir+sdate+'-00000.dat'
;
; loop over pressure levels and read /f77 binary to match MetO format
;
      qnew=fltarr(ncold,nr1old,nlvold)
      znew=fltarr(ncold,nr1old,nlvold)
      tnew=fltarr(ncold,nr1old,nlvold)
      nonew=fltarr(ncold,nr1old,nlvold)
      conew=fltarr(ncold,nr1old,nlvold)
      o3new=fltarr(ncold,nr1old,nlvold)
      unew=fltarr(ncold,nrold,nlvold)
      vnew=fltarr(ncold,nrold,nlvold)
      q=fltarr(ncold,nr1old)
      z=fltarr(ncold,nr1old)
      t=fltarr(ncold,nr1old)
      no=fltarr(ncold,nr1old)
      co=fltarr(ncold,nr1old)
      o3=fltarr(ncold,nr1old)
      u=fltarr(ncold,nrold)
      v=fltarr(ncold,nrold)
      p=fltarr(nlvold)
      close,1
      openr,1,ofile,/f77
      for kk=0L,nlvold-1L do begin
          plevel=p(kk)
          readu,1,plevel
          p(kk)=plevel
          readu,1,u,v,t,z,q,no,co,o3
          unew(*,*,kk)=u
          vnew(*,*,kk)=v
          tnew(*,*,kk)=t
          znew(*,*,kk)=z
          qnew(*,*,kk)=q
          nonew(*,*,kk)=no
          conew(*,*,kk)=co
          o3new(*,*,kk)=o3
          print,plevel,u(0,0),v(0,0),t(0,0),z(0,0),q(0,0),no(0,0),co(0,0),o3(0,0)

erase
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
slev=strtrim(string(p(kk)),2)+' hPa'
set_viewport,0.1,0.45,0.5,0.75
grd1=reform(znew(*,*,kk))
print,'Gp ',min(grd1),max(grd1)
index=where(grd1 ne 0.)
imin=min(grd1(index)) & imax=max(grd1)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
map_set,0,180,0,/contin,/grid,title='Gp '+slev,color=0,/noeras
contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
set_viewport,0.55,0.9,0.5,0.75
grd1=reform(tnew(*,*,kk))
print,'Tp ',min(grd1),max(grd1)
index=where(grd1 ne 0.)
imin=min(grd1(index)) & imax=max(grd1)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
map_set,0,180,0,/contin,/grid,title='Tp '+slev,color=0,/noeras
contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
set_viewport,0.1,0.45,0.15,0.4
grd1=reform(nonew(*,*,kk))
print,'NO ',min(grd1),max(grd1)
imin=min(grd1) & imax=max(grd1)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
map_set,0,180,0,/contin,/grid,title='NO '+slev,color=0,/noeras
contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
set_viewport,0.55,0.9,0.15,0.4
grd1=reform(o3new(*,*,kk))
print,'o3 ',min(grd1),max(grd1)
imin=min(grd1) & imax=max(grd1)
level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
map_set,0,180,0,/contin,/grid,title='o3 '+slev,color=0,/noeras
contour,grd1,alonold,alatold,levels=level,/fill,/cell_fill,c_color=col1,/noeras,/overplot
contour,grd1,alonold,alatold,levels=level,/follow,c_color=0,/overplot,/noeras
map_set,0,180,0,/contin,/grid,/noeras,color=mcolor
stop
     endfor		; loop over altitude
     close,1

thnew=tnew
for l=0L,nlvold-1L do THNEW(*,*,L)=TNEW(*,*,L)*(1000./P(L))^.286
goto,jump
end
