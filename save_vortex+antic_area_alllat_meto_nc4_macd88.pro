;
; nc4 files
; MetO zonal mean vortex and anticyclone area (all latitudes) for 18 year record
;
; VLH 9/26/09
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
!p.background=icolmax
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd

lstmn=9L & lstdy=28L & lstyr=1991L
;lstmn=3L & lstdy=12L & lstyr=2009L
ledmn=10L & leddy=6L & ledyr=2013L
;
; Get start and end dates
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd_all=lonarr(nfile)
syyyymmdd_all=strarr(nfile)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L
dir='/Volumes/earth/aura3/data/UKMO_data/Datfiles/ukmo_'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' Starting day outside range '
      if ndays gt ledday then goto,saveit
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      uyr=string(FORMAT='(I2.2)',iyr1)
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
      print,date
      yyyymmdd_all(icount)=long(syr+smn+sdy)
      syyyymmdd_all(icount)=date
      ifile=dir+mon(imn-1)+sdy+'_'+uyr+'.nc4'
;
;***Read UKMO marker field
;
      dum1=findfile(ifile)
      if dum1(0) eq '' then goto,skip
      ncid=ncdf_open(ifile)
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      thlev=fltarr(nth)
      marksf=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,mark2
      ncdf_close,ncid
      index=where(mark2 lt -1.)
      if index(0) ne -1L then mark2(index)=-1.*mark2(index)/mark2(index)	; highs are all -1
      index=where(mark2 gt 1.)
      if index(0) ne -1L then mark2(index)=mark2(index)/mark2(index)        ; lows are all +1
;
; declare time period array. compute area.
;
      if kcount eq 0L then begin
         amarkbar_all=-9999.+0.*fltarr(nfile,nr,nth)
         vmarkbar_all=-9999.+0.*fltarr(nfile,nr,nth)
         aarea_all=-9999.+0.*fltarr(nfile,nr,nth)
         varea_all=-9999.+0.*fltarr(nfile,nr,nth)
         dum=reform(mark2(*,*,0))
         area=0.*dum
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(j,*)=dx*hy    ; area of each grid point
         endfor
         kcount=1
      endif
;
; number of anticyclone points around each latitude circle
;
      azm=0./0.*fltarr(nr,nth)	; % longitudes with anticyclones
      aazm=0./0.*fltarr(nr,nth)	; anticyclone area
      vzm=0./0.*fltarr(nr,nth)	; % longitudes with vortex
      vazm=0./0.*fltarr(nr,nth)	; vortex area
      for k=0,nth-1 do begin
          for j=0,nr-1 do begin
              marklev=reform(mark2(j,*,k))
              index=where(marklev lt 0.0,nn)
              if index(0) ne -1 then azm(j,k)=100.*total(marklev(j,index))/float(nc)
              if index(0) ne -1 then aazm(j,k)=total(area(j,index))
              index=where(marklev gt 0.0,nn)
              if index(0) ne -1 then vzm(j,k)=100.*total(marklev(j,index))/float(nc)
              if index(0) ne -1 then vazm(j,k)=total(area(j,index))
if vzm(j,k) gt 100. then stop
          endfor
      endfor
;
; check
;
nlvls=21
col1=1+indgen(nlvls)*icolmax/nlvls
erase
!type=2^2+2^3
contour,vzm,alat,th,yrange=[min(th),max(th)],/noerase,levels=5.*findgen(nlvls),title=date,c_color=col1,/cell_fill,color=0
contour,azm,alat,th,/noerase,levels=-100.+5.*findgen(nlvls),/follow,/overplot,color=0
;
; retain all daily zonal means 
;
      amarkbar_all(icount,*,*)=azm	; % of latitude circles filled by anticyclones
      vmarkbar_all(icount,*,*)=vzm	; % of latitude circles filled by the vortex
      aarea_all(icount,*,*)=aazm      ; anticyclone area at each latitude
      varea_all(icount,*,*)=vazm      ; vortex area at each latitude
;
; skip if file is missing but still increment day
;
      skip:
      icount=icount+1L
goto,jump
saveit:

syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
sday=strmid(syyyymmdd_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nday=icount
save,file='MetO_antic_area_alllat_'+yearlab+'_nc4.sav',nth,nr,alat,nday,yyyymmdd_all,th,amarkbar_all,vmarkbar_all,aarea_all,varea_all
end
