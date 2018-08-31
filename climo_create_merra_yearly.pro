;
; create MERRA yearly files
;
@stddat
@kgmt
@ckday
@kdate
;
; get WACCM latitude grid
;
restore,'/Volumes/Data/WACCM/WACCM4/mee00fpl_FW2/mee00fpl_FW2.cam2.h3.Year12_1002_Q.sav'
nr=n_elements(lat)
latbin=LAT
dy=latbin(1)-latbin(0)
dirh='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
model_years=1979+indgen(36)
model_years=[2015]
model_years=string(FORMAT='(i4)',long(model_years))
nyears=n_elements(model_years)
;
; loop over model years
;
for iyear=0L,nyears-1L do begin

lstmn=1 & lstdy=1 & lstyr=1995 & lstday=0
ledmn=12 & leddy=31 & ledyr=1995 & ledday=0	; choose any non leap year
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
kday=ledday-lstday+1L
if kday ne 365L then stop,'check kday'
sdate_all=strarr(kday)
icount=0
kcount=0
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if icount eq kday then goto,saveit
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      sdate=smn+sdy
      sdate_all(icount)=sdate
      if sdate eq '0229' then stop,'check leap year?'	; there is no leap day in WACCM
;
; read data
;
      dum=findfile(dirh+model_years(iyear)+sdate+'.sav')
      if dum(0) eq '' then goto,jumpday
      restore,dum(0)
;
; declare yearly average and sigma arrays
;
      if kcount eq 0L then begin
         latitude=latitude_waccm
         nl=n_elements(pressure)
         UGRD_AVG=fltarr(nr,nl,kday)
         VGRD_AVG=fltarr(nr,nl,kday)
         GPGRD_AVG=fltarr(nr,nl,kday)
         QVGRD_AVG=fltarr(nr,nl,kday)
         TGRD_AVG=fltarr(nr,nl,kday)
         kcount=1L
      endif
;
; loop over years and retain all data
;
      print,'restored '+dum(0)
      UGRD_AVG(*,*,icount)=mean(UGRD,dimension=1)
      VGRD_AVG(*,*,icount)=mean(VGRD,dimension=1)
      GPGRD_AVG(*,*,icount)=mean(ZGRD,dimension=1)
      QVGRD_AVG(*,*,icount)=mean(QVGRD,dimension=1)
      TGRD_AVG(*,*,icount)=mean(TGRD,dimension=1)
print,max(TGRD_AVG(*,*,icount)),max(QVGRD_AVG(*,*,icount)),max(GPGRD_AVG(*,*,icount)),max(VGRD_AVG(*,*,icount)),max(UGRD_AVG(*,*,icount))
jumpday:
      icount=icount+1L
goto,jump
;
; save yearly file
;
saveit:
;
; convert specific humidity q to mixing ratio w (w=q/(1-q))
;
H2OGRD_AVG=QVGRD_AVG/(1.0-QVGRD_AVG)
;
; check
;
erase
!type=2^2+2^3
loadct,39
mcolor=byte(!p.color)
mcolor=fix(mcolor)
device,decompose=0
if mcolor eq 0 then mcolor=255
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
ilat=where(min(abs(latbin-0.94736842)) eq abs(latbin-0.94736842))
plotarray=transpose(reform(TGRD_AVG(ilat,*,*)))
plotarray2=transpose(reform(H2OGRD_AVG(ilat,*,*)))
omin=min(plotarray)
omax=max(plotarray)
level=omin+((omax-omin)/nlvls)*findgen(nlvls+1)
index=where(sdate_all  ne '')
contour,plotarray,index,pressure,levels=level,c_color=col1,/cell_fill,/noeras,yrange=[max(pressure),min(pressure)],/ylog,$
        xrange=[1.,kday],xticks=6,ytitle='Pressure',xtitle='DOY',charsize=1.5,min_value=-99.,title='MERRA Temp '+model_years(iyear)
contour,plotarray,index,pressure,levels=level,c_color=0,/follow,/noeras,/overplot
contour,plotarray2*1.e6,index,pressure,levels=findgen(10),c_color=mcolor,/follow,/noeras,/overplot,thick=3
;
; save yearly file
;
ofile=dirh+model_years(iyear)+'_daily_zm.sav'
save,file=ofile,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
endfor	; loop over years
end
