;
; Hovmoller of MLS water vapor 1 May to 1 Sepp
; overplot CIPS frequency
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.25]
yorig=[0.25]
xlen=0.5
ylen=0.7
cbaryoff=0.15
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
start_year=[2007,2008,2009,2010,2011,2012,2013,2014]
start_date=[-27, -21, -24, -24, -26, -27, -34, -28]
end_date=[66, 65, 61, 61, 64, 61, 64, 80]
nyear=n_elements(start_year)
for nn=start_year(0),start_year(nyear-1L) do begin

lstmn=5         ; NH
lstdy=1
lstyr=nn
ledmn=9         ; NH
leddy=1
ledyr=nn
lstday=0
ledday=0
;
; loop over years
;
for iyear=lstyr,ledyr do begin
for ilat=75,75 do begin ;50,80,5 do begin

syr=string(format='(i2.2)',iyear-2000)
slt=strcompress(ilat,/remove_all)
restore,'/Volumes/Data/CIPS_data/Pre_process/Save_files/xt_cips_3c_NH'+syr+'_DSC_'+slt+'Lat_2G_freq.sav
decfreq=XTFREQ
restore,'/Volumes/Data/CIPS_data/Pre_process/Save_files/xt_cips_3c_NH'+syr+'_ASC_'+slt+'Lat_2G_freq.sav
ascfreq=XTFREQ
cipsfreq=0.*decfreq
index=where(decfreq ne 0. and ascfreq eq 0.)
if index(0) ne -1L then cipsfreq(index)=decfreq(index)
index=where(decfreq eq 0. and ascfreq ne 0.)
if index(0) ne -1L then cipsfreq(index)=ascfreq(index)
index=where(decfreq ne 0. and ascfreq ne 0.)
if index(0) ne -1L then cipsfreq(index)=(decfreq(index)+ascfreq(index))/2.

kcount=0
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

z = stddat(lstmn,lstdy,iyear,lstday)
z = stddat(ledmn,leddy,iyear,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
;
; Compute initial Julian date
;
iyr = iyear
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; longitude grid
;
dx=15.
nc=long(360./dx)+1
longrid=dx*findgen(nc)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      sdate_all(icount)=sdate
      print,sdate
;
; restore MLS on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v3.3_'+sdate+'.sav'
;
; apply mask
;
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlsh2omix=mix
      good=where(mlsh2omix ne -99.)
      if good(0) ne -1L then mlsh2omix(good)=mlsh2omix(good)*1.e6
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
;        print,altitude
         ralt=83.
;        read,'Enter desired altitude ',ralt
         index=where(abs(altitude-ralt) eq min(abs(altitude-ralt)))
         ilev=index(0)
         salt=strcompress(long(ralt),/remove_all)
;        rlat=75.
;        read,'Enter desired latitude',rlat
;        slat=strcompress(long(rlat),/remove_all)
         mlspolarh2o_xt=fltarr(nc,kday)
         kcount=1
goto,quick
      endif
;
; compute polar water
;
      mlspolarh2o=fltarr(nc)
      mlsnh2oprof=lonarr(nc)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge rlat-2. and latitude(ii) le rlat+2. then begin
;         if latitude(ii) ge 80. then begin
             for i=0L,nc-2L do begin
                 if longitude(ii) ge longrid(i) and longitude(ii) le longrid(i+1) and mlsh2omix(ii,ilev) ne -99. then begin
                    mlspolarh2o(i)=mlspolarh2o(i)+mlsh2omix(ii,ilev)
                    mlsnh2oprof(i)=mlsnh2oprof(i)+1L
                 endif
             endfor
          endif
      endfor
      good=where(mlsnh2oprof gt 0L)
      if good(0) ne -1L then mlspolarh2o(good)=mlspolarh2o(good)/float(mlsnh2oprof(good))
      mlspolarh2o_xt(*,icount)=mlspolarh2o
skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; wrap around point
;
mlspolarh2o_xt(nc-1,*)=mlspolarh2o_xt(0,*)
;
; interpolate small gaps in time
;
for k=0,nc-1 do begin
    dlev=reform(mlspolarh2o_xt(k,*))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump2
               endif
           endfor
jump2:
        endif
    endfor
    mlspolarh2o_xt(k,*)=dlev
endfor
;
; year date label
;
syear=strmid(sdate_all,0,4)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
; save 
;
save,file='xt_mls_h2o_'+yearlab+'_'+slat+'_pmc_full_season.sav',mlspolarh2o_xt,kday,longrid,sdate_all
quick:
yearlab=string(format='(i4)',iyear)
restore,'xt_mls_h2o_'+yearlab+'_'+slat+'_pmc_full_season.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
;xindex=where(sday eq '01' or sday eq '05' or sday eq '10' or sday eq '15' or sday eq '20' or sday eq '25',nxticks)
xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='xt_mls_h2o+cips_freq_'+yearlab+'_'+slat+'_pmc_full_season.ps'
   !p.charsize=1.55
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean 
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolarh2o_xt eq 0.)
if index(0) ne -1L then mlspolarh2o_xt(index)=0./0.
mlspolarh2o_xt=smooth(mlspolarh2o_xt,3,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarh2o_xt(index)=0./0.

tlevel=0.5+.25*findgen(27)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlspolarh2o_xt,longrid,findgen(kday),/noeras,yrange=[kday-1,0.],xrange=[0.,360.],$
      charsize=2,color=0,xtitle='Longitude',/cell_fill,c_color=col1,charthick=2,$
      levels=tlevel,yticks=nxticks-1,ytickname=xlabs,ytickv=xindex,min_value=-99.
contour,mlspolarh2o_xt,longrid,findgen(kday),levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
print,min(mlspolarh2o_xt),max(mlspolarh2o_xt)
;contour,mlspolarh2o_xt,longrid,findgen(kday),levels=[150.],color=mcolor,/follow,/overplot,c_labels=[1],thick=4
;axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[-6,-51.],/save,charsize=2,color=0,charthick=2
;axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[-11,-42.],/save,charsize=2,color=0,charthick=2    ; May 10 to June 10
axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[74.,-51.],/save,charsize=2,color=0,charthick=2,ytitle='DFS'    ; May 1 to Sep 1
level=[1,5,10,20,30,40,50,60,70,80,90,100]
;contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=level,color=mcolor,/overplot,/noeras,$
;        c_charsize=2,c_labels=1+0*level,thick=10
index=where(iyr eq start_year)
loadct,0
plots,0.,start_date(index)
plots,360.,start_date(index),/continue,color=0,thick=20
xyouts,xmn+0.01,ymx-0.05,yearlab,/normal,charsize=3,color=mcolor,charthick=2
loadct,39

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Water Vapor (ppmv) '+slat+'N/'+salt+'km',charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim xt_mls_h2o+cips_freq_'+yearlab+'_'+slat+'_pmc_full_season.ps -rotate -90 xt_mls_h2o+cips_freq_'+yearlab+'_'+slat+'_pmc_full_season.png'
    endif
endfor  ; loop over latitudes
endfor  ; loop over years


endfor
end
