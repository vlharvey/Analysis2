;
; save daily zonal mean U, V, T on p from gridded data
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
nlvls=18
col1=(findgen(nlvls)/float(nlvls))*mcolor
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_Grid/'
;start_year=[2007,2008,2009,2010,2011,2012,2013,2014,2015]
;start_date=[-27, -21, -24, -24, -26, -27, -34, -28,-42]
;end_date=[66, 65, 61, 61, 64, 61, 64, 80]
;nyear=n_elements(start_year)

lstmn=8
lstdy=13
lstyr=2004
ledmn=3
leddy=10
ledyr=2018
lstday=0
ledday=0
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
;
; Compute initial Julian date
;
iyr = lstyr
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
nr=91L
latgrid=-90.+2.*findgen(nr)
dy=latgrid(1)-latgrid(0)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      if icount eq 0 then begin
         sdate_all=strarr(kday)
         dfs_all=fltarr(kday)
      endif
      sdate_all(icount)=sdate

      if iday ge 60 and iday le 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(6,21,long(syr))
      if iday lt 60 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr)-1L)
      if iday gt 258 then dfs=julday(long(smn),long(sdy),long(syr))-julday(12,21,long(syr))
      print,sdate,' ',dfs
      dfs_all(icount)=dfs
;
; restore MLS on this day
;
      dum=findfile(mdir+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'MLS_grid5_ALL_U_V_v4.2_'+sdate+'.sav'
;
; declare arrays on first day
;
      nlv=n_elements(pmls2)
      nr=n_elements(lat)
      if icount eq 0 then begin
         ubar_all=fltarr(kday,nr,nlv)-9999.
         vbar_all=fltarr(kday,nr,nlv)-9999.
         tbar_all=fltarr(kday,nr,nlv)-9999.
         zbar_all=fltarr(kday,nr,nlv)-9999.
         o3bar_all=fltarr(kday,nr,nlv)-9999.
         pressure=pmls2
      endif
      umean=mean(u,dim=4,/Nan)         ; mean over both nodes
      vmean=mean(v,dim=4,/Nan)
      tmean=mean(t,dim=4,/Nan)
      zmean=mean(gph,dim=4,/Nan)
      o3mean=mean(o3,dim=4,/Nan)
;
; compute zonal means
;
      ubar=mean(umean,dim=1,/nan)
      vbar=mean(vmean,dim=1,/nan)
      tbar=mean(tmean,dim=1,/nan)
      gpbar=mean(zmean,dim=1,/nan)/1000.
      o3bar=mean(o3mean,dim=1,/nan)
;
; convert geopotential height to geometric height
;
    ks=1.931853d-3
    ecc=0.081819
    gamma45=9.80
    rtd=double(180./!pi)
    dtr=1./rtd
    zbar=0.*gpbar
    for j=0L,nr-1L do begin
        sin2=sin( (lat(j)*dtr)^2.0 )
        numerator=1.0+ks*sin2
        denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
        gammas=gamma45*(numerator/denominator)
        r=6378.137/(1.006803-(0.006706*sin2))
        zbar(j,*)= (r*gpbar(j,*))/ ( (gammas/gamma45)*r - gpbar(j,*) )
    endfor
    index=where(finite(gpbar) ne 1)
    if index(0) ne -1L then zbar(index)=-9999.
;
;erase
;!type=2^2+2^3
;contour,tbar,lat,zbar,levels=160+10*findgen(nlvls),/noerase,c_color=col1,yrange=[10,90],xrange=[-90,90],/cell_fill,title=sdate
;contour,ubar,lat,zbar,levels=10+10*findgen(nlvls),/noerase,color=mcolor,/overplot
;contour,ubar,lat,zbar,levels=-200+10*findgen(nlvls),/overplot,c_linestyle=5
;contour,o3bar*1.e6,lat,zbar,levels=1+0.5*findgen(30),/overplot,c_linestyle=5
;
; retain all daily zonal means
;
    ubar_all(icount,*,*)=ubar
    vbar_all(icount,*,*)=vbar
    tbar_all(icount,*,*)=tbar
    zbar_all(icount,*,*)=zbar
    o3bar_all(icount,*,*)=o3bar

skipmls:
icount=icount+1L
goto,jump

saveit:
save,file='MLS_YZ_UVTO3_2004-2018.sav',sdate_all,dfs_all,vbar_all,ubar_all,tbar_all,zbar_all,o3bar_all,lat,pressure

end
