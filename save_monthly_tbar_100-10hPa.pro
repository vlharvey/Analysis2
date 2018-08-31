;
; save monthly mean UKMO T averaged poleward of 60 N and between 100 and 10 hPa.
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nwp
;
; Ask interactive questions- get starting/ending dates
;
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
imn=7
lstyr=2007
ledyr=2011
;read,' Enter month ',imn
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
smn=string(FORMAT='(i2.2)',imn)
syr0=string(FORMAT='(i4)',lstyr)
syr1=string(FORMAT='(i4)',ledyr)
if lstyr eq ledyr then yearlab=syr0
if lstyr ne ledyr then yearlab=syr0+'-'+syr1
ofile='/Volumes/earth/harvey/UKMO_data/Datfiles/ukmo_polar_Tavg_'+month(imn-1)+'_'+yearlab+'.sav'
;
; loop over days in imn over years from lstyr to ledyr
;
years=lonarr(ledyr-lstyr+1)
tmean=fltarr(ledyr-lstyr+1)
tmean_all=0.
kcount=0
for iyear=lstyr,ledyr do begin
    syr=string(FORMAT='(i4)',iyear)
    syr1=strmid(syr,2,2)
    years(iyear-lstyr)=iyear
    if iyear le 2006 then begin
       spawn,'ls /Volumes/earth/harvey/UKMO_data/Datfiles/ppassm_y'+syr1+'_m'+smn+'_d??_h12.pp.sav',ifiles
       nfiles=n_elements(ifiles)
    endif
    if iyear gt 2006 then begin
       spawn,'ls /Volumes/earth/harvey/UKMO_data/Datfiles/ukmo-nwp-strat_gbl-std_'+syr+smn+'??12_u-v-gph-t-w_uars.nc',ifiles
       nfiles=n_elements(ifiles)
    endif
;
; loop over days
;
    icount=0L
    for ifile=0,nfiles-1L do begin
        sfile=ifiles(ifile)
        print,sfile
        if iyear le 2006 then restore,sfile
        iflg=0
        if iyear gt 2006 then rd_ukmo_nwp,sfile,nc,nr,nc1,nr1,nlv,wlon,alon,wlat,alat,p,z3d,t3d,u3d,v3d,iflg
        if iflg ne 0 then goto,jumpday
;
; Average temperature from 100-10 hPa and from 40S-60S (July) and 60N-80N (Jan)
;
         if imn eq 1 then yindex=where(alat ge 60. and alat le 80.)
         if imn eq 7 then yindex=where(alat le -40. and alat ge -60.)
         zindex=where(p le 100. and p ge 10.)
         tmean_all=tmean_all+mean(t3d(*,yindex,zindex))
         tmean(iyear-lstyr)=tmean(iyear-lstyr)+mean(t3d(*,yindex,zindex))
         icount=icount+1L
         kcount=kcount+1L
      endfor
      tmean(iyear-lstyr)=tmean(iyear-lstyr)/float(icount)
      jumpday:
endfor
tmean_all=tmean_all/float(kcount)
tprime=tmean-tmean_all
;
; save file
;
ofile='/Volumes/earth/harvey/UKMO_data/Datfiles/ukmo_polar_Tavg_'+month(imn-1)+'_'+yearlab+'.sav'
save,filename=ofile,tmean_all,tmean,tprime,years
end
