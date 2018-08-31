;
; Save Average Composite SSW days -30 to +30
;
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
;
; SSW day zeros
;
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
kcount=0L
sswdates=[$
'20010206',$    ; for now choose six major SSW w/o ES from the most recent years for comparison to ES events since 2004
'20011229',$
'20020216',$
'20030117',$
'20070223',$
'20080222']
nevents=n_elements(sswdates)
for iES = 0L, nevents - 1L do begin
    sevent=string(format='(i2.2)',ies+1)
    icount=0L
    kdays=61

    sswdate0=sswdates(ies)
    iyr=long(strmid(sswdate0,0,4))
    imn=long(strmid(sswdate0,4,2))
    idy=long(strmid(sswdate0,6,2))
    jday = JULDAY(imn,idy,iyr)
;goto,plotit
koff=30
    jday0=jday-koff
    jday1=jday+koff
    CALDAT, jday0, lstmn ,lstdy , lstyr
    CALDAT, jday1, ledmn ,leddy , ledyr

lstday=0L & ledday=0L
if lstyr eq ledyr then yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
;goto,quick
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=long(ledday-lstday+1L)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
erase
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,colorbar
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
print,sdate,icount-30
;
; read daily file
;
        dum=findfile(dir+sdate+'.nc3')
        if dum ne '' then ncfile0=dir+sdate+'.nc3'
        rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        if iflag ne 0L then goto,jumpday
        tmp2=0.*p2
        for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286

        if kcount eq 0L then begin
           pv2_comp=fltarr(nr,nc,nth,kdays)
           p2_comp=fltarr(nr,nc,nth,kdays)
           u2_comp=fltarr(nr,nc,nth,kdays)
           v2_comp=fltarr(nr,nc,nth,kdays)
           qdf2_comp=fltarr(nr,nc,nth,kdays)
           q2_comp=fltarr(nr,nc,nth,kdays)
           qv2_comp=fltarr(nr,nc,nth,kdays)
           z2_comp=fltarr(nr,nc,nth,kdays)
           sf2_comp=fltarr(nr,nc,nth,kdays)
           mark2_comp=fltarr(nr,nc,nth,kdays)
           kcount=1L
        endif
;
; retain composite days
;
        pv2_comp(*,*,*,icount)=pv2_comp(*,*,*,icount)+pv2
        p2_comp(*,*,*,icount)=p2_comp(*,*,*,icount)+p2
        u2_comp(*,*,*,icount)=u2_comp(*,*,*,icount)+u2
        v2_comp(*,*,*,icount)=v2_comp(*,*,*,icount)+v2
        qdf2_comp(*,*,*,icount)=qdf2_comp(*,*,*,icount)+qdf2
        q2_comp(*,*,*,icount)=q2_comp(*,*,*,icount)+q2
        z2_comp(*,*,*,icount)=z2_comp(*,*,*,icount)+z2
        qv2_comp(*,*,*,icount)=qv2_comp(*,*,*,icount)+qv2
        sf2_comp(*,*,*,icount)=sf2_comp(*,*,*,icount)+sf2
        mark2_comp(*,*,*,icount)=mark2_comp(*,*,*,icount)+mark2
        icount=icount+1L

jumpday:

goto,jump
colorbar:
endfor          ; loop over SSW events
pv2_comp=pv2_comp/float(nevents)
p2_comp=p2_comp/float(nevents)
u2_comp=u2_comp/float(nevents)
v2_comp=v2_comp/float(nevents)
qdf2_comp=qdf2_comp/float(nevents)
q2_comp=q2_comp/float(nevents)
qv2_comp=qv2_comp/float(nevents)
z2_comp=z2_comp/float(nevents)
sf2_comp=sf2_comp/float(nevents)
mark2_comp=mark2_comp/float(nevents)
save,filename='/Users/harvey/Harvey_etal_2014/Post_process/SSW_composite_days_merra_dm30-dp30.sav',alon,alat,th,nc,nr,nth,pv2_comp,p2_comp,u2_comp,v2_comp,qdf2_comp,q2_comp,qv2_comp,z2_comp,sf2_comp,mark2_comp
end
