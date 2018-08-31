;
; Save Average Composite ES days -30 to +30
;
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
;
; Read ES day zeros
;
restore, '/Users/harvey/Desktop/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
dir1='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS520.MetO.'
kcount=0L

result=size(MAXHEIGHTTHETA)
nevents=result(1)
for iES = 0L, nevents - 1L do begin
    sevent=string(format='(i2.2)',ies+1)
    icount=0L
    ndays=61
    for iday =0L, ndays-1L do begin
        if iday lt 30L then sday=string(format='(I3.2)',iday-30)
        if iday ge 30L then sday=string(format='(I2.2)',iday-30)
        sdate=esdate(iday+60*ies)
        if sdate eq '' then goto,jumpday        ; missing day
print,iday-30,' ',sdate
;
; read daily file
;
        dum=findfile(dir1+sdate+'_AVG.V01.nc3')
        if dum ne '' then ncfile0=dir1+sdate+'_AVG.V01.nc3'
        if dum eq '' then ncfile0=dir+sdate+'_AVG.V01.nc3'
        rd_geos5_nc3_meto,ncfile0,nc,nr,nth,alon,alat,th,$
           pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
        if kcount eq 0L then begin
           pv2_comp=fltarr(nr,nc,nth,ndays)
           p2_comp=fltarr(nr,nc,nth,ndays)
           u2_comp=fltarr(nr,nc,nth,ndays)
           v2_comp=fltarr(nr,nc,nth,ndays)
           qdf2_comp=fltarr(nr,nc,nth,ndays)
           q2_comp=fltarr(nr,nc,nth,ndays)
           msf2_comp=fltarr(nr,nc,nth,ndays)
           vp2_comp=fltarr(nr,nc,nth,ndays)
           sf2_comp=fltarr(nr,nc,nth,ndays)
           mark2_comp=fltarr(nr,nc,nth,ndays)
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
        vp2_comp(*,*,*,icount)=vp2_comp(*,*,*,icount)+vp2
        msf2_comp(*,*,*,icount)=msf2_comp(*,*,*,icount)+msf2
        sf2_comp(*,*,*,icount)=sf2_comp(*,*,*,icount)+sf2
        mark2_comp(*,*,*,icount)=mark2_comp(*,*,*,icount)+mark2
        icount=icount+1L

jumpday:

endfor          ; loop over days -30 to +30
endfor          ; loop over ES events
pv2_comp=pv2_comp/float(nevents)
p2_comp=p2_comp/float(nevents)
u2_comp=u2_comp/float(nevents)
v2_comp=v2_comp/float(nevents)
qdf2_comp=qdf2_comp/float(nevents)
q2_comp=q2_comp/float(nevents)
msf2_comp=msf2_comp/float(nevents)
vp2_comp=vp2_comp/float(nevents)
sf2_comp=sf2_comp/float(nevents)
mark2_comp=mark2_comp/float(nevents)
save,filename='ES_composite_days_geos5_dm30-dp30.sav',alon,alat,th,nc,nr,nth,pv2_comp,p2_comp,u2_comp,v2_comp,qdf2_comp,q2_comp,msf2_comp,vp2_comp,sf2_comp,mark2_comp
end
