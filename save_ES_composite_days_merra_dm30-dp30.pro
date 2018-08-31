;
; Save Average Composite ES days -30 to +30
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
; Read ES day zeros
;
restore, '/Users/harvey/Harvey_etal_2014/Post_process/MLS_ES_daily_max_T_Z.sav'
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
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
        dum=findfile(dir+sdate+'.nc3')
        if dum ne '' then ncfile0=dir+sdate+'.nc3'
        rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
           u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
        if iflag ne 0L then goto,jumpday
        tmp2=0.*p2
        for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286

        if kcount eq 0L then begin
           pv2_comp=fltarr(nr,nc,nth,ndays)
           p2_comp=fltarr(nr,nc,nth,ndays)
           u2_comp=fltarr(nr,nc,nth,ndays)
           v2_comp=fltarr(nr,nc,nth,ndays)
           qdf2_comp=fltarr(nr,nc,nth,ndays)
           q2_comp=fltarr(nr,nc,nth,ndays)
           qv2_comp=fltarr(nr,nc,nth,ndays)
           z2_comp=fltarr(nr,nc,nth,ndays)
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
        z2_comp(*,*,*,icount)=z2_comp(*,*,*,icount)+z2
        qv2_comp(*,*,*,icount)=qv2_comp(*,*,*,icount)+qv2
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
qv2_comp=qv2_comp/float(nevents)
z2_comp=z2_comp/float(nevents)
sf2_comp=sf2_comp/float(nevents)
mark2_comp=mark2_comp/float(nevents)
save,filename='/Users/harvey/Harvey_etal_2014/Post_process/ES_composite_days_merra_dm30-dp30.sav',alon,alat,th,nc,nr,nth,pv2_comp,p2_comp,u2_comp,v2_comp,qdf2_comp,q2_comp,qv2_comp,z2_comp,sf2_comp,mark2_comp
end
