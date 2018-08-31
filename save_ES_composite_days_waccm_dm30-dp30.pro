;
; Save Average Composite ES days -30 to +30
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nc3

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
cbaryoff=0.055
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
!NOERAS=-1
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
;
; Read ES day zeros
;
restore,'/Users/harvey/Harvey_etal_2014/Post_process/elevated_strat.sav'
restore,'/Users/harvey/Harvey_etal_2014/Post_process/WACCM_ES_daily_max_T_Z.sav'
ndates = 30.*n_elements(DAYZERODATES)
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
nday = 0L
dayofES = 1L
niday = 0L
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
kcount=0L

for iES = 0L, n_elements(dayzerodates) - 1L do begin
    sevent=string(format='(i2.2)',ies+1)
    ydate = dayzerodates[iES]
    print,'Day Zero = ',ydate
    iyr=long(strmid(ydate,0,4))
    imn=long(strmid(ydate,4,2))
    idy=long(strmid(ydate,6,2))
    z = kgmt(imn,idy,iyr,kday)
    kday=kday-30
    if kday lt 0L then begin
       kday=kday+365
       iyr=iyr-1
    endif
    icount=0
    for iday=kday,kday+60L do begin
        iday0=iday
        if iday0 gt 366L then iday0=iday0-365L
        kdate,float(iday0),iyr,imn,idy
        ckday,iday0,iyr
        sdy=string(FORMAT='(i2.2)',idy)
        smn=string(FORMAT='(i2.2)',imn)
        syr=strtrim(string(iyr),2)

        if icount-30L lt 0L then sday=string(format='(i3.2)',icount-30L)
        if icount-30L ge 0L then sday=string(format='(i2.2)',icount-30L)
        ifile = syr+smn+sdy
        ifiles=file_search(dir+ifile+'_3D_dyn.nc3',count=nfile)
        if ifiles[0] eq '' then continue
        result=strsplit(ifiles(0),'.',/extract)
        result2=strsplit(result(4),'_',/extract)
        sdate=result2(0)
print,icount,' ',sday,' = ',sdate
;
; read daily file
;
        ncfile0=ifiles(0)
        ncid=ncdf_open(ncfile0)
        result0=ncdf_inquire(ncid)
        if kcount eq 0L then begin
           for idim=0,result0.ndims-1 do begin
               ncdf_diminq,ncid,idim,name,dim
               if name eq 'number_of_latitudes' then nr=dim
               if name eq 'number_of_longitudes' then nc=dim
               if name eq 'number_of_levels' then nth=dim
           endfor
           for ivar=0,result0.nvars-1 do begin
               result=ncdf_varinq(ncid,ivar)
               ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
               if result.name eq 'latitude' then alat=data
               if result.name eq 'longitude' then alon=data
               if result.name eq 'theta' then th=data
           endfor
           kcount=1L
           pv2_comp=fltarr(nr,nc,nth,61)
           p2_comp=fltarr(nr,nc,nth,61)
           u2_comp=fltarr(nr,nc,nth,61)
           v2_comp=fltarr(nr,nc,nth,61)
           qdf2_comp=fltarr(nr,nc,nth,61)
           q2_comp=fltarr(nr,nc,nth,61)
           gph2_comp=fltarr(nr,nc,nth,61)
           ttgw2_comp=fltarr(nr,nc,nth,61)
           sf2_comp=fltarr(nr,nc,nth,61)
           mark2_comp=fltarr(nr,nc,nth,61)
        endif
        for ivar=0,result0.nvars-1 do begin
            result=ncdf_varinq(ncid,ivar)
            ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
            if result.name eq 'IPV' then pv2=data
            if result.name eq 'P' then p2=data
            if result.name eq 'U' then u2=data
            if result.name eq 'V' then v2=data
            if result.name eq 'QDF' then qdf2=data
            if result.name eq 'Q' then q2=data
            if result.name eq 'GPH' then gph2=data
            if result.name eq 'TTGW' then ttgw2=data
            if result.name eq 'SF' then sf2=data
            if result.name eq 'MARK' then mark2=data
        endfor
        ncdf_close,ncid
;
; retain composite days
;
        pv2_comp(*,*,*,icount)=pv2_comp(*,*,*,icount)+pv2
        p2_comp(*,*,*,icount)=p2_comp(*,*,*,icount)+p2
        u2_comp(*,*,*,icount)=u2_comp(*,*,*,icount)+u2
        v2_comp(*,*,*,icount)=v2_comp(*,*,*,icount)+v2
        qdf2_comp(*,*,*,icount)=qdf2_comp(*,*,*,icount)+qdf2
        q2_comp(*,*,*,icount)=q2_comp(*,*,*,icount)+q2
        gph2_comp(*,*,*,icount)=gph2_comp(*,*,*,icount)+gph2
        ttgw2_comp(*,*,*,icount)=ttgw2_comp(*,*,*,icount)+ttgw2
        sf2_comp(*,*,*,icount)=sf2_comp(*,*,*,icount)+sf2
        mark2_comp(*,*,*,icount)=mark2_comp(*,*,*,icount)+mark2
        icount=icount+1L
endfor          ; loop over days -30 to +30
endfor          ; loop over ES events
pv2_comp=pv2_comp/float(n_elements(DAYZERODATES))
p2_comp=p2_comp/float(n_elements(DAYZERODATES))
u2_comp=u2_comp/float(n_elements(DAYZERODATES))
v2_comp=v2_comp/float(n_elements(DAYZERODATES))
qdf2_comp=qdf2_comp/float(n_elements(DAYZERODATES))
q2_comp=q2_comp/float(n_elements(DAYZERODATES))
gph2_comp=gph2_comp/float(n_elements(DAYZERODATES))
ttgw2_comp=ttgw2_comp/float(n_elements(DAYZERODATES))
sf2_comp=sf2_comp/float(n_elements(DAYZERODATES))
mark2_comp=mark2_comp/float(n_elements(DAYZERODATES))
save,filename='/Users/harvey/Harvey_etal_2014/Post_process/ES_composite_days_waccm_dm30-dp30.sav',alon,alat,th,nc,nr,nth,pv2_comp,p2_comp,u2_comp,v2_comp,qdf2_comp,q2_comp,gph2_comp,ttgw2_comp,sf2_comp,mark2_comp
end
