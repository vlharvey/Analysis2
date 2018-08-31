;
; save multi-year zonal mean MLS quantities
;
@stddat
@kgmt
@ckday
@kdate
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,0.8*cos(a),0.8*sin(a),/fill
setplot='x'
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=600 & nydim=600
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.6
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
dir='/aura6/data/MLS_data/Datfiles_SOSST/'
spawn,'ls '+dir+'o3_mls_v2.2_*sav',o3files
spawn,'ls '+dir+'h2o_mls_v2.2_*sav',h2ofiles
spawn,'ls '+dir+'tpd_mls_v2.2_*sav',tpfiles
spawn,'ls '+dir+'cat_mls_v2.2_*sav',catfiles
help,o3files,catfiles
;
; extract month numbers from filenames
;
wmon=strarr(n_elements(o3files))
for n=0L,n_elements(o3files)-1L do begin
    result=strsplit(o3files(n),'_',/extract)
    wmon(n)=strmid(result(5),4,2)
endfor
spawn,'/usr/bin/rm '+dir+'*-*'
th=[5000.,4800.,4600.,4400.,4200.,4000.,3800.,3600.,3400.,3200.,3000.,$
2800.,2600.,2400.,2200.,2000.,1800.,1600.,1400.,1200.,1000.,900.,800.,$
700.,600.,550.,500.,450.,400.,350.]
nth=n_elements(th)
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
for m=0,nmonths-1 do begin
;
; pull out days in this month from all years
;
    smon=string(format='(i2.2)',m+1L)
    index=where(wmon eq smon,nfile)
    ifiles=o3files(index)
    cfiles=catfiles(index)
    tfiles=tpfiles(index)
    hfiles=h2ofiles(index)
    for n=0L,nfile-1 do begin
        ifile=ifiles(n)
        cfile=cfiles(n)
        tfile=tfiles(n)
        hfile=hfiles(n)
        print,cfile
;
; read MLS data
;
        restore,cfile
        restore,tfile
        temp=temperature & tpmask=TEMPERATURE_MASK
        restore,ifile
        o3=mix & o3mask=mask
        restore,hfile
        h2o=mix & h2omask=mask
;help,temp,tpmask,h2o,h2omask,o3,o3mask
;
; eliminate data where mask is negative
;
        x=where(tpmask eq -99.)
        if x(0) ne -1L then temp(x)=-99.
        x=where(h2omask eq -99.)
        if x(0) ne -1L then h2o(x)=-99.
        x=where(o3mask eq -99.)
        if x(0) ne -1L then o3(x)=-99.
;
; bin MLS data in 5 degree latitude bins
;
        nprof=n_elements(time)
        nr=72L
        nz=n_elements(altitude)
        deltay=2.5
        latbin=-88.75+deltay*findgen(nr)
        alat=latbin
        tempyz=fltarr(nr,nz)
        ntempyz=fltarr(nr,nz)
        pyz=fltarr(nr,nz)
        npyz=fltarr(nr,nz)
        h2oyz=fltarr(nr,nz)
        nh2oyz=lonarr(nr,nz)
        o3yz=fltarr(nr,nz)
        no3yz=lonarr(nr,nz)
        for i=0L,nprof-1L do begin
            y0=latitude(i)
            tempprof=reform(temp(i,*))
            pprof=reform(pressure(i,*))
            h2oprof=reform(h2o(i,*))
            o3prof=reform(o3(i,*))
;help,tempprof,h2oprof,o3prof
            for j=0L,nr-1L do begin
                if latbin(j)-deltay/2. le y0 and latbin(j)+deltay/2. gt y0 then begin
        
                   index=where(tempprof ne -99.)
                   if index(0) ne -1L then begin
                      tempyz(j,index)=tempyz(j,index)+tempprof(index)
                      ntempyz(j,index)=ntempyz(j,index)+1L
                   endif
                   index=where(pprof ne -99.)
                   if index(0) ne -1L then begin
                      pyz(j,index)=pyz(j,index)+pprof(index)
                      npyz(j,index)=npyz(j,index)+1L
                   endif
                   index=where(h2oprof ne -99.)
                   if index(0) ne -1L then begin
                      h2oyz(j,index)=h2oyz(j,index)+h2oprof(index)
                      nh2oyz(j,index)=nh2oyz(j,index)+1L
                   endif
                   index=where(o3prof ne -99.)
                   if index(0) ne -1L then begin
                      o3yz(j,index)=o3yz(j,index)+o3prof(index)
                      no3yz(j,index)=no3yz(j,index)+1L
                   endif

                endif
            endfor
        endfor
;
; average contents of each bin
;
        index=where(ntempyz gt 0L)
        if index(0) ne -1L then tempyz(index)=tempyz(index)/float(ntempyz(index))
        index=where(ntempyz eq 0L)
        if index(0) ne -1L then tempyz(index)=-99.
        index=where(tempyz ne -99.)
;if index(0) ne -1L then print,'Temp ',min(tempyz(index)),max(tempyz(index)),max(ntempyz)

        index=where(npyz gt 0L)
        if index(0) ne -1L then pyz(index)=pyz(index)/float(npyz(index))
        index=where(npyz eq 0L)
        if index(0) ne -1L then pyz(index)=-99.
        index=where(pyz ne -99.)
;if index(0) ne -1L then print,'Pressure ',min(pyz(index)),max(pyz(index)),max(npyz)

        index=where(nh2oyz gt 0L)
        if index(0) ne -1L then h2oyz(index)=h2oyz(index)/float(nh2oyz(index))
        index=where(nh2oyz eq 0L)
        if index(0) ne -1L then h2oyz(index)=-99.
        index=where(h2oyz ne -99.)
;if index(0) ne -1L then print,'H2O ',min(h2oyz(index)),max(h2oyz(index)),max(nh2oyz)

        index=where(no3yz gt 0L)
        if index(0) ne -1L then o3yz(index)=o3yz(index)/float(no3yz(index))
        index=where(no3yz eq 0L)
        if index(0) ne -1L then o3yz(index)=-99.
        index=where(o3yz ne -99.)
;if index(0) ne -1L then print,'O3 ',min(o3yz(index)),max(o3yz(index)),max(no3yz)
;
; check daily zonal means
;
;!type=2^2+2^3
;set_viewport,.1,.45,.55,.9
;contour,tempyz,alat,altitude,nlevels=20,color=0,/noeras
;set_viewport,.55,.9,.55,.9
;contour,pyz,alat,altitude,nlevels=20,color=0,/noeras
;set_viewport,.1,.45,.1,.45
;contour,h2oyz*1.e6,alat,altitude,levels=findgen(15),color=0,/noeras
;set_viewport,.55,.9,.1,.45
;contour,o3yz*1.e6,alat,altitude,levels=findgen(15),color=0,/noeras
;stop
;
; summate
;
        if n eq 0 then begin
           h2oyz_mean=0.*h2oyz
           nh2oyz_mean=0L*nh2oyz
           o3yz_mean=0.*o3yz
           no3yz_mean=0L*no3yz
           tempyz_mean=0.*tempyz
           ntempyz_mean=0L*ntempyz
           pyz_mean=0.*pyz
           npyz_mean=0L*npyz
        endif
        index=where(h2oyz ne 99.)
        if index(0) ne -1L then begin
           h2oyz_mean(index)=h2oyz_mean(index)+h2oyz(index)
           nh2oyz_mean(index)=nh2oyz_mean(index)+1L
        endif
        index=where(o3yz ne 99.)
        if index(0) ne -1L then begin
           o3yz_mean(index)=o3yz_mean(index)+o3yz(index)
           no3yz_mean(index)=no3yz_mean(index)+1L
        endif
        index=where(tempyz ne 99.)
        if index(0) ne -1L then begin
           tempyz_mean(index)=tempyz_mean(index)+tempyz(index)
           ntempyz_mean(index)=ntempyz_mean(index)+1L
        endif
        index=where(pyz ne 99.)
        if index(0) ne -1L then begin
           pyz_mean(index)=pyz_mean(index)+pyz(index)
           npyz_mean(index)=npyz_mean(index)+1L
        endif
    endfor          ; loop over files
;
; divide by total number of days summed (except for PV)
;
    index=where(nh2oyz_mean gt 0L)
    if index(0) ne -1L then h2oyz_mean(index)=h2oyz_mean(index)/nh2oyz_mean(index)
    index=where(no3yz_mean gt 0L)
    if index(0) ne -1L then o3yz_mean(index)=o3yz_mean(index)/no3yz_mean(index)
    index=where(ntempyz_mean gt 0L)
    if index(0) ne -1L then tempyz_mean(index)=tempyz_mean(index)/ntempyz_mean(index)
    index=where(npyz_mean gt 0L)
    if index(0) ne -1L then pyz_mean(index)=pyz_mean(index)/npyz_mean(index)
;
; check monthly mean zonal means
;
erase
!type=2^2+2^3
set_viewport,.1,.45,.55,.9
contour,tempyz_mean,alat,altitude,nlevels=20,color=0,/noeras
set_viewport,.55,.9,.55,.9
contour,pyz_mean,alat,altitude,nlevels=20,color=0,/noeras
set_viewport,.1,.45,.1,.45
contour,h2oyz_mean*1.e6,alat,altitude,levels=findgen(15),color=0,/noeras
set_viewport,.55,.9,.1,.45
contour,o3yz_mean*1.e6,alat,altitude,levels=findgen(15),color=0,/noeras
;
; interpolate to theta surfaces
;
    h2oyth_mean=fltarr(nr,nth)
    o3yth_mean=fltarr(nr,nth)
    tempyth_mean=fltarr(nr,nth)
    pyth_mean=fltarr(nr,nth)
    for k=0L,nth-1L do begin
        thlev=th(k)
        for j=0L,nr-1L do begin
            o3prof=reform(o3yz_mean(j,*))
            h2oprof=reform(h2oyz_mean(j,*))
            tpprof=reform(tempyz_mean(j,*))
            pprof=reform(pyz_mean(j,*))
            thprof=0.*tpprof
            index=where(pprof ne -99.)
            if index(0) eq -1L then goto,jumplat
            thprof(index)=tpprof(index)*(1000./pprof(index))^0.286
            for kk=0L,nz-2L do begin
                if thprof(kk) ne -99. and thprof(kk+1) ne -99. and $
                   thprof(kk) lt thlev and thprof(kk+1) gt thlev then begin
                   zscale=(thlev-thprof(kk))/(thprof(kk+1)-thprof(kk))
                   if tpprof(kk) ne -99. and tpprof(kk+1) ne -99. then $
                      tempyth_mean(j,k)=tpprof(kk)+zscale*(tpprof(kk+1)-tpprof(kk))
                   if pprof(kk) ne -99. and pprof(kk+1) ne -99. then $
                      pyth_mean(j,k)=pprof(kk)+zscale*(pprof(kk+1)-pprof(kk))
                   if h2oprof(kk) ne -99. and h2oprof(kk+1) ne -99. then $
                      h2oyth_mean(j,k)=h2oprof(kk)+zscale*(h2oprof(kk+1)-h2oprof(kk))
                   if o3prof(kk) ne -99. and o3prof(kk+1) ne -99. then $
                      o3yth_mean(j,k)=o3prof(kk)+zscale*(o3prof(kk+1)-o3prof(kk))
                endif
            endfor
jumplat:
        endfor
    endfor
;
; write monthly means
;
    save,file=dir+mon(m)+'avg_yth.sav',nr,nth,alat,th,$
         pyth_mean,h2oyth_mean,o3yth_mean,tempyth_mean
endfor		; loop over months
end
