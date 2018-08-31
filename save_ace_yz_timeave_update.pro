;
; save ACE zonal mean quantities averaged over a specified time period
; filter out N2O values less than 1.5e-8 and N2O5 values less than 1.e-10
;
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,0.8*cos(a),0.8*sin(a),/fill
setplot='x'
;read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=600 & nydim=600
xorig=[0.15]
yorig=[0.45]
xlen=0.7
ylen=0.45
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
month='        '+['J','F','M','A','M','J','J','A','S','O','N','D',' ']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2005','2006']
lyear=long(syear)
;print,'2004   2005   2006'
kyear=2006L
;read,'Enter desired year ',kyear
x=where(kyear eq lyear)
if x(0) eq -1L then stop,'Invalid year'
syear0=syear(x)
;
; restore ACE SOSST data
;
restore,dira+'cat_ace_v2.2.'+syear0
restore,dira+'temp_ace_v2.2.'+syear0
temperature_save=t
temperatureprecision_save=mask
xx=where(t ne -99.)
print,'T ',min(t(xx)),max(t(xx))

restore,dira+'ch4_ace_v2.2.'+syear0
ch4_save=mix
ch4precision_save=mask
xx=where(mix ne -99.)
print,'CH4 ',min(mix(xx)),max(mix(xx))

restore,dira+'clono2_ace_v2.2.'+syear0
clono2_save=mix
clono2precision_save=mask
xx=where(mix ne -99.)
print,'ClONO2 ',min(mix(xx)),max(mix(xx))

restore,dira+'h2o_ace_v2.2.'+syear0
h2o_save=mix
h2oprecision_save=mask
xx=where(mix ne -99.)
print,'H2O ',min(mix(xx)),max(mix(xx))

restore,dira+'hno3_ace_v2.2.'+syear0
hno3_save=mix
hno3precision_save=mask
xx=where(mix ne -99.)
print,'HNO3 ',min(mix(xx)),max(mix(xx))

restore,dira+'n2o_ace_v2.2.'+syear0
n2o_save=mix
n2oprecision_save=mask
xx=where(mix ne -99.)
print,'N2O ',min(mix(xx)),max(mix(xx))

restore,dira+'n2o5_ace_v2.2.'+syear0
n2o5_save=mix
n2o5precision_save=mask
xx=where(mix ne -99.)
print,'N2O5 ',min(mix(xx)),max(mix(xx))

restore,dira+'no2_ace_v2.2.'+syear0
no2_save=mix
no2precision_save=mask
xx=where(mix ne -99.)
print,'NO2 ',min(mix(xx)),max(mix(xx))

restore,dira+'ccl3f_ace_v2.2.'+syear0
cfc11_save=mix
cfc11precision_save=mask
xx=where(mix ne -99.)
print,'CFC11 ',min(mix(xx)),max(mix(xx))

restore,dira+'ccl2f2_ace_v2.2.'+syear0
cfc12_save=mix
cfc12precision_save=mask
xx=where(mix ne -99.)
print,'CFC12 ',min(mix(xx)),max(mix(xx))

restore,dira+'o3_ace_v2.2.'+syear0
o3_save=mix
o3precision_save=mask
xx=where(mix ne -99.)
print,'O3 ',min(mix(xx)),max(mix(xx))

plotlat0:
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,fdoy,latitude,psym=8,color=0,yrange=[-90.,90.],xrange=[1.,366.],xticks=12,$
     xtickname=month,charsize=2,ytitle='Latitude',yticks=6,title=syear0
iday0=90L
;read,'Enter desired start day (1-365) ',iday0
plots,iday0,-90
plots,iday0,90.,/continue,thick=3,color=0
iday0logic='y'
;read,'Is start day correct? ',iday0logic
if iday0logic eq 'n' then goto,plotlat0

plotlat1:
erase
plot,fdoy,latitude,psym=8,color=0,yrange=[-90.,90.],xrange=[1.,366.],xticks=12,$
     xtickname=month,charsize=2,ytitle='Latitude',yticks=6,title=syear0
plots,iday0,-90
plots,iday0,90.,/continue,thick=3,color=0
iday1=152L
;read,'Enter desired end day (1-365) ',iday1
plots,iday1,-90
plots,iday1,90.,/continue,thick=3,color=0
iday1logic='y'
;read,'Is end day correct? ',iday1logic
if iday1logic eq 'n' then goto,plotlat1
;
; extract ACE data between sdate0 and sdate1
;
print,' '
print,'subset of dates'
xindex=where(fdoy ge iday0 and fdoy le iday1,nprof)
acedate=date
sdate0=strcompress(min(acedate(xindex)),/remove_all)
sdate1=strcompress(max(acedate(xindex)),/remove_all)
latitude_save=latitude(xindex)
temperature_save=reform(temperature_save(xindex,*))
temperatureprecision_save=reform(temperatureprecision_save(xindex,*))
x=where(temperature_save ne -99.)
print,'T ',min(temperature_save(x)),max(temperature_save(x))

ch4_save=reform(ch4_save(xindex,*))
ch4precision_save=reform(ch4precision_save(xindex,*))
x=where(ch4_save ne -99.)
print,'CH4 ',min(ch4_save(x)),max(ch4_save(x))

clono2_save=reform(clono2_save(xindex,*))
clono2precision_save=reform(clono2precision_save(xindex,*))
x=where(clono2_save ne -99.)
print,'ClONO2 ',min(clono2_save(x)),max(clono2_save(x))

h2o_save=reform(h2o_save(xindex,*))
h2oprecision_save=reform(h2oprecision_save(xindex,*))
x=where(h2o_save ne -99.)
print,'H2O ',min(h2o_save(x)),max(h2o_save(x))

hno3_save=reform(hno3_save(xindex,*))
hno3precision_save=reform(hno3precision_save(xindex,*))
x=where(hno3_save ne -99.)
print,'hno3 ',min(hno3_save(x)),max(hno3_save(x))

n2o_save=reform(n2o_save(xindex,*))
n2oprecision_save=reform(n2oprecision_save(xindex,*))
x=where(n2o_save ne -99.)
print,'N2O ',min(n2o_save(x)),max(n2o_save(x))

n2o5_save=reform(n2o5_save(xindex,*))
n2o5precision_save=reform(n2o5precision_save(xindex,*))
x=where(n2o5_save ne -99.)
print,'N2O5 ',min(n2o5_save(x)),max(n2o5_save(x))

no2_save=reform(no2_save(xindex,*))
no2precision_save=reform(no2precision_save(xindex,*))
x=where(no2_save ne -99.)
print,'NO2 ',min(no2_save(x)),max(no2_save(x))

cfc11_save=reform(cfc11_save(xindex,*))
cfc11precision_save=reform(cfc11precision_save(xindex,*))
x=where(cfc11_save ne -99.)
print,'CFC11 ',min(cfc11_save(x)),max(cfc11_save(x))

cfc12_save=reform(cfc12_save(xindex,*))
cfc12precision_save=reform(cfc12precision_save(xindex,*))
x=where(cfc12_save ne -99.)
print,'CFC12 ',min(cfc12_save(x)),max(cfc12_save(x))

o3_save=reform(o3_save(xindex,*))
o3precision_save=reform(o3precision_save(xindex,*))
x=where(o3_save ne -99.)
print,'O3 ',min(o3_save(x)),max(o3_save(x))
;
; filter masked data and data outside threshholds
;
print,' '
print,' no -99 precisions'
x=where(temperatureprecision_save eq -99.)
if x(0) ne -1L then temperature_save(x)=-99.
x=where(temperature_save ne -99.)
print,'T ',min(temperature_save(x)),max(temperature_save(x))

x=where(ch4precision_save eq -99. or ch4_save lt 1.e-10 or ch4_save gt 1.e-2)
if x(0) ne -1L then ch4_save(x)=-99.
x=where(ch4_save ne -99.)
print,'CH4 ',min(ch4_save(x)),max(ch4_save(x))

x=where(clono2precision_save eq -99. or clono2_save lt 1.e-12 or clono2_save gt 1.e-6)
if x(0) ne -1L then clono2_save(x)=-99.
x=where(clono2_save ne -99.)
print,'ClONO2 ',min(clono2_save(x)),max(clono2_save(x))

x=where(h2oprecision_save eq -99. or h2o_save lt 1.e-10 or h2o_save gt 1.e-2)
if x(0) ne -1L then h2o_save(x)=-99.
x=where(h2o_save ne -99.)
print,'H2O ',min(h2o_save(x)),max(h2o_save(x))

x=where(hno3precision_save eq -99. or hno3_save lt 1.e-15 or hno3_save gt 1.e-6)
if x(0) ne -1L then hno3_save(x)=-99.
x=where(hno3_save ne -99.)
print,'HNO3 ',min(hno3_save(x)),max(hno3_save(x))

x=where(n2oprecision_save eq -99 or n2o_save lt 1.5e-8 or n2o_save gt 1.e-5)
if x(0) ne -1L then n2o_save(x)=-99.
x=where(n2o_save ne -99.)
print,'N2O ',min(n2o_save(x)),max(n2o_save(x))

x=where(no2precision_save eq -99 or no2_save lt 1.e-12 or no2_save gt 1.e-5)
if x(0) ne -1L then no2_save(x)=-99.
x=where(no2_save ne -99.)
print,'NO2 ',min(no2_save(x)),max(no2_save(x))

x=where(cfc11precision_save eq -99)
if x(0) ne -1L then cfc11_save(x)=-99.
x=where(cfc11_save ne -99.)
print,'CFC11 ',min(cfc11_save(x)),max(cfc11_save(x))

x=where(cfc12precision_save eq -99)
if x(0) ne -1L then cfc12_save(x)=-99.
x=where(cfc12_save ne -99.)
print,'CFC12 ',min(cfc12_save(x)),max(cfc12_save(x))

x=where(n2o5precision_save eq -99 or n2o5_save lt 1.e-10 or n2o5_save gt 1.e-7)
if x(0) ne -1L then n2o5_save(x)=-99.
x=where(n2o5_save ne -99.)
print,'N2O5 ',min(n2o5_save(x)),max(n2o5_save(x))

x=where(o3precision_save eq -99 or o3_save lt 1.e-8 or o3_save gt 1.e-4)
if x(0) ne -1L then o3_save(x)=-99.
x=where(o3_save ne -99.)
print,'O3 ',min(o3_save(x)),max(o3_save(x))
;
; bin ACE data in latitude bins
;
nr=36L
nz=n_elements(altitude)
deltay=5.0
latbin=-87.5+deltay*findgen(nr)
tempyz=fltarr(nr,nz)
ntempyz=lonarr(nr,nz)
ch4yz=fltarr(nr,nz)
nch4yz=lonarr(nr,nz)
clono2yz=fltarr(nr,nz)
nclono2yz=lonarr(nr,nz)
h2oyz=fltarr(nr,nz)
nh2oyz=lonarr(nr,nz)
hno3yz=fltarr(nr,nz)
nhno3yz=lonarr(nr,nz)
n2oyz=fltarr(nr,nz)
nn2oyz=lonarr(nr,nz)
no2yz=fltarr(nr,nz)
nno2yz=lonarr(nr,nz)
cfc11yz=fltarr(nr,nz)
ncfc11yz=lonarr(nr,nz)
cfc12yz=fltarr(nr,nz)
ncfc12yz=lonarr(nr,nz)
n2o5yz=fltarr(nr,nz)
nn2o5yz=lonarr(nr,nz)
o3yz=fltarr(nr,nz)
no3yz=lonarr(nr,nz)
for i=0L,nprof-1L do begin
    y0=latitude_save(i)
    tempprof=reform(temperature_save(i,*))
    ch4prof=reform(ch4_save(i,*))
    clono2prof=reform(clono2_save(i,*))
    h2oprof=reform(h2o_save(i,*))
    hno3prof=reform(hno3_save(i,*))
    n2oprof=reform(n2o_save(i,*))
    no2prof=reform(no2_save(i,*))
    cfc11prof=reform(cfc11_save(i,*))
    cfc12prof=reform(cfc12_save(i,*))
    n2o5prof=reform(n2o5_save(i,*))
    o3prof=reform(o3_save(i,*))
    for j=0L,nr-1L do begin
        if latbin(j)-deltay/2. le y0 and latbin(j)+deltay/2. gt y0 then begin

           index=where(tempprof ne -99.)
           if index(0) ne -1L then begin
              tempyz(j,index)=tempyz(j,index)+tempprof(index)
              ntempyz(j,index)=ntempyz(j,index)+1L
           endif
           index=where(ch4prof ne -99.)
           if index(0) ne -1L then begin
              ch4yz(j,index)=ch4yz(j,index)+ch4prof(index)
              nch4yz(j,index)=nch4yz(j,index)+1L
           endif
           index=where(clono2prof ne -99.)
           if index(0) ne -1L then begin
              clono2yz(j,index)=clono2yz(j,index)+clono2prof(index)
              nclono2yz(j,index)=nclono2yz(j,index)+1L
           endif
           index=where(h2oprof ne -99.)
           if index(0) ne -1L then begin
              h2oyz(j,index)=h2oyz(j,index)+h2oprof(index)
              nh2oyz(j,index)=nh2oyz(j,index)+1L
           endif
           index=where(hno3prof ne -99.)
           if index(0) ne -1L then begin
              hno3yz(j,index)=hno3yz(j,index)+hno3prof(index)
              nhno3yz(j,index)=nhno3yz(j,index)+1L
           endif
           index=where(n2oprof ne -99. and abs(n2oprof) gt 1.e-8)
           if index(0) ne -1L then begin
              n2oyz(j,index)=n2oyz(j,index)+n2oprof(index)
              nn2oyz(j,index)=nn2oyz(j,index)+1L
           endif
           index=where(no2prof ne -99.)
           if index(0) ne -1L then begin
              no2yz(j,index)=no2yz(j,index)+no2prof(index)
              nno2yz(j,index)=nno2yz(j,index)+1L
           endif
           index=where(cfc11prof ne -99.)
           if index(0) ne -1L then begin
              cfc11yz(j,index)=cfc11yz(j,index)+cfc11prof(index)
              ncfc11yz(j,index)=ncfc11yz(j,index)+1L
           endif
           index=where(cfc12prof ne -99.)
           if index(0) ne -1L then begin
              cfc12yz(j,index)=cfc12yz(j,index)+cfc12prof(index)
              ncfc12yz(j,index)=ncfc12yz(j,index)+1L
           endif
           index=where(n2o5prof ne -99. and abs(n2o5prof) gt 1.e-10)
           if index(0) ne -1L then begin
              n2o5yz(j,index)=n2o5yz(j,index)+n2o5prof(index)
              nn2o5yz(j,index)=nn2o5yz(j,index)+1L
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
if index(0) ne -1L then print,'Temp ',min(tempyz(index)),max(tempyz(index)),max(ntempyz)

index=where(nch4yz gt 0L)
if index(0) ne -1L then ch4yz(index)=ch4yz(index)/float(nch4yz(index))
index=where(nch4yz eq 0L)
if index(0) ne -1L then ch4yz(index)=-99.
index=where(ch4yz ne -99.)
if index(0) ne -1L then print,'CH4 ',min(ch4yz(index)),max(ch4yz(index)),max(nch4yz)

index=where(nclono2yz gt 0L)
if index(0) ne -1L then clono2yz(index)=clono2yz(index)/float(nclono2yz(index))
index=where(nclono2yz eq 0L)
if index(0) ne -1L then clono2yz(index)=-99.
index=where(clono2yz ne -99.)
if index(0) ne -1L then print,'ClONO2 ',min(clono2yz(index)),max(clono2yz(index)),max(nclono2yz)

index=where(nh2oyz gt 0L)
if index(0) ne -1L then h2oyz(index)=h2oyz(index)/float(nh2oyz(index))
index=where(nh2oyz eq 0L)
if index(0) ne -1L then h2oyz(index)=-99.
index=where(h2oyz ne -99.)
if index(0) ne -1L then print,'H2O ',min(h2oyz(index)),max(h2oyz(index)),max(nh2oyz)

index=where(nhno3yz gt 0L)
if index(0) ne -1L then hno3yz(index)=hno3yz(index)/float(nhno3yz(index))
index=where(nhno3yz eq 0L)
if index(0) ne -1L then hno3yz(index)=-99.
index=where(hno3yz ne -99.)
if index(0) ne -1L then print,'HNO3 ',min(hno3yz(index)),max(hno3yz(index)),max(nhno3yz)

index=where(nn2oyz gt 0L)
if index(0) ne -1L then n2oyz(index)=n2oyz(index)/float(nn2oyz(index))
index=where(nn2oyz eq 0L)
if index(0) ne -1L then n2oyz(index)=-99.
index=where(n2oyz ne -99.)
if index(0) ne -1L then print,'N2O ',min(n2oyz(index)),max(n2oyz(index)),max(nn2oyz)

index=where(nno2yz gt 0L)
if index(0) ne -1L then no2yz(index)=no2yz(index)/float(nno2yz(index))
index=where(nno2yz eq 0L)
if index(0) ne -1L then no2yz(index)=-99.
index=where(no2yz ne -99.)
if index(0) ne -1L then print,'NO2 ',min(no2yz(index)),max(no2yz(index)),max(nno2yz)

index=where(ncfc11yz gt 0L)
if index(0) ne -1L then cfc11yz(index)=cfc11yz(index)/float(ncfc11yz(index))
index=where(ncfc11yz eq 0L)
if index(0) ne -1L then cfc11yz(index)=-99.
index=where(cfc11yz ne -99.)
if index(0) ne -1L then print,'CFC11 ',min(cfc11yz(index)),max(cfc11yz(index)),max(ncfc11yz)

index=where(ncfc12yz gt 0L)
if index(0) ne -1L then cfc12yz(index)=cfc12yz(index)/float(ncfc12yz(index))
index=where(ncfc12yz eq 0L)
if index(0) ne -1L then cfc12yz(index)=-99.
index=where(cfc12yz ne -99.)
if index(0) ne -1L then print,'CFC12 ',min(cfc12yz(index)),max(cfc12yz(index)),max(ncfc12yz)

index=where(nn2o5yz gt 0L)
if index(0) ne -1L then n2o5yz(index)=n2o5yz(index)/float(nn2o5yz(index))
index=where(nn2o5yz eq 0L)
if index(0) ne -1L then n2o5yz(index)=-99.
index=where(n2o5yz ne -99.)
if index(0) ne -1L then print,'N2O5 ',min(n2o5yz(index)),max(n2o5yz(index)),max(nn2o5yz)

index=where(no3yz gt 0L)
if index(0) ne -1L then o3yz(index)=o3yz(index)/float(no3yz(index))
index=where(no3yz eq 0L)
if index(0) ne -1L then o3yz(index)=-99.
index=where(o3yz ne -99.)
if index(0) ne -1L then print,'O3 ',min(o3yz(index)),max(o3yz(index)),max(no3yz)
;
; save time-averaged zonal mean ACE data
;
save,file='ACEupdate_'+sdate0+'-'+sdate1+'.sav',latbin,altitude,tempyz,ntempyz,ch4yz,nch4yz,$
     clono2yz,nclono2yz,h2oyz,nh2oyz,hno3yz,nhno3yz,n2oyz,nn2oyz,no2yz,nno2yz,n2o5yz,nn2o5yz,$
     o3yz,no3yz,cfc11yz,ncfc11yz,cfc12yz,ncfc12yz
;
; check zonal means
;
for ii=0L,10L do begin   ; loop over species
    if ii eq 0L then begin
       plotdata=tempyz & plottitle='Temperature'
    endif
    if ii eq 1L then begin
       plotdata=ch4yz & plottitle='Methane'
    endif
    if ii eq 2L then begin
       plotdata=clono2yz & plottitle='ClONO2'
    endif
    if ii eq 3L then begin
       plotdata=h2oyz & plottitle='H2O'
    endif
    if ii eq 4L then begin
       plotdata=hno3yz & plottitle='HNO3'
    endif
    if ii eq 5L then begin
       plotdata=n2oyz & plottitle='N2O'
    endif
    if ii eq 6L then begin
       plotdata=no2yz & plottitle='NO2'
    endif
    if ii eq 7L then begin
       plotdata=n2o5yz & plottitle='N2O5'
    endif
    if ii eq 8L then begin
       plotdata=o3yz & plottitle='Ozone'
    endif
    if ii eq 9L then begin
       plotdata=cfc11yz & plottitle='CFC11'
    endif
    if ii eq 10L then begin
       plotdata=cfc12yz & plottitle='CFC12'
    endif
;
; fill data voids regions?
;
    plotfilled=plotdata
    for k=0,n_elements(altitude)-1L do begin
        plotlev=reform(plotdata(*,k))
        index1=where(plotlev ne -99.,ngood)
        index2=where(plotlev eq -99.)
        if ngood gt 1 and index1(0) ne -1 and index2(0) ne -1 then begin
           filled=interpol(plotlev(index1),index1,index2)
           plotfilled(index2,k)=filled
        endif
    endfor
;   plotdata=plotfilled
;
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    nlvls=31
    index=where(plotdata ne -99.)
if index(0) eq -1L then goto,jumpspec
    imin=min(plotdata(index))
    imax=max(plotdata(index))
    iint=(imax-imin)/float(nlvls-1L)
    level=imin+iint*findgen(nlvls)
if max(level) eq 0. then goto,jumpspec
    col1=1+indgen(nlvls)*mcolor/nlvls
    contour,plotdata,latbin,altitude,xrange=[-90.,90.],yrange=[1.,120.],xticks=6,$
            charsize=2,xtitle='Latitude',ytitle='Altitude (km)',levels=level,/cell_fill,$
            title=sdate0+'-'+sdate1+' ACE '+plottitle,c_color=col1,color=0,min_value=-99.
    contour,plotdata,latbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
            c_labels=0*level

    xmnb=xmx+.07
    xmxb=xmnb+.01
    set_viewport,xmnb,xmxb,ymn,ymx
    !type=2^2+2^3+2^5
    plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0,charsize=2
    xbox=[0,10,10,0,0]
    y1=imin
    dy=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor

ymnb=ymn-.35
ymxb=ymnb+.25
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3
plot,fdoy,latitude,psym=8,color=0,yrange=[-90.,90.],xrange=[1.,366.],xticks=12,$
     xtickname=month,charsize=2,ytitle='Latitude',yticks=6
loadct,0
plots,iday0,-90
plots,iday0,90.,/continue,thick=8,color=200
plots,iday1,-90
plots,iday1,90.,/continue,thick=8,color=200
plots,iday0,90.
plots,iday1,90.,/continue,thick=8,color=200
plots,iday0,-90.
plots,iday1,-90.,/continue,thick=8,color=200
loadct,38

wait,1
jumpspec:
endfor          ; loop over species
end
