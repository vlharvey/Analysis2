;
; spring time period only
; plot time-altitude section of vortex NOx
;
@fillit
@smoothit

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=800 & nydim=800
xorig=[0.25,0.25,0.25]
yorig=[0.775,0.5,0.225]
xlen=0.5
ylen=0.2
cbaryoff=0.1
cbarydel=0.01
!NOERAS=-1
;!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
;  !p.font=0
   !p.charthick=3
   device,font_size=8
   device,/landscape,bits=8,filename='zt_ace_nox_spring_3pan+co.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2005','2006']
nyear=n_elements(syear)
;
; loop over years
;
for iyear=0L,nyear-1L do begin
;
; restore ACE SOSST data
;
    ex=findfile(dira+'cat_ace_v2.2.'+syear(iyear))
    if ex(0) eq '' then goto,jumpyear
    restore,dira+'dmps_meto_ace_v2.2.'+syear(iyear)
    dmpdate=date
    restore,dira+'cat_ace_v2.2.'+syear(iyear)
    restore,dira+'no_ace_v2.2.'+syear(iyear)
    nomix=mix & nomask=mask
    restore,dira+'no2_ace_v2.2.'+syear(iyear)
    no2mix=mix & no2mask=mask
help,no2mix,dmpdate
;
; NOx is NO + NO2 and take only good values
;
    noxmix=0.*nomix
    index=where(nomask ne -99. and no2mask ne -99.)	; both good
    noxmix(index)=nomix(index)+no2mix(index)
    index=where(nomask ne -99. and no2mask eq -99.)	; only NO good
    noxmix(index)=nomix(index)
    index=where(nomask eq -99. and no2mask ne -99.)	; only NO2 good
    noxmix(index)=no2mix(index)
;
; daily average NOx in the vortex
;
    ovelat_prof=reform(velat_prof(*,*,0))
    nday=long(max(fdoy))-long(min(fdoy))+1L
    nday=365.
    nz=n_elements(altitude)
    onox=fltarr(nday,nz)
    olat=fltarr(nday,nz)
    velatsave=0.*fltarr(nday)
    elatsave=0.*fltarr(nday)
    elatmin=99.+0.*fltarr(nday)
    elatmax=-99.+0.*fltarr(nday)
    num=lonarr(nday,nz)
    zz=where(altitude ge 40.)
    for iday=1L,nday do begin
        today=where(long(fdoy) eq iday and latitude gt 0.,nprof)	; Arctic
        if nprof le 1L then goto,skipday
        noxday=reform(noxmix(today,*))
        elatday=reform(elat_prof(today,*))
        velatday=reform(ovelat_prof(today,*))
index=where(elatday ne -99.)
;if index(0) ne -1L then print,iday,min(elatday(index)),max(elatday)
        for iprof=0L,nprof-1L do begin
            noxs=reform(noxday(iprof,*))
            elats=reform(elatday(iprof,*))
            velats=reform(velatday(iprof,*))
            index=where(elats ne -99. and velats ne -99.,npts)
;           if index(0) eq -1L then goto,jumpprof

            if latitude(today(iprof)) ge 50. then begin

            for k=nz-1L,0L,-1L do begin	; loop from the bottom up
;               if elats(k) ne -99. and velats(k) ne -99. and $
;                  elats(k) ge velats(k) and noxs(k) ne 0. then begin
                if noxs(k) ne 0. then begin
                   onox(iday-1L,k)=onox(iday-1L,k)+noxs(k)
                   olat(iday-1L,k)=olat(iday-1L,k)+latitude(today(iprof))
                   num(iday-1L,k)=num(iday-1L,k)+1L
                endif
            endfor
index=where(altitude ge 40. and elats gt 0.,nn)
if index(0) ne -1L then begin
   if min(elats(index)) lt elatmin(iday-1L) then elatmin(iday-1L)=min(elats(index))
   if max(elats(index)) gt elatmax(iday-1L) then elatmax(iday-1L)=max(elats(index))
   elatsave(iday-1L)=total(elats(index))/float(nn)
endif
index=where(altitude ge 40. and velats gt 0.,nn)
if index(0) ne -1L then begin
   velatsave(iday-1L)=total(velats(index))/float(nn)
;  print,iday,total(velats(index))/float(nn)
endif
            endif
;
; deal with mesosphere where there is no elat information
; assume the mesospheric profile is inside the vortex if it is inside above 40 km
;
            noxs=reform(noxday(iprof,zz))
            elats=reform(elatday(iprof,zz))
            velats=reform(velatday(iprof,zz))
            index=where(elats ne -99. and velats ne -99.,npts)
;           if index(0) eq -1L then goto,jumpprof
;           elat0=total(elats(index))/float(npts)
;           velat0=total(velats(index))/float(npts)
;           if elat0 ge velat0 then begin 	; in the upper stratospheric vortex 
            if latitude(today(iprof)) ge 50. then begin
               index=where(elats eq -99. or velats ne -99.,npts)   ; loop over unresolved points

               for kk=0L,npts-1L do begin
                   k=index(kk)
                   if altitude(k) ge 40. and velats(k) eq -99. and noxs(k) ne 0. then begin
                      onox(iday-1L,k)=onox(iday-1L,k)+noxs(k)
                      olat(iday-1L,k)=olat(iday-1L,k)+latitude(today(iprof))
                      num(iday-1L,k)=num(iday-1L,k)+1L
                   endif
               endfor
            endif
            jumpprof:
        endfor
        skipday:
    endfor
;
; daily average
;
    index=where(num gt 0L)
    if index(0) eq -1L then goto,jumpyear
    onox(index)=1.e6*onox(index)/float(num(index))
    olat(index)=olat(index)/float(num(index))
    oday=0.*olat
    for k=0L,nz-1L do oday(*,k)=1.+findgen(nday)
;
; fill
;
    onox_fill=onox
    for j=0,nz-1 do begin
        dummy=reform(onox(*,j))
        index1=where(dummy gt 0.,ngood)
        index2=where(dummy le 0.,nbad)
        if ngood gt 1L and nbad gt 1L then begin
           filled=interpol(dummy(index1),index1,index2)
           onox_fill(index2,j)=filled
        endif
    endfor
    if syear(iyear) eq '2004' then begin
       onox_fill(0:35,*)=0.      ; do not extrapolate data voids 
       index=where(fdoy lt 30.)
       latitude(index)=0.
       velatsave(0:35)=0.
       elatsave(0:35)=0.
       elatmin(0:35)=99.
       elatmax(0:35)=-99.
    endif
    onox=onox_fill
;
; smooth
;
    smoothit,onox,onoxsmooth
    onox=onoxsmooth
index=where(olat lt 50. and (oday lt 51. or oday gt 81.))
onox(index)=-9999.
;
; store save file
;
    save,file='vortex_nox_'+syear(iyear)+'.sav',nday,fdoy,altitude,latitude,onox
    restore,'vortex_co_'+syear(iyear)+'.sav'	;,nday,fdoy,altitude,latitude,oco
    smoothit,oco,ocosmooth
    oco=ocosmooth
;
; postscript file
;
    !type=2^2+2^3
    xmn=xorig(iyear)
    xmx=xorig(iyear)+xlen
    ymn=yorig(iyear)
    ymx=yorig(iyear)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    level=[1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,10000.]/1000.
    nlvls=n_elements(level)
    col1=1L+indgen(nlvls)*mcolor/float(nlvls)
    print,min(fdoy),max(fdoy)
    kday=365.
    leapdy=(long(syear(iyear)) mod 4)
    kday=kday+leapdy
;kday=long(max(fdoy))	; uncomment to "zoom" in on partial year
kday=91
    xlab=[' ',' ',' ',' ']
    contour,onox,1.+findgen(nday),altitude,/noeras,c_color=col1,/cell_fill,levels=level,$
         yrange=[30.,90.],xrange=[1.,kday],xticks=n_elements(xlab)-1L,xtickname=xlab,$
         charsize=2,ytitle='Altitude (km)',color=0,min_value=0.
    contour,onox,1.+findgen(nday),altitude,/noeras,color=0,/follow,levels=level,/overplot,$
            c_labels=0*level,min_value=0.

    colevel=[0.5,1.,2.5,5.,7.,10.]
    colevel=[0.5,5.]
    contour,oco,1.+findgen(nday),altitude,/noeras,color=mcolor,/follow,levels=colevel,/overplot,$
            c_labels=1L+0*colevel,min_value=0.,thick=6,c_charthick=4,c_charsize=2

    xyouts,5.,35.,syear(iyear),/data,color=0,charsize=3.2,charthick=8
    set_viewport,xorig(iyear),xorig(iyear)+xlen,ymn-0.05,ymn-0.005
    index=where(fdoy le kday)
    plot,fdoy(index),latitude(index),psym=8,color=0,yrange=[30.,90.],xrange=[1.,kday],xticks=n_elements(xlab)-1L,$
         xtickname=xlab,charsize=1.5,yticks=2,ytickv=[30.,60.,90.],symsize=0.75,ytitle='Lat',xticklen=0.02
;   xyouts,1.,34.,'ACE Latitudes',/data,charsize=1.6,color=0
loadct,0
ee=where(elatmin ne 99.)
;oplot,1.+ee,elatmin(ee),psym=8,color=220
ee=where(elatmax ne -99.)
;oplot,1.+ee,elatmax(ee),psym=8,color=220
ee=where(elatsave gt 0.)
;if ee(0) ne -1L then oplot,1.+ee,elatsave(ee),psym=8,color=175
ee=where(velatsave gt 0.)
;if ee(0) ne -1L then oplot,1.+ee,velatsave(ee),psym=8,color=125
oplot,fdoy(index),latitude(index),psym=8,symsize=0.5,color=0
loadct,38
    xyouts,15.,2.,'Jan',/data,color=0,charsize=2,alignment=0.5
    xyouts,46.,2.,'Feb',/data,color=0,charsize=2,alignment=0.5
    xyouts,74.,2.,'Mar',/data,color=0,charsize=2,alignment=0.5
    jumpyear:
endfor  ; loop over years
!type=2^2+2^3+2^6
xmnb=xorig(2)
xmxb=xmnb+xlen
ymnb=yorig(2)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmnb,xmxb,ymnb,ymxb
slab=' '+strarr(n_elements(level))
plot,[min(level),max(level)],[0,0],yrange=[0,10],color=0,$
     xticks=n_elements(level)-1L,xtickname=slab,$
     xrange=[min(level),max(level)],charsize=2,xtitle='NO!lx!n (ppmv)'
ybox=[0,10,10,0,0]
x1=min(level)
dx=(max(level)-min(level))/float(nlvls)
for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
endfor
slab=strcompress(string(format='(f7.3)',level),/remove_all)
x1=min(level)
for i=0L,n_elements(slab)-1L do begin
    slab0=slab(i)
    flab0=float(slab(i))
    if flab0 lt 0.01 then slab0=strcompress(string(format='(f5.3)',flab0),/remove_all)
    if flab0 ge 0.01 and flab0 lt 0.1 then slab0=strcompress(string(format='(f4.2)',flab0),/remove_all)
    if flab0 ge 0.1 and flab0 lt 1. then slab0=strcompress(string(format='(f3.1)',flab0),/remove_all)
    if flab0 ge 1. then slab0=strcompress(long(slab0),/remove_all)
    x1=x1+dx/2.
    xyouts,x1,-26.,slab0,orientation=90,charsize=1.5,/data,color=0,alignment=0.5
    x1=x1+dx/2.
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_ace_nox_spring_3pan+co.ps -rotate -90 zt_ace_nox_spring_3pan+co.jpg'
endif
end
