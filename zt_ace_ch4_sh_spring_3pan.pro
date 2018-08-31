;
; SH version
; spring time period only
; plot time-altitude section of vortex methane
;
@fillit
@smoothit

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,0.3*cos(a),0.3*sin(a),/fill
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=800 & nydim=1100
xorig=[0.25,0.25,0.25]
yorig=[0.775,0.5,0.225]
xlen=0.5
ylen=0.2
cbaryoff=0.08
cbarydel=0.01
!NOERAS=-1
!p.font=1
!x.ticklen=-0.05
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,font_size=9
   device,/landscape,bits=8,filename='zt_ace_ch4_sh_spring_3pan.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2005','2006']
nyear=n_elements(syear)
nlvls=21L
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
;
; loop over years
;
for iyear=0L,nyear-1L do begin
;
; restore ACE SOSST data
;
    ex=findfile(dira+'cat_ace_v2.2.'+syear(iyear))
    if ex(0) eq '' then goto,jumpyear
    restore,dira+'dmps_ace_v2.2.meto.'+syear(iyear)
    restore,dira+'cat_ace_v2.2.'+syear(iyear)
    restore,dira+'ch4_ace_v2.2.'+syear(iyear)
    nomix=mix & nomask=mask
;
; CH4
;
    noxmix=nomix
    index=where(nomask eq -99.)
    noxmix(index)=-99.
;
; daily average CH4 in the vortex
;
    ovelat_prof=reform(velat_prof(*,*,0))
    nday=long(max(fdoy))-long(min(fdoy))+1L
    nday=365.
    nz=n_elements(altitude)
    onox=fltarr(nday,nz)
    olat=fltarr(nday,nz)
    num=lonarr(nday,nz)
    zz=where(altitude ge 40.)
    for iday=1L,nday do begin
        today=where(long(fdoy) eq iday and latitude lt 0.,nprof)	; Antarctic
        if nprof le 1L then goto,skipday
        noxday=reform(noxmix(today,*))
        elatday=reform(elat_prof(today,*))
        velatday=reform(ovelat_prof(today,*))
        for iprof=0L,nprof-1L do begin
            noxs=reform(noxday(iprof,*))
            elats=reform(elatday(iprof,*))
            velats=reform(velatday(iprof,*))
            index=where(elats ne -99. and velats ne -99.,npts)
;           if index(0) eq -1L then goto,jumpprof

            if latitude(today(iprof)) le -50. then begin
            for k=nz-1L,0L,-1L do begin	; loop from the bottom up
                if noxs(k) ne 0. then begin
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
    endif
;   onox=onox_fill
;
; smooth
;
    smoothit,onox,onoxsmooth
;   onox=onoxsmooth
;index=where(olat gt -50. and (oday lt 51. or oday gt 81.))
;onox(index)=-9999.
;
; store save file
;
    save,file='vortex_ch4_sh_'+syear(iyear)+'.sav',nday,fdoy,altitude,latitude,onox
;
; postscript file
;
    !type=2^2+2^3
    xmn=xorig(iyear)
    xmx=xorig(iyear)+xlen
    ymn=yorig(iyear)
    ymx=yorig(iyear)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    level=[0.001,0.01,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0]
nlvls=n_elements(level)
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
    print,min(fdoy),max(fdoy)
    kday=365.
    leapdy=(long(syear(iyear)) mod 4)
    kday=kday+leapdy
;kday=long(max(fdoy))	; uncomment to "zoom" in on partial year
kday0=121
kday1=274
    xlab=[' ',' ',' ',' ',' ',' ']
    contour,onox,1.+findgen(nday),altitude,/noeras,c_color=col1,/cell_fill,levels=level,$
         yrange=[30.,90.],xrange=[kday0,kday1],xticks=n_elements(xlab)-1L,xtickname=xlab,$
         charsize=2,ytitle='Altitude (km)',color=0,min_value=0.
    contour,onox,1.+findgen(nday),altitude,/noeras,color=0,/follow,levels=level,/overplot,$
            c_labels=0*level,min_value=0.
    contour,onox,1.+findgen(nday),altitude,/noeras,color=mcolor,/follow,levels=[0.05,0.1],/overplot,$
            c_labels=[1,1],Min_value=0.,thick=2
    xyouts,kday0+5.,35.,syear(iyear),/data,color=0,charsize=3
    set_viewport,xorig(iyear),xorig(iyear)+xlen,ymn-0.05,ymn-0.005
    index=where(fdoy ge kday0 and fdoy le kday1)
    plot,fdoy(index),latitude(index),psym=1,color=0,yrange=[-90.,-30.],xrange=[kday0,kday1],xticks=n_elements(xlab)-1L,$
         xtickname=xlab,charsize=2,ytitle='Lat',yticks=2,$
         ytickv=[-30.,-60.,-90.],symsize=0.5
    xyouts,135.,-110.,'May',/data,color=0,charsize=2,alignment=0.5
    xyouts,166.,-110.,'Jun',/data,color=0,charsize=2,alignment=0.5
    xyouts,196.,-110.,'Jul',/data,color=0,charsize=2,alignment=0.5
    xyouts,227.,-110.,'Aug',/data,color=0,charsize=2,alignment=0.5
    xyouts,258.,-110.,'Sep',/data,color=0,charsize=2,alignment=0.5
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
     xrange=[min(level),max(level)],charsize=2,xtitle='CH!l4!n (ppmv)   '
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
    if flab0 lt 0.02 then slab0=strcompress(string(format='(f5.3)',flab0),/remove_all)
    if flab0 lt 1. and flab0 ge 0.02 then slab0=strcompress(string(format='(f4.2)',flab0),/remove_all)
    if flab0 ge 1. then slab0=strcompress(long(slab0),/remove_all)
    x1=x1+dx
    xyouts,x1,-25,slab0,orientation=90,charsize=1.5,/data,color=0,alignment=0.5
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_ace_ch4_sh_spring_3pan.ps -rotate -90 zt_ace_ch4_sh_spring_3pan.jpg'
endif
end
