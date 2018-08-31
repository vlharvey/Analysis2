;
; spring time period only
; plot time-altitude section of vortex CO
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
xorig=[0.25,0.25,0.25,0.25]
yorig=[0.8,0.6,0.4,0.2]
xlen=0.5
ylen=0.12
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
   device,/landscape,bits=8,filename='zt_ace_co_spring_3pan.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
erase
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2005','2006','2007','2008']
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
    restore,dira+'dmps_ace_v2.2.meto.'+syear(iyear)
    restore,dira+'cat_ace_v2.2.'+syear(iyear)
    restore,dira+'co_ace_v2.2.'+syear(iyear)
    comix=mix & comask=mask
altitude_save=altitude
;
; take only good values
;
    index=where(comask eq -99.)                         ; bad co data
    comix(index)=-99.
;
; daily average CO in the vortex
;
    nday=long(max(fdoy))-long(min(fdoy))+1L
    nday=365.
    nz=n_elements(altitude)
    oco=fltarr(nday,nz)
    olat=fltarr(nday,nz)
    num=lonarr(nday,nz)
    for iday=1L,nday do begin
        today=where(long(fdoy) eq iday and latitude gt 0.,nprof)	; Arctic
        if nprof le 1L then goto,skipday
        coday=reform(comix(today,*))
        for iprof=0L,nprof-1L do begin
            cos0=reform(coday(iprof,*))
            if latitude(today(iprof)) ge 50. then begin
            for k=0L,nz-1L do begin	; loop from the bottom up
                if cos0(k) ne -99. then begin
                   oco(iday-1L,k)=oco(iday-1L,k)+cos0(k)
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
    oco(index)=1.e6*oco(index)/float(num(index))
    olat(index)=olat(index)/float(num(index))
    oday=0.*olat
    for k=0L,nz-1L do oday(*,k)=1.+findgen(nday)
;
; fill
;
    oco_fill=oco
    for j=0,nz-1 do begin
        dummy=reform(oco(*,j))
        index1=where(dummy gt 0.,ngood)
        index2=where(dummy le 0.,nbad)
        if ngood gt 1L and nbad gt 1L then begin
           filled=interpol(dummy(index1),index1,index2)
           oco_fill(index2,j)=filled
        endif
    endfor
    if syear(iyear) eq '2004' then begin
       oco_fill(0:35,*)=0.      ; do not extrapolate data voids 
       index=where(fdoy lt 30.)
       latitude(index)=0.
    endif
    if syear(iyear) eq '2007' then begin
       if max(fdoy) lt 360. then oco_fill(long(max(fdoy)):long(nday)-1L,*)=0.      ; do not extrapolate data voids
    endif
    oco=oco_fill
;
; smooth
;
    smoothit,oco,ocosmooth
;   oco=ocosmooth
;index=where(olat lt 50. and (oday lt 51. or oday gt 81.))
;oco(index)=-9999.
;
; store save file
;
    save,file='vortex_co_'+syear(iyear)+'.sav',nday,fdoy,altitude,latitude,oco
;
; postscript file
;
    !type=2^2+2^3
    xmn=xorig(iyear)
    xmx=xorig(iyear)+xlen
    ymn=yorig(iyear)
    ymx=yorig(iyear)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    level=[0.001,0.01,0.1,0.5,1.0,2.,3.,5.,7.,10.,15.,20.,25.,30.,40.,50.]
nlvls=n_elements(level)
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
    print,min(fdoy),max(fdoy)
    kday=365.
    leapdy=(long(syear(iyear)) mod 4)
    kday=kday+leapdy
;kday=long(max(fdoy))	; uncomment to "zoom" in on partial year
kday=91
    xlab=[' ',' ',' ',' ']
    contour,oco,1.+findgen(nday),altitude_save,/noeras,c_color=col1,/cell_fill,levels=level,$
         yrange=[30.,90.],xrange=[1.,kday],xticks=n_elements(xlab)-1L,xtickname=xlab,$
         charsize=2,ytitle='Altitude (km)',color=0,min_value=0.
    contour,oco,1.+findgen(nday),altitude_save,/noeras,color=0,/follow,levels=level,/overplot,$
            c_labels=0*level,min_value=0.
    xyouts,5.,35.,syear(iyear),/data,color=0,charsize=3
    set_viewport,xorig(iyear),xorig(iyear)+xlen,ymn-0.05,ymn-0.005
    index=where(fdoy le kday)
    plot,fdoy(index),latitude(index),psym=1,color=0,yrange=[30.,90.],xrange=[1.,kday],xticks=n_elements(xlab)-1L,$
         xtickname=xlab,charsize=2,ytitle='Lat',yticks=2,$
         ytickv=[30.,60.,90.],symsize=0.5
    xyouts,15.,10.,'Jan',/data,color=0,charsize=2,alignment=0.5
    xyouts,46.,10.,'Feb',/data,color=0,charsize=2,alignment=0.5
    xyouts,74.,10.,'Mar',/data,color=0,charsize=2,alignment=0.5
    jumpyear:
endfor  ; loop over years
!type=2^2+2^3+2^6
xmnb=xorig(3)
xmxb=xmnb+xlen
ymnb=yorig(3)-cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmnb,xmxb,ymnb,ymxb
slab=' '+strarr(n_elements(level))
plot,[min(level),max(level)],[0,0],yrange=[0,10],color=0,$
     xticks=n_elements(level)-1L,xtickname=slab,$
     xrange=[min(level),max(level)],charsize=2,xtitle='CO (ppmv)   '
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
    if flab0 lt 1. and flab0 ge 0.01 then slab0=strcompress(string(format='(f4.2)',flab0),/remove_all)
    if flab0 ge 1. then slab0=strcompress(long(slab0),/remove_all)
    x1=x1+dx
    xyouts,x1,-25,slab0,orientation=90,charsize=1.5,/data,color=0,alignment=0.5
endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_ace_co_spring_3pan.ps -rotate -90 zt_ace_co_spring_3pan.jpg'
endif
end
