;
; spring time period only
; plot time-altitude section of vortex NOx 2004-2005 and 2006-2005
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
nxdim=600 & nydim=600
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2006']
;syear=['2006']
nyear=n_elements(syear)
nlvls=21L
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
;
; restore 2005 and save
;
restore,file='vortex_nox_2005.sav'
onox2005=onox
fdoy2005=fdoy
lat2005=latitude
;
; loop over years
;
for iyear=0L,nyear-1L do begin

restore,file='vortex_nox_'+syear(iyear)+'.sav'
onoxdiff=0.*onox2005
for i=0L,364L do begin
    index=where(onox(i,*) ne 0.)
    if index(0) ne -1L then onoxdiff(i,index)=onox(i,index)-onox2005(i,index)
endfor
if syear(iyear) eq '2004' then onoxdiff(0:35,*)=-9.e6      ; do not extrapolate data voids
;
; postscript file
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,font_size=9
       device,/landscape,bits=8,filename='zt_ace_noxdiff_spring_'+syear(iyear)+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    level=-1000.+100.*findgen(21)
;   level=-200.+20.*findgen(21)
    print,min(fdoy),max(fdoy)
    kday=365.
    leapdy=(long(syear(iyear)) mod 4)
    kday=kday+leapdy
;kday=long(max(fdoy))	; uncomment to "zoom" in on partial year
kday=91
    contour,onoxdiff,min(fdoy)+findgen(nday),altitude,/noeras,c_color=col1,/cell_fill,levels=level,$
         yrange=[30.,90.],xrange=[1.,kday],xticks=3,min_value=-9.e6,$
         xtickname=['                            Jan',$
                    '                            Feb',$
                    '                            Mar','  '],$
         charsize=1.75,ytitle='Alitude (km)',color=0,$
         title='ACE Arctic NOx '+syear(iyear)+' - 2005'
    index=where(level gt 0.)
    contour,onoxdiff,min(fdoy)+findgen(nday),altitude,/noeras,color=0,/follow,levels=level(index),/overplot,$
            c_labels=0*level(index),min_value=-9.e6
    index=where(level lt 0.)
    contour,onoxdiff,min(fdoy)+findgen(nday),altitude,/noeras,color=mcolor,/follow,levels=level(index),/overplot,$
            c_labels=0*level(index),min_value=-9.e6
        contour,onoxdiff,min(fdoy)+findgen(nday),altitude,/noeras,color=0,/follow,levels=[0.],/overplot,$
            c_labels=[0],thick=3,min_value=-9.e6
    !type=2^2+2^3+2^5
    xmnb=xorig(0)+xlen+cbaryoff
    xmxb=xmnb+cbarydel
    set_viewport,xmnb,xmxb,ymn,ymx
    slab=strcompress(long(level),/remove_all)
    plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
         yticks=n_elements(level)-1L,ytickname=slab,$
         yrange=[min(level),max(level)],charsize=1.5,title='ppbv'
    xbox=[0,10,10,0,0]
    y1=min(level)
    dy=(max(level)-min(level))/float(nlvls)
    for j=0,nlvls-1 do begin
        ybox=[y1,y1,y1+dy,y1+dy,y1]
        polyfill,xbox,ybox,color=col1(j)
        y1=y1+dy
    endfor
    !type=2^2+2^3
    set_viewport,xorig(0),xorig(0)+xlen,0.1,0.2
    plot,fdoy,latitude,psym=1,color=0,yrange=[30.,90.],xrange=[0.,kday],xticks=3,$
         xtickname=' '+strarr(12),charsize=1.75,ytitle='Latitude',yticks=2,$
         ytickv=[30.,60.,90.],symsize=0.5

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim zt_ace_noxdiff_spring_'+syear(iyear)+'.ps -rotate -90 zt_ace_noxdiff_spring_'+syear(iyear)+'.jpg'
    endif
    jumpyear:
endfor	; loop over years
end
