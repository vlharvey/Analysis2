;
; plot occultation locations from SOSST database
;
a=findgen(8)*(2*!pi/8.)
usersym,0.3*cos(a),0.3*sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=800 & nydim=800
xorig=[0.15]
yorig=[0.25]
cbaryoff=0.02
cbarydel=0.02
!NOERAS=-1
!p.font=0
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='sage2_loc.ps'
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dirp='/atmos/aura3/data/POAM_data/Datfiles_SOSST/'
dirh='/atmos/aura3/data/HALOE_data/Datfiles_SOSST/'
dirs2='/atmos/aura3/data/SAGE_II_data/Datfiles_SOSST/'
dirs3='/atmos/aura3/data/SAGE_III_data/Datfiles_SOSST/'
diri='/atmos/aura3/data/ILAS_data/Datfiles_SOSST/'
dira='/atmos/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['1984','1985','1986','1987','1989','1990','1991','1992','1993',$
       '1994','1995','1996','1997','1999','2000','2001','2002','2003',$
       '2004','2005']	;,'2006']
nyear=n_elements(syear)
col1=indgen(nyear)*mcolor/(nyear-1L)
col1(0)=1L
col1(nyear-1)=mcolor-5L
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+0.8
ymn=yorig(0)
ymx=yorig(0)+0.6
set_viewport,xmn,xmx,ymn,ymx
plot,findgen(366),-90.+findgen(181),/nodata,/noeras,color=0,charsize=1.5,charthick=2,$
     yrange=[-90.,90.],xrange=[1.,365.],yticks=6,xticks=11,xtickv=15+30*findgen(12),$
     xtickname=month,ytitle='Latitude',$
     title='SAGE II'
;    title='POAM II&III SAGE II&III ILAS I&II HALOE ACE'

loadct,23
xinc=0.8/float(nyear/2L)
for iyear=0L,nyear-1L do begin
    if iyear le 9L then xyouts,xmn+iyear*xinc,ymn-0.1,syear(iyear),color=col1(iyear),/normal,charsize=1.5
    if iyear gt 9L then xyouts,xmn+(iyear-10L)*xinc,ymn-0.125,syear(iyear),color=col1(iyear),/normal,charsize=1.5
;
; SAGE II
;
    ex=findfile(dirs2+'cat_sage2_v7.0.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirs2+'cat_sage2_v7.0.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; SAGE III
;
    ex=findfile(dirs3+'cat_sage3_v3.00.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirs3+'cat_sage3_v3.00.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; POAM II
;
    ex=findfile(dirp+'cat_poam2_v6.0.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirp+'cat_poam2_v6.0.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; POAM III
;
    ex=findfile(dirp+'cat_poam3_v4.0.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirp+'cat_poam3_v4.0.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; HALOE
;
    ex=findfile(dirh+'cat_haloe_v19.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirh+'cat_haloe_v19.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; ACE
;
    ex=findfile(dira+'cat_ace_v2.2.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dira+'cat_ace_v2.2.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; ILAS I
;
    ex=findfile(diri+'cat_ilas_v06.10.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,diri+'cat_ilas_v06.10.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; ILAS II
;
    ex=findfile(diri+'cat_ilas2_v1.4.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,diri+'cat_ilas2_v1.4.'+syear(iyear)
;      oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif

endfor
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim sage2_loc.ps -rotate -90 sage2_loc.jpg'
endif
end
