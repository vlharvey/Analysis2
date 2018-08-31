;
; plot occultation locations from SOSST database
;
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=600 & nydim=600
xorig=[0.15]
yorig=[0.25]
cbaryoff=0.02
cbarydel=0.02
!NOERAS=-1
!p.font=1
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
   device,/landscape,bits=8,filename='sosst_loc+mark_hirdls_days.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
month=['             Jan',$
       '             Feb',$
       '             Mar',$
       '             Apr',$
       '             May',$
       '             Jun',' ']
;       'Jul','Aug','Sep','Oct','Nov','Dec']
dirp='/aura3/data/POAM_data/Datfiles_SOSST/'
dirh='/aura3/data/HALOE_data/Datfiles_SOSST/'
dirs2='/aura3/data/SAGE_II_data/Datfiles_SOSST/'
dirs3='/aura3/data/SAGE_III_data/Datfiles_SOSST/'
diri='/aura3/data/ILAS_data/Datfiles_SOSST/'
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['1984','1985','1986','1987','1989','1990','1991','1992','1993',$
       '1994','1995','1996','1997','1999','2000','2001','2002','2003',$
       '2004','2005','2006']
syear=['2005']
nyear=n_elements(syear)
col1=indgen(nyear)*mcolor/(nyear-1L)
col1(0)=1L
col1(nyear-1)=mcolor-5L
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+0.7
ymn=yorig(0)
ymx=yorig(0)+0.5
set_viewport,xmn,xmx,ymn,ymx
plot,1.+findgen(365),-90.+findgen(181),/nodata,/noeras,color=0,$
     yrange=[-90.,90.],xrange=[1.,181.],yticks=6,xticks=6,$
     xtickname=month,charsize=1.75,ytitle='Latitude',$
     title=syear(0)
xyouts,xmx+0.01,ymx-0.15,'POAM III',color=130,/normal,charsize=2
xyouts,xmx+0.01,ymx-0.20,'SAGE II&III',color=70,/normal,charsize=2
xyouts,xmx+0.01,ymx-0.25,'HALOE',color=250,/normal,charsize=2
xyouts,xmx+0.01,ymx-0.30,'ACE',color=210,/normal,charsize=2

loadct,0
plots,21.,-90.,/data				; Jan 21
plots,21.,90.,/continue,/data,thick=2,color=200
xyouts,21.,-97.,'21',color=0,/data,charsize=1.25,alignment=0.5
plots,28.,-90.,/data				; Jan 28
plots,28.,90.,/continue,/data,thick=2,color=200
xyouts,28.,-97.,'28',color=0,/data,charsize=1.25,alignment=0.5
plots,38.,-90.,/data				; Feb 7
plots,38.,90.,/continue,/data,thick=2,color=200
xyouts,38.,-97.,'7',color=0,/data,charsize=1.25,alignment=0.5
plots,40.,-90.,/data				; Feb 9
plots,40.,90.,/continue,/data,thick=2,color=200
plots,48.,-90.,/data                            ; Feb 17
plots,48.,90.,/continue,/data,thick=2,color=200
xyouts,48.,-97.,'17',color=0,/data,charsize=1.25,alignment=0.5
plots,49.,-90.,/data                            ; Feb 18
plots,49.,90.,/continue,/data,thick=2,color=200
plots,50.,-90.,/data                            ; Feb 19
plots,50.,90.,/continue,/data,thick=2,color=200
plots,54.,-90.,/data                            ; Feb 23
plots,54.,90.,/continue,/data,thick=2,color=200
xyouts,54.,-97.,'23',color=0,/data,charsize=1.25,alignment=0.5
plots,60.,-90.,/data                            ; Mar 1
plots,60.,90.,/continue,/data,thick=2,color=200
xyouts,60.,-97.,'1',color=0,/data,charsize=1.25,alignment=0.5
plots,70.,-90.,/data                            ; Mar 11
plots,70.,90.,/continue,/data,thick=2,color=200
xyouts,70.,-97.,'11',color=0,/data,charsize=1.25,alignment=0.5
plots,81.,-90.,/data                            ; Mar 22
plots,81.,90.,/continue,/data,thick=2,color=200
xyouts,81.,-97.,'22',color=0,/data,charsize=1.25,alignment=0.5
plots,95.,-90.,/data                            ; Apr 5
plots,95.,90.,/continue,/data,thick=2,color=200
xyouts,95.,-97.,'5',color=0,/data,charsize=1.25,alignment=0.5
plots,110.,-90.,/data                            ; Apr 20
plots,110.,90.,/continue,/data,thick=2,color=200
xyouts,110.,-97.,'20',color=0,/data,charsize=1.25,alignment=0.5
plots,125.,-90.,/data                            ; May 5
plots,125.,90.,/continue,/data,thick=2,color=200
xyouts,125.,-97.,'5',color=0,/data,charsize=1.25,alignment=0.5
plots,126.,-90.,/data                            ; May 6
plots,126.,90.,/continue,/data,thick=2,color=200
plots,131.,-90.,/data                            ; May 11
plots,131.,90.,/continue,/data,thick=2,color=200
xyouts,131.,-97.,'11',color=0,/data,charsize=1.25,alignment=0.5
plots,162.,-90.,/data                            ; Jun 11
plots,162.,90.,/continue,/data,thick=2,color=200
xyouts,162.,-97.,'11',color=0,/data,charsize=1.25,alignment=0.5
plots,168.,-90.,/data                            ; Jun 17
plots,168.,90.,/continue,/data,thick=2,color=200
xyouts,168.,-97.,'17',color=0,/data,charsize=1.25,alignment=0.5
plots,170.,-90.,/data                            ; Jun 19
plots,170.,90.,/continue,/data,thick=2,color=200
plots,172.,-90.,/data                            ; Jun 21
plots,172.,90.,/continue,/data,thick=2,color=200
loadct,38

for iyear=0L,nyear-1L do begin
;
; read zonal mean marker
;
    close,10
    openr,10,'/aura2/harvey/Vortex_Ozone/Datfiles/markbar_'+syear(iyear)+'_v2.dat'
    nday=0L & nr=0L & nth=0L
    readu,10,nday,nr,nth
    days=fltarr(nday)
    alat=fltarr(nr)
    thlev=fltarr(nth)
    readu,10,days,alat,thlev
    markbar=fltarr(nday,nr,nth)
    marklow=fltarr(nday,nr,nth)
    markhigh=fltarr(nday,nr,nth)
    readu,10,markbar,marklow,markhigh
    close,10
    loadct,0
    for k=0,n_elements(thlev)-1L do begin
        if thlev(k) ge 600. then begin
           index=where(thlev eq thlev(k))
           markbarlev=reform(markbar(*,*,index(0)))
           contour,markbarlev,days,alat,/overplot,/cell_fill,level=[1.0],color=100,thick=2

           markbarlev=reform(markhigh(*,*,index(0)))
           index=where(markbarlev gt -1.0)
           markbarlev(index)=-9999.
           contour,markbarlev,days,alat,/overplot,/cell_fill,level=[-1],color=170,thick=2
        endif
    endfor
    loadct,38
;
; SAGE II
;
    ex=findfile(dirs2+'cat_sage2_v6.2.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirs2+'cat_sage2_v6.2.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=70
    endif
;
; SAGE III
;
    ex=findfile(dirs3+'cat_sage3_v3.00.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirs3+'cat_sage3_v3.00.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=70
    endif
;
; POAM II
;
    ex=findfile(dirp+'cat_poam2_v6.0.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirp+'cat_poam2_v6.0.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=130
    endif
;
; POAM III
;
    ex=findfile(dirp+'cat_poam3_v4.0.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirp+'cat_poam3_v4.0.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=130
    endif
;
; HALOE
;
    ex=findfile(dirh+'cat_haloe_v19.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dirh+'cat_haloe_v19.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=250
    endif
;
; ACE
;
    ex=findfile(dira+'cat_ace_v2.2.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,dira+'cat_ace_v2.2.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=210
    endif
;
; ILAS I
;
    ex=findfile(diri+'cat_ilas_v06.10.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,diri+'cat_ilas_v06.10.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif
;
; ILAS II
;
    ex=findfile(diri+'cat_ilas2_v1.4.'+syear(iyear))
    if ex(0) ne '' then begin
       restore,diri+'cat_ilas2_v1.4.'+syear(iyear)
       oplot,fdoy,latitude,psym=8,color=col1(iyear)
    endif

loadct,0
plots,21.,-90.,/data                            ; Jan 21
plots,21.,90.,/continue,/data,thick=2,color=0
plots,28.,-90.,/data                            ; Jan 28
plots,28.,90.,/continue,/data,thick=2,color=0
plots,38.,-90.,/data                            ; Feb 7
plots,38.,90.,/continue,/data,thick=2,color=0
plots,40.,-90.,/data                            ; Feb 9
plots,40.,90.,/continue,/data,thick=2,color=0
plots,48.,-90.,/data                            ; Feb 17
plots,48.,90.,/continue,/data,thick=2,color=0
plots,49.,-90.,/data                            ; Feb 18
plots,49.,90.,/continue,/data,thick=2,color=0
plots,50.,-90.,/data                            ; Feb 19
plots,50.,90.,/continue,/data,thick=2,color=0
plots,54.,-90.,/data                            ; Feb 23
plots,54.,90.,/continue,/data,thick=2,color=0
plots,60.,-90.,/data                            ; Mar 1
plots,60.,90.,/continue,/data,thick=2,color=0
plots,70.,-90.,/data                            ; Mar 11
plots,70.,90.,/continue,/data,thick=2,color=0
plots,81.,-90.,/data                            ; Mar 22
plots,81.,90.,/continue,/data,thick=2,color=0
plots,95.,-90.,/data                            ; Apr 5
plots,95.,90.,/continue,/data,thick=2,color=0
plots,110.,-90.,/data                            ; Apr 20
plots,110.,90.,/continue,/data,thick=2,color=0
plots,125.,-90.,/data                            ; May 5
plots,125.,90.,/continue,/data,thick=2,color=0
plots,126.,-90.,/data                            ; May 6
plots,126.,90.,/continue,/data,thick=2,color=0
plots,131.,-90.,/data                            ; May 11
plots,131.,90.,/continue,/data,thick=2,color=0
plots,162.,-90.,/data                            ; Jun 11
plots,162.,90.,/continue,/data,thick=2,color=0
plots,168.,-90.,/data                            ; Jun 17
plots,168.,90.,/continue,/data,thick=2,color=0
plots,170.,-90.,/data                            ; Jun 19
plots,170.,90.,/continue,/data,thick=2,color=0
plots,172.,-90.,/data                            ; Jun 21
plots,172.,90.,/continue,/data,thick=2,color=0
loadct,38

endfor
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim sosst_loc+mark_hirdls_days.ps -rotate -90 sosst_loc+mark_hirdls_days.jpg'
endif
end
