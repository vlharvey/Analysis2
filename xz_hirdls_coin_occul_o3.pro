;
; plot xz sections of HIRDLS at HALOE/SAGE/POAM latitudes
;
@aura2date
@rd_ukmo_nc3
@rd_sage3_o3_soundings
@rd_haloe_o3_soundings
@rd_poam3_o3_soundings
@rd_sage2_o3_soundings

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
nlvls=25
col1=1+indgen(nlvls)*mcolor/nlvls
device,decompose=0
icmm1=icolmax-1
icmm2=icolmax-2
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.15,0.15,0.15,0.15,0.15]
yorig=[0.80,0.65,0.50,0.35,0.20]
xlen=0.3
ylen=0.08
cbaryoff=0.08
cbarydel=0.02
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
dirh='/aura3/data/HALOE_data/Sound_data/haloe_'
dirs='/aura3/data/SAGE_II_data/Sound_data/sage2_'
dirs3='/aura3/data/SAGE_III_data/Sound_data/sage3_solar_'
dirp2='/aura3/data/POAM_data/Sound_data/poam2_'
dirp3='/aura3/data/POAM_data/Sound_data/poam3_'
SpeciesNames = ['Temperature',$
                'H2O', $
                'O3',  $
                'N2O', $
                'HNO3']
GeoLoc = ['Pressure',$
          'Time',$
          'Latitude',$
          'Longitude',$
          'SolarZenithAngle',$
          'LocalSolarTime']
hdir='/aura3/data/HIRDLS_data/Datfiles/'
Hfile=hdir+'HIRDLS2_2000d276_MZ3_c1.he5'

;; load HIRDLS data all at once
hirdls=LoadAuraData(Hfile, [GeoLoc, SpeciesNames])

;; file header and tail for MLS
Mfileh=hdir+'MLS-Aura_L2GP-'
Mfilet='_sAura2c--t_2000d276.he5'

;; loop over species
FOR is = 0,N_ELEMENTS(SpeciesNames)-1 DO BEGIN
    SpeciesName = SpeciesNames(is)
    Mfile=Mfileh + SpeciesName + Mfilet
    IF is EQ 0 THEN mls = LoadAuraData(Mfile, GeoLoc)
    mls=LoadAuraData(Mfile, SpeciesName, mls)
ENDFOR
;
; extract mls and hirdls variables
; time is elapsed seconds since midnight 1 Jan 1993
;
mpress=mls.p		; P               FLOAT     Array[37]
mlev=n_elements(mpress)
mtime=mls.time		; TIME            DOUBLE    Array[3495]
mlat=mls.lat		; LAT             FLOAT     Array[3495]
mlon=mls.lon		; LON             FLOAT     Array[3495]
msza=mls.sza		; SZA             FLOAT     Array[3495]
mlst=mls.lst		; LST             FLOAT     Array[3495]
mprof=n_elements(mlst)
mtemp=mls.t		; T               FLOAT     Array[37, 3495]
mh2o=mls.h2o		; H2O             FLOAT     Array[37, 3495]
mo3=mls.o3		; O3              FLOAT     Array[37, 3495]
mn2o=mls.n2o		; N2O             FLOAT     Array[37, 3495]
mhno3=mls.hno3		; HNO3            FLOAT     Array[37, 3495]

hpress=hirdls.p		;   P               FLOAT     Array[145]
hlev=n_elements(hpress)
htime=hirdls.time	;   TIME            DOUBLE    Array[7848]
hlat=hirdls.lat		;   LAT             FLOAT     Array[7848]
hlon=hirdls.lon		;   LON             FLOAT     Array[7848]
hsza=hirdls.sza		;   SZA             FLOAT     Array[7848]
hlst=hirdls.lst		;   LST             FLOAT     Array[7848]
hprof=n_elements(hlst)
htemp=hirdls.t		;   T               FLOAT     Array[145, 7848]
hh2o=hirdls.h2o		;   H2O             FLOAT     Array[145, 7848]
ho3=hirdls.o3		;   O3              FLOAT     Array[145, 7848]
hn2o=hirdls.n2o		;   N2O             FLOAT     Array[145, 7848]
hno3=hirdls.hno3	;   HNO3            FLOAT     Array[145, 7848]
;
; convert elapsed seconds to dates (yyyymmddhh)
;
aura2date,mdate,mtime
aura2date,hdate,htime
;
; make press,lat,lon 2d
;
mpress2=0.*mo3
mlat2=0.*mo3
mlon2=0.*mo3
for i=0L,mprof-1L do mpress2(*,i)=mpress
for i=0L,mlev-1L do begin
    mlat2(i,*)=mlat
    mlon2(i,*)=mlon
endfor
hpress2=0.*ho3
hlat2=0.*ho3
hlon2=0.*ho3
for i=0L,hprof-1L do hpress2(*,i)=hpress
for i=0L,hlev-1L do begin
    hlat2(i,*)=hlat
    hlon2(i,*)=hlon
endfor
htheta=htemp*(1000./hpress2)^0.286
;
; read UKMO data
;
sdate=strcompress(string(mdate(mprof-1)),/remove_all)
sdate=strmid(sdate,0,8)
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='xz_hirdls_sage_haloe_poam_o3_'+sdate+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
syr=strmid(sdate,0,4)
uyr=strmid(syr,2,2)
smn=strmid(sdate,4,2)
imn=fix(smn)
sdy=strmid(sdate,6,2)
ifile=mon(imn-1)+sdy+'_'+uyr
print,ifile
rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
            pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
      x2d=fltarr(nc,nth)
      y2d=fltarr(nc,nth)
      for k=0,nth-1 do x2d(*,k)=alon
      for j=0,nc-1 do y2d(j,*)=th
;
; read satellite ozone soundings
;
sfile=mon(imn-1)+sdy+'_'+syr
rd_sage3_o3_soundings,dirs3+sfile+'_o3.sound',norbits3,tsage3,$
   xsage3,ysage3,tropps3,tropzs3,tropths3,modes3,o3sage3,psage3,$
   thsage3,zsage3,clsage3,qo3sage3,nlevs3
rd_sage2_o3_soundings,dirs+sfile+'_o3.sound',norbits2,tsage2,$
   xsage2,ysage2,tropps2,tropzs2,tropths2,modes2,o3sage2,psage2,$
   thsage2,zsage2,clsage2,qo3sage2,nlevs2
rd_poam3_o3_soundings,dirp3+sfile+'_o3.sound',norbitp3,tpoam3,$
   xpoam3,ypoam3,troppp3,tropzp3,tropthp3,modep3,o3poam3,ppoam3,$
   thpoam3,zpoam3,clpoam3,qo3poam3,nlevp3
rd_haloe_o3_soundings,dirh+sfile+'_o3.sound',norbith,thal,$
   xhal,yhal,tropph,tropzh,tropthh,modeh,o3hal,phal,$
   thhal,zhal,clhal,qo3hal,nlevh
;
; remove missing and bad data
;
      yhal0=-999. & yhal1=-999.
      if norbith gt 0L then begin
         index=where(yhal ge -90. and yhal le 90.,norbith)
         if index(0) ne -1 then begin
            yhal=yhal(index)
            xhal=xhal(index)
            o3hal=o3hal(index,*)
            thhal=thhal(index,*)
         endif
         good=0.*fltarr(nlevh)
         for ilev=0,nlevh-1 do $
             if max(thhal(*,ilev)) lt 1.00000e+24 then good(ilev)=1.
         index=where(good eq 1.,nlevh)
         thhal=thhal(*,index)
         o3hal=o3hal(*,index)
         qo3hal=qo3hal(*,index)
         index=where(qo3hal/o3hal gt 0.5)
         if index(0) ne -1 then o3hal(index)=-999.
         index=where(modeh eq 0L)
         if index(0) ne -1 then begin
            yhalsr=yhal(index)
            xhalsr=xhal(index)
            o3halsr=o3hal(index,*)
            thhalsr=thhal(index,*)
            index=where(yhalsr ge -90. and yhalsr le 90.)
            if index(0) eq -1 then goto,jumphalsr
            yhal0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yhalsr(index))
               yhal0=result(0)
            endif
         endif
         jumphalsr:
         index=where(modeh eq 1L)
         if index(0) ne -1 then begin
            yhalss=yhal(index)
            xhalss=xhal(index)
            o3halss=o3hal(index,*)
            thhalss=thhal(index,*)
            index=where(yhalss ge -90. and yhalss le 90.)
            if index(0) eq -1 then goto,jumphalss
            yhal1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yhalss(index))
               yhal1=result(0)
            endif
         endif
         jumphalss:
      endif
      ysage30=-999. & ysage31=-999.
      if norbits3 gt 0L then begin
         index=where(qo3sage3/o3sage3 gt 0.5)
         if index(0) ne -1 then o3sage3(index)=-999.
         index=where(modes3 eq 0L)
         if index(0) ne -1 then begin
            ysage3sr=ysage3(index)
            xsage3sr=xsage3(index)
            o3sage3sr=o3sage3(index,*)
            thsage3sr=thsage3(index,*)
            index=where(ysage3sr ge -90. and ysage3sr le 90.)
            if index(0) eq -1 then goto,jumpsage3sr
            ysage30=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage3sr(index))
               ysage30=result(0)
            endif
         endif
         jumpsage3sr:
         index=where(modes3 eq 1L)
         if index(0) ne -1 then begin
            ysage3ss=ysage3(index)
            xsage3ss=xsage3(index)
            o3sage3ss=o3sage3(index,*)
            thsage3ss=thsage3(index,*)
            index=where(ysage3ss ge -90. and ysage3ss le 90.)
            if index(0) eq -1 then goto,jumpsage3ss
            ysage31=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage3ss(index))
               ysage31=result(0)
            endif
         endif
         jumpsage3ss:
      endif
      ysage20=-999. & ysage21=-999.
      if norbits2 gt 0L then begin
         index=where(qo3sage2/o3sage2 gt 0.5)
         if index(0) ne -1 then o3sage2(index)=-999.
         index=where(modes2 eq 0L)
         if index(0) ne -1 then begin
            ysage2sr=ysage2(index)
            xsage2sr=xsage2(index)
            o3sage2sr=o3sage2(index,*)
            thsage2sr=thsage2(index,*)
            index=where(ysage2sr ge -90. and ysage2sr le 90.)
            if index(0) eq -1 then goto,jumpsage2sr
            ysage20=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2sr(index))
               ysage20=result(0)
            endif
         endif
         jumpsage2sr:
         index=where(modes2 eq 1L)
         if index(0) ne -1 then begin
            ysage2ss=ysage2(index)
            xsage2ss=xsage2(index)
            o3sage2ss=o3sage2(index,*)
            thsage2ss=thsage2(index,*)
            index=where(ysage2ss ge -90. and ysage2ss le 90.)
            if index(0) eq -1 then goto,jumpsage2ss
            ysage21=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2ss(index))
               ysage21=result(0)
            endif
         endif
         jumpsage2ss:
      endif
      ypoam30=-999. & ypoam31=-999.
      if norbitp3 gt 0L then begin
         index=where(qo3poam3/o3poam3 gt 0.5)
         if index(0) ne -1 then o3poam3(index)=-999.
         index=where(modep3 eq 0L)
         if index(0) ne -1 then begin
            ypoam3sr=ypoam3(index)
            xpoam3sr=xpoam3(index)
            o3poam3sr=o3poam3(index,*)
            thpoam3sr=thpoam3(index,*)
            index=where(ypoam3sr ge -90. and ypoam3sr le 90.)
            if index(0) eq -1 then goto,jumppoam3sr
            ypoam30=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoam3sr(index))
               ypoam30=result(0)
            endif
         endif
         jumppoam3sr:
         index=where(modep3 eq 1L)
         if index(0) ne -1 then begin
            ypoam3ss=ypoam3(index)
            xpoam3ss=xpoam3(index)
            o3poam3ss=o3poam3(index,*)
            thpoam3ss=thpoam3(index,*)
            index=where(ypoam3ss ge -90. and ypoam3ss le 90.)
            if index(0) eq -1 then goto,jumppoam3ss
            ypoam31=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoam3ss(index))
               ypoam31=result(0)
            endif
         endif
         jumppoam3ss:
      endif
      yorder=[yhal0,yhal1,ypoam30,ypoam31,ysage30,ysage31,ysage20,ysage21]
      instorder=['HALOESR','HALOESS','POAM3SR','POAM3SS',$
                 'SAGE3SR','SAGE3SS','SAGE2SR','SAGE2SS']
;
; remove missing modes
;
      index=where(yorder gt -90.,npan)
      if index(0) eq -1 then goto,jump
      if index(0) ne -1 then begin
         yorder=yorder(index)
         instorder=instorder(index)
      endif
      index=reverse(sort(yorder))
      yorder=yorder(index)
      instorder=instorder(index)
      erase
      xyouts,.37,.95,sdate+' Ozone',charsize=2,/normal
      for ipan=0,npan-1 do begin
;
; plot longitude-altitude sections in correct order S-N
;
          case 1 of
            (instorder(ipan) eq 'HALOESR'): begin
             ydata=yhalsr
             xdata=xhalsr
             o3data=o3halsr
             thdata=thhalsr
            end
            (instorder(ipan) eq 'HALOESS'): begin
             ydata=yhalss
             xdata=xhalss
             o3data=o3halss
             thdata=thhalss
            end
            (instorder(ipan) eq 'POAM3SR'): begin
             ydata=ypoam3sr
             xdata=xpoam3sr
             o3data=o3poam3sr
             thdata=thpoam3sr
            end
            (instorder(ipan) eq 'POAM3SS'): begin
             ydata=ypoam3ss
             xdata=xpoam3ss
             o3data=o3poam3ss
             thdata=thpoam3ss
            end
            (instorder(ipan) eq 'SAGE3SR'): begin
             ydata=ysage3sr
             xdata=xsage3sr
             o3data=o3sage3sr
             thdata=thsage3sr
            end
            (instorder(ipan) eq 'SAGE3SS'): begin
             ydata=ysage3ss
             xdata=xsage3ss
             o3data=o3sage3ss
             thdata=thsage3ss
            end
            (instorder(ipan) eq 'SAGE2SR'): begin
             ydata=ysage2sr
             xdata=xsage2sr
             o3data=o3sage2sr
             thdata=thsage2sr
            end
            (instorder(ipan) eq 'SAGE2SS'): begin
             ydata=ysage2ss
             xdata=xsage2ss
             o3data=o3sage2ss
             thdata=thsage2ss
            end
            else: begin
            goto,noprof
            end
          endcase
          ydata=reform(ydata)
          xdata=reform(xdata)
          o3data=reform(o3data)
          thdata=reform(thdata)
;
; remove bad xdata and sort in longitude
;
          index=where(xdata ge 0. and xdata le 360.)
          if index(0) ne -1 then begin
             ydata=reform(ydata(index))
             xdata=reform(xdata(index))
             o3data=reform(o3data(index,*))
             thdata=reform(thdata(index,*))
          endif
          xsave=xdata
          xdata=0.*thdata
          result=size(thdata)
          ndim=result(0)
          nl=result(2)
          for i=0,nl-1 do begin
              sindex=sort(xsave)
              xdata=xsave(sindex)
              ydata=ydata(sindex)
              o3data(*,i)=o3data(sindex,i)
              thdata(*,i)=thdata(sindex,i)
          endfor
          y0=strcompress(string(FORMAT='(F5.1)',min(ydata)))
          y1=strcompress(string(FORMAT='(F5.1)',max(ydata)))
          level=0.5*findgen(nlvls)
          !type=2^2+2^3
          xmn=xorig(ipan)
          xmx=xorig(ipan)+xlen
          ymn=yorig(ipan)
          ymx=yorig(ipan)+ylen
          set_viewport,xmn,xmx,ymn,ymx
          contour,o3data*1.e6,xdata,thdata,levels=level,/cell_fill,$
                  title=instorder(ipan)+' '+y0+' to '+y1,c_color=col1,$
                  min_value=-999.,xticks=4,xrange=[0.,360.],yrange=[300.,2000.]
          contour,o3data*1.e6,xdata,thdata,levels=findgen(nlvls),/follow,$
                  /overplot,color=0
;
; closest longitude-altitude slice of UKMO data
;
;         result=moment(ydata)
;         y0=result(0)
;         if ndim eq 1 then y0=ydata(0)
;         index=where(abs(alat-y0) le 5.,nlat)
;         mark1=reform(mark2(index(0),*,*))
;         for j=1,nlat-1 do begin
;             llat=index(j)
;             mark1=mark1+reform(mark2(llat,*,*))
;         endfor
;         mark1=mark1/float(nlat)
;         contour,mark1,alon,th,levels=[0.1],/follow,$
;                 /overplot,color=0,thick=3,c_labels=[0]
;         index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=0.
;         contour,mark1,alon,th,levels=[-0.1],/follow,$
;                 /overplot,color=mcolor,thick=3,c_labels=[0]
;         index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=mcolor
;
; HIRDLS
;
          omin=0.
          omax=12.
          xmn=0.6
          xmx=0.6+xlen
          ymn=yorig(ipan)
          ymx=yorig(ipan)+ylen
          !type=2^2+2^3
          set_viewport,xmn,xmx,ymn,ymx
          index=where(hlat2(0,*) gt min(ydata) and hlat2(0,*) lt max(ydata),npt)
          y0=strcompress(string(FORMAT='(F5.1)',min(hlat2(0,index))))
          y1=strcompress(string(FORMAT='(F5.1)',max(hlat2(0,index))))
          xdata=transpose(hlon2(*,index)) & thdata=transpose(htheta(*,index)) 
          o3data=transpose(ho3(*,index)*1.e6)
          index=where(xdata lt 0.)
          if index(0) ne -1 then xdata(index)=xdata(index)+360.
          xsave=xdata
          xdata=0.*xdata
          result=size(thdata)
          nl=result(2)
          for i=0,nl-1 do begin
              sindex=sort(xsave(*,i))
              xdata(*,i)=xsave(sindex,i)
              o3data(*,i)=o3data(sindex,i)
              thdata(*,i)=thdata(sindex,i)
          endfor
          print,'hirdls ',min(o3data),max(o3data)
          contour,o3data,xdata,thdata,levels=level,/cell_fill,$
                  title='HIRDLS '+y0+' to '+y1,c_color=col1,$
                  min_value=-999.,xticks=4,xrange=[0.,360.],yrange=[300.,2000.]
          contour,o3data,xdata,thdata,levels=findgen(nlvls),/follow,$
                  /overplot,color=0

;         contour,mark1,alon,th,levels=[0.1],/follow,$
;                 /overplot,color=0,thick=3,c_labels=[0]
;         index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=0.
;         contour,mark1,alon,th,levels=[-0.1],/follow,$
;                 /overplot,color=mcolor,thick=3,c_labels=[0]
;         index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=mcolor

          noprof:
      endfor
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,0.15,0.9,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[omin,omax],[0,0],yrange=[0,10],xrange=[omin,omax],/noeras,$
          xtitle='(ppmv)',charsize=1.5
    ybox=[0,10,10,0,0]
    x2=omin
    dx=(omax-omin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert xz_hirdls_sage_haloe_poam_o3_'+sdate+$
         '.ps -rotate -90 xz_hirdls_sage_haloe_poam_o3_'+sdate+'.jpg'
endif
jump:

END
