
; check new nc data

@rd_WA3548T08CO_nc
@stddat
@kgmt
@ckday
@kdate

loadct,38
mcolor=byte(!p.color)
device,decompose=0
icmm1=mcolor-1B
icmm2=mcolor-2B
!noeras=1
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=750
nydim=750
xorig=[0.10]
yorig=[0.25]
xlen=0.8
ylen=0.5
cbaryoff=0.1
cbarydel=0.02
!NOERAS=-1
SETPLOT='ps'
read,'setplot',setplot
if setplot ne 'ps' then begin
   lc=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura7/harvey/WACCM_data/Datfiles/Datfiles_Liu/WA3548T08CO_2x.cam2.h2.'
uday=0L & lstday=0L & ledday=0L
lstmn=1L & lstdy=1L & lstyr=2002L 
ledmn=1L & leddy=1L & ledyr=2002L
mon=['jan_','feb_','mar_','apr_','may_','jun_','jul_','aug_','sep_','oct_','nov_','dec_']
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      ifile=dir+syr+'-'+smn+'-'+sdy+'-MetO.nc'
      rd_WA3548T08CO_nc,ifile,nc,nr,nth,alon,alat,th,pv1,p1,msf1,u1,v1,q1,qdf1,no1,co1,e1,iflag
;
; plot isentropic pressure at each level
;
      for ilev=0L,nth-1L do begin
;         grd1=transpose(pv1(*,*,ilev))
          grd1=transpose(e1(*,*,ilev))
          stheta=strtrim(string(long(th(ilev))),2)

          if setplot eq 'ps' then begin
             lc=0
             set_plot,'ps'
             xsize=nxdim/100.
             ysize=nydim/100.
             !psym=0
             !p.font=0
             device,font_size=9
             device,/landscape,bits=8,filename='merc_WA3548T08CO_pv_'+sdate+'_'+stheta+'K.ps'
             device,/color
             device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                    xsize=xsize,ysize=ysize
          endif

          nlvls=21
          index=where(grd1 ne 0. and grd1 ne 1.00000e+12)
          imin=min(grd1(index)) & imax=max(grd1(index))
          level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)

          col1=1+indgen(nlvls)*mcolor/nlvls
          !noeras=1
          erase
          xyouts,.4,.8,sdate,/normal,charsize=3
          xmn=xorig(0)
          xmx=xorig(0)+xlen
          ymn=yorig(0)
          ymx=yorig(0)+ylen
          set_viewport,xmn,xmx,ymn,ymx
          !type=2^2+2^3
          map_set,0,0,0,/grid,/contin,/noeras,color=mcolor,title='WA3548T08CO PV at '+stheta+' K',charsize=2
          contour,grd1,alon,alat,levels=level,/fill,/cell_fill,c_color=col1,/noeras,$
              xticks=6,yrange=[-90.,90.],xrange=[0.,360.],charsize=2,/overplot
index=where(level gt 0.)
if index(0) ne -1L then contour,grd1,alon,alat,levels=level(index),/follow,c_color=0,/overplot,/noeras
index=where(level lt 0.)
if index(0) ne -1L then contour,grd1,alon,alat,levels=level(index),/follow,c_color=mcolor,/overplot,/noeras,c_linestyle=5
          map_set,0,0,0,/grid,/contin,/noeras,color=mcolor
          map_set,0,0,0,/grid,/contin,/noeras 
          imin=min(level)
          imax=max(level)
          ymnb=yorig(0) -cbaryoff
          ymxb=ymnb  +cbarydel
          set_viewport,xmn,xmx,ymnb,ymxb
          !type=2^2+2^3+2^6
          plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],charsize=1.5,xtitle=''
          ybox=[0,10,10,0,0]
          x1=imin
          dx=(imax-imin)/float(nlvls)
          for j=0,nlvls-1 do begin
              xbox=[x1,x1,x1+dx,x1+dx,x1]
              polyfill,xbox,ybox,color=col1(j)
              x1=x1+dx
          endfor
      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim merc_WA3548T08CO_pv_'+sdate+'_'+stheta+'K.ps -rotate -90 '+$
               'merc_WA3548T08CO_pv_'+sdate+'_'+stheta+'K.jpg'
         spawn,'/usr/bin/rm merc_WA3548T08CO_pv_'+sdate+'_'+stheta+'K.ps'
      endif

      endfor
goto,jump
end
