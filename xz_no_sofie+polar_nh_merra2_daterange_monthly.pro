;
; loop over all months
;
; XZ sections of SOFIE NO measurements and vortex frequency and NH polar plots of MERRA-2 vortex
; average over some user specified date range
; VLH 12/20/17
;
@rd_merra2_nc3
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.10]
yorig=[0.10]
cbaryoff=0.03
cbarydel=0.05
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=0
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
longrid=30*findgen(12)
numx=n_elements(longrid)
dlon=longrid(1)-longrid(0)
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
dirs='/atmos/aura3/data/SOFIE_data/Datfiles_SOSST/'
dira='/atmos/aura3/data/ACE_data/Datfiles_SOSST/v3.5-6/'
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
ifile='                             '
mno=[31,28,31,30,31,30,31,31,30,31,30,31]

for iyear=2009,2009 do begin	;2008,2017 do begin
for imon=1,12 do begin

lstmn=imon & lstdy=1 & lstyr=iyear & lstday=0
ledmn=imon & leddy=mno(imon-1) & ledyr=iyear & ledday=0
if iyear mod 4 eq 0 and imon eq 2 then leddy=leddy+1

;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
icounts=0L
lcounta=0L
lcounts=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotavg
;
; read MERRA-2 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
if icount eq 0L then sdate0=sdate
      dum=findfile(dir+sdate+'00.nc3')
      if dum ne '' then ncfile0=dir+sdate+'00.nc3'
      rd_merra2_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag
      if iflag eq 1 then goto,jump
      if icount eq 0L then begin
         rtheta=4000.
;        print,th
;        read,' Enter theta level ',rtheta
         index=where(rtheta eq th)
         if index(0) eq -1 then stop
         itheta=index(0)
         stheta=strcompress(string(fix(th(itheta))),/remove_all)
         xz2d=fltarr(nc,nth)
         yz2d=fltarr(nc,nth)
         for k=0,nth-1 do xz2d(*,k)=alon
         for j=0,nc-1 do yz2d(j,*)=th
         x=fltarr(nc+1)
         x(0:nc-1)=alon
         x(nc)=alon(0)+360.
         x2d=fltarr(nc+1,nr)
         y2d=fltarr(nc+1,nr)
         for i=0,nc do y2d(i,*)=alat
         for j=0,nr-1 do x2d(*,j)=x
      endif
;
; restore SOSST yearly files on 1 January
;
      if iday eq 1L or icount eq 0L then begin
;
; SOFIE
         datesofie_all=[-99.]
         if iyr ge 2007 then begin
            restore,dirs+'cat_sofie_v1.3.'+syr+'.sav'
            restore,dirs+'dmps_sofie_v1.3.merra2.'+syr+'.sav'
            restore,dirs+'no_sofie_v1.3.'+syr+'.sav'
            restore,dirs+'p_sofie_v1.3.'+syr+'.sav'
            restore,dirs+'temp_sofie_v1.3.'+syr+'.sav'
            nz=n_elements(altitude)
            yyyydoysofie_all=date
            ysofie_all=latitude
            xsofie_all=longitude
            modes_all=sctype
            nosofie_all=mix
            thsofie_all=0.*t
            marksofie_all=MARK_PROF
            thmarksofie_all=0.*TP_PROF
            for k=0,nz-1 do thsofie_all(*,k)=t(*,k)*((1000./p(*,k))^(.286))
            for k=0,nz-1 do thmarksofie_all(*,k)=tp_prof(*,k)*((1000./p_prof(*,k))^(.286))
;
; convert YYYYDOY to YYYYMMDD
;
            nprof=n_elements(yyyydoysofie_all)
            syyyydoy=strcompress(yyyydoysofie_all,/r)
            doy=long(strmid(syyyydoy,4,3))
            yy=long(strmid(syyyydoy,0,4))
            datesofie_all=lonarr(nprof)
            for i=0L,nprof-1L do begin
                kdate,float(doy(i)),yy(i),imn,idy
                syr=string(FORMAT='(I4)',yy(i))
                smn=string(FORMAT='(I2.2)',imn)
                sdy=string(FORMAT='(I2.2)',idy)
                datesofie_all(i)=syr+smn+sdy
            endfor
         endif
;
; ACE
         dateace_all=[-99.]
         if iyr ge 2004 then begin
            restore,dira+'cat_ace_v3.5-6.'+syr
            restore,dira+'no_ace_v3.5-6.'+syr
            restore,dira+'temp_ace_v3.5-6.'+syr
            restore,dira+'p_ace_v3.5-6.'+syr
            dateace_all=date
            yace_all=latitude
            xace_all=longitude
            modea_all=sctype
            noace_all=mix
            thace_all=0.*t
            for k=0,nz-1 do thace_all(*,k)=t(*,k)*((1000./p(*,k))^(.286))
         endif

      endif
;
; extract daily SOSST data
;
      norbits=0L & norbita=0L

      sofieday=where(datesofie_all eq long(sdate),norbits)
      if norbits le 1L then goto,jumpsofie
      nosofie=reform(nosofie_all(sofieday,*))
      thsofie=reform(thsofie_all(sofieday,*))
      ysofie=reform(ysofie_all(sofieday))
      xsofie=reform(xsofie_all(sofieday))
      modes=reform(modes_all(sofieday))
      marksofie=reform(marksofie_all(sofieday,*))
      thmarksofie=reform(thmarksofie_all(sofieday,*))
jumpsofie:
       aceday=where(dateace_all eq long(sdate),norbita)
       if norbita le 1L then goto,jumpace
       noace=reform(noace_all(aceday,*))
       thace=reform(thace_all(aceday,*))
       yace=reform(yace_all(aceday))
       xace=reform(xace_all(aceday))
       modea=reform(modea_all(aceday))
jumpace:
;
; remove missing or bad data and separate sunrise and sunset 
;
      ysofie0=-999. & ysofie1=-999.
      if norbits gt 1L then begin
         index=where(ysofie ge 40. and ysofie le 90.,norbits)
         if index(0) ne -1 then begin
            ysofie=ysofie(index)
            xsofie=xsofie(index)
            nosofie=nosofie(index,*)*1.e6
            thsofie=thsofie(index,*)
            marksofie=marksofie(index,*)
            thmarksofie=thmarksofie(index,*)
         endif

         index=where(modes eq 'r')
         if index(0) ne -1 then begin
            ysofiesr=ysofie(index)
            xsofiesr=xsofie(index)
            nosofiesr=nosofie(index,*)
            thsofiesr=thsofie(index,*)
            marksofiesr=marksofie(index,*)
            thmarksofiesr=thmarksofie(index,*)
            index=where(ysofiesr ge 40. and ysofiesr le 90.)
            if index(0) eq -1 then goto,jumpsofiesr
            ysofie0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysofiesr(index))
               ysofie0=result(0)
            endif
         endif
         jumpsofiesr:
         index=where(modes eq 's')
         if index(0) ne -1 then begin
            ysofiess=ysofie(index)
            xsofiess=xsofie(index)
            nosofiess=nosofie(index,*)
            thsofiess=thsofie(index,*)
            marksofiess=marksofie(index,*)
            thmarksofiess=thmarksofie(index,*)
            index=where(ysofiess ge 40. and ysofiess le 90.)
            if index(0) eq -1 then goto,jumpsofiess
            ysofie1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysofiess(index))
               ysofie1=result(0)
            endif
         endif
         jumpsofiess:
      endif

      yace0=-999. & yace1=-999.
      if norbita gt 1L then begin
         index=where(yace ge 40. and yace le 90.,norbita)
         if index(0) ne -1 then begin
            yace=yace(index)
            xace=xace(index)
            noace=noace(index,*)*1.e6
            thace=thace(index,*)
         endif

         index=where(modea eq 'r')
         if index(0) ne -1 then begin
            yacesr=yace(index)
            xacesr=xace(index)
            noacesr=noace(index,*)
            thacesr=thace(index,*)
            index=where(yacesr ge 40. and yacesr le 90.)
            if index(0) eq -1 then goto,jumpacesr
            yace0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yacesr(index))
               yace0=result(0)
            endif
         endif
         jumpacesr:
         index=where(modea eq 's')
         if index(0) ne -1 then begin
            yacess=yace(index)
            xacess=xace(index)
            noacess=noace(index,*)
            thacess=thace(index,*)
            index=where(yacess ge 40. and yacess le 90.)
            if index(0) eq -1 then goto,jumpacess
            yace1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yacess(index))
               yace1=result(0)
            endif
         endif
         jumpacess:
      endif

      yorder=[ysofie0,ysofie1,yace0,yace1]
      instorder=['SOFIESR','SOFIESS','ACESR','ACESS']
;
; remove missing and SH modes
;
      index=where(yorder gt 30.,npan)
      if index(0) eq -1 then goto,jump
      if index(0) ne -1 then begin
         yorder=yorder(index)
         instorder=instorder(index)
      endif
      index=sort(yorder)
      yorder=yorder(index)
      instorder=instorder(index)
      npan0=4
      if npan gt npan0 then begin
         yadd=npan-npan0
         yorder=yorder(yadd:npan-1)
         instorder=instorder(yadd:npan-1)
      endif
      if npan lt npan0 then begin
         yadd=npan0-npan
         yorder=[fltarr(yadd),yorder]
         instorder=[strarr(yadd),instorder]
      endif
;
; postscript file
;
      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,$
                 filename='xz_no_sofie+polar_nh_merra2_'+sdate+'_'+stheta+'K.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; viewport location based on number of panels
;
      erase
      xmn=0.55 & ymn=0.1
      ywide=0.175
      xwide=0.35
      for ipan=0,npan0-1 do begin
          xmx=xmn+xwide
          ymx=ymn+ywide
          set_viewport,xmn,xmx,ymn,ymx
;
; plot longitude-altitude sections in correct order S-N
;
          case 1 of
            (instorder(ipan) eq 'SOFIESR'): begin
             ydata=ysofiesr
             xdata=xsofiesr
             nodata=nosofiesr
             thdata=thsofiesr
             markdata=marksofiesr
             thmarkdata=thmarksofiesr
            end
            (instorder(ipan) eq 'SOFIESS'): begin
             ydata=ysofiess
             xdata=xsofiess
             nodata=nosofiess
             thdata=thsofiess
             markdata=marksofiess
             thmarkdata=thmarksofiess
            end
            (instorder(ipan) eq 'ACESR'): begin
             ydata=yacesr
             xdata=xacesr
             nodata=noacesr
             thdata=thacesr
            end
            (instorder(ipan) eq 'ACESS'): begin
             ydata=yacess
             xdata=xacess
             nodata=noacess
             thdata=thacess
            end
            else: begin
            goto,noprof
            end
          endcase
          ydata=reform(ydata)
          xdata=reform(xdata)
          nodata=reform(nodata)
          thdata=reform(thdata)
; 
; remove bad xdata and sort in longitude
;
          index=where(xdata ge 0. and xdata le 360.)
          if index(0) ne -1 then begin
             ydata=reform(ydata(index))
             xdata=reform(xdata(index))
             nodata=reform(nodata(index,*))
             thdata=reform(thdata(index,*))
          endif
          xsave=xdata
          xdata=0.*thdata
          result=size(thdata)
          ndim=result(0)
          if ndim eq 1 then goto,oneprof
          nprof=result(1)
          nl=result(2)
          for i=0,nl-1 do begin
              sindex=sort(xsave)
              xdata(*,i)=xsave(sindex)
              nodata(*,i)=nodata(sindex,i)
              thdata(*,i)=thdata(sindex,i)
          endfor
          xlabels='Lat='+string(ydata(sindex),format='(f4.1)')
          level=[10.,20.,30.,40.,50.,70.,100.,150.,200.,500.,1000.,2000.,$
                 5000.,10000.,20000.,50000.,100000.,1.e6]/10000.
          nlvls=n_elements(level)
          col1=1+indgen(nlvls)*mcolor/nlvls
          !type=2^2+2^3
          contour,nodata,xdata,thdata,levels=level,/cell_fill,$
                  title=instorder(ipan)+'      '+xlabels(0),c_color=col1,$
                  min_value=-99.,xticks=4,xrange=[0.,360.],yrange=[500.,15000.],$
                  ytitle='Theta',color=0
          contour,nodata,xdata,thdata,levels=level,/follow,$
                  /overplot,color=0,min_value=-99.
;
; closest longitude-altitude slice of MERRA-2 data
;
          result=moment(ydata)
          y0=result(0)
          oneprof:
          if ndim eq 1 then y0=ydata(0)
          index=where(abs(alat-y0) le 7.5,nlat)
          pv1=reform(pv2(index(0),*,*))
          p1=reform(p2(index(0),*,*))
          mark1=reform(mark2(index(0),*,*))
          for j=1,nlat-1 do begin
          llat=index(j)
          pv1=pv1+reform(pv2(llat,*,*))
          p1=p1+reform(p2(llat,*,*))
          mark1=mark1+reform(mark2(llat,*,*))
          endfor
          pv1=pv1/float(nlat)
          p1=p1/float(nlat)
          mark1=mark1/float(nlat)
          t1=0.*mark1
          for k=0,nth-1 do t1(*,k)=th(k)*((p1(*,k)/1000.)^(.286))
          contour,smooth(mark1,3),alon,th,levels=[0.1],/follow,$
                  /overplot,color=0.9*mcolor,thick=6,c_labels=[0]
          index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=0.
          contour,smooth(mark1,3),alon,th,levels=[-0.1],/follow,$
                  /overplot,color=0,thick=6,c_labels=[0]
          index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=mcolor
;
; bin profiles in longitude for later averaging
;
            if instorder(ipan) eq 'SOFIESR' then begin
if lcounts eq 0L then begin
   nodata_avg_sofie=fltarr(numx,nz)
   thdata_avg_sofie=fltarr(numx,nz)
   nnodata_tot_sofie=fltarr(numx,nz)
   markdata_avg_sofie=fltarr(nc,nth)
endif
index=where(mark1 gt 0.)
if index(0) ne -1L then markdata_avg_sofie(index)=markdata_avg_sofie(index)+mark1(index)
lcounts=lcounts+1

for i=0L,nprof-1L do begin
    noprof=reform(nodata(i,*))
    thprof=reform(thdata(i,*))
    xpt=xdata(i,60)
    for ii=0L,numx-1L do begin
        if ii eq 0L then begin
           if xpt gt longrid(-1)+dlon/2. or xpt lt longrid(ii)+dlon/2. then begin
              index=where(noprof gt 0. and finite(thprof) eq 1)
              if index(0) ne -1L then begin
                 nnodata_tot_sofie(ii,index)=nnodata_tot_sofie(ii,index)+1.0
                 nodata_avg_sofie(ii,index)=nodata_avg_sofie(ii,index)+noprof(index)
                 thdata_avg_sofie(ii,index)=thdata_avg_sofie(ii,index)+thprof(index)
              endif
           endif
        endif

        if ii gt 0L then begin
           if xpt gt longrid(ii)-dlon/2. and xpt lt longrid(ii)+dlon/2. then begin
              index=where(noprof gt 0. and finite(thprof) eq 1)
              if index(0) ne -1L then begin
                 nnodata_tot_sofie(ii,index)=nnodata_tot_sofie(ii,index)+1.0
                 nodata_avg_sofie(ii,index)=nodata_avg_sofie(ii,index)+noprof(index)
                 thdata_avg_sofie(ii,index)=thdata_avg_sofie(ii,index)+thprof(index)
              endif
           endif
        endif

    endfor
endfor
            endif
            if instorder(ipan) eq 'ACESR' or instorder(ipan) eq 'ACESS' then begin
if lcounta eq 0L then begin
   nodata_avg_ace=fltarr(numx,nz)
   thdata_avg_ace=fltarr(numx,nz)
   nnodata_tot_ace=fltarr(numx,nz)
   markdata_avg_ace=fltarr(nc,nth)
endif
index=where(mark1 gt 0.)
if index(0) ne -1L then markdata_avg_ace(index)=markdata_avg_ace(index)+mark1(index)
lcounta=lcounta+1

for i=0L,nprof-1L do begin
    noprof=reform(nodata(i,*))
    thprof=reform(thdata(i,*))
    xpt=xdata(i,60)
    for ii=0L,numx-1L do begin
        if ii eq 0L then begin
           if xpt gt longrid(-1)+dlon/2. or xpt lt longrid(ii)+dlon/2. then begin
              index=where(noprof gt 0. and finite(thprof) eq 1)
              if index(0) ne -1L then begin
                 nnodata_tot_ace(ii,index)=nnodata_tot_ace(ii,index)+1.0
                 nodata_avg_ace(ii,index)=nodata_avg_ace(ii,index)+noprof(index)
                 thdata_avg_ace(ii,index)=thdata_avg_ace(ii,index)+thprof(index)
              endif
           endif
        endif

        if ii gt 0L then begin
           if xpt gt longrid(ii)-dlon/2. and xpt lt longrid(ii)+dlon/2. then begin
              index=where(noprof gt 0. and finite(thprof) eq 1)
              if index(0) ne -1L then begin
                 nnodata_tot_ace(ii,index)=nnodata_tot_ace(ii,index)+1.0
                 nodata_avg_ace(ii,index)=nodata_avg_ace(ii,index)+noprof(index)
                 thdata_avg_ace(ii,index)=thdata_avg_ace(ii,index)+thprof(index)
              endif
           endif
        endif

    endfor
endfor

            endif

          plots,0.,rtheta,/data
          plots,360.,rtheta,/data,/continue,color=0
          imin=min(level)
          imax=max(level)
          xmnb=xmx+.03
          xmxb=xmnb+.05
          set_viewport,xmnb,xmxb,ymn,ymx
          !type=2^2+2^3+2^5
slab=' '+strarr(n_elements(level))
plot,[0,0],[min(level),max(level)],xrange=[0,10],color=0,$
     yticks=n_elements(level)-1L,ytickname=slab,$
     yrange=[min(level),max(level)],charsize=1.2,title='NO (ppmv)',charthick=2
xbox=[0,10,10,0,0]
y1=min(level)
dy=(max(level)-min(level))/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
slab=strcompress(string(format='(f7.3)',level),/remove_all)
y1=min(level)+dy/2
for i=0L,n_elements(slab)-1L do begin
    slab0=slab(i)
    flab0=float(slab(i))
    if flab0 lt 0.01 then begin
       slab0=strcompress(string(format='(f5.3)',flab0),/remove_all)
       xyouts,xorig(0)+xwide+0.02,y1,slab0,charsize=1.2,/data,color=mcolor,charthick=2
    endif
    if flab0 lt 1. and flab0 ge 0.01 then begin
       slab0=strcompress(string(format='(f4.2)',flab0),/remove_all)
       xyouts,xorig(0)+xwide+0.02,y1,slab0,charsize=1.2,/data,color=0,charthick=2
    endif
    if flab0 ge 1. then begin
       slab0=strcompress(long(slab0),/remove_all)
       xyouts,xorig(0)+xwide+0.02,y1,slab0,charsize=1.2,/data,color=0,charthick=2
    endif
    y1=y1+dy
endfor
          noprof:
          ymn=ymx+0.05
      endfor		; loop over XZ panels
      pv1=transpose(pv2(*,*,itheta))
      p1=transpose(p2(*,*,itheta))
      sf1=transpose(sf2(*,*,itheta))
      mark1=transpose(mark2(*,*,itheta))
      pv=fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1
      pv(nc,*)=pv(0,*)
      p=fltarr(nc+1,nr)
      p(0:nc-1,0:nr-1)=p1
      p(nc,*)=p(0,*)
      sf=fltarr(nc+1,nr)
      sf(0:nc-1,0:nr-1)=sf1
      sf(nc,*)=sf(0,*)
      mark=fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
      !type=2^2+2^3
      set_viewport,.05,.45,.45,.85
      MAP_SET,90,0,-180,/stereo,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5,$
              title=sdate+'   '+stheta+' K'
      oplot,findgen(361),0.1+0.*findgen(361),psym=0,color=0
      contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark gt 0. and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark lt 0. and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      contour,mark,x,alat,levels=[.1],c_color=mcolor*.1,/overplot,/follow,c_labels=0,/noeras,thick=3
      contour,mark,x,alat,levels=[-.1],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=3

if icount eq 0L then begin
   mark_avg=0.*mark
   sf_avg=0.*sf
endif
mark_avg=mark_avg+mark
sf_avg=sf_avg+sf

      nomin=min(level)
      nomax=max(level)
      if norbits gt 1L then begin
         xsofie2d=0.*thsofie
         ysofie2d=0.*thsofie
         for k=0L,nz-1L do begin
             xsofie2d(*,k)=xsofie
             ysofie2d(*,k)=ysofie
         endfor
;        index=where(abs(thsofie-rtheta) le 100. and nosofie gt 0. and ysofie2d gt 0.,kcount)
         index=where(abs(thsofie-rtheta) le 100. and ysofie2d gt 0.,kcount)
         if kcount gt 0L then begin
            xday=xsofie2d(index)
            yday=ysofie2d(index)
            thday=thsofie(index)
            noday=nosofie(index)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,kcount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=250	;mcolor*(noday(i)-nomin)/(nomax-nomin)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
         if icounts eq 0L then begin
            xsave_sofie=xday
            ysave_sofie=yday
            icounts=1
         endif
         if icount gt 0L then begin
            xsave_sofie=[xsave_sofie,xday]
            ysave_sofie=[ysave_sofie,yday]
         endif
      endif
      if norbita gt 1L then begin
         xace2d=0.*thace
         yace2d=0.*thace
         for k=0L,nz-1L do begin
             xace2d(*,k)=xace
             yace2d(*,k)=yace
         endfor
         index=where(abs(thace-rtheta) le 100. and yace2d gt 0.,hcount)
         if hcount gt 0L then begin
            xday=xace2d(index)
            yday=yace2d(index)
            thday=thace(index)
            noday=noace(index)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,hcount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=150	;mcolor*(noday(i)-nomin)/(nomax-nomin)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
         if icount eq 0L then begin
            xsave_ace=xday
            ysave_ace=yday
         endif
         if icount gt 0L then begin
            xsave_ace=[xsave_ace,xday]
            ysave_ace=[ysave_ace,yday]
         endif
      endif

      imin=nomin
      imax=max(level)
      ymnb=.45-cbaryoff
      ymxb=ymnb+cbarydel
      set_viewport,.05,.45,ymnb,ymxb
      !type=2^2+2^3+2^6
      slab=' '+strarr(n_elements(level))
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
            xtitle='NO (ppmv)',charsize=1.2,color=0,xtickname=slab,charthick=2
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=col1(j)
          x1=x1+dx
      endfor
slab=strcompress(string(format='(f6.3)',level),/remove_all)
x1=min(level)+dx/2
for i=0L,n_elements(slab)-1L do begin
    slab0=slab(i)
    if i lt 5 then xyouts,x1,10,slab0,charsize=1.2,/data,color=mcolor,charthick=2,orientation=270
    if i ge 5 then xyouts,x1,10,slab0,charsize=1.2,/data,color=0,charthick=2,orientation=270
    x1=x1+dx
endfor

;     if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim xz_no_sofie+polar_nh_merra2_'+sdate+'_'+stheta+'K.ps -rotate -90 '+$
               ' xz_no_sofie+polar_nh_merra2_'+sdate+'_'+stheta+'K.jpg'
;        spawn,'/usr/bin/rm xz_no_sofie+polar_nh_merra2_'+sdate+'_'+stheta+'K.ps'
      endif
      icount=icount+1L
      goto,jump

plotavg:
sdate1=sdate
daterange=sdate0+'-'+sdate1
mark_avg=mark_avg/float(icount)
sf_avg=sf_avg/float(icount)

nodata_avg_sofie=nodata_avg_sofie/nnodata_tot_sofie
thdata_avg_sofie=thdata_avg_sofie/nnodata_tot_sofie
markdata_avg_sofie=markdata_avg_sofie/float(lcounts)

if lcounta gt 0L then begin
nodata_avg_ace=nodata_avg_ace/nnodata_tot_ace
thdata_avg_ace=thdata_avg_ace/nnodata_tot_ace
markdata_avg_ace=markdata_avg_ace/float(lcounta)
endif
erase
      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,$
                 filename='xz_no_sofie+polar_nh_merra2_'+daterange+'_'+stheta+'K.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif

      xmn=0.55 & ymn=0.5
      ywide=0.175
      xwide=0.35
      xmx=xmn+xwide
      ymx=ymn+ywide
      set_viewport,xmn,xmx,ymn,ymx

;     xlabels='Lat='+string(ydata(sindex),format='(f4.1)')
      level=[10.,20.,30.,40.,50.,70.,100.,150.,200.,$
             500.,1000.,2000.,5000.,10000.,20000.,50000.,100000.,1.e6]/10000.
      nlvls=n_elements(level)
      col1=1+indgen(nlvls)*mcolor/nlvls
      !type=2^2+2^3
;
; wrap around point in longitude
;
longrid1=fltarr(numx+1)
longrid1(0:numx-1)=longrid
longrid1(numx)=longrid(0)+360.

if lcounta gt 0L then begin

nodata_avg_ace1=fltarr(numx+1,nz)
thdata_avg_ace1=fltarr(numx+1,nz)
nodata_avg_ace1(0:numx-1,*)=nodata_avg_ace
thdata_avg_ace1(0:numx-1,*)=thdata_avg_ace
nodata_avg_ace1(numx,*)=nodata_avg_ace(0,*)
thdata_avg_ace1(numx,*)=thdata_avg_ace(0,*)

      contour,nodata_avg_ace1,longrid1,thdata_avg_ace1,levels=level,/cell_fill,$
              title='ACE '+strcompress(mean(ysave_ace),/r)+' '+daterange,c_color=col1,$
              min_value=-99.,xticks=4,xrange=[0.,360.],yrange=[1000.,15000.],$
              ytitle='Theta',color=0
      contour,nodata_avg_ace1,longrid1,thdata_avg_ace1,levels=[0.02,0.1,0.5,5],/follow,/overplot,color=0,min_value=-99.,thick=2
      contour,markdata_avg_ace,alon,th,levels=[0.5],/follow,$
                  /overplot,c_color=[250],thick=3,c_labels=[0]
      index=where(markdata_avg_ace gt .5)
      if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=0,symsize=0.5
thprof=mean(thdata_avg_ace,dim=1,/nan)
yindex=where(thprof ge 1000. and thprof le 15000.)
axis,yrange=[min(altitude(yindex)),max(altitude(yindex))],/save,color=0,yaxis=1,ytitle='                 Approx. Altitude (km)',charsize=1.5,charthick=2
endif

      ymn=ymx+0.05
      ymx=ymn+ywide
      set_viewport,xmn,xmx,ymn,ymx
nodata_avg_sofie1=fltarr(numx+1,nz)
thdata_avg_sofie1=fltarr(numx+1,nz)
nodata_avg_sofie1(0:numx-1,*)=nodata_avg_sofie
thdata_avg_sofie1(0:numx-1,*)=thdata_avg_sofie
nodata_avg_sofie1(numx,*)=nodata_avg_sofie(0,*)
thdata_avg_sofie1(numx,*)=thdata_avg_sofie(0,*)

      contour,nodata_avg_sofie1,longrid1,thdata_avg_sofie1,levels=level,/cell_fill,$
              title='SOFIESR '+strcompress(mean(ysave_sofie),/r)+' '+daterange,c_color=col1,$
              min_value=-99.,xticks=4,xrange=[0.,360.],yrange=[1000.,15000.],$
              ytitle='Theta',color=0
      contour,nodata_avg_sofie1,longrid1,thdata_avg_sofie1,levels=[0.02,0.1,0.5,5],/follow,/overplot,color=0,min_value=-99.,thick=2
      contour,markdata_avg_sofie,alon,th,levels=[0.5],/follow,$
                  /overplot,c_color=[250],thick=3,c_labels=[0]
      index=where(markdata_avg_sofie gt .5)
      if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=0,symsize=0.5
thprof=mean(thdata_avg_sofie,dim=1,/nan)
yindex=where(thprof ge 1000. and thprof le 15000.)
axis,yrange=[min(altitude(yindex)),max(altitude(yindex))],/save,color=0,yaxis=1,charsize=1.5,charthick=2


      !type=2^2+2^3
      set_viewport,.05,.45,.45,.85
      MAP_SET,90,0,-180,/stereo,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5,$
              title=month(long(smn-1))+' '+syr+'   '+stheta+' K'
      oplot,findgen(361),0.1+0.*findgen(361),psym=0,color=0
      contour,sf_avg,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark_avg gt 0.1 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark_avg lt -0.1 and y2d gt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      contour,mark_avg,x,alat,levels=[.1],c_color=mcolor*.1,/overplot,/follow,c_labels=0,/noeras,thick=3
      contour,mark_avg,x,alat,levels=[-.1],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=3

      a=findgen(6)*(2*!pi/5.)
      usersym,2.*cos(a),2.*sin(a),/fill
      oplot,xsave_sofie,ysave_sofie,psym=8,color=250  ;mcolor*(noday(i)-nomin)/(nomax-nomin)
      a=findgen(6)*(2*!pi/5.)
      usersym,2.*cos(a),2.*sin(a)
;     oplot,xsave_sofie,ysave_sofie,psym=8,color=lc

      a=findgen(6)*(2*!pi/5.)
      usersym,2.*cos(a),2.*sin(a),/fill
      if lcounta gt 0L then oplot,xsave_ace,ysave_ace,psym=8,color=150  ;mcolor*(noday(i)-nomin)/(nomax-nomin)
      a=findgen(6)*(2*!pi/5.)
      usersym,2.*cos(a),2.*sin(a)
      if lcounta gt 0L then oplot,xsave_ace,ysave_ace,psym=8,color=lc

      imin=nomin
      imax=max(level)
      ymnb=.45-cbaryoff
      ymxb=ymnb+cbarydel
      set_viewport,.45,.95,ymnb,ymxb
      !type=2^2+2^3+2^6
      slab=' '+strarr(n_elements(level))
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
            xtitle='NO (ppmv)',charsize=1.2,color=0,xtickname=slab,charthick=2
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=col1(j)
          x1=x1+dx
      endfor
slab=strcompress(string(format='(f7.3)',level),/remove_all)
x1=min(level)+dx/2
for i=0L,n_elements(slab)-1L do begin
    slab0=slab(i)
    if level(i) ge 1. then slab0=strcompress(string(format='(i3)',level(i)),/remove_all)
    if i lt 5 then xyouts,x1,10,slab0,charsize=1.2,/data,color=mcolor,charthick=2,orientation=270
    if i ge 5 and level(i) lt 1. then xyouts,x1,10,slab0,charsize=1.2,/data,color=0,charthick=2,orientation=270
    if i ge 5 and level(i) ge 1. then xyouts,x1,7,slab0,charsize=1.2,/data,color=0,charthick=2,orientation=270
    x1=x1+dx
endfor

      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim xz_no_sofie+polar_nh_merra2_'+daterange+'_'+stheta+'K.ps -rotate -90 '+$
               ' xz_no_sofie+polar_nh_merra2_'+daterange+'_'+stheta+'K.jpg'
;        spawn,'/usr/bin/rm xz_no_sofie+polar_nh_merra2_'+daterange+'_'+stheta+'K.ps'
      endif

skipmonth:
endfor	; loop over months
endfor	; loop over years
end
