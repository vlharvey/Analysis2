; 
; SAVE WA3 LOPs
; polar and xz sections of WACCM 3 ozone + anticyclones and Arctic vortex
;
@stddat
@kgmt
@ckday
@kdate
@rd_waccm3_nc3
@range_ring

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
!NOERAS=-1
device,decompose=0
nxdim=700
nydim=700
xorig=[0.1]
yorig=[0.35]
xlen=0.8
ylen=0.8
cbaryoff=0.03
cbarydel=0.02
!NOERAS=-1
lstmn=1
lstdy=1
lstyr=1991
ledmn=1
leddy=1
ledyr=2004
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
print, ' '
print, '      WACCM Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 57 then lstyr=lstyr+2000
if ledyr lt 57 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1957 then stop,'Year out of range '
if ledyr lt 1957 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
RADG = !PI / 180.
FAC20 = 1.0 / TAN(45.*RADG)
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
!noeras=1
dir='/aura7/harvey/WACCM_data/Datfiles/wa3_tnv3_'

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
nday=ledday-lstday+1L
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then begin
         save,file='WA3_LOP_ZT_nh.sav',date_zt,lop_zt,ho3mean_zt,ao3mean_zt,ho3sigma_zt,ao3sigma_zt
         stop,' Normal termination condition '
      endif
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
      date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
                              month(imn-1),' ',idy,', ',iyr))
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      print,syr+smn+sdy
      uyr=strmid(syr,2,2)
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      ifile=string(FORMAT='(i4.4,i2.2,i2.2,a4)',iyr,imn,idy,'.nc3')
      rd_waccm3_nc3,dir+ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,sf2,o32,ch42,no22,h2o2,iflag
;     rd_waccm3_nc,dir+ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
;        u2,v2,qdf2,o32,ch42,no22,h2o2,iflag
      if iflag eq 1 then goto,jump
;
; read new marker field
;
;     file1=dir+string(FORMAT='(i4.4,i2.2,i2.2,a4)',iyr,imn,idy,'.nc4')
;     dum1=findfile(file1)
;     if dum1(0) ne '' then begin
;        ncid=ncdf_open(file1)
;        print,'opened '+file1,ncid
;        mark2=fltarr(nr,nc,nth)
;        ncdf_varget,ncid,3,mark2
;        ncdf_close,ncid
;     endif
;     if dum1(0) eq '' then stop,'cannot find nc4 file'
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.

; select theta level
    if icount eq 0L then begin
       rlev=1000.
;      print,th
;      read,'Enter theta surface ',rlev
       zindex=where(th eq rlev)
       ilev=zindex(0)
       slev=strcompress(string(fix(th(ilev))),/remove_all)+'K'
       x=fltarr(nc+1)
       x(0:nc-1)=alon
       x(nc)=alon(0)+360.
       x2d=fltarr(nc+1,nr)
       y2d=fltarr(nc+1,nr)
       for i=0,nc do y2d(i,*)=alat
       for j=0,nr-1 do x2d(*,j)=x
       xz2d=fltarr(nc,nth)
       yz2d=fltarr(nc,nth)
       for k=0,nth-1 do xz2d(*,k)=alon
       for j=0,nc-1 do yz2d(j,*)=th
;
; zt arrays
;
       lop_zt=fltarr(nday,nth)-99.
       ho3mean_zt=fltarr(nday,nth)-99.
       ao3mean_zt=fltarr(nday,nth)-99.
       ho3sigma_zt=fltarr(nday,nth)-99.
       ao3sigma_zt=fltarr(nday,nth)-99.
       hthmean_zt=fltarr(nday,nth)-99.
       athmean_zt=fltarr(nday,nth)-99.
       hthsigma_zt=fltarr(nday,nth)-99.
       athsigma_zt=fltarr(nday,nth)-99.
       date_zt=strarr(nday)
;
; area of each gridpoint
;
         adata=fltarr(nc,nr)
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         hy=re*deltay*dtr
         for j=0L,nr-1L do begin
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             for i=0L,nc-1L do begin
                 adata(i,j)=dx*hy
             endfor
         endfor
    endif
    date_zt(icount)=lfile
;
; build WACCM ozone, mask, lon, lat, theta profile arrays like MLS (nprof,nth)
;
    iprof=0L
    nprof=nr*nc
    wo3=fltarr(nprof,nth)
    wmark2=fltarr(nprof,nth)
    wlon2=fltarr(nprof,nth)
    wlat2=fltarr(nprof,nth)
    wtheta2=fltarr(nprof,nth)
    for i=0L,nc-1L do begin
    for j=0L,nr-1L do begin
        wo3(iprof,*)=o32(j,i,*)
        wmark2(iprof,*)=mark2(j,i,*)
        wlon2(iprof,*)=alon(i)
        wlat2(iprof,*)=alat(j)
        wtheta2(iprof,*)=th
        iprof=iprof+1L
    endfor
    endfor

; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='Figures/lop_o3_wa3_nhDL_'+ifile+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif

    mark1=transpose(mark2(*,*,ilev))
;   sf1=transpose(sf2(*,*,ilev))
;sfmax=max(sf1(*,nr/2:nr-1))
;sfmax=sfmax-0.1*sfmax
    o31=transpose(o32(*,*,ilev))
;   sf=0.*fltarr(nc+1,nr)
;   sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
;   sf(nc,*)=sf(0,*)
    mark=0.*fltarr(nc+1,nr)
    mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
    mark(nc,*)=mark(0,*)
    o3=0.*fltarr(nc+1,nr)
    o3(0:nc-1,0:nr-1)=o31(0:nc-1,0:nr-1)
    o3(nc,*)=o3(0,*)
    erase
    xmn=0.1
    xmx=0.5
    ymn=0.5
    ymx=0.9
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    MAP_SET,90,-90,0,/stereo,/noeras,/grid,/contin,title='WA3 O3',charsize=1.5,color=0,/noborder
    oplot,findgen(361),0.1+0.*findgen(361)
    nlvls=19
    level=1.0+0.5*findgen(nlvls)
    o3min=min(level)
    o3max=max(level)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,o3,x,alat,/overplot,levels=level,c_color=col1,/cell_fill,/noeras
    contour,o3,x,alat,/overplot,levels=level,/follow,c_labels=0*level,/noeras,color=0
    loadct,0
    contour,mark,x,alat,/overplot,levels=[0.1],thick=7,color=150,c_labels=[0]
    index=where(mark gt 0. and y2d gt 0.)
    if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=1,color=150,symsize=0.5
    loadct,38
    contour,mark,x,alat,/overplot,levels=[-0.1],thick=7,color=0,c_labels=[0]
    index=where(mark lt 0. and y2d gt 0.)
    if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=1,color=0,symsize=0.5
    MAP_SET,90,-90,0,/stereo,/noeras,/grid,/contin,charsize=1.5,color=0,/noborder
    oplot,findgen(361),0.2+0.*findgen(361),color=0,psym=8,symsize=0.5
    imin=min(level)
    imax=max(level)
    ymnb=ymn-cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.05,xmx-0.05,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle='Ozone (ppmv)',charsize=1.5,color=0
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(j)
        x1=x1+dx
    endfor
;
; ***ADD MLS LOP ALGORITHM HERE ***
;
; loop over theta levels ranging from 2000-500K
; save anticyclone points and nearby ambient points (t,x,y,th,o3)
;
    !type=2^2+2^3
    set_viewport,.5,.9,.5,.9
    map_set,90,-90,0,/stereo,/contin,/grid,color=0,title=syr+smn+sdy,charsize=2,/noeras,/noborder
    oplot,findgen(361),0.2+0.*findgen(361),color=0,psym=8,symsize=0.5
    hareamean_save=fltarr(nth)
    ho3mean_save=fltarr(nth)
    ao3mean_save=fltarr(nth)
    hthmean_save=fltarr(nth)
    athmean_save=fltarr(nth)
    ho3sigma_save=fltarr(nth)
    ao3sigma_save=fltarr(nth)
    hthsigma_save=fltarr(nth)
    athsigma_save=fltarr(nth)
    ao3sigma_save=fltarr(nth)
    ho3_save=-99.
    ao3_save=-99.
    hth_save=-99.
    ath_save=-99.
    hx_save=-99.
    ax_save=-99.
    hy_save=-99.
    ay_save=-99.
    nlvls=nth
    col1=reverse(1+indgen(nlvls)*mcolor/(float(nlvls)))
    for k=26L,15L,-1L do begin		; 500 K to 2000 K
;
; extract MetO anticyclone info
;
        rtheta2=th(k)
        mark1=transpose(mark2(*,*,k))
        mark=fltarr(nc+1,nr)
        mark(0:nc-1,0:nr-1)=mark1
        mark(nc,*)=mark(0,*)
        if rtheta2 eq th(ilev) then begin
        loadct,0
        contour,mark,x,alat,levels=[.05],c_color=150,/overplot,/follow,c_labels=0,/noeras,thick=5
        contour,mark,x,alat,levels=[-.05],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=5
        loadct,38
        endif

        rtheta2=th(k)
        index=where(rtheta2 eq th)
        itheta=index(0)
;
; look in NH for anticyclones
;
        index=where(wtheta2 eq rtheta2 and wlat2 gt 20.,npt)
        if index(0) eq -1L then goto,jumplev
        xdata=wlon2(index) & ydata=wlat2(index) 
        o3data=wo3(index)
        o3mark=wmark2(index)
        o3theta=wtheta2(index)
        index=where(xdata lt 0.)
        if index(0) ne -1 then xdata(index)=xdata(index)+360.
;
; if there are no anticyclones on this level then jumplev
;
        if min(o3mark) eq 0. then goto,jumplev
        mmin=min(o3mark)
;
; find highest latitude anticyclone
;
        if mmin lt -1L then begin
           ymax=0.
           for imin=mmin,-1L do begin
               index=where(o3mark eq imin)
               if index(0) ne -1L then begin
               if max(ydata(index)) gt ymax then ymax=max(ydata(index))
               endif
           endfor
           index=where(o3mark lt 0. and ydata eq ymax)
           mmin=min(o3mark(index))
        endif
        hindex=where(o3mark eq mmin)
        if n_elements(hindex) lt 2L then goto,jumplev
        resultx=moment(xdata(hindex))
;
; GM logic
;
        if min(xdata(hindex)) lt alon(1) and max(xdata(hindex)) gt alon(nc-2) then begin
           xsave=xdata(hindex)
           xindex=where(xsave gt 180.)
           xsave(xindex)=xsave(xindex)-360.
           resultx=moment(xsave)
        endif
        x0=resultx(0)
        if x0 lt 0. then x0=x0+360.
        resulty=moment(ydata(hindex))
        y0=resulty(0)
        if y0 le 20. then goto,jumplev                ; ignore tropical anticyclones
if rtheta2 eq th(ilev) then begin
   oplot,xdata(hindex),ydata(hindex),psym=1,color=230
   oplot,[x0,x0],[y0,y0],psym=8,symsize=3,color=0
endif
;
; calculate area of the anticyclone
;
hareamean_save(k)=total(adata(hindex))
;
; calculate distance from anticyclone center of mass to all NH data points
;
; GM logic
;
        dxf=re*abs(x0-xdata)*dtr*cos(y0*dtr)
        if min(xdata) eq alon(0) and max(xdata) eq alon(nc-1) then begin
           xsave=xdata
           if x0 gt 180. then begin
              index=where(xsave lt 180.)
              xsave(index)=xsave(index)+360.
           endif
           if x0 lt 180. then begin
              index=where(xsave gt 180.)
              xsave(index)=xsave(index)-360.
           endif
           dxx=abs(x0-xsave)
           index=where(dxx gt 360.)
           if index(0) ne -1L then dxx(index)=dxx(index)-360.
           dxf=re*dxx*dtr*cos(y0*dtr)
        endif
        dyf=re*abs(y0-ydata)*dtr
        dist=sqrt(dxf*dxf+dyf*dyf)
        dist0=max(dist)+1000.
        if dist0 gt 5000. then dist0=5000.
;
; draw range ring 1000 km beyond anticyclone edge
;
        range_ring,y0,x0,dist0,360,bear,lats,lons
        print,'th ',rtheta2,' CM ',x0,y0,' d ',dist0
        oplot,[x0,x0],[y0,y0],psym=8,color=col1(k),symsize=3
        oplot,lons,lats,color=col1(k),thick=3
;
; extract ambient MLS ozone within dist0 of anticyclone and plot
;
        agood=where(o3mark eq 0. and dist le dist0 and ydata ge 20.)
        hgood=where(o3mark eq mmin and ydata ge 30.)
;
; save ambient and anticyclone ozone and theta
;
        if hgood(0) ne -1L then begin
if rtheta2 eq th(ilev) then oplot,xdata(hgood),ydata(hgood),psym=1,color=0
           ho3_save=[ho3_save,o3data(hgood)]
           hth_save=[hth_save,o3theta(hgood)]
           hx_save=[hx_save,xdata(hgood)]
           hy_save=[hy_save,ydata(hgood)]
           if n_elements(hgood) ge 3L then begin
              result=moment(o3data(hgood))
              ho3mean_save(k)=result(0)
              ho3sigma_save(k)=sqrt(result(1))
              result=moment(o3theta(hgood))
              hthmean_save(k)=result(0)
              hthsigma_save(k)=sqrt(result(1))
           endif
        endif
        if agood(0) ne -1L then begin
if rtheta2 eq th(ilev) then oplot,xdata(agood),ydata(agood),psym=1,color=70
           ao3_save=[ao3_save,o3data(agood)]
           ath_save=[ath_save,o3theta(agood)]
           ax_save=[ax_save,xdata(agood)]
           ay_save=[ay_save,ydata(agood)]
           if n_elements(agood) ge 3L then begin
              result=moment(o3data(agood))
              ao3mean_save(k)=result(0)
              ao3sigma_save(k)=sqrt(result(1))
              result=moment(o3theta(agood))
              athmean_save(k)=result(0)
              athsigma_save(k)=sqrt(result(1))
           endif
        endif
        jumplev:
      endfor    ; end loop over theta
;
; save mean ozone profiles on this day
;
      ho3mean_zt(icount,*)=ho3mean_save
      ao3mean_zt(icount,*)=ao3mean_save
      ho3sigma_zt(icount,*)=ho3sigma_save
      ao3sigma_zt(icount,*)=ao3sigma_save
      hthmean_zt(icount,*)=hthmean_save
      athmean_zt(icount,*)=athmean_save
      hthsigma_zt(icount,*)=hthsigma_save
      athsigma_zt(icount,*)=athsigma_save

      thmin=600. & thmax=1600.
      set_viewport,.3,.7,.1,.4
      index=where(ho3mean_save gt 0.)
      if index(0) ne -1L then begin
         plot,ho3mean_save(index),th(index),color=0,/nodata,thick=5,psym=0,$
              xrange=[2.,8.],yrange=[thmin,thmax],charsize=1.5,ytitle='Theta'
         oplot,ho3mean_save(index),th(index),color=mcolor*.9,thick=5,psym=0
         oplot,ho3mean_save(index)-ho3sigma_save(index),th(index),color=mcolor*.9,thick=5,psym=0,linestyle=5
         oplot,ho3mean_save(index)+ho3sigma_save(index),th(index),color=mcolor*.9,thick=5,psym=0,linestyle=5
         index=where(ao3mean_save gt 0.)
         if index(0) ne -1L then begin
            oplot,ao3mean_save(index),th(index),color=0,thick=5,psym=0
            oplot,ao3mean_save(index)-ao3sigma_save(index),th(index),color=0,thick=5,psym=0,linestyle=5
            oplot,ao3mean_save(index)+ao3sigma_save(index),th(index),color=0,thick=5,psym=0,linestyle=5
         endif
      endif
;
; loop over theta again
; mark anticyclone points with ozone less than ambient-1.5 sigma and 1 ppmv lower than ambient mean
; lop arrays
;
      tlop_save=[-99.]
      xlop_save=[-99.]
      ylop_save=[-99.]
      o3lop_save=[-99.]
      thlop_save=[-99.]
      o3lopmean_save=fltarr(nth)
      nlvls=27
      col1=reverse(1+indgen(nlvls)*mcolor/(float(nlvls)))
      for k=26L,15L,-1L do begin		; 500 K to 2000 K
          rtheta2=th(k)
          index=where(rtheta2 eq th)
          itheta=index(0)
          zindex=where(athmean_save eq rtheta2)
          if zindex(0) eq -1L then goto,jumplev2
          zindex2=where(hthmean_save eq rtheta2)
          if zindex2(0) eq -1L then goto,jumplev2
;
; require ao3mean_save > ho3mean_save by at least 0.25 ppmv
;
          if ao3mean_save(zindex(0))-ho3mean_save(zindex2(0)) lt 0.25 then goto,jumplev2
;
; 1.5 sigma threshold
;
          ao3thresh=ao3mean_save(zindex(0))-1.*ao3sigma_save(zindex(0))
          index=where(hth_save eq rtheta2 and ho3_save lt ao3thresh,npt)
          if index(0) eq -1L then goto,jumplev2
          oplot,ho3_save(index),hth_save(index),psym=8,color=mcolor*.3,symsize=1.5
;
; save LOP location and ozone
;
          xlop_save=[xlop_save,hx_save(index)]
          ylop_save=[ylop_save,hy_save(index)]
          thlop_save=[thlop_save,hth_save(index)]
          o3lop_save=[o3lop_save,ho3_save(index)]
;
; mean LOP profiles
;
          o3lopmean_save(k)=total(ho3_save(index))/float(npt)
          oplot,[o3lopmean_save(k),o3lopmean_save(k)],[rtheta2,rtheta2],psym=8,color=75,symsize=3
;
; retain LOP altitude
;
          lop_zt(icount,k)=1.
          jumplev2:
      endfor
;
; match dimension to ho3mean_save
;
      for k=0L,nth-1L do begin
          index=where(hthmean_save eq th(k))
          if index(0) eq -1L then o3lopmean_save(k)=-99
      endfor
      good=where(o3lopmean_save ne -99.)
      if good(0) ne -1L then o3lopmean_save=o3lopmean_save(good)
;
; remove first point
;
      index=where(xlop_save ne -99.,nlop)
      if index(0) ne -1L then begin
         tlop_save=tlop_save(index)
         xlop_save=xlop_save(index)
         ylop_save=ylop_save(index)
         thlop_save=thlop_save(index)
         o3lop_save=o3lop_save(index)
;
; remove low-level false positives
;
         if max(thlop_save) lt 900. then begin
            tlop_save=[-99.]
            xlop_save=[-99.]
            ylop_save=[-99.]
            o3lop_save=[-99.]
            thlop_save=[-99.]
            nlop=0L
         endif
      endif
      xyouts,2.2,620,strcompress(nlop,/remove_all)+' LOP Points',color=0,/data,charsize=1.5
;
; overplot LOP locations on first polar plot
;
;   index=where(thlop_save ne -99.,klop)
    index=where(thlop_save eq th(ilev),klop)
    if index(0) ne -1L then begin
    xmn=0.1
    xmx=0.5
    ymn=0.5
    ymx=0.9
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3
    MAP_SET,90,-90,0,/stereo,/noeras,/grid,/contin,charsize=1.5,color=0,/noborder
    for ii=0L,klop-1L do begin
        oplot,[xlop_save(index(ii)),xlop_save(index(ii))],$
              [ylop_save(index(ii)),ylop_save(index(ii))],psym=8,symsize=3,color=75
;             color=((o3lop_save(index(ii)-o3min)/(o3max-o3min))*mcolor,psym=8,symsize=2
;             color=((thlop_save(index(ii)-thmin)/(thmax-thmin))*mcolor,psym=8,symsize=2
    endfor
    endif
;
; save marked LOP points along with all anticyclone and nearby ambient t,x,y,th,ozone/theta points
; and mean/sigma profiles binned by theta
;
      save,file='/aura7/harvey/WACCM_data/Datfiles/wa3_nhDLlop_'+syr+smn+sdy+'.sav',xlop_save,ylop_save,$
           thlop_save,o3lop_save,o3lopmean_save,th,hx_save,hy_save,hth_save,ho3_save,ho3mean_save,$
           ho3sigma_save,hthmean_save,hthsigma_save,ax_save,ay_save,ath_save,ao3_save,ao3mean_save,$
           ao3sigma_save,athmean_save,athsigma_save,hareamean_save

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim Figures/lop_o3_wa3_nhDL_'+ifile+'.ps -rotate -90 Figures/lop_o3_wa3_nhDL_'+ifile+'.jpg'
    endif
    icount=icount+1L
goto, jump

end
