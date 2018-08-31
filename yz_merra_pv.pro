;
; yz MERRA PV
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

sver='v3.3'

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15]
yorig=[.2]
xlen=0.8
ylen=0.6
cbaryoff=0.04
cbarydel=0.01
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
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
mon=['jan','feb','mar','apr','may','jun',$
     'jul','aug','sep','oct','nov','dec']
month=['January','February','March','April','May','June',$
       'July','August','September','October','November','December']
stimes=[$
;'_0000.V01.',$
;'_0600.V01.',$
;'_1200.V01.',$
'_1800.V01.']
slabs=['00Z','06Z','12Z','18Z']
slabs=['18Z']
ntimes=n_elements(stimes)
!noeras=1
dirm='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
lstmn=3L & lstdy=1L & lstyr=2005L
ledmn=6L & leddy=30L & ledyr=2005L
lstday=0L & ledday=0L
;
; get date range
;
print, ' '
print, '      MERRA Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
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
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal termination condition '
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
print,sdate
;
; read MLS data
;
      dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) eq '' then begin
         mprof=0L
         goto,jumpmerra
      endif
      restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
      restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
      restore,dirm+'n2o_mls_'+sver+'_'+sdate+'.sav'              ; mix
      nz=n_elements(altitude)
      nthlev=n_elements(thlev)
      mprof=n_elements(longitude)
      mlev=n_elements(altitude)
      muttime=time
      mlat=latitude
      mlon=longitude
      bad=where(mask eq -99.)
      if bad(0) ne -1L then mix(bad)=-99.
      good=where(mix ne -99.)
;     if good(0) eq -1L then goto,jumpmerra
      mco=mix
      restore,dirm+'h2o_mls_'+sver+'_'+sdate+'.sav'              ; water vapor mix
      bad=where(mask eq -99.)
      if bad(0) ne -1L then mix(bad)=-99.
      good=where(mix ne -99.)
;     if good(0) eq -1L then goto,jumpmerra
      mh2o=mix
      mtemp=temperature
      mpress=pressure
;
; eliminate bad uttimes and SH
;
      index=where(muttime gt 0. and mlat gt 10.)
      if index(0) eq -1L then goto,jumpmerra
      muttime=reform(muttime(index))
      mlat=reform(mlat(index))
      mlon=reform(mlon(index))
      mtemp=reform(mtemp(index,*))
      mpress=reform(mpress(index,*))
      mco=reform(mco(index,*))
      mh2o=reform(mh2o(index,*))
      mtheta=mtemp*(1000./mpress)^0.286
      index=where(mtemp lt 0.)
      if index(0) ne -1L then mtheta(index)=-99.
;
; construct 2d MLS arrays to match N2O
;
      mpress2=mpress
      mtime2=0.*mco
      mlat2=0.*mco
      mlon2=0.*mco
      for i=0L,mlev-1L do begin
          mtime2(*,i)=muttime
          mlat2(*,i)=mlat
          mlon2(*,i)=mlon
      endfor
;
; extract MLS near rlev
;
      rlev=700.
      kindex=where(abs(mtheta-rlev) le 25. and mco ne -99. and mh2o ne -99.,mprof)
      if kindex(0) eq -1L then goto,jumpmerra
      codata=mco(kindex)*1.e7
      h2odata=mh2o(kindex)*1.e6
      xdata=mlon2(kindex)
      ydata=mlat2(kindex)
;
; loop over daily output times
;
;     for itime=0L,ntimes-1L do begin
;
; read MERRA data
;
jumpmerra:
    ncfile0=dir+sdate+'.nc3'
    rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
       u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
      if iflag eq 1 then goto,jump
    ncid=ncdf_open(dir+sdate+'.nc4')
    marknew2=fltarr(nr,nc,nth)
    ncdf_varget,ncid,3,marknew2
    ncdf_close,ncid

      x2d=fltarr(nr,nth)
      for k=0,nr-1 do x2d(k,*)=th
;
; select theta level to plot
;
      if icount eq 0 then begin
         print,th
;        read,' Enter desired theta surface ',rlev
         zindex=where(th eq rlev)
         ilev=zindex(0)
         slev=strcompress(long(th(ilev)),/remove_all)+'K'
      endif
      speed2=sqrt(u2^2.+v2^2.)
pvzm=mean(pv2,dim=2)

      mark1=transpose(mark2(*,*,ilev))
      marknew1=transpose(marknew2(*,*,ilev))
      sf1=transpose(sf2(*,*,ilev))
      pv1=transpose(pv2(*,*,ilev))
      p1=transpose(p2(*,*,ilev))
      q1=transpose(q2(*,*,ilev))
      u1=transpose(u2(*,*,ilev))
      speed1=transpose(speed2(*,*,ilev))
      mark=fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1
      mark(nc,*)=mark(0,*)
      marknew=fltarr(nc+1,nr)
      marknew(0:nc-1,0:nr-1)=marknew1
      marknew(nc,*)=marknew(0,*)
      sf=fltarr(nc+1,nr)
      sf(0:nc-1,0:nr-1)=sf1
      sf(nc,*)=sf(0,*)
      pv=fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1
      pv(nc,*)=pv(0,*)
      p=fltarr(nc+1,nr)
      p(0:nc-1,0:nr-1)=p1
      p(nc,*)=p(0,*)
      q=fltarr(nc+1,nr)
      q(0:nc-1,0:nr-1)=q1
      q(nc,*)=q(0,*)
      u=fltarr(nc+1,nr)
      u(0:nc-1,0:nr-1)=u1
      u(nc,*)=u(0,*)
      speed=fltarr(nc+1,nr)
      speed(0:nc-1,0:nr-1)=speed1
      speed(nc,*)=speed(0,*)
      temp=th(ilev)*((p/1000.)^(.286))
;
; save postscript version
;
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename='yz_merra_pv_'+sdate+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.thick=2.0                   ;Plotted lines twice as thick
      endif
;
; PV
;
      erase
      !type=2^2+2^3
      xyouts,.35,.05,sdate+' '+slev,/normal,color=0,charsize=2,charthick=2
      !type=2^2+2^3
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      index=where(x2d le 500.)
;      index=where(abs(pvzm) lt 10000)
      smin=min(pvzm(index))
      smax=max(pvzm(index))
      nlvls=20L
      if icount eq 0L then pvlevel=smin+((smax-smin)/float(nlvls))*indgen(nlvls+1)
      col1=1+indgen(nlvls)*mcolor/nlvls
      contour,pvzm,alat,th,levels=pvlevel,c_color=col1,/cell_fill,/noeras,xrange=[-90,90],yrange=[350,500],color=0,ytitle='Theta (K)',xtitle='Latitude'
      contour,pvzm,alat,th,/overplot,levels=pvlevel,color=mcolor,/noeras,/follow
      set_viewport,xmx+cbaryoff,xmx+cbaryoff+cbarydel,ymn,ymx
      !type=2^2+2^3+2^5
      omin=min(pvlevel)
      omax=max(pvlevel)
      plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax],title='(PVU)',color=0
      xbox=[0,10,10,0,0]
      y1=omin
      dy=(omax-omin)/float(nlvls)
      for j=0,nlvls-1 do begin
          ybox=[y1,y1,y1+dy,y1+dy,y1]
          polyfill,xbox,ybox,color=col1(j)
          y1=y1+dy
      endfor
;
      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device,/close
         spawn,'convert -trim yz_merra_pv_'+sdate+'.ps -rotate -90 yz_merra_pv_'+sdate+'.jpg'
      endif
;     endfor	; loop over 4 daily times
      icount=icount+1L
goto,jump
end
