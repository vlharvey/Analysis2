;
; interpolate the edge of the Arctic vortex from GEOS-4
; to each MLS CO swath

@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

sver='v2.2'
sver='v1.52'

loadct,38
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=1000
nydim=800
xorig=[.1,.3,.5,.7,.1,.3,.5,.7,.1,.3,.5,.7,.1,.3,.5,.7]
yorig=[.7,.7,.7,.7,.5,.5,.5,.5,.3,.3,.3,.3,.1,.1,.1,.1]
xlen=0.15
ylen=0.15
cbaryoff=0.085
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
!noeras=1
dirm='/aura6/data/MLS_data/Datfiles_SOSST/'
idir='/aura6/data/GEOS4_data/Analysis/Save_files/'
dir='/aura7/harvey/GEOS4_data/Datfiles/DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'
lstmn=5L & lstdy=1L & lstyr=2006L
ledmn=8L & leddy=31L & ledyr=2006L
lstday=0L & ledday=0L
;
; get date range
;
print, ' '
print, '      GEOS-4 Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
if lstyr lt 2004 then stop,'Year out of range '
if ledyr lt 2004 then stop,'Year out of range '
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
      rd_geos5_nc3_meto,dir+sdate+'_1200.V01.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.

; select theta levels to plot
      if icount eq 0 then begin
         index=where(th ge 500. and th le 5000.,nth2)
         th2=reverse(th(index))
         thlevs=strcompress(string(fix(th2)))+' K'
         thlw=min(th2)
         thup=max(th2)
      endif
;
; read MLS.  ozone and MetO mark are on altitude.  elat is on theta
;
      dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      dum=findfile(dirm+'dmps_mls_'+sver+'.meto.'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      dum=findfile(dirm+'tpd_mls_'+sver+'_'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      dum=findfile(dirm+'mark_mls_'+sver+'.meto.'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
      restore,dirm+'dmps_mls_'+sver+'.meto.'+sdate+'.sav'       ; elat_prof
      restore,dirm+'tpd_mls_'+sver+'_'+sdate+'.sav'             ; temperature, pressure
      restore,dirm+'co_mls_'+sver+'_'+sdate+'.sav'              ; mix
      restore,dirm+'mark_mls_'+sver+'.meto.'+sdate+'.sav'       ; mark_prof
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
      if good(0) eq -1L then goto,jump
      mco=mix*1.e6
      mtemp=temperature
      mpress=pressure
      metomark=MARK_PROF
;
; eliminate bad uttimes and NH
;
      index=where(muttime gt 0. and mlat lt -25.,mprof)
      if index(0) eq -1L then goto,jump
      muttime=reform(muttime(index))
      mlat=reform(mlat(index))
      mlon=reform(mlon(index))
      mtemp=reform(mtemp(index,*))
      mpress=reform(mpress(index,*))
      mco=reform(mco(index,*))
      metomark=reform(metomark(index,*))
      mtheta=mtemp*(1000./mpress)^0.286
      index=where(mtemp lt 0.)
      if index(0) ne -1L then mtheta(index)=-99.
;
; construct 2d MLS arrays to match ozone
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
; determine number of SH swaths and print time of swaths on polar plot
;
    index=where(mlat le -81.7,npt)
    tswath=muttime(index)
    flag=0.*tswath
    nswath=50L
    tsave=fltarr(nswath)
    mcount=0L
    for i=0L,npt-1L do begin
        index=where(abs(tswath(i)-tswath) lt 1. and flag eq 0.)
        if index(0) ne -1L then begin
           flag(index)=1.0
           kindex=where(abs(muttime-tswath(index(0))) le 0.5 and mlat lt 0.)
           stime=string(FORMAT='(I2.2)',long(tswath(index(0))))
           tsave(mcount)=tswath(index(0))
           mcount=mcount+1L
           xtmp=mlon(kindex)
           ytmp=mlat(kindex)
           index=where(ytmp eq min(ytmp))
        endif
    endfor
    nswath=mcount
    tsave=tsave(0:mcount-1L)
;
; check for marker file
;
;   dum=findfile(idir+'geos4_mark_at_mls_'+sdate+'_sh.sav')
;   if dum(0) ne '' then begin
;      restore,idir+'geos4_mark_at_mls_'+sdate+'_sh.sav'
;      goto,jumpint
;   endif
;
; interpolate GEOS-4 marker to 2d MLS array (mco)
; but retain field of integers
;
    gmark=0.*mco
    for ii=0L,mprof-1L do begin
        for kk=0L,mlev-1L do begin
            slon=mlon2(ii,kk)
            slat=mlat2(ii,kk)
            slev=mtheta(ii,kk)
            if slon lt alon(0) then slon=slon+360.
            for i=0L,nc-1L do begin
                ip1=i+1
                if i eq nc-1L then ip1=0L
                xlon=alon(i)
                xlonp1=alon(ip1)
                if i eq nc-1L then xlonp1=360.+alon(ip1)
                if slon ge xlon and slon le xlonp1 then begin
                   xscale=(slon-xlon)/(xlonp1-xlon)
                   goto,jumpx
                endif
            endfor
jumpx:
            for j=0L,nr-2L do begin
                jp1=j+1
                xlat=alat(j)
                xlatp1=alat(jp1)
                if slat ge xlat and slat le xlatp1 then begin
                    yscale=(slat-xlat)/(xlatp1-xlat)
                    goto,jumpy
                endif
            endfor
jumpy:
            for k=1L,nth-1L do begin
                kp1=k-1
                uz=th(k)
                uzp1=th(kp1)
                if slev ge uz and slev le uzp1 then begin
                   zscale=(slev-uz)/(uzp1-uz)
                   p1=mark2(j,i,k)
                   pp2=mark2(jp1,i,k)
                   p3=mark2(j,ip1,k)
                   p4=mark2(jp1,ip1,k)
                   p5=mark2(j,i,kp1)
                   p6=mark2(jp1,i,kp1)
                   p7=mark2(j,ip1,kp1)
                   p8=mark2(jp1,ip1,kp1)
                   if p1 eq 1. or pp2 eq 1. or p3 eq 1. or p4 eq 1. or $
                      p5 eq 1. or p6 eq 1. or p7 eq 1. or p8 eq 1. then $
                      gmark(ii,kk)=1.0
                   if p1 lt 0. or pp2 lt 0. or p3 lt 0. or p4 lt 0. or $
                      p5 lt 0. or p6 lt 0. or p7 lt 0. or p8 lt 0. then $
                      gmark(ii,kk)=min([p1,pp2,p3,p4,p5,p6,p7,p8])
                   goto,jumpz1
                endif
            endfor
jumpz1:
        endfor
    endfor
    save,file=idir+'geos4_mark_at_mls_'+sdate+'_sh.sav',gmark,altitude
jumpint:
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
       device,/landscape,bits=8,filename='sh_geos4_vortex+mls_coswaths_'+sdate+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
       !p.thick=2.0                   ;Plotted lines twice as thick
       !p.charsize=2.0
    endif
;
; plot 
;
    erase
    xyouts,0.2,0.9,'MLS CO + MetO and GEOS-4 Antarctic Vortex '+sdate,/normal,charsize=2.0,color=0,charthick=2
    !type=2^2+2^3
    for iswath=0L,nswath-1L do begin
        tplot=tsave(iswath)
        kindex=where(abs(muttime-tplot) le 0.5 and mlat lt 0.,mprofs)
        stime=strcompress(string(FORMAT='(F4.1)',muttime(kindex(0))),/remove_all)
        coswath=smooth(mco(kindex,*),3,/edge_truncate)
        index=where(coswath lt 0.)
        if index(0) ne -1L then coswath(index)=0.1
        copr=0.*coswath
        for k=0L,nz-1L do begin
            index=where(coswath(*,k) ne 0.1)
            if n_elements(index) ge 3L then begin
               result=moment(coswath(index,k))
               copr(index,k)=coswath(index,k)-result(0)
            endif
        endfor
        mmarkswath=metomark(kindex,*)
        gmarkswath=gmark(kindex,*)
        xswath=findgen(mprofs)
        ylabels=string(format='(f4.1)',mlat(kindex))
        xlabels=string(format='(f5.1)',mlon(kindex))

        xmn=xorig(iswath)
        xmx=xorig(iswath)+xlen
        ymn=yorig(iswath)
        ymx=yorig(iswath)+ylen
        set_viewport,xmn,xmx,ymn,ymx
        level=[0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,30.,40.,50.]
;       level=-10.+0.5*findgen(41)	; for copr
        nlvls=n_elements(level)
        imin=min(level) & imax=max(level)
        col1=1+indgen(nlvls)*icolmax/nlvls
        contour,coswath,xswath,altitude,title='t='+stime,levels=level,/fill,c_color=col1,charsize=1.5,$
                yrange=[10.,85.],/noeras,color=0
        contour,coswath,xswath,altitude,levels=[0.5,5.],color=0,/follow,/overplot,thick=2
        contour,mmarkswath,xswath,altitude,levels=[0.2],color=mcolor,/follow,/overplot,thick=6,c_linestyle=5,c_labels=0
;       contour,mmarkswath,xswath,altitude,levels=[-0.2],color=0,/follow,/overplot,thick=6,c_linestyle=5,c_labels=0
        contour,gmarkswath,xswath,altitude,levels=[0.2],color=mcolor,/follow,/overplot,thick=3,c_labels=0
;       contour,gmarkswath,xswath,altitude,levels=[-0.2],color=0,/follow,/overplot,thick=3,c_labels=0
    endfor	; loop over MLS swaths
    xmn=xorig(nswath)
    xmx=xorig(nswath)+xlen
    ymn=yorig(nswath)
    ymx=yorig(nswath)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    map_set,-90,0,-90,/stereo,/contin,/grid,color=0,/noeras
    oplot,mlon,mlat,psym=1,color=0,symsize=0.5
    marklev=reform(gmark(*,40))
    index=where(marklev ge 0.2)
    if index(0) ne -1L then oplot,mlon2(index),mlat2(index),psym=1,color=75,symsize=0.5

    set_viewport,.2,.78,.14-cbaryoff,.14-cbaryoff+cbarydel
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],$
          xrange=[imin,imax],xtitle='Carbon Monoxide (ppmv)',/noeras,$
          xtickname=strcompress(string(format='(f5.2)',level),/remove_all),$
          xstyle=1,charsize=1.5,color=0,charthick=2,xticks=nlvls-1
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
      xbox=[x1,x1,x1+dx,x1+dx,x1]
      polyfill,xbox,ybox,color=col1(j)
      x1=x1+dx
    endfor
    !p.charthick=1.
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim sh_geos4_vortex+mls_coswaths_'+sdate+'.ps -rotate -90 '+$
             'sh_geos4_vortex+mls_coswaths_'+sdate+'.jpg'
       spawn,'/usr/bin/rm -f sh_geos4_vortex+mls_coswaths_'+sdate+'.ps'
    endif
    icount=icount+1L
goto,jump
end
