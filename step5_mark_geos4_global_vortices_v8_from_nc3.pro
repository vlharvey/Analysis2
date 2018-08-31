;
; reads in .nc3 files.  marks new polar vortex
; and writes out .nc5 files.
; v8 polar vortices
;
; this version adds SAO and STJ logic.  poleward most jet above 2000 and below 500 K
;
; this version adds the "subvortex modification" below 500 K
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3
@compvort
@calcelat2d
@marker_lows_v8
@marker_highs_v6
@write_ukmo_nc4
@drawvectors

loadct,38
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.15]
xlen=0.7
ylen=0.7
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dir='/aura7/harvey/GEOS4_data/Datfiles/'
lstmn=3 & lstdy=31 & lstyr=4 & lstday=0
ledmn=5 & leddy=15 & ledyr=7 & ledday=0
lstmn=9 & lstdy=1 & lstyr=4 & lstday=0
ledmn=9 & leddy=15 & ledyr=4 & ledday=0
;
; Ask interactive questions- get starting/ending date
;
;print, ' '
;print, '      GEOS-4 Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 2000 then lstyr=lstyr+2000
if ledyr lt 2000 then ledyr=ledyr+2000
if lstyr lt 2004 then stop,'Year out of range '
if ledyr lt 2004 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
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
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition'
;
; read UKMO data
;
      syr=strtrim(string(iyr),2)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      date=syr+smn+sdy
      ifile=dir+'DAS.flk.asm.tavg3d_mis_e.GEOS403.MetO.'+date
      print,ifile
      rd_geos5_nc3,ifile+'_1200.V01.nc3',nc,nr,nth,xlon,xlat,th,$
                 pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
mark2=0.*mark2
      x=fltarr(nc+1)
      x(0:nc-1)=xlon(0:nc-1)
      x(nc)=xlon(0)+360.
      mark2=0.*qdf2

; loop over theta
      mpvnh=999. & mpvsh=-999.
;     for thlev=0,nth-1 do begin
      for thlev=0,13 do begin		; only down to 400 K to save time
          theta=th(thlev)
          pv1=transpose(pv2(*,*,thlev))
;         elat1=calcelat2d(pv1,xlon,xlat)
          p1=transpose(p2(*,*,thlev))
          mpv1=pv1*((th(thlev)/300.))^(-9./2.)
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          zeta1=u1*0.0
          compvort,u1,v1,zeta1,xlon,xlat,nc,nr
          mpv=fltarr(nc+1,nr)
          mpv(0:nc-1,0:nr-1)=mpv1(0:nc-1,0:nr-1)
          mpv(nc,*)=mpv(0,*)
          pv=fltarr(nc+1,nr)
          pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
          pv(nc,*)=pv(0,*)
;         elat=fltarr(nc+1,nr)
;         elat(0:nc-1,0:nr-1)=elat1(0:nc-1,0:nr-1)
;         elat(nc,*)=elat(0,*)
          p=fltarr(nc+1,nr)
          p(0:nc-1,0:nr-1)=p1(0:nc-1,0:nr-1)
          p(nc,*)=p(0,*)
          tp=theta*(p/1000.)^0.286
          zeta=fltarr(nc+1,nr)
          zeta(0:nc-1,0:nr-1)=zeta1(0:nc-1,0:nr-1)
          zeta(nc,*)=zeta(0,*)
          u=fltarr(nc+1,nr)
          u(0:nc-1,0:nr-1)=u1(0:nc-1,0:nr-1)
          u(nc,*)=u(0,*)
          v=fltarr(nc+1,nr)
          v(0:nc-1,0:nr-1)=v1(0:nc-1,0:nr-1)
          v(nc,*)=v(0,*)
          sp=sqrt(u^2.+v^2.)
          qdf=fltarr(nc+1,nr)
          qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
          qdf(nc,*)=qdf(0,*)
          sf=0.*fltarr(nc+1,nr)
          sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
          sf(nc,*)=sf(0,*)
          x2d=0.*sf
          y2d=0.*sf
          for i=0,nc do y2d(i,*)=xlat
          for j=0,nr-1 do x2d(*,j)=x
  
; streamfunction based polar vortex marker
          markl=0.*qdf
          marker_lows_v8,sf,markl,qdf,zeta,u,v,x,xlat,theta,pv
;
; sub-vortex modification
;
        if theta eq 500. then begin
           nhindex=where(markl gt 0. and y2d gt 0.)
           shindex=where(markl gt 0. and y2d lt 0.)
           if nhindex(0) ne -1 then begin
              mpvnh=total(mpv(nhindex))/n_elements(nhindex)
;             elatnh=min(elat(nhindex))
           endif
           if shindex(0) ne -1 then begin
              mpvsh=total(mpv(shindex))/n_elements(shindex)
;             elatsh=max(elat(shindex))
           endif
           print,'AVG MPV ',mpvnh,mpvsh
;          print,'ELAT EDGE ',elatnh,elatsh
        endif
        if theta lt 500. then begin
           index=where(y2d gt 20. and mpv gt mpvnh and tp lt 230.)
;          if index(0) ne -1 then markl(index)=1.
index=where(y2d gt 20. and tp gt 240.)
if index(0) ne -1 then markl(index)=0.
;index=where(y2d gt 20. and mpv lt mpvnh)
;if index(0) ne -1 then markl(index)=0.

           index=where(y2d lt -20. and mpv lt mpvsh and tp lt 230.)
;          if index(0) ne -1 then markl(index)=1.
index=where(y2d lt -20. and tp gt 240.)
if index(0) ne -1 then markl(index)=0.
;index=where(y2d lt -20. and mpv gt mpvsh)
;if index(0) ne -1 then markl(index)=0.
        endif
;  
; streamfunction based anticyclone marker
;
          markh=0.*qdf
;         marker_highs_v6,sf,markh,qdf,zeta,u,v,x,xlat,pv
;
; check for vortex and anticyclone overlap
;
          lindex=where(markl gt 0. and markh lt 0.)
          if lindex(0) ne -1 then begin
;            if min(y2d(lindex)) lt 0. and max(y2d(lindex)) gt 0. then stop
;
; NH
; 
             if min(y2d(lindex)) gt 0. then begin
                s0=min(sf(lindex))
                kindex=where(markl gt 0. and y2d gt 0.)
                s1=min(sf(kindex))
                index=where(sf ge (s0+s1)/2.0 and y2d gt 0.)
                markl(index)=0.
             endif
;
; SH
; 
             if max(y2d(lindex)) lt 0. then begin
                s0=max(sf(lindex))
                kindex=where(markl gt 0. and y2d lt 0.)
                s1=max(sf(kindex))
                index=where(sf le (s0+s1)/2.0 and y2d lt 0.)
                markl(index)=0.
             endif
          endif
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
index=where(pv ne 1.e12)
if index(0) eq -1L then goto,jumplev
imin=min(pv)
imax=max(pv(index))
level=imin+((imax-imin)/nlvls)*findgen(nlvls)
contour,pv,x,xlat,levels=level,/fill,c_color=col1,/noeras,title=date+'  '+string(theta)+' K',charsize=2
contour,pv,x,xlat,levels=level,/follow,/overplot,color=icolmax,thick=2,max_value=1.e12
drawvectors,nc+1,nr,x,xlat,u,v,4,0
loadct,0
contour,markl,x,xlat,levels=[0.1],/follow,/overplot,thick=5,color=100
loadct,38
;index=where(markl gt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=.2*mcolor
;index=where(markh lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=.9*mcolor
;if max(markl) gt 0. then stop
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markh1=0.*qdf1
          markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)
jumplev:

      ENDFOR	; loop over theta

; Write UKMO isentropic data in netCDF format
      write_ukmo_nc4,ifile+'_1200.V01.nc5',nc,nr,nth,xlon,xlat,th,mark2
      goto,jump
end
