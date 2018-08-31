;
; reads in .nc2 files.  marks polar vortex and anticyclones
; and writes out .nc3 files.
; v7 polar vortices
;
; this version uses the poleward-most jet to define the vortex
; instead of the strongest
;
@stddat
@kgmt
@ckday
@kdate
@compvort
@calcelat2d
@marker_lows_v8
@marker_highs_v6
@write_waccm3_nc4

loadct,38
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dirh='/aura7/harvey/WACCM_data/Datfiles/Aurora/waccm3_'
lstmn=1 & lstdy=1 & lstyr=6 & lstday=0
ledmn=1 & leddy=1 & ledyr=6 & ledday=0
;
; Ask interactive questions- get starting/ending date
;
print, ' '
print, '      WACCM Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 50 then lstyr=lstyr+2000
if ledyr lt 50 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
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
      syr=strtrim(string(iyr),2)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
;
;***Read WACCM data
      dum=findfile(dirh+sdate+'.nc3')
      if dum(0) ne '' then begin
         ifile=dirh+sdate+'.nc3'
         ofile=dirh+sdate+'.nc5'
      endif
      ncid=ncdf_open(ifile)
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      alon=fltarr(nc)
      alat=fltarr(nr)
      th=fltarr(nth)
      pv2=fltarr(nr,nc,nth)
      p2=fltarr(nr,nc,nth)
      msf2=fltarr(nr,nc,nth)
      u2=fltarr(nr,nc,nth)
      v2=fltarr(nr,nc,nth)
      qdf2=fltarr(nr,nc,nth)
      mark2=fltarr(nr,nc,nth)
      vp2=fltarr(nr,nc,nth)
      sf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,alon
      ncdf_varget,ncid,1,alat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,pv2
      ncdf_varget,ncid,4,p2
      ncdf_varget,ncid,5,msf2
      ncdf_varget,ncid,6,u2
      ncdf_varget,ncid,7,v2
      ncdf_varget,ncid,8,qdf2
      ncdf_varget,ncid,9,mark2
      ncdf_varget,ncid,10,vp2
      ncdf_varget,ncid,11,sf2
      ncdf_close,ncid
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.
      mark2=0.*qdf2
;
; compute 3d relative vorticity (this version of compvort requires dimensions lon,lat,lev
;
      zz=fltarr(nc,nr,nth)
      uu=fltarr(nc,nr,nth)
      vv=fltarr(nc,nr,nth)
      for k=0,nth-1 do begin
          uu(*,*,k)=transpose(u2(*,*,k))
          vv(*,*,k)=transpose(v2(*,*,k))
      endfor
      relvort,uu,vv,zz,alon,alat,nc,nr
      zeta2=0.*qdf2
      for k=0,nth-1 do zeta2(*,*,k)=transpose(zz(*,*,k))

; loop over theta
      mpvnh=999. & mpvsh=-999.
      for thlev=0,nth-1 do begin
          theta=th(thlev)
          pv1=transpose(pv2(*,*,thlev))
;         elat1=calcelat2d(pv1,alon,alat)
          p1=transpose(p2(*,*,thlev))
          mpv1=pv1*((th(thlev)/300.))^(-9./2.)
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          zeta1=transpose(zeta2(*,*,thlev))
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
          qdf=fltarr(nc+1,nr)
          qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
          qdf(nc,*)=qdf(0,*)
          sf=0.*fltarr(nc+1,nr)
          sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
          sf(nc,*)=sf(0,*)
          x2d=0.*sf
          y2d=0.*sf
          for i=0,nc do y2d(i,*)=alat
          for j=0,nr-1 do x2d(*,j)=x
  
; streamfunction based polar vortex marker
          markl=0.*qdf
          marker_lows_v8,sf,markl,qdf,zeta,u,v,x,alat,theta,pv
;
; sub-vortex modification
;
        if theta lt 200. then begin
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
;         marker_highs_v6,sf,markh,qdf,zeta,u,v,x,alat,pv
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
contour,sf,x,alat,nlevels=20,/noeras,xrange=[0.,360.],yrange=[-90.,90],$
        title=sdate+' '+string(theta)+' K',charsize=2
map_set,0,180,0,/contin,/grid,/noeras
index=where(markl gt 0.)
if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=.2*mcolor
index=where(markh lt 0.)
if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=.9*mcolor

          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markh1=0.*qdf1
          markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)
jumplev:

      ENDFOR	; loop over theta

; Write UKMO isentropic data in netCDF format
      write_waccm3_nc4,ofile,nc,nr,nth,alon,alat,th,mark2

      goto,jump
end
