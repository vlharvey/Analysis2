;
; reads in .nc2 files.  marks polar vortex and anticyclones
; and writes out .nc3 files.
; v7 polar vortices
;
@stddat
@kgmt
@ckday
@kdate
@relvort
@marker_lows_v7
@marker_highs_v6

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
lstmn=1 & lstdy=1 & lstyr=1995 & lstday=0
ledmn=1 & leddy=1 & ledyr=1995 & ledday=0
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
      if ndays gt ledday then stop,' Normal termination condition '
      sday=string(FORMAT='(i3.3)',iday-1)
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      sdate=syr+smn+sdy
;
;***Read WACCM data
      dum=findfile(dirh+sdate+'.nc2')
      if dum(0) ne '' then begin
         ifile=dirh+sdate+'.nc2'
         ofile=dirh+sdate+'.nc3'
      endif
      ncid=ncdf_open(ifile)
      ncdf_diminq,ncid,0,name,nr
      ncdf_diminq,ncid,1,name,nc
      ncdf_diminq,ncid,2,name,nth
      xlon=fltarr(nc)
      xlat=fltarr(nr)
      th=fltarr(nth)
      pv2=fltarr(nr,nc,nth)
      p2=fltarr(nr,nc,nth)
      msf2=fltarr(nr,nc,nth)
      u2=fltarr(nr,nc,nth)
      v2=fltarr(nr,nc,nth)
      qdf2=fltarr(nr,nc,nth)
      vp2=fltarr(nr,nc,nth)
      sf2=fltarr(nr,nc,nth)
      ncdf_varget,ncid,0,xlon
      ncdf_varget,ncid,1,xlat
      ncdf_varget,ncid,2,th
      ncdf_varget,ncid,3,pv2
      ncdf_varget,ncid,4,p2
      ncdf_varget,ncid,5,msf2
      ncdf_varget,ncid,6,u2
      ncdf_varget,ncid,7,v2
      ncdf_varget,ncid,8,qdf2
      ncdf_varget,ncid,9,vp2
      ncdf_varget,ncid,10,sf2
      ncdf_close,ncid
      print,'read ',ifile
      x=fltarr(nc+1)
      x(0:nc-1)=xlon(0:nc-1)
      x(nc)=xlon(0)+360.
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
      relvort,uu,vv,zz,xlon,xlat,nc,nr
      zeta2=0.*qdf2
      for k=0,nth-1 do zeta2(*,*,k)=transpose(zz(*,*,k))

; loop over theta
      for thlev=0,nth-1 do begin
          theta=th(thlev)
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          zeta1=transpose(zeta2(*,*,thlev))
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
          for i=0,nc do y2d(i,*)=xlat
          for j=0,nr-1 do x2d(*,j)=x
;erase
;!type=2^2+2^3
;imin=min(zeta) & imax=max(zeta)
;nlvls=30
;level=imin+((imax-imin)/nlvls)*findgen(nlvls)
;map_set,0,0,0,/contin,/grid,/noeras
;contour,zeta,x,xlat,levels=level,/overplot,/noeras,title=string(th(thlev)),c_linestyle=(level lt 0.)
  
; streamfunction based polar vortex marker
          markl=0.*qdf
          marker_lows_v7,sf,markl,qdf,zeta,u,v,x,xlat,theta
  
; streamfunction based anticyclone marker
          markh=0.*qdf
          marker_highs_v6,sf,markh,qdf,zeta,u,v,x,xlat
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
;erase
;!type=2^2+2^3
;contour,sf,x,xlat,nlevels=20,/noeras,title=string(th(thlev))
;index=where(qdf lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=1,color=mcolor,symsize=0.5
;index=where(zeta lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=mcolor*.3,symsize=0.5
;index=where(markl gt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=.1*mcolor
;index=where(markh lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=.5*mcolor
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markh1=0.*qdf1
          markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)

      ENDFOR	; loop over theta

; Write isentropic data in netCDF format
      print,'writing ',ofile
      ncid = ncdf_create(ofile,/CLOBBER)
      latdimid=ncdf_dimdef(ncid, 'number_of_latitudes' , nr)
      londimid=ncdf_dimdef(ncid, 'number_of_longitudes',  nc)
      levdimid=ncdf_dimdef(ncid, 'number_of_levels'    ,  nth)
      lonsid = ncdf_vardef(ncid, 'longitudes',  londimid)
      latsid = ncdf_vardef(ncid, 'latitudes' ,  latdimid)
      levsid = ncdf_vardef(ncid, 'th_levels' ,  levdimid)
      ipvid  = ncdf_vardef(ncid, 'ipv'       , [latdimid,londimid,levdimid])
      prsid  = ncdf_vardef(ncid, 'press'     , [latdimid,londimid,levdimid])
      msfid  = ncdf_vardef(ncid, 'msf'       , [latdimid,londimid,levdimid])
      uuuid  = ncdf_vardef(ncid, 'u_wind'    , [latdimid,londimid,levdimid])
      vvvid  = ncdf_vardef(ncid, 'v_wind'    , [latdimid,londimid,levdimid])
      qdfid  = ncdf_vardef(ncid, 'qdf'       , [latdimid,londimid,levdimid])
      mksid  = ncdf_vardef(ncid, 'mark'      , [latdimid,londimid,levdimid])
      vptid  = ncdf_vardef(ncid, 'vel_pot'   , [latdimid,londimid,levdimid])
      stfid  = ncdf_vardef(ncid, 'strm'      , [latdimid,londimid,levdimid])
      ncdf_control,ncid,/ENDEF
      ncdf_varput, ncid, lonsid, xlon , COUNT=[nc]
      ncdf_varput, ncid, latsid, xlat , COUNT=[nr]
      ncdf_varput, ncid, levsid, th   , COUNT=[nth]
      ncdf_varput, ncid, ipvid , pv2  , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, prsid , p2   , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, msfid , msf2 , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, uuuid , u2   , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, vvvid , v2   , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, qdfid , qdf2 , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, mksid , mark2, COUNT=[nr,nc,nth]
      ncdf_varput, ncid, vptid , vp2  , COUNT=[nr,nc,nth]
      ncdf_varput, ncid, stfid , sf2  , COUNT=[nr,nc,nth]
      ncdf_close,ncid
      goto,jump
end
