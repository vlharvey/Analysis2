;
; run circumpolar anticyclone algorithm
;
@stddat
@kgmt
@ckday
@kdate
@relvort
@marker_lows_v7
@marker_highs_v6
@marker_circhighs_v7
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
lstmn=1 & lstdy=1 & lstyr=95 & lstday=0
ledmn=1 & leddy=1 & ledyr=95 & ledday=0
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
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      sdate=syr+smn+sdy
;
;***Read WACCM data
      dum=findfile(dirh+sdate+'.nc3')
      if dum(0) ne '' then begin
         ifile=dirh+sdate+'.nc3'
         ofile=dirh+sdate+'.nc4'
      endif
      print,ifile
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
      mark2=fltarr(nr,nc,nth)
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
      ncdf_varget,ncid,9,mark2
      ncdf_varget,ncid,10,vp2
      ncdf_varget,ncid,11,sf2
      ncdf_close,ncid

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
;     for thlev=10,20 do begin
          print,th(thlev)
          theta=th(thlev)
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          zeta1=transpose(zeta2(*,*,thlev))
          mark1=transpose(mark2(*,*,thlev))
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
          mark=fltarr(nc+1,nr)
          mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
          mark(nc,*)=mark(0,*)
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
;stop
  
; streamfunction based polar vortex marker
          markl=mark
;         marker_lows_v7,sf,markl,qdf,zeta,u,v,x,xlat,theta
  
; streamfunction based anticyclone marker
          markh=mark
;         marker_highs_v6,sf,markh,qdf,zeta,u,v,x,xlat

; like polar vortex marker routine modified to find circumpolar highs
          markh2=0.*qdf
          marker_circhighs_v7,sf,markh2,qdf,zeta,u,v,x,xlat,theta
;
; check for overlap between circumpolar high and non circumpolar highs
;
          lindex=where(markh2 lt 0. and markh lt 0.)
          if lindex(0) ne -1 then begin
;print,'highs overlapping'
;erase
;!type=2^2+2^3
;contour,sf,x,xlat,nlevels=20,/noeras,title=string(th(thlev))
;index=where(markh lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=8,color=mcolor
;index=where(markh2 lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=mcolor*.3,symsize=0.5
;stop
             if min(y2d(lindex)) gt 0. then begin
                kindex=where(markh lt 0. and y2d gt 0.,kpt)
                jindex=where(markh2 lt 0. and y2d gt 0.,jpt)
                if max(y2d(jindex)) eq max(y2d) then markh(lindex)=0.		; keep circumpolar
                if max(y2d(jindex)) lt max(y2d) then markh2(lindex)=0.
;               if kpt gt jpt then markh2(lindex)=0.		; keep bigger
;               if kpt lt jpt then markh(lindex)=0.
             endif
             if max(y2d(lindex)) lt 0. then begin
                kindex=where(markh lt 0. and y2d lt 0.,kpt)
                jindex=where(markh2 lt 0. and y2d lt 0.,jpt)
                if min(y2d(jindex)) eq min(y2d) then markh(lindex)=0.           ; keep circumpolar
                if min(y2d(jindex)) gt min(y2d) then markh2(lindex)=0.
;               if kpt gt jpt then markh2(lindex)=0.		; keep bigger
;               if kpt lt jpt then markh(lindex)=0.
             endif
;print,'are the highs still overlapping ?'
;erase
;!type=2^2+2^3
;contour,sf,x,xlat,nlevels=20,/noeras,title=string(th(thlev))
;index=where(markh lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=8,color=mcolor
;index=where(markh2 lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=mcolor*.3,symsize=0.5
;stop
          endif
;
; merge 2 anticyclone marker fields
;
         index=where(markh2 lt 0.)
         if index(0) ne -1L then markh(index)=min(markh)-1.0
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
;stop
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markh1=0.*qdf1
          markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)

      ENDFOR	; loop over theta

; Write isentropic data in netCDF format
      write_waccm3_nc4,ofile,nc,nr,nth,xlon,xlat,th,mark2

      goto,jump
end
