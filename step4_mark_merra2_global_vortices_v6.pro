;
; reads in .nc2 files.  marks polar vortex and anticyclones
; and writes out .nc3 files.
; v7 polar vortices
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc2
@relvort
@marker_lows_v7
@marker_highs_v6
@write_merra2_nc3

;loadct,39
;mcolor=byte(!p.color)
;device,decompose=0
;a=findgen(8)*(2*!pi/8.)
;usersym,cos(a),sin(a),/fill
;nxdim=700
;nydim=700
;xorig=[0.15]
;yorig=[0.25]
;xlen=0.7
;ylen=0.5
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dirw='/atmos/harvey/MERRA2_data/Datfiles/'
ifiles=file_search(dirw+'MERRA2-on-WACCM_theta_2018????*.nc2',count=nfile)
print,'starting'
spawn,'date'
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN
    result=strsplit(ifiles(n),'_',/extract)
    result2=strsplit(result(3),'.',/extract)
    sdate=result2(0)
;   print,sdate
;
; skip if nc3 file already exists
;
    dum=findfile(dirw+'MERRA2-on-WACCM_theta_'+sdate+'.nc3')
    if dum(0) ne '' then goto,jumpfile
;
;***Read data
;
      ifile=dirw+'MERRA2-on-WACCM_theta_'+sdate+'.nc2'
      print,'reading ',ifile
      rd_merra2_nc2,ifile,nc,nr,nth,xlon,xlat,th,pv2,p2,u2,v2,qdf2,qv2,z2,sf2,q2,o32,iflg
      if iflg ne 0L then goto,jumpfile
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
          pv1=transpose(pv2(*,*,thlev))
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          pv=fltarr(nc+1,nr)
          pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
          pv(nc,*)=pv(0,*)
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
;stop
  
; streamfunction based polar vortex marker
          markl=0.*qdf
          marker_lows_v7,sf,markl,qdf,zeta,u,v,x,xlat,theta
  
; streamfunction based anticyclone marker
          markh=0.*qdf
          marker_highs_v6,sf,markh,qdf,zeta,u,v,x,xlat,pv
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
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=.1*mcolor
;index=where(markh lt 0.)
;if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=.9*mcolor
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markh1=0.*qdf1
          markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)

      ENDFOR	; loop over theta

; Write isentropic data in netCDF format
      ofile=dirw+'MERRA2-on-WACCM_theta_'+sdate+'.nc3'
      print,'writing ',ofile
      write_merra2_nc3,ofile,nc,nr,nth,xlon,xlat,th,pv2,p2,u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32

jumpfile:
endfor
print,'ending'
spawn,'date'
end
