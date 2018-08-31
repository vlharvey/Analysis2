;
; reads in .nc2 files.  marks polar vortex and anticyclones
; and writes out .nc3 files.
; v7 polar vortices
;
@stddat
@kgmt
@ckday
@kdate
@rd_mls_nc2
@relvort
@marker_lows_v7
@marker_lows_pv_v8
@marker_highs_v6
@write_mls_nc3

loadct,39
mcolor=byte(!p.color)
nlvls=30L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.1]
yorig=[0.1]
xlen=0.8
ylen=0.8
setplot='ps'
;read,'setplot?',setplot
if setplot ne 'ps' then begin
   lc=0
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dirw='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/'
ifiles=file_search(dirw+'MLS_grid_theta_*.nc2',count=nfile)
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN
    result=strsplit(ifiles(n),'_',/extract)
    result2=strsplit(result(-1),'.',/extract)
    sdate=result2(0)
    print,sdate
;
; skip if nc3 file already exists
;
    dum=findfile(dirw+'MLS_grid_theta_'+sdate+'.nc3')
    if dum(0) ne '' then goto,jumpfile
;
;***Read data
;
      ifile=dirw+'MLS_grid_theta_'+sdate+'.nc2'
      print,'reading ',ifile
      rd_mls_nc2,ifile,nc,nr,nth,xlon,xlat,th,pv2,p2,u2,v2,qdf2,co2,z2,sf2,h2o2,iflg
      if iflg ne 0L then goto,jumpfile
      x=fltarr(nc+1)
      x(0:nc-1)=xlon(0:nc-1)
      x(nc)=xlon(0)+360.
      mark2=0.*qdf2
      markco2=0.*qdf2
;
; meridional CO gradient calculation
;
co2=smooth(co2,3,/edge_truncate,/Nan)
cograd2=0*co2
cograd4=0*co2
alat=xlat
alon=xlon
dy2=2.*(alat(1)-alat(0))*!pi/180.
dy4=12.*(alat(1)-alat(0))*!pi/180.
for k=0,nth-1 do begin
for j=2,nr-3 do begin
    jm1=j-1
    jp1=j+1
    jm2=j-2
    jp2=j+2
    for i=0,nc-1 do begin
        cograd2(j,i,k) = (co2(jp1,i,k)-co2(jm1,i,k))/dy2        ; 2nd order
        cograd4(j,i,k) = (-1.*co2(jp2,i,k)+8.*co2(jp1,i,k) $
                          -8.*co2(jm1,i,k)+1.*co2(jm2,i,k))/dy4 ; 4th order
    endfor
endfor
cograd2(0,*,k)=cograd2(2,*,k)
cograd2(1,*,k)=cograd2(2,*,k)
cograd4(nr-2,*,k)=cograd2(nr-3,*,k)
cograd4(nr-1,*,k)=cograd2(nr-3,*,k)
endfor
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
      for thlev=0,nth-1 do begin	;nth-1 do begin	; 5000 K and below
          theta=th(thlev)
          pv1=transpose(pv2(*,*,thlev))
          u1=transpose(u2(*,*,thlev))
          v1=transpose(v2(*,*,thlev))
          qdf1=transpose(qdf2(*,*,thlev))
          sf1=transpose(sf2(*,*,thlev))
          z1=transpose(z2(*,*,thlev))
          co1=transpose(co2(*,*,thlev))
          cograd1=transpose(cograd2(*,*,thlev))
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
speed=sqrt(u^2.0+v^2.0)
          z=fltarr(nc+1,nr)
          z(0:nc-1,0:nr-1)=z1(0:nc-1,0:nr-1)
          z(nc,*)=z(0,*)
          qdf=fltarr(nc+1,nr)
          qdf(0:nc-1,0:nr-1)=qdf1(0:nc-1,0:nr-1)
          qdf(nc,*)=qdf(0,*)
          sf=0.*fltarr(nc+1,nr)
          sf(0:nc-1,0:nr-1)=sf1(0:nc-1,0:nr-1)
          sf(nc,*)=sf(0,*)
          co=0.*fltarr(nc+1,nr)
          co(0:nc-1,0:nr-1)=co1(0:nc-1,0:nr-1)
          co(nc,*)=co(0,*)
co=co*1.e6
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
;
; pv based polar vortex marker
          markcol=0.*qdf
          marker_lows_pv_v8,co,markcol,qdf,zeta,u,v,x,xlat,theta
;
; co gradient marker
; calculate horizontal co gradient
          markcogradl=0.*qdf
;         marker_lows_pv_v9,co,markcogradl,cograd1,zeta,u,v,x,xlat,theta,sfval
  
; streamfunction based anticyclone marker
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
      !type=2^2+3^2
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
index=where(z1 gt 0.)
if index(0) ne -1L then meanz=mean(z1(index))
erase
!type=2^2+2^3
;MAP_SET,0,0,0,/noeras,/grid,/contin,title=sdate+' '+string(th(thlev))+' K ('+strcompress(round(meanz/1000.),/remove_all)+' km)',charsize=2.0
index=where(finite(co) eq 1)
if index(0) ne -1L then if max(co(index)) gt 0. then begin
   imin=min(co(index))
   imax=max(co(index))
   if imin eq imax then imin=imax/2.
   level=imin+((imax-imin)/float(nlvls))*findgen(nlvls)
   contour,co,x,xlat,levels=level,/cell_fill,c_color=col1,/noeras,title=sdate+' '+string(th(thlev))+' K ('+strcompress(round(meanz/1000.),/remove_all)+' km)',charsize=2.0
   z=z/1000.
   contour,z,x,xlat,levels=min(z)+((max(z)-min(z))/20.)*findgen(20),/noeras,thick=3,/overplot,color=mcolor
endif
loadct,0
contour,markl,x,xlat,levels=[0.1],/noeras,thick=10,/overplot,color=0
contour,markcol,x,xlat,levels=[0.1],/noeras,thick=10,/overplot,color=mcolor*.9
;contour,markcogradl,x,xlat,levels=[0.1],/noeras,thick=10,/overplot,color=mcolor*.5
loadct,39
;contour,speed,x,xlat,levels=[50],/noeras,thick=3,/overplot,color=mcolor*.3
;contour,speed,x,xlat,levels=[70],/noeras,thick=3,/overplot,color=mcolor*.4
;contour,speed,x,xlat,levels=[90],/noeras,thick=3,/overplot,color=mcolor*.5
;contour,speed,x,xlat,levels=[110],/noeras,thick=3,/overplot,color=mcolor*.6
;contour,speed,x,xlat,levels=[130],/noeras,thick=3,/overplot,color=mcolor*.7
;contour,speed,x,xlat,levels=[150],/noeras,thick=3,/overplot,color=mcolor*.8
;contour,speed,x,xlat,levels=[170],/noeras,thick=3,/overplot,color=mcolor*.9
index=where(finite(co) eq 1)
if index(0) ne -1L and max(co(index)) gt 0. then wait,.1
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markl(0:nc-1,0:nr-1)
          mark2(*,*,thlev)=transpose(markl1)
          markl1=0.*qdf1
          markl1(0:nc-1,0:nr-1)=markcol(0:nc-1,0:nr-1)
          markco2(*,*,thlev)=transpose(markl1)

;         markh1=0.*qdf1
;         markh1(0:nc-1,0:nr-1)=markh(0:nc-1,0:nr-1)
;         mark2(*,*,thlev)=mark2(*,*,thlev)+transpose(markh1)

      ENDFOR	; loop over theta

; Write isentropic data in netCDF format
      ofile=dirw+'MLS_grid_theta_'+sdate+'.nc3'
      print,'writing ',ofile
      write_mls_nc3,ofile,nc,nr,nth,xlon,xlat,th,pv2,p2,u2,v2,qdf2,mark2,co2,z2,sf2,h2o2,markco2
jumpfile:
endfor
end
