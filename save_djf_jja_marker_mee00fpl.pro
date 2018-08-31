;
; save DJF and JJA 40-year average marker for Cora
; 6/3/2012
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
cbaryoff=0.065
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
!noeras=1
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_WACCM4/mee00fpl_FW2.cam2.h3.dyns.'
lstmn=1L & lstdy=1L & lstyr=2003L
ledmn=12L & leddy=31L & ledyr=2042L
lstday=0L & ledday=0L
;
; get date range
;
print, ' '
print, '      WACCM Version '
print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
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
      if ndays gt ledday then goto,saveit
;
; construct date string
;
      syr=strcompress(iyr,/remove_all)
      smn=string(FORMAT='(i2.2)',imn)
      sdy=string(FORMAT='(i2.2)',idy)
      sdate=syr+smn+sdy
;
; only read DJF and JJA
;
      if smn ne '12' and smn ne '01' and smn ne '02' and smn ne '06' and smn ne '07' and smn ne '08' then goto,jump
;
; read daily data
;
    ncfile0=dir+sdate+'_3D_dyn.nc3'
    dum=findfile(ncfile0)
    if dum(0) eq '' then goto,jump
    print,ncfile0
    ncid=ncdf_open(ncfile0)
    result0=ncdf_inquire(ncid)
    for idim=0,result0.ndims-1 do begin
        ncdf_diminq,ncid,idim,name,dim
        if name eq 'number_of_latitudes' then nr=dim
        if name eq 'number_of_longitudes' then nc=dim
        if name eq 'number_of_levels' then nth=dim
;       print,'read ',name,' dimension ',dim
    endfor
    for ivar=0,result0.nvars-1 do begin
        result=ncdf_varinq(ncid,ivar)
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if result.name eq 'latitude' then alat=data
        if result.name eq 'longitude' then alon=data
        if result.name eq 'theta' then th=data
;       if result.name eq 'IPV' then pv2=data
        if result.name eq 'P' then p2=data
;       if result.name eq 'U' then u2=data
;       if result.name eq 'V' then v2=data
;       if result.name eq 'QDF' then qdf2=data
;       if result.name eq 'Q' then q2=data
;       if result.name eq 'GPH' then gph2=data
;       if result.name eq 'SF' then sf2=data
        if result.name eq 'MARK' then mark2=data
;       print,ivar,result.name,min(data),max(data)
    endfor
    ncdf_close,ncid
;
; first time declare mean arrays
;
    if icount eq 0L then begin
       djf_mark=0.*mark2
       jja_mark=0.*mark2
       jja_pressure=0.*p2
       djf_pressure=0.*p2
       ndjf_mark=0.*mark2
       njja_mark=0.*mark2
       icount=1L
    endif
;
; DJF
;
    if smn eq '12' or smn eq '01' or smn eq '02' then begin
       djf_mark=djf_mark+mark2
       djf_pressure=djf_pressure+p2
       ndjf_mark=ndjf_mark+1.
    endif
;
; JJA
;
    if smn eq '06' or smn eq '07' or smn eq '08' then begin
       jja_mark=jja_mark+mark2
       jja_pressure=jja_pressure+p2
       njja_mark=njja_mark+1.
    endif
    goto,jump
;
; check
;
;   if smn ne '12' and smn ne '01' and smn ne '02' then goto,jump
;
;    x=fltarr(nc+1)
;    x(0:nc-1)=alon(0:nc-1)
;    x(nc)=alon(0)+360.
;
;; select theta levels to plot
;;   if l eq 0 then begin
;       zindex=where(th ge 500. and th le 30000.,nth2)
;       thlevs=strcompress(string(fix(th(zindex))))
;;   endif
;
;    erase
;    set_viewport,.1,.9,.1,.9
;    !type=2^2+2^3     ; suppress x and y axes
;    MAP_SET,90,0,0,/stereo,/contin,/grid,/noeras,londel=90.,$
;            label=1,lonlab=1,charsize=2,latdel=180.,color=0
;    for kk=0,nth2-1 do begin
;        lev=zindex(nth2-1-kk)
;        mark1=transpose(djf_mark(*,*,lev))/transpose(ndjf_mark(*,*,lev))
;        sf1=transpose(sf2(*,*,lev))
;        pv1=transpose(pv2(*,*,lev))
;        p1=transpose(p2(*,*,lev))
;        mpv1=pv1*((th(lev)/300.))^(-9./2.)
;
;; temperature
;        temp1=transpose(t2(*,*,lev))
;        temp=fltarr(nc+1,nr)
;        temp(0:nc-1,0:nr-1)=temp1(0:nc-1,0:nr-1)    ; NH
;        temp(nc,*)=temp(0,*)
;
;; height of theta surface
;        zth=transpose(z2(*,*,lev))
;        zth=reform(zth(*,nr/2:nr-1))	; NH
;	result=moment(zth)
;	avgz=result(0)
;        savgz=strcompress(string(fix(avgz)),/remove_all)
;;
;        dum=fltarr(nc+1,nr)
;        dum(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)    ; NH
;        dum(nc,*)=dum(0,*)
;        contour,dum,x,alat,levels=[0.1],color=(float(kk)/nth2)*mcolor,/overplot,$
;                c_labels=0,thick=10,max_value=1.e15
;      endfor
;    !type=2^2+2^3+2^5
;nlvls=nth2-1
;level=th(zindex)
;col1=1+(findgen(nlvls)/nlvls)*mcolor
;imin=min(level)
;imax=max(level)
;    set_viewport,.975,.985,.1,.9
;    plot,[0,0],[imin,imax],xrange=[0,10],$
;          yrange=[imin,imax],title='Theta (K)       ',/noeras,$
;          ytickname=thlevs,ytickv=level,yticks=nlvls-1,charsize=1,color=0,charthick=2
;    xbox=[0,10,10,0,0]
;    y1=imin
;    dy=(imax-imin)/float(nlvls)
;    for j=0,nlvls-1 do begin
;      ybox=[y1,y1,y1+dy,y1+dy,y1]
;      polyfill,xbox,ybox,color=col1(j)
;      y1=y1+dy
;    endfor
;goto, jump
;
; save IDL save file
;
saveit:
djf_mark=djf_mark/ndjf_mark
djf_pressure=djf_pressure/ndjf_mark
jja_mark=jja_mark/njja_mark
jja_pressure=jja_pressure/njja_mark
djf_press=fltarr(nth)
jja_press=fltarr(nth)
for kk=0,nth-1 do djf_press(kk)=mean(djf_pressure(nr/2:nr-1,*,kk))  ; NH
for kk=0,nth-1 do jja_press(kk)=mean(jja_pressure(0:nr/2-1,*,kk))  ; SH
save,file='mee00fpl_40-year_mark_djf_jja.sav',nc,nr,nth,alon,alat,th,djf_mark,jja_mark,djf_press,jja_press
end
