;
; plot vortex max wind speed as a function of altitude and time
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo_nc3

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
device,decompose=0
!noeras=1
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.35]
xlen=0.8
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
set_plot,'x'
setplot='x'
read,'setplot= ',setplot
if setplot ne 'ps' then $
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
smon=['J','F','M','A','M','J','J','A','S','O','N','D']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
;
; Ask interactive questions- get starting/ending date
;
;lstyr=1992L
lstyr=2005L
ledyr=2005L
;read,' Enter starting date (month, day, year) ',lstyr
;read,' Enter ending date   (month, day, year) ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
;
; loop over years
;
for kyr=lstyr,ledyr do begin
    lstmn=10 & lstdy=1
    ledmn=5 & leddy=1
    z = stddat(lstmn,lstdy,kyr,lstday)
    z = stddat(ledmn,leddy,kyr+1L,ledday)
    if ledday lt lstday then stop,' Wrong dates! '
    ndays=ledday-lstday+1L

    icount=0L
    idy = lstdy
    imn = lstmn
    iyr = kyr
    z = kgmt(imn,idy,iyr,iday)
    iday = iday - 1
;
; loop over days
;
    jump: iday = iday + 1
          kdate,float(iday),iyr,imn,idy
          ckday,iday,iyr

; --- Test for end condition and close windows.
          z = stddat(imn,idy,iyr,kdays)
          if kdays lt lstday then stop,' starting day outside range '
          if kdays gt ledday then goto,plotit

          syr=strtrim(string(iyr),2)
          uyr=strmid(syr,2,2)
          smn=string(FORMAT='(i2.2)',imn)
          sdy=string(FORMAT='(i2.2)',idy)

          ifile=mon(imn-1)+sdy+'_'+uyr
          dum1=findfile(diru+ifile+'.nc3')
          if dum1(0) ne '' then ncid=ncdf_open(diru+ifile+'.nc3')
          if dum1(0) eq '' then goto,jumpday
          ncdf_diminq,ncid,0,name,nr
          ncdf_diminq,ncid,1,name,nc
          ncdf_diminq,ncid,2,name,nth
          alon=fltarr(nc)
          alat=fltarr(nr)
          th=fltarr(nth)
          u2=fltarr(nr,nc,nth)
          v2=fltarr(nr,nc,nth)
          ncdf_varget,ncid,0,alon
          ncdf_varget,ncid,1,alat
          ncdf_varget,ncid,2,th
          ncdf_varget,ncid,6,u2
          ncdf_varget,ncid,7,v2
          ncdf_close,ncid
          s2=sqrt(u2^2.0+v2^2.0)
          print,ifile,max(s2)
          if icount eq 0L then begin
             wzt=fltarr(ndays,nth)
             latzt=fltarr(ndays,nth)
             sfile=strarr(ndays)
             dum=reform(u2(*,*,0))
             y2d=0.*dum
             for i=0,nc-1 do y2d(*,i)=alat
          endif

          for thlev=0L,nth-1L do begin
              s1=reform(s2(*,*,thlev))
              index=where(y2d gt 30.)
              wzt(icount,thlev)=max(s1(index))
;
; at what latitude did max(s1(index)) occur?
;
              index=where(y2d le 30.)
              s1(index)=0.
              index=where(s1 eq max(s1))
              latzt(icount,thlev)=y2d(index(0))
          endfor
          jumpday:
          sfile(icount)='y'+uyr+'_m'+smn+'_d'+sdy
          icount=icount+1L
    goto,jump
;
; plot
;
plotit:
    yy=strmid(sfile(0),1,2)
    if long(yy) lt 90L then y1='20'+yy
    if long(yy) gt 90L then y1='19'+yy
    yy=strmid(sfile(ndays-1),1,2)
    if long(yy) lt 90L then y2='20'+yy
    if long(yy) gt 90L then y2='19'+yy
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='ukmo_wmax_'+y1+'-'+y2+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3+2^7
    plot,[1,ndays,ndays,1,1],[500.,500.,2000.,2000.,500.],min_value=0.,$
          xrange=[1,ndays],yrange=[500.,2000.],/nodata,charsize=2,$
          ytitle='Theta (K)',title='MetO Analyses '+y1+'-'+y2,xtickname=[' ',' '],xticks=1
    kindex=where(strmid(sfile,9,2) eq '15',nxtick)
    xmon=long(strmid(sfile(kindex),5,2))
    for i=0,nxtick-1 do begin
        xlab=smon(xmon(i)-1)
        plots,kindex(i)+1,450.
        plots,kindex(i)+1,500.,/continue,/data
        xyouts,kindex(i)+1,325.,xlab,/data,alignment=0.5,charsize=3
    endfor
    nlvls=21
    level=10.*findgen(nlvls)
    col1=1+indgen(nlvls)*icolmax/nlvls
    contour,wzt,1.+findgen(ndays),th,levels=level,/fill,$
            /cell_fill,/overplot,c_color=col1,min_value=0.
    contour,wzt,1.+findgen(ndays),th,levels=level,c_color=0,$
            /follow,/overplot,min_value=0.
    imin=min(level)
    imax=max(level)
    ymnb=yorig(0) -cbaryoff
    ymxb=ymnb  +cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
         xtitle='Maximum Vortex Wind Speed (m/s)',charsize=2
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlvls)
    for j=0,nlvls-1 do begin
    xbox=[x1,x1,x1+dx,x1+dx,x1]
    polyfill,xbox,ybox,color=col1(j)
    x1=x1+dx
    endfor
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim ukmo_wmax_'+y1+'-'+y2+'.ps -rotate -90 ukmo_wmax_'+y1+'-'+y2+'.jpg'
    endif
endfor	; loop over years
end
