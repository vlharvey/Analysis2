;
; plot hovmoller of temperature
;
@rd_ukmo_nc3
loadct,38
mcolor=byte(fix(!p.color))
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=800
nydim=800
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.12
cbarydel=0.02
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   lc=icolmax
   set_plot,'x'
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
!noeras=1
dir='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
ifiles=[$
'ukmo_files_92.fil',$
'ukmo_files_93.fil',$
'ukmo_files_94.fil',$
'ukmo_files_95.fil',$
'ukmo_files_96.fil',$
'ukmo_files_97.fil',$
'ukmo_files_98.fil',$
'ukmo_files_99.fil',$
'ukmo_files_00.fil',$
'ukmo_files_01.fil',$
'ukmo_files_02.fil',$
'ukmo_files_03.fil',$
'ukmo_files_04.fil',$
'ukmo_files_05.fil']
nyear=n_elements(ifiles)
for n=0,nyear-1L do begin 
    close,1
    openr,1,'Fil_files/'+ifiles(n)
    nfile=0L
    readf,1,nfile
    sfile=strarr(nfile)
    for l=0,nfile-1 do begin
        readf,1,ifile
        print,ifile
        sfile(l)=ifile
        dum1=findfile(dir+ifile+'.nc3')
        if dum1(0) eq '' then goto,jump
        ncid=ncdf_open(dir+ifile+'.nc3')
        ncdf_diminq,ncid,0,name,nr
        ncdf_diminq,ncid,1,name,nc
        ncdf_diminq,ncid,2,name,nth
        alon=fltarr(nc)
        alat=fltarr(nr)
        th=fltarr(nth)
        p2=fltarr(nr,nc,nth)
        ncdf_varget,ncid,0,alon
        ncdf_varget,ncid,1,alat
        ncdf_varget,ncid,2,th
        ncdf_varget,ncid,4,p2
        ncdf_close,ncid
        t2=0.*p2
        for k=0,nth-1 do t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
        if l eq 0L then xt=fltarr(nc,nfile)
;
; user prompted input at very beginning
;
        if n eq 0L and l eq 0L then begin
           rlat=0.
           print,alat
           read,'Enter desired latitude ',rlat
           index=where(alat eq rlat)
           ilat=index(0)
           ralt=0.
           print,th
           read,'Enter desired altitude ',ralt
           index=where(th eq ralt)
           ialt=index(0)
        endif
;
; save hovmoller
;
        for i=0,nc-1 do xt(*,l)=t2(ilat,*,ialt)
jump:
    endfor
;
; plot Hov
;
    sth=strcompress(string(fix(ralt)),/remove_all)
    slat=strcompress(string(rlat),/remove_all)
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='hov_T_'+sth+'K_'+slat+'_'+ifile+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    !type=2^2+2^3
    set_viewport,xmn,xmx,ymn,ymx
    if rlat gt 0. then mtitle='MetO Temperature  '+sth+' K '+slat+' N'
    if rlat lt 0. then mtitle='MetO Temperature  '+sth+' K '+slat+' S'
    if n eq 0L then begin
       imin=fix(min(xt))
       imax=round(max(xt))
imin=200.
imax=220.
    endif
    nlev=11
    col1=1+(findgen(nlev)/nlev)*mcolor
    cint=(imax-imin)/(nlev-1)
    level=imin+cint*findgen(nlev)
    index=where(strmid(sfile,4,2) eq '01',nytick)
    contour,xt,alon,findgen(nfile),xrange=[0.,360.],/fill,$
            /cell_fill,yrange=[0,nfile-1L],xstyle=1,ystyle=1,xticks=6,$
            yticks=nytick-1,ytickv=index,ytickname=sfile(index),$
            xtitle='Longitude',ytitle='Time',c_color=col1,$
            title=mtitle,/noeras,levels=level,charsize=2
    contour,xt,alon,findgen(nfile),/overplot,/follow,levels=level,$
            c_color=0,c_labels=0*level,/noeras
    contour,xt,alon,findgen(nfile),/overplot,/follow,levels=0,$
            c_color=0,c_labels=0,thick=3,/noeras
    ymnb=yorig-cbaryoff
    ymxb=ymnb +cbarydel
    set_viewport,xmn,xmx,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],charsize=1.5
    ybox=[0,10,10,0,0]
    x1=imin
    dx=(imax-imin)/float(nlev)
    for jj=0,nlev-1 do begin
        xbox=[x1,x1,x1+dx,x1+dx,x1]
        polyfill,xbox,ybox,color=col1(jj)
        x1=x1+dx
    endfor
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim hov_T_'+sth+'K_'+slat+'_'+ifile+'.ps -rotate -90 hov_T_'+sth+'K_'+slat+'_'+ifile+'.jpg'
    endif
    if setplot ne 'ps' then stop
endfor
end
