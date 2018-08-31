;
; save daytime zonal means
; store multi-year monthly mean WACCM3 in all fields 
;
@rd_waccm_nc3

loadct,38
mcolor=!p.color
mcolor=byte(!p.color)
device,decompose=0
dir='/aura7/harvey/WACCM_data/Datfiles/wa3_tnv3_'
spawn,'ls '+dir+'*nc3',allfiles
;
; extract month numbers from filenames
;
wmon=strarr(n_elements(allfiles))
for n=0L,n_elements(allfiles)-1L do begin
    result=strsplit(allfiles(n),'_',/extract)
    wmon(n)=strmid(result(3),4,2)
endfor
spawn,'/usr/bin/rm '+dir+'*-*'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
nmonths=n_elements(mon)
for m=0,nmonths-1 do begin
;
; pull out days in this month from all years
;
    smon=string(format='(i2.2)',m+1L)
    index=where(wmon eq smon,nfile)
    ifiles=allfiles(index)
    kfile=0L
    for n=0L,nfile-1 do begin
        ifile=ifiles(n)
        print,ifile
;
; read data
;
        rd_waccm3_nc3,ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
             u2,v2,qdf2,mark2,sf2,o3,ch4,no2,h2o,iflag

        if n eq 0 then no2_mean=fltarr(nr,nth)
;
; compute solar zenith angle at each gridpoint
;
doy=
      pi=3.14159265
      dtor=pi/180.
      earinc=23.5
      zangle=fltarr(nr,nc)
      for jj=0L,nr-1 do begin
      for ii=0L,nc-1 do begin
          rlat=alat(jj)
          rlon=alon(ii)
          gmt=12.
          sinlat=sin(rlat*dtor)
          coslat=sqrt(1.-sinlat^2.)
          sinlon=sin(rlon*dtor)
          coslon=cos(rlon*dtor)
          soya=(doy-81.25)*pi/182.5           ; day angle
          soha=2.*pi*(gmt-12.)/24.            ; hour angle
          soha=-soha
          sininc=sin(earinc*dtor)
          sindec=sininc*sin(soya)
          cosdec= sqrt(1.-sindec^2.)
          coszen=cos(soha)*coslon+sin(soha)*sinlon
          coszen=coszen*cosdec*coslat
          coszen=sindec*sinlat+coszen
          coszen=min([max([coszen,-1.]),1.])
          chi = acos(coszen)
          zangle(jj,ii) = chi/dtor
      endfor
      endfor
;
; summate zonal mean daytime data
;
      for k=0,nth-1L do begin
          no2lev=reform(no2lev(*,*,k))
          for j=0,nr-1L do begin
              zanglelat=reform(zangle(jj,*))
              index=where(zanglelat lt 80.,nn)
              no2_mean(j,k)=no2_mean(j,k)+total(no2lev(index,j))/float(nn)
          endfor
      endfor

      kfile=kfile+1L
    endfor          ; loop over files
;
; divide by total number of days summed
;
    no2_mean=no2_mean/float(kfile)

; write monthly means 
    save,file=dir+mon(m)+'avg.sav',nc,nr,nth,alon,alat,th,no2_mean

    thlev=20
    plt1=transpose(reform(no2_mean(*,*,thlev),nr,nc))
    plt=fltarr(nc+1,nr)
    plt(0:nc-1,0:nr-1)=plt1
    plt(nc,*)=plt(0,*)
    x=fltarr(nc+1)
    x(0:nc-1)=alon
    x(nc)=x(0)
    !type=2^2+2^3
    nlvls=30
    erase
    col1=3+indgen(nlvls)*mcolor/nlvls
    MAP_SET,0,0,0,/contin,/grid,/noeras,title=mon(m)
    contour,plt,x,alat,nlevels=30,/fill,/cell_fill,/overplot,c_color=col1,/noeras
    contour,plt,x,alat,nlevels=30,/follow,/overplot,c_color=0,/noeras,c_labels=0*indgen(30)
    MAP_SET,0,0,0,/contin,/noeras
endfor  ; loop over months
end

