;
; No plotting- just store.
;
; Calculate GSFC forecast surface and mean sea level pressure.
; Read in 2.5x2 topograpy data and use geopotential height
; to calculate the pressure at the station height.  MSLP
; is computed hydrostatically.
;
; GSFC fields from 1000 hPa up to 500 hPa are discontinuous and 
; contain special values of 1.e+15 where the pressure surface is
; below the ground.  Therefore, cannot hydrostatically compute MSLP
; where there is no Z and T at the ground.  To compute Psfc find
; lowest good data.  Where this data is above the topography height
; set Psfc = MSLP.

dir='/usr71/users/ukmo/'
nam1='    '
nam2='    '
nam3='    '
gup=0.
glw=0.
gtp=0.
prs=[1000.00,850.000,700.000,500.000,400.000,300.000,250.000,$
     200.000,150.000,100.000,70.0000,50.0000,30.0000,10.0000,$
     5.00000,2.00000,1.00000,0.400000]

nc=0L
nr=0L
nl=0L
nfld=0
fld=' '
ifile='                            '
qfile='                                 '
close,1
openr,1,'ukmo_prs2theta_nc.fil'
readf,1,nfile
for nn=0,nfile-1 do begin
    readf,1,ifile
    close,2
    openr,2,dir+ifile,/f77
    print,'opening ',ifile
    readu,2,nfld
    for m=0,nfld-1 do begin
        readu,2,fld
        readu,2,nc,nr,nl
        lon=fltarr(nc)
        lat=fltarr(nr)
        plev=fltarr(nl)
        readu,2,lon,lat,plev
        data=fltarr(nc,nr,nl)
        readu,2,data
        if m eq 0 then t3d=data
        if m eq 1 then u3d=data
        if m eq 2 then v3d=data
        if m eq 3 then g3d=data
    endfor
   index=where(g3d lt 0.)
   g3d(index)=1.e15
   q3d=0.*t3d
   if nn eq 0 then begin
; read global topography (geopotential meters)
; topg3.75 matches daily UKMO z,t format 90- -90 and 0-360
topg=fltarr(nc,nr)
close,10
openr,10,'topg3.75',/f77
readu,10,topg
close,10
nlvls=21
icolmax=!p.color
; convert meters to geopotential meters
       z=topg/9.806
       z2=fltarr(nc+1,nr)
       z2(0:nc-1,0:nr-1)=z(0:nc-1,0:nr-1)
       z2(nc,*)=z2(0,*)
       zlevel=min(z)+((max(z)-min(z))/nlvls)*findgen(nlvls)
       lon2=fltarr(nc+1)
       lon2(0:nc-1)=lon(0:nc-1)
       lon2(nc)=lon2(0)+360.
       col1=1+indgen(nlvls)*icolmax/nlvls
    endif

    nlg=nc
    nlat=nr

; calculate MSLP (in mb) from Z at 1000mb (100000Pa)
; g1000 and t1000 have some 1.e+15 values
; set mslp at these values=1000mb
    g1000=fltarr(nlg,nlat)
    t1000=fltarr(nlg,nlat)
    mslp=fltarr(nlg,nlat)
    g1000(0:nlg-1,0:nlat-1)=g3d(0:nlg-1,0:nlat-1,0)
    t1000(0:nlg-1,0:nlat-1)=t3d(0:nlg-1,0:nlat-1,0)
    mslp=(100000.+((100000.*g1000*9.806)/(287.*t1000)))/100.
    index=where(g1000 eq 1.e15) 
    mslp(index)=1000.00

; surface arrays do not include ipvsfc,qdfsfc
    psfc=fltarr(nlg,nlat)
    tsfc=fltarr(nlg,nlat)
    usfc=fltarr(nlg,nlat)
    vsfc=fltarr(nlg,nlat)
    qsfc=fltarr(nlg,nlat)
    thsfc=fltarr(nlg,nlat)
    msfsfc=fltarr(nlg,nlat)

; calculate theta and msf
    p3d=0.*t3d
    for k=0,nl-1 do begin
        p3d(*,*,k)=plev(k)
    endfor
    th3d=t3d*(1000./p3d)^.286
    msf3d=1004.*t3d+9.806*g3d
    index=where(g3d eq 1.e15)
    t3d(index)=1.e15
    th3d(index)=1.e15
    msf3d(index)=1.e15

; obtain Ps by finding GSFC Z bounding topg Z and interpolating pressure
    for i=0,nlg-1 do begin
        for j=0,nlat-1 do begin
            gtp=z(i,j)
            for k=0,n_elements(prs)-2 do begin
                gup=g3d(i,j,k+1)
                glw=g3d(i,j,k)
                scale=(gup-gtp)/(gup-glw)
                if gup ne 1.e15 then begin			; lowest good data
                if gup ge gtp AND glw le gtp then begin		; Zs bounding topg
                   psfc(i,j)=prs(k+1)-scale*(prs(k+1)-prs(k))
                   tsfc(i,j)=t3d(i,j,k+1)-scale*(t3d(i,j,k+1)-t3d(i,j,k))
                   usfc(i,j)=u3d(i,j,k+1)-scale*(u3d(i,j,k+1)-u3d(i,j,k))
                   vsfc(i,j)=v3d(i,j,k+1)-scale*(v3d(i,j,k+1)-v3d(i,j,k))
                   thsfc(i,j)=th3d(i,j,k+1)-scale*(th3d(i,j,k+1)-th3d(i,j,k))
                   msfsfc(i,j)=msf3d(i,j,k+1)-scale*(msf3d(i,j,k+1)-msf3d(i,j,k))
                   qsfc(i,j)=q3d(i,j,k+1)-scale*(q3d(i,j,k+1)-q3d(i,j,k))
; diabatic heating array has more special values than analyses
                   if q3d(i,j,k) eq 1.e15 or q3d(i,j,k+1) eq 1.e15 then begin
                      index=where(q3d(i,j,*) ne 1.e15)
                      if index(0) ne -1 then qsfc(i,j)=q3d(i,j,index(0))
                   endif
                   goto,jumpout
                endif
                if glw eq 1.e15 then begin  ; none or 1 level w/gup ge gtp
                   psfc(i,j)=prs(k+1)		; sfc first good val
                   tsfc(i,j)=t3d(i,j,k+1)
                   usfc(i,j)=u3d(i,j,k+1)
                   vsfc(i,j)=v3d(i,j,k+1)
                   qsfc(i,j)=q3d(i,j,k+1)
                   thsfc(i,j)=th3d(i,j,k+1)
                   msfsfc(i,j)=msf3d(i,j,k+1)
; diabatic heating array has more special values than analyses
                   if q3d(i,j,k+1) eq 1.e15 then begin
                      index=where(q3d(i,j,*) ne 1.e15)  
                      if index(0) ne -1 then qsfc(i,j)=q3d(i,j,index(0))
                   endif
                   goto,jumpout
                endif
                endif
                if glw ne 1.e15 and glw gt gtp then begin ; glw exceeds gtp
                   psfc(i,j)=mslp(i,j)
                   tsfc(i,j)=t3d(i,j,k)
                   usfc(i,j)=u3d(i,j,k)
                   vsfc(i,j)=v3d(i,j,k)
                   qsfc(i,j)=q3d(i,j,k)
                   thsfc(i,j)=th3d(i,j,k)
                   msfsfc(i,j)=msf3d(i,j,k)
; diabatic heating array has more special values than analyses
                   if q3d(i,j,k) eq 1.e15 then begin
                      index=where(q3d(i,j,*) ne 1.e15)  
                      if index(0) ne -1 then qsfc(i,j)=q3d(i,j,index(0))
                   endif
                   goto,jumpout
                endif
                jumplev:
            endfor
            jumpout:
;           if psfc(i,j) eq 0. then print,gup,gtp,glw,mslp(i,j),i,j
         endfor
     endfor

; store surface quantities
     close,6
     openw,6,dir+ifile+'.psfc',/f77
     writeu,6,nc,nr
     writeu,6,lon,lat
     writeu,6,mslp,psfc,tsfc,usfc,vsfc,qsfc,thsfc,msfsfc
     close,6
endfor
end
