;
; Read in daily .nc3 files and write out monthly averages.
; Programmed by VLH Wed Mar  5 14:36:17 EST 2003
;
@rd_ukmo_nc3
@write_ukmo_nc3

dir='/usr72/users/ukmo/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'store_monthly_aves.fil'
;
; loop over months
;
nmonth=0L
readf,1,nmonth
for m=0,nmonth-1 do begin
;
; loop over days
;
    nfile=0L
    readf,1,nfile
    for n=0,nfile-1 do begin
;
; read UKMO isentropic data
;
        readf,1,ifile
        rd_ukmo_nc3,dir+ifile+'.nc3',nc,nr,nth,xlon,xlat,th,$
                   pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
;
; declare monthly average arrays
;
        if n eq 0 then begin
           pvave=0.*pv2
           pave=0.*pv2
           msfave=0.*pv2
           uave=0.*pv2
           vave=0.*pv2
           qave=0.*pv2
           qdfave=0.*pv2
           markave=0.*pv2
           vpave=0.*pv2
           sfave=0.*pv2
        endif
;
; compute monthly averages
;
        pvave=pvave+pv2
        pave=pave+p2
        msfave=msfave+msf2
        uave=uave+u2
        vave=vave+v2
        qave=qave+q2
        qdfave=qdfave+qdf2
        markave=markave+mark2
        vpave=vpave+vp2
        sfave=sfave+sf2

    endfor		; loop over days

    pvave=pvave/float(nfile)
    pave=pave/float(nfile)
    msfave=msfave/float(nfile)
    uave=uave/float(nfile)
    vave=vave/float(nfile)
    qave=qave/float(nfile)
    qdfave=qdfave/float(nfile)
    markave=markave/float(nfile)
    vpave=vpave/float(nfile)
    sfave=sfave/float(nfile)
;
; write monthly averages
;
    write_ukmo_nc3,dir+ifile+'_ave.nc3',nc,nr,nth,xlon,xlat,th,$
          pvave,pave,msfave,uave,vave,qave,qdfave,markave,vpave,sfave

endfor			; loop over months
end
