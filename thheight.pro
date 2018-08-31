;
; print UKMO vertical profiles of hemispheric averaged theta, height, pressure, temp, u, pvu

@rd_ukmo_nc3

nam1='    '
nam2='    '
nam3='    '
th0=0
tmin=0.0

nc=0L
nr=0L
nl=0L
nfld=0
fld=' '
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                                       '
close,1
openr,1,'thheight.fil'
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    if iflag eq 1 then goto,jump

; Height of isentropic surface = (msf - cp*T)/g
z2=0.*qdf2
t2=0.*qdf2
zavg=fltarr(nth)
pavg=fltarr(nth)
for k=0,nth-1 do begin
    t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
    z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
    res1=moment(z2(*,*,k))
    res2=moment(p2(*,*,k))
    zavg(k)=res1(0)
    pavg(k)=res2(0)
    if k eq 0 then print,fix(th(k)),res1(0),res2(0)
    if k gt 0 then print,fix(th(k)),res1(0),res2(0),$
       zavg(k-1)-zavg(k),pavg(k-1)-pavg(k)
endfor

print,'    Theta     Height      Pressure     Temperature    Uwind        PVU'
for k=0,nth-1 do begin
    t2(*,*,k) = th(k)*( (p2(*,*,k)/1000.)^(.286) )
    z2(*,*,k) = (msf2(*,*,k) - 1004.*t2(*,*,k))/(9.86*1000.)
    res1=moment(z2((nr/2):nr-1,*,k))
    res2=moment(p2((nr/2):nr-1,*,k))
    res3=moment(t2((nr/2):nr-1,*,k))
    res4=moment(u2((nr/2):nr-1,*,k))
    res5=moment(pv2((nr/2):nr-1,*,k))
    print,fix(th(k)),res1(0),res2(0),res3(0),res4(0),res5(0)*1.e4
endfor
;print,' '
;print,'    Theta   Height      Pressure   Temperature     Uwind          PVU'
;    for k=0,nth-1 do begin
;        zstat=strcompress(string(FORMAT='(F4.1,A1,F4.1)',$
;              min(z2((nr/2):nr-1,*,k)),'-',max(z2((nr/2):nr-1,*,k))))
;        pstat=strcompress(string(FORMAT='(F6.1,A1,F6.1)',$
;              min(p2((nr/2):nr-1,*,k)),'-',max(p2((nr/2):nr-1,*,k))))
;        tstat=strcompress(string(FORMAT='(F5.1,A1,F5.1)',$
;              min(t2((nr/2):nr-1,*,k)),'-',max(t2((nr/2):nr-1,*,k))))
;        ustat=strcompress(string(FORMAT='(F5.1,A1,F5.1)',$
;              min(u2((nr/2):nr-1,*,k)),'-',max(u2((nr/2):nr-1,*,k))))
;        pvstat=strcompress(string(FORMAT='(F7.1,A1,F7.1)',$
;              1.e4*min(pv2((nr/2):nr-1,*,k)),'-',1.e4*max(pv2((nr/2):nr-1,*,k))))
;        print,fix(th(k)),'   ',zstat,'   ',pstat,'   ',tstat,'   ',ustat,'   ',pvstat
;    endfor
stop
jump:
endfor		; loop over days
end
