;
; check smoothed edge
;
@rd_ukmo_nc3
vortexmin=1.
highmax=-1.
outmin=0.
outmax=0.
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
close,1
openr,1,'ukmo_sh_files_00.fil'
nfile=0L
readf,1,nfile
for n=0,nfile-1 do begin
    readf,1,ifile
    iflag=0
    rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                pv2,p2,msf2,u2,v2,q2,qdf2,marksf2,vp2,sf2,iflag
    index=where(marksf2 ne 0.)
    if index(0) ne -1 then $
       marksf2(index)=marksf2(index)/abs(marksf2(index))
    savemark=marksf2
    marksmooth2=0.*marksf2
;
; smooth marker field
;
    for i=0,nc-1 do begin
        ip1=i+1
        im1=i-1
        if i eq 0 then im1=nc-1
        if i eq nc-1 then ip1=0
        for j=0,nr-1 do begin
            jp1=j+1
            jm1=j-1
            if j eq 0 then jm1=j
            if j eq nr-1 then jp1=j
            for k=0,nth-1 do begin
                kp1=k+1
                km1=k-1
                if k eq 0 then km1=k
                if k eq nth-1 then kp1=k
                marksmooth2(j,i,k)=( (1./4.)*(savemark(j,im1,k) + $
                               2.0*savemark(j,i,k)+savemark(j,ip1,k)) + $
                                  (1./4.)*(savemark(jm1,i,k) + $
                               2.0*savemark(j,i,k)+savemark(jp1,i,k)) + $
                                  (1./4.)*(savemark(j,i,km1) + $
                               2.0*savemark(j,i,k)+savemark(j,i,kp1)))/3.0
             endfor
        endfor
    endfor

index=where(marksf2 eq 1.)
if index(0) ne -1 then begin
   if min(marksmooth2(index)) lt vortexmin then begin
      vortexmin=min(marksmooth2(index))
      print,'Min in Vortex ',min(marksmooth2(index))
   endif
endif
index=where(marksf2 eq -1.) 
if index(0) ne -1 then begin
   if max(marksmooth2(index)) gt highmax then begin
      highmax=max(marksmooth2(index))
      print,'Max in Highs  ',max(marksmooth2(index))
   endif
endif
index=where(marksf2 eq 0.)
if index(0) ne -1 then begin
   if min(marksmooth2(index)) lt outmin then begin
      outmin=min(marksmooth2(index))
      print,'Outside min ',min(marksmooth2(index))
   endif
   if max(marksmooth2(index)) gt outmax then begin
      outmax=max(marksmooth2(index))
      print,'Outside max ',max(marksmooth2(index))
   endif
endif
endfor		; loop over files
end
