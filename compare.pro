;compare.pro -- For the same rev, compare the # clouds detected in v3.22 vs. v3.20.

restore,'/Volumes/earth/harvey/CIPS_data/Summary_files/cips_l4_summary_v3.22_north_07_01G_all.sav'
rev22=rev
ncld22=num_cld
nobs22=num_obs
doy22=doy
dfs22=dfs
alb22=alb
sza22=sza

restore,'/Volumes/earth/harvey/CIPS_data/Summary_files/cips_l4_summary_v03.20_north_07_01G_all.sav'
rev20=rev
ncld20=num_cld
nobs20=num_obs
doy20=doy
dfs20=dfs
alb20=alb
sza20=sza

n=n_elements(rev22)
index20=lonarr(n)-99
index22=index20
for i=0,n-1 do begin
   x=where(rev20 eq rev22(i),nx)
   if nx eq 1 then begin
      index20(i)=x(0)
      index22(i)=i
   endif
   if nx gt 1 then stop
endfor

good=where(index20 ne -99,norbit)
index20=index20(good)
r20=rev20(index20)
ncld20=ncld20(index20,*)
nobs20=nobs20(index20,*)
doy20=doy20(index20)
dfs20=dfs20(index20)
alb20=alb20(index20,*)
sza20=sza20(index20,*)
;
; make 2d rev number array to identify rev number where v3.20 has more clouds than v3.22
;
rev20_2d=0*nobs20
for i=0,nbin-1L do rev20_2d(*,i)=r20
;
; 2d latitude array
;
lat_2d=0*nobs20
for i=0,norbit-1L do lat_2d(i,*)=latlo
;
; 2d DOY array
;
doy_2d=0*nobs20
for i=0,nbin-1L do doy_2d(*,i)=doy20

good=where(index22 ne -99)
index22=index22(good)
r22=rev22(index22)
ncld22=ncld22(index22,*)
nobs22=nobs22(index22,*)
doy22=doy22(index22)
dfs22=dfs22(index22)
alb22=alb22(index22,*)
sza22=sza22(index22,*)

ma=max([ncld20,ncld22])

plot,r20,ncld20(*,0),psym=1,yrange=[0,ma]
for i=0,13 do oplot,r20,ncld20(*,i),psym=1
for i=0,13 do oplot,r22,ncld22(*,i),psym=1,color=240
print,'Number of Cloud Observations'
wait,1.

erase
npts=n_elements(good)
nfrac20=0.*ncld20
nfrac22=0.*ncld22
index=where(nobs20 ne 0.)
nfrac20(index)=float(ncld20(index))/float(nobs20(index))
index=where(nobs22 ne 0.)
nfrac22(index)=float(ncld22(index))/float(nobs22(index))
plot,r20,nfrac20(*,0),psym=1,yrange=[0.,1.]
for i=0,13 do oplot,r20,nfrac20(*,i),psym=1
for i=0,13 do oplot,r22,nfrac22(*,i),psym=1,color=240
print,'Fraction of Clouds'
wait,1.

erase
!type=2^2+2^3
set_viewport,0.1,0.9,0.1,0.9
plot,nfrac20,nfrac22,psym=1,xtitle='Cloud Fraction v3.20',ytitle='Cloud Fraction v3.22',title='north_07_01G_all'
oplot,findgen(10),color=200,thick=5
print,'Fraction of Clouds'

index=where(nfrac20 gt nfrac22)
latbad=reform(lat_2d(index))
doybad=reform(doy_2d(index))
revbad=reform(rev20_2d(index))
nfracbad=reform(nfrac22(index)-nfrac20(index))
ncld22bad=ncld22(index)
ncld20bad=ncld20(index)
nobs22bad=nobs22(index)
nobs20bad=nobs20(index)

x=histogram(doy_2d(index))
dmax=max(doy_2d(index))
dmin=min(doy_2d(index))
nday=n_elements(x)
set_viewport,.6,.85,.175,.45
plot,dmin+1.+findgen(nday),x,psym=1,/noeras,xtitle='DOY',ytitle='# Bad'
oplot,dmin+1.+findgen(nday),x,psym=1,color=240

index=where(doybad eq 148 and revbad eq 487)
print,'DOY = ',doybad(index) 
print,'ORBIT = ',revbad(index) 
print,'LAT = ',latbad(index) 
print,'NCLD v3.20 = ',ncld20bad(index) 
print,'NCLD v3.22 = ',ncld22bad(index) 
print,'NOBS v3.20 = ',nobs20bad(index)
print,'NOBS v3.22 = ',nobs22bad(index)

stop

d=ncld22-ncld20

plot,dfs22,ncld22(*,0),/nodata,yrange=[-50,10],xrange=[-40,80]
index=ncld22*0.0-99
for i=0,13 do begin
   x=where((ncld22(*,i) lt ncld20(*,i)) and (nobs22(*,i) eq nobs20(*,i)),nx)
   if nx gt 0 then begin
      index(0:nx-1,i)=x
      oplot,dfs22(x),d(x,i),psym=1
   endif
endfor

stop

plot,sza22,ncld22(*,0),/nodata,yrange=[-50,0],xrange=[-10,30]
index=ncld22*0.0-99
for i=0,13 do begin
   x=where((ncld22(*,i) lt ncld20(*,i)) and (nobs22(*,i) eq nobs20(*,i)),nx)
   if nx gt 0 then begin
      index(0:nx-1,i)=x
      oplot,alb22(x,i),d(x,i),psym=1
      print,min(alb22(x,i))
   endif
endfor

end
