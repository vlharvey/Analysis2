;
; create MERRA2 3D monthly means 4 UT
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra2_nc3

loadct,39
mcolor=byte(!p.color)
mcolor=fix(mcolor)
device,decompose=0
if mcolor eq 0 then mcolor=255
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
;model_years=1990+indgen(26)	; 1990-2015
model_years=[2017]
nyear=n_elements(model_years)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmonth=n_elements(smon)
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
stime=['00','06','12','18']
ntime=n_elements(stime)
;
; loop over model years
;
for iyear=model_years(0),model_years(nyear-1L) do begin
for imon=0L,nmonth-1L do begin
    if iyear mod 4 eq 0L then mday(1)=29
    if iyear mod 4 ne 0L then mday(1)=28
for itime=0L,ntime-1L do begin
    for iday=0L,mday(imon)-1L do begin

        sdate=strcompress(string(FORMAT='(I4,I2.2,I2.2,A2)',iyear,imon+1,iday+1,stime(itime)))
        rd_merra2_nc3,dir+sdate+'.nc3',nc,nr,nth,alon,alat,th,pv2,p2,u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag
        if iflag ne 0L then stop,'missing file on '+sdate
; 
; re-declare arrays at the beginning of each month
;
        if iday eq 0L then begin
           pv2avg=0.*pv2
           p2avg=0.*pv2
           u2avg=0.*pv2
           v2avg=0.*pv2
           qdf2avg=0.*pv2
           mark2avg=0.*pv2
           qv2avg=0.*pv2
           z2avg=0.*pv2
           sf2avg=0.*pv2
           q2avg=0.*pv2
           o32avg=0.*pv2
           n2avg=0.*pv2
        endif
        good=where(p2 ne 0.)
        pv2avg(good)=pv2avg(good)+pv2(good)
        p2avg(good)=p2avg(good)+p2(good)
        u2avg(good)=u2avg(good)+u2(good)
        v2avg(good)=v2avg(good)+v2(good)
        qdf2avg(good)=qdf2avg(good)+qdf2(good)
        mark2avg(good)=mark2avg(good)+mark2(good)
        qv2avg(good)=qv2avg(good)+qv2(good)
        z2avg(good)=z2avg(good)+z2(good)
        sf2avg(good)=sf2avg(good)+sf2(good)
        q2avg(good)=q2avg(good)+q2(good)
        o32avg(good)=o32avg(good)+o32(good)
        n2avg(good)=n2avg(good)+1.
;
; check
;
;erase
;!type=2^2+2^3
;nlvls=20
;col1=1+mcolor*findgen(20)/nlvls
;plotarray=mean(u2,dim=2)
;omin=-100.
;omax=100.
;level=omin+((omax-omin)/nlvls)*findgen(nlvls+1)
;contour,plotarray,alat,th,levels=level,c_color=col1,/cell_fill,/noeras,$
;        xrange=[-90,90],xticks=6,ytitle='Theta (K)',yrange=[min(th),max(th)],charsize=1.5,min_value=-99.,title='MERRA Ubar '+sdate
;index=where(level gt 0)
;contour,plotarray,alat,th,levels=level(index),c_color=0,/follow,/noeras,/overplot
;index=where(level lt 0)
;contour,plotarray,alat,th,levels=level(index),c_color=mcolor,c_linestyle=5,/follow,/noeras,/overplot

    endfor	; loop over days
;
; average
;
    good=where(n2avg gt 0.)
    pv2avg(good)=pv2avg(good)/n2avg(good)
    p2avg(good)=p2avg(good)/n2avg(good)
    u2avg(good)=u2avg(good)/n2avg(good)
    v2avg(good)=v2avg(good)/n2avg(good)
    qdf2avg(good)=qdf2avg(good)/n2avg(good)
    mark2avg(good)=mark2avg(good)/n2avg(good)
    qv2avg(good)=qv2avg(good)/n2avg(good)
    z2avg(good)=z2avg(good)/n2avg(good)
    sf2avg(good)=sf2avg(good)/n2avg(good)
    q2avg(good)=q2avg(good)/n2avg(good)
    o32avg(good)=o32avg(good)/n2avg(good)
erase
!type=2^2+2^3
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
plotarray=mean(u2,dim=2)
omin=-100.
omax=100.
level=omin+((omax-omin)/nlvls)*findgen(nlvls+1)
plotarray=mean(u2avg,dim=2)
contour,plotarray,alat,th,levels=level,c_color=col1,/cell_fill,/noeras,$
        xrange=[-90,90],xticks=6,ytitle='Theta (K)',yrange=[min(th),max(th)],charsize=1.5,min_value=-99.,title='MERRA Ubar '+string(imon+1)+'/'+string(iyear)+' '+stime(itime)+'Z'
index=where(level gt 0)
contour,plotarray,alat,th,levels=level(index),c_color=0,/follow,/noeras,/overplot
index=where(level lt 0)
contour,plotarray,alat,th,levels=level(index),c_color=mcolor,c_linestyle=5,/follow,/noeras,/overplot
dum=mean(mark2avg,dim=2)
contour,dum,alat,th,levels=0.1+0.1*findgen(10),color=0,thick=8,/follow,/noeras,/overplot
;
; save monthly file
;
    ofile=dir+'AVG_'+strcompress(string(FORMAT='(I4,I2.2,A1,A2)',iyear,imon+1,'_',stime(itime)))+'.sav'
    save,file=ofile,alat,alon,th,pv2avg,p2avg,u2avg,v2avg,qdf2avg,mark2avg,qv2avg,z2avg,sf2avg,q2avg,o32avg
endfor	; loop over months
endfor	; loop over years
endfor	; loop over UT times
end
