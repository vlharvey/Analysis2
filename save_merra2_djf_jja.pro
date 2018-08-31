;
; save MERRA2 DJF and JJA for all years
;
loadct,39
mcolor=byte(!p.color)
mcolor=fix(mcolor)
device,decompose=0
if mcolor eq 0 then mcolor=255
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_AVG_2'
;
; file listing of MLS CO marker files
;
spawn,'ls '+dir+'*.sav',ifiles
nfile=n_elements(ifiles)
;
; loop over files
;
icount=0L
for ifile=0L,nfile-1L do begin
    result=strsplit(ifiles(ifile),'_',/extract)
    yyyymm=result(-2)
    smon=strmid(yyyymm,4,2)
    imon=long(smon)

    if long(yyyymm) lt 200408L then goto,skipmonth
    if imon ne 1 and imon ne 2 and imon ne 12 and imon ne 6 and imon ne 7 and imon ne 8 then goto,skipmonth             ; DJF and JJA only

    print,'reading ',ifiles(ifile)
    restore,ifiles(ifile)
;
; declare DJF and JJA arrays
;
        if icount eq 0L then begin
           djf_pv2avg=0.*pv2avg
           djf_p2avg=0.*pv2avg
           djf_u2avg=0.*pv2avg
           djf_v2avg=0.*pv2avg
           djf_qdf2avg=0.*pv2avg
           djf_mark2avg=0.*pv2avg
           djf_qv2avg=0.*pv2avg
           djf_z2avg=0.*pv2avg
           djf_sf2avg=0.*pv2avg
           djf_q2avg=0.*pv2avg
           djf_o32avg=0.*pv2avg
           djf_n2avg=0.*pv2avg

           jja_pv2avg=0.*pv2avg
           jja_p2avg=0.*pv2avg
           jja_u2avg=0.*pv2avg
           jja_v2avg=0.*pv2avg
           jja_qdf2avg=0.*pv2avg
           jja_mark2avg=0.*pv2avg
           jja_qv2avg=0.*pv2avg
           jja_z2avg=0.*pv2avg
           jja_sf2avg=0.*pv2avg
           jja_q2avg=0.*pv2avg
           jja_o32avg=0.*pv2avg
           jja_n2avg=0.*pv2avg
           icount=1L
        endif
;
; DJF
;
    if smon eq '12' or smon eq '01' or smon eq '02' then begin
        good=where(p2avg ne 0.)
        djf_pv2avg(good)=djf_pv2avg(good)+pv2avg(good)
        djf_p2avg(good)=djf_p2avg(good)+p2avg(good)
        djf_u2avg(good)=djf_u2avg(good)+u2avg(good)
        djf_v2avg(good)=djf_v2avg(good)+v2avg(good)
        djf_qdf2avg(good)=djf_qdf2avg(good)+qdf2avg(good)
        djf_mark2avg(good)=djf_mark2avg(good)+mark2avg(good)
        djf_qv2avg(good)=djf_qv2avg(good)+qv2avg(good)
        djf_z2avg(good)=djf_z2avg(good)+z2avg(good)
        djf_sf2avg(good)=djf_sf2avg(good)+sf2avg(good)
        djf_q2avg(good)=djf_q2avg(good)+q2avg(good)
        djf_o32avg(good)=djf_o32avg(good)+o32avg(good)
        djf_n2avg(good)=djf_n2avg(good)+1.
    endif
;
; JJA
;
    if smon eq '06' or smon eq '07' or smon eq '08' then begin
        good=where(p2avg ne 0.)
        jja_pv2avg(good)=jja_pv2avg(good)+pv2avg(good)
        jja_p2avg(good)=jja_p2avg(good)+p2avg(good)
        jja_u2avg(good)=jja_u2avg(good)+u2avg(good)
        jja_v2avg(good)=jja_v2avg(good)+v2avg(good)
        jja_qdf2avg(good)=jja_qdf2avg(good)+qdf2avg(good)
        jja_mark2avg(good)=jja_mark2avg(good)+mark2avg(good)
        jja_qv2avg(good)=jja_qv2avg(good)+qv2avg(good)
        jja_z2avg(good)=jja_z2avg(good)+z2avg(good)
        jja_sf2avg(good)=jja_sf2avg(good)+sf2avg(good)
        jja_q2avg(good)=jja_q2avg(good)+q2avg(good)
        jja_o32avg(good)=jja_o32avg(good)+o32avg(good)
        jja_n2avg(good)=jja_n2avg(good)+1.
    endif

skipmonth:
endfor  ; loop over files
;
; average
;
good=where(djf_n2avg gt 0.)
djf_pv2avg(good)=djf_pv2avg(good)/djf_n2avg(good)
djf_p2avg(good)=djf_p2avg(good)/djf_n2avg(good)
djf_u2avg(good)=djf_u2avg(good)/djf_n2avg(good)
djf_v2avg(good)=djf_v2avg(good)/djf_n2avg(good)
djf_qdf2avg(good)=djf_qdf2avg(good)/djf_n2avg(good)
djf_mark2avg(good)=djf_mark2avg(good)/djf_n2avg(good)
djf_qv2avg(good)=djf_qv2avg(good)/djf_n2avg(good)
djf_z2avg(good)=djf_z2avg(good)/djf_n2avg(good)
djf_sf2avg(good)=djf_sf2avg(good)/djf_n2avg(good)
djf_q2avg(good)=djf_q2avg(good)/djf_n2avg(good)
djf_o32avg(good)=djf_o32avg(good)/djf_n2avg(good)

good=where(jja_n2avg gt 0.)
jja_pv2avg(good)=jja_pv2avg(good)/jja_n2avg(good)
jja_p2avg(good)=jja_p2avg(good)/jja_n2avg(good)
jja_u2avg(good)=jja_u2avg(good)/jja_n2avg(good)
jja_v2avg(good)=jja_v2avg(good)/jja_n2avg(good)
jja_qdf2avg(good)=jja_qdf2avg(good)/jja_n2avg(good)
jja_mark2avg(good)=jja_mark2avg(good)/jja_n2avg(good)
jja_qv2avg(good)=jja_qv2avg(good)/jja_n2avg(good)
jja_z2avg(good)=jja_z2avg(good)/jja_n2avg(good)
jja_sf2avg(good)=jja_sf2avg(good)/jja_n2avg(good)
jja_q2avg(good)=jja_q2avg(good)/jja_n2avg(good)
jja_o32avg(good)=jja_o32avg(good)/jja_n2avg(good)

erase
!type=2^2+2^3
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
omin=-100.
omax=100.
level=omin+((omax-omin)/nlvls)*findgen(nlvls+1)
plotarray=mean(jja_u2avg,dim=2)
contour,plotarray,alat,th,levels=level,c_color=col1,/cell_fill,/noeras,position = [0.1,0.25,0.45,0.75],$	;x0,y0,x1,y1],$
        xrange=[-90,90],xticks=6,ytitle='Theta (K)',yrange=[min(th),max(th)],charsize=1.5,min_value=-99.,title='JJA'
index=where(level gt 0)
contour,plotarray,alat,th,levels=level(index),c_color=0,/follow,/noeras,/overplot
index=where(level lt 0)
contour,plotarray,alat,th,levels=level(index),c_color=mcolor,c_linestyle=5,/follow,/noeras,/overplot
dum=mean(jja_mark2avg,dim=2)
contour,dum,alat,th,levels=0.1+0.1*findgen(10),color=0,thick=8,/follow,/noeras,/overplot

plotarray=mean(djf_u2avg,dim=2)
contour,plotarray,alat,th,levels=level,c_color=col1,/cell_fill,/noeras,position = [0.55,0.25,0.9,0.75],$ 
        xrange=[-90,90],xticks=6,ytitle='Theta (K)',yrange=[min(th),max(th)],charsize=1.5,min_value=-99.,title='DJF'
index=where(level gt 0)
contour,plotarray,alat,th,levels=level(index),c_color=0,/follow,/noeras,/overplot
index=where(level lt 0)
contour,plotarray,alat,th,levels=level(index),c_color=mcolor,c_linestyle=5,/follow,/noeras,/overplot
dum=mean(djf_mark2avg,dim=2)
contour,dum,alat,th,levels=0.1+0.1*findgen(10),color=0,thick=8,/follow,/noeras,/overplot
;
; save djf/jja file
;
ofile='MERRA2_djf_jja.sav'
save,file=ofile,alat,alon,th,djf_pv2avg,djf_p2avg,djf_u2avg,djf_v2avg,djf_qdf2avg,djf_mark2avg,djf_qv2avg,$
     djf_z2avg,djf_sf2avg,djf_q2avg,djf_o32avg,jja_pv2avg,jja_p2avg,jja_u2avg,jja_v2avg,jja_qdf2avg,$
     jja_mark2avg,jja_qv2avg,jja_z2avg,jja_sf2avg,jja_q2avg,jja_o32avg
end
