;
; plot daily median species in each hem to make zt 
; i.e. daily hemispheric averages
;
loadct,38
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=700
nydim=700
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   lc=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
ifile='               '
close,1
openr,1,'sage+haloe+poam+ilas.fil'
nmonth=0L
readf,1,nmonth
for imonth=0,nmonth-1 do begin
    readf,1,ifile
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='daily_o3_medians_'+ifile+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
    close,2
    openr,2,'../Datfiles/'+ifile
    nday=0L & nth=0L
    readu,2,nday,nth
    days=fltarr(nday)
    th=fltarr(nth)
    readu,2,days,th
    nh_all_o3_median=9999+fltarr(nday,nth)
    nh_all_h2o_median=9999+fltarr(nday,nth)
    nh_all_ch4_median=9999+fltarr(nday,nth)
    nh_all_no2_median=9999+fltarr(nday,nth)
    nh_all_sad_median=9999+fltarr(nday,nth)
    sh_all_o3_median=9999+fltarr(nday,nth)
    sh_all_h2o_median=9999+fltarr(nday,nth)
    sh_all_ch4_median=9999+fltarr(nday,nth)
    sh_all_no2_median=9999+fltarr(nday,nth)
    sh_all_sad_median=9999+fltarr(nday,nth)
    nh_vortex_o3_median=9999+fltarr(nday,nth)
    nh_vortex_h2o_median=9999+fltarr(nday,nth)
    nh_vortex_ch4_median=9999+fltarr(nday,nth)
    nh_vortex_no2_median=9999+fltarr(nday,nth)
    nh_vortex_sad_median=9999+fltarr(nday,nth)
    sh_vortex_o3_median=9999+fltarr(nday,nth)
    sh_vortex_h2o_median=9999+fltarr(nday,nth)
    sh_vortex_ch4_median=9999+fltarr(nday,nth)
    sh_vortex_no2_median=9999+fltarr(nday,nth)
    sh_vortex_sad_median=9999+fltarr(nday,nth)
    nh_high_o3_median=9999+fltarr(nday,nth)
    nh_high_h2o_median=9999+fltarr(nday,nth)
    nh_high_ch4_median=9999+fltarr(nday,nth)
    nh_high_no2_median=9999+fltarr(nday,nth)
    nh_high_sad_median=9999+fltarr(nday,nth)
    sh_high_o3_median=9999+fltarr(nday,nth)
    sh_high_h2o_median=9999+fltarr(nday,nth)
    sh_high_ch4_median=9999+fltarr(nday,nth)
    sh_high_no2_median=9999+fltarr(nday,nth)
    sh_high_sad_median=9999+fltarr(nday,nth)
    nh_out_o3_median=9999+fltarr(nday,nth)
    nh_out_h2o_median=9999+fltarr(nday,nth)
    nh_out_ch4_median=9999+fltarr(nday,nth)
    nh_out_no2_median=9999+fltarr(nday,nth)
    nh_out_sad_median=9999+fltarr(nday,nth)
    sh_out_o3_median=9999+fltarr(nday,nth)
    sh_out_h2o_median=9999+fltarr(nday,nth)
    sh_out_ch4_median=9999+fltarr(nday,nth)
    sh_out_no2_median=9999+fltarr(nday,nth)
    sh_out_sad_median=9999+fltarr(nday,nth)
    readu,2,nh_all_o3_median, nh_all_h2o_median, nh_all_ch4_median,$
       nh_all_no2_median, nh_all_sad_median, sh_all_o3_median,$
       sh_all_h2o_median, sh_all_ch4_median, sh_all_no2_median,$
       sh_all_sad_median
    readu,2,nh_vortex_o3_median, nh_vortex_h2o_median, nh_vortex_ch4_median,$
       nh_vortex_no2_median, nh_vortex_sad_median, sh_vortex_o3_median,$
       sh_vortex_h2o_median, sh_vortex_ch4_median, sh_vortex_no2_median,$
       sh_vortex_sad_median
    readu,2,nh_high_o3_median, nh_high_h2o_median, nh_high_ch4_median,$
       nh_high_no2_median, nh_high_sad_median, sh_high_o3_median,$
       sh_high_h2o_median, sh_high_ch4_median, sh_high_no2_median,$
       sh_high_sad_median
    readu,2,nh_out_o3_median, nh_out_h2o_median, nh_out_ch4_median,$
       nh_out_no2_median, nh_out_sad_median, sh_out_o3_median,$
       sh_out_h2o_median, sh_out_ch4_median, sh_out_no2_median,$
       sh_out_sad_median
    nh_all_o3_sigma=9999+fltarr(nday,nth)
    nh_all_h2o_sigma=9999+fltarr(nday,nth)
    nh_all_ch4_sigma=9999+fltarr(nday,nth)
    nh_all_no2_sigma=9999+fltarr(nday,nth)
    nh_all_sad_sigma=9999+fltarr(nday,nth)
    sh_all_o3_sigma=9999+fltarr(nday,nth)
    sh_all_h2o_sigma=9999+fltarr(nday,nth)
    sh_all_ch4_sigma=9999+fltarr(nday,nth)
    sh_all_no2_sigma=9999+fltarr(nday,nth)
    sh_all_sad_sigma=9999+fltarr(nday,nth)
    nh_vortex_o3_sigma=9999+fltarr(nday,nth)
    nh_vortex_h2o_sigma=9999+fltarr(nday,nth)
    nh_vortex_ch4_sigma=9999+fltarr(nday,nth)
    nh_vortex_no2_sigma=9999+fltarr(nday,nth)
    nh_vortex_sad_sigma=9999+fltarr(nday,nth)
    sh_vortex_o3_sigma=9999+fltarr(nday,nth)
    sh_vortex_h2o_sigma=9999+fltarr(nday,nth)
    sh_vortex_ch4_sigma=9999+fltarr(nday,nth)
    sh_vortex_no2_sigma=9999+fltarr(nday,nth)
    sh_vortex_sad_sigma=9999+fltarr(nday,nth)
    nh_high_o3_sigma=9999+fltarr(nday,nth)
    nh_high_h2o_sigma=9999+fltarr(nday,nth)
    nh_high_ch4_sigma=9999+fltarr(nday,nth)
    nh_high_no2_sigma=9999+fltarr(nday,nth)
    nh_high_sad_sigma=9999+fltarr(nday,nth)
    sh_high_o3_sigma=9999+fltarr(nday,nth)
    sh_high_h2o_sigma=9999+fltarr(nday,nth)
    sh_high_ch4_sigma=9999+fltarr(nday,nth)
    sh_high_no2_sigma=9999+fltarr(nday,nth)
    sh_high_sad_sigma=9999+fltarr(nday,nth)
    nh_out_o3_sigma=9999+fltarr(nday,nth)
    nh_out_h2o_sigma=9999+fltarr(nday,nth)
    nh_out_ch4_sigma=9999+fltarr(nday,nth)
    nh_out_no2_sigma=9999+fltarr(nday,nth)
    nh_out_sad_sigma=9999+fltarr(nday,nth)
    sh_out_o3_sigma=9999+fltarr(nday,nth)
    sh_out_h2o_sigma=9999+fltarr(nday,nth)
    sh_out_ch4_sigma=9999+fltarr(nday,nth)
    sh_out_no2_sigma=9999+fltarr(nday,nth)
    sh_out_sad_sigma=9999+fltarr(nday,nth)
    readu,2,nh_all_o3_sigma, nh_all_h2o_sigma, nh_all_ch4_sigma,$
       nh_all_no2_sigma, nh_all_sad_sigma, sh_all_o3_sigma,$
       sh_all_h2o_sigma, sh_all_ch4_sigma, sh_all_no2_sigma,$
       sh_all_sad_sigma
    readu,2,nh_vortex_o3_sigma, nh_vortex_h2o_sigma, nh_vortex_ch4_sigma,$
       nh_vortex_no2_sigma, nh_vortex_sad_sigma, sh_vortex_o3_sigma,$
       sh_vortex_h2o_sigma, sh_vortex_ch4_sigma, sh_vortex_no2_sigma,$
       sh_vortex_sad_sigma
    readu,2,nh_high_o3_sigma, nh_high_h2o_sigma, nh_high_ch4_sigma,$
       nh_high_no2_sigma, nh_high_sad_sigma, sh_high_o3_sigma,$
       sh_high_h2o_sigma, sh_high_ch4_sigma, sh_high_no2_sigma,$
       sh_high_sad_sigma
    readu,2,nh_out_o3_sigma, nh_out_h2o_sigma, nh_out_ch4_sigma,$
       nh_out_no2_sigma, nh_out_sad_sigma, sh_out_o3_sigma,$
       sh_out_h2o_sigma, sh_out_ch4_sigma, sh_out_no2_sigma,$
       sh_out_sad_sigma
    close,2
;
; plot
;
    erase
    !type=2^2+2^3
    level=-2.+.2*findgen(20)
    xyouts,.45,.9,strmid(ifile,24,6),/normal,charsize=2
    set_viewport,.15,.85,.55,.85
    tmp=0.*nh_high_o3_median
    index=where(nh_high_o3_median ne 0.)
    tmp(index)=(nh_high_o3_median(index)-nh_out_o3_median(index))*1.e6
    contour,tmp,days,th,/noeras,title='NH High-Ambient Ozone',$
         ytitle='Theta',xrange=[0.,nday-1],yrange=[th(nth-1),th(0)],$
         levels=level,max_value=9999.,c_labels=findgen(20) ge 0,$
         c_linestyle=level lt 0
    set_viewport,.15,.85,.15,.45
    tmp=0.*sh_high_o3_median
    index=where(sh_high_o3_median ne 0.)
    tmp(index)=(sh_high_o3_median(index)-sh_out_o3_median(index))*1.e6
    contour,tmp,days,th,/noeras,title='SH High-Ambient Ozone',$
         ytitle='Theta',xrange=[0.,nday-1],yrange=[th(nth-1),th(0)],$
         levels=level,max_value=9999.,c_labels=findgen(20) ge 0,$
         c_linestyle=level lt 0

    if setplot eq 'ps' then device,/close
stop
endfor		; loop over months
end
