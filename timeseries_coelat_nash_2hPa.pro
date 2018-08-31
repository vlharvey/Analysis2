;
; timeseries of Elat edge values based on CO and PV
;
loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,0.15,0.15,0.45,0.75]
yorig=[0.65,0.4,0.1,0.1,0.1]
cbaryoff=0.02
cbarydel=0.01
xlen=0.8
ylen=0.3
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
nyear=n_elements(syear)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; get file listing
;
spawn,'ls elat_time_co_edges_????-????.sav',ifiles
restore,ifiles(0)
lowlat_elatedge_time_all=LOWLAT_ELATEDGE_TIME
nashedge_time_all=nashedge_time
sdate_time_all=sdate_time

for ifile=1L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    LOWLAT_ELATEDGE_TIME_all=[LOWLAT_ELATEDGE_TIME_all,LOWLAT_ELATEDGE_TIME]
    nashedge_time_all=[nashedge_time_all,nashedge_time]
    sdate_time_all=[sdate_time_all,sdate_time]
endfor

; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_coelat_nash_'+spress+'.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; DJF
;
syear=strmid(sdate_time_all,0,4)
smon=strmid(sdate_time_all,4,2)
sday=strmid(sdate_time_all,6,2)
index=where(smon ne '11' and smon ne '03')
LOWLAT_ELATEDGE_TIME_all=LOWLAT_ELATEDGE_TIME_all(index)
NASHEDGE_TIME_all=NASHEDGE_TIME_all(index)
sdate_time_all=sdate_time_all(index)
syear=strmid(sdate_time_all,0,4)
smon=strmid(sdate_time_all,4,2)
sday=strmid(sdate_time_all,6,2)

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
xlabs0=strmid(sdate_time_all(xindex),2,2)
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
xlabs1=strmid(sdate_time_all(xindex),2,2)
xlabs=xlabs0+xlabs1
plot,findgen(n_elements(sdate_time_all)),LOWLAT_ELATEDGE_TIME_ALL,color=0,/noeras,charsize=1.25,ytitle='Equivalent Latitude',yrange=[20,90],$
     title='DJF '+spress,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,psym=8,/nodata
index=where(LOWLAT_ELATEDGE_TIME_ALL gt 10.)
oplot,index,LOWLAT_ELATEDGE_TIME_ALL(index),psym=8,color=0,symsize=0.75
oplot,findgen(n_elements(sdate_time_all)),nashedge_time_all,color=mcolor*.9,psym=8,symsize=0.5
xyouts,xmn+0.03,ymn+0.05,'CO',/normal,color=0,charsize=2,charthick=4
xyouts,xmn+0.03,ymn+0.01,'PV',/normal,color=mcolor*.9,charsize=2,charthick=4
;
; indicate new December
;
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
for i=0L,nxticks-1L do begin
    plots,xindex(i),20
    plots,xindex(i),90,/continue,color=0
endfor
;
; Difference between CO Elat and PV Elat
;
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+0.2
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
xlabs0=strmid(sdate_time_all(xindex),2,2)
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
xlabs1=strmid(sdate_time_all(xindex),2,2)
xlabs=xlabs0+xlabs1
diff=LOWLAT_ELATEDGE_TIME_ALL-nashedge_time_all
ymax=60
ymin=-20
plot,findgen(n_elements(sdate_time_all)),diff,color=0,/noeras,charsize=1.25,ytitle='!7D!3Elat',yrange=[ymin,ymax],psym=8,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,symsize=0.8,/nodata
index=where(LOWLAT_ELATEDGE_TIME_ALL gt 10.)
oplot,index,diff,symsize=0.5,color=0,psym=8
oplot,findgen(n_elements(sdate_time_all)),0.*diff,color=0
;
; indicate new December
;
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
for i=0L,nxticks-1L do begin
    plots,xindex(i),ymin
    plots,xindex(i),ymax,/continue,color=0,/data
endfor
;
; monthly PDFs
;
index=where(smon eq '12')
imin=-20.
imax=60.
iinc=4.
y2=histogram(diff(index),min=imin,max=imax,binsize=iinc)/float(n_elements(index))
xmn=xorig(2)
xmx=xorig(2)+0.2
ymn=yorig(2)
ymx=yorig(2)+0.2
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=long((imax-imin)/iinc + 1)
level=imin+iinc*findgen(nlvls)
plot,level,smooth(y2,3),color=0,title='DEC',xtitle='!7D!3Elat',ytitle='Frequency',thick=5
xyouts,xmx-0.06,ymx-0.03,string(format='(f4.1)',median(diff(index))),color=0,/normal

index=where(smon eq '01')
y2=histogram(diff(index),min=imin,max=imax,binsize=iinc)/float(n_elements(index))
xmn=xorig(3)
xmx=xorig(3)+0.2
ymn=yorig(3)
ymx=yorig(3)+0.2
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=long((imax-imin)/iinc + 1)
level=imin+iinc*findgen(nlvls)
plot,level,smooth(y2,3),color=0,title='JAN',xtitle='!7D!3Elat',thick=5
xyouts,xmx-0.06,ymx-0.03,string(format='(f4.1)',median(diff(index))),color=0,/normal

index=where(smon eq '02')
y2=histogram(diff(index),min=imin,max=imax,binsize=iinc)/float(n_elements(index))
xmn=xorig(4)
xmx=xorig(4)+0.2
ymn=yorig(4) 
ymx=yorig(4)+0.2
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=long((imax-imin)/iinc + 1)
level=imin+iinc*findgen(nlvls)
plot,level,smooth(y2,3),color=0,title='FEB',xtitle='!7D!3Elat',thick=5
xyouts,xmx-0.06,ymx-0.03,string(format='(f4.1)',median(diff(index))),color=0,/normal

; Close PostScript file and return control to X-windows
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_coelat_nash_'+spress+'.ps -rotate -90 '+$
                       'timeseries_coelat_nash_'+spress+'.jpg'
;  spawn,'rm -f timeseries_coelat_nash_'+spress+'.ps'
endif
end
