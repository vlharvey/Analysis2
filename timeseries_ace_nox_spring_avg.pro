;
; spring time period only
; plot timeseries of vortex NOx at different altitudes for each year 2004, 2005, 2006
;
@fillit
@smoothit

loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='x'
read,'setplot=',setplot
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
nxdim=600 & nydim=600
xorig=[0.20]
yorig=[0.30]
xlen=0.7
ylen=0.4
cbaryoff=0.08
cbarydel=0.02
!NOERAS=-1
!p.font=1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
syear=['2004','2005','2006']
nyear=n_elements(syear)
nlvls=21L
col1=1L+indgen(nlvls)*mcolor/float(nlvls)
;
; restore NOx files
;
restore,'vortex_nox_2004.sav
nday4=nday
FDOY4=FDOY
LATITUDE4=LATITUDE
ONOX4=ONOX
onox4(0:35,*)=0.
restore,'vortex_nox_2005.sav
nday5=nday
FDOY5=FDOY
LATITUDE5=LATITUDE
ONOX5=ONOX
restore,'vortex_nox_2006.sav
nday6=nday
FDOY6=FDOY
LATITUDE6=LATITUDE
ONOX6=ONOX
nz=n_elements(ALTITUDE)
;
; loop over altitudes
;
for kk=0L,nz-1L do begin
    salt=strcompress(long(altitude(kk)),/remove_all)
;
; postscript file
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       device,font_size=9
       device,/landscape,bits=8,filename='timeseries_ace_nox_spring_'+salt+'km_avg.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif
;
; extract vortex NOx at this altitude
;
    nox2004=reform(ONOX4(*,kk))
    nox2005=reform(ONOX5(*,kk))
    nox2006=reform(ONOX6(*,kk))
;
; plot timeseries
;
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
;   kday=365.
;kday=long(max(fdoy))	; uncomment to "zoom" in on partial year
kday=91
    nmin=0.
    nmax=max([nox2004,nox2005,nox2006])
    if nmax le 0. then goto,jumplev
;
; make same scale as plot with all data
;
      if altitude(kk) ge 71.0000 then nmax=     5000.00
      if altitude(kk) ge 61.0000 and altitude(kk) lt 71. then nmax=     1000.00
      if altitude(kk) eq 60.0000 then nmax=     980.000
      if altitude(kk) eq 59.0000 then nmax=     920.000
      if altitude(kk) eq 58.0000 then nmax=     912.000
      if altitude(kk) eq 57.0000 then nmax=     1000.00
      if altitude(kk) eq 56.0000 then nmax=     1000.00
      if altitude(kk) eq 55.0000 then nmax=     1000.00
      if altitude(kk) eq 54.0000 then nmax=     1000.00
      if altitude(kk) eq 53.0000 then nmax=     1000.00
      if altitude(kk) eq 52.0000 then nmax=     1000.00
      if altitude(kk) eq 51.0000 then nmax=     1000.00
      if altitude(kk) eq 50.0000 then nmax=     874.285
      if altitude(kk) eq 49.0000 then nmax=     785.160
      if altitude(kk) eq 48.0000 then nmax=     671.100
      if altitude(kk) eq 47.0000 then nmax=     554.125
      if altitude(kk) eq 46.0000 then nmax=     463.795
      if altitude(kk) eq 45.0000 then nmax=     388.305
      if altitude(kk) eq 44.0000 then nmax=     294.350
      if altitude(kk) eq 43.0000 then nmax=     236.185
      if altitude(kk) eq 42.0000 then nmax=     186.270
      if altitude(kk) eq 41.0000 then nmax=     137.070
      if altitude(kk) eq 40.0000 then nmax=     96.1750
      if altitude(kk) eq 39.0000 then nmax=     60.8450
      if altitude(kk) eq 38.0000 then nmax=     31.3850
      if altitude(kk) eq 37.0000 then nmax=     39.5350
      if altitude(kk) eq 36.0000 then nmax=     21.9850
      if altitude(kk) eq 35.0000 then nmax=     19.3000
      if altitude(kk) eq 34.0000 then nmax=     17.4650
      if altitude(kk) eq 33.0000 then nmax=     15.2550
      if altitude(kk) eq 32.0000 then nmax=     14.8250
      if altitude(kk) eq 31.0000 then nmax=     15.6150
      if altitude(kk) eq 30.0000 then nmax=     13.3300

    if altitude(kk) gt 80. or altitude(kk) lt 30. then goto,jumplev
    plot,1.+findgen(nday4),nox2004,psym=8,/noeras,yrange=[nmin,nmax],xticks=3,$
         xtickname=['                            Jan',$
                    '                            Feb',$
                    '                            Mar','  '],$
         charsize=1.75,ytitle='NOx (ppbv)',color=0,xrange=[1.,kday],$
         title='ACE Daily Avg NOx < 50 N '+salt+' km'
    oplot,1.+findgen(nday4),nox2004,psym=8,color=mcolor*.9
    oplot,1.+findgen(nday4),nox2004,psym=0,color=mcolor*.9,thick=2
    oplot,1.+findgen(nday5),nox2005,psym=8,color=mcolor*.45
    oplot,1.+findgen(nday5),nox2005,psym=0,color=mcolor*.45,thick=2
    oplot,1.+findgen(nday6),nox2006,psym=8,color=0
    oplot,1.+findgen(nday6),nox2006,psym=0,color=0,thick=2
    xyouts,xmn+0.02,ymx-0.03,'2004',/normal,charsize=2,color=mcolor*.9
    xyouts,xmn+0.02,ymx-0.06,'2005',/normal,charsize=2,color=mcolor*.45
    xyouts,xmn+0.02,ymx-0.09,'2006',/normal,charsize=2,color=0
    set_viewport,xorig(0),xorig(0)+xlen,0.15,0.25
    plot,fdoy,latitude,psym=1,color=0,yrange=[30.,90.],xrange=[0.,kday],xticks=3,$
         xtickname=' '+strarr(12),charsize=1.75,ytitle='Latitude',yticks=2,$
         ytickv=[30.,60.,90.],symsize=0.5
    oplot,fdoy4,latitude4,psym=1,color=mcolor*.9
    oplot,fdoy5,latitude5,psym=1,color=mcolor*.45
    oplot,fdoy6,latitude6,psym=1,color=0

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_ace_nox_spring_'+salt+'km_avg.ps -rotate -90 '+$
             'timeseries_ace_nox_spring_'+salt+'km_avg.jpg'
       spawn,'/usr/bin/rm timeseries_ace_nox_spring_'+salt+'km_avg.ps'
    endif
    jumplev:
endfor	; loop over years
end
