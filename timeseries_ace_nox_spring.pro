;
; not daily averages, all data
; spring time period only
; plot timeseries of vortex NOx at different altitudes for each year 2004, 2005, 2006
;
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
; restore ACE SOSST data
;
restore,dira+'dmps_ace_v2.2.meto.2004
restore,dira+'cat_ace_v2.2.2004
restore,dira+'no_ace_v2.2.2004
nomix=mix & nomask=mask
restore,dira+'no2_ace_v2.2.2004
no2mix=mix & no2mask=mask
noxmix=0.*nomix
index=where(nomask ne -99. and no2mask ne -99.)     ; both good
if index(0) ne -1L then noxmix(index)=nomix(index)+no2mix(index)
index=where(nomask ne -99. and no2mask eq -99.)     ; only NO good
if index(0) ne -1L then noxmix(index)=nomix(index)
index=where(nomask eq -99. and no2mask ne -99.)     ; only NO2 good
if index(0) ne -1L then noxmix(index)=no2mix(index)
onox4=noxmix
nday4=long(max(fdoy))-long(min(fdoy))+1L
fdoy4=fdoy
latitude4=latitude
index=where(latitude4 gt 50.)
if index(0) ne -1L then begin
   onox4=onox4(index,*)
   fdoy4=fdoy4(index)
   latitude4=latitude4(index)
endif

restore,dira+'dmps_ace_v2.2.meto.2005
restore,dira+'cat_ace_v2.2.2005
restore,dira+'no_ace_v2.2.2005
nomix=mix & nomask=mask
restore,dira+'no2_ace_v2.2.2005
no2mix=mix & no2mask=mask
noxmix=0.*nomix
index=where(nomask ne -99. and no2mask ne -99.)     ; both good
if index(0) ne -1L then noxmix(index)=nomix(index)+no2mix(index)
index=where(nomask ne -99. and no2mask eq -99.)     ; only NO good
if index(0) ne -1L then noxmix(index)=nomix(index)
index=where(nomask eq -99. and no2mask ne -99.)     ; only NO2 good
if index(0) ne -1L then noxmix(index)=no2mix(index)
onox5=noxmix
nday5=long(max(fdoy))-long(min(fdoy))+1L
fdoy5=fdoy
latitude5=latitude
index=where(latitude5 gt 50.)
if index(0) ne -1L then begin
   onox5=onox5(index,*)
   fdoy5=fdoy5(index)
   latitude5=latitude5(index)
endif

restore,dira+'dmps_ace_v2.2.meto.2006
restore,dira+'cat_ace_v2.2.2006
restore,dira+'no_ace_v2.2.2006
nomix=mix & nomask=mask
restore,dira+'no2_ace_v2.2.2006
no2mix=mix & no2mask=mask
noxmix=0.*nomix
index=where(nomask ne -99. and no2mask ne -99.)     ; both good
if index(0) ne -1L then noxmix(index)=nomix(index)+no2mix(index)
index=where(nomask ne -99. and no2mask eq -99.)     ; only NO good
if index(0) ne -1L then noxmix(index)=nomix(index)
index=where(nomask eq -99. and no2mask ne -99.)     ; only NO2 good
if index(0) ne -1L then noxmix(index)=no2mix(index)
onox6=noxmix
nday6=long(max(fdoy))-long(min(fdoy))+1L
fdoy6=fdoy
latitude6=latitude
index=where(latitude6 gt 50.)
if index(0) ne -1L then begin
   onox6=onox6(index,*)
   fdoy6=fdoy6(index)
   latitude6=latitude6(index)
endif

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
       device,/landscape,bits=8,filename='timeseries_ace_nox_spring_'+salt+'km.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
    endif
;
; extract vortex NOx at this altitude
;
    nox2004=reform(ONOX4(*,kk))*1.e9
    nox2005=reform(ONOX5(*,kk))*1.e9
    nox2006=reform(ONOX6(*,kk))*1.e9
    index=where(nox2004 gt 0.)
    if index(0) ne -1L then begin
       nox2004=nox2004(index)
       fdoy2004=fdoy4(index)
    endif
    index=where(nox2005 gt 0.)
    if index(0) ne -1L then begin
       nox2005=nox2005(index)
       fdoy2005=fdoy5(index)
    endif
    index=where(nox2006 gt 0.)
    if index(0) ne -1L then begin
       nox2006=nox2006(index)
       fdoy2006=fdoy6(index)
    endif
;print,altitude(kk),min(nox2004),max(nox2004),min(nox2005),max(nox2005),min(nox2006),max(nox2006)
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
kday=91
    nmin=0.
    nmax=max([nox2004,nox2005,nox2006])
    if nmax le 0. then goto,jumplev
if nmax gt 5000. then nmax=5000.
if altitude(kk) le 70. and nmax gt 1000. then nmax=1000.
print,altitude(kk),nmax
    if altitude(kk) gt 80. or altitude(kk) lt 30. then goto,jumplev
    plot,fdoy2004,nox2004,psym=8,/noeras,yrange=[nmin,nmax],xticks=3,$
         xtickname=['                            Jan',$
                    '                            Feb',$
                    '                            Mar','  '],$
         charsize=1.75,ytitle='NOx (ppbv)',color=0,xrange=[1.,kday],$
         title='ACE NOx > 50 N '+salt+' km'
    oplot,fdoy2004,nox2004,psym=8,color=mcolor*.9,symsize=2
;   oplot,fdoy2004,nox2004,psym=0,color=mcolor*.9,thick=2
    oplot,fdoy2005,nox2005,psym=8,color=mcolor*.45,symsize=2
;   oplot,fdoy2005,nox2005,psym=0,color=mcolor*.45,thick=2
    oplot,fdoy2006,nox2006,psym=8,color=0,symsize=2
;   oplot,fdoy2006,nox2006,psym=0,color=0,thick=2
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
       spawn,'convert -trim timeseries_ace_nox_spring_'+salt+'km.ps -rotate -90 '+$
             'timeseries_ace_nox_spring_'+salt+'km.jpg'
       spawn,'/usr/bin/rm timeseries_ace_nox_spring_'+salt+'km.ps'
    endif
    jumplev:
endfor	; loop over years
end
