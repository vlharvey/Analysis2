loadct,39
nlvls=30
mcolor=byte(!p.color)
col1=1+indgen(nlvls+1)*mcolor/float(nlvls)
!type=2^2+2^3
device,decompose=0

R = 287.
g0 = 9.81
Cp = 1004.
Pressure = [.015,0.02,.03,.04,.05,.06,.07,.08,.1,.15,.2,.3,.4,.5,.6,.7,.8,1.,1.5,2.,3.,4.,5.,6.,7.,8.,10.,15.,20.,30.,40.,50.,60.,70.,80.,100.,150.,200.,300.,400.,500.]
np=n_elements(pressure)

dirw='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_'
dum1=file_search(dirw+'????????.sav',count=count1)

for ifile=0L,count1-1L do begin
;
; extract date
;
    result=strsplit(dum1(ifile),'/',/extract)
    result2=strsplit(result(4),'.',/extract)
    result3=strsplit(result2(0),'_',/extract)
    sdate=result3(1)
    print,sdate
;
; skip if press file already exists
;
    dum=findfile(dirw+'press_'+sdate+'.sav')
    if dum(0) ne '' then goto,jumpfile
;
; restore T, U, V, DP, QV, PS
;
    restore,dum1(ifile)
    nc=n_elements(tnew[*,0,0])
    nr=n_elements(tnew[0,*,0])
    nz=n_elements(tnew[0,0,*])
    Tsum = fltarr(nc,nr,nz)
    lnPsum = fltarr(nc,nr,nz)
    tlayermean = fltarr(nc,nr,nz)
    GPHT_on_P = fltarr(nc, nr, np)
    QVGRD=0.*GPHT_on_P
    TGRD=0.*GPHT_on_P
    UGRD=0.*GPHT_on_P
    VGRD=0.*GPHT_on_P

    lnP = alog(levels_assim)	; 1d Pressure array
    lnPs = alog(1000.)		; 2d Surface Pressure
;
; arrays go top-down. nz-1 is surface
;
;tsum[*,*,nz-1L] = tnew[*,*,nz- 1L]*(alog(1000.) - lnP[nz- 1L])
;lnPsum[*,*,nz- 1L] = alog(1000.) - lnP[nz - 1L]

    lndp = alog(dpnew/100.)
    dp = dpnew/100.

    press = fltarr(nc,nr,nz+1L)
    gridboxThick = fltarr(nc,nr,nz)
    GPHT = fltarr(nc,nr,nz)
    GPHT1 = fltarr(nc,nr,nz+1L)
    GPHT1[*,*,nz] = 0.
    P = fltarr(nc,nr,nz)
    P1 = fltarr(nc,nr,nz+1L)
    P1[*,*,nz] = 1013.
    press[*,*,1L] = exp((alog(levels_assim[0L]) + alog(levels_assim[1L]))/2.)
    press[*,*,0L] = press[*,*,1L] - dp[*,*,0L]
    P[*,*,0L] = Pressure(0)	;levels_assim[0L]
    for kk = 2L, nz do begin
	press[*,*,kk] = press[*,*,kk-1L] + dp[*,*,kk-1L] 
	P[*,*,kk-1L] = exp((alog(press[*,*,kk-1L]) + alog(press[*,*,kk]))/2.)
    endfor

    SurfacedP = psnew - reform(press[*,*,nz-1L])

    theta = reform(tnew[*,*,nz-1L])*((1013./reform(press[*,*,nz]))^(R/Cp))
    tlayer = (reform(tnew[*,*,nz-1L]) * alog(reform(p[*,*,nz-1L])) + (theta * alog(1013.))) / (alog(reform(p[*,*,nz-1L])) + alog(1013.))
    bottomthick = ((R/g0)*tlayer*alog(1013./reform(press[*,*,nz])))
    for ii = 0L, nc - 1L do begin
	for jj = 0L, nr - 1L do begin
		for kk = nz-1L, 0L, -1L do begin
			gridboxThick[ii,jj,kk] = (R/g0)*tnew[ii,jj,kk]*alog(press[ii,jj,kk+1L]/press[ii,jj,kk])
			if psnew[ii,jj] ge 1013. then GPHT[ii,jj,kk] = total(gridboxThick[ii,jj,kk:nz-1L])
			if psnew[ii,jj] lt 1013. then GPHT[ii,jj,kk] = total(gridboxThick[ii,jj,kk:nz-1L]) + bottomthick[ii,jj]
		endfor

		P1[ii,jj,0:nz-1L] = P[ii,jj,*]
                GPHT1[ii,jj,0:nz-1L] = GPHT[ii,jj,*]
		GPHT_on_P[ii,jj,*] = interpol(GPHT1[ii,jj,*],alog(P1[ii,jj,*]),alog(Pressure))

                TGRD[ii,jj,*] = interpol(tnew[ii,jj,*],alog(P[ii,jj,*]),alog(Pressure))
                UGRD[ii,jj,*] = interpol(unew[ii,jj,*],alog(P[ii,jj,*]),alog(Pressure))
                VGRD[ii,jj,*] = interpol(vnew[ii,jj,*],alog(P[ii,jj,*]),alog(Pressure))
                QVGRD[ii,jj,*] = interpol(qvnew[ii,jj,*],alog(P[ii,jj,*]),alog(Pressure))

	endfor
    endfor
    psgrd=psnew
    ZGRD = GPHT_on_P/1000.
;
; check
;
;   for k = 0, np - 1L do begin
;      erase
;      z2d=reform(ZGRD(*,*,k))
;      t2d=reform(TGRD(*,*,k))
;      qv2d=reform(QVGRD(*,*,k))
;      map_set, 0,-180,/contin, title=sdate+' '+string(Pressure(k))+'hPa'
;      contour,z2d,LONGITUDE_WACCM,LATITUDE_WACCM,nlevels=30,/cell_fill,c_color=col1,/overplot
;      contour,t2d,LONGITUDE_WACCM,LATITUDE_WACCM,nlevels=20,thick=2,/follow,color=mcolor,/overplot
;      contour,qv2d,LONGITUDE_WACCM,LATITUDE_WACCM,nlevels=20,thick=2,/follow,color=0,/overplot
;      MAP_CONTINENTS, COLOR=white,/overplot
;      stop
;   endfor

save,filename=dirw+'press_'+sdate+'.sav',LATITUDE_WACCM,LONGITUDE_WACCM,PRESSURE,PSGRD,QVGRD,TGRD,UGRD,VGRD,ZGRD
jumpfile:
endfor	; loop over files

end
