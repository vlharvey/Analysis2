;
; scatter plot and PDFs of physically coincident 
; HIRDLS and HALOE/SAGE II/SAGE III/POAM III ozone
; VLH 11/12/2003
;
@aura2date
@loadauradata
@rd_sage3_o3_soundings
@rd_haloe_o3_soundings
@rd_poam3_o3_soundings
@rd_sage2_o3_soundings
@rd_ukmo_nc3
@stddat
@kgmt
@ckday
@kdate
@ks_stats

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
nlvls=25
col1=1+indgen(nlvls)*mcolor/nlvls
icmm1=icolmax-1
icmm2=icolmax-2
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.1,0.6,0.1,0.6,0.1,0.6,0.1,0.6]
yorig=[0.75,0.75,0.525,0.525,0.3,0.3,0.075,0.075]
xlen=0.30
ylen=0.14
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
SpeciesNames = ['Temperature',$
                'H2O', $
                'O3',  $
                'N2O', $
                'HNO3']
GeoLoc = ['Pressure',$
          'Time',$
          'Latitude',$
          'Longitude',$
          'SolarZenithAngle',$
          'LocalSolarTime']
hdir='/aura3/data/HIRDLS_data/Datfiles/'
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
dirh='/aura3/data/HALOE_data/Sound_data/haloe_'
dirs='/aura3/data/SAGE_II_data/Sound_data/sage2_'
dirs3='/aura3/data/SAGE_III_data/Sound_data/sage3_solar_'
dirp='/aura3/data/POAM_data/Sound_data/poam3_'
ifile='                             '
lstmn=10 & lstdy=2 & lstyr=0 & lstday=0
ledmn=10 & leddy=2 & ledyr=0 & ledday=0
;
; Ask interactive questions- get starting/ending date
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; --- Loop here --------
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
;
; read satellite ozone soundings
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      sfile=mon(imn-1)+sdy+'_'+syr
      rd_sage3_o3_soundings,dirs3+sfile+'_o3.sound',norbits3,tsage3,$
         xsage3,ysage3,tropps3,tropzs3,tropths3,modes3,o3sage3,psage3,$
         thsage3,zsage3,clsage3,qo3sage3,nlevs3
      print,norbits3,' SAGE III'
      rd_sage2_o3_soundings,dirs+sfile+'_o3.sound',norbits2,tsage2,$
         xsage2,ysage2,tropps2,tropzs2,tropths2,modes2,o3sage2,psage2,$
         thsage2,zsage2,clsage2,qo3sage2,nlevs2
      print,norbits2,' SAGE II'
      rd_poam3_o3_soundings,dirp+sfile+'_o3.sound',norbitp3,tpoam3,$
         xpoam3,ypoam3,troppp3,tropzp3,tropthp3,modep3,o3poam3,ppoam3,$
         thpoam3,zpoam3,clpoam3,qo3poam3,nlevp3
      print,norbitp3,' POAM III'
      rd_haloe_o3_soundings,dirh+sfile+'_o3.sound',norbith,thal,$
         xhal,yhal,tropph,tropzh,tropthh,modeh,o3hal,phal,$
         thhal,zhal,clhal,qo3hal,nlevh
      print,norbith,' HALOE'
;
; read HIRDLS and MLS data
;
      sday=strcompress(string(iday),/remove_all)
      Hfile=hdir+'HIRDLS2_'+syr+'d'+sday+'_MZ3_c1.he5'
      print,hfile

;; load HIRDLS data all at once
      hirdls=LoadAuraData(Hfile, [GeoLoc, SpeciesNames])

;; file header and tail for MLS
      Mfileh=hdir+'MLS-Aura_L2GP-'
      Mfilet='_sAura2c--t_'+syr+'d'+sday+'.he5'
      print,mfilet

;; loop over species
      FOR is = 0,N_ELEMENTS(SpeciesNames)-1 DO BEGIN
          SpeciesName = SpeciesNames(is)
          Mfile=Mfileh + SpeciesName + Mfilet
          IF is EQ 0 THEN mls = LoadAuraData(Mfile, GeoLoc)
          mls=LoadAuraData(Mfile, SpeciesName, mls)
      ENDFOR
;
; extract mls and hirdls variables
; time is elapsed seconds since midnight 1 Jan 1993
;
      mpress=mls.p            ; P               FLOAT     Array[37]
      mlev=n_elements(mpress)
      mtime=mls.time          ; TIME            DOUBLE    Array[3495]
      mlat=mls.lat            ; LAT             FLOAT     Array[3495]
      mlon=mls.lon            ; LON             FLOAT     Array[3495]
      msza=mls.sza            ; SZA             FLOAT     Array[3495]
      mlst=mls.lst            ; LST             FLOAT     Array[3495]
      mprof=n_elements(mlst)
      mtemp=mls.t             ; T               FLOAT     Array[37, 3495]
      mh2o=mls.h2o            ; H2O             FLOAT     Array[37, 3495]
      mo3=mls.o3              ; O3              FLOAT     Array[37, 3495]
      mn2o=mls.n2o            ; N2O             FLOAT     Array[37, 3495]
      mhno3=mls.hno3          ; HNO3            FLOAT     Array[37, 3495]

      hpress=hirdls.p         ;   P               FLOAT     Array[145]
      hlev=n_elements(hpress)
      htime=hirdls.time       ;   TIME            DOUBLE    Array[7848]
      hlat=hirdls.lat         ;   LAT             FLOAT     Array[7848]
      hlon=hirdls.lon         ;   LON             FLOAT     Array[7848]
      hsza=hirdls.sza         ;   SZA             FLOAT     Array[7848]
      hlst=hirdls.lst         ;   LST             FLOAT     Array[7848]
      hprof=n_elements(hlst)
      htemp=hirdls.t          ;   T               FLOAT     Array[145, 7848]
      hh2o=hirdls.h2o         ;   H2O             FLOAT     Array[145, 7848]
      ho3=hirdls.o3           ;   O3              FLOAT     Array[145, 7848]
      hn2o=hirdls.n2o         ;   N2O             FLOAT     Array[145, 7848]
      hno3=hirdls.hno3        ;   HNO3            FLOAT     Array[145, 7848]
;
; make press,lat,lon 2d
;
      mpress2=0.*mo3
      mlat2=0.*mo3
      mlon2=0.*mo3
      for i=0L,mprof-1L do mpress2(*,i)=mpress
      for i=0L,mlev-1L do begin
          mlat2(i,*)=mlat
          mlon2(i,*)=mlon
      endfor
      hpress2=0.*ho3
      hlat2=0.*ho3
      hlon2=0.*ho3
      for i=0L,hprof-1L do hpress2(*,i)=hpress
      for i=0L,hlev-1L do begin
          hlat2(i,*)=hlat
          hlon2(i,*)=hlon
      endfor
      mtheta=mtemp*(1000./mpress2)^0.286
      htheta=htemp*(1000./hpress2)^0.286
;
; retain coincident soundings
;
      if icount eq 0L then begin
         ncoin=1000L & nlev=300L & dxc=100.
         xcoinhal=-9999.+fltarr(ncoin,nlev)
         ycoinhal=-9999.+fltarr(ncoin,nlev)
         pcoinhal=-9999.+fltarr(ncoin,nlev)
         thcoinhal=-9999.+fltarr(ncoin,nlev)
         o3coinhal=-9999.+fltarr(ncoin,nlev)
         xcoinhirdlshal=-9999.+fltarr(ncoin)
         ycoinhirdlshal=-9999.+fltarr(ncoin)
         pcoinhirdlshal=-9999.+fltarr(ncoin,nlev)
         thcoinhirdlshal=-9999.+fltarr(ncoin,nlev)
         o3coinhirdlshal=-9999.+fltarr(ncoin,nlev)
         xcoinsage3=-9999.+fltarr(ncoin,nlev)
         ycoinsage3=-9999.+fltarr(ncoin,nlev)
         pcoinsage3=-9999.+fltarr(ncoin,nlev)
         thcoinsage3=-9999.+fltarr(ncoin,nlev)
         o3coinsage3=-9999.+fltarr(ncoin,nlev)
         xcoinhirdls3=-9999.+fltarr(ncoin)
         ycoinhirdls3=-9999.+fltarr(ncoin)
         pcoinhirdls3=-9999.+fltarr(ncoin,nlev)
         thcoinhirdls3=-9999.+fltarr(ncoin,nlev)
         o3coinhirdls3=-9999.+fltarr(ncoin,nlev)
         xcoinsage2=-9999.+fltarr(ncoin,nlev)
         ycoinsage2=-9999.+fltarr(ncoin,nlev)
         pcoinsage2=-9999.+fltarr(ncoin,nlev)
         thcoinsage2=-9999.+fltarr(ncoin,nlev)
         o3coinsage2=-9999.+fltarr(ncoin,nlev)
         xcoinhirdls2=-9999.+fltarr(ncoin)
         ycoinhirdls2=-9999.+fltarr(ncoin)
         pcoinhirdls2=-9999.+fltarr(ncoin,nlev)
         thcoinhirdls2=-9999.+fltarr(ncoin,nlev)
         o3coinhirdls2=-9999.+fltarr(ncoin,nlev)
         xcoinpoam3=-9999.+fltarr(ncoin,nlev)
         ycoinpoam3=-9999.+fltarr(ncoin,nlev)
         pcoinpoam3=-9999.+fltarr(ncoin,nlev)
         thcoinpoam3=-9999.+fltarr(ncoin,nlev)
         o3coinpoam3=-9999.+fltarr(ncoin,nlev)
         xcoinhirdlsp3=-9999.+fltarr(ncoin)
         ycoinhirdlsp3=-9999.+fltarr(ncoin)
         pcoinhirdlsp3=-9999.+fltarr(ncoin,nlev)
         thcoinhirdlsp3=-9999.+fltarr(ncoin,nlev)
         o3coinhirdlsp3=-9999.+fltarr(ncoin,nlev)
         hcoin=0L & scoin2=0L & scoin3=0L & pcoin3=0L
      endif
;
; find HIRDLS soundings within dxc km of HALOE/SAGE II/SAGE III/POAM III soundings 
;
; HIRDLS/HALOE
;
      if norbith gt 0L then begin
         for i=0,norbith-1L do begin
             index=where(thhal(i,*) gt 0. and thhal(i,*) ne 1.00000e+24)
             thmin=min(thhal(i,index)) & thmax=max(thhal(i,index))
             xh=xhal(i) & yh=yhal(i)
             dxf=re*abs(xh-hlon)*dtr*cos(yh*dtr)
             dyf=re*abs(yh-hlat)*dtr
             dist=sqrt(dxf*dxf+dyf*dyf)
             hindex=where(dist le dxc,ncoin0)
             if hindex(0) ne -1 then begin
                for icoin=0L,ncoin0-1L do begin
                    ii=hindex(icoin)
                    xcoinhirdlshal(hcoin,*)=hlon(ii)
                    ycoinhirdlshal(hcoin,*)=hlat(ii)
                    pcoinhirdlshal(hcoin,0:hlev-1L)=hpress2(*,ii)
                    thcoinhirdlshal(hcoin,0:hlev-1L)=htheta(*,ii)
                    o3coinhirdlshal(hcoin,0:hlev-1L)=ho3(*,ii)
                    xcoinhal(hcoin,*)=xh
                    ycoinhal(hcoin,*)=yh
                    for kk=0,hlev-1 do begin
                        thlev=htheta(kk,ii)
                        if thlev lt thmin or thlev gt thmax then goto,jumpkk1
                        for k=0,nlevh-2 do begin
                            if thhal(i,k) ne 1.00000e+24 and thlev lt thhal(i,k) then goto,jumpkk1
                            if o3hal(i,k) lt 0. or o3hal(i,k+1) lt 0. then goto,jumplev1
                            if o3hal(i,k) gt 1. or o3hal(i,k+1) gt 1. then goto,jumplev1
                            if thhal(i,k) gt 3000. then goto,jumpkk1
                            p0=phal(i,k) & p1=phal(i,k+1)
                            th0=thhal(i,k) & th1=thhal(i,k+1)
                            o30=o3hal(i,k) & o31=o3hal(i,k+1)
                            if th0 le thlev and th1 gt thlev then begin
                               scale=(th1-thlev)/(th1-th0)
                               thcoinhal(hcoin,kk)=th1-scale*(th1-th0)
                               o3coinhal(hcoin,kk)=o31-scale*(o31-o30)
                               pcoinhal(hcoin,kk)=p1^.286-scale*(p1^.286-p0^.286)
                               pcoinhal(hcoin,kk)=pcoinhal(hcoin,kk)^(1./.286)
                               goto,jumpkk1
                            endif
                            jumplev1:
                        endfor
                        jumpkk1:
                    endfor
                    hcoin=hcoin+1L
                    if hcoin ge ncoin then stop,'increase ncoin'
                endfor
             endif
         endfor
         print,'HALOE done'
      endif
;
; HIRDLS/SAGE III
;
      if norbits3 gt 0L then begin
         for i=0,norbits3-1L do begin
             index=where(thsage3(i,*) gt 0. and thsage3(i,*) ne 1.00000e+24)
             thmin=min(thsage3(i,index)) & thmax=max(thsage3(i,index))
             xh=xsage3(i) & yh=ysage3(i)
             dxf=re*abs(xh-hlon)*dtr*cos(yh*dtr)
             dyf=re*abs(yh-hlat)*dtr
             dist=sqrt(dxf*dxf+dyf*dyf)
             hindex=where(dist le dxc,ncoin0)
             if hindex(0) ne -1 then begin
                for icoin=0L,ncoin0-1L do begin
                    ii=hindex(icoin)
                    xcoinhirdls3(scoin3,*)=hlon(ii)
                    ycoinhirdls3(scoin3,*)=hlat(ii)
                    pcoinhirdls3(scoin3,0:hlev-1L)=hpress2(*,ii)
                    thcoinhirdls3(scoin3,0:hlev-1L)=htheta(*,ii)
                    o3coinhirdls3(scoin3,0:hlev-1L)=ho3(*,ii)
                    xcoinsage3(scoin3,*)=xh
                    ycoinsage3(scoin3,*)=yh
                    for kk=0,hlev-1 do begin
                        thlev=htheta(kk,ii)
                        if thlev lt thmin or thlev gt thmax then goto,jumpkk2
                        for k=0,nlevs3-2 do begin
                            if thsage3(i,k) ne 1.00000e+24 and thlev lt thsage3(i,k) then goto,jumpkk2
                            if o3sage3(i,k) lt 0. or o3sage3(i,k+1) lt 0. then goto,jumplev2
                            if o3sage3(i,k) gt 1. or o3sage3(i,k+1) gt 1. then goto,jumplev2
                            if thsage3(i,k) gt 3000. then goto,jumpkk2

                            p0=psage3(i,k) & p1=psage3(i,k+1)
                            th0=thsage3(i,k) & th1=thsage3(i,k+1)
                            o30=o3sage3(i,k) & o31=o3sage3(i,k+1)
                            if th0 le thlev and th1 gt thlev then begin
                               scale=(th1-thlev)/(th1-th0)
                               thcoinsage3(scoin3,kk)=th1-scale*(th1-th0)
                               o3coinsage3(scoin3,kk)=o31-scale*(o31-o30)
                               pcoinsage3(scoin3,kk)=p1^.286-scale*(p1^.286-p0^.286)
                               pcoinsage3(scoin3,kk)=pcoinsage3(scoin3,kk)^(1./.286)
                               goto,jumpkk2
                            endif
                            jumplev2:
                        endfor
                        jumpkk2:
                    endfor
                    scoin3=scoin3+1L
                    if scoin3 ge ncoin then stop,'increase ncoin'
                endfor
             endif
         endfor
         print,'SAGE III done'
      endif
;
; HIRDLS/SAGE II
;
      if norbits2 gt 0L then begin
         for i=0,norbits2-1L do begin
             index=where(thsage2(i,*) gt 0. and thsage2(i,*) ne 1.00000e+24)
             thmin=min(thsage2(i,index)) & thmax=max(thsage2(i,index))
             xh=xsage2(i) & yh=ysage2(i)
             dxf=re*abs(xh-hlon)*dtr*cos(yh*dtr)
             dyf=re*abs(yh-hlat)*dtr
             dist=sqrt(dxf*dxf+dyf*dyf)
             hindex=where(dist le dxc,ncoin0)
             if hindex(0) ne -1 then begin
                for icoin=0L,ncoin0-1L do begin
                    ii=hindex(icoin)
                    xcoinhirdls2(scoin2,*)=hlon(ii)
                    ycoinhirdls2(scoin2,*)=hlat(ii)
                    pcoinhirdls2(scoin2,0:hlev-1L)=hpress2(*,ii)
                    thcoinhirdls2(scoin2,0:hlev-1L)=htheta(*,ii)
                    o3coinhirdls2(scoin2,0:hlev-1L)=ho3(*,ii)
                    xcoinsage2(scoin2,*)=xh
                    ycoinsage2(scoin2,*)=yh
                    for kk=0,hlev-1 do begin
                        thlev=htheta(kk,ii)
                        if thlev lt thmin or thlev gt thmax then goto,jumpkk3
                        for k=0,nlevs2-2 do begin
                            if thsage2(i,k) ne 1.00000e+24 and thlev lt thsage2(i,k) then goto,jumpkk3
                            if o3sage2(i,k) lt 0. or o3sage2(i,k+1) lt 0. then goto,jumplev3
                            if o3sage2(i,k) gt 1. or o3sage2(i,k+1) gt 1. then goto,jumplev3
                            if thsage2(i,k) gt 3000. then goto,jumpkk3

                            p0=psage2(i,k) & p1=psage2(i,k+1)
                            th0=thsage2(i,k) & th1=thsage2(i,k+1)
                            o30=o3sage2(i,k) & o31=o3sage2(i,k+1)
                            if th0 le thlev and th1 gt thlev then begin
                               scale=(th1-thlev)/(th1-th0)
                               thcoinsage2(scoin2,kk)=th1-scale*(th1-th0)
                               o3coinsage2(scoin2,kk)=o31-scale*(o31-o30)
                               pcoinsage2(scoin2,kk)=p1^.286-scale*(p1^.286-p0^.286)
                               pcoinsage2(scoin2,kk)=pcoinsage2(scoin2,kk)^(1./.286)
                               goto,jumpkk3
                            endif
                            jumplev3:
                        endfor
                        jumpkk3:
                    endfor
                    scoin2=scoin2+1L
                    if scoin2 ge ncoin then stop,'increase ncoin'
                endfor
             endif
         endfor
         print,'SAGE II done'
      endif
;
; HIRDLS/POAM III
;
      if norbitp3 gt 0L then begin
         for i=0,norbitp3-1L do begin
             index=where(thpoam3(i,*) gt 0. and thpoam3(i,*) ne 1.00000e+24)
             thmin=min(thpoam3(i,index)) & thmax=max(thpoam3(i,index))
             xh=xpoam3(i) & yh=ypoam3(i)
             dxf=re*abs(xh-hlon)*dtr*cos(yh*dtr)
             dyf=re*abs(yh-hlat)*dtr
             dist=sqrt(dxf*dxf+dyf*dyf)
             hindex=where(dist le dxc,ncoin0)
             if hindex(0) ne -1 then begin
                for icoin=0L,ncoin0-1L do begin
                    ii=hindex(icoin)
                    xcoinhirdlsp3(pcoin3)=hlon(ii)
                    ycoinhirdlsp3(pcoin3)=hlat(ii)
                    pcoinhirdlsp3(pcoin3,0:hlev-1L)=hpress2(*,ii)
                    thcoinhirdlsp3(pcoin3,0:hlev-1L)=htheta(*,ii)
                    o3coinhirdlsp3(pcoin3,0:hlev-1L)=ho3(*,ii)
                    xcoinpoam3(pcoin3,*)=xh
                    ycoinpoam3(pcoin3,*)=yh
                    for kk=0,hlev-1 do begin
                        thlev=htheta(kk,ii)
                        if thlev lt thmin or thlev gt thmax then goto,jumpkk4
                        for k=0,nlevp3-2 do begin
                            if thpoam3(i,k) ne 1.00000e+24 and thlev lt thpoam3(i,k) then goto,jumpkk4
                            if o3poam3(i,k) lt 0. or o3poam3(i,k+1) lt 0. then goto,jumplev4
                            if o3poam3(i,k) gt 1. or o3poam3(i,k+1) gt 1. then goto,jumplev4
                            if thpoam3(i,k) gt 3000. then goto,jumpkk4

                            p0=ppoam3(i,k) & p1=ppoam3(i,k+1)
                            th0=thpoam3(i,k) & th1=thpoam3(i,k+1)
                            o30=o3poam3(i,k) & o31=o3poam3(i,k+1)
                            if th0 le thlev and th1 gt thlev then begin
                               scale=(th1-thlev)/(th1-th0)
                               thcoinpoam3(pcoin3,kk)=th1-scale*(th1-th0)
                               o3coinpoam3(pcoin3,kk)=o31-scale*(o31-o30)
                               pcoinpoam3(pcoin3,kk)=p1^.286-scale*(p1^.286-p0^.286)
                               pcoinpoam3(pcoin3,kk)=pcoinpoam3(pcoin3,kk)^(1./.286)
                               goto,jumpkk4
                            endif
                            jumplev4:
                        endfor
                        jumpkk4:
                    endfor
                    pcoin3=pcoin3+1L
                    if pcoin3 ge ncoin then stop,'increase ncoin'
                endfor
             endif
         endfor
         print,'POAM III done'
      endif
      icount=icount+1L
goto,jump
           
plotit:
;
; remove data voids
;
if hcoin gt 0L then begin
   index=where(o3coinhal gt 0. and o3coinhirdlshal gt 0.,hcoin) 
   xcoinhirdlshal=xcoinhirdlshal(index)
   ycoinhirdlshal=ycoinhirdlshal(index)
   pcoinhirdlshal=pcoinhirdlshal(index)
   thcoinhirdlshal=thcoinhirdlshal(index)
   o3coinhirdlshal=o3coinhirdlshal(index)*1.e6
   xcoinhal=xcoinhal(index)
   ycoinhal=ycoinhal(index)
   pcoinhal=pcoinhal(index)
   thcoinhal=thcoinhal(index)
   o3coinhal=o3coinhal(index)*1.e6
endif

if scoin3 gt 0L then begin
   index=where(o3coinsage3 gt 0. and o3coinhirdls3 gt 0.,scoin3)
   xcoinhirdls3=xcoinhirdls3(index)
   ycoinhirdls3=ycoinhirdls3(index)
   pcoinhirdls3=pcoinhirdls3(index)
   thcoinhirdls3=thcoinhirdls3(index)
   o3coinhirdls3=o3coinhirdls3(index)*1.e6
   xcoinsage3=xcoinsage3(index)
   ycoinsage3=ycoinsage3(index)
   pcoinsage3=pcoinsage3(index)
   thcoinsage3=thcoinsage3(index)
   o3coinsage3=o3coinsage3(index)*1.e6
endif

if scoin2 gt 0L then begin
   index=where(o3coinsage2 gt 0. and o3coinhirdls2 gt 0.,scoin2)
   xcoinhirdls2=xcoinhirdls2(index)
   ycoinhirdls2=ycoinhirdls2(index)
   pcoinhirdls2=pcoinhirdls2(index)
   thcoinhirdls2=thcoinhirdls2(index)
   o3coinhirdls2=o3coinhirdls2(index)*1.e6
   xcoinsage2=xcoinsage2(index)
   ycoinsage2=ycoinsage2(index)
   pcoinsage2=pcoinsage2(index)
   thcoinsage2=thcoinsage2(index)
   o3coinsage2=o3coinsage2(index)*1.e6
endif

if pcoin3 gt 0L then begin
   index=where(o3coinpoam3 gt 0. and o3coinhirdlsp3 gt 0.,pcoin3)
   xcoinhirdlsp3=xcoinhirdlsp3(index)
   ycoinhirdlsp3=ycoinhirdlsp3(index)
   pcoinhirdlsp3=pcoinhirdlsp3(index)
   thcoinhirdlsp3=thcoinhirdlsp3(index)
   o3coinhirdlsp3=o3coinhirdlsp3(index)*1.e6
   xcoinpoam3=xcoinpoam3(index)
   ycoinpoam3=ycoinpoam3(index)
   pcoinpoam3=pcoinpoam3(index)
   thcoinpoam3=thcoinpoam3(index)
   o3coinpoam3=o3coinpoam3(index)*1.e6
endif

daterange=strcompress(string(FORMAT='(A3,A1,I2,A2,I4,A3,A3,A1,I2,A2,I4,A3)',$
 month(lstmn-1),' ',lstdy,', ',lstyr,' - ',month(ledmn-1),' ',leddy,', ',ledyr))
datelab=strcompress(string(FORMAT='(I4,I2.2,I2.2,A1,I4,I2.2,I2.2)',$
 lstyr,lstmn,lstdy,'-',ledyr,ledmn,leddy))

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,$
           filename='scatter_hirdls_coin_occul_o3_'+datelab+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
           xsize=xsize,ysize=ysize
endif

; Set plot boundaries
erase
!type=2^2+2^3
xyouts,.25,.95,'Ozone '+daterange,/normal,charsize=2
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,findgen(13),findgen(13),xrange=[0.,12.],yrange=[0.,12.],charsize=1.2,$
     ytitle='HALOE',xtitle='HIRDLS',title='Scatter Plot'
if hcoin gt 0. then begin
result=correlate(o3coinhirdlshal,o3coinhal)
r=result(0)
xyouts,8.,3.,'N ='+strcompress(string(hcoin)),/data,charsize=1.2
xyouts,8.,1.,'r = '+strcompress(string(format='(f6.3)',r)),/data,charsize=1.2
thmax=3000. & thmin=200.
for icoin=0L,hcoin-1L do begin
   a=findgen(8)*(2*!pi/8.)
   usersym,cos(a),sin(a),/fill
   xx=o3coinhirdlshal(icoin)
   yy=o3coinhal(icoin)
   oplot,[xx,xx],[yy,yy],psym=8,color=((thcoinhirdlshal(icoin)-thmin)/(thmax-thmin))*icolmax,symsize=.5
   a=findgen(9)*(2*!pi/8.)
   usersym,cos(a),sin(a)
;  oplot,[xx,xx],[yy,yy],psym=8,color=lc,symsize=.5
endfor
endif
omin=thmin & omax=thmax
xmnb=xmn+xlen+0.01
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax]
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
axis,10,omin,0,YAX=1,/DATA,charsize=1.2,/ynozero
xyouts,xmxb,ymx+0.01,'Theta',/normal,charsize=1.2

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
x=.5*findgen(nlvls)
if hcoin gt 0L then begin
y1=histogram(o3coinhirdlshal,min=0,max=12.,binsize=.5)/(1.*hcoin)
y1=smooth(y1,3)     ; mapped PDF
y2=histogram(o3coinhal,min=0,max=12.,binsize=.5)/(1.*hcoin)
y2=smooth(y2,3)     ; flight PDF
ymax=max(y1,y2)+0.2*max(y1,y2)
!linetype=0
plot,x,y1,xtitle='Ozone (ppmv)',ytitle='Frequency',charsize=1.2,$
     title='PDFs',xrange=[0.,12.],yrange=[0.,ymax]
y2=histogram(o3coinhal,min=0,max=12.,binsize=.5)/(1.*hcoin)
y2=smooth(y2,3)     ; flight PDF
!linetype=1
oplot,x,y2
plots,[7.5,.3*ymax],/data
plots,[8.5,.3*ymax],/continue,/data
!linetype=0
xyouts,8.75,.3*ymax,'HALOE',/data,charsize=1.2
y2=histogram(o3coinhirdlshal,min=0,max=12.,binsize=.5)/(1.*hcoin)
y2=smooth(y2,3)
plots,[7.5,.1*ymax],/data
plots,[8.5,.1*ymax],/continue,/data
xyouts,8.75,.1*ymax,'HIRDLS',/data,charsize=1.2
ks_stats,o3coinhal,o3coinhirdlshal,kstest,cprob
xyouts,7.5,.9*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.8*ymax,'KS sig='+$
       strcompress(string(format='(f5.3)',100.*cprob),/remove_all)+'%',/data
print,'KS=',kstest
print,'KS significance=',100.*cprob,' %'
so3bar=total(o3coinhirdlshal)/n_elements(o3coinhirdlshal)
ho3bar=total(o3coinhal)/n_elements(o3coinhal)
o3coinhirdlshal2=o3coinhirdlshal-so3bar+ho3bar
ks_stats,o3coinhal,o3coinhirdlshal2,kstest,cprob
xyouts,7.5,.7*ymax,'w/o Mean Bias:',/data
xyouts,7.5,.6*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.5*ymax,'KS sig='+strcompress(string(format='(f5.3)',$
        100.*cprob),/remove_all)+'%',/data
endif

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
plot,findgen(13),findgen(13),xrange=[0.,12.],yrange=[0.,12.],charsize=1.2,$
     ytitle='SAGE III',xtitle='HIRDLS',title='Scatter Plot'
thmax=3000. & thmin=200.
if scoin3 gt 0. then begin
result=correlate(o3coinhirdls3,o3coinsage3)
r=result(0)
xyouts,8.,3.,'N ='+strcompress(string(scoin3)),/data,charsize=1.2
xyouts,8.,1.,'r = '+strcompress(string(format='(f6.3)',r)),/data,charsize=1.2
for icoin=0L,scoin3-1L do begin
   a=findgen(8)*(2*!pi/8.)
   usersym,cos(a),sin(a),/fill
   xx=o3coinhirdls3(icoin)
   yy=o3coinsage3(icoin)
   oplot,[xx,xx],[yy,yy],psym=8,color=((thcoinhirdls3(icoin)-thmin)/(thmax-thmin))*icolmax,symsize=.5
   a=findgen(9)*(2*!pi/8.)
   usersym,cos(a),sin(a)
;  oplot,[xx,xx],[yy,yy],psym=8,color=lc,symsize=.5
endfor
endif
omin=thmin & omax=thmax
xmnb=xmn+xlen+0.01
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax]
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
axis,10,omin,0,YAX=1,/DATA,charsize=1.2,/ynozero
xyouts,xmxb,ymx+0.01,'Theta',/normal,charsize=1.2

xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
x=.5*findgen(nlvls)
if scoin3 gt 0L then begin
y1=histogram(o3coinhirdls3,min=0,max=12.,binsize=.5)/(1.*scoin3)
y1=smooth(y1,3)     ; mapped PDF
y2=histogram(o3coinsage3,min=0,max=12.,binsize=.5)/(1.*scoin3)
y2=smooth(y2,3)     ; flight PDF
ymax=max(y1,y2)+0.2*max(y1,y2)
!linetype=0
plot,x,y1,xtitle='Ozone (ppmv)',ytitle='Frequency',charsize=1.2,$
     title='PDFs',xrange=[0.,12.],yrange=[0.,ymax]
y2=histogram(o3coinsage3,min=0,max=12.,binsize=.5)/(1.*scoin3)
y2=smooth(y2,3)     ; flight PDF
!linetype=1
oplot,x,y2
plots,[7.5,.3*ymax],/data
plots,[8.5,.3*ymax],/continue,/data
!linetype=0
xyouts,8.75,.3*ymax,'SAGE III',/data,charsize=1.2
y2=histogram(o3coinhirdls3,min=0,max=12.,binsize=.5)/(1.*scoin3)
y2=smooth(y2,3)
plots,[7.5,.1*ymax],/data
plots,[8.5,.1*ymax],/continue,/data
xyouts,8.75,.1*ymax,'HIRDLS',/data,charsize=1.2
ks_stats,o3coinsage3,o3coinhirdls3,kstest,cprob
xyouts,7.5,.9*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.8*ymax,'KS sig='+$
       strcompress(string(format='(f5.3)',100.*cprob),/remove_all)+'%',/data
print,'KS=',kstest
print,'KS significance=',100.*cprob,' %'
so3bar=total(o3coinhirdls3)/n_elements(o3coinhirdls3)
ho3bar=total(o3coinsage3)/n_elements(o3coinsage3)
o3coinhirdls32=o3coinhirdls3-so3bar+ho3bar
ks_stats,o3coinsage3,o3coinhirdls32,kstest,cprob
xyouts,7.5,.7*ymax,'w/o Mean Bias:',/data
xyouts,7.5,.6*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.5*ymax,'KS sig='+strcompress(string(format='(f5.3)',$
        100.*cprob),/remove_all)+'%',/data
endif

xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
plot,findgen(13),findgen(13),xrange=[0.,12.],yrange=[0.,12.],charsize=1.2,$
     ytitle='SAGE II',xtitle='HIRDLS',title='Scatter Plot'
thmax=3000. & thmin=200.
if scoin2 gt 0. then begin
result=correlate(o3coinhirdls2,o3coinsage2)
r=result(0)
xyouts,8.,3.,'N ='+strcompress(string(scoin2)),/data,charsize=1.2
xyouts,8.,1.,'r = '+strcompress(string(format='(f6.3)',r)),/data,charsize=1.2
for icoin=0L,scoin2-1L do begin
   a=findgen(8)*(2*!pi/8.)
   usersym,cos(a),sin(a),/fill
   xx=o3coinhirdls2(icoin)
   yy=o3coinsage2(icoin)
   oplot,[xx,xx],[yy,yy],psym=8,color=((thcoinhirdls2(icoin)-thmin)/(thmax-thmin))*icolmax,symsize=.5
   a=findgen(9)*(2*!pi/8.)
   usersym,cos(a),sin(a)
;  oplot,[xx,xx],[yy,yy],psym=8,color=lc,symsize=.5
endfor
endif
omin=thmin & omax=thmax
xmnb=xmn+xlen+0.01
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,10],yrange=[omin,omax]
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
axis,10,omin,0,YAX=1,/DATA,charsize=1.2,/ynozero
xyouts,xmxb,ymx+0.01,'Theta',/normal,charsize=1.2

xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
x=.5*findgen(nlvls)
if scoin2 gt 0L then begin
y1=histogram(o3coinhirdls2,min=0,max=12.,binsize=.5)/(1.*scoin2)
y1=smooth(y1,3)     ; mapped PDF
y2=histogram(o3coinsage2,min=0,max=12.,binsize=.5)/(1.*scoin2)
y2=smooth(y2,3)     ; flight PDF
ymax=max(y1,y2)+0.2*max(y1,y2)
!linetype=0
plot,x,y1,xtitle='Ozone (ppmv)',ytitle='Frequency',charsize=1.2,$
     title='PDFs',xrange=[0.,12.],yrange=[0.,ymax]
y2=histogram(o3coinsage2,min=0,max=12.,binsize=.5)/(1.*scoin2)
y2=smooth(y2,3)     ; flight PDF
!linetype=1
oplot,x,y2
plots,[7.5,.3*ymax],/data
plots,[8.5,.3*ymax],/continue,/data
!linetype=0
xyouts,8.75,.3*ymax,'SAGE II',/data,charsize=1.2
y2=histogram(o3coinhirdls2,min=0,max=12.,binsize=.5)/(1.*scoin2)
y2=smooth(y2,3)
plots,[7.5,.1*ymax],/data
plots,[8.5,.1*ymax],/continue,/data
xyouts,8.75,.1*ymax,'HIRDLS',/data,charsize=1.2
ks_stats,o3coinsage2,o3coinhirdls2,kstest,cprob
xyouts,7.5,.9*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.8*ymax,'KS sig='+$
       strcompress(string(format='(f5.3)',100.*cprob),/remove_all)+'%',/data
print,'KS=',kstest
print,'KS significance=',100.*cprob,' %'
so3bar=total(o3coinhirdls2)/n_elements(o3coinhirdls2)
ho3bar=total(o3coinsage2)/n_elements(o3coinsage2)
o3coinhirdls22=o3coinhirdls2-so3bar+ho3bar
ks_stats,o3coinsage2,o3coinhirdls22,kstest,cprob
xyouts,7.5,.7*ymax,'w/o Mean Bias:',/data
xyouts,7.5,.6*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.5*ymax,'KS sig='+strcompress(string(format='(f5.3)',$
        100.*cprob),/remove_all)+'%',/data
endif

xmn=xorig(6)
xmx=xorig(6)+xlen
ymn=yorig(6)
ymx=yorig(6)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
plot,findgen(13),findgen(13),xrange=[0.,12.],yrange=[0.,12.],charsize=1.2,$
     ytitle='POAM III',xtitle='HIRDLS',title='Scatter Plot'
thmax=3000. & thmin=200.
if pcoin3 gt 0. then begin
result=correlate(o3coinhirdlsp3,o3coinpoam3)
r=result(0)
xyouts,8.,3.,'N ='+strcompress(string(pcoin3)),/data,charsize=1.2
xyouts,8.,1.,'r = '+strcompress(string(format='(f6.3)',r)),/data,charsize=1.2
for icoin=0L,pcoin3-1L do begin
   a=findgen(8)*(2*!pi/8.)
   usersym,cos(a),sin(a),/fill
   xx=o3coinhirdlsp3(icoin)
   yy=o3coinpoam3(icoin)
   oplot,[xx,xx],[yy,yy],psym=8,color=((thcoinhirdlsp3(icoin)-thmin)/(thmax-thmin))*icolmax,symsize=.5
   a=findgen(9)*(2*!pi/8.)
   usersym,cos(a),sin(a)
;  oplot,[xx,xx],[yy,yy],psym=8,color=lc,symsize=.5
endfor
endif
omin=thmin & omax=thmax
xmnb=xmn+xlen+0.01
xmxb=xmnb+cbarydel
set_viewport,xmnb,xmxb,ymn,ymx
!type=2^2+2^3+2^5+2^6
plot,[0,0],[omin,omax],xrange=[0,12],yrange=[omin,omax]
xbox=[0,10,10,0,0]
y1=omin
dy=(omax-omin)/float(nlvls)
for j=0,nlvls-1 do begin
    ybox=[y1,y1,y1+dy,y1+dy,y1]
    polyfill,xbox,ybox,color=col1(j)
    y1=y1+dy
endfor
!type=2^2+2^3+2^5
axis,10,omin,0,YAX=1,/DATA,charsize=1.2,/ynozero
xyouts,xmxb,ymx+0.01,'Theta',/normal,charsize=1.2

xmn=xorig(7)
xmx=xorig(7)+xlen
ymn=yorig(7)
ymx=yorig(7)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
x=.5*findgen(nlvls)
if pcoin3 gt 0L then begin
y1=histogram(o3coinhirdlsp3,min=0,max=12.,binsize=.5)/(1.*pcoin3)
y1=smooth(y1,3)     ; mapped PDF
y2=histogram(o3coinpoam3,min=0,max=12.,binsize=.5)/(1.*pcoin3)
y2=smooth(y2,3)     ; flight PDF
ymax=max(y1,y2)+0.2*max(y1,y2)
!linetype=0
plot,x,y1,xtitle='Ozone (ppmv)',ytitle='Frequency',charsize=1.2,$
     title='PDFs',xrange=[0.,12.],yrange=[0.,ymax]
y2=histogram(o3coinpoam3,min=0,max=12.,binsize=.5)/(1.*pcoin3)
y2=smooth(y2,3)     ; flight PDF
!linetype=1
oplot,x,y2
plots,[7.5,.3*ymax],/data
plots,[8.5,.3*ymax],/continue,/data
!linetype=0
xyouts,8.75,.3*ymax,'POAM III',/data,charsize=1.2
y2=histogram(o3coinhirdlsp3,min=0,max=12.,binsize=.5)/(1.*pcoin3)
y2=smooth(y2,3)
plots,[7.5,.1*ymax],/data
plots,[8.5,.1*ymax],/continue,/data
xyouts,8.75,.1*ymax,'HIRDLS',/data,charsize=1.2
ks_stats,o3coinpoam3,o3coinhirdlsp3,kstest,cprob
xyouts,7.5,.9*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.8*ymax,'KS sig='+$
       strcompress(string(format='(f5.3)',100.*cprob),/remove_all)+'%',/data
print,'KS=',kstest
print,'KS significance=',100.*cprob,' %'
so3bar=total(o3coinhirdlsp3)/n_elements(o3coinhirdlsp3)
ho3bar=total(o3coinpoam3)/n_elements(o3coinpoam3)
o3coinhirdlsp32=o3coinhirdlsp3-so3bar+ho3bar
ks_stats,o3coinpoam3,o3coinhirdlsp32,kstest,cprob
xyouts,7.5,.7*ymax,'w/o Mean Bias:',/data
xyouts,7.5,.6*ymax,'KS='+strmid(string(kstest),5,4),/data
xyouts,7.5,.5*ymax,'KS sig='+strcompress(string(format='(f5.3)',$
        100.*cprob),/remove_all)+'%',/data
endif

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert scatter_hirdls_coin_occul_o3_'+datelab+$
         '.ps -rotate -90 scatter_hirdls_coin_occul_o3_'+datelab+'.jpg'
endif
end
