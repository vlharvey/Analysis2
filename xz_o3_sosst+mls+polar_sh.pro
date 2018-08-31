;
; SH
; plot SOSST ozone in longitude altitude section
; add polar projection at a user prompted theta
; VLH 2/6/05
;
; swaths of MLS added	VLH 5/19/06
;
@aura2date
@rd_ukmo_nc3
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,38
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='x'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.10]
yorig=[0.10]
cbaryoff=0.015
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=0
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
nmon=['01','02','03','04','05','06','07','08','09','10','11','12']
dirs2='/aura3/data/SAGE_II_data/Datfiles_SOSST/'
dirs3='/aura3/data/SAGE_III_data/Datfiles_SOSST/'
diri='/aura3/data/ILAS_data/Datfiles_SOSST/'
dirh='/aura3/data/HALOE_data/Datfiles_SOSST/'
dirp='/aura3/data/POAM_data/Datfiles_SOSST/'
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
dirm='/aura6/data/MLS_data/Datfiles/'
ifile='                             '
lstmn=7 & lstdy=16 & lstyr=2005 & lstday=0
ledmn=5 & leddy=31 & ledyr=2006 & ledday=0
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
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
      if ndays gt ledday then stop,' Normal termination condition '
;
; read UKMO data
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      smn=string(FORMAT='(i2.2)',imn)
      idate=long(smn+sdy)
      uyr=strmid(syr,2,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
      if iflag eq 1 then goto,jump
      if icount eq 0L then begin
         rtheta=1000.
;        print,th
;        read,' Enter theta level ',rtheta
         index=where(rtheta eq th)
         if index(0) eq -1 then stop
         itheta=index(0)
         stheta=strcompress(string(fix(th(itheta))),/remove_all)
         xz2d=fltarr(nc,nth)
         yz2d=fltarr(nc,nth)
         for k=0,nth-1 do xz2d(*,k)=alon
         for j=0,nc-1 do yz2d(j,*)=th
         x=fltarr(nc+1)
         x(0:nc-1)=alon
         x(nc)=alon(0)+360.
         x2d=fltarr(nc+1,nr)
         y2d=fltarr(nc+1,nr)
         for i=0,nc do y2d(i,*)=alat
         for j=0,nr-1 do x2d(*,j)=x
      endif
;
; read MLS data
;
dum=findfile(dirm+'MLS-Aura_L2GP_'+syr+smn+sdy+'.sav')
if dum(0) eq '' then begin
   print,'No MLS data today'
   goto,jump
endif
restore,dirm+'MLS-Aura_L2GP_'+syr+smn+sdy+'.sav'
mpress=o3.pressure
mlev=n_elements(mpress)
mprof=o3.ntimes
mtime=o3.TIME
mlat=o3.latitude
mlon=o3.longitude
index=where(mlon lt 0.)
mlon(index)=mlon(index)+360.
mtemp=tp.L2GPVALUE
mo3=o3.L2GPVALUE*1.e6
mo3prec=o3.L2GPPRECISION
mo3stat=o3.STATUS
mo3qual=o3.QUALITY
mo3mask=0.*mo3
;
; remove profiles with negative times
;
xx=where(mtime gt 0.,mprof)
mtime=mtime(xx)
mlat=mlat(xx)
mlon=mlon(xx)
mo3=mo3(*,xx)
mo3stat=mo3stat(xx)
mo3qual=mo3qual(xx)
mo3prec=mo3prec(*,xx)
mo3mask=mo3mask(*,xx)
mtemp=mtemp(*,xx)
;
; use quality, status, and precision flags to remove suspect data
;
o3bad=where(mo3prec lt 0.)
if o3bad(0) ne -1L then mo3mask(o3bad)=-99.
o3bad=where(mo3stat mod 2 ne 0L)                ; o3status=0 is good, all odd values are bad
if o3bad(0) ne -1L then mo3mask(*,o3bad)=-99.
o3bad=where(mo3qual lt 0.1)                     ; do not use if quality < 0.1
if o3bad(0) ne -1L then mo3mask(*,o3bad)=-99.
;
; convert elapsed seconds to dates (yyyymmddhh)
;
aura2date,mdate,mtime
;
; time is elapsed seconds since midnight 1 Jan 1993 to hours today
; convert to daily UT time (0-24 hours)
;
muttime=mtime
istime=1993010100L
ehr=muttime/60./60.       ; convert time from seconds to hours
hh2=0.d*muttime
for n=0L,mprof-1L do begin
    yy1=istime/1000000
    if yy1 mod 4 eq 0 then mno(1)=29L
    if yy1 mod 4 ne 0 then mno(1)=28L
    mm1=istime/10000L-yy1*100L
    dd1=istime/100L-yy1*10000L-mm1*100L
    dd2=dd1+long(ehr(n))/24L
    hh1=istime-yy1*1000000L-mm1*10000L-dd1*100L
    yy2=yy1 & mm2=mm1
    while dd2 gt mno(mm2-1) do begin
          dd2=dd2-mno(mm2-1)
          mm2=mm2+1L
          if mm2 gt 12L then begin
             mm2=mm2-12L
             yy2=yy2+1L
             if yy2 mod 4 eq 0 then mno(1)=29
             if yy2 mod 4 ne 0 then mno(1)=28
          endif
    endwhile
    hh2(n)=ehr(n) mod 24
    if hh2(n) ge 24. then begin
       hh2(n)=hh2(n)-24.
       dd2=dd2+1L
       if dd2 gt mno(mm2-1L) then begin
          dd2=dd2-mno(mm2-1L)
          mm2=mm2+1L
          if mm2 gt 12L then begin
             mm2=mm2-12L
             yy2=yy2+1L
          endif
       endif
    endif
endfor
muttime=hh2
;
; calculate potential temperature
;
mtheta=0.*mtemp
for i=0L,mlev-1L do mtheta(i,*)=mtemp(i,*)*(1000./mpress(i))^0.286
index=where(mtemp lt 0.)
if index(0) ne -1L then mtheta(index)=-999.
mpress2=0.*mo3
mtime2=0.*mo3
mlat2=0.*mo3
mlon2=0.*mo3
for i=0L,mprof-1L do mpress2(*,i)=mpress
for i=0L,mlev-1L do begin
    mtime2(i,*)=muttime
    mlat2(i,*)=mlat
    mlon2(i,*)=mlon
endfor
;
; restore SOSST yearly files on 1 January
;
      if iday eq 1L or icount eq 0L then begin
;
; SAGE II
         datesage2_all=[-99.]
         if iyr lt 2006 then begin
            restore,dirs2+'cat_sage2_v6.2.'+syr
            restore,dirs2+'Theta/o3_sage2_v6.2_theta.'+syr
            datesage2_all=date
            ysage2_all=latitude
            xsage2_all=longitude
            modes2_all=sctype
            o3sage2_all=mix
         endif
;
; SAGE III
         datesage3_all=[-99.]
         if iyr ge 2002 and iyr lt 2006 then begin
            restore,dirs3+'cat_sage3_v3.00.'+syr
            restore,dirs3+'Theta/o3mlr_sage3_v3.00_theta.'+syr
            datesage3_all=date
            ysage3_all=latitude
            xsage3_all=longitude
            modes3_all=sctype
            o3sage3_all=mix
         endif
;
; HALOE
         datehal_all=[-99.]
         if iyr lt 2006 then begin
            restore,dirh+'cat_haloe_v19.'+syr
            restore,dirh+'Theta/o3_haloe_v19_theta.'+syr
            datehal_all=date
            yhal_all=latitude
            xhal_all=longitude
            modeh_all=sctype
            o3hal_all=mix
         endif
;
; POAM
         datepoam_all=[-99.]
         if iyr ge 1993 and iyr le 1996 then begin
            restore,dirp+'cat_poam2_v6.0.'+syr
            restore,dirp+'Theta/o3_poam2_v6.0_theta.'+syr
            datepoam_all=date
            ypoam_all=latitude
            xpoam_all=longitude
            modep_all=sctype
            o3poam_all=mix
         endif
         if iyr ge 1998 and iyr lt 2006 then begin
            restore,dirp+'cat_poam3_v4.0.'+syr
            restore,dirp+'Theta/o3_poam3_v4.0_theta.'+syr
            datepoam_all=date
            ypoam_all=latitude
            xpoam_all=longitude
            modep_all=sctype
            o3poam_all=mix
         endif
;
; ILAS
         dateilas_all=[-99.]
         if iyr ge 1996 and iyr le 1997 then begin
            restore,diri+'cat_ilas_v06.10.'+syr
            restore,diri+'Theta/o3_ilas_v06.10_theta.'+syr
            dateilas_all=date
            yilas_all=latitude
            xilas_all=longitude
            modei_all=sctype
            o3ilas_all=mix
         endif
         if iyr eq 2003 then begin
            restore,diri+'cat_ilas2_v1.4.'+syr
            restore,diri+'Theta/o3_ilas2_v1.4_theta.'+syr
            dateilas_all=date
            yilas_all=latitude
            xilas_all=longitude
            modei_all=sctype
            o3ilas_all=mix
         endif
;
; ACE
         dateace_all=[-99.]
         if iyr ge 2004 then begin
            restore,dira+'cat_ace_v2.2.'+syr
            restore,dira+'Theta/o3_ace_v2.2_theta.'+syr
            dateace_all=date
            yace_all=latitude
            xace_all=longitude
            modea_all=sctype
            o3ace_all=mix
         endif

      endif
;
; extract daily SOSST data
;
      norbits3=0L & norbits2=0L & norbitp=0L & norbith=0L & norbiti=0L & norbita=0L

      sage2day=where(datesage2_all eq iyr*10000L+idate,norbits2)
      if norbits2 le 1L then goto,jumpsage2
      o3sage2=reform(o3sage2_all(sage2day,*))
      thsage2=0.*o3sage2
      nlevs2=n_elements(theta)
      for k=0L,nlevs2-1L do thsage2(*,k)=theta(k)
      ysage2=reform(ysage2_all(sage2day))
      xsage2=reform(xsage2_all(sage2day))
      modes2=reform(modes2_all(sage2day))
jumpsage2:
      sage3day=where(datesage3_all eq iyr*10000L+idate,norbits3)
      if norbits3 le 1L then goto,jumpsage3
      o3sage3=reform(o3sage3_all(sage3day,*))
      thsage3=0.*o3sage3
      nlevs3=n_elements(theta)
      for k=0L,nlevs3-1L do thsage3(*,k)=theta(k)
      ysage3=reform(ysage3_all(sage3day))
      xsage3=reform(xsage3_all(sage3day))
      modes3=reform(modes3_all(sage3day))
jumpsage3:
      halday=where(datehal_all eq iyr*10000L+idate,norbith)
      if norbith le 1L then goto,jumphal
      o3hal=reform(o3hal_all(halday,*))
      thhal=0.*o3hal
      nlevh=n_elements(theta)
      for k=0L,nlevh-1L do thhal(*,k)=theta(k)
      yhal=reform(yhal_all(halday))
      xhal=reform(xhal_all(halday))
      modeh=reform(modeh_all(halday))
jumphal:
      poamday=where(datepoam_all eq iyr*10000L+idate,norbitp)
      if norbitp le 1L then goto,jumppoam
      o3poam=reform(o3poam_all(poamday,*))
      thpoam=0.*o3poam
      nlevp=n_elements(theta)
      for k=0L,nlevp-1L do thpoam(*,k)=theta(k)
      ypoam=reform(ypoam_all(poamday))
      xpoam=reform(xpoam_all(poamday))
      modep=reform(modep_all(poamday))
jumppoam:
      ilasday=where(dateilas_all eq iyr*10000L+idate,norbiti)
      if norbiti le 1L then goto,jumpilas
      o3ilas=reform(o3ilas_all(ilasday,*))
      thilas=0.*o3ilas
      nlevi=n_elements(theta)
      for k=0L,nlevi-1L do thilas(*,k)=theta(k)
      yilas=reform(yilas_all(ilasday))
      xilas=reform(xilas_all(ilasday))
      modei=reform(modei_all(ilasday))
jumpilas:
      aceday=where(dateace_all eq iyr*10000L+idate,norbita)
      if norbita le 10L then begin
         norbita=0L
         goto,jumpace
      endif
      o3ace=reform(o3ace_all(aceday,*))
      thace=0.*o3ace
      nleva=n_elements(theta)
      for k=0L,nleva-1L do thace(*,k)=theta(k)
      yace=reform(yace_all(aceday))
      xace=reform(xace_all(aceday))
      modea=reform(modea_all(aceday))
jumpace:

print,iyr,imn,idy,norbiti
;
; remove missing or bad data and separate sunrise and sunset 
;
      yhal0=-999. & yhal1=-999.
      if norbith gt 1L then begin
         index=where(yhal ge -90. and yhal le 0.,norbith)
         if index(0) ne -1 then begin
            modeh=modeh(index)
            yhal=yhal(index)
            xhal=xhal(index)
            o3hal=o3hal(index,*)*1.e6
            thhal=thhal(index,*)
         endif

         index=where(modeh eq 'r')
         if index(0) ne -1 then begin
            yhalsr=yhal(index)
            xhalsr=xhal(index)
            o3halsr=o3hal(index,*)
            thhalsr=thhal(index,*)
            index=where(yhalsr ge -90. and yhalsr le 0.)
            if index(0) eq -1 then goto,jumphalsr
            yhal0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yhalsr(index))
               yhal0=result(0)
            endif
         endif
         jumphalsr:
         index=where(modeh eq 's')
         if index(0) ne -1 then begin
            yhalss=yhal(index)
            xhalss=xhal(index)
            o3halss=o3hal(index,*)
            thhalss=thhal(index,*)
            index=where(yhalss ge -90. and yhalss le 0.)
            if index(0) eq -1 then goto,jumphalss
            yhal1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yhalss(index))
               yhal1=result(0)
            endif
         endif
         jumphalss:
      endif
      ysage30=-999. & ysage31=-999.
      if norbits3 gt 1L then begin
         index=where(ysage3 ge -90. and ysage3 le 0.,norbits3)
         if index(0) ne -1 then begin
            modes3=modes3(index)
            ysage3=ysage3(index)
            xsage3=xsage3(index)
            o3sage3=o3sage3(index,*)*1.e6
            thsage3=thsage3(index,*)
         endif

         index=where(modes3 eq 'r')
         if index(0) ne -1 then begin
            ysage3sr=ysage3(index)
            xsage3sr=xsage3(index)
            o3sage3sr=o3sage3(index,*)
            thsage3sr=thsage3(index,*)
            index=where(ysage3sr ge -90. and ysage3sr le 0.)
            if index(0) eq -1 then goto,jumpsage3sr
            ysage30=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage3sr(index))
               ysage30=result(0)
            endif
         endif
         jumpsage3sr:
         index=where(modes3 eq 's')
         if index(0) ne -1 then begin
            ysage3ss=ysage3(index)
            xsage3ss=xsage3(index)
            o3sage3ss=o3sage3(index,*)
            thsage3ss=thsage3(index,*)
            index=where(ysage3ss ge -90. and ysage3ss le 0.)
            if index(0) eq -1 then goto,jumpsage3ss
            ysage31=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage3ss(index))
               ysage31=result(0)
            endif
         endif
         jumpsage3ss:
      endif
      ysage20=-999. & ysage21=-999.
      if norbits2 gt 1L then begin
         index=where(ysage2 ge -90. and ysage2 le 0.,norbits2)
         if index(0) ne -1 then begin
            modes2=modes2(index)
            ysage2=ysage2(index)
            xsage2=xsage2(index)
            o3sage2=o3sage2(index,*)*1.e6
            thsage2=thsage2(index,*)
         endif

         index=where(modes2 eq 'r')
         if index(0) ne -1 then begin
            ysage2sr=ysage2(index)
            xsage2sr=xsage2(index)
            o3sage2sr=o3sage2(index,*)
            thsage2sr=thsage2(index,*)
            index=where(ysage2sr ge -90. and ysage2sr le 0.)
            if index(0) eq -1 then goto,jumpsage2sr
            ysage20=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2sr(index))
               ysage20=result(0)
            endif
         endif
         jumpsage2sr:
         index=where(modes2 eq 's')
         if index(0) ne -1 then begin
            ysage2ss=ysage2(index)
            xsage2ss=xsage2(index)
            o3sage2ss=o3sage2(index,*)
            thsage2ss=thsage2(index,*)
            index=where(ysage2ss ge -90. and ysage2ss le 0.)
            if index(0) eq -1 then goto,jumpsage2ss
            ysage21=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2ss(index))
               ysage21=result(0)
            endif
         endif
         jumpsage2ss:
      endif
      ypoam0=-999. & ypoam1=-999.
      if norbitp gt 1L then begin
         index=where(ypoam ge -90. and ypoam le 0.,norbitp)
         if index(0) ne -1 then begin
            modep=modep(index)
            ypoam=ypoam(index)
            xpoam=xpoam(index)
            o3poam=o3poam(index,*)*1.e6
            thpoam=thpoam(index,*)
         endif

         index=where(modep eq 'r')
         if index(0) ne -1 then begin
            ypoamsr=ypoam(index)
            xpoamsr=xpoam(index)
            o3poamsr=o3poam(index,*)
            thpoamsr=thpoam(index,*)
            index=where(ypoamsr ge -90. and ypoamsr le 0.)
            if index(0) eq -1 then goto,jumppoamsr
            ypoam0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoamsr(index))
               ypoam0=result(0)
            endif
         endif
         jumppoamsr:
         index=where(modep eq 's')
         if index(0) ne -1 then begin
            ypoamss=ypoam(index)
            xpoamss=xpoam(index)
            o3poamss=o3poam(index,*)
            thpoamss=thpoam(index,*)
            index=where(ypoamss ge -90. and ypoamss le 0.)
            if index(0) eq -1 then goto,jumppoamss
            ypoam1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoamss(index))
               ypoam1=result(0)
            endif
         endif
         jumppoamss:
      endif
      yilas0=-999. & yilas1=-999.
      if norbiti gt 1L then begin
         index=where(yilas ge -90. and yilas le 0.,norbiti)
         if index(0) ne -1 then begin
            modei=modei(index)
            yilas=yilas(index)
            xilas=xilas(index)
            o3ilas=o3ilas(index,*)*1.e6
            thilas=thilas(index,*)
         endif

         index=where(modei eq 'r')
         if index(0) ne -1 then begin
            yilassr=yilas(index)
            xilassr=xilas(index)
            o3ilassr=o3ilas(index,*)
            thilassr=thilas(index,*)
            index=where(yilassr ge -90. and yilassr le 0.)
            if index(0) eq -1 then goto,jumpilassr
            yilas0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yilassr(index))
               yilas0=result(0)
            endif
         endif
         jumpilassr:
         index=where(modei eq 's')
         if index(0) ne -1 then begin
            yilasss=yilas(index)
            xilasss=xilas(index)
            o3ilasss=o3ilas(index,*)
            thilasss=thilas(index,*)
            index=where(yilasss ge -90. and yilasss le 0.)
            if index(0) eq -1 then goto,jumpilasss
            yilas1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yilasss(index))
               yilas1=result(0)
            endif
         endif
         jumpilasss:
      endif

      yace0=-999. & yace1=-999.
      if norbita gt 1L then begin
         index=where(yace ge -90. and yace le 0.,norbita)
         if index(0) ne -1 then begin
            modea=modea(index)
            yace=yace(index)
            xace=xace(index)
            o3ace=o3ace(index,*)*1.e6
            thace=thace(index,*)
         endif

         index=where(modea eq 'r')
         if index(0) ne -1 then begin
            yacesr=yace(index)
            xacesr=xace(index)
            o3acesr=o3ace(index,*)
            thacesr=thace(index,*)
            index=where(yacesr ge -90. and yacesr le 0.)
            if index(0) eq -1 then goto,jumpacesr
            yace0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yacesr(index))
               yace0=result(0)
            endif
         endif
         jumpacesr:
         index=where(modea eq 's')
         if index(0) ne -1 then begin
            yacess=yace(index)
            xacess=xace(index)
            o3acess=o3ace(index,*)
            thacess=thace(index,*)
            index=where(yacess ge -90. and yacess le 0.)
            if index(0) eq -1 then goto,jumpacess
            yace1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yacess(index))
               yace1=result(0)
            endif
         endif
         jumpacess:
      endif

      yorder=[yhal0,yhal1,ypoam0,ypoam1,ysage30,ysage31,$
              ysage20,ysage21,yilas0,yilas1,yace0,yace1]
      instorder=['HALOESR','HALOESS','POAMSR','POAMSS','SAGE3SR','SAGE3SS',$
                 'SAGE2SR','SAGE2SS','ILASSR','ILASSS','ACESR','ACESS']
;
; remove missing and NH modes
;
      index=where(yorder ge -90. and yorder lt -30.,npan)
      if index(0) eq -1 then goto,jump
      if index(0) ne -1 then begin
         yorder=yorder(index)
         instorder=instorder(index)
      endif
      index=sort(yorder)
      yorder=yorder(index)
      instorder=instorder(index)
      npan0=4
      if npan gt npan0 then begin
         yadd=npan-npan0
         yorder=yorder(yadd:npan-1)
         instorder=instorder(yadd:npan-1)
      endif
      if npan lt npan0 then begin
         yadd=npan0-npan
         yorder=[fltarr(yadd),yorder]
         instorder=[strarr(yadd),instorder]
      endif
;
; postscript file
;
      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,$
                 filename='XZ_Ozone/xz_o3_sosst+mls+polar_sh_'+lfile+'_'+stheta+'K.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; viewport location based on number of panels
;
      erase
      xmn=0.55 & ymn=0.1
      ywide=0.175
      xwide=0.35
      for ipan=0,npan0-1 do begin
          xmx=xmn+xwide
          ymx=ymn+ywide
          set_viewport,xmn,xmx,ymn,ymx
;
; plot longitude-altitude sections in correct order S-N
;
          case 1 of
            (instorder(ipan) eq 'HALOESR'): begin
             ydata=yhalsr
             xdata=xhalsr
             o3data=o3halsr
             thdata=thhalsr
            end
            (instorder(ipan) eq 'HALOESS'): begin
             ydata=yhalss
             xdata=xhalss
             o3data=o3halss
             thdata=thhalss
            end
            (instorder(ipan) eq 'POAMSR'): begin
             ydata=ypoamsr
             xdata=xpoamsr
             o3data=o3poamsr
             thdata=thpoamsr
            end
            (instorder(ipan) eq 'POAMSS'): begin
             ydata=ypoamss
             xdata=xpoamss
             o3data=o3poamss
             thdata=thpoamss
            end
            (instorder(ipan) eq 'SAGE3SR'): begin
             ydata=ysage3sr
             xdata=xsage3sr
             o3data=o3sage3sr
             thdata=thsage3sr
            end
            (instorder(ipan) eq 'SAGE3SS'): begin
             ydata=ysage3ss
             xdata=xsage3ss
             o3data=o3sage3ss
             thdata=thsage3ss
            end
            (instorder(ipan) eq 'SAGE2SR'): begin
             ydata=ysage2sr
             xdata=xsage2sr
             o3data=o3sage2sr
             thdata=thsage2sr
            end
            (instorder(ipan) eq 'SAGE2SS'): begin
             ydata=ysage2ss
             xdata=xsage2ss
             o3data=o3sage2ss
             thdata=thsage2ss
            end
            (instorder(ipan) eq 'ILASSR'): begin
             ydata=yilassr
             xdata=xilassr
             o3data=o3ilassr
             thdata=thilassr
            end
            (instorder(ipan) eq 'ILASSS'): begin
             ydata=yilasss
             xdata=xilasss
             o3data=o3ilasss
             thdata=thilasss
            end
            (instorder(ipan) eq 'ACESR'): begin
             ydata=yacesr
             xdata=xacesr
             o3data=o3acesr
             thdata=thacesr
            end
            (instorder(ipan) eq 'ACESS'): begin
             ydata=yacess
             xdata=xacess
             o3data=o3acess
             thdata=thacess
            end
            else: begin
            goto,noprof
            end
          endcase
          ydata=reform(ydata)
          xdata=reform(xdata)
          o3data=reform(o3data)
          thdata=reform(thdata)
; 
; remove bad xdata and sort in longitude
;
          index=where(xdata ge 0. and xdata le 360.)
          if index(0) ne -1 then begin
             ydata=reform(ydata(index))
             xdata=reform(xdata(index))
             o3data=reform(o3data(index,*))
             thdata=reform(thdata(index,*))
          endif
          xsave=xdata
          xdata=0.*thdata
          result=size(thdata)
          ndim=result(0)
          if ndim eq 1 then goto,oneprof
          nl=result(2)
          for i=0,nl-1 do begin
              sindex=sort(xsave)
              xdata(*,i)=xsave(sindex)
              o3data(*,i)=o3data(sindex,i)
              thdata(*,i)=thdata(sindex,i)
          endfor
          xlabels='Lat='+string(ydata(sindex),format='(f5.1)')
          level=2.0+0.5*findgen(17)
          nlvls=n_elements(level)
          col1=1+indgen(nlvls)*mcolor/nlvls
          !type=2^2+2^3
          contour,o3data,xdata,thdata,levels=level,/cell_fill,$
                  title=instorder(ipan)+'      '+xlabels(0),c_color=col1,$
                  min_value=-99.,xticks=4,xrange=[0.,360.],yrange=[500.,1600.],$
                  color=0
          contour,o3data,xdata,thdata,levels=findgen(nlvls),/follow,$
                  /overplot,color=0,min_value=-99.
;
; closest longitude-altitude slice of UKMO data
;
          result=moment(ydata)
          y0=result(0)
          oneprof:
          if ndim eq 1 then y0=ydata(0)
          index=where(abs(alat-y0) le 7.5,nlat)
          pv1=reform(pv2(index(0),*,*))
          p1=reform(p2(index(0),*,*))
          mark1=reform(mark2(index(0),*,*))
          for j=1,nlat-1 do begin
          llat=index(j)
          pv1=pv1+reform(pv2(llat,*,*))
          p1=p1+reform(p2(llat,*,*))
          mark1=mark1+reform(mark2(llat,*,*))
          endfor
          pv1=pv1/float(nlat)
          p1=p1/float(nlat)
          mark1=mark1/float(nlat)
          t1=0.*mark1
          for k=0,nth-1 do t1(*,k)=th(k)*((p1(*,k)/1000.)^(.286))
          contour,smooth(mark1,3),alon,th,levels=[0.1],/follow,$
                  /overplot,color=mcolor,thick=3,c_labels=[0]
          index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=0.
          contour,smooth(mark1,3),alon,th,levels=[-0.1],/follow,$
                  /overplot,color=0,thick=6,c_labels=[0]
          index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=mcolor

          plots,0.,rtheta,/data
          plots,360.,rtheta,/data,/continue,color=0
          imin=min(level)
          imax=max(level)
          xmnb=xmx+.03
          xmxb=xmnb+.01
          set_viewport,xmnb,xmxb,ymn,ymx
          !type=2^2+2^3+2^5
          plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],color=0
          xbox=[0,10,10,0,0]
          y1=imin
          dy=(imax-imin)/float(nlvls)
          for j=0,nlvls-1 do begin
              ybox=[y1,y1,y1+dy,y1+dy,y1]
              polyfill,xbox,ybox,color=col1(j)
              y1=y1+dy
          endfor
          noprof:
          ymn=ymx+0.05
      endfor
      pv1=transpose(pv2(*,*,itheta))
      p1=transpose(p2(*,*,itheta))
      sf1=transpose(sf2(*,*,itheta))
      mark1=transpose(mark2(*,*,itheta))
      pv=fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1
      pv(nc,*)=pv(0,*)
      p=fltarr(nc+1,nr)
      p(0:nc-1,0:nr-1)=p1
      p(nc,*)=p(0,*)
      sf=fltarr(nc+1,nr)
      sf(0:nc-1,0:nr-1)=sf1
      sf(nc,*)=sf(0,*)
      mark=fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
      !type=2^2+2^3
      xyouts,.15,.95,syr+smn+sdy+'  '+stheta+' K',/normal,color=0,charsize=2
      set_viewport,.05,.45,.5,.9
      MAP_SET,-90,0,0,/stereo,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5
      oplot,findgen(361),-0.1+0.*findgen(361),psym=0,color=0
      contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark gt 0. and y2d lt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark lt 0. and y2d lt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      contour,mark,x,alat,levels=[.1],c_color=mcolor*.1,/overplot,/follow,c_labels=0,/noeras,thick=3
      contour,mark,x,alat,levels=[-.1],c_color=0,/overplot,/follow,c_labels=0,/noeras,thick=3
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
    index=where(mtheta gt rtheta-80. and mtheta le rtheta+80. and mlat2 ge -90. and mlat2 lt 0.,npt)
    th0=string(FORMAT='(I4)',itheta)
    xdata=mlon2(index) & ydata=mlat2(index)
    o3data=mo3(index)
    o3mask=mo3mask(index)
    tdata=mtemp(index)
    index=where(xdata lt 0.)
    if index(0) ne -1 then xdata(index)=xdata(index)+360.
;   omin=-2.0
;   omax=max(o3data)
omin=min(level)
omax=max(level)
    for i=0L,n_elements(o3data)-1L do begin
        if o3mask(i) ne -99. then $
           oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                 color=((o3data(i)-omin)/(omax-omin))*mcolor,symsize=1.5
        if o3data(i) gt omax then $
           oplot,[xdata(i),xdata(i)],[ydata(i),ydata(i)],psym=8,$
                 color=0.95*mcolor,symsize=1.5
    endfor
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
index=where(mlat lt -81.5,npt)
if index(0) eq -1L then goto,jump
tswath=muttime(index)
flag=0.*tswath
nswath=50L
tsave=fltarr(nswath)
icount=0L
    for i=0L,npt-1L do begin
        index=where(abs(tswath(i)-tswath) lt 1. and flag eq 0.)
        if index(0) ne -1L then begin
           flag(index)=1.0
           kindex=where(abs(muttime-tswath(index(0))) le 0.5 and mlat lt 0.,nn)
;          oplot,mlon(kindex),mlat(kindex),psym=8,symsize=0.5,color=(float(i+1)/float(npt))*mcolor
           stime=string(FORMAT='(F4.1)',tswath(index(0)))
           tsave(icount)=tswath(index(0))
           icount=icount+1L
           xtmp=mlon(kindex)
           ytmp=mlat(kindex)
;          index=where(ytmp eq max(ytmp))
;          xyouts,xtmp(index(0)),5.,stime,/data,charsize=1.25,alignment=0.5,color=0
;
; put time labels on both endpoints of swath
;
x0=xtmp(0) & y0=ytmp(0)
x1=xtmp(nn-1) & y1=ytmp(nn-1)
           xyouts,x0,5.,stime,/data,charsize=1.25,alignment=0.5,color=0
           xyouts,x1,5.,stime,/data,charsize=1.25,alignment=0.5,color=0
        endif
    endfor
tsave=tsave(0:icount-1L)
tplot=19.9791
print,tsave
read,'Central Swath times ',tplot
kindex=where(abs(muttime-tplot) le 0.5 and mlat lt 0.,mprof)
stime=string(FORMAT='(F4.1)',tplot)
for kk=0,mprof-1L,4L do $
    oplot,[mlon(kindex(kk)),mlon(kindex(kk))],[mlat(kindex(kk)),mlat(kindex(kk))],$
           psym=2,symsize=1.75,color=lc
;
; superimpose SOSST locations colored by ozone mixing ratio
;
      o3min=min(level)
      o3max=max(level)
      if norbith gt 1L then begin
         xhal2d=0.*thhal
         yhal2d=0.*thhal
         for k=0L,long(nlevh)-1L do begin
             xhal2d(*,k)=xhal
             yhal2d(*,k)=yhal
         endfor
         index=where(abs(thhal-rtheta) le 50. and o3hal gt 0. and yhal2d lt 0.,hcount)
         if hcount gt 0L then begin
            xday=xhal2d(index)
            yday=yhal2d(index)
            thday=thhal(index)
            o3day=o3hal(index)
            a=findgen(9)*(2*!pi/8.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,hcount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(9)*(2*!pi/8.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
      endif
      if norbits2 gt 1L then begin
         xsage2d=0.*thsage2
         ysage2d=0.*thsage2
         for k=0L,long(nlevs2)-1L do begin
             xsage2d(*,k)=xsage2
             ysage2d(*,k)=ysage2
         endfor
         index=where(abs(thsage2-rtheta) le 50. and o3sage2 gt 0. and ysage2d lt 0.,scount)
         if scount gt 0L then begin
            xday=xsage2d(index)
            yday=ysage2d(index)
            thday=thsage2(index)
            o3day=o3sage2(index)
            a=findgen(5)*(2*!pi/4.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,scount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(5)*(2*!pi/4.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
      endif
      if norbits3 gt 1L then begin
         xsage3d=0.*thsage3
         ysage3d=0.*thsage3
         for k=0L,long(nlevs3)-1L do begin
             xsage3d(*,k)=xsage3
             ysage3d(*,k)=ysage3
         endfor
         index=where(abs(thsage3-rtheta) le 50. and o3sage3 gt 0. and ysage3d lt 0.,s3count)
         if s3count gt 0L then begin
            xday=xsage3d(index)
            yday=ysage3d(index)
            thday=thsage3(index)
            o3day=o3sage3(index)
            a=findgen(5)*(2*!pi/4.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,s3count-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(5)*(2*!pi/4.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
      endif
      if norbitp gt 1L then begin
         xpoam2d=0.*thpoam
         ypoam2d=0.*thpoam
         for k=0L,long(nlevp)-1L do begin
             xpoam2d(*,k)=xpoam
             ypoam2d(*,k)=ypoam
         endfor
         index=where(abs(thpoam-rtheta) le 50. and o3poam gt 0. and ypoam2d lt 0.,pcount)
         if pcount gt 0L then begin
            xday=xpoam2d(index)
            yday=ypoam2d(index)
            thday=thpoam(index)
            o3day=o3poam(index)
            a=findgen(4)*(2*!pi/3.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,pcount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],symsize=1.5,$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(4)*(2*!pi/3.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc,symsize=1.5
         endif
      endif
      if norbiti gt 1L then begin
         xilas2d=0.*thilas
         yilas2d=0.*thilas
         for k=0L,long(nlevi)-1L do begin
             xilas2d(*,k)=xilas
             yilas2d(*,k)=yilas
         endfor
         index=where(abs(thilas-rtheta) le 50. and o3ilas gt 0. and yilas2d lt 0.,icount)
         if icount gt 0L then begin
            xday=xilas2d(index)
            yday=yilas2d(index)
            thday=thilas(index)
            o3day=o3ilas(index)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,icount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
      endif
      if norbita gt 1L then begin
         xace2d=0.*thace
         yace2d=0.*thace
         for k=0L,long(nleva)-1L do begin
             xace2d(*,k)=xace
             yace2d(*,k)=yace
         endfor
         index=where(abs(thace-rtheta) le 50. and o3ace gt 0. and yace2d lt 0.,hcount)
         if hcount gt 0L then begin
            xday=xace2d(index)
            yday=yace2d(index)
            thday=thace(index)
            o3day=o3ace(index)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,hcount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(6)*(2*!pi/5.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc
         endif
      endif
;
; plot MLS swath
;
      !type=2^2+2^3
      set_viewport,.075,.475,.1,.4
o3swath=transpose(smooth(mo3(*,kindex),3,/edge_truncate))
o3mask=transpose(mo3mask(*,kindex))
index=where(o3swath lt 0.)
if index(0) ne -1L then o3swath(index)=0.1
o3pr=0.*o3swath
index=where(abs(o3swath) gt 100. or o3mask eq -99.)
if index(0) ne -1L then o3swath(index)=0.
result=size(o3swath)
nz=result(2)
for k=0L,nz-1L do begin
    result=moment(o3swath(*,k))
    o3pr(*,k)=o3swath(*,k)-result(0)
endfor
thswath=transpose(mtheta(*,kindex))
lonswath=transpose(mlon2(*,kindex))
latswath=transpose(mlat2(*,kindex))
markswath=0.*latswath
xswath=0.*thswath
ylabels=string(format='(f5.1)',mlat(kindex))
xlabels=string(format='(f5.1)',mlon(kindex))
for i=0L,mlev-1L do xswath(*,i)=findgen(mprof)
;
; interpolate MetO marker to MLS swath
;
for ii=0L,mprof-1L do begin
    for kk=0L,mlev-1L do begin
        slon=lonswath(ii,kk)
        slat=latswath(ii,kk)
        slev=thswath(ii,kk)

        if slon lt alon(0) then slon=slon+360.
        for i=0L,nc-1L do begin
            ip1=i+1
            if i eq nc-1L then ip1=0L
            xlon=alon(i)
            xlonp1=alon(ip1)
            if i eq nc-1L then xlonp1=360.+alon(ip1)
            if slon ge xlon and slon le xlonp1 then begin
               xscale=(slon-xlon)/(xlonp1-xlon)
               goto,jumpx
            endif
        endfor
jumpx:
        for j=0L,nr-2L do begin
            jp1=j+1
            xlat=alat(j)
            xlatp1=alat(jp1)
            if slat ge xlat and slat le xlatp1 then begin
                yscale=(slat-xlat)/(xlatp1-xlat)
                goto,jumpy
            endif
        endfor
jumpy:
        for k=1L,nth-1L do begin
            kp1=k-1             ; UKMO data is "top down"
            uz=th(k)
            uzp1=th(kp1)
            if slev ge uz and slev le uzp1 then begin
               zscale=(slev-uz)/(uzp1-uz)
               pj1=mark2(j,i,k)+xscale*(mark2(j,ip1,k)-mark2(j,i,k))
               pjp1=mark2(jp1,i,k)+xscale*(mark2(jp1,ip1,k)-mark2(jp1,i,k))
               pj2=mark2(j,i,kp1)+xscale*(mark2(j,ip1,kp1)-mark2(j,i,kp1))
               pjp2=mark2(jp1,i,kp1)+xscale*(mark2(jp1,ip1,kp1)-mark2(jp1,i,kp1))
               p1=pj1+yscale*(pjp1-pj1)
               p2=pj2+yscale*(pjp2-pj2)
               markswath(ii,kk)=p1+zscale*(p2-p1)
               goto,jumpz
            endif
        endfor
jumpz:
    endfor
endfor
;
loadct,38
    !type=2^2+2^3
;o3swath=o3pr
;    omin=min(o3swath)
;    omax=max(o3swath)
;omin=0.
;omax=10.
;    nlvls=21
;    level=omin+((omax-omin)/nlvls)*findgen(nlvls)
;    col1=1+indgen(nlvls)*mcolor/nlvls
          level=2.0+0.5*findgen(17)
omin=min(level)
omax=max(level)
          nlvls=n_elements(level)
          col1=1+indgen(nlvls)*mcolor/nlvls

    contour,o3swath,xswath,thswath,levels=level,/cell_fill,title='MLS Swath at '+stime+' UT',$
            c_color=col1,min_value=-999.,xticks=1,xrange=[0.,mprof-1L],xtickname=[' ',' '],$
            yrange=[500.,1600.],charsize=1.25,color=0
    index=where(level gt 0.)
    contour,o3swath,xswath,thswath,levels=level(index),/follow,color=0,min_value=-999.,/overplot,c_labels=0*level
;   index=where(level lt 0.)
;   contour,o3swath,xswath,thswath,levels=level(index),/follow,color=0,min_value=-999.,/overplot,c_labels=0*level
    contour,o3swath,xswath,thswath,levels=[0.0],/follow,color=0,thick=3,min_value=-999.,/overplot,c_labels=0*level
    contour,markswath,xswath,thswath,levels=[-0.1],/follow,color=0,thick=5,min_value=-999.,/overplot,c_labels=0*level
    contour,markswath,xswath,thswath,levels=[0.1],/follow,color=mcolor,thick=5,min_value=-999.,/overplot,c_labels=0*level
for i=0L,mprof-1L,20 do begin
    plots,i,500
    plots,i,450,/data,/continue,color=0
    xyouts,i,400.,xlabels(i),/data,charsize=1.25,alignment=0.5,color=0
    xyouts,i,300.,ylabels(i),/data,charsize=1.25,alignment=0.5,color=0
endfor
plots,mprof-1,500
plots,mprof-1,450,/data,/continue,color=0
xyouts,mprof-1,400.,xlabels(mprof-1),/data,charsize=1.25,alignment=0.5,color=0
xyouts,mprof-1,300.,ylabels(mprof-1),/data,charsize=1.25,alignment=0.5,color=0

;   imin=omin
;   imax=omax
;   xmnb=0.45+.03
;   xmxb=xmnb+.01
;   set_viewport,xmnb,xmxb,0.15,0.45
;   !type=2^2+2^3+2^5
;   plot,[0,0],[imin,imax],xrange=[0,10],yrange=[imin,imax],title='O3 ppmv',color=0
;   xbox=[0,10,10,0,0]
;   y1=imin
;   dy=(imax-imin)/float(nlvls)
;   for j=0,nlvls-1 do begin
;       ybox=[y1,y1,y1+dy,y1+dy,y1]
;       polyfill,xbox,ybox,color=col1(j)
;       y1=y1+dy
;   endfor

      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert -trim XZ_Ozone/xz_o3_sosst+mls+polar_sh_'+lfile+'_'+stheta+'K.ps -rotate -90 '+$
               ' XZ_Ozone/xz_o3_sosst+mls+polar_sh_'+lfile+'_'+stheta+'K.jpg'
         spawn,'/usr/bin/rm XZ_Ozone/xz_o3_sosst+mls+polar_sh_'+lfile+'_'+stheta+'K.ps'
      endif
      icount=icount+1L
      goto,jump
end
