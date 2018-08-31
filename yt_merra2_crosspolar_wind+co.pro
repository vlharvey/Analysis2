;
; MERRA2 version
; cross polar Hovmoller of the Arctic vortex and MLS CO
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra2_nc3

loadct,39
mcolor=byte(!p.color)
icolmax=byte(!p.color)
icolmax=fix(icolmax)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,2*cos(a),2*sin(a),/fill
nxdim=800
nydim=800
xorig=[0.15]
yorig=[0.25]
xlen=0.8
ylen=0.6
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
lstmn=12
lstdy=1
ledmn=4
leddy=1
lstday=0
ledday=0
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
;
; loop over years
;
;for lstyr=1979,2013 do begin
for lstyr=2010,2010 do begin
   ledyr=lstyr+1
;
; Ask interactive questions- get starting/ending date
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
!noeras=1
dir='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
dir2='/atmos/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_U_V_v4.2_'
dirm='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit

      sdate=strcompress(string(FORMAT='(I4,I2.2,I2.2)',iyr,imn,idy))
;
; check for sufficient number of daily profiles
;
      dum=findfile(dirm+'cat_mls_v4.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmerra
      restore,dirm+'cat_mls_v4.2_'+sdate+'.sav'
      print,sdate,' MLS profs ',n_elements(id)
      if n_elements(id) lt 2000L then goto,skipmerra
;
; BRO             FLOAT     = Array[144, 96, 37, 2]
; CLO             FLOAT     = Array[144, 96, 37, 2]
; CO              FLOAT     = Array[144, 96, 37, 2]
; GPH             FLOAT     = Array[144, 96, 55, 2]
; H2O             FLOAT     = Array[144, 96, 55, 2]
; HCL             FLOAT     = Array[144, 96, 37, 2]
; HNO3            FLOAT     = Array[144, 96, 37, 2]
; HO2             FLOAT     = Array[144, 96, 49, 2]
; LAT             DOUBLE    = Array[96]
; LON             DOUBLE    = Array[144]
; N2O             FLOAT     = Array[144, 96, 37, 2]
; NODE            STRING    = Array[2]
; O3              FLOAT     = Array[144, 96, 55, 2]
; OH              FLOAT     = Array[144, 96, 49, 2]
; PMLS            FLOAT     = Array[37]
; PMLS2           FLOAT     = Array[55]
; PMLS3           FLOAT     = Array[49]
; T               FLOAT     = Array[144, 96, 55, 2]
; U               FLOAT     = Array[144, 96, 55, 2]
; V               FLOAT     = Array[144, 96, 55, 2]
;
      dum=findfile(dir2+sdate+'.sav')
      if dum(0) eq '' then goto,skipmerra
      restore,dir2+sdate+'.sav'
      mu2=mean(u,dim=4,/Nan)
      mv2=mean(v,dim=4,/Nan)
      mz2=mean(gph,dim=4,/Nan)
      mco2=mean(co,dim=4,/Nan)
      msp2=sqrt(mu2^2+mv2^2)

      ncfile0=dir+sdate+'00.nc3'
      rd_merra2_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
         u2,v2,qdf2,mark2,qv2,z2,sf2,q2,o32,iflag
      if iflag ne 0L then goto,skipmerra
      tmp2=0.*p2
      for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
      index=where(mark2 lt 0.)
      if index(0) ne -1L then mark2(index)=-1.*(mark2(index)/mark2(index))
      sp2=sqrt(u2^2+v2^2)
;
; ask for level and lon on first day
;
      if kcount eq 0L then begin
         merravortex_yt=fltarr(kday,nr)
         merratemp_yt=fltarr(kday,nr)
         merrasp_yt=fltarr(kday,nr)
         mlssp_yt=fltarr(kday,nr)
         mlsco_yt=fltarr(kday,nr)
         sdate_all=strarr(kday)
;        print,th
         rtheta=4000.
;        read,'Enter desired theta surface ',rtheta
         index=where(th eq rtheta)
         itheta=index(0)
         stheta=strcompress(long(rtheta),/remove_all)
;
; find closest pressure level for 55 element and 37 element arrays
;
t1=mean(t,dim=4,/Nan)
tzm=mean(t1,dim=1,/Nan)
tprof=mean(tzm,dim=1,/Nan)
theta_mls55=tprof*(1000./PMLS2)^0.286
index=where(abs(theta_mls55-rtheta) eq min(abs(theta_mls55-rtheta),/Nan))
ilev=index(0)
print,'does MERRA2 theta=',rtheta,' jibe with this pressure level? ',PMLS2(ilev)
index2=where(abs(pmls-PMLS2(ilev)) eq min(abs(pmls-PMLS2(ilev)),/Nan))
print,'closest CO pressure level is ',PMLS(index2)
ilev2=index2(0)

;        print,alon
         rlon=180.
;        read,'Enter desired longitude ',rlon
         index=where(rlon eq alon)
         ilon=index(0)
         if ilon lt nc/2 then ilon2=ilon+nc/2
         if ilon ge nc/2 then ilon2=ilon-nc/2
         slon=strcompress(long(rlon),/remove_all)
         slon2=strcompress(long(alon(ilon2)),/remove_all)
         kcount=1
goto,quick
      endif
      sdate_all(icount)=sdate
;
; extract cross polar swath on this day
;
      for j=nr/2L,nr-1L do begin
          merravortex_yt(icount,j-nr/2)=mark2(j,ilon,itheta)
          merratemp_yt(icount,j-nr/2)=tmp2(j,ilon,itheta)
          merrasp_yt(icount,j-nr/2)=sp2(j,ilon,itheta)
          mlssp_yt(icount,j-nr/2)=msp2(ilon,j,ilev)
          mlsco_yt(icount,j-nr/2)=mco2(ilon,j,ilev2)
      endfor
      jj=nr/2
      for j=nr-1,nr/2,-1L do begin
          merravortex_yt(icount,jj)=mark2(j,ilon2,itheta)
          merratemp_yt(icount,jj)=tmp2(j,ilon2,itheta)
          merrasp_yt(icount,jj)=sp2(j,ilon2,itheta)
          mlssp_yt(icount,jj)=msp2(ilon2,j,ilev)
          mlsco_yt(icount,jj)=mco2(ilon2,j,ilev2)
          jj=jj+1
      endfor
skipmerra:
      icount=icount+1L

goto, jump

plotit:
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)      ;+'/'+sday(xindex)
index=where(syear ne 0)
minyear=long(min(long(syear(index))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)

save,file='yt_merra2_crosspolar_vortex_'+yearlab+'_'+stheta+'_'+slon+'.sav',merravortex_yt,merratemp_yt,merrasp_yt,mlssp_yt,mlsco_yt,yearlab,stheta,slon,sdate_all,kday,alat,slon2
quick:
restore,'yt_merra2_crosspolar_vortex_2010-2011_4000_180.sav
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)      ;+'/'+sday(xindex)
index=where(syear ne 0)
minyear=long(min(long(syear(index))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)


if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='yt_merra_crosspolar_vortex_'+yearlab+'_'+stheta+'_'+slon+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif

;
; plot latitude-time temperature and the vortex/anticyclones
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(merravortex_yt eq 0.)
if index(0) ne -1L then merravortex_yt(index)=0./0.
if index(0) ne -1L then merratemp_yt(index)=0./0.
merravortex_yt=smooth(merravortex_yt,3,/NaN)
level=.5*findgen(21)
nlvls=n_elements(level)
tlevel=180.+5*findgen(nlvls)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,merravortex_yt,1.+findgen(kday),alat,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
      charsize=1.5,color=0,ytitle='Latitude',/fill,c_color=0.*col1,title=yearlab+' MERRA2 Arctic Vortex '+stheta+' K',$
      levels=[0.1],yticks=6,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.,$
      ytickname=['Eq','30N','60N','NP','60N','30N','Eq']
;colevel=[0.01,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.]
;nlvls=n_elements(colevel)
;col1=1+indgen(nlvls)*icolmax/nlvls
contour,smooth(mlsco_yt,5,/Nan,/edge_truncate)*1.e6,1.+findgen(kday),alat,/fill,levels=level,c_color=col1,/overplot
contour,smooth(mlsco_yt,5,/Nan,/edge_truncate)*1.e6,1.+findgen(kday),alat,/foll,levels=level,color=0,/overplot,thick=3,c_labels=0*level

loadct,0
;contour,merratemp_yt,1.+findgen(kday),alat,/cell_fill,levels=tlevel,c_color=col1,/overplot
;contour,merratemp_yt,1.+findgen(kday),alat,levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
index=where(finite(merravortex_yt) ne 1) 
merravortex_yt(index)=0.
contour,merravortex_yt,1.+findgen(kday),alat,levels=[0.1,0.5,0.9],color=mcolor,/follow,/overplot,thick=10,c_labels=[0]
merrasp_yt=smooth(merrasp_yt,3,/edge_truncate)
loadct,10
contour,merrasp_yt,1.+findgen(kday),alat,levels=30+20*findgen(6),c_color=[100,170,180,200,230],/follow,/overplot,thick=5
xyouts,xmn-0.1,ymn,slon,/normal,color=0,charsize=2
xyouts,xmn-0.1,ymx,slon2,/normal,color=0,charsize=2
loadct,39

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS CO (ppmv)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yt_merra_crosspolar_vortex_'+yearlab+'_'+stheta+'_'+slon+'.ps -rotate -90 yt_merra_crosspolar_vortex_'+yearlab+'_'+stheta+'_'+slon+'.png'
endif
endfor
end
