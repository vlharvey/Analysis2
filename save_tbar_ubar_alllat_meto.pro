;
; MetO zonal mean temperature and zonal wind at all latitudes for 18 year record
;
; VLH 9/10/09
;
@stddat
@kgmt
@ckday
@kdate
@rd_ukmo

lstmn=9L & lstdy=28L & lstyr=1991L
lstmn=3L & lstdy=12L & lstyr=2009L
ledmn=8L & leddy=1L & ledyr=2009L
;
; Get start and end dates
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 1991 or lstyr gt 2009 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 2009 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
nfile=kday
yyyymmdd_all=lonarr(nfile)
syyyymmdd_all=strarr(nfile)
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L

; --- Loop here over days --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' Starting day outside range '
      if ndays gt ledday then goto,saveit
if iyr ge 2000L then iyr1=iyr-2000L
if iyr lt 2000L then iyr1=iyr-1900L

      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      date=syr+smn+sdy
print,date
      yyyymmdd_all(icount)=long(syr+smn+sdy)
      syyyymmdd_all(icount)=date
;
; read MetO pressure data, i.e. /aura7/harvey/UKMO_data/Datfiles/ppassm_yYY_mMM_dDD_h12.pp.dat
;
;***Read UKMO data
      file='/aura7/harvey/UKMO_data/Datfiles/ppassm_y'+$
            string(FORMAT='(i2.2,a2,i2.2,a2,i2.2,a11)',$
            iyr1,'_m',imn,'_d',idy,'_h12.pp.dat')
      rd_ukmo,file,iflg,nc,nr,nl,alon,alat,wlon,wlat,press,gp,tp,uu,vv
      if iflg eq 1 or max(gp) eq -10000.0 or max(gp) eq 0. then goto,skip
;
; declare time period arrays
;
      if kcount eq 0L then begin
         tbar_all=-9999.+0.*fltarr(nfile,nr,nl)
         ubar_all=-9999.+0.*fltarr(nfile,nr,nl)
         gbar_all=-9999.+0.*fltarr(nfile,nr,nl)
         kcount=1
      endif
;
; interpolate MetO winds to temperature latitudes
;
; t,gph go from 90   N to -90   S by 2.5  (73 lats)
;           and 0      to 360     by 3.75 (96 lons)
; winds go from 88.75N to -88.75S by 2.5  (72 lats)
;           and 1.875  to 358.125 by 3.75 (96 lons)
;
      up=0.*tp
      vp=0.*tp
      for l=0,nl-1 do begin
          for i=0,nc-2L do begin
            im1=i-1
            ip1=i
            if i eq 0 then ip1=0
            if i eq 0 then im1=nc-1
            for j=1,nr-2 do begin
                jm1=j-1
                jp1=j
                uim1=.5*(uu(im1,j-1,l)+uu(im1,j,l))
                uip1=.5*(uu(ip1,j-1,l)+uu(ip1,j,l))
                up(i,j-1,l)=.5*(uim1+uip1)

                vim1=.5*(vv(im1,j-1,l)+vv(im1,j,l))
                vip1=.5*(vv(ip1,j-1,l)+vv(ip1,j,l))
                vp(i,j-1,l)=.5*(vim1+vip1)
;print,wlat(jm1),alat(j),wlat(jp1)
;print,wlon(im1),alon(i),wlon(ip1)
;print,' '
            endfor
;
; assign poles equal to 88.75 values
;
            up(*,0,l)=up(*,1,l)
            vp(*,0,l)=vp(*,1,l)
            up(*,nr-1,l)=up(*,nr-2,l)
            vp(*,nr-1,l)=vp(*,nr-2,l)
          endfor
      endfor
;
; reverse latitudes to go from SP to NP
;
      alat=reverse(alat)
      upr=0.*alat
      vpr=0.*alat
      tpr=0.*alat
      gpr=0.*alat
      for l=0,nl-1 do begin
          for i=0,nc-1L do begin
              upr=reverse(reform(up(i,*,l)))
              vpr=reverse(reform(vp(i,*,l)))
              tpr=reverse(reform(tp(i,*,l)))
              gpr=reverse(reform(gp(i,*,l)))
              up(i,*,l)=upr
              vp(i,*,l)=vpr
              tp(i,*,l)=tpr
              gp(i,*,l)=gpr
          endfor
      endfor
;
; calculate zonal mean temperature and zonal wind
;
      uzm=-9999.+0.*fltarr(nr,nl)
      tzm=-9999.+0.*fltarr(nr,nl)
      gzm=-9999.+0.*fltarr(nr,nl)
      for k=0,nl-1 do begin
          for j=0,nr-1 do begin
              tzm(j,k)=total(tp(*,j,k))/float(nc)
              uzm(j,k)=total(up(*,j,k))/float(nc)
              gzm(j,k)=total(gp(*,j,k))/float(nc)
              if tzm(j,k) gt 400. or tzm(j,k) lt 100. then stop,'check temperature'
              if uzm(j,k) lt -200. or uzm(j,k) gt 200. then stop,'check zonal wind'
          endfor
      endfor
;
; check
;
erase
contour,uzm,alat,press,/ylog,yrange=[1000.,min(press)],/noerase,levels=5.*findgen(100),title=date
contour,uzm,alat,press,/ylog,/noerase,levels=-100.+5.*findgen(20),c_linestyle=5,/overplot
;contour,tzm,alat,press,/ylog,yrange=[1000.,1],/noerase
;stop
;
; retain all daily zonal means 
;
      tbar_all(icount,*,*)=tzm
      ubar_all(icount,*,*)=uzm
      gbar_all(icount,*,*)=gzm
;
; skip if file is missing but still increment day
;
      skip:
      icount=icount+1L
goto,jump
saveit:

syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
sday=strmid(syyyymmdd_all,6,2)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
nday=icount
save,file='MetO_tbar_ubar_alllat_'+yearlab+'.sav',nl,nr,alat,nday,yyyymmdd_all,press,tbar_all,ubar_all,gbar_all
end
