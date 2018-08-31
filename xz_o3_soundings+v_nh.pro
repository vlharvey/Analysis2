;
; NH
; plot satellite data in longitude altitude section
; VLH 10/5/2003
;
@rd_sage3_o3_soundings
@rd_haloe_o3_soundings
@rd_poam3_o3_soundings
@rd_sage2_o3_soundings
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
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
icmm1=icolmax-1
icmm2=icolmax-2
setplot='x'
;read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.10]
yorig=[0.10]
cbaryoff=0.02
cbarydel=0.02
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
dirs='/aura3/data/SAGE_II_data/Sound_data/sage2_'
dirs3='/aura3/data/SAGE_III_data/Sound_data/sage3_solar_'
diri='/aura3/data/ILAS_data/Sound_data/ilas_'
dirh='/aura3/data/HALOE_data/Sound_data/haloe_'
dirp3='/aura3/data/POAM_data/Sound_data/poam3_'
dirp2='/aura3/data/POAM_data/Sound_data/poam2_'
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
lstmn=2 & lstdy=10 & lstyr=3 & lstday=0
ledmn=2 & leddy=10 & ledyr=3 & ledday=0
;
; Ask interactive questions- get starting/ending date
;
print, ' '
print, '      UKMO Version '
print, ' '
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
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
      uyr=strmid(syr,2,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
      if iflag eq 1 then goto,jump
      x2d=fltarr(nc,nth)
      y2d=fltarr(nc,nth)
      for k=0,nth-1 do x2d(*,k)=alon
      for j=0,nc-1 do y2d(j,*)=th
;
; read satellite ozone soundings
;
      sfile=mon(imn-1)+sdy+'_'+syr
      rd_sage3_o3_soundings,dirs3+sfile+'_o3.sound',norbits3,tsage3,$
         xsage3,ysage3,tropps3,tropzs3,tropths3,modes3,o3sage3,psage3,$
         thsage3,zsage3,clsage3,qo3sage3,nlevs3
      rd_sage2_o3_soundings,dirs+sfile+'_o3.sound',norbits2,tsage2,$
         xsage2,ysage2,tropps2,tropzs2,tropths2,modes2,o3sage2,psage2,$
         thsage2,zsage2,clsage2,qo3sage2,nlevs2
      if iyr lt 1998 then begin
      rd_poam3_o3_soundings,dirp2+sfile+'_o3.sound',norbitp3,tpoam3,$
         xpoam3,ypoam3,troppp3,tropzp3,tropthp3,modep3,o3poam3,ppoam3,$
         thpoam3,zpoam3,clpoam3,qo3poam3,nlevp3
      endif
      if iyr ge 1998 then begin
      rd_poam3_o3_soundings,dirp3+sfile+'_o3.sound',norbitp3,tpoam3,$
         xpoam3,ypoam3,troppp3,tropzp3,tropthp3,modep3,o3poam3,ppoam3,$
         thpoam3,zpoam3,clpoam3,qo3poam3,nlevp3
      endif
      rd_haloe_o3_soundings,dirh+sfile+'_o3.sound',norbith,thal,$
         xhal,yhal,tropph,tropzh,tropthh,modeh,o3hal,phal,$
         thhal,zhal,clhal,qo3hal,nlevh
;
; remove missing and bad data
;
      yhal0=-999. & yhal1=-999.
      if norbith gt 0L then begin
         index=where(yhal ge -90. and yhal le 90.,norbith)
         if index(0) ne -1 then begin
            yhal=yhal(index)
            xhal=xhal(index)
            o3hal=o3hal(index,*)
            thhal=thhal(index,*)
         endif
         good=0.*fltarr(nlevh)
         for ilev=0,nlevh-1 do $
             if max(thhal(*,ilev)) lt 1.00000e+24 then good(ilev)=1.
         index=where(good eq 1.,nlevh)
         thhal=thhal(*,index)
         o3hal=o3hal(*,index)
         qo3hal=qo3hal(*,index)
         index=where(qo3hal/o3hal gt 0.5)
         if index(0) ne -1 then o3hal(index)=-999.
         index=where(modeh eq 0L)
         if index(0) ne -1 then begin
            yhalsr=yhal(index)
            xhalsr=xhal(index)
            o3halsr=o3hal(index,*)
            thhalsr=thhal(index,*)
            index=where(yhalsr ge -90. and yhalsr le 90.)
            if index(0) eq -1 then goto,jumphalsr
            yhal0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yhalsr(index))
               yhal0=result(0)
            endif
         endif
         jumphalsr:
         index=where(modeh eq 1L)
         if index(0) ne -1 then begin
            yhalss=yhal(index)
            xhalss=xhal(index)
            o3halss=o3hal(index,*)
            thhalss=thhal(index,*)
            index=where(yhalss ge -90. and yhalss le 90.)
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
      if norbits3 gt 0L then begin
         index=where(qo3sage3/o3sage3 gt 0.5)
         if index(0) ne -1 then o3sage3(index)=-999.
         index=where(modes3 eq 0L)
         if index(0) ne -1 then begin
            ysage3sr=ysage3(index)
            xsage3sr=xsage3(index)
            o3sage3sr=o3sage3(index,*)
            thsage3sr=thsage3(index,*)
            index=where(ysage3sr ge -90. and ysage3sr le 90.)
            if index(0) eq -1 then goto,jumpsage3sr
            ysage30=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage3sr(index))
               ysage30=result(0)
            endif
         endif
         jumpsage3sr:
         index=where(modes3 eq 1L)
         if index(0) ne -1 then begin
            ysage3ss=ysage3(index)
            xsage3ss=xsage3(index)
            o3sage3ss=o3sage3(index,*)
            thsage3ss=thsage3(index,*)
            index=where(ysage3ss ge -90. and ysage3ss le 90.)
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
      if norbits2 gt 0L then begin
         index=where(qo3sage2/o3sage2 gt 0.5)
         if index(0) ne -1 then o3sage2(index)=-999.
         index=where(modes2 eq 0L)
         if index(0) ne -1 then begin
            ysage2sr=ysage2(index)
            xsage2sr=xsage2(index)
            o3sage2sr=o3sage2(index,*)
            thsage2sr=thsage2(index,*)
            index=where(ysage2sr ge -90. and ysage2sr le 90.)
            if index(0) eq -1 then goto,jumpsage2sr
            ysage20=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2sr(index))
               ysage20=result(0)
            endif
         endif
         jumpsage2sr:
         index=where(modes2 eq 1L)
         if index(0) ne -1 then begin
            ysage2ss=ysage2(index)
            xsage2ss=xsage2(index)
            o3sage2ss=o3sage2(index,*)
            thsage2ss=thsage2(index,*)
            index=where(ysage2ss ge -90. and ysage2ss le 90.)
            if index(0) eq -1 then goto,jumpsage2ss
            ysage21=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ysage2ss(index))
               ysage21=result(0)
            endif
         endif
         jumpsage2ss:
      endif
      ypoam30=-999. & ypoam31=-999.
      if norbitp3 gt 0L then begin
         index=where(qo3poam3/o3poam3 gt 0.5)
         if index(0) ne -1 then o3poam3(index)=-999.
         index=where(modep3 eq 0L)
         if index(0) ne -1 then begin
            ypoam3sr=ypoam3(index)
            xpoam3sr=xpoam3(index)
            o3poam3sr=o3poam3(index,*)
            thpoam3sr=thpoam3(index,*)
            index=where(ypoam3sr ge -90. and ypoam3sr le 90.)
            if index(0) eq -1 then goto,jumppoam3sr
            ypoam30=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoam3sr(index))
               ypoam30=result(0)
            endif
         endif
         jumppoam3sr:
         index=where(modep3 eq 1L)
         if index(0) ne -1 then begin
            ypoam3ss=ypoam3(index)
            xpoam3ss=xpoam3(index)
            o3poam3ss=o3poam3(index,*)
            thpoam3ss=thpoam3(index,*)
            index=where(ypoam3ss ge -90. and ypoam3ss le 90.)
            if index(0) eq -1 then goto,jumppoam3ss
            ypoam31=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(ypoam3ss(index))
               ypoam31=result(0)
            endif
         endif
         jumppoam3ss:
      endif
      yorder=[yhal0,yhal1,ypoam30,ypoam31,ysage30,ysage31,ysage20,ysage21]
      instorder=['HALOESR','HALOESS','POAM3SR','POAM3SS',$
                 'SAGE3SR','SAGE3SS','SAGE2SR','SAGE2SS']
;
; remove missing and SH modes
;
      index=where(yorder gt 0.,npan)
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
      if npan lt 4 then begin
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
                 filename='xz_o3_soundings+v_nh_'+ifile+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
      erase
;
; polar orthographic of PV, vortex/anticyclone boundaries, and occultations
;
      !type=2^2+2^3
      set_viewport,0.05,0.35,0.55,0.85
      date=sdy+' '+month(imn-1)+' '+syr
      xyouts,.08,.93,date,/normal,charsize=3
      theta=1200.
      index=where(theta eq th)
      thlev=index(0)
      stheta=strcompress(string(fix(theta)),/remove_all)
      x=fltarr(nc+1)
      x(0:nc-1)=alon
      x(nc)=alon(0)+360.
      lon=fltarr(nc+1,nr)
      lat=fltarr(nc+1,nr)
      for i=0,nc   do lat(i,*)=alat
      for j=0,nr-1 do lon(*,j)=x
      pv1=transpose(pv2(*,*,thlev))
      mark1=transpose(mark2(*,*,thlev))
      pv=0.*fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
      pv(nc,*)=pv(0,*)
      mark=0.*fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
      MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,/noborder,$
              title=stheta+' K',charsize=2.0
      index=where(lat gt 0.)
      pvmin=0.
      pvmax=max(pv(index))+0.05*max(pv(index))
      pvint=(pvmax-pvmin)/nlvls
      pvlevel=pvmin+pvint*findgen(nlvls)
      col1=1+indgen(nlvls)*icolmax/float(nlvls)
      contour,pv,x,alat,/overplot,levels=pvlevel,c_color=col1,$
             /cell_fill,/noeras
      contour,pv,x,alat,/overplot,levels=pvlevel,/follow,$
              c_labels=0*pvlevel,/noeras,color=0
      contour,mark,x,alat,/overplot,levels=[0.1],thick=5,color=0
      contour,mark,x,alat,/overplot,levels=[-0.1],thick=5,color=mcolor
      MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,/noborder,$
              charsize=2.0,latdel=10,color=0
      for j=0,npan0-1 do begin
          oplot,findgen(361),yorder(j)+0.*findgen(361),psym=0,color=mcolor
      endfor
;
; occultation points colored by ozone
;
          omin=0.
          omax=10.
          if norbits3 gt 0 then begin
             norbit=norbits3
             for i=0,norbit-1 do begin
                 o3prof=reform(o3sage3(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jump1s3
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thsage3(i,index))
                 xx=xsage3(i)
                 yy=ysage3(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(8)*(2*!pi/8.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(9)*(2*!pi/8.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jump1s3:
             endfor
          endif
          if norbitp3 gt 0 then begin
             norbit=norbitp3
             for i=0,norbit-1 do begin
                 o3prof=reform(o3poam3(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jump1p
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thpoam3(i,index))
                 xx=xpoam3(i)
                 yy=ypoam3(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(4)*(2*!pi/4.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(5)*(2*!pi/4.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jump1p:
             endfor
          endif
          if norbits2 gt 0 then begin
             norbit=norbits2
             for i=0,norbit-1 do begin
                 o3prof=reform(o3sage2(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jump1s
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thsage2(i,index))
                 xx=xsage2(i)
                 yy=ysage2(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(5)*(2*!pi/5.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(6)*(2*!pi/5.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jump1s:
             endfor
          endif
          if norbith gt 0 then begin
             norbit=norbith
             for i=0,norbit-1 do begin
                 o3prof=reform(o3hal(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jump1h
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thhal(i,index))
                 xx=xhal(i)
                 yy=yhal(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(3)*(2*!pi/3.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(4)*(2*!pi/3.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jump1h:
             endfor
          endif
; horizontal PV color bar
          ymnb=0.55-cbaryoff
          ymxb=ymnb+cbarydel
          set_viewport,0.05,.35,ymnb,ymxb
          !type=2^2+2^3+2^6
          plot,[min(pvlevel),max(pvlevel)],[0,0],yrange=[0,10],$
               xrange=[min(pvlevel),max(pvlevel)],$
               xtitle='Potential Vorticity (PVU)'
          ybox=[0,10,10,0,0]
          x1=min(pvlevel)
          dx=(max(pvlevel)-min(pvlevel))/float(nlvls)
          for j=0,nlvls-1 do begin
              xbox=[x1,x1,x1+dx,x1+dx,x1]
              polyfill,xbox,ybox,color=col1(j)
              x1=x1+dx
          endfor

; vertical ozone color bar
          xmnb=0.35+0.01
          xmxb=xmnb+cbarydel
          set_viewport,xmnb,xmxb,0.55,0.85
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
          axis,10,omin,0,YAX=1,/DATA,charsize=1.5,/ynozero
          xyouts,xmxb+0.03,.75,'Ozone (ppmv)',$
                 /normal,charsize=1.5,orientation=-90.

      !type=2^2+2^3
      set_viewport,0.05,0.35,0.1,0.4
      date=sdy+' '+month(imn-1)+' '+syr
      theta=800.
      index=where(theta eq th)
      thlev=index(0)
      stheta=strcompress(string(fix(theta)),/remove_all)
      x=fltarr(nc+1)
      x(0:nc-1)=alon
      x(nc)=alon(0)+360.
      pv1=transpose(pv2(*,*,thlev))
      mark1=transpose(mark2(*,*,thlev))
      pv=0.*fltarr(nc+1,nr)
      pv(0:nc-1,0:nr-1)=pv1(0:nc-1,0:nr-1)
      pv(nc,*)=pv(0,*)
      mark=0.*fltarr(nc+1,nr)
      mark(0:nc-1,0:nr-1)=mark1(0:nc-1,0:nr-1)
      mark(nc,*)=mark(0,*)
      MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,/noborder,$
              title=stheta+' K',charsize=2.0
      index=where(lat gt 0.)
      pvmin=0.
      pvmax=max(pv(index))+0.05*max(pv(index))
      pvint=(pvmax-pvmin)/nlvls
      pvlevel=pvmin+pvint*findgen(nlvls)
      col1=1+indgen(nlvls)*icolmax/float(nlvls)
      contour,pv,x,alat,/overplot,levels=pvlevel,c_color=col1,$
             /cell_fill,/noeras
      contour,pv,x,alat,/overplot,levels=pvlevel,/follow,$
              c_labels=0*pvlevel,/noeras,color=0
      contour,mark,x,alat,/overplot,levels=[0.1],thick=5,color=0
      contour,mark,x,alat,/overplot,levels=[-0.1],thick=5,color=mcolor
      MAP_SET,90,0,-90,/ortho,/noeras,/grid,/contin,/noborder,$
              charsize=2.0,latdel=10,color=0
      for j=0,npan0-1 do begin
          oplot,findgen(361),yorder(j)+0.*findgen(361),psym=0,color=mcolor
      endfor
          omin=0.
          omax=10.
          if norbits3 gt 0 then begin
             norbit=norbits3
             for i=0,norbit-1 do begin
                 o3prof=reform(o3sage3(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jumps3
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thsage3(i,index))
                 xx=xsage3(i)
                 yy=ysage3(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(8)*(2*!pi/8.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(9)*(2*!pi/8.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jumps3:
             endfor
          endif
          if norbitp3 gt 0 then begin
             norbit=norbitp3
             for i=0,norbit-1 do begin
                 o3prof=reform(o3poam3(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jumpp
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thpoam3(i,index))
                 xx=xpoam3(i)
                 yy=ypoam3(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(4)*(2*!pi/4.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(5)*(2*!pi/4.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jumpp:
             endfor
          endif
          if norbits2 gt 0 then begin
             norbit=norbits2
             for i=0,norbit-1 do begin
                 o3prof=reform(o3sage2(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jumps
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thsage2(i,index))
                 xx=xsage2(i)
                 yy=ysage2(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(5)*(2*!pi/5.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(6)*(2*!pi/5.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jumps:
             endfor
          endif
          if norbith gt 0 then begin
             norbit=norbith
             for i=0,norbit-1 do begin
                 o3prof=reform(o3hal(i,*))
                 index=where(o3prof gt 0.)
                 if index(0) eq -1 then goto,jumph
                 o3prof=o3prof(index)*1.e6
                 thprof=reform(thhal(i,index))
                 xx=xhal(i)
                 yy=yhal(i)
                 dth=min(abs(thprof-theta))
                 kindex=where(abs(thprof-theta) eq dth)
                 a=findgen(3)*(2*!pi/3.)
                 usersym,cos(a),sin(a),/fill
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,$
                        color=(o3prof(kindex(0))/omax)*icolmax
                 a=findgen(4)*(2*!pi/3.)
                 usersym,cos(a),sin(a)
                 oplot,[xx,xx],[yy,yy],psym=8,symsize=2,color=mcolor
                 jumph:
             endfor
          endif
; horizontal PV color bar
          ymnb=0.1-cbaryoff
          ymxb=ymnb+cbarydel
          set_viewport,0.05,.35,ymnb,ymxb
          !type=2^2+2^3+2^6
          plot,[min(pvlevel),max(pvlevel)],[0,0],yrange=[0,10],$
               xrange=[min(pvlevel),max(pvlevel)],$
               xtitle='Potential Vorticity (PVU)'
          ybox=[
