;
; SH
; plot satellite data in longitude altitude section
; add polar projection at a user prompted theta
; VLH 1/27/04
;
@rd_sage3_o3_soundings
@rd_haloe_o3_soundings
@rd_poam3_o3_soundings
@rd_sage2_o3_soundings
@rd_ilas_o3_soundings
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
setplot='ps'
read,'setplot=',setplot
nxdim=750 & nydim=750
xorig=[0.10]
yorig=[0.10]
cbaryoff=0.015
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
nmon=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
dirs='/aura3/data/SAGE_II_data/Sound_data/sage2_'
dirs3='/aura3/data/SAGE_III_data/Sound_data/sage3_solar_'
diri='/aura3/data/ILAS_data/Sound_data/ilas_'
dirh='/aura3/data/HALOE_data/Sound_data/haloe_'
dirp3='/aura3/data/POAM_data/Sound_data/poam3_'
dirp2='/aura3/data/POAM_data/Sound_data/poam2_'
diru='/aura3/data/UKMO_data/Datfiles/ukmo_'
ifile='                             '
lstmn=6 & lstdy=1 & lstyr=4 & lstday=0
ledmn=12 & leddy=1 & ledyr=5 & ledday=0
;
; Ask interactive questions- get starting/ending date
;
;print, ' '
;print, '      UKMO Version '
;print, ' '
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
      if ndays gt ledday then stop,' Normal termination condition '
;
; read UKMO data
;
      syr=strtrim(string(iyr),2)
      sdy=string(FORMAT='(i2.2)',idy)
      uyr=strmid(syr,2,2)
      ifile=mon(imn-1)+sdy+'_'+uyr
      lfile=nmon(imn-1)+'_'+sdy+'_'+uyr
      rd_ukmo_nc3,diru+ifile+'.nc3',nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,vp2,sf2,iflag
      if iflag eq 1 then goto,jump
      if icount eq 0L then begin
         rtheta=1000.
         print,th
         read,' Enter theta level ',rtheta
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
; read satellite ozone soundings
;
      sfile=mon(imn-1)+sdy+'_'+syr
      norbits3=0L & norbits2=0L & norbitp3=0L & norbith=0L & norbiti=0L
;     rd_sage3_o3_soundings,dirs3+sfile+'_o3.sound',norbits3,tsage3,$
;        xsage3,ysage3,tropps3,tropzs3,tropths3,modes3,o3sage3,psage3,$
;        thsage3,zsage3,clsage3,qo3sage3,nlevs3
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
         thhal,zhal,qo3hal,nlevh
      rd_ilas_o3_soundings,diri+sfile+'_o3.sound',norbiti,tilas,$
         xilas,yilas,troppi,tropzi,tropthi,modei,o3ilas,pilas,$
         thilas,zilas,clilas,qo3ilas,nlevi
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
             if max(thhal(*,ilev)) gt 0. then good(ilev)=1.
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
      yilas0=-999. & yilas1=-999.
      if norbiti gt 0L then begin
         index=where(yilas ge -90. and yilas le 90.,norbith)
         if index(0) ne -1 then begin
            yilas=yilas(index)
            xilas=xilas(index)
            o3ilas=o3ilas(index,*)
            thilas=thilas(index,*)
         endif
         good=0.*fltarr(nlevi)
         for ilev=0,nlevi-1 do $
             if max(thilas(*,ilev)) lt 1.00000e+24 then good(ilev)=1.
         index=where(good eq 1.,nlevi)
         thilas=thilas(*,index)
         o3ilas=o3ilas(*,index)
         qo3ilas=qo3ilas(*,index)
         index=where(qo3ilas/o3ilas gt 0.5)
         if index(0) ne -1 then o3ilas(index)=-999.
         index=where(modei eq 0L)
         if index(0) ne -1 then begin
            yilassr=yilas(index)
            xilassr=xilas(index)
            o3ilassr=o3ilas(index,*)
            thilassr=thilas(index,*)
            index=where(yilassr ge -90. and yilassr le 90.)
            if index(0) eq -1 then goto,jumpilassr
            yilas0=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yilassr(index))
               yilas0=result(0)
            endif
         endif
         jumpilassr:
         index=where(modei eq 1L)
         if index(0) ne -1 then begin
            yilasss=yilas(index)
            xilasss=xilas(index)
            o3ilasss=o3ilas(index,*)
            thilasss=thilas(index,*)
            index=where(yilasss ge -90. and yilasss le 90.)
            if index(0) eq -1 then goto,jumpilasss
            yilas1=index(0)
            if n_elements(index) gt 1 then begin
               result=moment(yilasss(index))
               yilas1=result(0)
            endif
         endif
         jumpilasss:
      endif
      yorder=[yhal0,yhal1,ypoam30,ypoam31,ysage30,ysage31,$
              ysage20,ysage21,yilas0,yilas1]
      instorder=['HALOESR','HALOESS','POAMSR','POAMSS','SAGE3SR','SAGE3SS',$
                 'SAGE2SR','SAGE2SS','ILASSR','ILASSS']
;
; remove missing and NH modes
;
      index=where(yorder ge -90. and yorder lt -10.,npan)
      if index(0) eq -1 then goto,jump
      if index(0) ne -1 then begin
         yorder=yorder(index)
         instorder=instorder(index)
      endif
      index=sort(yorder)
      yorder=yorder(index)
      instorder=instorder(index)
      npan0=3
      if npan gt npan0 then begin
         yadd=npan-npan0
         yorder=yorder(0:npan0-1)
         instorder=instorder(0:npan0-1)
      endif
      if npan lt 3 then begin
         yadd=npan0-npan
         yorder=[yorder,fltarr(yadd)]
         instorder=[instorder,strarr(yadd)]
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
                 filename='xz_o3_soundings+polar_sh_'+lfile+'_'+stheta+'K.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; viewport location based on number of panels
;
      erase
      xmn=0.55 & ymn=0.2
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
             o3data=o3halsr*1.e-6
             thdata=thhalsr
            end
            (instorder(ipan) eq 'HALOESS'): begin
             ydata=yhalss
             xdata=xhalss
             o3data=o3halss*1.e-6
             thdata=thhalss
            end
            (instorder(ipan) eq 'POAMSR'): begin
             ydata=ypoam3sr
             xdata=xpoam3sr
             o3data=o3poam3sr
             thdata=thpoam3sr
            end
            (instorder(ipan) eq 'POAMSS'): begin
             ydata=ypoam3ss
             xdata=xpoam3ss
             o3data=o3poam3ss
             thdata=thpoam3ss
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
          level=0.5+0.5*findgen(20)
          nlvls=n_elements(level)
          col1=1+indgen(nlvls)*mcolor/nlvls
          !type=2^2+2^3
          contour,o3data*1.e6,xdata,thdata,levels=level,/cell_fill,$
                  title=instorder(ipan)+'      '+xlabels(0),c_color=col1,$
                  min_value=-999.,xticks=4,xrange=[0.,360.],yrange=[300.,2000.],$
                  ytitle='Theta',xtitle='Longitude'
          contour,o3data*1.e6,xdata,thdata,levels=findgen(nlvls),/follow,$
                  /overplot,color=0
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
                  /overplot,color=0,thick=3,c_labels=[0]
          index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=0.
          contour,smooth(mark1,3),alon,th,levels=[-0.1],/follow,$
                  /overplot,color=mcolor,thick=6,c_labels=[0]
          index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,xz2d(index),yz2d(index),psym=3,color=mcolor

          plots,0.,rtheta,/data
          plots,360.,rtheta,/data,/continue,color=mcolor
          imin=0.
          imax=max(level)
          xmnb=xmx+.03
          xmxb=xmnb+.01
          set_viewport,xmnb,xmxb,ymn,ymx
          !type=2^2+2^3+2^5
          plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax]
          xbox=[0,10,10,0,0]
          y1=imin
          dy=(imax-imin)/float(nlvls)
          for j=0,nlvls-1 do begin
              ybox=[y1,y1,y1+dy,y1+dy,y1]
              polyfill,xbox,ybox,color=col1(j)
              y1=y1+dy
          endfor
          noprof:
          ymn=ymx+0.075
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
      set_viewport,.05,.45,.45,.85
      MAP_SET,-90,0,-180,/stereo,/contin,/grid,/noeras,color=lc,/noborder,charsize=1.5,$
              title=lfile+'   '+stheta+' K'
      oplot,findgen(361),-0.1+0.*findgen(361),psym=0
      contour,sf,x,alat,nlevels=30,c_color=lc,/overplot,/follow,c_labels=0,/noeras
      index=where(mark gt 0. and y2d lt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=2,color=lc
      index=where(mark lt 0. and y2d lt 0.)
      if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=4,color=lc
      contour,mark,x,alat,levels=[.1],c_color=mcolor*.1,/overplot,/follow,c_labels=0,/noeras,thick=3
      contour,mark,x,alat,levels=[-.1],c_color=lc,/overplot,/follow,c_labels=0,/noeras,thick=3
      o3min=2.0
      o3max=max(level)
      if norbith gt 0L then begin
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
      if norbits2 gt 0L then begin
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
            o3day=o3sage2(index)*1.e6
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
      if norbits3 gt 0L then begin
         xsage3d=0.*thsage3
         ysage3d=0.*thsage3
         for k=0L,long(nlevs3)-1L do begin
             xsage3d(*,k)=xsage3
             ysage3d(*,k)=ysage3
         endfor
         index=where(abs(thsage3-rtheta) le 50. and o3sage3 gt 0. and ysage3d lt 0.,scount)
         if scount gt 0L then begin
            xday=xsage3d(index)
            yday=ysage3d(index)
            thday=thsage3(index)
            o3day=o3sage3(index)*1.e6
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
      if norbitp3 gt 0L then begin
         xpoam2d=0.*thpoam3
         ypoam2d=0.*thpoam3
         for k=0L,long(nlevp3)-1L do begin
             xpoam2d(*,k)=xpoam3
             ypoam2d(*,k)=ypoam3
         endfor
         index=where(abs(thpoam3-rtheta) le 50. and o3poam3 gt 0. and ypoam2d lt 0.,pcount)
         if pcount gt 0L then begin
            xday=xpoam2d(index)
            yday=ypoam2d(index)
            thday=thpoam3(index)
            o3day=o3poam3(index)*1.e6
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
      if norbiti gt 0L then begin
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
            o3day=o3ilas(index)*1.e6
            a=findgen(4)*(2*!pi/3.)
            usersym,2.*cos(a),2.*sin(a),/fill
            for i=0,icount-1 do $
                oplot,[xday(i),xday(i)],[yday(i),yday(i)],symsize=1.5,$
                      psym=8,color=mcolor*(o3day(i)-o3min)/(o3max-o3min)
            a=findgen(4)*(2*!pi/3.)
            usersym,2.*cos(a),2.*sin(a)
            oplot,xday,yday,psym=8,color=lc,symsize=1.5
         endif
      endif

      imin=o3min
      imax=max(level)
      ymnb=.45-cbaryoff
      ymxb=ymnb+cbarydel
      set_viewport,.05,.45,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],$
            xtitle='Ozone (ppmv)',charsize=1.5
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=col1(j)
          x1=x1+dx
      endfor

      if setplot ne 'ps' then stop
      if setplot eq 'ps' then begin
         device, /close
         spawn,'convert xz_o3_soundings+polar_sh_'+lfile+'_'+stheta+'K.ps -rotate -90 '+$
               ' XZ_Ozone/xz_o3_soundings+polar_sh_'+lfile+'_'+stheta+'K.jpg'
      endif
      icount=icount+1L
      goto,jump
end
