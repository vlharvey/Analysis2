;
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
nlvls=20
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
; remove missing modes
;
      index=where(yorder ne -999.,npan)
      if index(0) eq -1 then goto,jump
      if index(0) ne -1 then begin
         yorder=yorder(index)
         instorder=instorder(index)
      endif
      index=sort(yorder)
      yorder=yorder(index)
      instorder=instorder(index)
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
                 filename='xz_o3_soundings+v_'+ifile+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; viewport location based on number of panels
;
      erase
      xmn=0.1 & ymn=0.1
      ywide=(0.9-0.05*float(npan))/float(npan)
      xwide=ywide*2.
      if npan lt 4 then xwide=ywide*1.5
      for ipan=0,npan-2 do begin
          xmx=xmn+xwide
          ymx=ymn+ywide
          xmn=(xmx+xmn)/2.
          ymn=ymx+0.05
      endfor
      plots,0.1,0.1,/normal
      plots,xmn,ymn,/normal,/continue
      plots,0.1+xwide,0.1,/normal
      plots,xmn+xwide,ymn,/normal,/continue
      xmn=0.1 & ymn=0.1
      ywide=(0.9-0.05*float(npan))/float(npan)
      xwide=ywide*2.
      if npan lt 4 then xwide=ywide*1.5
      for ipan=0,npan-1 do begin
          xmx=xmn+xwide
          ymx=ymn+ywide
          set_viewport,xmn,xmx,ymn,ymx
          xmn=(xmx+xmn)/2.
          ymn=ymx+0.05
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
            (instorder(ipan) eq 'POAM3SR'): begin
             ydata=ypoam3sr
             xdata=xpoam3sr
             o3data=o3poam3sr
             thdata=thpoam3sr
            end
            (instorder(ipan) eq 'POAM3SS'): begin
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
          xlabels=strcompress(string(ydata(sindex)),/remove_all)
          level=0.5*findgen(nlvls)
          !type=2^2+2^3
          contour,o3data*1.e6,xdata,thdata,levels=level,/cell_fill,$
                  title=instorder(ipan)+'      '+xlabels(0),c_color=col1,$
                  min_value=-999.,xticks=4,xrange=[0.,360.],yrange=[300.,2000.]
          contour,o3data*1.e6,xdata,thdata,levels=findgen(nlvls),/follow,$
                  /overplot,color=0
;
; closest longitude-altitude slice of UKMO data
;
          result=moment(ydata)
          y0=result(0)
          oneprof:
          if ndim eq 1 then y0=ydata(0)
          index=where(abs(alat-y0) le 5.,nlat)
          pv1=reform(pv2(index(0),*,*))
          p1=reform(p2(index(0),*,*))
          v1=reform(v2(index(0),*,*))
          mark1=reform(mark2(index(0),*,*))
          for j=1,nlat-1 do begin
          llat=index(j)
          pv1=pv1+reform(pv2(llat,*,*))
          p1=p1+reform(p2(llat,*,*))
          v1=v1+reform(v2(llat,*,*))
          mark1=mark1+reform(mark2(llat,*,*))
          endfor
          pv1=pv1/float(nlat)
          p1=p1/float(nlat)
          v1=v1/float(nlat)
          mark1=mark1/float(nlat)
          t1=0.*mark1
          for k=0,nth-1 do t1(*,k)=th(k)*((p1(*,k)/1000.)^(.286))
          if ydata(0) gt 0. then begin
          level=[5.,15.,30.,50.,75.,100.]
          contour,v1,alon,th,levels=level,/follow,$
                  /overplot,color=mcolor,c_labels=[0],thick=2
          index=where(v1 gt 5.)
          if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=mcolor
          level=[-100.,-75.,-50.,-30.,-15.,-5.]
          contour,v1,alon,th,levels=level,/follow,$
                  /overplot,color=0,c_labels=[0],thick=2
          index=where(v1 lt -5.)
          if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=0
          endif
          if ydata(0) lt 0. then begin
          level=[5.,15.,30.,50.,75.,100.]
          contour,v1,alon,th,levels=level,/follow,$
                  /overplot,color=0,c_labels=[0],thick=2
          index=where(v1 gt 5.)
          if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=0
          level=[-100.,-75.,-50.,-30.,-15.,-5.]
          contour,v1,alon,th,levels=level,/follow,$
                  /overplot,color=mcolor,c_labels=[0],thick=2
          index=where(v1 lt -5.)
          if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=mcolor
          endif

;         contour,mark1,alon,th,levels=[0.1],/follow,$
;                 /overplot,color=0,thick=3,c_labels=[0]
;         index=where(mark1 gt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=0.
;         contour,mark1,alon,th,levels=[-0.1],/follow,$
;                 /overplot,color=mcolor,thick=3,c_labels=[0]
;         index=where(mark1 lt 0.)
;         if index(0) ne -1 then oplot,x2d(index),y2d(index),psym=3,color=mcolor
      endfor
      xyouts,.1,.8,sdy+' '+month(imn-1)+' '+syr,charsize=3,/normal
      if setplot eq 'x' then begin
         save=assoc(3,bytarr(nxdim,nydim))
         img=bytarr(nxdim,nydim)
         img(0,0)=TVRD(0,0,nxdim,nydim)
         write_gif,'xz_o3_soundings+v_'+ifile+'.gif',img
      endif
      if setplot eq 'ps' then device, /close
      icount=icount+1L
;     stop
      goto,jump
end
