;
; zt of chemical distinction value d
;
@rd_mls_nc3
@write_mls_nc3

PRO Bimodal, X, A, F, pder

        F = A[0] * exp(-((x-A[1])^2.)/(2.*(A[2]^2.)))+ A[3] * exp(-((x-A[4])^2.)/(2.*(A[5]^2.)))

        IF N_PARAMS() GE 4 THEN begin ; If the procedure is called with four parameters, calculate the partial derivatives.
                pder = FLTARR(N_ELEMENTS(X), 6)
                ; Compute the partial derivatives with respect to A
                pder[*, 0] = exp(-((x-A[1])^2)/(2*A[2]^2.))
                pder[*, 3] = exp(-((x-A[4])^2)/(2*A[5]^2.))

                pder[*, 1] = A[0]*(x-A[1])*exp(-(x-A[1])^2./(2*A[2]^2.))/A[2]^2.
                pder[*, 4] = A[3]*(x-A[4])*exp(-(x-A[4])^2./(2*A[5]^2.))/A[5]^2.

                pder[*, 2] = A[0]*(A[1]-x)^2.*exp(-(A[1]-x)^2./(2.*A[2]^2.))/(2.*A[2]^3)
                pder[*, 5] = A[3]*(A[4]-x)^2.*exp(-(A[4]-x)^2./(2.*A[5]^2.))/(2.*A[5]^3)
        endif
end

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.20]
yorig=[0.25]
xlen=0.6
ylen=0.6
device,decompose=0
mcolor=byte(!p.color)
nlvls=26L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
syear=['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
nyear=n_elements(syear)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; get file listing
;
dir='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid_theta_'

for iyear=0L,nyear-1L do begin
for imon=0L,nmon-1L do begin
ifiles=file_search(dir+syear(iyear)+smon(imon)+'??.nc3',count=nfile)
print,syear(iyear),smon(imon),nfile
if nfile eq 0L then goto,jumpmon
;
; loop over files
;
icount=0L
FOR n=0l,nfile-1l DO BEGIN
    result=strsplit(ifiles(n),'.',/extract)
    result2=result(0)
    result3=strsplit(result2,'_',/extract)
    sdate=result3(-1)
    print,sdate
    ifile=dir+sdate+'.nc3'
    rd_mls_nc3,ifile,nc,nr,nth,alon,alat,th,pv2,p2,$
       u2,v2,qdf2,mark2,co2,z2,sf2,h2o2,markco2,iflag
    if iflag ne 0L then goto,jumpstep
    tmp2=0.*p2
    for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
;
; monthly means
;
if icount eq 0L then begin
;  rth=2400.
;  rth=3000.
;  print,th
;  read,'Enter desired theta ',rth
;  index=where(abs(rth-th) eq min(abs(rth-th)))
;  ith=index(0)
;  sth=strcompress(long(th(ith)),/remove_all)+'K'
   x=fltarr(nc+1)
   x(0:nc-1)=alon(0:nc-1)
   x(nc)=alon(0)+360.
   x2d=fltarr(nc,nr)
   y2d=fltarr(nc,nr)
   for i=0,nc-1 do y2d(i,*)=alat
   for j=0,nr-1 do x2d(*,j)=alon
   zt_chemd=fltarr(nfile,nth)
endif
;
; loop over theta
;
for ith=0L,nth-1L do begin
    rth=th(ith)
;
; extract level
;
    h2o=transpose(h2o2(*,*,ith))*1.e6
    co=transpose(co2(*,*,ith))*1.e6
;
; hemispheric PDF
;
nhindex=where(y2d lt 0.,nn)
imin=min(co(nhindex))
imax=max(co(nhindex))
nlvls=26L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
level=-0.5+0.1*findgen(nlvls)
if imax gt 2.0 then level=imin+((2.*imax-imin)/float(nlvls))*findgen(nlvls)
y2=histogram(co(nhindex),min=min(level),max=max(level),binsize=level(1)-level(0))/float(nn)
result=size(y2)	; not always nlvls in dimension
nlvls=result(1)
level=level(0:nlvls-1)
y2hem=smooth(100.*y2,3)
;y2hem=100.*y2
;
; code from Jeff France
;
;                       X = level
                        Y = y2hem
                        weights = Y*0. + 1. ; Define a vector of weights.
                        Aguess = fltarr(6)
                        m1 = where(y eq max(y))
                        Aguess[0] = y[m1]*1.
                        Aguess[1] = level[m1]*1.
                        Aguess[2] = .3

                        m2 = 0L
;                       smoothcopdf = smooth(y2hem,3)
smoothcopdf=y2hem
                        for i = 1L, n_elements(Y) - 2L do begin
if i ne m1(0) then begin
                                if smoothcopdf[i] gt smoothcopdf[i-1L] and $
                                   smoothcopdf[i] gt smoothcopdf[i+1L] then m2 = i ;selects next local maximum as initial guess for 2nd gaussian
endif
                        endfor

                        Aguess[3] = y[m2]*1.
                        Aguess[4] = level[m2]*1.
                        Aguess[5] = .3
                        A = aguess

                        ; Provide an initial guess of the function's parameters.
                        yfit = CURVEFIT(level, Y, weights, A, SIGMA, FUNCTION_NAME='bimodal')

tmp1=A[0] * exp((-((level-A[1])^2)/(2*A[2]^2)))
tmp2=A[3] * exp((-((level-A[4])^2)/(2*A[5]^2)))

copdflow=tmp1
copdfhigh=tmp2

index1=where(tmp1 eq max(tmp1))
index2=where(tmp2 eq max(tmp2))
if index2(0) ne 0 and index1(0) ne 0 and index1 gt index2 then begin
copdfhigh=tmp1
copdflow=tmp2
endif
index=where(copdfhigh eq max(copdfhigh))
comean=level(index(0))
;
; integrate area under high lat curve and intersecting area
;
higharea=total(copdfhigh)
intarea=0.
for ilvl=0L,nlvls-1L do begin
    if copdflow(ilvl) ne 0. and copdfhigh(ilvl) ne 0. then begin
       intarea=intarea+min([copdflow(ilvl),copdfhigh(ilvl)])
;print,ilvl,min([copdflow(ilvl),copdfhigh(ilvl)])
    endif
endfor
d=intarea/higharea
zt_chemd(n,ith)=d
print,th(ith),d
;
; end code
;
;nhindex=where(y2d lt 60.)
;y2=histogram(co(nhindex),min=min(level),max=max(level),binsize=level(1)-level(0))/float(nn)
;y2s=smooth(100.*y2,3)
;;;y2s=100.*y2
;y2slow=max(y2hem)*y2s/max(y2s)
;;oplot,level,y2slow,thick=3,color=mcolor*.3
;;comean1=mean(co(nhindex))
;;cosig1=stdev(co(nhindex))
;yfit = GAUSSFIT(level,y2slow, coeff, nterms=5)
;OPLOT, level, yfit, THICK=3,color=mcolor*.3
;;
;nhindex=where(y2d gt 60.)
;y2=histogram(co(nhindex),min=min(level),max=max(level),binsize=level(1)-level(0))/float(nn)
;y2shigh=smooth(100.*y2,3)
;;;oplot,level,y2shigh,thick=3,color=mcolor*.9
;yfit = GAUSSFIT(level,y2shigh, coeff, nterms=5)
;OPLOT, level, yfit, THICK=3,color=mcolor*.9
index=where(yfit eq max(yfit))
;comax=level(index(0))
;;comean=mean(co(nhindex))
;;cosig=stdev(co(nhindex))

endfor	; loop over theta
stop
;
icount=icount+1L
jumpstep:
endfor	; loop over files
;
; save postscript version
;
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename='zt_chemical_dist_'+syear(iyear)+smon(imon)+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
         !p.thick=2.0                   ;Plotted lines twice as thick
         !p.charsize=2.0
      endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=21L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
level=0.05*findgen(nlvls)
if imax gt 2.0 then level=imin+((2.*imax-imin)/float(nlvls))*findgen(nlvls)
contour,zt_chemd,1+findgen(nfile),th,levels=level,/noeras,charsize=2,c_color=col1,/cell_fill,title=syear(iyear)+smon(imon),ytitle='Theta (K)',color=0
contour,zt_chemd,1+findgen(nfile),th,levels=level,/noeras,charsize=2,color=0,/foll,/overplot

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim zt_chemical_dist_'+syear(iyear)+smon(imon)+'.ps -rotate -90 zt_chemical_dist_'+syear(iyear)+smon(imon)+'.jpg'
endif

jumpmon:
endfor	; loop over months
endfor	; loop over years
end
