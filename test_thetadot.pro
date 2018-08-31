

@stddat
@kgmt
@ckday
@kdate

num=0L
nlvls=20
loadct,38
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
col1=1+indgen(nlvls)*icolmax/nlvls
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.10,0.10]
yorig=[0.60,0.20]
xlen=0.8
ylen=0.3
cbaryoff=0.03
cbarydel=0.01
mon=strarr(12)*4
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
lstmn=0l
lstdy=0l
lstyr=0l
ledmn=0l
leddy=0l
ledyr=0l
thlev=0l
lstday=0l
ledday=0l
date=''
print, ' '
print, '      UKMO Version '
print, ' '
read,'Enter date (month, day, year) ',lstmn,lstdy,lstyr
read,'Enter ending date (month, day, year) ',ledmn,leddy,ledyr
print, ' '
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 or lstyr gt 1999 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 1999 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   wdelete,1
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='ukmo_marks.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
   xsize=xsize,ysize=ysize
endif

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z=kgmt(imn,idy,iyr,iday)
iday = iday - 1

; loop over days
jump: iday = iday + 1
kdate,float(iday),iyr,imn,idy
ckday,iday,iyr
iyr1=iyr-1900
print,imn,idy,iyr1

; Test for end condition and close windows.
z=stddat(imn,idy,iyr,ndays)
if ndays lt lstday then stop,' starting day outside range '
if ndays gt ledday then goto, jumpout 

date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
     month(imn-1),' ',idy,', ',iyr))
if idy ge 10 then datelab=strcompress(string(FORMAT='(A4,I2,A1,I4)',$
     mon(imn-1),idy,'_',iyr))
if idy lt 10 then datelab=strcompress(string(FORMAT='(A4,A1,I1,A1,I4)',$
     mon(imn-1),'0',idy,'_',iyr))
;
; read input
;
nc=0l
nr=0l 
nth=0l 
close,10
openr,10,datelab+'_theta_dot.dat'
readu,10,nr,nc,nth
xlat=fltarr(nr)
xlon=fltarr(nc)
th=fltarr(nth)
readu,10,xlat,xlon,th
thd2=fltarr(nr,nc,nth)
q2=fltarr(nr,nc,nth)
p2=fltarr(nr,nc,nth)
readu,10,thd2,q2,p2
if ndays eq lstday then begin
icount=0
nmax=ledday-lstday+1
thd=fltarr(nr,nc,nth,nmax)
p=fltarr(nr,nc,nth,nmax)
endif
q2=0.
thd(*,*,*,icount)=thd2
thd2=0.
p(*,*,*,icount)=p2
icount=icount+1
goto, jump	; loop over days
jumpout:
;
; smooth in time
;
for l=0,nth-1 do begin
for i=0,nc-1 do begin
for j=0,nr-1 do begin
dummy1=reform(p(j,i,l,*),nmax)
dummy2=reform(thd(j,i,l,*),nmax)
index=where(dummy1 ge 70.)
if index(0) ne -1 then begin
dummy2=smooth(dummy2,3,/edge_truncate)
dummy2=smooth(dummy2,3,/edge_truncate)
dummy2=smooth(dummy2,3,/edge_truncate)
dummy2=smooth(dummy2,3,/edge_truncate)
dummy2=smooth(dummy2,3,/edge_truncate)
thd(j,i,l,*)=dummy2
endif
endfor
endfor
endfor

end
