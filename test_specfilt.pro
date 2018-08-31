@specfilt
nr=72
nc=96
lat=-88.75+180.*findgen(nr)/nr
lon=1.875+360.*findgen(nc)/nc
data=fltarr(nc,nr)
for j=0,nr-1 do begin
data(*,j)=cos(4.*lon*!pi/180.)
endfor
data0=data
specfilt,data,lon,lat,nr,nc
data1=data
 for i=0,nc-1 do begin
  dummy1=reform(data(i,*),nr)
  dummy2=dummy1
  for j=1,nr-2 do begin
   dummy2(j)=(dummy1(j-1)+2.*dummy1(j)+dummy1(j+1))/3.
  endfor
;  dummy1=dummy2
;  for j=1,nr-2 do begin
;   dummy2(j)=(-1.*dummy1(j-1)+5.*dummy1(j)-1.*dummy1(j+1))/5.
;  endfor
  data(i,*)=dummy2
 endfor
stop
;
; smooth in longitude
;
 for j=0,nr-1 do begin
  dummy1=reform(data(*,j),nc)
  dummy1=reform([dummy1(nc-2:nc-1),dummy1,dummy1(0:1)],nc+4)
  dummy2=dummy1
  for i=1,nc+2 do begin
   dummy2(i)=(dummy1(i-1)+2.*dummy1(i)+dummy1(i+1))/3.
  endfor
;  dummy1=dummy2
;  for i=1,nc+2 do begin
;   dummy2(i)=(-1.*dummy1(i-1)+5.*dummy1(i)-1.*dummy1(i+1))/5.
;  endfor
  data(*,j)=dummy2(2:nc+1)
 endfor
stop
end
