nr=96
dummy1=cos(2.*!pi*6.*findgen(96)/96.)
dummy2=dummy1
plot,dummy1
stop
  for j=1,nr-2 do begin
   dummy2(j)=(dummy1(j-1)+2.*dummy1(j)+dummy1(j+1))/5.
  endfor
  dummy1=dummy2
  oplot,dummy1
stop
  for j=1,nr-2 do begin
   dummy2(j)=(-1.*dummy1(j-1)+5.*dummy1(j)-1.*dummy1(j+1))/3.
  endfor
  dummy1=dummy2
  oplot,dummy1
stop
end
