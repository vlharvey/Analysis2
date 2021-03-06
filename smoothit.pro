pro smoothit,cn,cnsmooth
result=size(cn)
nc=result(1)
nr=result(2)
cnsm1=cn
for j=0,nr-1 do begin
    dummy=reform(cnsm1(*,j),nc)
    index=where(dummy ne -9999.,ngood)
    if ngood gt 1 then $
       cnsm1(index,j)=smooth(dummy(index),3)
endfor
cnsm2=cn
for i=0,nc-1 do begin
    dummy=reform(cnsm2(i,*),nr)
    index=where(dummy ne -9999.,ngood)
    if ngood gt 1 then $
       cnsm2(i,index)=smooth(dummy(index),3)
endfor
cnsmooth=.5*(cnsm1+cnsm2)
;
; remove bad values introduced by interpol
;
index=where(cnsm1 lt 0. or cnsm1 gt 15. or $
            cnsm2 lt 0. or cnsm2 gt 15.)
if index(0) ne -1 then cnsmooth(index)=-9999.
return
end
