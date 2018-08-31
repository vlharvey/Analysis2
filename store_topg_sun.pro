;
; Read ASCII 2.5x3.75 topograpy data and store in binary on a Sun.
;
close,10
openr,10,'topg3.75.ascii'
nc=0L & nr=0L
readf,10,nc,nr
topg=fltarr(nc,nr)
readf,10,topg
close,10
openw,10,'topg3.75'
writeu,10,topg
close,10
end
