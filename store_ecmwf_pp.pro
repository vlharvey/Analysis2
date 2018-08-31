pro store_ecmwf_pp,filename,nlg,nlat,alon,alat,mslp,sfcp

print, 'opening '+filename

close,10
openw,10,filename,/f77
writeu,10,nlg,nlat
writeu,10,alon,alat
writeu,10,mslp,sfcp
close,10
end

