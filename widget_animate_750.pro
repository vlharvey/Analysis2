@xinteranimate.pro
nfram=220
ifile=' '
read,'input file?',ifile
read,'input number of frames ',nframe
OPENR, 1, ifile
                H = BYTARR(750,750, nframe)
                READU, 1, H
                CLOSE, 1

print, 'prepare to initailize XINTERANIMATE'
                XinterANIMATE, SET=[750,750, nframe], /Showload

print, ' XinterANIMATE initailized'
print, 'prepare to load the images into XinterANIMATE'

                FOR I=0,nframe-1 DO begin
                XinterANIMATE, FRAME = I, IMAGE = H(*,*,I)
                endfor
xinteranimate,1
end


