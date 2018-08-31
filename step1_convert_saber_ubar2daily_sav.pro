;
; read SABER temperature, geopotential height, and gradient winds 
; (use Ubar, U', and V' from Ruth Lieberman to reconstruct 3-D data
; and store in IDL save format 
;
; ftp ftp.cora.nwra.com
; login: anonymous
; password: Your email
; cd pub/ruth
; get GRID_PHI_WINDS_2003_016_040.sav
;
; restore SABER data
;
dir='/aura6/data/SABER_data/Datfiles_winds/'
restore,dir+'GRID_PHI_WINDS_2003_016_040.sav'
;
; The dimensions and coordinates are:
;
SOUTH    = -85
NORTH    =  85
LATINC   =   5.
NUMLAT   = FIX( (NORTH-SOUTH)/LATINC) + 1
nr=numlat
LATITUDE = SOUTH + FINDGEN(NUMLAT)*LATINC
alat=latitude
LONINC   =  30.
NUMLON    = FIX( 360./LONINC)
nc=numlon
LONGITUDE = FINDGEN(NUMLON)*LONINC
alon=longitude
MAXDAYS    = 30
NUMDAYS    = 25
MAXWAVES   = 7
NUMZ      = 120   ; Maximum number of levels on a fixed altitude grid
nz=numz
ZINC      = .15
ZBASE     = 2.
ZGRID     = ZBASE + FINDGEN(NUMZ)*ZINC	; Vertical coordinate is x = -(log(p/p0) so P = P0 * exp(-zgrid)
pzero=1000.
press=pzero*exp(-1.*zgrid)
;
; zm_atemp_rec(*,*,0) and uwind(*,*,0) correspond to January 16th 
; zm_atemp_rec(*,*,24) and uwind(*,*,24) correspond to February 9th 
;
; ZM_ATEMP_REC=FLTARR(NUMLAT,NUMZ,numdays)         Tbar
; TPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     T in K
; UPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     u' in m/s
; VPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     v' in m/s
; GPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     geop height in m
; UWIND    =FLTARR(NUMLAT,NUMZ,MAXDAYS)            Ubar in m/s
;
yyyymmdd=[$
'20030116',$
'20030117',$
'20030118',$
'20030119',$
'20030120',$
'20030121',$
'20030122',$
'20030123',$
'20030124',$
'20030125',$
'20030126',$
'20030127',$
'20030128',$
'20030129',$
'20030130',$
'20030131',$
'20030201',$
'20030202',$
'20030203',$
'20030204',$
'20030205',$
'20030206',$
'20030207',$
'20030208',$
'20030209'$
]
;
; loop over days
;
for iday=0L,numdays-1L do begin
    date=yyyymmdd(iday)
;
; extract day of SABER data from these arrays
;
; ZM_ATEMP_REC=FLTARR(NUMLAT,NUMZ,numdays)         Tbar
; TPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     T in K
; UPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     u' in m/s
; VPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     v' in m/s
; GPROFILES=FLTARR(NUMZ,NUMLON,NUMLAT,MAXDAYS)     geop height in m
; UWIND    =FLTARR(NUMLAT,NUMZ,MAXDAYS)            Ubar in m/s

    v3da=reform(VPROFILES(*,*,*,iday),numz,numlon,numlat)
    z3da=reform(GPROFILES(*,*,*,iday),numz,numlon,numlat)
    upr3d=reform(UPROFILES(*,*,*,iday),numz,numlon,numlat)
    ubar=reform(UWIND(*,*,iday),numlat,numz)
    tpr3d=reform(TPROFILES(*,*,*,iday),numz,numlon,numlat)
    tbar=reform(ZM_ATEMP_REC(*,*,iday),numlat,numz)
;
; build 3d T and U based on zonal mean and eddy components (T=Tbar+T' and U=Ubar+U')
;
    u3d=fltarr(numlon,numlat,numz)
    t3d=fltarr(numlon,numlat,numz)
    v3d=fltarr(numlon,numlat,numz)
    z3d=fltarr(numlon,numlat,numz)
    for k=0L,numz-1L do begin
    for j=0L,numlat-1L do begin
    for i=0L,numlon-1L do begin
        u3d(i,j,k)=ubar(j,k)+upr3d(k,i,j)
        t3d(i,j,k)=tbar(j,k)+tpr3d(k,i,j)
        v3d(i,j,k)=v3da(k,i,j)
        z3d(i,j,k)=z3da(k,i,j)
    endfor
    endfor
    endfor
    print,'Read SABER data on ',date,' ',min(u3d),max(u3d),min(v3d),max(v3d),min(t3d),max(t3d),min(z3d),max(z3d)
;
; save 3d SABER data
;
    ofile=dir+'GRID_PHI_WINDS.'+date+'.sav'
    save,file=ofile,alon,alat,press,u3d,v3d,z3d,t3d
endfor
end
