;
; save JJA and DJF zonal mean T, U, V, Z, CO, N2O, H2O, and O3 for MLS record
;
; KDAY            FLOAT     =       4664.00
; NLV             LONG      =           55
; NLV2            LONG      =           37
; NR              LONG      =           96
; LAT             DOUBLE    = Array[96]
; PRESS37         FLOAT     = Array[37]
; PRESS55         FLOAT     = Array[55]
; SDATE_ALL       STRING    = Array[4664]
; TBAR            FLOAT     = Array[96, 55, 4664]
; UBAR            FLOAT     = Array[96, 55, 4664]
; VBAR            FLOAT     = Array[96, 55, 4664]
; ZBAR            FLOAT     = Array[96, 55, 4664]
; H2OBAR          FLOAT     = Array[96, 55, 4664]
; O3BAR           FLOAT     = Array[96, 55, 4664]
; N2OBAR          FLOAT     = Array[96, 37, 4664]
; COBAR           FLOAT     = Array[96, 37, 4664]
;
restore,'/atmos/aura6/data/MLS_data/Pre_process/mls_daily_zonal_means_from_gridded_20040808-20170515.sav
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2) 
jja=where(smon eq '06' or smon eq '07' or smon eq '08')
djf=where(smon eq '12' or smon eq '01' or smon eq '02')
;
; JJA (1128 days)
;
TBAR_JJA=mean(TBAR(*,*,jja),dim=3,/Nan)
UBAR_JJA=mean(UBAR(*,*,jja),dim=3,/Nan)
VBAR_JJA=mean(VBAR(*,*,jja),dim=3,/Nan)
ZBAR_JJA=mean(ZBAR(*,*,jja),dim=3,/Nan)
H2OBAR_JJA=mean(H2OBAR(*,*,jja),dim=3,/Nan)
O3BAR_JJA=mean(O3BAR(*,*,jja),dim=3,/Nan)
N2OBAR_JJA=mean(N2OBAR(*,*,jja),dim=3,/Nan)
COBAR_JJA=mean(COBAR(*,*,jja),dim=3,/Nan)
;
; DJF (1173 days)
;
TBAR_DJF=mean(TBAR(*,*,djf),dim=3,/Nan)       
UBAR_DJF=mean(UBAR(*,*,djf),dim=3,/Nan)       
VBAR_DJF=mean(VBAR(*,*,djf),dim=3,/Nan)       
ZBAR_DJF=mean(ZBAR(*,*,djf),dim=3,/Nan)       
H2OBAR_DJF=mean(H2OBAR(*,*,djf),dim=3,/Nan)
O3BAR_DJF=mean(O3BAR(*,*,djf),dim=3,/Nan)
N2OBAR_DJF=mean(N2OBAR(*,*,djf),dim=3,/Nan)
COBAR_DJF=mean(COBAR(*,*,djf),dim=3,/Nan)
;
; save seasonal means
;
save,file='mls_djf_jja.sav',lat,PRESS37,PRESS55,TBAR_JJA,UBAR_JJA,VBAR_JJA,ZBAR_JJA,H2OBAR_JJA,O3BAR_JJA,$
     N2OBAR_JJA,COBAR_JJA,TBAR_DJF,UBAR_DJF,VBAR_DJF,ZBAR_DJF,H2OBAR_DJF,O3BAR_DJF,N2OBAR_DJF,COBAR_DJF
end
