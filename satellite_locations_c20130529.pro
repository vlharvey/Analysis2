;produce satellite position/time data file
output_path='Z:\Brakebusch\WACCM\SAT_HIST\'

;;;;; available data (as of 20130227) - IMPORTANT: UPDATE THOSE LINES BEFORE RUNNING !!!!!

;live instruments (or instruments that need a last update)
MLS_period = [20040808,20130224] ;uses original MLS he5 O3 files [last update 20130226]
ACE_FTS_period = [20040110,20120318] ;uses yearly ACE-FTS nc files (produced using ACE_raw2netcdf.pro) [last updated 20130226]
SABER_period = [20020125,20130201] ;uses original SABER data coincidence tables [last update 20130226]
SOFIE_period = [20070515,20130128] ;uses origianl SOFIE Event Date & Location Files [last update 20130326]
MIPAS_period = [20020605,20120408] ;uses IDL save files created by Laura [last update 20130227]
SBUV_period = [19700410,20111231] ;uses catalog files created from original v8.6 BUV, SBUV, and SBUV2 files [last update 20130227]
HIRDLS_period = [20050122,20080317] ;uses original HIRDLS he5 files [last update 20130227]
SMR_period = [20010619,20130224] ;uses ascii dump from SMR database [last update 20130227]
GOMOS_period = [20020826,20110524] ;uses IDL save files created by Lynn [last update 20130226, new data version is more recent but not available yet]

;dead instruments that don't need an update (look for new data versions though!)
CONCORDIASI_PSC16_period = [20100910,20101003] ;uses original CONCORDIASI ascii files ;NO UPDATE NECESSARY
CONCORDIASI_PSC17_period = [20100913,20101014] ;uses original CONCORDIASI ascii files ;NO UPDATE NECESSARY
CONCORDIASI_PSC19_period = [20101003,20101122] ;uses original CONCORDIASI ascii files ;NO UPDATE NECESSARY
POAM2_period = [19931015,19961113] ;uses IDL save files (cat files) created by Lynn ;NO UPDATE NECESSARY
POAM3_period = [19980422,20051205] ;uses IDL save files (cat files) created by Lynn ;NO UPDATE NECESSARY
SAGE1_period = [19790221,19811118] ;uses SABER1 catalog file (produced using read_sage1_ascii.pro) ;NO UPDATE NECESSARY 
SAGE2_period = [19841024,20050822] ;uses IDL save files (cat files) created by Lynn ;NO UPDATE NECESSARY
SAGE3_period = [20020227,20051231] ;uses IDL save files (cat files) created by Lynn ;NO UPDATE NECESSARY
HALOE_period = [19911011,20051121] ;uses IDL save files (cat files) created by Lynn ;NO UPDATE NECESSARY
PEX2005_period = [20050307,20050307] ;uses ascii table created by Tobias (one file for PEX2005,REC2010,PEX2010) ;NO UPDATE NECESSARY
REC2010_period = [20100117,20100305] ;uses ascii table created by Tobias (one file for PEX2005,REC2010,PEX2010) ;NO UPDATE NECESSARY
PEX2010_period = [20100310,20100310] ;uses ascii table created by Tobias (one file for PEX2005,REC2010,PEX2010) ;NO UPDATE NECESSARY
UARS_MLS_period = [19910919,19990727] ;uses original V0004, V0005, V0006 3AT files ;NO UPDATE NECESSARY
SMILES_period = [20091012,20100421] ;uses ASCII file created by Doug ;NO UPDATE NECESSARY
LIMS_period = [19781025,19790528] ;uses original V006 files ;NO UPDATE NECESSARY
CONCORDIASI_PSC14_period = [20080713,20101221] ;uses original CONCORDIASI files (french data) ;NO UPDATE NECESSARY
CONCORDIASI_PSC15_period = [20090101,20100914] ;uses original CONCORDIASI files (french data) ;NO UPDATE NECESSARY
CONCORDIASI_PSC18_period = [20090101,20101130] ;uses original CONCORDIASI files (french data) ;NO UPDATE NECESSARY
SME_period = [19811215,19861121] ;uses NetCDF files prepared by Aimee Merkel ;NO UPDATE NECESSARY

mls_obs_path='Z:\Brakebusch\MLS3\original he5\'
ace_obs_path='Z:\Brakebusch\SOSST\ACE original nc\'
;hirdls_obs_path='Z:\Brakebusch\HIRDLS\original he5\'
;hirdls_obs_path='Z:\Brakebusch\HIRDLS\v5.05\'
hirdls_obs_path='Z:\Brakebusch\HIRDLS\v6\'
saber_obs_path='Z:\Brakebusch\SABER\original ascii\'
sofie_obs_path='Z:\Brakebusch\SOFIE\'
concordiasi_obs_path = 'Z:\Brakebusch\CONCORDIASI\'
mipas_obs_path='Z:\Brakebusch\MIPAS\'
poam2_obs_path='Z:\Brakebusch\SOSST\POAM_cat\' ;this path contains cat !!!
poam3_obs_path='Z:\Brakebusch\SOSST\POAM_cat\' ;this path contains cat !!!
sage1_obs_path='Z:\Brakebusch\SOSST\SAGE1\original (HDF)\'
sage2_obs_path='Z:\Brakebusch\SOSST\SAGE2\'
sage3_obs_path='Z:\Brakebusch\SOSST\SAGE3_cat\' ;this path contains cat !!!
haloe_obs_path='Z:\Brakebusch\SOSST\HALOE\'
juelich_obs_path='Z:\Brakebusch\JUELICH\'
uars_mls_obs_path='Z:\Brakebusch\UARS_MLS\UARML3AT\'
smiles_obs_path='Z:\Brakebusch\SMILES\'
lims_obs_path='Z:\Brakebusch\LIMS\'
sbuv_obs_path='Z:\Brakebusch\SBUV\catalog_files\'
smr_obs_path='Z:\Brakebusch\ODIN\'
sme_obs_path='Z:\Brakebusch\SME\'
gomos_obs_path='Z:\Brakebusch\GOMOS\'

;WACCM grid
waccm_lon = indgen(144)/144.*360. & waccm_lon_step = 2.5
waccm_lat = indgen(96)/95.*180.-90. & waccm_lat_step = 1.89475

period_array = [mls_period, ace_fts_period, hirdls_period, saber_period, sofie_period, concordiasi_psc16_period, concordiasi_psc17_period,$
                concordiasi_psc19_period, mipas_period, poam2_period, poam3_period, sage2_period, sage3_period, haloe_period, pex2005_period,$
                rec2010_period, pex2010_period, concordiasi_psc14_period, concordiasi_psc15_period, concordiasi_psc18_period, sage1_period,$
                uars_mls_period, smiles_period, lims_period, sbuv_period] 
start_date = strcompress(min(period_array),/r)
end_date   = strcompress(max(period_array),/r)
;start_date = '20070514'
;end_date = '20130129'


saber_tpl = {VERSION:1.00000, DATASTART:3, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:8,FIELDTYPES:[3,3,7,7,3,4,4,4],$
             FIELDNAMES:['orbit','event','date','time','doy','latitude','longitude','sza'],$
             FIELDLOCATIONS:[0,6,10,23,34,50,66,82],FIELDGROUPS:[0,1,2,3,4,5,6,7]}
             
sofie_tpl = {VERSION:1.00000, DATASTART:1, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:6,FIELDTYPES:[3,3,3,4,4,4],$
             FIELDNAMES:['event','orbit','date','time','latitude','longitude'],$
             FIELDLOCATIONS:[0,4,8,17,23,32],FIELDGROUPS:[0,1,2,3,4,5]}             
             
concordiasi_tpl = {VERSION:1.00000, DATASTART:1, DELIMITER:44b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:7,FIELDTYPES:[7,4,4,4,4,4,4],$
             FIELDNAMES:['julian','ozone','latitude','longitude','altitude','temperature','pressure'],$
             FIELDLOCATIONS:[0,13,23,34,45,55,63],FIELDGROUPS:[0,1,2,3,4,5,6]}
             
juelich_tpl = {VERSION:1.00000, DATASTART:1, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:5,FIELDTYPES:[3,4,4,4,7],$
             FIELDNAMES:['date','time','lon','lat','tag'],$
             FIELDLOCATIONS:[0,9,17,25,33],FIELDGROUPS:[0,1,2,3,4]}                          
             
;aura_tpl =  {VERSION:1.00000, DATASTART:0, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
;             COMMENTSYMBOL:'',FIELDCOUNT:2,FIELDTYPES:[7,3],$
;             FIELDNAMES:['date','TAI93At0zOfGranule'],$
;             FIELDLOCATIONS:[0,9],FIELDGROUPS:[0,1]}

smiles_tpl = {VERSION:1.00000, DATASTART:1, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:8,FIELDTYPES:[3,7,4,3,4,4,4,3],$
             FIELDNAMES:['obs_nr','date_time','ltime','time','lat','lon','sza','asc_dsc'],$
             FIELDLOCATIONS:[0,2,26,32,43,50,58,64],FIELDGROUPS:[0,1,2,3,4,5,6,7]}

smiles_tpl = {VERSION:1.00000, DATASTART:1, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:8,FIELDTYPES:[3,7,4,3,4,4,4,3],$
             FIELDNAMES:['obs_nr','date_time','ltime','time','lat','lon','sza','asc_dsc'],$
             FIELDLOCATIONS:[0,2,26,32,43,50,58,64],FIELDGROUPS:[0,1,2,3,4,5,6,7]}

concordiasi2_tpl = {VERSION:1.00000, DATASTART:0, DELIMITER:32b, MISSINGVALUE:!values.f_nan,$
             COMMENTSYMBOL:'',FIELDCOUNT:31,FIELDTYPES:[4,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,3],$
             FIELDNAMES:['decimal_day','year','month','day','hour_UTC','minute_UTC','second_UTC','lon','lat','alt_m',$
             'pressure_Pa','temperature1_K','temperature2_K','SZA','u','v','corr_T1_K','corr_T2_K','GPS_mode','ISBA_press_hPa',$
             'ISBA_T1_K','ISBA_T2_K','u_GPS','v_GPS','internal_T_TSEN_C','volt_TSEN_power_V','He_press_balloon_C','He_T_balloon_C',$
             'internal_T_ISBA_C','volt_ISBA_power_V','volt_ISBA_heat_V'],$
             FIELDLOCATIONS:[1,13,19,21,24,27,31,34,46,58,65,75,87,100,110,122,140,144,153,160,165,173,181,191,194,208,214,223,229,235,242],$
             FIELDGROUPS:[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]}


ndays=date2julday(end_date) - date2julday(start_date) +1
      
start_julday = date2julday(start_date)
start_doy = strmid(start_date,0,4)+strcompress(string(date2doy(start_date),format='(i3.3)'),/r)
end_doy = strmid(end_date,0,4)+strcompress(string(date2doy(end_date),format='(i3.3)'),/r)

date_reference = 19700101
julian_reference = date2julday(date_reference)
if julian_reference gt start_julday then stop

;;read Aura time correction file (TAI93At0zOfGranule)
;file='Z:\Brakebusch\MLS\TAIat0z.txt'
;data=read_ascii(file,template=aura_tpl)
;aura_strdatecode = strmid(data.date,0,4)+'d'+strmid(data.date,5,3)
;aura_time_corr = data.TAI93At0zOfGranule

str_mls_period=strcompress(string(mls_period,format='(i8)'),/r) & str_mls_period = str_mls_period[0]+' - '+str_mls_period[1]
str_ace_fts_period=strcompress(string(ace_fts_period,format='(i8)'),/r) & str_ace_fts_period = str_ace_fts_period[0]+' - '+str_ace_fts_period[1]
str_hirdls_period=strcompress(string(hirdls_period,format='(i8)'),/r) & str_hirdls_period = str_hirdls_period[0]+' - '+str_hirdls_period[1]
str_saber_period=strcompress(string(saber_period,format='(i8)'),/r) & str_saber_period = str_saber_period[0]+' - '+str_saber_period[1]
str_sofie_period=strcompress(string(sofie_period,format='(i8)'),/r) & str_sofie_period = str_sofie_period[0]+' - '+str_sofie_period[1]
str_mipas_period=strcompress(string(mipas_period,format='(i8)'),/r) & str_mipas_period = str_mipas_period[0]+' - '+str_mipas_period[1]
str_concordiasi_psc16_period=strcompress(string(concordiasi_psc16_period,format='(i8)'),/r) & str_concordiasi_psc16_period = str_concordiasi_psc16_period[0]+' - '+str_concordiasi_psc16_period[1]
str_concordiasi_psc17_period=strcompress(string(concordiasi_psc17_period,format='(i8)'),/r) & str_concordiasi_psc17_period = str_concordiasi_psc17_period[0]+' - '+str_concordiasi_psc17_period[1]
str_concordiasi_psc19_period=strcompress(string(concordiasi_psc19_period,format='(i8)'),/r) & str_concordiasi_psc19_period = str_concordiasi_psc19_period[0]+' - '+str_concordiasi_psc19_period[1]
str_poam2_period=strcompress(string(poam2_period,format='(i8)'),/r) & str_poam2_period = str_poam2_period[0]+' - '+str_poam2_period[1]
str_poam3_period=strcompress(string(poam3_period,format='(i8)'),/r) & str_poam3_period = str_poam3_period[0]+' - '+str_poam3_period[1]
str_sage2_period=strcompress(string(sage2_period,format='(i8)'),/r) & str_sage2_period = str_sage2_period[0]+' - '+str_sage2_period[1]
str_sage3_period=strcompress(string(sage3_period,format='(i8)'),/r) & str_sage3_period = str_sage3_period[0]+' - '+str_sage3_period[1]
str_haloe_period=strcompress(string(haloe_period,format='(i8)'),/r) & str_haloe_period = str_haloe_period[0]+' - '+str_haloe_period[1]
str_pex2005_period=strcompress(string(pex2005_period,format='(i8)'),/r) & str_pex2005_period = str_pex2005_period[0]+' - '+str_pex2005_period[1]
str_rec2010_period=strcompress(string(rec2010_period,format='(i8)'),/r) & str_rec2010_period = str_rec2010_period[0]+' - '+str_rec2010_period[1]
str_pex2010_period=strcompress(string(pex2010_period,format='(i8)'),/r) & str_pex2010_period = str_pex2010_period[0]+' - '+str_pex2010_period[1]
str_concordiasi_psc14_period=strcompress(string(concordiasi_psc14_period,format='(i8)'),/r) & str_concordiasi_psc14_period = str_concordiasi_psc14_period[0]+' - '+str_concordiasi_psc14_period[1]
str_concordiasi_psc15_period=strcompress(string(concordiasi_psc15_period,format='(i8)'),/r) & str_concordiasi_psc15_period = str_concordiasi_psc15_period[0]+' - '+str_concordiasi_psc15_period[1]
str_concordiasi_psc18_period=strcompress(string(concordiasi_psc18_period,format='(i8)'),/r) & str_concordiasi_psc18_period = str_concordiasi_psc18_period[0]+' - '+str_concordiasi_psc18_period[1]
str_sage1_period=strcompress(string(sage1_period,format='(i8)'),/r) & str_sage1_period = str_sage1_period[0]+' - '+str_sage1_period[1]
str_uars_mls_period=strcompress(string(uars_mls_period,format='(i8)'),/r) & str_uars_mls_period = str_uars_mls_period[0]+' - '+str_uars_mls_period[1]
str_smiles_period=strcompress(string(smiles_period,format='(i8)'),/r) & str_smiles_period = str_smiles_period[0]+' - '+str_smiles_period[1]
str_lims_period=strcompress(string(lims_period,format='(i8)'),/r) & str_lims_period = str_lims_period[0]+' - '+str_lims_period[1]
str_sbuv_period=strcompress(string(sbuv_period,format='(i8)'),/r) & str_sbuv_period = str_sbuv_period[0]+' - '+str_sbuv_period[1]
str_smr_period=strcompress(string(smr_period,format='(i8)'),/r) & str_smr_period = str_smr_period[0]+' - '+str_smr_period[1]
str_sme_period=strcompress(string(sme_period,format='(i8)'),/r) & str_sme_period = str_sme_period[0]+' - '+str_sme_period[1]
str_gomos_period=strcompress(string(gomos_period,format='(i8)'),/r) & str_gomos_period = str_gomos_period[0]+' - '+str_gomos_period[1]


;open NetCDF file for sequential writing
ncid = ncdf_create(output_path+'sat_hist_'+start_date+'-'+end_date+'_c'+today()+'.nc',/clobber)
  ncdf_attput, ncid, 'Description', 'Satellite profile positions and times from '+start_date+' through '+end_date, /global
  ncdf_attput, ncid, 'Author', 'Matthias Brakebusch, file was created on '+today()+' using satellite_locations_c20130326.pro', /global
  ncdf_attput, ncid, 'Comment', '(1) EOS MLS O3 he5 files version v03-30 were used (leap seconds not corrected!), '+str_mls_period+'; '+$
                                '(2) ACE-FTS ascii files version v3.0 were used, '+str_ace_fts_period+'; '+$
                                '(3) EOS HIRDLS he5 files version v5.05 were used, '+str_hirdls_period+'; '+$
                                '(4) SABER ascii files version v1.07 were used, '+str_saber_period+'; '+$
                                '(5) SOFIE eventLocation files v1.2 were used, '+str_sofie_period+'; '+$
                                '(6) CONCORDIASI PSC 16 ascii file version ? was used, '+str_concordiasi_psc16_period+'; '+$
                                '(7) CONCORDIASI PSC 17 ascii file version ? was used, '+str_concordiasi_psc17_period+'; '+$
                                '(8) CONCORDIASI PSC 19 ascii file version ? was used, '+str_concordiasi_psc19_period+'; '+$
                                '(9) MIPAS E_IMK O3 ascii files version 30, 40, 5R were used, '+str_mipas_period+'; '+$
                                '(10) POAM2 files version 6.0 were used, '+str_poam2_period+'; '+$
                                '(11) POAM3 files version 4.0 were used, '+str_poam3_period+'; '+$
                                '(12) SAGE2 files version 6.2 were used, '+str_sage2_period+'; '+$
                                '(13) SAGE3 files version 3.00 were used, '+str_sage3_period+'; '+$
                                '(14) HALOE files version 19 were used, '+str_haloe_period+'; '+$
                                '(15) Juelich PEX2005 file version ? was used, '+str_pex2005_period+'; '+$
                                '(16) RECONCILE REC2010 file version ? was used, '+str_rec2010_period+'; '+$
                                '(17) Juelich PEX2010 file version ? was used, '+str_pex2010_period+'; '+$
                                '(18) SAGE1 files version LaRC 3/4/6 were used, '+str_sage1_period+'; '+$
                                '(19) UARS MLS L3AT files version V0004, V0005, V0006 were used, '+str_uars_mls_period+'; '+$
                                '(20) SMILES file version ? was used, '+str_smiles_period+'; '+$
                                '(21) LIMS files version 6 were used, '+str_lims_period+'; '+$
                                '(22) CONCORDIASI PSC 14 ascii file version 1105 was used, '+str_concordiasi_psc14_period+'; '+$
                                '(23) CONCORDIASI PSC 15 ascii file version 1105 was used, '+str_concordiasi_psc15_period+'; '+$
                                '(24) CONCORDIASI PSC 18 ascii file version 1105 was used, '+str_concordiasi_psc18_period+'; '+$
                                '(25) BUV on Nimbus 4 catalog file based on profile data v8.6 was used, 19700410 - 19760430; '+$
                                '(26) SBUV on Nimbus 7 catalog file based on profile data v8.6 was used, 19781031 - 19900531; '+$
                                '(27) SBUV2 on NOAA 9 catalog file based on profile data v8.6 was used, 19850202 - 19980131; '+$
                                '(28) SBUV2 on NOAA 11 catalog file based on profile data v8.6 was used, 19881201 - 20010327; '+$
                                '(29) SBUV2 on NOAA 14 catalog file based on profile data v8.6 was used, 19950205 - 20060928; '+$
                                '(30) SBUV2 on NOAA 16 catalog file based on profile data v8.6 was used, 20001003 - 20111231; '+$ ;update this line
                                '(31) SBUV2 on NOAA 17 catalog file based on profile data v8.6 was used, 20020711 - 20111231; '+$ ;update this line
                                '(32) SBUV2 on NOAA 18 catalog file based on profile data v8.6 was used, 20050605 - 20111231; '+$ ;update this line
                                '(33) SBUV2 on NOAA 19 catalog file based on profile data v8.6 was used, 20090223 - 20111231; '+$ ;update this line
                                '(34) SMR geolocations from ODIN/SMR datbase were used, '+str_smr_period+'; '+$
                                '(35) SME AG (airglow) observations from March 2002 were used, '+str_sme_period+'; '+$
                                '(36) SME UV (ultraviolet) observations from October 2001 were used, '+str_sme_period+'; '+$
                                '(37) GOMOS files version ??? were used, '+str_gomos_period+'.', /global
  ncdf_attput, ncid, 'prof_num', '(1) MLS: OrbitGeodeticAngle*10; '+$
                                 '(2) ACE-FTS: OccultationNumber; '+$
                                 '(3) HIRDLS: Profile_ID; '+$
                                 '(4) SABER: Event; '+$
                                 '(5) SOFIE: event; '+$
                                 '(6) CONCORDIASI PSC 16: 16*100000+measurement*100+neighbor; '+$
                                 '(7) CONCORDIASI PSC 17: 17*100000+measurement*100+neighbor; '+$
                                 '(8) CONCORDIASI PSC 19: 19*100000+measurement*100+neighbor; '+$
                                 '(9) MIPAS: days_since_20020101*100000+sid; '+$
                                 '(10) POAM2: ID as in Lynn Harvey''s cat files; '+$
                                 '(11) POAM3: ID as in Lynn Harvey''s cat files; '+$
                                 '(12) SAGE2: ID as in Lynn Harvey''s cat files; '+$
                                 '(13) SAGE3: ID as in Lynn Harvey''s cat files; '+$
                                 '(14) HALOE: ID as in Lynn Harvey''s cat files; '+$
                                 '(15) PEX2005: year*100000+flight*1000+observation; '+$
                                 '(16) REC2010: year*100000+flight*1000+observation; '+$
                                 '(17) PEX2010: year*100000+flight*1000+observation; '+$
                                 '(18) SAGE1: number of observation; '+$
                                 '(19) UARS MLS: number of observation; '+$
                                 '(20) SMILES: observation number; '+$
                                 '(21) LIMS: number of observation; '+$
                                 '(22) CONCORDIASI PSC 14: 14*100000+measurement*100+neighbor; '+$
                                 '(23) CONCORDIASI PSC 15: 15*100000+measurement*100+neighbor; '+$
                                 '(24) CONCORDIASI PSC 18: 18*100000+measurement*100+neighbor; '+$
                                 '(25) BUV on Nimbus 4: 04*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(26) SBUV on Nimbus 7: 07*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(27) SBUV2 on NOAA 9: 09*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(28) SBUV2 on NOAA 11: 11*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(29) SBUV2 on NOAA 14: 14*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(30) SBUV2 on NOAA 16: 16*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(31) SBUV2 on NOAA 17: 17*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(32) SBUV2 on NOAA 18: 18*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(33) SBUV2 on NOAA 19: 19*100000000+day_of_instrument*10000+obs_of_day; '+$
                                 '(34) SMR: days_since_20010619*100000+seconds_in_day; '+$
                                 '(35) SME AG: days_since_19811215*100000+seconds_in_day; '+$
                                 '(36) SME UV: days_since_19811215*100000+seconds_in_day; '+$
                                 '(37) GOMOS: days_since_20020101*100000+seconds_in_day.', /global
  tid = ncdf_dimdef(ncid, 'profs', /unlimited)
  
  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'orbit_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'prof_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'instr_num', tid, /LONG)
  
  vid = ncdf_vardef(ncid, 'doy', tid, /LONG)
  vid = ncdf_vardef(ncid, 'local_time', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'occ_type', tid, /SHORT)
  vid = ncdf_vardef(ncid, 'instr_sza', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'julian', tid, /DOUBLE)
  
  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  ncdf_attput, ncid, 'orbit_num', 'long_name', 'orbit number' & ncdf_attput, ncid, 'orbit_num', 'units', '1'
  ncdf_attput, ncid, 'prof_num', 'long_name', 'profile number, ... see ''prof_num'' comment' & ncdf_attput, ncid, 'prof_num', 'units', '1'
  ncdf_attput, ncid, 'instr_num', 'long_name', 'MLS=1, ACE-FTS=2, HIRDLS=3, ... see global comment' & ncdf_attput, ncid, 'instr_num', 'units', '1'
  
  ncdf_attput, ncid, 'doy', 'long_name', 'year, day of year' & ncdf_attput, ncid, 'doy', 'units', 'yyyyddd'
  ncdf_attput, ncid, 'local_time', 'long_name', 'local solar time' & ncdf_attput, ncid, 'local_time', 'units', 'hours'
  ncdf_attput, ncid, 'occ_type', 'long_name', 'type of occultation' & ncdf_attput, ncid, 'occ_type', 'units', '1 = sunrise, -1 = sunset, 0 = N/A'
  ncdf_attput, ncid, 'instr_sza', 'long_name', 'solar zenith angle' & ncdf_attput, ncid, 'instr_sza', 'units', 'degrees'
  ncdf_attput, ncid, 'julian', 'long_name', 'julian date' & ncdf_attput, ncid, 'julian', 'units', 'days since '+strcompress(date_reference,/r)
  
  ncdf_control, ncid, /endef

  ;sequentially read data and write to NetCDF
  offset=0l
  ace_old_file=''
  saber_old_file=''
  hirdls_old_file=''
  sofie_sr_old_file=''
  sofie_ss_old_file=''
  mipas_old_file=''
  concordiasi_psc16_old_file=''
  concordiasi_psc17_old_file=''
  concordiasi_psc19_old_file=''
  poam2_old_file=''
  poam3_old_file=''
  sage2_old_file=''
  sage3_old_file=''
  haloe_old_file=''
  juelich_old_file=''
  sage1_old_file=''
  uars_mls_old_file=''
  smiles_old_file=''
  lims_old_file=''
  concordiasi_psc14_old_file=''
  concordiasi_psc15_old_file=''
  concordiasi_psc18_old_file=''
  sbuv_old_file=''
  smr_old_file=''
  sme_ag_old_file=''
  sme_uv_old_file=''
  gomos_old_file=''
  openw,2,output_path+'satellite_locations_per_day_c'+today()+'.txt'
  printf,2,'date  nprofs_orig nprofs_incl MLS_orig MLS_incl ACE-FTS_orig ACE-FTS_incl HIRDLS_orig HIRDLS_incl '+$
           'SABER_orig SABER_incl SOFIE_orig SOFIE_incl CONCORDIASI_PSC16_orig CONCORDIASI_PSC16_incl CONCORDIASI_PSC17_orig '+ $
           'CONCORDIASI_PSC17_incl CONCORDIASI_PSC19_orig CONCORDIASI_PSC19_incl MIPAS_orig MIPAS_incl POAM2_orig '+$
           'POAM2_incl POAM3_orig POAM3_incl SAGE2_orig SAGE2_incl SAGE3_orig SAGE3_incl HALOE_orig HALOE_incl '+ $
           'PEX2005_orig PEX2005_incl REC2010_orig REC2010_incl PEX2010_orig PEX2010_incl SAGE1_orig SAGE1_incl '+ $
           'UARS_MLS_orig UARS_MLS_incl SMILES_orig SMILES_incl LIMS_orig LIMS_incl CONCORDIASI_PSC14_orig '+$
           'CONCORDIASI_PSC14_incl CONCORDIASI_PSC15_orig CONCORDIASI_PSC15_incl CONCORDIASI_PSC18_orig CONCORDIASI_PSC18_incl '+ $
           'BUV_orig BUV_incl SBUV_orig SBUV_incl SBUV2_orig SBUV2_incl SMR_orig SMR_incl SME_AG_orig SME_AG_incl SME_UV_orig SME_UV_incl '+ $
           'GOMOS_orig GOMOS_incl'

  for iday=0l,ndays-1l do begin
    print, iday+1, ndays
  
    ;produce date variables
    current_julday = start_julday + iday
    caldat, current_julday, month, day, year
    current_doy = current_julday-julday(1,1,year)+1
    current_date = long(year*10000l+month*100l+day)
    str_year = strcompress(string(year),/r)
    str_date_code = date2datecode(julday2date(start_julday + iday))
    seconds_since_19930101=86400*(current_julday-julday(1,1,1993)) ;for MLS
    ;print, str_date_code
        
    ;read POAM2
    file=file_search(poam2_obs_path+'cat_poam2_v6.0.'+str_year,count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with POAM2 data anyway
    if nfiles ne 0 then begin
      if file ne poam2_old_file then begin 
        restore, file[0]
        poam2_orig_date = date
        poam2_orig_id = id
        poam2_orig_latitudes = latitude
        poam2_orig_longitudes = longitude
        poam2_orig_occ_type = fix(1*(sctype eq 'r')-1*(sctype eq 's'))
        poam2_orig_time = time
        poam2_old_file=file
        undefine, date & undefine, id & undefine, latitude & undefine, longitude & undefine, sctype & undefine, time 
        print,'starting POAM2 '+str_year
      endif
      index=where(poam2_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        poam2_nprofs = nindex
        poam2_longitudes = poam2_orig_longitudes[index]
        poam2_latitudes = poam2_orig_latitudes[index]
        poam2_time = poam2_orig_time[index]
        poam2_occ_num = poam2_orig_id[index] 
        poam2_occ_type = poam2_orig_occ_type[index]
      endif else poam2_nprofs = 0
    endif else poam2_nprofs = 0
    
    ;read POAM3
    file=file_search(poam3_obs_path+'cat_poam3_v4.0.'+str_year,count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with POAM3 data anyway
    if nfiles ne 0 then begin
      if file ne poam3_old_file then begin 
        restore, file[0]
        poam3_orig_date = date
        poam3_orig_id = id
        poam3_orig_latitudes = latitude
        poam3_orig_longitudes = longitude
        poam3_orig_occ_type = fix(1*(sctype eq 'r')-1*(sctype eq 's'))
        poam3_orig_time = time
        poam3_old_file=file
        undefine, date & undefine, id & undefine, latitude & undefine, longitude & undefine, sctype & undefine, time 
        print,'starting POAM3 '+str_year
      endif
      index=where(poam3_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        poam3_nprofs = nindex
        poam3_longitudes = poam3_orig_longitudes[index]
        poam3_latitudes = poam3_orig_latitudes[index]
        poam3_time = poam3_orig_time[index]
        poam3_occ_num = poam3_orig_id[index] 
        poam3_occ_type = poam3_orig_occ_type[index]
      endif else poam3_nprofs = 0
    endif else poam3_nprofs = 0
    
    ;read SAGE2
    file=file_search(sage2_obs_path+'cat_sage2_v6.2.'+str_year,count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with SAGE2 data anyway
    if nfiles ne 0 then begin
      if file ne sage2_old_file then begin 
        restore, file[0]
        sage2_orig_date = date
        sage2_orig_id = id
        sage2_orig_latitudes = latitude
        sage2_orig_longitudes = longitude
        sage2_orig_occ_type = fix(1*(sctype eq 'r')-1*(sctype eq 's'))
        sage2_orig_time = time
        sage2_old_file=file
        undefine, date & undefine, id & undefine, latitude & undefine, longitude & undefine, sctype & undefine, time 
        print,'starting SAGE2 '+str_year
      endif
      index=where(sage2_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sage2_nprofs = nindex
        sage2_longitudes = sage2_orig_longitudes[index]
        sage2_latitudes = sage2_orig_latitudes[index]
        sage2_time = sage2_orig_time[index]
        sage2_occ_num = sage2_orig_id[index] 
        sage2_occ_type = sage2_orig_occ_type[index]
      endif else sage2_nprofs = 0
    endif else sage2_nprofs = 0
    
    ;read SAGE3
    file=file_search(sage3_obs_path+'cat_sage3_v3.00.'+str_year,count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with SAGE3 data anyway
    if nfiles ne 0 then begin
      if file ne sage3_old_file then begin 
        restore, file[0]
        sage3_orig_date = date
        sage3_orig_id = id
        sage3_orig_latitudes = latitude
        sage3_orig_longitudes = longitude
        sage3_orig_occ_type = fix(1*(sctype eq 'r')-1*(sctype eq 's'))
        sage3_orig_time = time
        sage3_old_file=file
        undefine, date & undefine, id & undefine, latitude & undefine, longitude & undefine, sctype & undefine, time 
        print,'starting SAGE3 '+str_year
      endif
      index=where(sage3_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sage3_nprofs = nindex
        sage3_longitudes = sage3_orig_longitudes[index]
        sage3_latitudes = sage3_orig_latitudes[index]
        sage3_time = sage3_orig_time[index]
        sage3_occ_num = sage3_orig_id[index] 
        sage3_occ_type = sage3_orig_occ_type[index]
      endif else sage3_nprofs = 0
    endif else sage3_nprofs = 0
    
    ;read HALOE
    file=file_search(haloe_obs_path+'cat_haloe_v19.'+str_year,count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with HALOE data anyway
    if nfiles ne 0 then begin
      if file ne haloe_old_file then begin 
        restore, file[0]
        haloe_orig_date = date
        haloe_orig_id = id
        haloe_orig_latitudes = latitude
        haloe_orig_longitudes = longitude
        haloe_orig_occ_type = fix(1*(sctype eq 'r')-1*(sctype eq 's'))
        haloe_orig_time = time
        haloe_old_file=file
        undefine, date & undefine, id & undefine, latitude & undefine, longitude & undefine, sctype & undefine, time 
        print,'starting HALOE '+str_year
      endif
      index=where(haloe_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        haloe_nprofs = nindex
        haloe_longitudes = haloe_orig_longitudes[index]
        haloe_latitudes = haloe_orig_latitudes[index]
        haloe_time = haloe_orig_time[index]
        haloe_occ_num = haloe_orig_id[index] 
        haloe_occ_type = haloe_orig_occ_type[index]
      endif else haloe_nprofs = 0
    endif else haloe_nprofs = 0
    
    ;read MLS
    file=file_search(mls_obs_path+'O3\'+str_year+'\MLS-Aura_L2GP-O3_v??-??-c??_'+str_date_code+'.he5',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version
    if nfiles ne 0 then begin
;      data=readl2gp_std(file[0])
      data=read_mls_v3(file[0])
      mls_nprofs = data.ntimes
      mls_latitudes = data.latitude
      mls_longitudes = data.longitude
      mls_time = data.time
      mls_ltime = data.localsolartime
      mls_sza = data.solarzenithangle
      mls_oga = (data.orbitgeodeticangle)*10
      
;      ;read and apply time correction
;      hdfid = h5f_open(file[0])
;        groupid = h5g_open(hdfid, '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES')
;          attributeid = h5a_open_name(groupid,'TAI93At0zOfGranule')
;          tai93atozofgranule = (h5a_read(attributeid))[0]
;          h5a_close, attributeid
;        h5g_close,groupid
;      h5f_close, hdfid
;      mls_time = mls_time - tai93atozofgranule
    endif else mls_nprofs = 0
    
    ;read ACE-FTS
    file=file_search(ace_obs_path+'ACE-FTS_v3.0_'+str_year+'.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with ACE data anyway
    if nfiles ne 0 then begin
      if file ne ace_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'longitude'), ace_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'latitude'), ace_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), ace_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), ace_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'occultation_number'), ace_orig_occ_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'occ_type'), ace_orig_occ_type
        ncdf_close,ncid2
        ace_old_file=file
        ace_orig_nprofs=n_elements(ace_orig_longitudes)
        print,'starting ACE-FTS '+str_year
      endif
      index=where(ace_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        ace_nprofs = nindex
        ace_longitudes = ace_orig_longitudes[index]
        ace_latitudes = ace_orig_latitudes[index]
        ace_time = ace_orig_time[index]
        ace_occ_num = ace_orig_occ_num[index] 
        ace_occ_type = ace_orig_occ_type[index]
        ace_occ_type = fix(1*(ace_occ_type eq 114)-1*(ace_occ_type eq 115)) ;change r (114) to 1 and s (115) to -1
      endif else ace_nprofs = 0
    endif else ace_nprofs = 0
    
;    ;read HIRDLS
;    ;file=file_search(hirdls_obs_path+str_year+'\HIRDLS-Aura_L2_v??-??-??-c??_'+str_date_code+'.he5',count=nfiles)
;;    file=file_search(hirdls_obs_path+'\HIRDLS2ALL_v05-05-00-c03_'+str_date_code+'.he5',count=nfiles)
;    file=file_search(hirdls_obs_path+str_year+'\HIRDLS-Aura_L2_v06-??-??-c??_'+str_date_code+'.he5',count=nfiles)
;    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version
;    if nfiles ne 0 then begin
;      data=read_hirdls_v5(file[0],'O3') ;seems to work for v6 as well
;      hirdls_nprofs = data.ntimes
;      hirdls_latitudes = data.latitude
;      hirdls_longitudes = data.longitude
;      hirdls_daytime = data.secondsinday
;;      hirdls_ltime = data.localsolartime ;since original ltime exceeds 24h on many occasions, it's inferred from SID and lon instead
;      hirdls_ltime = utc2ltime(hirdls_daytime/3600.,hirdls_longitudes)
;      hirdls_profileid = data.profile_id
;      hirdls_sza = data.solarzenithangle
;    endif else hirdls_nprofs = 0
    
    ;read HIRDLS (from HIRDLS location file)
    file=file_search(hirdls_obs_path+'\HIRDLS_sat_hist_20050122-20080317_c20121019.nc',count=nfiles)
    if nfiles ne 0 then begin
      if file ne hirdls_old_file then begin
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), hirdls_date_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), hirdls_time_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), hirdls_lat_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), hirdls_lon_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'prof_num'), hirdls_id_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'local_time'), hirdls_ltime_all
          ncdf_varget,ncid2,ncdf_varid(ncid2,'instr_sza'), hirdls_sza_all          
        ncdf_close,ncid2
        hirdls_old_file = file[0]
      endif
      index=where(hirdls_date_all eq current_date,nindex)
      if nindex ne 0 then begin
        hirdls_nprofs = nindex
        hirdls_latitudes = hirdls_lat_all[index]
        hirdls_longitudes = hirdls_lon_all[index]
        hirdls_daytime = hirdls_time_all[index]
        hirdls_ltime = hirdls_ltime_all[index]
        hirdls_profileid = hirdls_id_all[index]
        hirdls_sza = hirdls_sza_all[index]
      endif else hirdls_nprofs = 0
    endif else hirdls_nprofs = 0
    
    ;read SABER
    file=file_search(saber_obs_path+'SABER'+strmid(str_date_code,0,4)+'.txt',count=nfiles)
    if nfiles ne 0 then begin
      if file ne saber_old_file then begin
        data=read_ascii(file[0],template=saber_tpl)
        saber_orig_orbit = data.orbit
        saber_orig_event = fix(data.event)
        saber_orig_date = long(strmid(data.date,0,4)+strmid(data.date,5,2)+strmid(data.date,8,2))
        saber_orig_daytime = long(strmid(data.time,0,2)*3600.+strmid(data.time,3,2)*60.+strmid(data.time,6,2))
        saber_orig_doy = data.doy
        saber_orig_lat = data.latitude
        saber_orig_lon = data.longitude
        saber_orig_sza = data.sza
        saber_orig_profiles = n_elements(data.latitude)
        saber_old_file=file[0]
      endif
      index=where(saber_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        saber_nprofs = nindex
        saber_longitudes = saber_orig_lon[index]
        saber_latitudes = saber_orig_lat[index]
        saber_daytime = saber_orig_daytime[index]
        
        ;fix seconds gt 86400 issue (SABER switches date but not seconds timer)
        tindex = where(saber_daytime ge 86460,ntindex) ;86400-86459 are left in to be filtered later, rest is corrected
        if ntindex ne 0 then saber_daytime[tindex]=saber_daytime[tindex]-86400
        
        saber_orbit = saber_orig_orbit[index]
        saber_event = saber_orig_event[index]
        saber_ltime = utc2ltime(saber_daytime/3600.,saber_longitudes)
        saber_sza = saber_orig_sza[index]
      endif else saber_nprofs = 0
    endif else saber_nprofs=0
    
;    ;read SOFIE
;    file=file_search(sofie_obs_path+'SOFIE_Level2_1.2_ic_????????.sav',count=nfiles)
;    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version
;    if nfiles ne 0 then begin
;      if sofie_old_file ne file[0] then restore, file[0]
;      index=where(sofie_orig_date eq current_date,nindex)
;      if nindex ne 0 then begin
;        sofie_nprofs = nindex
;        sofie_longitudes = sofie_orig_lon[index]
;        sofie_latitudes = sofie_orig_lat[index]
;        sofie_daytime = sofie_orig_sid[index]
;        sofie_orbit = sofie_orig_orbit[index]
;        sofie_event = sofie_orig_event[index]
;        sofie_occ = sofie_orig_mode[index]
;      endif else sofie_nprofs=0
;      sofie_old_file=file[0]      
;    endif else sofie_nprofs=0
    
    ;read SOFIE
    file_sr=file_search(sofie_obs_path+'sofie_eventLocation_V1.2_rise.txt',count=nfiles_sr)
    file_ss=file_search(sofie_obs_path+'sofie_eventLocation_V1.2_set.txt',count=nfiles_ss)
    if nfiles_sr gt 1 then file_sr=(file_sr[sort(file_sr)])[nfiles_sr-1] ;pick highest version
    if nfiles_ss gt 1 then file_ss=(file_ss[sort(file_ss)])[nfiles_ss-1] ;pick highest version
    if nfiles_sr*nfiles_ss ne 0 and sofie_sr_old_file ne file_sr[0] and sofie_ss_old_file ne file_ss[0] then begin
      data_sr=read_ascii(file_sr[0],template=sofie_tpl)
      data_ss=read_ascii(file_ss[0],template=sofie_tpl)
      sofie_orig_event = [data_sr.event,data_ss.event]
      sofie_orig_orbit = [data_sr.orbit,data_ss.orbit]
      sofie_orig_date = [data_sr.date,data_ss.date]
      sofie_orig_time = [data_sr.time,data_ss.time]
      sofie_orig_latitude = [data_sr.latitude,data_ss.latitude]
      sofie_orig_longitude = [data_sr.longitude,data_ss.longitude]
      sofie_orig_occtype = [replicate(1,n_elements(data_sr.event)),replicate(-1,n_elements(data_ss.event))]
      isort = sort(sofie_orig_event)
      sofie_orig_event = sofie_orig_event[isort]
      sofie_orig_orbit = sofie_orig_orbit[isort]
      sofie_orig_date = sofie_orig_date[isort]
      sofie_orig_time = sofie_orig_time[isort]
      sofie_orig_lat = sofie_orig_latitude[isort]
      sofie_orig_lon = sofie_orig_longitude[isort]
      sofie_orig_occtype = sofie_orig_occtype[isort]
      tmp_year = floor(sofie_orig_date/1000)
      tmp_doy = sofie_orig_date - tmp_year*1000
      sofie_orig_date = doy2date(tmp_doy,tmp_year)
      undefine, tmp_year & undefine, tmp_doy      
      sofie_orig_sid = long(sofie_orig_time*3600.)
      sofie_sr_old_file = file_sr[0]
      sofie_ss_old_file = file_ss[0]
    endif else sofie_nprofs=0
    index=where(sofie_orig_date eq current_date,nindex)
    if nindex ne 0 then begin
      sofie_nprofs = nindex
      sofie_longitudes = sofie_orig_lon[index]
      sofie_latitudes = sofie_orig_lat[index]
      sofie_daytime = sofie_orig_sid[index]
      sofie_orbit = sofie_orig_orbit[index]
      sofie_event = sofie_orig_event[index]
      sofie_occ = sofie_orig_occtype[index]
    endif else sofie_nprofs=0
    
        
    ;read CONCORDIASI PSC16
    file=file_search(concordiasi_obs_path+'PSC16_UCOZ.txt',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc16_old_file then begin
        data=read_ascii(file[0],template=concordiasi_tpl)
        caldat,data.julian,mm,dd,yyyy,h,m,s
        concordiasi_psc16_orig_date = long(strcompress(string(yyyy)+string(mm,format='(i2.2)')+string(dd,format='(i2.2)'),/r))
        concordiasi_psc16_orig_daytime = long(round(h*3600.+m*60.+s))
        concordiasi_psc16_orig_doy = date2doy(concordiasi_psc16_orig_date)
        concordiasi_psc16_orig_profiles = n_elements(data.latitude)
        concordiasi_psc16_orig_ident = indgen(concordiasi_psc16_orig_profiles)
        concordiasi_psc16_orig_lat = data.latitude
        concordiasi_psc16_orig_lon = data.longitude
        concordiasi_psc16_old_file=file[0]
      endif
      index=where(concordiasi_psc16_orig_date eq current_date and $
                  concordiasi_psc16_orig_lat ge -90 and concordiasi_psc16_orig_lat le 90 and $
                  concordiasi_psc16_orig_lon ge -180 and concordiasi_psc16_orig_lon le 180 and $
                  concordiasi_psc16_orig_daytime ge 0 and concordiasi_psc16_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc16_nprofs = nindex
        concordiasi_psc16_longitudes = concordiasi_psc16_orig_lon[index]
        concordiasi_psc16_latitudes = concordiasi_psc16_orig_lat[index]
        concordiasi_psc16_daytime = concordiasi_psc16_orig_daytime[index]
        concordiasi_psc16_date = concordiasi_psc16_orig_date[index]
        concordiasi_psc16_doy = concordiasi_psc16_orig_doy[index]
        concordiasi_psc16_ident = concordiasi_psc16_orig_ident[index]
        
        concordiasi_psc16_ltime = utc2ltime(concordiasi_psc16_daytime/3600.,concordiasi_psc16_longitudes)
      endif else concordiasi_psc16_nprofs = 0
    endif else concordiasi_psc16_nprofs=0
    
    ;read CONCORDIASI PSC17
    file=file_search(concordiasi_obs_path+'PSC17_UCOZ.txt',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc17_old_file then begin
        data=read_ascii(file[0],template=concordiasi_tpl)
        caldat,data.julian,mm,dd,yyyy,h,m,s
        concordiasi_psc17_orig_date = long(strcompress(string(yyyy)+string(mm,format='(i2.2)')+string(dd,format='(i2.2)'),/r))
        concordiasi_psc17_orig_daytime = long(round(h*3600.+m*60.+s))
        concordiasi_psc17_orig_doy = date2doy(concordiasi_psc17_orig_date)
        concordiasi_psc17_orig_profiles = n_elements(data.latitude)
        concordiasi_psc17_orig_ident = indgen(concordiasi_psc17_orig_profiles)
        concordiasi_psc17_orig_lat = data.latitude
        concordiasi_psc17_orig_lon = data.longitude
        concordiasi_psc17_old_file=file[0]
      endif
      index=where(concordiasi_psc17_orig_date eq current_date and $
                  concordiasi_psc17_orig_lat ge -90 and concordiasi_psc17_orig_lat le 90 and $
                  concordiasi_psc17_orig_lon ge -180 and concordiasi_psc17_orig_lon le 180 and $
                  concordiasi_psc17_orig_daytime ge 0 and concordiasi_psc17_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc17_nprofs = nindex
        concordiasi_psc17_longitudes = concordiasi_psc17_orig_lon[index]
        concordiasi_psc17_latitudes = concordiasi_psc17_orig_lat[index]
        concordiasi_psc17_daytime = concordiasi_psc17_orig_daytime[index]
        concordiasi_psc17_date = concordiasi_psc17_orig_date[index]
        concordiasi_psc17_doy = concordiasi_psc17_orig_doy[index]
        concordiasi_psc17_ident = concordiasi_psc17_orig_ident[index]
        
        concordiasi_psc17_ltime = utc2ltime(concordiasi_psc17_daytime/3600.,concordiasi_psc17_longitudes)
      endif else concordiasi_psc17_nprofs = 0
    endif else concordiasi_psc17_nprofs=0
    
    ;read CONCORDIASI PSC19
    file=file_search(concordiasi_obs_path+'PSC19_UCOZ.txt',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc19_old_file then begin
        data=read_ascii(file[0],template=concordiasi_tpl)
        caldat,data.julian,mm,dd,yyyy,h,m,s
        concordiasi_psc19_orig_date = long(strcompress(string(yyyy)+string(mm,format='(i2.2)')+string(dd,format='(i2.2)'),/r))
        concordiasi_psc19_orig_daytime = long(round(h*3600.+m*60.+s))
        concordiasi_psc19_orig_doy = date2doy(concordiasi_psc19_orig_date)
        concordiasi_psc19_orig_profiles = n_elements(data.latitude)
        concordiasi_psc19_orig_ident = indgen(concordiasi_psc19_orig_profiles)
        concordiasi_psc19_orig_lat = data.latitude
        concordiasi_psc19_orig_lon = data.longitude
        concordiasi_psc19_old_file=file[0]
      endif
      index=where(concordiasi_psc19_orig_date eq current_date and $
                  concordiasi_psc19_orig_lat ge -90 and concordiasi_psc19_orig_lat le 90 and $
                  concordiasi_psc19_orig_lon ge -180 and concordiasi_psc19_orig_lon le 180 and $
                  concordiasi_psc19_orig_daytime ge 0 and concordiasi_psc19_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc19_nprofs = nindex
        concordiasi_psc19_longitudes = concordiasi_psc19_orig_lon[index]
        concordiasi_psc19_latitudes = concordiasi_psc19_orig_lat[index]
        concordiasi_psc19_daytime = concordiasi_psc19_orig_daytime[index]
        concordiasi_psc19_date = concordiasi_psc19_orig_date[index]
        concordiasi_psc19_doy = concordiasi_psc19_orig_doy[index]
        concordiasi_psc19_ident = concordiasi_psc19_orig_ident[index]
        
        concordiasi_psc19_ltime = utc2ltime(concordiasi_psc19_daytime/3600.,concordiasi_psc19_longitudes)
      endif else concordiasi_psc19_nprofs = 0
    endif else concordiasi_psc19_nprofs=0
    
    ;read MIPAS
    file=file_search(mipas_obs_path+'\mipas_locations_????????.sav',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version
    if nfiles ne 0 then begin
      if mipas_old_file ne file[0] then begin
        restore, file[0]
        mipas_date = long(mipas_date)
      endif
      index=where(mipas_date eq current_date,nindex)
      if nindex ne 0 then begin
        mipas_nprofs = nindex
        mipas_jul = mipas_jday[index]
        mipas_lat = mipas_latitude[index]
        mipas_lon = mipas_longitude[index]
        mipas_sid = mipas_time[index]
        mipas_localtime = mipas_ltime[index]
        mipas_sza_tmp = mipas_sza[index]
        mipas_identifier = (mipas_jul-date2julday(20020101))*100000+round(mipas_sid) 
      endif else mipas_nprofs=0
      mipas_old_file=file[0]
    endif else mipas_nprofs = 0
    
    ;read JUELICH (PEX2005, REC2010, PEX2010)
    file=file_search(juelich_obs_path+'out.dat',count=nfiles)
    if nfiles ne 0 then begin
      if file ne juelich_old_file then begin
        data=read_ascii(file[0],template=juelich_tpl)
        juelich_orig_date = data.date
        juelich_orig_daytime = 60.*data.time
        juelich_orig_doy = date2doy(juelich_orig_date)
        juelich_orig_profiles = n_elements(data.date)
        
        juelich_orig_instr = strmid(data.tag,0,3)
        tmp_year = fix(strmid(data.tag,3,4))
        tmp_flightnr = fix(strmid(data.tag,8,2))
        tmp_n_pex2005 = total(strmid(data.tag,0,7) eq 'PEX2005')
        tmp_n_rec2010 = total(strmid(data.tag,0,7) eq 'REC2010')
        tmp_n_pex2010 = total(strmid(data.tag,0,7) eq 'PEX2010')
        tmp_obsnr = 1+[indgen(tmp_n_pex2005),indgen(tmp_n_rec2010),indgen(tmp_n_pex2010)]
        if n_elements(tmp_obsnr) ne juelich_orig_profiles then stop
                
        juelich_orig_ident = tmp_year*100000+tmp_flightnr*1000+tmp_obsnr
        juelich_orig_lat = data.lat
        juelich_orig_lon = data.lon
        juelich_old_file=file[0]
        print,'starting JUELICH (PEX2005, REC2010, PEX2010) '+str_year
      endif
      index=where(juelich_orig_date eq current_date and $
                  juelich_orig_lat ge -90 and juelich_orig_lat le 90 and $
                  juelich_orig_lon ge 0 and juelich_orig_lon le 360 and $
                  juelich_orig_daytime ge 0 and juelich_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        juelich_nprofs = nindex
        juelich_longitudes = juelich_orig_lon[index]
        juelich_latitudes = juelich_orig_lat[index]
        juelich_daytime = juelich_orig_daytime[index]
        juelich_date = juelich_orig_date[index]
        juelich_doy = juelich_orig_doy[index]
        juelich_ident = juelich_orig_ident[index]
        juelich_instr = juelich_orig_instr[index]
        
        juelich_ltime = utc2ltime(juelich_daytime/3600.,juelich_longitudes)
      endif else juelich_nprofs = 0
    endif else juelich_nprofs=0
    
    ;read SAGE1
    file=file_search(sage1_obs_path+'SAGE1_catalogfile_c20120917.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with SAGE1 data anyway
    if nfiles ne 0 then begin
      if file ne sage1_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), sage1_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), sage1_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), sage1_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), sage1_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'occ_type'), sage1_orig_occ_type
        ncdf_close,ncid2
        sage1_orig_id = lindgen(n_elements(sage1_orig_longitudes))
        sage1_old_file=file
      endif
      index=where(sage1_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sage1_nprofs = nindex
        sage1_longitudes = sage1_orig_longitudes[index]
        sage1_latitudes = sage1_orig_latitudes[index]
        sage1_time = sage1_orig_time[index]
        sage1_ident = sage1_orig_id[index] 
        sage1_occ_type = sage1_orig_occ_type[index]
      endif else sage1_nprofs = 0
    endif else sage1_nprofs = 0
    
    ;read UARS MLS
    file=file_search(uars_mls_obs_path+'UARS_MLS_catalogfile_c20120918.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with UARS MLS data anyway
    if nfiles ne 0 then begin
      if file ne uars_mls_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), uars_mls_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), uars_mls_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), uars_mls_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), uars_mls_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'sza'), uars_mls_orig_sza
        ncdf_close,ncid2
        uars_mls_orig_id = lindgen(n_elements(uars_mls_orig_longitudes))
        uars_mls_old_file=file
      endif
      index=where(uars_mls_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        uars_mls_nprofs = nindex
        uars_mls_longitudes = uars_mls_orig_longitudes[index]
        uars_mls_latitudes = uars_mls_orig_latitudes[index]
        uars_mls_time = uars_mls_orig_time[index]
        uars_mls_ident = uars_mls_orig_id[index] 
        uars_mls_sza = uars_mls_orig_sza[index]
      endif else uars_mls_nprofs = 0
    endif else uars_mls_nprofs = 0
    
    ;read LIMS
    file=file_search(lims_obs_path+'LIMS_catalogfile_c20120918.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with UARS MLS data anyway
    if nfiles ne 0 then begin
      if file ne lims_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), lims_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), lims_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), lims_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), lims_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'sza'), lims_orig_sza
        ncdf_close,ncid2
        lims_orig_id = lindgen(n_elements(lims_orig_longitudes))
        lims_old_file=file
      endif
      index=where(lims_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        lims_nprofs = nindex
        lims_longitudes = lims_orig_longitudes[index]
        lims_latitudes = lims_orig_latitudes[index]
        lims_time = lims_orig_time[index]
        lims_ident = lims_orig_id[index] 
        lims_sza = lims_orig_sza[index]
      endif else lims_nprofs = 0
    endif else lims_nprofs = 0
    
    ;read SMILES
    file=file_search(smiles_obs_path+'obs_info_20120405.txt',count=nfiles)
    if nfiles ne 0 then begin
      if file ne smiles_old_file then begin
        data=read_ascii(file[0],template=smiles_tpl)
        smiles_orig_obs_nr = fix(data.obs_nr)
        smiles_orig_date = long(strmid(data.date_time,0,4)+strmid(data.date_time,5,2)+strmid(data.date_time,8,2))
        smiles_orig_sid = long(strmid(data.date_time,11,2))*3600.+long(strmid(data.date_time,14,2))*60.+long(strmid(data.date_time,17,2))
        smiles_orig_ltime = data.ltime
        smiles_orig_lat = data.lat
        smiles_orig_lon = data.lon
        smiles_orig_sza = data.sza
        smiles_orig_ascdsc = data.asc_dsc
        smiles_orig_profiles = n_elements(data.lat)
        smiles_old_file=file[0]
      endif
      index=where(smiles_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        smiles_nprofs = nindex
        smiles_longitudes = smiles_orig_lon[index]
        smiles_latitudes = smiles_orig_lat[index]
        smiles_time = smiles_orig_sid[index]
        
        smiles_ident = smiles_orig_obs_nr[index]
        smiles_ltime = smiles_orig_ltime[index]
        smiles_sza = smiles_orig_sza[index]
        smiles_node = smiles_orig_ascdsc[index]
      endif else smiles_nprofs = 0
    endif else smiles_nprofs=0
    
    ;read CONCORDIASI PSC14
    file=file_search(concordiasi_obs_path+'french data\bb10V14N42.dat',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc19_old_file then begin
        data=read_ascii(file[0],template=concordiasi2_tpl)
        ;'decimal_day','year','month','day','hour_UTC','minute_UTC','second_UTC','lon','lat','SZA'
        concordiasi_psc14_orig_date = data.year*10000l+data.month*100l+data.day*1l
        concordiasi_psc14_orig_daytime = long(round(data.hour_UTC*3600l+data.minute_UTC*60l+data.second_UTC*1l))
        concordiasi_psc14_orig_doy = date2doy(concordiasi_psc14_orig_date)
        concordiasi_psc14_orig_profiles = n_elements(data.lat)
        concordiasi_psc14_orig_ident = indgen(concordiasi_psc14_orig_profiles)
        concordiasi_psc14_orig_lat = data.lat
        concordiasi_psc14_orig_lon = data.lon
        concordiasi_psc14_orig_sza = data.SZA
        concordiasi_psc14_old_file=file[0]
      endif
      index=where(concordiasi_psc14_orig_date eq current_date and $
                  concordiasi_psc14_orig_lat ge -90 and concordiasi_psc14_orig_lat le 90 and $
                  concordiasi_psc14_orig_lon ge -180 and concordiasi_psc14_orig_lon le 180 and $
                  concordiasi_psc14_orig_daytime ge 0 and concordiasi_psc14_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc14_nprofs = nindex
        concordiasi_psc14_longitudes = concordiasi_psc14_orig_lon[index]
        concordiasi_psc14_latitudes = concordiasi_psc14_orig_lat[index]
        concordiasi_psc14_daytime = concordiasi_psc14_orig_daytime[index]
        concordiasi_psc14_date = concordiasi_psc14_orig_date[index]
        concordiasi_psc14_doy = concordiasi_psc14_orig_doy[index]
        concordiasi_psc14_ident = concordiasi_psc14_orig_ident[index]
        concordiasi_psc14_sza = concordiasi_psc14_orig_sza[index]
        
        concordiasi_psc14_ltime = utc2ltime(concordiasi_psc14_daytime/3600.,concordiasi_psc14_longitudes)
      endif else concordiasi_psc14_nprofs = 0
    endif else concordiasi_psc14_nprofs=0
    
    ;read CONCORDIASI PSC15
    file=file_search(concordiasi_obs_path+'french data\bb10V15N32.dat',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc19_old_file then begin
        data=read_ascii(file[0],template=concordiasi2_tpl)
        ;'decimal_day','year','month','day','hour_UTC','minute_UTC','second_UTC','lon','lat','SZA'
        concordiasi_psc15_orig_date = data.year*10000l+data.month*100l+data.day*1l
        concordiasi_psc15_orig_daytime = long(round(data.hour_UTC*3600l+data.minute_UTC*60l+data.second_UTC*1l))
        concordiasi_psc15_orig_doy = date2doy(concordiasi_psc15_orig_date)
        concordiasi_psc15_orig_profiles = n_elements(data.lat)
        concordiasi_psc15_orig_ident = indgen(concordiasi_psc15_orig_profiles)
        concordiasi_psc15_orig_lat = data.lat
        concordiasi_psc15_orig_lon = data.lon
        concordiasi_psc15_orig_sza = data.SZA
        concordiasi_psc15_old_file=file[0]
      endif
      index=where(concordiasi_psc15_orig_date eq current_date and $
                  concordiasi_psc15_orig_lat ge -90 and concordiasi_psc15_orig_lat le 90 and $
                  concordiasi_psc15_orig_lon ge -180 and concordiasi_psc15_orig_lon le 180 and $
                  concordiasi_psc15_orig_daytime ge 0 and concordiasi_psc15_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc15_nprofs = nindex
        concordiasi_psc15_longitudes = concordiasi_psc15_orig_lon[index]
        concordiasi_psc15_latitudes = concordiasi_psc15_orig_lat[index]
        concordiasi_psc15_daytime = concordiasi_psc15_orig_daytime[index]
        concordiasi_psc15_date = concordiasi_psc15_orig_date[index]
        concordiasi_psc15_doy = concordiasi_psc15_orig_doy[index]
        concordiasi_psc15_ident = concordiasi_psc15_orig_ident[index]
        concordiasi_psc15_sza = concordiasi_psc15_orig_sza[index]
        
        concordiasi_psc15_ltime = utc2ltime(concordiasi_psc15_daytime/3600.,concordiasi_psc15_longitudes)
      endif else concordiasi_psc15_nprofs = 0
    endif else concordiasi_psc15_nprofs=0
    
    ;read CONCORDIASI PSC18
    file=file_search(concordiasi_obs_path+'french data\bb10V18N43.dat',count=nfiles)
    if nfiles ne 0 then begin
      if file ne concordiasi_psc19_old_file then begin
        data=read_ascii(file[0],template=concordiasi2_tpl)
        ;'decimal_day','year','month','day','hour_UTC','minute_UTC','second_UTC','lon','lat','SZA'
        concordiasi_psc18_orig_date = data.year*10000l+data.month*100l+data.day*1l
        concordiasi_psc18_orig_daytime = long(round(data.hour_UTC*3600l+data.minute_UTC*60l+data.second_UTC*1l))
        concordiasi_psc18_orig_doy = date2doy(concordiasi_psc18_orig_date)
        concordiasi_psc18_orig_profiles = n_elements(data.lat)
        concordiasi_psc18_orig_ident = indgen(concordiasi_psc18_orig_profiles)
        concordiasi_psc18_orig_lat = data.lat
        concordiasi_psc18_orig_lon = data.lon
        concordiasi_psc18_orig_sza = data.SZA
        concordiasi_psc18_old_file=file[0]
      endif
      index=where(concordiasi_psc18_orig_date eq current_date and $
                  concordiasi_psc18_orig_lat ge -90 and concordiasi_psc18_orig_lat le 90 and $
                  concordiasi_psc18_orig_lon ge -180 and concordiasi_psc18_orig_lon le 180 and $
                  concordiasi_psc18_orig_daytime ge 0 and concordiasi_psc18_orig_daytime le 86400, nindex)
      if nindex ne 0 then begin
        concordiasi_psc18_nprofs = nindex
        concordiasi_psc18_longitudes = concordiasi_psc18_orig_lon[index]
        concordiasi_psc18_latitudes = concordiasi_psc18_orig_lat[index]
        concordiasi_psc18_daytime = concordiasi_psc18_orig_daytime[index]
        concordiasi_psc18_date = concordiasi_psc18_orig_date[index]
        concordiasi_psc18_doy = concordiasi_psc18_orig_doy[index]
        concordiasi_psc18_ident = concordiasi_psc18_orig_ident[index]
        concordiasi_psc18_sza = concordiasi_psc18_orig_sza[index]
        
        concordiasi_psc18_ltime = utc2ltime(concordiasi_psc18_daytime/3600.,concordiasi_psc18_longitudes)
      endif else concordiasi_psc18_nprofs = 0
    endif else concordiasi_psc18_nprofs=0
    
    ;read BUV/SBUV/SBUV2
    file=file_search(sbuv_obs_path+'BUV_SBUV_SBUV2_'+str_year+'_catalog_c????????.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick latest processing of catalog files
    if nfiles ne 0 then begin
      if file ne sbuv_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), sbuv_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), sbuv_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), sbuv_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), sbuv_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'orbit_num'), sbuv_orig_orbit_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'sat_num'), sbuv_orig_sat_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'local_time'), sbuv_orig_ltime
          ncdf_varget,ncid2,ncdf_varid(ncid2,'instr_sza'), sbuv_orig_sza
          ncdf_varget,ncid2,ncdf_varid(ncid2,'julian'), sbuv_orig_julian
        ncdf_close,ncid2
        sbuv_old_file=file
        sbuv_orig_nprofs=n_elements(sbuv_orig_longitudes)
        print,'starting BUV/SBUV/SBUV2 '+str_year
      endif
      index=where(sbuv_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sbuv_nprofs = nindex
        sbuv_time = sbuv_orig_time[index]
        sbuv_latitudes = sbuv_orig_latitudes[index]
        sbuv_longitudes = sbuv_orig_longitudes[index]
        sbuv_orbit_num = sbuv_orig_orbit_num[index]
        sbuv_sat_num = sbuv_orig_sat_num[index]
        sbuv_ltime = sbuv_orig_ltime[index]
        sbuv_sza = sbuv_orig_sza[index]
        sbuv_julian = sbuv_orig_julian[index]
      endif else sbuv_nprofs = 0
    endif else sbuv_nprofs = 0
    if sbuv_nprofs eq 0 then sbuv_sat_num=0 ; needed for log file
    
    ;read SMR
    file=file_search(smr_obs_path+'smr_geolocation_c????????.sav',count=nfiles)
    if nfiles ne 0 then begin
      if nfiles gt 1 then file = (reverse(file[sort(file)]))[0] else file = file[0]
      if file ne smr_old_file then begin
        restore, file
        smr_orig_date = smr_date
        smr_orig_sid = smr_sid
        smr_orig_lat = smr_lat
        smr_orig_lon = smr_lon
        smr_orig_profiles = n_elements(smr_orig_lat)
        smr_old_file=file
        undefine, smr_date & undefine, smr_lat & undefine, smr_lon & undefine, smr_sid
      endif
      index=where(smr_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        smr_nprofs = nindex
        smr_longitudes = smr_orig_lon[index]
        smr_latitudes = smr_orig_lat[index]
        smr_time = smr_orig_sid[index]
        
        smr_ident = (date2julday(current_date)-date2julday(20010619))*100000l + smr_time 
        smr_ltime = utc2ltime(smr_time/86400.,smr_longitudes)
      endif else smr_nprofs = 0
    endif else smr_nprofs=0
    
    ;read SME AG
    file=file_search(sme_obs_path+'SME_AG_profilelist_c.nc',count=nfiles)
    if nfiles ne 0 then begin
      if file ne sme_ag_old_file then begin
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), sme_ag_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), sme_ag_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), sme_ag_orig_lat
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), sme_ag_orig_lon
          ncdf_varget,ncid2,ncdf_varid(ncid2,'orbit_num'), sme_ag_orig_orbit_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'instr_sza'), sme_ag_orig_sza
        ncdf_close,ncid2
        sme_ag_orig_profiles = n_elements(sme_ag_orig_lat)
        sme_ag_old_file=file
      endif
      index=where(sme_ag_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sme_ag_nprofs = nindex
        sme_ag_longitudes = sme_ag_orig_lon[index]
        sme_ag_latitudes = sme_ag_orig_lat[index]
        sme_ag_time = sme_ag_orig_time[index]
        sme_ag_orbit_num = sme_ag_orig_orbit_num[index]
        sme_ag_sza = sme_ag_orig_sza[index]
        
        sme_ag_ident = (date2julday(current_date)-date2julday(19811215))*100000l + sme_ag_time 
        sme_ag_ltime = utc2ltime(sme_ag_time/86400.,sme_ag_longitudes)
      endif else sme_ag_nprofs = 0
    endif else sme_ag_nprofs=0
    
    ;read SME UV
    file=file_search(sme_obs_path+'SME_UV_profilelist_c.nc',count=nfiles)
    if nfiles ne 0 then begin
      if file ne sme_uv_old_file then begin
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), sme_uv_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), sme_uv_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lat'), sme_uv_orig_lat
          ncdf_varget,ncid2,ncdf_varid(ncid2,'lon'), sme_uv_orig_lon
          ncdf_varget,ncid2,ncdf_varid(ncid2,'orbit_num'), sme_uv_orig_orbit_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'instr_sza'), sme_uv_orig_sza
        ncdf_close,ncid2
        sme_uv_orig_profiles = n_elements(sme_uv_orig_lat)
        sme_uv_old_file=file
      endif
      index=where(sme_uv_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        sme_uv_nprofs = nindex
        sme_uv_longitudes = sme_uv_orig_lon[index]
        sme_uv_latitudes = sme_uv_orig_lat[index]
        sme_uv_time = sme_uv_orig_time[index]
        sme_uv_orbit_num = sme_uv_orig_orbit_num[index]
        sme_uv_sza = sme_uv_orig_sza[index]
        
        sme_uv_ident = (date2julday(current_date)-date2julday(19811215))*100000l + sme_uv_time 
        sme_uv_ltime = utc2ltime(sme_uv_time/86400.,sme_uv_longitudes)
      endif else sme_uv_nprofs = 0
    endif else sme_uv_nprofs=0
    
    ;read GOMOS
    file=file_search(gomos_obs_path+'gomos_geolocations_c20130226.sav',count=nfiles)
    if nfiles ne 0 then begin
      if file ne gomos_old_file then begin
        restore, file[0]
        gomos_orig_profiles = n_elements(gomos_date)
        gomos_old_file=file
      endif
      index=where(gomos_date eq current_date,nindex)
      if nindex ne 0 then begin
        gomos_nprofs = nindex
        gomos_latitudes = [gomos_latitude_bottom[index],gomos_latitude_top[index]]
        gomos_longitudes = [gomos_longitude_bottom[index],gomos_longitude_top[index]]
        gomos_sid = [gomos_sid_bottom[index],gomos_sid_top[index]]
        gomos_sza = [gomos_sza_bottom[index],gomos_sza_top[index]]
                
        gomos_ident = (date2julday(current_date)-date2julday(20020101))*100000l + gomos_sid 
        gomos_ltime = utc2ltime(gomos_sid/86400.,gomos_longitudes)
      endif else gomos_nprofs = 0
    endif else gomos_nprofs=0
           
    ;concatenate
    instrument=0 & latitudes=0 & longitudes=0 & daytime=0 & ltime=0
    orbit=0l & event=0 & occ_type=0 & sza=0 & julian=0
    nprofs = long(mls_nprofs*1l + ace_nprofs + hirdls_nprofs + saber_nprofs + sofie_nprofs + $
             concordiasi_psc16_nprofs*4 + concordiasi_psc17_nprofs*4 + concordiasi_psc19_nprofs*4 + $
             mipas_nprofs + poam2_nprofs + poam3_nprofs + sage2_nprofs + sage3_nprofs + haloe_nprofs + $
             juelich_nprofs + sage1_nprofs + uars_mls_nprofs + smiles_nprofs + lims_nprofs + $
             concordiasi_psc14_nprofs*4 + concordiasi_psc15_nprofs*4 + concordiasi_psc18_nprofs*4 + $
             sbuv_nprofs + smr_nprofs + sme_ag_nprofs + sme_uv_nprofs + gomos_nprofs*2)
    if nprofs eq 0 then goto, nodata
    if mls_nprofs ne 0 then begin
      instrument=[instrument,replicate(1,mls_nprofs)]
      latitudes=[latitudes,mls_latitudes]
      longitudes=[longitudes,mls_longitudes]
      mls_sid=(mls_time-seconds_since_19930101)
      daytime=[daytime,mls_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,mls_ltime]
      orbit=[orbit,replicate(-1l,mls_nprofs)]
      event=[event,mls_oga]
      occ_type=[occ_type,replicate(0,mls_nprofs)]
      sza=[sza,mls_sza]
      mls_h=floor(mls_sid/3600.)
      mls_min=floor((mls_sid/3600.-mls_h)*60.)
      mls_s=round(mls_sid-mls_h*3600.-mls_min*60.)
      julian=[julian,julday(month,day,year,mls_h,mls_min,mls_s)-julian_reference]
    endif
    if ace_nprofs ne 0 then begin
      instrument=[instrument,replicate(2,ace_nprofs)]
      latitudes=[latitudes,ace_latitudes]
      longitudes=[longitudes,ace_longitudes]
      ace_sid=round(ace_time*3600)
      daytime=[daytime,ace_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(ace_time,ace_longitudes)]
      orbit=[orbit,replicate(-1l,ace_nprofs)]
      event=[event,ace_occ_num]
      occ_type=[occ_type,ace_occ_type]
      sza=[sza,replicate(-1,ace_nprofs)]
      ace_h=floor(ace_sid/3600.)
      ace_min=floor((ace_sid/3600.-ace_h)*60.)
      ace_s=round(ace_sid-ace_h*3600.-ace_min*60.)
      julian=[julian,julday(month,day,year,ace_h,ace_min,ace_s)-julian_reference]
    endif
    if hirdls_nprofs ne 0 then begin
      instrument=[instrument,replicate(3,hirdls_nprofs)]
      latitudes=[latitudes,hirdls_latitudes]
      longitudes=[longitudes,hirdls_longitudes]
      hirdls_sid=round(hirdls_daytime)
      daytime=[daytime,hirdls_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,hirdls_ltime]
      orbit=[orbit,replicate(-1l,hirdls_nprofs)]
      event=[event,hirdls_profileid]
      occ_type=[occ_type,replicate(0,hirdls_nprofs)]
      sza=[sza,hirdls_sza]
      hirdls_h=floor(hirdls_sid/3600.)
      hirdls_min=floor((hirdls_sid/3600.-hirdls_h)*60.)
      hirdls_s=round(hirdls_sid-hirdls_h*3600.-hirdls_min*60.)
      julian=[julian,julday(month,day,year,hirdls_h,hirdls_min,hirdls_s)-julian_reference]
    endif
    if saber_nprofs ne 0 then begin
      instrument=[instrument,replicate(4,saber_nprofs)]
      latitudes=[latitudes,saber_latitudes]
      longitudes=[longitudes,saber_longitudes]
      saber_sid=saber_daytime
      daytime=[daytime,saber_sid]
      ltime=[ltime,saber_ltime]
      orbit=[orbit,long(saber_orbit)]
      event=[event,saber_event]
      occ_type=[occ_type,replicate(0,saber_nprofs)]
      sza=[sza,saber_sza]
      saber_h=floor(saber_sid/3600.)
      saber_min=floor((saber_sid/3600.-saber_h)*60.)
      saber_s=round(saber_sid-saber_h*3600.-saber_min*60.)
      julian=[julian,julday(month,day,year,saber_h,saber_min,saber_s)-julian_reference]
    endif
    if sofie_nprofs ne 0 then begin
      instrument=[instrument,replicate(5,sofie_nprofs)]
      latitudes=[latitudes,sofie_latitudes]
      longitudes=[longitudes,sofie_longitudes]
      sofie_sid=sofie_daytime
      daytime=[daytime,sofie_sid]
      ltime=[ltime,utc2ltime(sofie_daytime/3600.,sofie_longitudes)]
      orbit=[orbit,long(sofie_orbit)]
      event=[event,sofie_event]
      occ_type=[occ_type,sofie_occ]
      sza=[sza,replicate(-1,sofie_nprofs)]
      sofie_h=floor(sofie_sid/3600.)
      sofie_min=floor((sofie_sid/3600.-sofie_h)*60.)
      sofie_s=round(sofie_sid-sofie_h*3600.-sofie_min*60.)
      julian=[julian,julday(month,day,year,sofie_h,sofie_min,sofie_s)-julian_reference]
    endif
    if concordiasi_psc16_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc16_longitudes = (concordiasi_psc16_longitudes lt 0)*(360+concordiasi_psc16_longitudes) + (concordiasi_psc16_longitudes ge 0)*concordiasi_psc16_longitudes
      concordiasi_psc16_new_nprofs=0 & concordiasi_psc16_new_lat=0 & concordiasi_psc16_new_lon=0 & concordiasi_psc16_new_sid=0 & concordiasi_psc16_new_ltime=0 & concordiasi_psc16_new_ident=0 & concordiasi_psc16_new_jul=0
      for iprof=0,concordiasi_psc16_nprofs-1 do begin
        west_lon = (concordiasi_psc16_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc16_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc16_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc16_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlats eq 0 then stop
        concordiasi_psc16_new_nprofs+=4
        concordiasi_psc16_new_lat=[concordiasi_psc16_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc16_new_lon=[concordiasi_psc16_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc16_new_sid=[concordiasi_psc16_new_sid,replicate(concordiasi_psc16_daytime[iprof],4)]
        concordiasi_psc16_new_ltime=[concordiasi_psc16_new_ltime,replicate(concordiasi_psc16_ltime[iprof],4)]
        concordiasi_psc16_new_ident=[concordiasi_psc16_new_ident,16000000+concordiasi_psc16_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc16_h=floor(concordiasi_psc16_daytime[iprof]/3600.)
        concordiasi_psc16_min=floor((concordiasi_psc16_daytime[iprof]/3600.-concordiasi_psc16_h)*60.)
        concordiasi_psc16_s=round(concordiasi_psc16_daytime[iprof]-concordiasi_psc16_h*3600.-concordiasi_psc16_min*60.)
        concordiasi_psc16_new_jul=[concordiasi_psc16_new_jul,replicate((julday(month,day,year,concordiasi_psc16_h,concordiasi_psc16_min,concordiasi_psc16_s)-julian_reference),4)]
      endfor
      concordiasi_psc16_new_lat=concordiasi_psc16_new_lat[1:*]
      concordiasi_psc16_new_lon=concordiasi_psc16_new_lon[1:*] & concordiasi_psc16_new_sid=concordiasi_psc16_new_sid[1:*]
      concordiasi_psc16_new_ltime=concordiasi_psc16_new_ltime[1:*] & concordiasi_psc16_new_ident=concordiasi_psc16_new_ident[1:*]
      concordiasi_psc16_new_jul=concordiasi_psc16_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(6,concordiasi_psc16_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc16_new_lat]
      longitudes=[longitudes,concordiasi_psc16_new_lon]
      daytime=[daytime,concordiasi_psc16_new_sid]
      ltime=[ltime,concordiasi_psc16_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc16_new_nprofs)]
      event=[event,concordiasi_psc16_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc16_new_nprofs)]
      sza=[sza,replicate(-1,concordiasi_psc16_new_nprofs)]
      julian=[julian,concordiasi_psc16_new_jul]
    endif
    if concordiasi_psc17_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc17_longitudes = (concordiasi_psc17_longitudes lt 0)*(360+concordiasi_psc17_longitudes) + (concordiasi_psc17_longitudes ge 0)*concordiasi_psc17_longitudes
      concordiasi_psc17_new_nprofs=0 & concordiasi_psc17_new_lat=0 & concordiasi_psc17_new_lon=0 & concordiasi_psc17_new_sid=0 & concordiasi_psc17_new_ltime=0 & concordiasi_psc17_new_ident=0 & concordiasi_psc17_new_jul=0
      for iprof=0,concordiasi_psc17_nprofs-1 do begin
        west_lon = (concordiasi_psc17_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc17_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc17_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc17_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlats eq 0 then stop
        concordiasi_psc17_new_nprofs+=4
        concordiasi_psc17_new_lat=[concordiasi_psc17_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc17_new_lon=[concordiasi_psc17_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc17_new_sid=[concordiasi_psc17_new_sid,replicate(concordiasi_psc17_daytime[iprof],4)]
        concordiasi_psc17_new_ltime=[concordiasi_psc17_new_ltime,replicate(concordiasi_psc17_ltime[iprof],4)]
        concordiasi_psc17_new_ident=[concordiasi_psc17_new_ident,17000000+concordiasi_psc17_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc17_h=floor(concordiasi_psc17_daytime[iprof]/3600.)
        concordiasi_psc17_min=floor((concordiasi_psc17_daytime[iprof]/3600.-concordiasi_psc17_h)*60.)
        concordiasi_psc17_s=round(concordiasi_psc17_daytime[iprof]-concordiasi_psc17_h*3600.-concordiasi_psc17_min*60.)
        concordiasi_psc17_new_jul=[concordiasi_psc17_new_jul,replicate((julday(month,day,year,concordiasi_psc17_h,concordiasi_psc17_min,concordiasi_psc17_s)-julian_reference),4)]
      endfor
      concordiasi_psc17_new_lat=concordiasi_psc17_new_lat[1:*]
      concordiasi_psc17_new_lon=concordiasi_psc17_new_lon[1:*] & concordiasi_psc17_new_sid=concordiasi_psc17_new_sid[1:*]
      concordiasi_psc17_new_ltime=concordiasi_psc17_new_ltime[1:*] & concordiasi_psc17_new_ident=concordiasi_psc17_new_ident[1:*]
      concordiasi_psc17_new_jul=concordiasi_psc17_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(7,concordiasi_psc17_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc17_new_lat]
      longitudes=[longitudes,concordiasi_psc17_new_lon]
      daytime=[daytime,concordiasi_psc17_new_sid]
      ltime=[ltime,concordiasi_psc17_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc17_new_nprofs)]
      event=[event,concordiasi_psc17_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc17_new_nprofs)]
      sza=[sza,replicate(-1,concordiasi_psc17_new_nprofs)]
      julian=[julian,concordiasi_psc17_new_jul]
    endif
    if concordiasi_psc19_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc19_longitudes = (concordiasi_psc19_longitudes lt 0)*(360+concordiasi_psc19_longitudes) + (concordiasi_psc19_longitudes ge 0)*concordiasi_psc19_longitudes
      concordiasi_psc19_new_nprofs=0 & concordiasi_psc19_new_lat=0 & concordiasi_psc19_new_lon=0 & concordiasi_psc19_new_sid=0 & concordiasi_psc19_new_ltime=0 & concordiasi_psc19_new_ident=0 & concordiasi_psc19_new_jul=0
      for iprof=0,concordiasi_psc19_nprofs-1 do begin
        west_lon = (concordiasi_psc19_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc19_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc19_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc19_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlons eq 0 then stop
        concordiasi_psc19_new_nprofs+=4
        concordiasi_psc19_new_lat=[concordiasi_psc19_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc19_new_lon=[concordiasi_psc19_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc19_new_sid=[concordiasi_psc19_new_sid,replicate(concordiasi_psc19_daytime[iprof],4)]
        concordiasi_psc19_new_ltime=[concordiasi_psc19_new_ltime,replicate(concordiasi_psc19_ltime[iprof],4)]
        concordiasi_psc19_new_ident=[concordiasi_psc19_new_ident,19000000+concordiasi_psc19_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc19_h=floor(concordiasi_psc19_daytime[iprof]/3600.)
        concordiasi_psc19_min=floor((concordiasi_psc19_daytime[iprof]/3600.-concordiasi_psc19_h)*60.)
        concordiasi_psc19_s=round(concordiasi_psc19_daytime[iprof]-concordiasi_psc19_h*3600.-concordiasi_psc19_min*60.)
        concordiasi_psc19_new_jul=[concordiasi_psc19_new_jul,replicate((julday(month,day,year,concordiasi_psc19_h,concordiasi_psc19_min,concordiasi_psc19_s)-julian_reference),4)]
      endfor
      concordiasi_psc19_new_lat=concordiasi_psc19_new_lat[1:*]
      concordiasi_psc19_new_lon=concordiasi_psc19_new_lon[1:*] & concordiasi_psc19_new_sid=concordiasi_psc19_new_sid[1:*]
      concordiasi_psc19_new_ltime=concordiasi_psc19_new_ltime[1:*] & concordiasi_psc19_new_ident=concordiasi_psc19_new_ident[1:*]
      concordiasi_psc19_new_jul=concordiasi_psc19_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(8,concordiasi_psc19_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc19_new_lat]
      longitudes=[longitudes,concordiasi_psc19_new_lon]
      daytime=[daytime,concordiasi_psc19_new_sid]
      ltime=[ltime,concordiasi_psc19_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc19_new_nprofs)]
      event=[event,concordiasi_psc19_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc19_new_nprofs)]
      sza=[sza,replicate(-1,concordiasi_psc19_new_nprofs)]
      julian=[julian,concordiasi_psc19_new_jul]
    endif
    if mipas_nprofs ne 0 then begin
      instrument=[instrument,replicate(9,mipas_nprofs)]
      latitudes=[latitudes,mipas_lat]
      longitudes=[longitudes,mipas_lon]
      daytime=[daytime,mipas_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,mipas_localtime]
      orbit=[orbit,replicate(-1l,mipas_nprofs)] ;this information *is* available in the MIPAS data, update this !!!
      event=[event,mipas_identifier]
      occ_type=[occ_type,replicate(0,mipas_nprofs)]
      sza=[sza,mipas_sza_tmp]
      mipas_h=floor(mipas_sid/3600.)
      mipas_min=floor((mipas_sid/3600.-mipas_h)*60.)
      mipas_s=round(mipas_sid-mipas_h*3600.-mipas_min*60.)
      julian=[julian,julday(month,day,year,mipas_h,mipas_min,mipas_s)-julian_reference]
    endif
    if poam2_nprofs ne 0 then begin
      instrument=[instrument,replicate(10,poam2_nprofs)]
      latitudes=[latitudes,poam2_latitudes]
      longitudes=[longitudes,poam2_longitudes]
      poam2_sid=round(poam2_time*3600)
      daytime=[daytime,poam2_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(poam2_time,poam2_longitudes)]
      orbit=[orbit,replicate(-1l,poam2_nprofs)]
      event=[event,poam2_occ_num]
      occ_type=[occ_type,poam2_occ_type]
      sza=[sza,replicate(-1,poam2_nprofs)]
      poam2_h=floor(poam2_sid/3600.)
      poam2_min=floor((poam2_sid/3600.-poam2_h)*60.)
      poam2_s=round(poam2_sid-poam2_h*3600.-poam2_min*60.)
      julian=[julian,julday(month,day,year,poam2_h,poam2_min,poam2_s)-julian_reference]
    endif
    if poam3_nprofs ne 0 then begin
      instrument=[instrument,replicate(11,poam3_nprofs)]
      latitudes=[latitudes,poam3_latitudes]
      longitudes=[longitudes,poam3_longitudes]
      poam3_sid=round(poam3_time*3600)
      daytime=[daytime,poam3_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(poam3_time,poam3_longitudes)]
      orbit=[orbit,replicate(-1l,poam3_nprofs)]
      event=[event,poam3_occ_num]
      occ_type=[occ_type,poam3_occ_type]
      sza=[sza,replicate(-1,poam3_nprofs)]
      poam3_h=floor(poam3_sid/3600.)
      poam3_min=floor((poam3_sid/3600.-poam3_h)*60.)
      poam3_s=round(poam3_sid-poam3_h*3600.-poam3_min*60.)
      julian=[julian,julday(month,day,year,poam3_h,poam3_min,poam3_s)-julian_reference]
    endif
    if sage2_nprofs ne 0 then begin
      instrument=[instrument,replicate(12,sage2_nprofs)]
      latitudes=[latitudes,sage2_latitudes]
      longitudes=[longitudes,sage2_longitudes]
      sage2_sid=round(sage2_time*3600)
      daytime=[daytime,sage2_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(sage2_time,sage2_longitudes)]
      orbit=[orbit,replicate(-1l,sage2_nprofs)]
      event=[event,sage2_occ_num]
      occ_type=[occ_type,sage2_occ_type]
      sza=[sza,replicate(-1,sage2_nprofs)]
      sage2_h=floor(sage2_sid/3600.)
      sage2_min=floor((sage2_sid/3600.-sage2_h)*60.)
      sage2_s=round(sage2_sid-sage2_h*3600.-sage2_min*60.)
      julian=[julian,julday(month,day,year,sage2_h,sage2_min,sage2_s)-julian_reference]
    endif
    if sage3_nprofs ne 0 then begin
      instrument=[instrument,replicate(13,sage3_nprofs)]
      latitudes=[latitudes,sage3_latitudes]
      longitudes=[longitudes,sage3_longitudes]
      sage3_sid=round(sage3_time*3600)
      daytime=[daytime,sage3_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(sage3_time,sage3_longitudes)]
      orbit=[orbit,replicate(-1l,sage3_nprofs)]
      event=[event,sage3_occ_num]
      occ_type=[occ_type,sage3_occ_type]
      sza=[sza,replicate(-1,sage3_nprofs)]
      sage3_h=floor(sage3_sid/3600.)
      sage3_min=floor((sage3_sid/3600.-sage3_h)*60.)
      sage3_s=round(sage3_sid-sage3_h*3600.-sage3_min*60.)
      julian=[julian,julday(month,day,year,sage3_h,sage3_min,sage3_s)-julian_reference]
    endif
    if haloe_nprofs ne 0 then begin
      instrument=[instrument,replicate(14,haloe_nprofs)]
      latitudes=[latitudes,haloe_latitudes]
      longitudes=[longitudes,haloe_longitudes]
      haloe_sid=round(haloe_time*3600)
      daytime=[daytime,haloe_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(haloe_time,haloe_longitudes)]
      orbit=[orbit,replicate(-1l,haloe_nprofs)]
      event=[event,haloe_occ_num]
      occ_type=[occ_type,haloe_occ_type]
      sza=[sza,replicate(-1,haloe_nprofs)]
      haloe_h=floor(haloe_sid/3600.)
      haloe_min=floor((haloe_sid/3600.-haloe_h)*60.)
      haloe_s=round(haloe_sid-haloe_h*3600.-haloe_min*60.)
      julian=[julian,julday(month,day,year,haloe_h,haloe_min,haloe_s)-julian_reference]
    endif
    if juelich_nprofs ne 0 then begin
      tmp_instr = 15*(juelich_instr+str_year eq 'PEX2005')+16*(juelich_instr+str_year eq 'REC2010')+17*(juelich_instr+str_year eq 'PEX2010')
      instrument=[instrument,tmp_instr]
      latitudes=[latitudes,juelich_latitudes]
      longitudes=[longitudes,juelich_longitudes]
      daytime=[daytime,juelich_daytime] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,juelich_ltime]
      orbit=[orbit,replicate(-1l,juelich_nprofs)]
      event=[event,juelich_ident]
      occ_type=[occ_type,replicate(0,juelich_nprofs)]
      sza=[sza,replicate(-1l,juelich_nprofs)]
      juelich_h=floor(juelich_daytime/3600.)
      juelich_min=floor((juelich_daytime/3600.-juelich_h)*60.)
      juelich_s=round(juelich_daytime-juelich_h*3600.-juelich_min*60.)
      julian=[julian,julday(month,day,year,juelich_h,juelich_min,juelich_s)-julian_reference]
    endif
    if sage1_nprofs ne 0 then begin
      instrument=[instrument,replicate(18,sage1_nprofs)]
      latitudes=[latitudes,sage1_latitudes]
      longitudes=[longitudes,sage1_longitudes]
      daytime=[daytime,sage1_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(sage1_time/24.,sage1_longitudes)]
      orbit=[orbit,replicate(-1l,sage1_nprofs)]
      event=[event,sage1_ident]
      occ_type=[occ_type,sage1_occ_type]
      sza=[sza,replicate(-1l,sage1_nprofs)]
      sage1_h=floor(sage1_time/3600.)
      sage1_min=floor((sage1_time/3600.-sage1_h)*60.)
      sage1_s=round(sage1_time-sage1_h*3600.-sage1_min*60.)
      julian=[julian,julday(month,day,year,sage1_h,sage1_min,sage1_s)-julian_reference]
    endif
    if uars_mls_nprofs ne 0 then begin
      instrument=[instrument,replicate(19,uars_mls_nprofs)]
      latitudes=[latitudes,uars_mls_latitudes]
      longitudes=[longitudes,uars_mls_longitudes]
      daytime=[daytime,uars_mls_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(uars_mls_time/24.,uars_mls_longitudes)]
      orbit=[orbit,replicate(-1l,uars_mls_nprofs)]
      event=[event,uars_mls_ident]
      occ_type=[occ_type,replicate(0l,uars_mls_nprofs)]
      sza=[sza,uars_mls_sza]
      uars_mls_h=floor(uars_mls_time/3600.)
      uars_mls_min=floor((uars_mls_time/3600.-uars_mls_h)*60.)
      uars_mls_s=round(uars_mls_time-uars_mls_h*3600.-uars_mls_min*60.)
      julian=[julian,julday(month,day,year,uars_mls_h,uars_mls_min,uars_mls_s)-julian_reference]
    endif
    if smiles_nprofs ne 0 then begin
      instrument=[instrument,replicate(20,smiles_nprofs)]
      latitudes=[latitudes,smiles_latitudes]
      longitudes=[longitudes,smiles_longitudes]
      daytime=[daytime,smiles_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,smiles_ltime]
      orbit=[orbit,replicate(-1l,smiles_nprofs)]
      event=[event,smiles_ident]
      occ_type=[occ_type,replicate(0l,smiles_nprofs)]
      sza=[sza,smiles_sza]
      smiles_h=floor(smiles_time/3600.)
      smiles_min=floor((smiles_time/3600.-smiles_h)*60.)
      smiles_s=round(smiles_time-smiles_h*3600.-smiles_min*60.)
      julian=[julian,julday(month,day,year,smiles_h,smiles_min,smiles_s)-julian_reference]
    endif
    if lims_nprofs ne 0 then begin
      instrument=[instrument,replicate(21,lims_nprofs)]
      latitudes=[latitudes,lims_latitudes]
      longitudes=[longitudes,lims_longitudes]
      daytime=[daytime,lims_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(lims_time/24.,lims_longitudes)]
      orbit=[orbit,replicate(-1l,lims_nprofs)]
      event=[event,lims_ident]
      occ_type=[occ_type,replicate(0l,lims_nprofs)]
      sza=[sza,lims_sza]
      lims_h=floor(lims_time/3600.)
      lims_min=floor((lims_time/3600.-lims_h)*60.)
      lims_s=round(lims_time-lims_h*3600.-lims_min*60.)
      julian=[julian,julday(month,day,year,lims_h,lims_min,lims_s)-julian_reference]
    endif
    if concordiasi_psc14_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc14_longitudes = (concordiasi_psc14_longitudes lt 0)*(360+concordiasi_psc14_longitudes) + (concordiasi_psc14_longitudes ge 0)*concordiasi_psc14_longitudes
      concordiasi_psc14_new_nprofs=0 & concordiasi_psc14_new_lat=0 & concordiasi_psc14_new_lon=0 & concordiasi_psc14_new_sid=0 & concordiasi_psc14_new_ltime=0 & concordiasi_psc14_new_ident=0 & concordiasi_psc14_new_jul=0
      for iprof=0,concordiasi_psc14_nprofs-1 do begin
        west_lon = (concordiasi_psc14_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc14_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc14_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc14_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlats eq 0 then stop
        concordiasi_psc14_new_nprofs+=4
        concordiasi_psc14_new_lat=[concordiasi_psc14_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc14_new_lon=[concordiasi_psc14_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc14_new_sid=[concordiasi_psc14_new_sid,replicate(concordiasi_psc14_daytime[iprof],4)]
        concordiasi_psc14_new_ltime=[concordiasi_psc14_new_ltime,replicate(concordiasi_psc14_ltime[iprof],4)]
        concordiasi_psc14_new_ident=[concordiasi_psc14_new_ident,14000000+concordiasi_psc14_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc14_h=floor(concordiasi_psc14_daytime[iprof]/3600.)
        concordiasi_psc14_min=floor((concordiasi_psc14_daytime[iprof]/3600.-concordiasi_psc14_h)*60.)
        concordiasi_psc14_s=round(concordiasi_psc14_daytime[iprof]-concordiasi_psc14_h*3600.-concordiasi_psc14_min*60.)
        concordiasi_psc14_new_jul=[concordiasi_psc14_new_jul,replicate((julday(month,day,year,concordiasi_psc14_h,concordiasi_psc14_min,concordiasi_psc14_s)-julian_reference),4)]
      endfor
      concordiasi_psc14_new_lat=concordiasi_psc14_new_lat[1:*]
      concordiasi_psc14_new_lon=concordiasi_psc14_new_lon[1:*] & concordiasi_psc14_new_sid=concordiasi_psc14_new_sid[1:*]
      concordiasi_psc14_new_ltime=concordiasi_psc14_new_ltime[1:*] & concordiasi_psc14_new_ident=concordiasi_psc14_new_ident[1:*]
      concordiasi_psc14_new_jul=concordiasi_psc14_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(22,concordiasi_psc14_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc14_new_lat]
      longitudes=[longitudes,concordiasi_psc14_new_lon]
      daytime=[daytime,concordiasi_psc14_new_sid]
      ltime=[ltime,concordiasi_psc14_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc14_new_nprofs)]
      event=[event,concordiasi_psc14_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc14_new_nprofs)]
      sza=[sza,(transpose(rebin(concordiasi_psc14_sza,concordiasi_psc14_nprofs,4)))[*]] ; Ferkelzeile: a=[234,78,54] & print, (transpose(rebin(a,n_elements(a),4)))[*]
      julian=[julian,concordiasi_psc14_new_jul]
    endif
    if concordiasi_psc15_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc15_longitudes = (concordiasi_psc15_longitudes lt 0)*(360+concordiasi_psc15_longitudes) + (concordiasi_psc15_longitudes ge 0)*concordiasi_psc15_longitudes
      concordiasi_psc15_new_nprofs=0 & concordiasi_psc15_new_lat=0 & concordiasi_psc15_new_lon=0 & concordiasi_psc15_new_sid=0 & concordiasi_psc15_new_ltime=0 & concordiasi_psc15_new_ident=0 & concordiasi_psc15_new_jul=0
      for iprof=0,concordiasi_psc15_nprofs-1 do begin
        west_lon = (concordiasi_psc15_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc15_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc15_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc15_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlats eq 0 then stop
        concordiasi_psc15_new_nprofs+=4
        concordiasi_psc15_new_lat=[concordiasi_psc15_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc15_new_lon=[concordiasi_psc15_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc15_new_sid=[concordiasi_psc15_new_sid,replicate(concordiasi_psc15_daytime[iprof],4)]
        concordiasi_psc15_new_ltime=[concordiasi_psc15_new_ltime,replicate(concordiasi_psc15_ltime[iprof],4)]
        concordiasi_psc15_new_ident=[concordiasi_psc15_new_ident,15000000+concordiasi_psc15_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc15_h=floor(concordiasi_psc15_daytime[iprof]/3600.)
        concordiasi_psc15_min=floor((concordiasi_psc15_daytime[iprof]/3600.-concordiasi_psc15_h)*60.)
        concordiasi_psc15_s=round(concordiasi_psc15_daytime[iprof]-concordiasi_psc15_h*3600.-concordiasi_psc15_min*60.)
        concordiasi_psc15_new_jul=[concordiasi_psc15_new_jul,replicate((julday(month,day,year,concordiasi_psc15_h,concordiasi_psc15_min,concordiasi_psc15_s)-julian_reference),4)]
      endfor
      concordiasi_psc15_new_lat=concordiasi_psc15_new_lat[1:*]
      concordiasi_psc15_new_lon=concordiasi_psc15_new_lon[1:*] & concordiasi_psc15_new_sid=concordiasi_psc15_new_sid[1:*]
      concordiasi_psc15_new_ltime=concordiasi_psc15_new_ltime[1:*] & concordiasi_psc15_new_ident=concordiasi_psc15_new_ident[1:*]
      concordiasi_psc15_new_jul=concordiasi_psc15_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(23,concordiasi_psc15_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc15_new_lat]
      longitudes=[longitudes,concordiasi_psc15_new_lon]
      daytime=[daytime,concordiasi_psc15_new_sid]
      ltime=[ltime,concordiasi_psc15_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc15_new_nprofs)]
      event=[event,concordiasi_psc15_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc15_new_nprofs)]
      sza=[sza,(transpose(rebin(concordiasi_psc15_sza,concordiasi_psc15_nprofs,4)))[*]] ; Ferkelzeile: a=[234,78,54] & print, (transpose(rebin(a,n_elements(a),4)))[*]
      julian=[julian,concordiasi_psc15_new_jul]
    endif
    if concordiasi_psc18_nprofs ne 0 then begin
      ;add neigbors (1. determine lat and lon neighbors, keep closest time, 2. add profiles and replicate constant data 4 times)
      concordiasi_psc18_longitudes = (concordiasi_psc18_longitudes lt 0)*(360+concordiasi_psc18_longitudes) + (concordiasi_psc18_longitudes ge 0)*concordiasi_psc18_longitudes
      concordiasi_psc18_new_nprofs=0 & concordiasi_psc18_new_lat=0 & concordiasi_psc18_new_lon=0 & concordiasi_psc18_new_sid=0 & concordiasi_psc18_new_ltime=0 & concordiasi_psc18_new_ident=0 & concordiasi_psc18_new_jul=0
      for iprof=0,concordiasi_psc18_nprofs-1 do begin
        west_lon = (concordiasi_psc18_longitudes[iprof]-waccm_lon_step) mod 360
        east_lon = (concordiasi_psc18_longitudes[iprof]+waccm_lon_step) mod 360
        if west_lon lt east_lon then new_lons = waccm_lon[where(waccm_lon ge west_lon and waccm_lon le east_lon,nnewlons)] & if nnewlons eq 0 then stop
        if west_lon gt east_lon then new_lons = reverse(waccm_lon[where(waccm_lon ge west_lon or waccm_lon le east_lon,nnewlons)]) & if nnewlons eq 0 then stop
        new_lats = waccm_lat[where(waccm_lat ge concordiasi_psc18_latitudes[iprof]-waccm_lat_step and waccm_lat le concordiasi_psc18_latitudes[iprof]+waccm_lat_step,nnewlats)] & if nnewlats eq 0 then stop
        concordiasi_psc18_new_nprofs+=4
        concordiasi_psc18_new_lat=[concordiasi_psc18_new_lat,new_lats[1],new_lats[0],new_lats[0],new_lats[1]]
        concordiasi_psc18_new_lon=[concordiasi_psc18_new_lon,new_lons[0],new_lons[0],new_lons[1],new_lons[1]]
        concordiasi_psc18_new_sid=[concordiasi_psc18_new_sid,replicate(concordiasi_psc18_daytime[iprof],4)]
        concordiasi_psc18_new_ltime=[concordiasi_psc18_new_ltime,replicate(concordiasi_psc18_ltime[iprof],4)]
        concordiasi_psc18_new_ident=[concordiasi_psc18_new_ident,18000000+concordiasi_psc18_ident[iprof]*100+[1,2,3,4]]
        concordiasi_psc18_h=floor(concordiasi_psc18_daytime[iprof]/3600.)
        concordiasi_psc18_min=floor((concordiasi_psc18_daytime[iprof]/3600.-concordiasi_psc18_h)*60.)
        concordiasi_psc18_s=round(concordiasi_psc18_daytime[iprof]-concordiasi_psc18_h*3600.-concordiasi_psc18_min*60.)
        concordiasi_psc18_new_jul=[concordiasi_psc18_new_jul,replicate((julday(month,day,year,concordiasi_psc18_h,concordiasi_psc18_min,concordiasi_psc18_s)-julian_reference),4)]
      endfor
      concordiasi_psc18_new_lat=concordiasi_psc18_new_lat[1:*]
      concordiasi_psc18_new_lon=concordiasi_psc18_new_lon[1:*] & concordiasi_psc18_new_sid=concordiasi_psc18_new_sid[1:*]
      concordiasi_psc18_new_ltime=concordiasi_psc18_new_ltime[1:*] & concordiasi_psc18_new_ident=concordiasi_psc18_new_ident[1:*]
      concordiasi_psc18_new_jul=concordiasi_psc18_new_jul[1:*]
      ;concatenate
      instrument=[instrument,replicate(24,concordiasi_psc18_new_nprofs)]
      latitudes=[latitudes,concordiasi_psc18_new_lat]
      longitudes=[longitudes,concordiasi_psc18_new_lon]
      daytime=[daytime,concordiasi_psc18_new_sid]
      ltime=[ltime,concordiasi_psc18_new_ltime]
      orbit=[orbit,replicate(-1l,concordiasi_psc18_new_nprofs)]
      event=[event,concordiasi_psc18_new_ident]
      occ_type=[occ_type,replicate(0,concordiasi_psc18_new_nprofs)]
      sza=[sza,(transpose(rebin(concordiasi_psc18_sza,concordiasi_psc18_nprofs,4)))[*]] ; Ferkelzeile: a=[234,78,54] & print, (transpose(rebin(a,n_elements(a),4)))[*]
      julian=[julian,concordiasi_psc18_new_jul]
    endif
    if sbuv_nprofs ne 0 then begin
      tmp_instr_num = 24+sbuv_sat_num
      instrument=[instrument,tmp_instr_num]
      latitudes=[latitudes,sbuv_latitudes]
      longitudes=[longitudes,sbuv_longitudes]
      daytime=[daytime,sbuv_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,sbuv_ltime]
      orbit=[orbit,sbuv_orbit_num]
      tmp_instr_day = current_julday - date2julday(19700410l*(sbuv_sat_num eq 1) + 19781031l*(sbuv_sat_num eq 2) + 19850202l*(sbuv_sat_num eq 3) + $
                                                   19881201l*(sbuv_sat_num eq 4) + 19950205l*(sbuv_sat_num eq 5) + 20001003l*(sbuv_sat_num eq 6) + $
                                                   20020711l*(sbuv_sat_num eq 7) + 20050605l*(sbuv_sat_num eq 8) + 20090223l*(sbuv_sat_num eq 9)) +1
      if min(tmp_instr_day) lt 1 then stop
      tmp_event_num = tmp_instr_num*100000000l + tmp_instr_day*10000l + lindgen(sbuv_nprofs) 
      event=[event,tmp_event_num]
      occ_type=[occ_type,replicate(0,sbuv_nprofs)]
      sza=[sza,sbuv_sza]
      sbuv_h=floor(sbuv_time/3600.)
      sbuv_min=floor((sbuv_time/3600.-sbuv_h)*60.)
      sbuv_s=round(sbuv_time-sbuv_h*3600.-sbuv_min*60.)
      julian=[julian,julday(month,day,year,sbuv_h,sbuv_min,sbuv_s)-julian_reference]
    endif
    if smr_nprofs ne 0 then begin
      instrument=[instrument,replicate(34,smr_nprofs)]
      latitudes=[latitudes,smr_latitudes]
      longitudes=[longitudes,smr_longitudes]
      daytime=[daytime,smr_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,smr_ltime]
      orbit=[orbit,replicate(-1l,smr_nprofs)]
      event=[event,smr_ident]
      occ_type=[occ_type,replicate(0,smr_nprofs)]
      sza=[sza,replicate(-1l,smr_nprofs)]
      smr_h=floor(smr_time/3600.)
      smr_min=floor((smr_time/3600.-smr_h)*60.)
      smr_s=round(smr_time-smr_h*3600.-smr_min*60.)
      julian=[julian,julday(month,day,year,smr_h,smr_min,smr_s)-julian_reference]
    endif
    if sme_ag_nprofs ne 0 then begin
      instrument=[instrument,replicate(35,sme_ag_nprofs)]
      latitudes=[latitudes,sme_ag_latitudes]
      longitudes=[longitudes,sme_ag_longitudes]
      daytime=[daytime,sme_ag_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,sme_ag_ltime]
      orbit=[orbit,sme_ag_orbit_num]
      event=[event,sme_ag_ident]
      occ_type=[occ_type,replicate(0,sme_ag_nprofs)]
      sza=[sza,sme_ag_sza]
      sme_ag_h=floor(sme_ag_time/3600.)
      sme_ag_min=floor((sme_ag_time/3600.-sme_ag_h)*60.)
      sme_ag_s=round(sme_ag_time-sme_ag_h*3600.-sme_ag_min*60.)
      julian=[julian,julday(month,day,year,sme_ag_h,sme_ag_min,sme_ag_s)-julian_reference]
    endif
    if sme_uv_nprofs ne 0 then begin
      instrument=[instrument,replicate(36,sme_uv_nprofs)]
      latitudes=[latitudes,sme_uv_latitudes]
      longitudes=[longitudes,sme_uv_longitudes]
      daytime=[daytime,sme_uv_time] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,sme_uv_ltime]
      orbit=[orbit,sme_uv_orbit_num]
      event=[event,sme_uv_ident]
      occ_type=[occ_type,replicate(0,sme_uv_nprofs)]
      sza=[sza,sme_uv_sza]
      sme_uv_h=floor(sme_uv_time/3600.)
      sme_uv_min=floor((sme_uv_time/3600.-sme_uv_h)*60.)
      sme_uv_s=round(sme_uv_time-sme_uv_h*3600.-sme_uv_min*60.)
      julian=[julian,julday(month,day,year,sme_uv_h,sme_uv_min,sme_uv_s)-julian_reference]
    endif
    if gomos_nprofs ne 0 then begin
      instrument=[instrument,replicate(37,gomos_nprofs)]
      latitudes=[latitudes,gomos_latitudes]
      longitudes=[longitudes,gomos_longitudes]
      daytime=[daytime,gomos_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,gomos_ltime]
      orbit=[orbit,replicate(-1l,gomos_nprofs)]
      event=[event,gomos_ident]
      occ_type=[occ_type,replicate(0,gomos_nprofs)]
      sza=[sza,gomos_sza]
      gomos_h=floor(gomos_sid/3600.)
      gomos_min=floor((gomos_sid/3600.-gomos_h)*60.)
      gomos_s=round(gomos_sid-gomos_h*3600.-gomos_min*60.)
      julian=[julian,julday(month,day,year,gomos_h,gomos_min,gomos_s)-julian_reference]
    endif

    instrument=instrument[1:*]
    latitudes=latitudes[1:*]
    longitudes=longitudes[1:*]
    daytime=daytime[1:*]
    ltime=ltime[1:*]
    orbit=orbit[1:*]
    event=event[1:*]
    occ_type=occ_type[1:*]
    sza=sza[1:*]
    julian=julian[1:*]
        
    ;data preparation
    yyyymmdd = replicate(current_date,nprofs)
    yyyyddd = replicate(long(year*1000l+current_doy),nprofs)
    longitudes = (longitudes lt 0)*(360+longitudes) + (longitudes ge 0)*longitudes
    index = where(longitudes ge 359.99, nindex) & if nindex ne 0 then longitudes[index]=0.
    isort=sort(daytime)
    instrument=instrument[isort]
    latitudes=latitudes[isort]
    longitudes=longitudes[isort]
    daytime=daytime[isort]
    ltime=ltime[isort]
    orbit=orbit[isort]
    event=event[isort]
    occ_type=occ_type[isort]
    yyyymmdd=yyyymmdd[isort]
    yyyyddd=yyyyddd[isort]
    sza=sza[isort]
    julian=julian[isort]
        
    ;intercept error codes (-999.99) and values out of range   
    nprofs_old = nprofs
    good = where(yyyymmdd ge start_date and yyyymmdd le end_date and $
                 daytime ge 0 and daytime le 86400 and $
                 latitudes ge -90 and latitudes le 90 and $
                 longitudes ge 0 and longitudes lt 360 and $
                 yyyyddd ge start_doy and yyyyddd le end_doy and $
                 ltime ge 0 and ltime le 24.5 and $
                 event ne -999.99,ngood,complement=bad,ncomplement=nbad)
    if ngood eq 0 then goto, nodata else begin
      yyyymmdd = yyyymmdd[good]
      daytime = daytime[good]
      latitudes = latitudes[good]
      longitudes = longitudes[good]
      yyyyddd = yyyyddd[good]
      ltime = ltime[good]
      instrument = instrument[good]
      orbit=orbit[good]
      event=event[good]
      occ_type=occ_type[good]
      sza=sza[good]
      julian=julian[good]
    endelse
    nprofs=long(ngood)
    
    mls_incl_nprofs = total(instrument eq 1,/int)
    ace_incl_nprofs = total(instrument eq 2,/int)
    hirdls_incl_nprofs = total(instrument eq 3,/int)
    saber_incl_nprofs = total(instrument eq 4,/int)
    sofie_incl_nprofs = total(instrument eq 5,/int)
    concordiasi_psc16_incl_nprofs = total(instrument eq 6,/int)/4
    concordiasi_psc17_incl_nprofs = total(instrument eq 7,/int)/4
    concordiasi_psc19_incl_nprofs = total(instrument eq 8,/int)/4
    mipas_incl_nprofs = total(instrument eq 9,/int)
    poam2_incl_nprofs = total(instrument eq 10,/int)
    poam3_incl_nprofs = total(instrument eq 11,/int)
    sage2_incl_nprofs = total(instrument eq 12,/int)
    sage3_incl_nprofs = total(instrument eq 13,/int)
    haloe_incl_nprofs = total(instrument eq 14,/int)
    pex2005_incl_nprofs = total(instrument eq 15,/int)
    rec2010_incl_nprofs = total(instrument eq 16,/int)
    pex2010_incl_nprofs = total(instrument eq 17,/int)
    sage1_incl_nprofs = total(instrument eq 18,/int)
    uars_mls_incl_nprofs = total(instrument eq 19,/int)
    smiles_incl_nprofs = total(instrument eq 20,/int)
    lims_incl_nprofs = total(instrument eq 21,/int)
    concordiasi_psc14_incl_nprofs = total(instrument eq 22,/int)/4
    concordiasi_psc15_incl_nprofs = total(instrument eq 23,/int)/4
    concordiasi_psc18_incl_nprofs = total(instrument eq 24,/int)/4
    buv_incl_nprofs = total(instrument eq 25,/int)
    sbuv_incl_nprofs = total(instrument eq 26,/int)
    sbuv2_incl_nprofs = total(instrument ge 27 and instrument le 33,/int)
    smr_incl_nprofs = total(instrument eq 34,/int)
    sme_ag_incl_nprofs = total(instrument eq 35,/int)
    sme_uv_incl_nprofs = total(instrument eq 36,/int)
    gomos_incl_nprofs = total(instrument eq 37,/int)
        
    ;write satellite_locations_per_day.txt
    if juelich_nprofs eq 0 then juelich_instr=''
    printf,2,str_date_code+$
             '  '+strcompress(string(nprofs_old),/r)+$
             '  '+strcompress(string(nprofs),/r)+$
             '  '+strcompress(string(mls_nprofs),/r)+$
             '  '+strcompress(string(mls_incl_nprofs),/r)+$
             '  '+strcompress(string(ace_nprofs),/r)+$
             '  '+strcompress(string(ace_incl_nprofs),/r)+$
             '  '+strcompress(string(hirdls_nprofs),/r)+$
             '  '+strcompress(string(hirdls_incl_nprofs),/r)+$
             '  '+strcompress(string(saber_nprofs),/r)+$
             '  '+strcompress(string(saber_incl_nprofs),/r)+$
             '  '+strcompress(string(sofie_nprofs),/r)+$
             '  '+strcompress(string(sofie_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc16_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc16_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc17_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc17_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc19_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc19_incl_nprofs),/r)+$
             '  '+strcompress(string(mipas_nprofs),/r)+$
             '  '+strcompress(string(mipas_incl_nprofs),/r)+$
             '  '+strcompress(string(poam2_nprofs),/r)+$
             '  '+strcompress(string(poam2_incl_nprofs),/r)+$
             '  '+strcompress(string(poam3_nprofs),/r)+$
             '  '+strcompress(string(poam3_incl_nprofs),/r)+$
             '  '+strcompress(string(sage2_nprofs),/r)+$
             '  '+strcompress(string(sage2_incl_nprofs),/r)+$
             '  '+strcompress(string(sage3_nprofs),/r)+$
             '  '+strcompress(string(sage3_incl_nprofs),/r)+$
             '  '+strcompress(string(haloe_nprofs),/r)+$
             '  '+strcompress(string(haloe_incl_nprofs),/r)+$
             '  '+strcompress(string(total(juelich_instr+str_year eq 'PEX2005',/int)),/r)+$
             '  '+strcompress(string(pex2005_incl_nprofs),/r)+$
             '  '+strcompress(string(total(juelich_instr+str_year eq 'REC2010',/int)),/r)+$
             '  '+strcompress(string(rec2010_incl_nprofs),/r)+$
             '  '+strcompress(string(total(juelich_instr+str_year eq 'PEX2010',/int)),/r)+$
             '  '+strcompress(string(pex2010_incl_nprofs),/r)+$
             '  '+strcompress(string(sage1_nprofs),/r)+$
             '  '+strcompress(string(sage1_incl_nprofs),/r)+$
             '  '+strcompress(string(uars_mls_nprofs),/r)+$
             '  '+strcompress(string(uars_mls_incl_nprofs),/r)+$
             '  '+strcompress(string(smiles_nprofs),/r)+$
             '  '+strcompress(string(smiles_incl_nprofs),/r)+$
             '  '+strcompress(string(lims_nprofs),/r)+$
             '  '+strcompress(string(lims_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc14_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc14_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc15_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc15_incl_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc18_nprofs),/r)+$
             '  '+strcompress(string(concordiasi_psc18_incl_nprofs),/r)+$
             '  '+strcompress(string(total(sbuv_sat_num eq 1,/int)),/r)+$ ;BUV
             '  '+strcompress(string(buv_incl_nprofs),/r)+$ ;BUV
             '  '+strcompress(string(total(sbuv_sat_num eq 2,/int)),/r)+$ ;SBUV
             '  '+strcompress(string(sbuv_incl_nprofs),/r)+$ ;SBUV
             '  '+strcompress(string(total(sbuv_sat_num ge 3,/int)),/r)+$ ;all SBUV2
             '  '+strcompress(string(sbuv2_incl_nprofs),/r)+$ ;all SBUV2
             '  '+strcompress(string(smr_nprofs),/r)+$
             '  '+strcompress(string(smr_incl_nprofs),/r)+$
             '  '+strcompress(string(sme_ag_nprofs),/r)+$
             '  '+strcompress(string(sme_ag_incl_nprofs),/r)+$
             '  '+strcompress(string(sme_uv_nprofs),/r)+$
             '  '+strcompress(string(sme_uv_incl_nprofs),/r)+$
             '  '+strcompress(string(gomos_nprofs),/r)+$
             '  '+strcompress(string(gomos_incl_nprofs),/r)
    
             
;    print,str_date_code+$
;         '  '+strcompress(string(mls_nprofs),/r)+$
;         '  '+strcompress(string(ace_nprofs),/r)+$
;         '  '+strcompress(string(hirdls_nprofs),/r)+$
;         '  '+strcompress(string(saber_nprofs),/r)+' done.' ;MLS, ACE, HIRDLS, SABER           
    
    ;write NetCDF file
    ncdf_varput, ncid, 'date', yyyymmdd, offset=offset, count=nprofs
    ncdf_varput, ncid, 'time', daytime, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lat', latitudes, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lon', longitudes, offset=offset, count=nprofs
    ncdf_varput, ncid, 'orbit_num', orbit, offset=offset, count=nprofs
    ncdf_varput, ncid, 'prof_num', event, offset=offset, count=nprofs
    ncdf_varput, ncid, 'instr_num', instrument, offset=offset, count=nprofs
    
    ncdf_varput, ncid, 'doy', yyyyddd, offset=offset, count=nprofs
    ncdf_varput, ncid, 'local_time', ltime, offset=offset, count=nprofs
    ncdf_varput, ncid, 'occ_type', occ_type, offset=offset, count=nprofs
    ncdf_varput, ncid, 'instr_sza', sza, offset=offset, count=nprofs
    ncdf_varput, ncid, 'julian', julian, offset=offset, count=nprofs
    
    ;increase offset for next day
    offset+=nprofs
    
    ;jump here if no data on this day
    nodata: 
  endfor

ncdf_close, ncid
close,2
print,'...done.'
end