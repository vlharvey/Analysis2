;produce ACE-FTS raw data files
;one per year with all species
data_path='D:\SOSST\ACE original\v3.0\'
out_path='D:\SOSST\ACE original nc\'

start_occ=long(2221)
end_occ=long(48557) ;UPDATE THIS NUMBER BEFORE RUNNING !!! ;last update on 09/26/2012
str_year='2004' ;this is the year the code will start with

;;use this ASCII template for ACE-FTS v2.2 
;tpl = {VERSION:        1.00000, $ ;always the same
;      DATASTART:      12, $ ;line in which the data starts, below header and column labels
;      DELIMITER:      32b,$ ;delimiter between columns, note that the 'b' is necessary but won't be displayed in the output in 1.
;      MISSINGVALUE:   !values.f_nan, $ ;error code: for NaNs use !values.f_nan, or -99,...
;      COMMENTSYMBOL:  '', $ ;the symbol that indicates a comment, like the semicolon in IDL
;      FIELDCOUNT:     71, $ ;number of columns
;      FIELDTYPES:     [4,4,3,replicate(4,68)], $ ;variable type of each column, 7=string, 4=float, 2=integer (comprehensive table will all types can be found in the help of SIZE)
;      FIELDNAMES:     ['z','T','T_fit','P_atm','dens','H2O','H2O_err','O3','O3_err','N2O','N2O_err','CO','CO_err','CH4','CH4_err','NO','NO_err',$
;                       'NO2','NO2_err','HNO3','HNO3_err','HF','HF_err','HCl','HCl_err','OCS','OCS_err','N2O5','N2O5_err','ClONO2','ClONO2_err',$
;                       'HCN','HCN_err','CH3Cl','CH3Cl_err','CF4','CF4_err','CCl2F2','CCl2F2_err','CCl3F','CCl3F_err','COF2','COF2_err','C2H6',$
;                       'C2H6_err','C2H2','C2H2_err','CHF2Cl','CHF2Cl_err','HCOOH','HCOOH_err','SF6','SF6_err','ClO','ClO_err','HO2NO2',$
;                       'HO2NO2_err','H2O2','H2O2_err','HOCl','HOCl_err','H2CO','H2CO_err','CCl4','CCl4_err','N2','N2_err','CFC13','CFC13_err',$
;                       'HCFC142b','HCFC142b_err'], $ ;column names (those will be the variable names in IDL)
;      FIELDLOCATIONS: [3,8,17,22,33,43,55,67,79,91,103,115,127,139,151,163,175,187,199,211,223,235,247,259,271,283,295,307,319,331,343,355,367,$
;                       379,391,403,415,427,439,451,463,475,487,499,511,523,535,547,559,571,583,595,607,619,631,643,655,667,679,691,703,715,727,$
;                       739,751,763,775,787,799,811,823], $ ;this indicates the beginning of the next column, just copy from the output in 1.
;      FIELDGROUPS:    [indgen(71)]} ;sequential number, almost always equal to indgen(FIELDCOUNT)
      
;use this ASCII template for ACE-FTS v3.0 
tpl = {VERSION:        1.00000, $ ;always the same
      DATASTART:      12, $ ;line in which the data starts, below header and column labels
      DELIMITER:      32b,$ ;delimiter between columns, note that the 'b' is necessary but won't be displayed in the output in 1.
      MISSINGVALUE:   !values.f_nan, $ ;error code: for NaNs use !values.f_nan, or -99,...
      COMMENTSYMBOL:  '', $ ;the symbol that indicates a comment, like the semicolon in IDL
      FIELDCOUNT:     81, $ ;number of columns
      FIELDTYPES:     [4,4,3,replicate(4,78)], $ ;variable type of each column, 7=string, 4=float, 2=integer (comprehensive table will all types can be found in the help of SIZE)
      FIELDNAMES:     ['z','T','T_fit','P_atm','dens','H2O','H2O_err','O3','O3_err','N2O','N2O_err','CO','CO_err','CH4','CH4_err','NO','NO_err',$
                       'NO2','NO2_err','HNO3','HNO3_err','HF','HF_err','HCl','HCl_err','OCS','OCS_err','N2O5','N2O5_err','ClONO2','ClONO2_err',$
                       'HCN','HCN_err','CH3Cl','CH3Cl_err','CF4','CF4_err','CCl2F2','CCl2F2_err','CCl3F','CCl3F_err','COF2','COF2_err',$
                       'COCl2','COCl2_err','COClF','COClF_err',$ ;new in v3
                       'C2H6',$
                       'C2H6_err','C2H2','C2H2_err','CHF2Cl','CHF2Cl_err','HCOOH','HCOOH_err','SF6','SF6_err',$
                       ;ClO and ClO_err are not in ACE-FTS v3
                       'HO2NO2',$
                       'HO2NO2_err','H2O2','H2O2_err',$
                       ;HOCl and HOCl_err are not in ACE-FTS v3
                       'H2CO','H2CO_err',$
                       'CH3OH','CH3OH_err',$ ;new in ACE-FTS v3
                       'CCl4','CCl4_err','N2','N2_err',$
                       'O2','O2_err',$ ;new in ACE-FTS v3
                       'CFC13','CFC13_err',$ ;CFC113 is the correct label, in v2.2 it is sometimes erroneously named CFC13
                       'HCFC141b','HCFC141b_err',$ ;new in ACE-FTS v3
                       'HCFC142b','HCFC142b_err',$
                       'HFC134a','HFC134a_err','CO2','CO2_err'], $ ;new in ACE-FTS v3 ;column names (those will be the variable names in IDL)
      FIELDLOCATIONS: [3,8,17,22,33,43,55,67,79,91,103,115,127,139,151,163,175,187,199,211,223,235,247,259,271,283,295,307,319,331,$
                       343,355,367,379,391,403,415,427,439,451,463,475,487,499,511,523,535,547,559,571,583,595,607,619,631,643,655,$
                       667,679,691,703,715,727,739,751,763,775,787,799,811,823,835,847,859,871,883,895,907,919,931,943], $ ;this indicates the beginning of the next column, just copy from the output in 1.
      FIELDGROUPS:    [indgen(81)]} ;sequential number, almost always equal to indgen(FIELDCOUNT)      

;array of occultations flagged as bad measurements, see https://databace.uwaterloo.ca/validation/data_issues_table.php
bad_nocc = [1439+indgen(1454-1439+1),2206+indgen(2549-2206+1),2221+indgen(2255-2221+1),2471,2551+indgen(2830-2551+1),2968,4108+indgen(4281-4108+1),40578] ;updated on 20110919
bad_n2o_nocc = (read_ascii(data_path+'data_issues_v3.0_due_to_N2O_20110919.txt',template=ascii_template(data_path+'data_issues_v3.0_due_to_N2O_20110919.txt'))).field1
bad_n2o_nocc_type = strmid(bad_n2o_nocc,1,1)
bad_n2o_nocc = long(strmid(bad_n2o_nocc,2))

ioffset=0
start_new_file=1

openw,1,out_path+'log_ACE-FTS_v3.0_netcdf_c'+today()+'.txt'

for iocc=start_occ,end_occ-1l do begin
  filename='s?'+strcompress(string(iocc),/r)+'v3.0.asc'
  files=file_search(data_path+'\'+str_year+'*\'+filename,count=nfiles)
  if nfiles eq 0 then begin
    str_newyear = strcompress(str_year+1,/r)
    files=file_search(data_path+'\'+str_newyear+'*\'+filename,count=nfiles)
    if nfiles ne 0 then begin
      str_year = str_newyear
      ncdf_close, ncid
      start_new_file=1
    endif
  endif
  
  ;if strmid(files[0],strpos(files[0],'v3.0')+5,4) ne '2005' then continue ;USE THIS LINE TO DO ONE SPECIFIC YEAR ONLY !!!!
  if nfiles eq 0 then continue
  
  print, 'Starting ',iocc, ' in year ', str_year
  printf,1, 'Starting ',iocc, ' in year ', str_year
  
  for ifile=0l,nfiles-1l do begin
    data = read_ascii(files[ifile],template=tpl,header=header)
    ioffset+=1
    
    ;check consistency (are the fieldnames in the header the same as in the ascii template?)
    fieldnames=header[11]
    strput, fieldnames, 'P_atm  ', strpos(fieldnames,'P (atm)')
    fieldnames=strsplit(fieldnames,/extract)
    comparison=strcmp(fieldnames,tpl.fieldnames)
    if total(comparison) ne n_elements(fieldnames) and total(comparison[71:72]) eq 2 then stop
        
    ;extract information from header
    name = strmid(header[0],strpos(header[0],'|')+2)
    date = long(strmid(header[5],strpos(header[5],'|')+2,4)+strmid(header[5],strpos(header[5],'|')+2+5,2)+strmid(header[5],strpos(header[5],'|')+2+5+3,2))
    time = float(strmid(header[5],strpos(header[5],' ',/reverse_search)+1,2))+$
           float(strmid(header[5],strpos(header[5],' ',/reverse_search)+4,2))/60.+$
           float(strmid(header[5],strpos(header[5],' ',/reverse_search)+7,5))/3600.
    latitude = float(strmid(header[6],strpos(header[6],'|')+2))
    longitude = float(strmid(header[7],strpos(header[7],'|')+2))
    beta_angle = float(strmid(header[8],strpos(header[8],'|')+2))
    occ_type = strmid(files[ifile],strpos(files[ifile],'\',/reverse_search)+1,2)
    
    ;check consistency (is occ_type the same in filename and header?)
    if occ_type ne strmid(name,strpos(name,'.')+1,2) then stop
    occ_type = strmid(occ_type,1)
    
    if fix(date/10000) lt fix(str_year) then continue
    
    ;setup new output file
    if start_new_file then begin
      ncid = ncdf_create(out_path+'ACE-FTS_v3.0_'+str_year+'.nc',/clobber)
        ncdf_attput, ncid,'Description','ACE-FTS v3.0 dataset for '+str_year, /global
        ncdf_attput, ncid,'Author','file created by Matthias Brakebusch on '+today()+' using ACE_raw2netcdf.pro', /global
        
        nid = ncdf_dimdef(ncid, 'occultation_number', /unlimited)
        zid = ncdf_dimdef(ncid, 'Z', 150)
        
        vid = ncdf_vardef(ncid, 'occultation_number', [nid], /long)
        vid = ncdf_vardef(ncid, 'flag', [nid], /long)
        vid = ncdf_vardef(ncid, 'Z', [zid], /float)
        
        ;vid = ncdf_vardef(ncid, 'name', [nid], /char)
        vid = ncdf_vardef(ncid, 'date', [nid], /long)
        vid = ncdf_vardef(ncid, 'latitude', [nid], /float)
        vid = ncdf_vardef(ncid, 'longitude', [nid], /float)
        vid = ncdf_vardef(ncid, 'beta_angle', [nid], /float)
        
        vid = ncdf_vardef(ncid, 'occ_type', [nid], /char)
        vid = ncdf_vardef(ncid, 'time', [nid], /float)
        
        vid = ncdf_vardef(ncid, 'T', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'T_FIT', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'P_ATM', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'DENS', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'H2O', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'H2O_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'O3', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'O3_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2O', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2O_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CO', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CO_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CH4', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CH4_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'NO', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'NO_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'NO2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'NO2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HNO3', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HNO3_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HF', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HF_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCL', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCL_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'OCS', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'OCS_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2O5', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2O5_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CLONO2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CLONO2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCN', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCN_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CH3CL', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CH3CL_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CF4', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CF4_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CCL2F2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CCL2F2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CCL3F', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CCL3F_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'COF2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'COF2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'COCL2', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'COCL2_ERR', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'COCLF', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'COCLF_ERR', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'C2H6', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'C2H6_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'C2H2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'C2H2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CHF2CL', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CHF2CL_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCOOH', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCOOH_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'SF6', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'SF6_ERR', [zid,nid], /float)
;        vid = ncdf_vardef(ncid, 'CLO', [zid,nid], /float) ;not in v3.0
;        vid = ncdf_vardef(ncid, 'CLO_ERR', [zid,nid], /float) ;not in v3.0
        vid = ncdf_vardef(ncid, 'HO2NO2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HO2NO2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'H2O2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'H2O2_ERR', [zid,nid], /float)
;        vid = ncdf_vardef(ncid, 'HOCL', [zid,nid], /float) ;not in v3.0
;        vid = ncdf_vardef(ncid, 'HOCL_ERR', [zid,nid], /float) ;not in v3.0
        vid = ncdf_vardef(ncid, 'H2CO', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'H2CO_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CH3OH', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'CH3OH_ERR', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'CCL4', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CCL4_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'N2_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'O2', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'O2_ERR', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'CFC113', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'CFC113_ERR', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCFC141B', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'HCFC141B_ERR', [zid,nid], /float) ;new in v3.0
        vid = ncdf_vardef(ncid, 'HCFC142B', [zid,nid], /float)
        vid = ncdf_vardef(ncid, 'HCFC142B_ERR', [zid,nid], /float)
  
        ncdf_attput,ncid,'occultation_number','longname','# of occultation (first='+strcompress(string(start_occ),/r)+')'
          ncdf_attput,ncid,'occultation_number','units','1'
        ncdf_attput,ncid,'flag','longname','0 = good, 1 = do not use at all, 2 = do not use N2O'
          ncdf_attput,ncid,'flag','units','bitmask'
        ncdf_attput,ncid,'Z','longname','altitude' & ncdf_attput,ncid,'Z','units','km'
        ;ncdf_attput,ncid,'name','longname','name of event' & ncdf_attput,ncid,'name','units','N/A'
        ncdf_attput,ncid,'date','longname','date' & ncdf_attput,ncid,'date','units','yyyymmdd'
        ncdf_attput,ncid,'latitude','longname','latitude' & ncdf_attput,ncid,'latitude','units','deg'
        ncdf_attput,ncid,'longitude','longname','longitude' & ncdf_attput,ncid,'longitude','units','deg'
        ncdf_attput,ncid,'beta_angle','longname','beta angle of occultation' & ncdf_attput,ncid,'beta_angle','units','deg'        
        ncdf_attput,ncid,'occ_type','longname','type of occultation' & ncdf_attput,ncid,'occ_type','units','r = sunrise, s = sunset'
        ncdf_attput,ncid,'time','longname','UTC 24h time (float value)' & ncdf_attput,ncid,'time','units','h'
        ncdf_attput,ncid,'T','longname','temperature' & ncdf_attput,ncid,'T','units','K'
        ncdf_attput,ncid,'T_FIT','longname','temperature origin' & ncdf_attput,ncid,'T_FIT','units','1 = measured, 0 = MSIS model'
        ncdf_attput,ncid,'P_ATM','longname','pressure (1 atm = 1.01325 bar)' & ncdf_attput,ncid,'P_ATM','units','atm'
        ncdf_attput,ncid,'DENS','longname','atmospheric density' & ncdf_attput,ncid,'DENS','units','1/cm3'
        ncdf_attput,ncid,'H2O','longname','water vapor' & ncdf_attput,ncid,'H2O','units','ppv'
        ncdf_attput,ncid,'H2O_ERR','longname','water vapor error, -888 = not retrieved' & ncdf_attput,ncid,'H2O_ERR','units','Ppv'
        ncdf_attput,ncid,'O3','longname','ozone' & ncdf_attput,ncid,'O3','units','ppv'
        ncdf_attput,ncid,'O3_ERR','longname','ozone error, -888 = not retrieved' & ncdf_attput,ncid,'O3_ERR','units','Ppv'
        ncdf_attput,ncid,'N2O','longname','nitrous oxide' & ncdf_attput,ncid,'N2O','units','ppv'
        ncdf_attput,ncid,'N2O_ERR','longname','nitrous oxide error, -888 = not retrieved' & ncdf_attput,ncid,'N2O_ERR','units','Ppv'
        ncdf_attput,ncid,'CO','longname','carbon monoxide' & ncdf_attput,ncid,'CO','units','ppv'
        ncdf_attput,ncid,'CO_ERR','longname','carbon monoxide, -888 = not retrieved' & ncdf_attput,ncid,'CO_ERR','units','Ppv'
        ncdf_attput,ncid,'CH4','longname','methane' & ncdf_attput,ncid,'CH4','units','ppv'
        ncdf_attput,ncid,'CH4_ERR','longname','methane error, -888 = not retrieved' & ncdf_attput,ncid,'CH4_ERR','units','Ppv'
        ncdf_attput,ncid,'NO','longname','nitric oxide radical' & ncdf_attput,ncid,'NO','units','ppv'
        ncdf_attput,ncid,'NO_ERR','longname','ozone error, -888 = not retrieved' & ncdf_attput,ncid,'NO_ERR','units','Ppv'
        ncdf_attput,ncid,'NO2','longname','nitrogen dioxide' & ncdf_attput,ncid,'NO2','units','ppv'
        ncdf_attput,ncid,'NO2_ERR','longname','nitrogen dioxide error, -888 = not retrieved' & ncdf_attput,ncid,'NO2_ERR','units','Ppv'
        ncdf_attput,ncid,'HNO3','longname','nitric acid' & ncdf_attput,ncid,'HNO3','units','ppv'
        ncdf_attput,ncid,'HNO3_ERR','longname','nitric acid error, -888 = not retrieved' & ncdf_attput,ncid,'HNO3_ERR','units','Ppv'
        ncdf_attput,ncid,'HF','longname','hydrogen fluoride' & ncdf_attput,ncid,'HF','units','ppv'
        ncdf_attput,ncid,'HF_ERR','longname','hydrogen fluoride error, -888 = not retrieved' & ncdf_attput,ncid,'HF_ERR','units','Ppv'
        ncdf_attput,ncid,'HCL','longname','hydrogen chloride' & ncdf_attput,ncid,'HCL','units','ppv'
        ncdf_attput,ncid,'HCL_ERR','longname','hydrogen chloride error, -888 = not retrieved' & ncdf_attput,ncid,'HCL_ERR','units','Ppv'
        ncdf_attput,ncid,'OCS','longname','carbonyl sulfide' & ncdf_attput,ncid,'OCS','units','ppv'
        ncdf_attput,ncid,'OCS_ERR','longname','carbonyl sulfide error, -888 = not retrieved' & ncdf_attput,ncid,'OCS_ERR','units','Ppv'
        ncdf_attput,ncid,'N2O5','longname','nitrogen pentoxide' & ncdf_attput,ncid,'N2O5','units','ppv'
        ncdf_attput,ncid,'N2O5_ERR','longname','nitrogen pentoxide error, -888 = not retrieved' & ncdf_attput,ncid,'N2O5_ERR','units','Ppv'
        ncdf_attput,ncid,'CLONO2','longname','chlorine nitrate' & ncdf_attput,ncid,'CLONO2','units','ppv'
        ncdf_attput,ncid,'CLONO2_ERR','longname','chlorine nitrate error, -888 = not retrieved' & ncdf_attput,ncid,'CLONO2_ERR','units','Ppv'
        ncdf_attput,ncid,'HCN','longname','hydrogen cyanid' & ncdf_attput,ncid,'HCN','units','ppv'
        ncdf_attput,ncid,'HCN_ERR','longname','hydrogen cyanid error, -888 = not retrieved' & ncdf_attput,ncid,'HCN_ERR','units','Ppv'
        ncdf_attput,ncid,'CH3CL','longname','methyl chloride' & ncdf_attput,ncid,'CH3CL','units','ppv'
        ncdf_attput,ncid,'CH3CL_ERR','longname','methyl chloride error, -888 = not retrieved' & ncdf_attput,ncid,'CH3CL_ERR','units','Ppv'
        ncdf_attput,ncid,'CF4','longname','carbon tetrafluoride' & ncdf_attput,ncid,'CF4','units','ppv'
        ncdf_attput,ncid,'CF4_ERR','longname','carbon tetrafluoride error, -888 = not retrieved' & ncdf_attput,ncid,'CF4_ERR','units','Ppv'
        ncdf_attput,ncid,'CCL2F2','longname','dichlorodifluoromethane (CFC-12)' & ncdf_attput,ncid,'CCL2F2','units','ppv'
        ncdf_attput,ncid,'CCL2F2_ERR','longname','dichlorodifluoromethane error, -888 = not retrieved' & ncdf_attput,ncid,'CCL2F2_ERR','units','Ppv'
        ncdf_attput,ncid,'CCL3F','longname','trichlorofluoromethane (CFC-11)' & ncdf_attput,ncid,'CCL3F','units','ppv'
        ncdf_attput,ncid,'CCL3F_ERR','longname','trichlorofluoromethane error, -888 = not retrieved' & ncdf_attput,ncid,'CCL3F_ERR','units','Ppv'
        ncdf_attput,ncid,'COF2','longname','carbonyl fluoride' & ncdf_attput,ncid,'COF2','units','ppv'
        ncdf_attput,ncid,'COF2_ERR','longname','carbonyl fluoride error, -888 = not retrieved' & ncdf_attput,ncid,'COF2_ERR','units','Ppv'
        ncdf_attput,ncid,'COCL2','longname','carbonyl chloride' & ncdf_attput,ncid,'COCL2','units','ppv' ;new in v3.0
        ncdf_attput,ncid,'COCL2_ERR','longname','carbonyl chloride error, -888 = not retrieved' & ncdf_attput,ncid,'COCL2_ERR','units','Ppv' ;new in v3.0
        ncdf_attput,ncid,'COCLF','longname','carbonyl chlorofluoride' & ncdf_attput,ncid,'COCLF','units','ppv' ;new in v3.0
        ncdf_attput,ncid,'COCLF_ERR','longname','carbonyl chlorofluoride error, -888 = not retrieved' & ncdf_attput,ncid,'COCLF_ERR','units','Ppv' ;new in v3.0
        ncdf_attput,ncid,'C2H6','longname','ethane' & ncdf_attput,ncid,'C2H6','units','ppv'
        ncdf_attput,ncid,'C2H6_ERR','longname','ethane error, -888 = not retrieved' & ncdf_attput,ncid,'C2H6_ERR','units','Ppv'
        ncdf_attput,ncid,'C2H2','longname','acetylene' & ncdf_attput,ncid,'C2H2','units','ppv'
        ncdf_attput,ncid,'C2H2_ERR','longname','acetylene error, -888 = not retrieved' & ncdf_attput,ncid,'C2H2_ERR','units','Ppv'
        ncdf_attput,ncid,'CHF2CL','longname','chlorodifluoromethane (HCFC-22)' & ncdf_attput,ncid,'CHF2CL','units','ppv'
        ncdf_attput,ncid,'CHF2CL_ERR','longname','chlorodifluoromethane error, -888 = not retrieved' & ncdf_attput,ncid,'CHF2CL_ERR','units','Ppv'
        ncdf_attput,ncid,'HCOOH','longname','formic acid' & ncdf_attput,ncid,'HCOOH','units','ppv'
        ncdf_attput,ncid,'HCOOH_ERR','longname','formic acid error, -888 = not retrieved' & ncdf_attput,ncid,'HCOOH_ERR','units','Ppv'
        ncdf_attput,ncid,'SF6','longname','sulfur hexafluoride' & ncdf_attput,ncid,'SF6','units','ppv'
        ncdf_attput,ncid,'SF6_ERR','longname','sulfur hexafluoride error, -888 = not retrieved' & ncdf_attput,ncid,'SF6_ERR','units','Ppv'
;        ncdf_attput,ncid,'CLO','longname','chlorine monoxide' & ncdf_attput,ncid,'CLO','units','ppv' ;not in v3.0
;        ncdf_attput,ncid,'CLO_ERR','longname','chlorine monoxide error, -888 = not retrieved' & ncdf_attput,ncid,'CLO_ERR','units','Ppv' ;not in v3.0
        ncdf_attput,ncid,'HO2NO2','longname','peroxynitric acid' & ncdf_attput,ncid,'HO2NO2','units','ppv'
        ncdf_attput,ncid,'HO2NO2_ERR','longname','peroxynitric acid error, -888 = not retrieved' & ncdf_attput,ncid,'HO2NO2_ERR','units','Ppv'
        ncdf_attput,ncid,'H2O2','longname','hydrogen peroxide' & ncdf_attput,ncid,'H2O2','units','ppv'
        ncdf_attput,ncid,'H2O2_ERR','longname','hydrogen peroxide error, -888 = not retrieved' & ncdf_attput,ncid,'H2O2_ERR','units','Ppv'
;        ncdf_attput,ncid,'HOCL','longname','hypochlorous acid' & ncdf_attput,ncid,'HOCL','units','ppv' ;not in v3.0
;        ncdf_attput,ncid,'HOCL_ERR','longname','hypochlorous acid error, -888 = not retrieved' & ncdf_attput,ncid,'HOCL_ERR','units','Ppv' ;not in v3.0
        ncdf_attput,ncid,'H2CO','longname','formaldehyde' & ncdf_attput,ncid,'H2CO','units','ppv'
        ncdf_attput,ncid,'H2CO_ERR','longname','formaldehyde error, -888 = not retrieved' & ncdf_attput,ncid,'H2CO_ERR','units','Ppv'
        ncdf_attput,ncid,'CH3OH','longname','methanol' & ncdf_attput,ncid,'CH3OH','units','ppv' ;new in v3.0
        ncdf_attput,ncid,'CH3OH_ERR','longname','methanol error, -888 = not retrieved' & ncdf_attput,ncid,'CH3OH_ERR','units','Ppv' ;new in v3.0
        ncdf_attput,ncid,'CCL4','longname','carbon tetrachloride' & ncdf_attput,ncid,'CCL4','units','ppv'
        ncdf_attput,ncid,'CCL4_ERR','longname','carbon tetrachloride error, -888 = not retrieved' & ncdf_attput,ncid,'CCL4_ERR','units','Ppv'
        ncdf_attput,ncid,'N2','longname','nitrogen' & ncdf_attput,ncid,'N2','units','ppv'
        ncdf_attput,ncid,'N2_ERR','longname','nitrogen error, -888 = not retrieved' & ncdf_attput,ncid,'N2_ERR','units','Ppv'
        ncdf_attput,ncid,'O2','longname','oxygen' & ncdf_attput,ncid,'O2','units','ppv' ;new in v3.0
        ncdf_attput,ncid,'O2_ERR','longname','oxygen error, -888 = not retrieved' & ncdf_attput,ncid,'O2_ERR','units','Ppv' ;new in v3.0
        ncdf_attput,ncid,'CFC113','longname','trichlorotrifluoroethane (CFC-113)' & ncdf_attput,ncid,'CFC113','units','ppv'
        ncdf_attput,ncid,'CFC113_ERR','longname','trichlorotrifluoroethane error, -888 = not retrieved' & ncdf_attput,ncid,'CFC113_ERR','units','Ppv'
        ncdf_attput,ncid,'HCFC141B','longname','dichlorofluoroethane (HCFC141b)' & ncdf_attput,ncid,'HCFC141B','units','ppv'
        ncdf_attput,ncid,'HCFC141B_ERR','longname','dichlorofluoroethane error, -888 = not retrieved' & ncdf_attput,ncid,'HCFC141B_ERR','units','Ppv' ;new in v3.0
        ncdf_attput,ncid,'HCFC142B','longname','monochlorodifluoroethane (HCFC142b)' & ncdf_attput,ncid,'HCFC142B','units','ppv' ;new in v3.0
        ncdf_attput,ncid,'HCFC142B_ERR','longname','monochlorodifluoroethane error, -888 = not retrieved' & ncdf_attput,ncid,'HCFC142B_ERR','units','Ppv'
          
        ncdf_control, ncid, /endef
        
        ncdf_varput, ncid, 'Z', data.z
        z_reference = data.z
        start_new_file=0
        ioffset=0
    endif
    
    index=where(bad_nocc eq iocc,nindex)
    n2oindex=where(bad_n2o_nocc eq iocc and bad_n2o_nocc_type eq occ_type,nn2oindex)
    flag = 0 + 1*(nindex ne 0) + 2*(nn2oindex ne 0)
    
    if array_equal(z_reference,data.z) ne 1 then continue ;those are actually wrong temperatures overlapping into vertical axis column
    
    ;write data to NetCDF file
    ncdf_varput, ncid, 'occultation_number', iocc, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'flag', flag, offset=[ioffset], count=[1]
    ;ncdf_varput, ncid, 'name', name, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'date', date, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'latitude', latitude, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'longitude', longitude, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'beta_angle', beta_angle, offset=[ioffset], count=[1]    
    ncdf_varput, ncid, 'occ_type', occ_type, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'time', time, offset=[ioffset], count=[1]
    ncdf_varput, ncid, 'T', data.T, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'T_FIT', data.T_FIT, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'P_ATM', data.P_ATM, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'DENS', data.DENS, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'H2O', data.H2O, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'H2O_ERR', data.H2O_ERR, offset=[0,ioffset], count=[150,1]    
    ncdf_varput, ncid, 'O3', data.O3, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'O3_ERR', data.O3_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2O', data.N2O, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2O_ERR', data.N2O_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CO', data.CO, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CO_ERR', data.CO_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CH4', data.CH4, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CH4_ERR', data.CH4_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'NO', data.NO, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'NO_ERR', data.NO_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'NO2', data.NO2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'NO2_ERR', data.NO2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HNO3', data.HNO3, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HNO3_ERR', data.HNO3_ERR, offset=[0,ioffset], count=[150,1]    
    ncdf_varput, ncid, 'HF', data.HF, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HF_ERR', data.HF_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCL', data.HCL, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCL_ERR', data.HCL_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'OCS', data.OCS, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'OCS_ERR', data.OCS_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2O5', data.N2O5, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2O5_ERR', data.N2O5_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CLONO2', data.CLONO2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CLONO2_ERR', data.CLONO2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCN', data.HCN, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCN_ERR', data.HCN_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CH3CL', data.CH3CL, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CH3CL_ERR', data.CH3CL_ERR, offset=[0,ioffset], count=[150,1]    
    ncdf_varput, ncid, 'CF4', data.CF4, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CF4_ERR', data.CF4_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL2F2', data.CCL2F2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL2F2_ERR', data.CCL2F2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL3F', data.CCL3F, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL3F_ERR', data.CCL3F_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'COF2', data.COF2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'COF2_ERR', data.COF2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'C2H6', data.C2H6, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'C2H6_ERR', data.C2H6_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'C2H2', data.C2H2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'C2H2_ERR', data.C2H2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CHF2CL', data.CHF2CL, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CHF2CL_ERR', data.CHF2CL_ERR, offset=[0,ioffset], count=[150,1]    
    ncdf_varput, ncid, 'HCOOH', data.HCOOH, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCOOH_ERR', data.HCOOH_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'SF6', data.SF6, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'SF6_ERR', data.SF6_ERR, offset=[0,ioffset], count=[150,1]
;    ncdf_varput, ncid, 'CLO', data.CLO, offset=[0,ioffset], count=[150,1] ;not in v3.0
;    ncdf_varput, ncid, 'CLO_ERR', data.CLO_ERR, offset=[0,ioffset], count=[150,1] ;not in v3.0
    ncdf_varput, ncid, 'HO2NO2', data.HO2NO2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HO2NO2_ERR', data.HO2NO2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'H2O2', data.H2O2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'H2O2_ERR', data.H2O2_ERR, offset=[0,ioffset], count=[150,1]
;    ncdf_varput, ncid, 'HOCL', data.HOCL, offset=[0,ioffset], count=[150,1] ;not in v3.0
;    ncdf_varput, ncid, 'HOCL_ERR', data.HOCL_ERR, offset=[0,ioffset], count=[150,1] ;not in v3.0
    ncdf_varput, ncid, 'H2CO', data.H2CO, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'H2CO_ERR', data.H2CO_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL4', data.CCL4, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CCL4_ERR', data.CCL4_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2', data.N2, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'N2_ERR', data.N2_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CFC113', data.CFC13, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'CFC113_ERR', data.CFC13_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCFC142B', data.HCFC142B, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'HCFC142B_ERR', data.HCFC142B_ERR, offset=[0,ioffset], count=[150,1]
    ncdf_varput, ncid, 'COCL2', data.COCL2, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'COCL2_ERR', data.COCL2_ERR, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'COCLF', data.COCLF, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'COCLF_ERR', data.COCLF_ERR, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'CH3OH', data.CH3OH, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'CH3OH_ERR', data.CH3OH_ERR, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'O2', data.O2, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'O2_ERR', data.O2_ERR, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'HCFC141B', data.HCFC141B, offset=[0,ioffset], count=[150,1] ;new in v3.0
    ncdf_varput, ncid, 'HCFC141B_ERR', data.HCFC141B_ERR, offset=[0,ioffset], count=[150,1] ;new in v3.0
  endfor
    
endfor
ncdf_close, ncid
close,1

print,'...done.'
end