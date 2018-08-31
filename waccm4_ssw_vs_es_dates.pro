;
; WACCM SSW dates vs ES dates
;
loadct,39
device,decompose=0
mcolor=255
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill

restore,'Save_files/WACCM4_dTdy_Ubar_SSW_Climo.sav
index=where(th eq 900.)
NH_UBAR900=reform(NH_UBAR(*,index(0)))
smon=strmid(strcompress(yyyymmdd,/remove_all),4,2)
index=where(NH_UBAR900 lt 0. and (smon eq '01' or smon eq '02' or smon eq '12'))
esdayzeros=[20030224, 20041218, 20060204, 20081226, 20151221, 20231218, 20260224, 20261211, 20300123, 20320207, 20331221, 20370224, 20390226, 20410101, 20420104]
majorsswdays=YYYYMMDD(index)       
esmon=strmid(strcompress(esdayzeros,/remove_all),4,2)
esyr=strmid(strcompress(esdayzeros,/remove_all),0,4)
esdy=strmid(strcompress(esdayzeros,/remove_all),6,2)
sswmon=strmid(strcompress(majorsswdays,/remove_all),4,2)
sswyr=strmid(strcompress(majorsswdays,/remove_all),0,4)
sswdy=strmid(strcompress(majorsswdays,/remove_all),6,2)
sswjday=julday(long(sswmon),long(sswdy),long(sswyr))
esjday=julday(long(esmon),long(esdy),long(esyr))
;
; plot
;
plot,sswjday,psym=2,yrange=[min(sswjday)-10.,max(sswjday)+10.]
for ies=0,n_elements(esdayzeros)-1L do begin
    index=where(abs(esjday(ies)-sswjday) eq min(abs(esjday(ies)-sswjday)))
    oplot,[index(0),index(0)],[esjday(ies),esjday(ies)],psym=8,color=mcolor*.9
print,esjday(ies)-sswjday
stop
endfor
end
