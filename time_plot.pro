pro time_plot,setplot,mode,ifile0,ifile1,nday,theta,thal,xhal,yhal,xsat,ysat,$
phal,mhal,ch4hal,hfhal,h2ohal,o3hal,hclhal,no2hal,nohal,aerhal,here1,here2,$
here3,here4,halcomp,haldens,halmedr,haldisw,halconc,halsurf,halvolu,haleffr
;
!p.multi=[0,4,2]
;
if mode eq 'both' then begin
!psym=1
index=where(thal lt 1e10 and ch4hal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,ch4hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' CH4 vs Time'

!psym=2
index=where(thal lt 1e10 and ch4hal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,ch4hal(index) 
;
!psym=1
index=where(thal lt 1e10 and hfhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,hfhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HF vs Time'
!psym=2
index=where(thal lt 1e10 and hfhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,hfhal(index) 
;
!psym=1
index=where(thal lt 1e10 and h2ohal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,h2ohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' H2O vs Time'
!psym=2
index=where(thal lt 1e10 and h2ohal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,h2ohal(index) 
;
!psym=1
index=where(thal lt 1e10 and o3hal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,o3hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' O3 vs Time'
!psym=2
index=where(thal lt 1e10 and o3hal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,o3hal(index) 
;
!psym=1
index=where(thal lt 1e10 and hclhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,hclhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HCl vs Time'
!psym=2
index=where(thal lt 1e10 and hclhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,hclhal(index) 
;
!psym=1
index=where(thal lt 1e10 and no2hal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,no2hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO2 vs Time'
!psym=2
index=where(thal lt 1e10 and no2hal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,no2hal(index) 
;
!psym=1
index=where(thal lt 1e10 and nohal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,nohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO vs Time'
!psym=2
index=where(thal lt 1e10 and nohal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,nohal(index) 
;
!psym=1
index=where(thal lt 1e10 and aerhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,aerhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' AER vs Time'
!psym=2
index=where(thal lt 1e10 and aerhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,aerhal(index) 
if setplot eq 'x' then stop
;
!psym=1
index=where(thal lt 1e10 and halcomp lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halcomp(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-COMP vs Time'
!psym=2
index=where(thal lt 1e10 and halcomp lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,halcomp(index) 
;
!psym=1
index=where(thal lt 1e10 and haldens lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,haldens(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DENS vs Time'
!psym=2
index=where(thal lt 1e10 and haldens lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,haldens(index) 
;
!psym=1
index=where(thal lt 1e10 and halmedr lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halmedr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-MEDR vs Time'
!psym=2
index=where(thal lt 1e10 and halmedr lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,halmedr(index) 
;
!psym=1
index=where(thal lt 1e10 and haldisw lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,haldisw(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DISW vs Time'
!psym=2
index=where(thal lt 1e10 and haldisw lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,haldisw(index) 
;
!psym=1
index=where(thal lt 1e10 and halconc lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halconc(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-CONC vs Time'
!psym=2
index=where(thal lt 1e10 and halconc lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,halconc(index) 
;
!psym=1
index=where(thal lt 1e10 and halsurf lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halsurf(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-SURF vs Time'
!psym=2
index=where(thal lt 1e10 and halsurf lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,halsurf(index) 
;
!psym=1
index=where(thal lt 1e10 and halvolu lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halvolu(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-VOLU vs Time'
!psym=2
index=where(thal lt 1e10 and halvolu lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,halvolu(index) 
;
!psym=1
index=where(thal lt 1e10 and haleffr lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,haleffr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-EFFR vs Time'
!psym=2
index=where(thal lt 1e10 and haleffr lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,thal(index)/24.,haleffr(index) 
if setplot eq 'x' then stop
;
endif  ; mode endif
;
if mode eq 'ris' or mode eq 'set' then begin
if mode eq 'ris' then sym=1
if mode eq 'set' then sym=2
!psym=sym
index=where(thal lt 1e10 and ch4hal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,ch4hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' CH4 vs Time'
;
!psym=sym
index=where(thal lt 1e10 and hfhal lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,hfhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HF vs Time'
;
!psym=sym
index=where(thal lt 1e10 and h2ohal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,h2ohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' H2O vs Time'
;
!psym=sym
index=where(thal lt 1e10 and o3hal lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,o3hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' O3 vs Time'
;
!psym=sym
index=where(thal lt 1e10 and hclhal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,hclhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HCl vs Time'
;
!psym=sym
index=where(thal lt 1e10 and no2hal lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,no2hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO2 vs Time'
;
!psym=sym
index=where(thal lt 1e10 and nohal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,nohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO vs Time'
;
!psym=sym
index=where(thal lt 1e10 and aerhal lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,aerhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' AER vs Time'
if setplot eq 'x' then stop
;
!psym=sym
index=where(thal lt 1e10 and halcomp lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halcomp(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-COMP vs Time'
;
!psym=sym
index=where(thal lt 1e10 and haldens lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,haldens(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DENS vs Time'
;
!psym=sym
index=where(thal lt 1e10 and halmedr lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halmedr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-MEDR vs Time'
;
!psym=sym
index=where(thal lt 1e10 and haldisw lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,haldisw(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DISW vs Time'
;
!psym=sym
index=where(thal lt 1e10 and halconc lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,thal(index)/24.,halconc(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-CONC vs Time'
;
!psym=sym
index=where(thal lt 1e10 and halsurf lt 1e10)
if index(0) ne -1 then $
plot,thal(index)/24.,halsurf(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-SURF vs Time'
;
!psym=sym
index=where(thal lt 1e10 and halvolu lt 1e10)
plot,thal(index)/24.,halvolu(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-VOLU vs Time'
;
!psym=sym
index=where(thal lt 1e10 and haleffr lt 1e10)
plot,thal(index)/24.,haleffr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-EFFR vs Time'
if setplot eq 'x' then stop
;
endif   ; mode endif
return
end
