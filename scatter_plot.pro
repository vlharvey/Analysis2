pro scatter_plot,setplot,mode,ifile0,ifile1,theta,thal,xhal,yhal,xsat,ysat,$
phal,mhal,ch4hal,hfhal,h2ohal,o3hal,hclhal,no2hal,nohal,aerhal,here1,here2,$
here3,here4,halcomp,haldens,halmedr,haldisw,halconc,halsurf,halvolu,haleffr
;
!p.multi=[0,4,2]
;
if mode eq 'both' then begin
!psym=1
index=where(ch4hal lt 1e10 and yhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),yhal(index),title=ifile0+'-'+ifile1+' '+theta+' Lat vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and yhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),yhal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and hfhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),hfhal(index),title=ifile0+'-'+ifile1+' '+theta+' HF vs CH4' 
!psym=2
index=where(ch4hal lt 1e10 and hfhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),hfhal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and h2ohal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),h2ohal(index),title=ifile0+'-'+ifile1+' '+theta+' H2O vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and h2ohal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),h2ohal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and o3hal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),o3hal(index),title=ifile0+'-'+ifile1+' '+theta+' O3 vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and o3hal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),o3hal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and hclhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),hclhal(index),title=ifile0+'-'+ifile1+' '+theta+' HCl vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and hclhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),hclhal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and no2hal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),no2hal(index),title=ifile0+'-'+ifile1+' '+theta+' NO2 vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and no2hal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),no2hal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and nohal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),nohal(index),title=ifile0+'-'+ifile1+' '+theta+' NO vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and nohal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),nohal(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and aerhal lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),aerhal(index),title=ifile0+'-'+ifile1+' '+theta+' AER vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and aerhal lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),aerhal(index) 
if setplot eq 'x' then stop
;
!psym=1
index=where(ch4hal lt 1e10 and halcomp lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halcomp(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-COMP vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and halcomp lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),halcomp(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and haldens lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),haldens(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-DENS vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and haldens lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),haldens(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and halmedr lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halmedr(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-MEDR vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and halmedr lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),halmedr(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and haldisw lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),haldisw(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-DISW vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and haldisw lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),haldisw(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and halconc lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halconc(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-CONC vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and halconc lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),halconc(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and halsurf lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halsurf(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-SURF vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and halsurf lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),halsurf(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and halvolu lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halvolu(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-VOLU vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and halvolu lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),halvolu(index) 
;
!psym=1
index=where(ch4hal lt 1e10 and haleffr lt 1e10 and mhal eq 0)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),haleffr(index),title=ifile0+'-'+ifile1+' '+theta+' HAL-EFFR vs CH4'
!psym=2
index=where(ch4hal lt 1e10 and haleffr lt 1e10 and mhal eq 1)
if index(0) ne -1 then $
oplot,ch4hal(index),haleffr(index) 
if setplot eq 'x' then stop
;
endif  ; mode endif
;
if mode eq 'ris' or mode eq 'set' then begin
if mode eq 'ris' then sym=1
if mode eq 'set' then sym=2
!psym=sym
index=where(ch4hal lt 1e10 and yhal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),yhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' Lat vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and hfhal lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),hfhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HF vs CH4' 
;
!psym=sym
index=where(ch4hal lt 1e10 and h2ohal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),h2ohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' H2O vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and o3hal lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),o3hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' O3 vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and hclhal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),hclhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HCl vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and no2hal lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),no2hal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO2 vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and nohal lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),nohal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' NO vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and aerhal lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),aerhal(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' AER vs CH4'
if setplot eq 'x' then stop
;
!psym=sym
index=where(ch4hal lt 1e10 and halcomp lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halcomp(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-COMP vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and haldens lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),haldens(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DENS vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and halmedr lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halmedr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-MEDR vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and haldisw lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),haldisw(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-DISW vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and halconc lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halconc(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-CONC vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and halsurf lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),halsurf(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-SURF vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and halvolu lt 1e10)
if index(0) eq -1 then begin
plot,[0,1,1,0,0],[0,0,1,1,0],/nodata 
xyouts,.25,.5,'NO DATA'
endif
if index(0) ne -1 then $
plot,ch4hal(index),halvolu(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-VOLU vs CH4'
;
!psym=sym
index=where(ch4hal lt 1e10 and haleffr lt 1e10)
if index(0) ne -1 then $
plot,ch4hal(index),haleffr(index),xrange=[0,nday],title=ifile0+'-'+ifile1+' '+theta+' HAL-EFFR vs CH4'
if setplot eq 'x' then stop
;
endif  ; mode endif
stop
return
end
