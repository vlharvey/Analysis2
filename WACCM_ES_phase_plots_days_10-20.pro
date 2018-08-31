SETPLOT='ps'
;read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.1,0.6,0.1,0.6,0.1,0.6]
yorig=[0.7,0.7,0.4,0.4,0.1,0.1]

xlen=0.25
ylen=0.25
cbaryoff=0.02
cbarydel=0.01

a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill


if setplot eq 'ps' then begin
  xsize=nxdim/100.
  ysize=nydim/100.
  set_plot,'ps'
  device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
         /bold,/color,bits_per_pixel=8,/helvetica,filename='idl.ps'
  !p.charsize=1.25    ; test with 1.5
  !p.thick=2
  !p.charthick=5
  !p.charthick=5
  !y.thick=2
  !x.thick=2
endif

loadct,39
mcolor=!p.color
icolmax=255
mcolor=icolmax
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
;

inday = 0L

restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/Post_process_1_p.sav'
restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'
restore, '/Volumes/MacD68-1/france/ES_paper/Post_process/WACCM_ES_daily_max_T_Z.sav'
restore,'/Volumes/MacD68-1/france/ES_paper/Post_process/WACCM_ES_wave-1_phase_pre_process_5km.sav' ; ,amp60,amp70,amp75,phase60,phase70,phase75,Date

phase60 = phase60 - 180.
x = where(phase60 gt 360.,nx)
if nx gt 0L then phase60[x] = phase60[x]-360.

for iES = 0L, n_elements(dayzeros) - 1L do begin
	for iday = 0L, 29L do begin
	
		sdateindex = where(strmatch(dates, date[inday]) eq 1L)
		print, dates[sdateindex], sdateindex
		press = reform(waccmp[sdateindex,*,*])
		x = where(press lt 0.,nx)
		if nx gt 0L then press[x] = !values.f_nan
		stratopause = fltarr(n_elements(lon))
		x = where(abs(lat-60.) eq min(abs(lat-60.)))
		for jj = 0L, n_elements(lon) - 1L do stratopause[jj] = press[jj,x]
		stratopause = smooth(stratopause,5,/edge,/nan)
		
     
    ;----------------------plot code--------------------- 


	if iday ge 9L and iday lt 19L then begin

		if iday eq 9L then plot,/ylog,phase60[inday,*], pressure_grid,position=[.001,.1,.999,.999],xrange = [-180.,180],yrange = [100,.01], thick = 3,$
		title = 'ES days 10-20 event ' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2), color = 0,xstyle = 1, ystyle = 1
	
		if iday gt 9L then oplot,phase60[inday,*], pressure_grid, thick = 3, color = (iday-9L)/10.*255.
	
		oplot,  lon-180.,stratopause, linestyle = 1, thick = 2,  color = (iday-9L)/10.*255.
	endif
    
    inday = inday+1L
endfor

    ; -----------------plot the color bar-----------------------
      level1 = findgen(10)+10.
      nlvls  = n_elements(level1)
      col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors

  	 loadct, 39
    ;print, max(meanGEOS5strats)
      ;plot the color bar
      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level1)
      imax=max(level1)
      slab=' '+strarr(n_elements(level1))

      !p.title = ' '
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xticks=n_elements(level1)-1L,$
      position = [.2,.03,.8,.05],xstyle=1,xtickname=slab, xtitle = 'Day of Event'


      ybox=[0,10,10,0,0]

      x2=imin
      for j=1,n_elements(col1)-1 do begin
        dx= level1[j] - level1[j-1]
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

;;;;;;;;      
     slab=strcompress(string(format='(f7.3)',level1),/remove_all)
		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255

		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:3] = 255

    x1=min(level1)+dx/2 + dx

for i=1L,n_elements(slab)-1L do begin
   slab0=slab(i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(I5.2)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=.8,/data,color=slabcolor[i],charthick=1, orientation= 90.
   x1=x1+dx
endfor


  if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'
	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/60N_phase_plot_days_10-20_event_' + strtrim(strcompress(string(format='(I2.2)',ies + 1L)),2)+ '.png'
	spawn,'rm idl001.ppm idl.ps'
    endif
endfor
end
