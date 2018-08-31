;
; GEOS-5 version
; plot polar projections and yz cross polar sections

@/Volumes/MacD68-1/france/idl_files/stddat			; Determines the number of days since Jan 1, 1956
@/Volumes/MacD68-1/france/idl_files/kgmt			; This function computes the Julian day number (GMT) from the
								;    day, month, and year information.
@/Volumes/MacD68-1/france/idl_files/ckday			; This routine changes the Julian day from 365(6 if leap yr)
								;    to 1 and increases the year, if necessary.
@/Volumes/MacD68-1/france/idl_files/kdate			; gives back kmn,kdy information from the Julian day #.
@/Volumes/MacD68-1/france/idl_files/rd_ukmo_nc3			; reads the data from nc3 files
@/Volumes/MacD68-1/france/idl_files/date2uars			; This code returns the UARS day given (jday,year) information.
@/Volumes/MacD68-1/france/idl_files/plotPosition		; defines positions for n number of plots
@/Volumes/MacD68-1/france/idl_files/rd_GEOS5_nc3




SETPLOT='ps'
read,'setplot',setplot
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
icolmax=250
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

pxa = [.1,.1,.1,.1,.1,.4,.4,.4,.4,.4,.7,.7,.7,.7,.7]
pxb = [.3,.3,.3,.3,.3,.6,.6,.6,.6,.6,.9,.9,.9,.9,.9]

pya = [.82,.64,.46,.28,.1, .82,.64,.46,.28,.1, .82,.64,.46,.28,.1]
pyb = [.97,.79,.61,.43,.25,.97,.79,.61,.43,.25,.97,.79,.61,.43,.25]

npan = n_elements(pxa)
; --- NH --------
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
MONTH = ['01','02','03','04','05','06','07','08','09','10','11','12']
restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/latband_pre_process_NH.sav'

restore, '/Volumes/MacD68-1/france/ES_paper/Post_process/WACCM_ES_75latband_pre_process.sav'

x = where(temp le 0., nx)
if nx gt 0L then temp[x] = !values.f_nan

antimarkcount = fltarr(n_elements(alon), n_elements(altitude))
vortexmarkcount = fltarr(n_elements(alon), n_elements(altitude))
WACCMtmean = fltarr(n_elements(alon), n_elements(altitude))

for ii = 0L, n_elements(alon) - 1L do begin
	for jj = 0L, n_elements(altitude) - 1L do begin
		
		b = where(mark[*,ii,jj] lt 0., nx)
		if nx ge 1L then antimarkcount[ii,jj] = nx

		b = where(mark[*,ii,jj] gt 0., nx)
		if nx ge 1L then vortexmarkcount[ii,jj] = nx
		
		WACCMtmean[ii,jj] = mean(temp[*,ii,jj],/nan)
	endfor
endfor
antimarkcount = antimarkcount /max(antimarkcount,/nan)
vortexmarkcount = vortexmarkcount /max(vortexmarkcount,/nan)



;--------------------------- Plot Code ---------------

    erase
loadct,39


         !type=2^2+2^3
          level1=200.+5.*findgen(13)
	nlvls = n_elements(level1)
      nlvls  = n_elements(level1)
      col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors

          	 contour,WACCMtmean,lon,altitude,levels=level1,/cell_fill,c_color=col1,color=0,charsize = .8,$
                 	yticks=7,xticks=6,xtickname = ['0','60','120','180','240','300','360'],yrange=[20,100.], $
			position = [.15,.5,.5,.9], ytitle = 'Altitude (km)', xtitle = 'Longitude', title = '55-65!E!4'+string("157B)+ '!X!N N!CLatitude DJF'
  		contour,WACCMtmean,lon,altitude, levels=level1,color=0,charsize=1.5,$
      		  min_value = -2., max_value =40000.,/overplot,  thick=1



loadct, 0

    contour,smooth(vortexmarkcount,7,/nan),lon,altitude,levels=[.4,.8],/overplot,color=0, thick = 10
  	contour,smooth(antimarkcount,7,/nan),lon,altitude,levels=[.1,.5],/overplot,color=250, thick = 10

  	
for i = 0L, n_elements(lon) - 1L do begin
	if i eq 0L then strat = fltarr(n_elements(lon))
	z = where(altitude gt 25. and altitude lt 90.)
		strat[i] = altitude[where(reform(waccmtmean[i,*]) eq max(waccmtmean[i,z],/nan))]
endfor

  	oplot, lon, smooth(strat, 7,/nan), linestyle = 0, thick = 12, color = 100


loadct, 39
    ; -----------------plot the color bar -----------------------

      !type=2^2+2^3+2^6			; no y title or ticsks
      imin=min(level1)
      imax=max(level1)
      slab=' '+strarr(n_elements(level1))


     !p.title = ' '
     plot,[0,1],[0,0],yrange=[0,10],xrange=[0,1],/noeras,xticks=n_elements(level1)-1L,$
     position = [.25,.4,.8,.42],xstyle=1,xtickname=slab,$
          xtitle = 'Mean Temperature (K)', charsize = 1.2


      ybox=[0,10,10,0,0]

      x2=0.
      for j=1,n_elements(col1)-1 do begin
        dx= 1./(nlvls-1.)
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j-1)
        x2=x2+dx
      endfor

loadct,0
     slab=strcompress(string(format='(f7.3)',level1),/remove_all)
		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255

		slabcolor = fltarr(n_elements(level1))*0.
		slabcolor[0:2] = 255
      x1=0.+dx/2
for i=0L,n_elements(slab)-2L do begin
   slab0=slab(i)
   flab0=float(slab(i))
      slab0=strcompress(string(format='(i5.1)',flab0),/remove_all)
      xyouts,x1-dx/2,.76,slab0,charsize=1.,/data,color=slabcolor[i],charthick=2
   x1=x1+dx
endfor


   if setplot eq 'ps' then begin
        device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'
	spawn,$
	'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/WACCM_T_75N_Z_lon.png'
	spawn,'rm idl001.ppm idl.ps'
    endif

end