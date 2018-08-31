;-----------------------------------------------------------------------------------------------------------------------------
;WACCM AND MLS ELEVATED STRATOPAUSE COMPOSITE FIGURE FOR WACCM STRATOPAUSE CLIMATOLOGY
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 08/28/2012      |
;	 -------------------------------
;
;
;
;-----------------------------------------------------------------------------------------------------------------------------
;

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
@/Volumes/MacD68-1/france/idl_files/rd_waccm_nc3
@/Volumes/MacD68-1/france/idl_files/frac_index


;--------------------------




px1a = .01
px1b = 0.21
px2a = .23
px2b = .43

py1a = .54
py1b = .94
py2a = .11
py2b = .51


SETPLOT='ps'
read,'setplot',setplot
nxdim=1500
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
  	!p.charsize=2.    ; test with 1.5
  	!p.thick=3
  	!p.charthick=8
  	!y.thick=3
  	!x.thick=3
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

days = 0

months = ['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
MONTH = ['04','05','06','07','08','09','10','11']
titles1 = ['-(60-30)','-(30-0)','0-30','30-60','i)','k) Sep','m) Oct','Nov']
titles2 = ['b)','d)','f)','h)','j)','l)','n)','p)','r)']
titles3 = ['a)','c)','e)','g)','i)','k)','m)','o)','q)']




restore, '/Volumes/MacD68-1/france/ES_paper/Post_proces/WACCM_ES_daily_min_T_Z.sav'
for iday = 30, 60L do begin




	plot,[0,0],[0,0], position=[.001,.001,.4999,.999], xstyle = 4, ystyle = 4

	; ------------- TEMPERATURE PLOTS -----------------------




	xyouts, .03,.97,'Temperature (K)'
	xyouts, .56,.97,'Height (km)'
;xyouts, .3,.97,'WACCM', charsize = 4., charthick = 10

    	;----------------------plot code---------------------      ;x = where(lons eq 0.00)
      	;lons[x] = 360.
      	level1 = findgen(13)*3.+240.
      	level1a = findgen(26)*3.+222.
	level3 = findgen(1)+.2
	nlvls  = n_elements(level1)
	col1 = (1 + indgen(nlvls)) * icolmax / nlvls	; define colors



	dailymintemplon[*,iday] = mean(dailymintemplon[*,iday],/nan);
	dailymintemplat[*,iday] = mean(dailymintemplat[*,iday],/nan)
	dailyminheightlon[*,iday] = mean(dailyminheightlon[*,iday],/nan)
	dailyminheightlat[*,iday] = mean(dailyminheightlat[*,iday],/nan)
;	

  	if iday eq 30L then map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,	position = [px1a,py1a,px1b,py1b]
 	if iday gt 30L then map_set,90.,-90.,0,/ortho, /noeras,/noborder,position = [px1a,py1a,px1b,py1b]
 	
	if iday eq 30L then oplot, dailymintemplon[*,iday],dailymintemplat[*,iday], psym = 8, color = (iday-30.)*8., symsize = .3
	if iday gt 30L then oplot, dailymintemplon[*,iday],dailymintemplat[*,iday], psym = 8, color = (iday-30.)*8., symsize = .3


  	if iday eq 30L then map_set,90.,-90.,0,/ortho, /grid,/noeras,/noborder,/contin,position = [px2a,py1a,px2b,py1b]
 	if iday gt 30L then map_set,90.,-90.,0,/ortho, /noeras,/noborder,position = [px2a,py1a,px2b,py1b]


	if iday eq 30L then oplot, dailyminheightlon[*,iday],dailyminheightlat[*,iday], psym = 8, color = (iday-30.)*8., symsize = .3
	if iday gt 30L then oplot, dailyminheightlon[*,iday],dailyminheightlat[*,iday], psym = 8, color = (iday-30.)*8., symsize = .3

ENDFOR; iday





;--------------------------------  CONVERT .PS TO .PNG  ---------------------------------------


if setplot eq 'ps' then begin
	device, /close
	spawn,'pstopnm -dpi=300 -landscape idl.ps'	spawn,'pnmtopng idl001.ppm > /Volumes/MacD68-1/france/ES_paper/Figures/WACCM_ES_min_daily_evolution.png'
	spawn,'rm idl001.ppm idl.ps'
endif
;-------------------------------- END: CONVERT .PS TO .PNG  ---------------------------------------

end
