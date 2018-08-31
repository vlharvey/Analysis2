;plt_l4_summary.pro -- Make example plots to check out the level 4 summary files.

print,'MUST DEFINE DFS FOR NH11 SEASON AND BEYOND!!!'
PRINT,' '
stop

;FOR FREQUENCY, WE ONLY NEED THE "ALL" FILES.

;THIS CODE WORKS WITH THE 1-DEGREE LATITUDE BIN FILES.
;LATITUDES GO FROM 50-51 TO 84-85 AND 85-84 TO 51-50, IN 1-DEGREE INCREMENTS.

pth = 'C:\CIPS\data\v4.20\level3c\'
ptho = 'C:\CIPS\analysis\v4.20\contour_plots\'

season='north_07'
read,'Input the season (nh07, nh08, nh09, nh10, nh11, sh07, sh08, sh09, sh10) >> ',season
season=strupcase(season)
case season of
    'NH07': sea='north_07'
    'NH08': sea='north_08'
    'NH09': sea='north_09'
    'NH10': sea='north_10'
    'NH11': sea='north_11'
    'SH07': sea='south_0708'
    'SH08': sea='south_0809'
    'SH09': sea='south_0910'
    'SH10': sea='south_1011'
endcase

ver='04.20'

pre=['cips_3c_'+sea+'_v'+ver+'_r04_']
pre = replicate(pre,9)


Gs = ['1','2','5']
ng=n_elements(gs)

;***************************** BEGIN LOOP OVER TEST CASES ************************************
FOR TESTCASE = 0,ng-1 DO BEGIN

g = Gs(testcase)

case g of
   '1': fn=4
   '2': fn=5
   '5': fn=6
   '10': fn=7
endcase

SUBTI=season+', v'+VER+', Rev04, '+G+'G'

suf = ['.sav']
suf = replicate(suf,9)

fnames=['01G_all','02G_all','05G_all','10G_all','01G_cld','02G_cld','05G_cld','10G_cld']
fnames=[fnames,'01G_nocld','02G_nocld','05G_nocld','10G_nocld']

fnames=pre+fnames+suf

FILL=-999

if sea eq 'north_07' then hem='NH'
if sea eq 'north_08' then hem='NH'
if sea eq 'north_09' then hem='NH'
if sea eq 'north_10' then hem='NH'
if sea eq 'north_11' then hem='NH'
if sea eq 'south_0708' then hem='SH'
if sea eq 'south_0809' then hem='SH'
if sea eq 'south_0910' then hem='SH'
if sea eq 'south_1011' then hem='SH'

   restore,PTH+fnames(fn)
   lat=1.0*lathi
   nlat=n_elements(lathi)
   for i=0,nlat-1 do lat(i)=(lathi(i)+latlo(i))/2.0
   ndays=max(dfs)-min(dfs)+1
   ddd = indgen(ndays)+min(dfs)
   array=fltarr(ndays,nlat)-99
   for i=0,ndays-1 do begin
      x=where(dfs eq ddd(i),nx)
      if nx gt 0 then begin
         for j=0,nlat-1 do begin
            good=where(iwc(x,j) ne fill,ngood)
            if ngood gt 2 then begin
               array(i,j)=mean(iwc(x(good),j))
            endif
         endfor
      endif
   endfor
;endfor


SET_PLOT,'PS'
DEVICE,/COLOR,BITS_PER_PIXEL=8
DEVICE,/BOLD
DEVICE, FILENAME=PTHo+'IWC_V4.20_'+season+'_'+G+'G.ps'

restore,'c:\idl_cora\c11.tbl'
tvlct,c1,c2,c3

!p.charsize=1.8
!P.FONT=0
!P.THICK=3

!Y.THICK=3
!y.style=1
IF HEM EQ 'NH' THEN !y.range=[60,85] ELSE !Y.RANGE=[-60,-85]
!y.title='Latitude'
if hem eq 'SH' then lat = (-1.0)*lat

!X.THICK=3
!x.style=1
!x.range=[-40,80]
;!x.range=[20,33]
!x.title='Days From Solstice'

FLEV = FINDGEN(11)*15

;MUST USE /CELL_FILL, NOT FILL, WHEN SETTING ARRAY TO 0./0
BAD=WHERE(array EQ -99)
array(BAD)=0./0

col=indgen(11)+1

contour,array(*,0:34),Ddd,lat(0:34),c_colors=col,levels=flev,/cell_fill, $
SUBTITLE=SUBTI,title='Ascending IWC'

POS=[0.95,0.20,1.00,0.775]
COLOR_PANEL_INTEGER,POS,FLEV,COL,-0.00,1.4
XYOUTS,0.975,0.855,'!mm!xg m!u-2!n',/NORMAL,CHARSIZE=1.5,ALIGNMENT=0.5

contour,array(*,35:69),Ddd,lat(35:69),c_colors=col,levels=flev,/cell_fill, $
SUBTITLE=SUBTI,title='Descending IWC'

POS=[0.95,0.20,1.00,0.775]
COLOR_PANEL_INTEGER,POS,FLEV,COL,-0.00,1.4
XYOUTS,0.975,0.855,'!mm!xg m!u-2!n',/NORMAL,CHARSIZE=1.5,ALIGNMENT=0.5

restore,'c:\idl_cora\c11_rb.tbl'
tvlct,c1,c2,c3
COL=[1,2,3,4,5,7,8,9,10,11] ;BLUE < 0, RED > 0
FLEV = [-15,-12,-9,-6,-3,0,3,6,9,12]
diff=ARRAY(*,0:34)-reverse(ARRAY(*,35:69),2)
contour,diff,Ddd,lat(0:34),c_colors=col,levels=flev,/cell_fill, $
title='Asc Minus Desc IWC',SUBTITLE=SUBTI
COLOR_PANEL_INTEGER,POS,FLEV,COL,-0.00,1.4
XYOUTS,0.975,0.855,'!mm!xg m!u-2!n',/NORMAL,CHARSIZE=1.5,ALIGNMENT=0.5

DEVICE,/CLOSE
SET_PLOT,'WIN'

ENDFOR
; ******************************* END LOOP OVER TEST CASES *****************************

end
