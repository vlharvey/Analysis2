FUNCTION YYYYMMDD_TO_DFS_NH,DATE

;CALL THIS WITH:  RESULT=YYYYMMDD_TO_DFS_NH(DT)
;   WHERE DT IS THE DATE ARRAY.

;CONVERT DATE ARRAY OF FORM YYYYMMDD TO DAYS FROM SOLSTICE
;ASSUMES ARRAY (NOT ONE NUMBER) IS PASSED.
;CALCULATES DAYS FROM 21 JUNE.
;IF DATES IN MULTIPLE YEARS ARE INCLUDED,
;   IT ASSUMES 21 JUNE OF *EACH* YEAR SEPARATELY.
;   DAYS FROM 22 JUNE TO 21 DEC ARE POSITIVE DFS. DAYS FROM
;   22 DEC TO 20 JUN ARE NEGATIVE DFS.

   jday=yyyymmdd_to_jday(date)
   yyyymmdd=strcompress(date,/remove_all)
   yyyy=strmid(yyyymmdd,0,4)
   mmdd=strmid(yyyymmdd,4,4)

   dfs=fix(jday)	;initialize dfs array

   ;Define the correct day of solstice to subtract.
   pos=where(mmdd ge 1222,npos)	;Dec 22 to 31
   if npos gt 0 then begin
      for i=0,npos-1 do begin
         index=pos(i)
         dfs(index) = jday(index)-julday(6,21,yyyy(index)+1)
      endfor
   endif
   pos=where(mmdd ge 101 and mmdd le 621,npos)	;Jan 1 to June 21
   if npos gt 0 then begin
      for i=0,npos-1 do begin
         index=pos(i)
         dfs(index)=jday(index)-julday(6,21,yyyy(index))
      endfor
   endif
   pos=where(mmdd gt 621 and mmdd lt 1222,npos)	;June 22 to Dec 21
   if npos gt 0 then begin
      for i=0,npos-1 do begin
         index=pos(i)
         dfs(index)=jday(index)-julday(6,21,yyyy(index))
      endfor
   endif


RETURN,DFS

END