;This will return three arrays that are calendars.  The first is the regular Gregorian calendar.
;The second has all the orbit numbers that begin that day, and the last one give the day number
;as a running total.
;Input
;  month_number - the number of the month for which you want the calendar.
;Output
;  prints a trio of calendars
Function Calendar, month_number
    ;December causes a problem because the program looks to month 13
	If month_number eq 12 then begin
	   number_of_days=31
	Endif Else begin
	   number_of_days=(ymd2yd([2007,month_number+1,1])-ymd2yd([2007,month_number,1]))[0]
	Endelse
	start_day=(ymd2yd([2007,month_number,1]))[0]
	set_cips_production_vars
	month_day=(indgen(number_of_days+1)+start_day)
	print, month_day
	time=yd2usec(month_day)
	print, time
	;read the data in
	orbit=get_orbit_info(time[0],time[number_of_days])
       orbit_start_time=orbit[*].start_time
       orbit_number=orbit[*].orbit_number
	first_orbit_number=intarr(number_of_days)
    ;This gives the first orbit on that day of the calendar.
	For i=0,number_of_days-1 do begin
       today=time[i]
       today_first_orbit=orbit_number[min(where(orbit_start_time gt today))]
	   first_orbit_number[i]=today_first_orbit
    end

	;Convert January 7th, 2007 (a Sunday) to Julian calendar.
	first_sunday=2454107
	;Convert the first day of the month to Julian calendar.
	day_of_week=yd2jd(start_day)
	;Subtract the two dates and find the remainder.
	d=day_of_week-first_sunday
	day=d mod 7

	;Place this data in an array to look like a calendar.
	;Create a 6x7 array
	calendar=intarr(7,6)
	calendar=calendar-1
	;The day is the position where the calendar starts and the end position is the number of days.
	calendar[day:day+number_of_days-1]=indgen(number_of_days)
	;Add 1 to find the real date
	calendar=calendar+1
	print,calendar

    ;Print the calendar with orbit numbers.
	orbit_calendar=intarr(7,6)
	orbit_calendar[day:day+number_of_days-1]=first_orbit_number
	print, orbit_calendar

    ;Print the calendar with a running total of the days.
	year_calendar=intarr(7,6)
	year_calendar[day:day+number_of_days-1]=indgen(number_of_days)
	year_calendar[day:day+number_of_days-1]=year_calendar[day:day+number_of_days-1]+start_day-2007000
	print, year_calendar
End
