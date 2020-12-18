/*
 
Copyright (c)  1995-2000   The Regents of the University of California
All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this PDB
software for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear
in all copies.
 
Those desiring to incorporate this PDB Software into commercial products
or use for commercial purposes should contact the Technology Transfer
Office, University of California, San Diego, 9500 Gilman Drive, La Jolla,
CA 92093-0910, Ph: (619) 534-5815, FAX: (619) 534-7345.
 
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS PDB SOFTWARE, EVEN IF THE
UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
 
THE PDB SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
THE USE OF THE PDB SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.
 
*/
/* $Revision: 2.3 $ $Date: 2000/01/05 18:00:30 $ */
#include "jdate.h"
#include <string.h>
#include <stdlib.h>

const Weekdays JDate::weekdays[] = {
	"Sun", "Sunday",
	"Mon", "Monday",
	"Tue", "Tuesday",
	"Wed", "Wednesday",
	"Thu", "Thursday",
	"Fri", "Friday",
	"Sat", "Saturday"
};

const Months JDate::months[] = {
	"Jan", "January",
	"Feb", "February",
	"Mar", "March",
	"Apr", "April",
	"May", "May",
	"Jun", "June",
	"Jul", "July",
	"Aug", "August",
	"Sep", "September",
	"Oct", "October",
	"Nov", "November",
	"Dec", "December"
};

int JDate::jday(int month, int day, int year) {

	// my second Y2K bug
	if (year < 100 && year >= 72) year += 1900;
	// this will be a Y2072 bug
	else if (year < 72) year += 2000;

	if (month > 2) {
		month -= 3;
	} else {
		month += 9;
		year--;
	}

	int century = abs(year/100);
	int anni = year - 100*century;

	int jd = abs((146097 * century)/4)
		 + abs((1461 * anni)/4)
		 + abs((153 * month + 2)/5)
		 + day + 1721119;

	return jd;
}

int JDate::jday(char *month, int day, int year) {

	for (int y = 0; y < 12; y++)
		if (!strncasecmp(months[y].shortName, month, 3))
			return jday(y+1, day, year);
	
	return -1;
}

int *JDate::jdate(int julian) {

	int *result = new int[4];

	int weekday = (julian + 1)%7;

	int jtmp = julian - 1721119;
	int year = (int) (4 * jtmp - 1)/146097;

	jtmp = 4 * jtmp - 1 - 146097 * year;
	int day = (int) jtmp/4;
	jtmp = (int) (4 * day + 3)/1461;
	day = 4 * day + 3 - 1461 * jtmp;
	day = (int) (day + 4)/4;

	int month = (int) (5 * day - 3)/153;
	day = 5 * day - 3 - 153 * month;
	day = (int) (day + 5)/5;

	year = 100 * year + jtmp;

	if (month < 10) {
		month += 3;
	} else {
		month -= 9;
		year++;
	}

	result[0] = month;
	result[1] = day;
	result[2] = year;
	result[3] = weekday;

	return result;
}
