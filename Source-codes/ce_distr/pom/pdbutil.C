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
/* $Revision: 2.1 $ $Date: 2000/01/21 20:38:25 $ */
#include "pdbutil.h"
#include "jdate.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int toJday(char *date) {
	int day, year;
	char *month = new char[4];

	char *myDate = new char[strlen(date)+1];
	strcpy(myDate, date);

	myDate[2] = '\0';
	myDate[6] = '\0';

	day = atoi(myDate);
	strcpy(month, &myDate[3]);
	year = atoi(&myDate[7]);

	int result = JDate::jday(month, day, year);

	delete [] myDate;
	delete [] month;

	return result;
}

char *fromJday(int jday) {
	if (!jday) return "???";
	int *result = new int[4];
	result = JDate::jdate(jday);
	result[2] -= 1900;

	// fixing a Y2K with a Y3K problem!!
	if (result[2] >= 100) result[2] -= 100;

	char *date = new char[9];

	sprintf(date, "%2.2d-%s-%2.2d", result[1],
				    JDate::months[result[0]-1].shortName,
				    result[2]);

	date[4] -= 32; date[5] -= 32;

	delete [] result;

	return date;
}

int isId(char *id) {
	if (!id) return 0;
	int l = strlen(id);
	if (l != 4) return 0;
	if (id[0] < 48 || id[0] > 57) return 0;
	for (int n = 1; n < l; n++) {
		if (id[n] < 48 ||
		    (id[n] > 57 && id[n] < 65) ||
		    (id[n] > 91 && id[n] < 97) ||
		    id[n] > 122) return 0;
	}

	return 1;
}
