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
#include "miscutil.h"
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void bubbleSort(float **array, int nArray) {
	int notSorted = 1;
	float tmp;
	while (notSorted) {
		notSorted = 0;
		for (int n = 0; n < nArray - 1; n++) {
			if (*(*array + n) > *(*array + n + 1)) {
				notSorted = 1;
				tmp = *(*array + n);
				*(*array + n) = *(*array + n + 1);
				*(*array + n + 1) = tmp;
			}
		}
	}
}

void bubbleSort(double **array, int nArray) {
	int notSorted = 1;
	double tmp;
	while (notSorted) {
		notSorted = 0;
		for (int n = 0; n < nArray - 1; n++) {
			if (*(*array + n) > *(*array + n + 1)) {
				notSorted = 1;
				tmp = *(*array + n);
				*(*array + n) = *(*array + n + 1);
				*(*array + n + 1) = tmp;
			}
		}
	}
}

void bubbleSort(int **array, int nArray) {
	int tmp, notSorted = 1;
	while (notSorted) {
		notSorted = 0;
		for (int n = 0; n < nArray - 1; n++) {
			if (*(*array + n) > *(*array + n + 1)) {
				notSorted = 1;
				tmp = *(*array + n);
				*(*array + n) = *(*array + n + 1);
				*(*array + n + 1) = tmp;
			}
		}
	}
}

void andArrays(int *array1, int nA1, int *array2, int nA2, int **result,
	       int &nResult) {
	
	if (nA1 > nA2) {
		int tmp = nA1; nA1 = nA2; nA2 = tmp;
		int *tmp2 = array1; array1 = array2; array2 = tmp2;
	}

	int n2 = 0;
	(*result) = new int[nA1]; nResult = 0;

	for (int n = 0; n < nA1; n++) {
		while (array1[n] > array2[n2] && n2 < nA2) n2++;

		if (n2 >= nA2) break;

		if (array1[n] == array2[n2]) {
			*(*result + nResult) = array1[n];
			nResult++;
		}
	}
}

int sortedArrayContains(int *array, int number, int nArray) {
	if (array[nArray-1] < number) return 0;
	
	int high = nArray - 1, low = 0, mid;


	while (low < high && (low + 1) != high) {
		if (number == array[high] || number == array[low]) return 1;
		mid = low + (int) ((high - low)/2);

		if (number == array[mid]) return 1;

		if (number < array[mid]) {
			high = mid;
		} else {
			low = mid;
		}
	}

	if (number == array[low]) return 1;

	return 0;
}

void stddev(int *values, int num, double *ave, double *stdDev) {
	*stdDev = *ave = 0.;
        if (num) {
                double ssq = 0.0, diff;

                *ave = average(values, num);
 
                for (int i = 0; i < num; i++) {
                        diff = values[i] - *ave;
                        *stdDev += diff*diff;
                }
 
                *stdDev /= num;
		*stdDev = sqrt(*stdDev);
        }
}

void stddev(double *values, int num, double *ave, double *stdDev) {
	*stdDev = *ave = 0.;
        if (num) {
                double ssq = 0.0, diff;

                *ave = average(values, num);
 
                for (int i = 0; i < num; i++) {
                        diff = values[i] - *ave;
                        *stdDev += diff*diff;
                }
 
                *stdDev /= num;
		*stdDev = sqrt(*stdDev);
        }
}

void stddev(float *values, int num, double *ave, double *stdDev) {
	*stdDev = *ave = 0.;
        if (num) {
                double ssq = 0.0, diff;

                *ave = average(values, num);
 
                for (int i = 0; i < num; i++) {
                        diff = (double) (values[i] - *ave);
                        *stdDev += diff*diff;
                }
 
                *stdDev /= num;
		*stdDev = sqrt(*stdDev);
        }
}

double average(double *value, int count) {
	double total = 0.;
	for (int i = 0; i < count; i++) total += value[i];
	return (double) total/count;
}

double average(float *value, int count) {
	float total = 0.;
	for (int i = 0; i < count; i++) total += value[i];
	return (double) total/count;
}

double average(int *value, int count) {
	int total = 0;
	for (int i = 0; i < count; i++) total += value[i];
	return (double) total/count;
}

float leftOver(float nom, float denom) {
	if (denom) {
		return ((nom/denom) - (int) (nom/denom)) * denom;
	} else {
		return 0.;
	}
}

int leftOver(int nominator, int denominator) {
	if (denominator) {
		return nominator%denominator;
	} else {
		return 0;
	}
}

const char *itoa(int num) {
	static char *out = new char[255];

	if (!num) return "0";

	out[0] = '\0';
	int pos = 0;

	if (num < 0) {
		num *= -1;
		out[pos++] = '-';
	}

	int nDigs = (int) log10((double) num);
	int digit, tenPow;

	for (int i = nDigs; i > 0; i--) {
		tenPow = (int) pow(10., (double) i);
		digit = (int) (num/tenPow);
		out[pos++] = digit + 48;
		num -= tenPow*digit;
	}

	out[pos++] = num + 48;
	out[pos] = '\0';
	return out;
}

int file_cmp(char *file_test, char *file_tmpl)
{
  char *star=strchr(file_tmpl, '*');
  if(!star) return(strcmp(file_test, file_tmpl));
  if(strncmp(file_test, file_tmpl, strlen(file_tmpl)-strlen(star))) return(1);
  return(strcmp(&file_test[strlen(file_test)-strlen(star)+1], star+1));
}

void compressText(char ** text)
{
  int i, space=0, cLen=0;
  if(!(*text)) return;
  if(!strlen(*text)) return;
  for(i=0; i<strlen(*text); i++)
    {
      if(*(*text+i)==' ' && space==1) continue;
      if(*(*text+i)==' ' && space==0) space=1;
      if(*(*text+i)!=' ') space=0;
      cLen++;
    }

  char * text_=new char [cLen+1];
  text_[cLen]='\0';

  cLen=0;
  for(i=0; i<strlen(*text); i++)
    {
      if(*(*text+i)==' ' && space==1) continue;
      if(*(*text+i)==' ' && space==0) space=1;
      if(*(*text+i)!=' ') space=0;
      text_[cLen]=*(*text+i);
      cLen++;
    }

  if (text_[cLen-1] == ' ') {
	text_[cLen-1] = '\0'; // delete [] (text_+cLen);
	// char *tmp = new char[cLen]; tmp[cLen-1] = '\0';
	// strncpy(tmp, text_, cLen-1);
	// delete [] text_;
	// text_ = tmp;
  }

  delete [] (*text);
  *text=text_;
}

double ranBase_ = 1.0;
double ran_() {
  static double ranNew, const1=470001.0, const2=999563.0;
  ranBase_=fmod(const1*ranBase_, const2);
  ranNew=ranBase_/const2;
  return(ranNew);
}

void AddBuffer(char *** array,
               char * buffer,
               int Narray) {
 
        char ** array_;
        array_=new char* [Narray+1];
        for(int i=0; i<Narray; i++) array_[i]=*(*array+i);
        array_[Narray]=new char[strlen(buffer)+1];
        strcpy(array_[Narray], buffer);
        if(Narray) delete [] (*array);
        (*array)=array_;
}

void setText(char **text, char *buffer) {
  *text=new char [(buffer?strlen(buffer)+1:1)];
  if(buffer) strcpy(*text, buffer);
}

void addText(char ** text, char * buffer) {
        char * text_=new char [(*text?strlen(*text):0)+
                         (buffer?strlen(buffer)+1:1)];
        if(*text) {
                if(*text) strcpy(text_, *text);
                if(buffer) strcat(text_, buffer);
                delete [] (*text);
                *text=text_;
        } else {
                if(buffer) strcpy(text_, buffer);
                *text=text_;
        }
}

void addChar(char ** text, char c) {
	int l = (*text?strlen(*text):0);
        char * text_=new char [l+2];

        if(*text) {
                strcpy(text_, *text);
                delete [] (*text);
        }

	text_[l] = c; text_[l+1] = '\0';
	*text=text_;
}

void StringArrayToBuffer(int nItems, char **array, char ***names) {
  char bar[2]; bar[0]=1; bar[1]='\0';
  *array=NULL;
  for(int i=0; i<nItems; i++) {
    addText(array, (*names)[i]);
    addText(array, bar);
  }
}

void BufferToStringArray(int nItems, char **array, char ***names) {
  int nItems_, i, j;
  *names= new char* [nItems];
 
  for(i=0, j=0, nItems_=0; i<strlen(*array); i++) {
    if((*array)[i]==1) {
      (*names)[nItems_]=new char [i-j+1];
      if(i>j)
	strncpy((*names)[nItems_], *array+j, i-j);
      *((*names)[nItems_]+i-j)='\0';
      j=i+1; nItems_++;
    }
  }
  delete [] *array;
}

char *getPath(int argc, char **argv) {
	int c;
	char *dbpath;
	int i;

	dbpath = NULL;

	for (i = 0; i < argc; i++) {
		if (!strcmp("-D", argv[i]) && i < argc - 1) {
			dbpath = argv[i+1];
			break;
		}
	}

	if (!dbpath || strlen(dbpath) == 0) dbpath = DBPATH;

	return dbpath;
}
