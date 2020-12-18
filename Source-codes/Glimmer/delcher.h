//  A. L. Delcher
//
//  File:  ANNOTATION/Glimmer/delcher.h
//  Version:  2.03  3 Dec 2002
//
//  Copyright (c) 2002 by Arthur Delcher, Steven Salzberg, Simon
//  Kasif, and Owen White.  All rights reserved.  Redistribution
//  is not permitted without the express written permission of
//  the authors.
//
//  Common generic routines.
//


#ifndef  __DELCHER_H_INCLUDED
#define  __DELCHER_H_INCLUDED


#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream.h>
#include  <iomanip.h>
#include  <fstream.h>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <errno.h>
#include  <unistd.h>


#define  TRUE  1
#define  FALSE  0
#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif


const int  FASTA_WIDTH = 60;
  // Max number of characters to print on a FASTA data line


const char *  Commatize
    (long int  n);
void  Fasta_Print
    (FILE * fp, char * s, char * hdr);
void  Fasta_Print_N
    (FILE * fp, char * s, int n, char * hdr);
FILE *  File_Open
    (const char *, const char *);
double  Percent
    (double a, double b);
void *  Safe_calloc
    (size_t N, size_t Len, size_t line_num = 0);
void *  Safe_malloc
    (size_t Len, size_t line_num = 0);
void *  Safe_realloc
    (void * Q, size_t Len, size_t line_num = 0);


template <class DT>
  DT  Max  (DT, DT);
template <class DT>
  DT  Min  (DT, DT);
template <class DT>
  void  Swap  (DT &, DT &);


template <class DT>
DT  Max  (DT A, DT B)

//  Return the larger of  A  and  B .

  {
   if  (A > B)
       return  A;
     else
       return  B;
  }



template <class DT>
DT  Min  (DT A, DT B)

//  Return the smaller of  A  and  B .

  {
   if  (A < B)
       return  A;
     else
       return  B;
  }



template <class DT>
void  Swap  (DT & A, DT & B)

//  Swap the values in  A  and  B .

  {
   DT  Save;

   Save = A;
   A = B;
   B = Save;

   return;
  }



#endif
