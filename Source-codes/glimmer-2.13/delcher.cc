//  A. L. Delcher
//
//  File:  ANNOTATION/Glimmer/delcher.cc
//  Version:  2.03  3 Dec 2002
//
//  Copyright (c) 2002 by Arthur Delcher, Steven Salzberg, Simon
//  Kasif, and Owen White.  All rights reserved.  Redistribution
//  is not permitted without the express written permission of
//  the authors.
//
//  Common generic routines.
//


#include  "delcher.h"


const char *  Commatize
    (long int  n)

//  Return a string representing the value of  n  with commas
//  every three digits.

  {
   static char  buff [50];
   bool  is_negative = false;
   int  i, ct;

   buff [49] = '\0';

   if  (n == 0)
       {
        buff [48] = '0';
        return  buff + 48;
       }

   i = 48;
   if  (n < 0)
       {
        is_negative = true;
        n *= -1;
       }

   for  (ct = 0;  n > 0;  ct ++)
     {
      if  (ct == 3)
          {
           buff [i --] = ',';
           ct = 0;
          }
      buff [i --] = char ('0' + n % 10);
      n /= 10;
     }

   if  (is_negative)
       buff [i --] = '-';

   return  buff + i + 1;
  }



void  Fasta_Print
    (FILE * fp, char * s, char * hdr)

//  Print string  s  in fasta format to  fp .  Put string  hdr
//  on header line.

  {
   int  ct = 0;

   if  (hdr != NULL)
       fprintf (fp, ">%s\n", hdr);

   while  (* s != '\0')
     {
      if  (ct == FASTA_WIDTH)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (* s, fp);
      s ++;
      ct ++;
     }

   fputc ('\n', fp);

   return;
  }



void  Fasta_Print_N
    (FILE * fp, char * s, int n, char * hdr)

//  Print first  n  bytes of  string  s  in fasta format to  fp .
//  Put string  hdr  on header line.

  {
   int  ct = 0, i;

   if  (hdr != NULL)
       fprintf (fp, ">%s\n", hdr);

   for  (i = 0;  i < n;  i ++)
     {
      if  (ct == FASTA_WIDTH)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (s [i], fp);
      ct ++;
     }

   fputc ('\n', fp);

   return;
  }



FILE *  File_Open  (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }



double  Percent
    (double a, double b)

//  Return  a / b  as a percentage.  Return  0.0  if  b = 0.0 .

  {
   if  (b == 0.0)
       return  0.0;

   return  100.0 * (a / b);
  }



void *  Safe_calloc
    (size_t N, size_t Len, size_t line_num)

/* Allocate and return a pointer to enough memory to hold an
*  array with  N  entries of  Len  bytes each.  All memory is
*  cleared to 0.  Exit if fail. */

  {
   void  * P;

   P = calloc (N, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  calloc failed");
        if  (line_num != 0)
            fprintf (stderr, " at line %lu", (long unsigned) line_num);
        fputc ('\n', stderr);
        fprintf (stderr, "  %lu elements of %lu bytes each\n",
            (long unsigned) N, (long unsigned) Len);
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_malloc
    (size_t Len, size_t line_num)

/* Allocate and return a pointer to  Len  bytes of memory.
*  Exit if fail. */

  {
   void  * P;

   P = malloc (Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  malloc failed");
        if  (line_num != 0)
            fprintf (stderr, " at line %lu", (long unsigned) line_num);
        fputc ('\n', stderr);
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_realloc
    (void * Q, size_t Len, size_t line_num)

/* Reallocate memory for  Q  to  Len  bytes and return a
*  pointer to the new memory.  Exit if fail. */

  {
   void  * P;

   P = realloc (Q, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  realloc failed");
        if  (line_num != 0)
            fprintf (stderr, " at line %lu", (long unsigned) line_num);
        fputc ('\n', stderr);
        exit (EXIT_FAILURE);
       }

   return  P;
  }



