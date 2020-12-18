//
//  A. L. Delcher
//
//     File:  Glimmer/get-len.c
//
//  Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Arthur Delcher, Steven Salzberg,
//    Simon Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program reads the fasta file specified on the command line
//  and prints the number of letters in it.


#include  "delcher.h"


#define  INIT_SIZE  10000
#define  INCR_SIZE  10000
#define  MAX_LINE  300


int  Read_String  (FILE *, char * &, long int &, char [], int);



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   char  * Data;
   char  Name [MAX_LINE];
   double  GC_Percent, AT_Fraction, GC_Fraction;
   double  Lambda, Indep_Len;
   long int  i, GC_Count, Input_Size, Len;
   
   if  (argc < 2)
       {
        fprintf (stderr, "USAGE:  %s <fasta-file> <coord-file> \n",
            argv [0]);
        exit (EXIT_FAILURE);
       }

   fp = File_Open (argv [1], "r");

   Data = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Read_String (fp, Data, Input_Size, Name, FALSE);

   fclose (fp);

   Len = strlen (Data + 1);

   printf ("Length of sequence in %s is %ld\n", argv [1], Len);

   GC_Count = 0;
   for  (i = 1;  i <= Len;  i ++)
     if  (Data [i] == 'g' || Data [i] == 'c')
         GC_Count ++;

   GC_Percent = (100.0 * GC_Count) / Len;
   printf ("  AT = %.2f%%  GC = %.2f%%\n", 100.0 - GC_Percent,
           GC_Percent);

   GC_Fraction = GC_Percent / 200.0;
   AT_Fraction = 0.5 - GC_Fraction;
   Lambda = AT_Fraction * AT_Fraction * (2.0 * GC_Fraction + AT_Fraction);
   Indep_Len = 3.0 * log (2.0 * 1e6 * Lambda) / Lambda;

   printf ("  Expected longest orf in 1 million bases = %.0f\n", Indep_Len);

   return  0;
  }



int  Read_String  (FILE * fp, char * & T, long int & Size, char Name [],
                   int Partial)

/* Read next string from  fp  (assuming FASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch == EOF)
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
           break;
         case  's' :
           Ch = 'c';
           break;
         case  'w' :
           if  ((Ct - Lo) % 3 == 1)
               Ch = 'a';
             else
               Ch = 't';
           break;
         case  'r' :
           Ch = 'g';
           break;
         case  'y' :
           Ch = 'c';
           break;
         case  'm' :
           Ch = 'c';
           break;
         case  'k' :
           if  ((Ct - Lo) % 3 == 1)
               Ch = 'g';
             else
               Ch = 't';
           break;
         case  'b' :
           Ch = 'c';
           break;
         case  'd' :
           if  ((Ct - Lo) % 3 == 1)
               Ch = 'a';
             else
               Ch = 't';
           break;
         case  'h' :
           Ch = 'c';
           break;
         case  'v' :
           Ch = 'c';
           break;
         case  'n' :
           Ch = 'c';
           break;
          
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'c';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }



