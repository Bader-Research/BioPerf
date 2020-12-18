//
//  A. L. Delcher
//
//     File:  Glimmer/extract.c
//
//  Version:  1.02  25 Feb 98
//
//  Version:  1.03   9 Aug 2000
//            Ignore lines not beginning with a tag followed by 2 numbers
//
//     Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Arthur Delcher, Steven Salzberg,
//    Simon Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program reads the bases in the first command-line file
//  and then extracts the substrings from it that are
//  specified by the start and end positions in
//  the second command-line file.
//
//  Two command-line options are allowed:
//    -skip  will make the output omit the first 3 characters of
//           each string produced
//    -l n   will make the output omit any string shorter than  n
//           characters (n includes any skipped characters)
//


#include  "delcher.h"


#define  INIT_SIZE  10000
#define  INCR_SIZE  10000
#define  MAX_LINE  300
#define  DEFAULT_MIN_GENE_LEN  1


char  Complement  (char);
int  Read_String  (FILE *, char * &, long int &, char [], int);


long int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;


int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   char  * P;
   char  * Data, Line [MAX_LINE], Name [MAX_LINE];
   long int  i, End, Input_Size, Len, Start, Gene_Len;
   int  Skip_Start_Codon = FALSE;
   
   if  (argc < 3)
       {
        fprintf (stderr, "USAGE:  %s <dna-file> <coord-file> \n",
            argv [0]);
        exit (EXIT_FAILURE);
       }

   for  (i = 3;  i < argc;  i ++)
     if  (strcmp (argv [i], "-skip") == 0)
         Skip_Start_Codon = TRUE;
     else if  (strncmp (argv [i], "-l", 2) == 0)
         {
          if  (strlen (argv [i]) == 2)
              P = argv [++ i];
            else
              P = argv [i] + 2;
          Min_Gene_Len = strtol (P, NULL, 10);
          if  (Min_Gene_Len < 1)
              {
               fprintf (stderr, "ERROR:  Bad minimum length = %s\n", P);
               return  -1;
              }
         }
       else
         {
          fprintf (stderr, "ERROR:  Unrecognized option  %s\n", argv [i]);
          return  -1;
         }

   fp = File_Open (argv [1], "r");

   Data = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Read_String (fp, Data, Input_Size, Name, FALSE);

   fclose (fp);

   Len = strlen (Data + 1);

   fp = File_Open (argv [2], "r");

   while  (fgets (Line, MAX_LINE, fp) != NULL)
     {
      if  (sscanf (Line, "%s %ld %ld", Name, & Start, & End) != 3)
          continue;
      if  ((Start < End && End - Start <= Len / 2) || Start - End > Len / 2)
          {
           Gene_Len = 1 + End - Start;
           if  (Gene_Len < 0)
               Gene_Len += Len;
           if  (Gene_Len < Min_Gene_Len)
               continue;
           printf ("%-10s ", Name);
           if  (Skip_Start_Codon)
               i = Start + 2;
             else
               i = Start - 1;
           do
             {
              i ++;
              if  (i > Len)
                  i -= Len;
              putchar (Data [i]);
             }  while  (i != End);
           putchar ('\n');
          }
        else
          {
           Gene_Len = 1 + Start - End;
           if  (Gene_Len < 0)
               Gene_Len += Len;
           if  (Gene_Len < Min_Gene_Len)
               continue;
           printf ("%-10s ", Name);
           if  (Skip_Start_Codon)
               i = Start - 2;
             else
               i = Start + 1;
           do
             {
              i --;
              if  (i < 1)
                  i += Len;
              putchar (Complement (Data [i]));
             }  while  (i != End);
           putchar ('\n');
          }
        
     }

   return  0;
  }



char  Complement  (char Ch)

/* Returns the DNA complement of  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      case  'r' :          // a or g
        return  'y';
      case  'y' :          // c or t
        return  'r';
      case  's' :          // c or g
        return  's';
      case  'w' :          // a or t
        return  'w';
      case  'm' :          // a or c
        return  'k';
      case  'k' :          // g or t
        return  'm';
      case  'b' :          // c, g or t
        return  'v';
      case  'd' :          // a, g or t
        return  'h';
      case  'h' :          // a, c or t
        return  'd';
      case  'v' :          // a, c or g
        return  'b';
      default :            // anything
        return  'n';
     }
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



