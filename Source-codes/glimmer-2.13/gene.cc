//  A. L. Delcher
//
//     File:  ~adelcher/Glimmer/gene.cc
//  Version:  1.03   8 Feb 99
//            Make easier to modify start/stop codons
//
//    Copyright (c) 1997-99 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  Common routines for DNA sequences.



#include "delcher.h"
#include "gene.h"


unsigned  Ch_Mask
    (char Ch)

/* Returns a bit mask representing character  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  0x1;
      case  'c' :
        return  0x2;
      case  'g' :
        return  0x4;
      case  't' :
        return  0x8;
      case  'r' :     // a or g
        return  0x5;
      case  'y' :     // c or t
        return  0xA;
      case  's' :     // c or g
        return  0x6;
      case  'w' :     // a or t
        return  0x9;
      case  'm' :     // a or c
        return  0x3;
      case  'k' :     // g or t
        return  0xC;
      case  'b' :     // c, g or t
        return  0xE;
      case  'd' :     // a, g or t
        return  0xD;
      case  'h' :     // a, c or t
        return  0xB;
      case  'v' :     // a, c or g
        return  0x7;
      default :       // anything
        return  0xF;
     }
  }



int  Codon_To_Subscript
    (char * s)

//  Return the subscript that corresponds to the first 3 nucleotide
//  letters in string  s .  Return  -1  if any of these positions
//  is not  a, c, g, or t.

  {
   int  i, j, sub = 0;

   for  (i = 0;  i < 3;  i ++)
     {
      sub *= 4;
      j = Nucleotide_To_Subscript (s [i]);
      if  (j < 0)
          return  -1;
        else
          sub += j;
     }

   return  sub;
  }



char  Complement
    (char Ch)

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



char  Filter
    (char Ch)

//  Return a single  a, c, g or t  for  Ch .

  {
   switch  (tolower (Ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  Ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        return  'c';
    }
  }



void  Find_Stop_Codons
    (char X [], int T, int Stop [])

//  Set  Stop [0 .. 6]  TRUE  or  FALSE   according to whether
//  X [1 .. T] has a stop codon in the corresponding reading frame.
//  Stop [6]  is always set  FALSE .

  {
   unsigned  Codon;
   int  i;

   for  (i = 0;  i < 7;  i ++)
     Stop [i] = 0;

   if  (T < 3)
       return;

   Codon = Ch_Mask (X [1]) << 4 | Ch_Mask (X [2]);

   for  (i = 3;  i <= T;  i ++)
     {
      Codon = (Codon & SHIFT_MASK) << 4;
      Codon |= Ch_Mask (X [i]);
      
      if  (Is_Forward_Stop (Codon))
          Stop [i % 3] = TRUE;
      if  (Is_Reverse_Stop (Codon))
          Stop [3 + i % 3] = TRUE;
     }

   return;
  }



int  Is_Forward_Start
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  forward direction.

  {
   return  (
            (Codon & ATG_MASK) == Codon
//               || (Codon & CTG_MASK) == Codon
               || (Codon & GTG_MASK) == Codon
               || (Codon & TTG_MASK) == Codon
           );
  }



int  Is_Forward_Stop
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  forward direction.

  {
   return  (
            (Codon & TAA_MASK) == Codon
               || (Codon & TAG_MASK) == Codon
               || (Codon & TGA_MASK) == Codon
           );
  }



int  Is_Reverse_Start
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  reverse direction.

  {
   return  (
            (Codon & CAT_MASK) == Codon
//               || (Codon & CAG_MASK) == Codon
               || (Codon & CAC_MASK) == Codon
               || (Codon & CAA_MASK) == Codon
           );
  }



int  Is_Reverse_Stop
    (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  reverse direction.

  {
   return  (
            (Codon & TTA_MASK) == Codon
               || (Codon & CTA_MASK) == Codon
               || (Codon & TCA_MASK) == Codon
           );
  }



int  Is_Start
    (char * S)

/* Return  TRUE  iff  S  is a start codon. */

  {
   return  (
            strncmp (S, "atg", 3) == 0
//               || strncmp (S, "ctg", 3) == 0
               || strncmp (S, "gtg", 3) == 0
               || strncmp (S, "ttg", 3) == 0
           );
  }



int  Is_Stop
    (char * S)

/* Return  TRUE  iff  S  is a stop codon. */

  {
   return  (
            strncmp (S, "taa", 3) == 0
               || strncmp (S, "tag", 3) == 0
               || strncmp (S, "tga", 3) == 0
           );
  }



int  Nucleotide_To_Subscript
    (char ch)

//  Return the subscript that corresponds to nucleotide  ch .
//  Return  -1  if  ch  is not a, c, g or t.

  {
   switch  (tolower (ch))
     {
      case  'a' :
        return  0;
      case  'c' :
        return  1;
      case  'g' :
        return  2;
      case  't' :
        return  3;
      default :
        return -1;
     }
  }



int  Read_String
    (FILE * fp, char * & T, long int & Size, char Name [],
     int Partial)

/* Read next string from  fp  (assuming FASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read.
*  Put first string on fasta header line into  Name .
*/

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
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }



int  Rev_Codon_To_Subscript
    (char * s)

//  Return the subscript that corresponds to the reverse complement
//  of the first 3 nucleotide
//  letters in string  s .  Return  -1  if any of these positions
//  is not  a, c, g, or t.

  {
   int  i, j, sub = 0;

   for  (i = 2;  i >= 0;  i --)
     {
      sub *= 4;
      j = Rev_Nucleotide_To_Subscript (s [i]);
      if  (j < 0)
          return  -1;
        else
          sub += j;
     }

   return  sub;
  }



void  Reverse_Complement
  (char S [], long int T)

//  Set  S [1 .. T]  to its DNA reverse complement.

  {
   char  Ch;
   long int  i, j;

   for  (i = 1, j = T;  i < j;  i ++, j --)
     {
      Ch = S [j];
      S [j] = Complement (S [i]);
      S [i] = Complement (Ch);
     }

   if  (i == j)
       S [i] = Complement (S [i]);
  }



int  Rev_Nucleotide_To_Subscript
    (char ch)

//  Return the subscript that corresponds to the reverse complement
//  of nucleotide  ch .
//  Return  -1  if  ch  is not a, c, g or t.

  {
   switch  (tolower (ch))
     {
      case  'a' :
        return  3;
      case  'c' :
        return  2;
      case  'g' :
        return  1;
      case  't' :
        return  0;
      default :
        return -1;
     }
  }



char *  Subscript_To_Codon
    (int sub)

//  Return the nucleotide triple that corresponds to subscript
//   sub .  Return  "bad"  if  sub < 0  or  sub > 63 .

  {
   char  * alphabet = "acgt";
   static char  ans [4];
   int  i;

   if  (sub < 0 || sub > 63)
       {
        strcpy (ans, "bad");
        return  ans;
       }

   ans [3] = '\0';
   for  (i = 2;  i >= 0;  i --)
     {
      ans [i] = alphabet [sub % 4];
      sub /= 4;
     }

   return  ans;
  }



