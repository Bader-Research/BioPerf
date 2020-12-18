//  Programmer:  A. Delcher
//
//        File:  Glimmer/anomaly.cc
//
//     Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Arthur Delcher, Steven Salzberg,
//    Simon Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program reads the bases in the first command-line file
//  and then takes the start and end positions specified in
//  the second command-line file and checks for anomalous
//  start/stop codons and frame shifts.
//


#include  "delcher.h"
#include  "gene.h"


#define  INIT_SIZE  10000
#define  INCR_SIZE  10000
#define  MAX_LINE  3000



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   char  * Buffer, * Data, Line [MAX_LINE], Name [MAX_LINE];
   char  Codon [4] = "tag";
   int  Frame_Shift, Check_Start_Codon = TRUE;
   int  Check_Previous_Stop = FALSE, Direction;
   long int  Buffer_Len, Gene_Len;
   long int  i, j, Begin, End, Input_Size, Len, Start, Stop;
   
   if  (argc < 3)
       {
        fprintf (stderr, "USAGE:  %s <dna-file> <coord-file> \n",
            argv [0]);
        exit (EXIT_FAILURE);
       }

   fp = File_Open (argv [1], "r");

   for  (i = 3;  i < argc;  i ++)
     if  (strcmp (argv [i], "-start") == 0)
         Check_Start_Codon = FALSE;
     else if  (strcmp (argv [i], "+orf") == 0)
         Check_Previous_Stop = TRUE;
       else
         {
          fprintf (stderr, "ERROR:  Unrecognized option  %s\n", argv [i]);
          return  -1;
         }
     
   Buffer = (char *) Safe_malloc (INIT_SIZE);
   Buffer_Len = INIT_SIZE;

   Data = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Read_String (fp, Data, Input_Size, Name, FALSE);

   fclose (fp);

   Len = strlen (Data + 1);

   fp = File_Open (argv [2], "r");

   while  (fgets (Line, MAX_LINE, fp) != NULL)
     {
      sscanf (Line, "%s %ld %ld", Name, & Start, & End);
      if  (Start < End && End - Start <= Len / 2 || Start - End > Len / 2)
          {
           Direction = +1;
           Gene_Len = 4 + End - Start;
           if  (Gene_Len < 0)
               Gene_Len += Len;

           if  (Buffer_Len < Gene_Len + 1)
               Buffer = (char *) Safe_realloc (Buffer, 1 + Gene_Len);
           Buffer_Len = 1 + Gene_Len;
           for  (i = 0;  i < Gene_Len;  i ++)
             {
              if  (Start + i <= Len)
                  j = Start + i;
                else
                  j = Start + i - Len;
              Buffer [i] = tolower (Data [j]);
             }
           Buffer [i] = '\0';
          }
        else
          {
           Direction = -1;
           Gene_Len = 4 + Start - End;
           if  (Gene_Len < 0)
               Gene_Len += Len;

           if  (Buffer_Len < Gene_Len + 1)
               Buffer = (char *) Safe_realloc (Buffer, 1 + Gene_Len);
           Buffer_Len = 1 + Gene_Len;
           for  (i = 0;  i < Gene_Len;  i ++)
             {
              if  (Start - i >= 1)
                  j = Start - i;
                else
                  j = Start - i + Len;
              Buffer [i] = Complement (tolower (Data [j]));
             }
           Buffer [i] = '\0';
          }

      if  (Check_Previous_Stop)
          {
           if  (Direction == +1)
               {
                for  (i = 3;  i > 0;  i --)
                  if  (Start - i < 1)
                      Codon [i] = tolower (Data [Start - i + Len]);
                    else
                      Codon [i] = tolower (Data [Start - i]);
               }
             else
               {
                for  (i = 3;  i > 0;  i --)
                  if  (Start + i > Len)
                      Codon [i] = Complement (tolower (Data [Start + i - Len]));
                    else
                      Codon [i] = Complement (tolower (Data [Start + i]));
               }
           if  (strncmp (Codon, "taa", 3) != 0
                   && strncmp (Codon, "tag", 3) != 0
                   && strncmp (Codon, "tga", 3) != 0)
               printf ("%-10s %8ld %8ld no stop before start\n",
                              Name, Start, End);
          }
      if  (Check_Start_Codon
             && strncmp (Buffer, "atg", 3) != 0
             && strncmp (Buffer, "ctg", 3) != 0
             && strncmp (Buffer, "gtg", 3) != 0
             && strncmp (Buffer, "ttg", 3) != 0
          )
          printf ("%-10s has bad start codon = %.3s\n", Name, Buffer);
      if  (strcmp (Buffer + Gene_Len - 3, "taa") != 0
             && strcmp (Buffer + Gene_Len - 3, "tag") != 0
             && strcmp (Buffer + Gene_Len - 3, "tga") != 0)
          {
           printf ("%-10s has bad stop codon = %s\n", Name, Buffer + Gene_Len - 3);
           for  (j = Gene_Len;  j < Len;  j += 3)
             {
              for  (i = 0;  i < 3;  i ++)
                if  (Direction == +1)
                    {
                     if  (Start + i + j > Len)
                         Codon [i] = tolower (Data [Start + i + j - Len]);
                       else
                         Codon [i] = tolower (Data [Start + i + j]);
                    }
                  else
                    {
                     if  (Start - i - j < 1)
                         Codon [i] = Complement (tolower (Data [Start - i - j + Len]));
                       else
                         Codon [i] = Complement (tolower (Data [Start - i - j]));
                    }
              if  (strncmp (Codon, "taa", 3) == 0
                     || strncmp (Codon, "tag", 3) == 0
                     || strncmp (Codon, "tga", 3) == 0)
                  break;
             }
           assert (j < Len);
           printf ("           next stop occurs at offset %ld   Gene_Len = %ld\n",
                       j, Gene_Len);
          }

      Frame_Shift = (Gene_Len % 3);
      if  (Frame_Shift)
          {
           printf ("%-10s %8ld %8ld has %+d frame shift\n",
                           Name, Start, End, Frame_Shift);

           for  (i = 0;  i < Gene_Len - 3;  i += 3)
             if  (strncmp (Buffer + i, "taa", 3) == 0
                     || strncmp (Buffer + i, "tag", 3) == 0
                     || strncmp (Buffer + i, "tga", 3) == 0)
                 break;
           if  (i < Gene_Len - 3)
               {
                Stop = Start + Direction * (i - 1);
                if  (Stop < 1)
                    Stop += Len;
                else if  (Stop > Len)
                    Stop -= Len;
                printf ("   Best prefix is %8ld %8ld   Len = %ld\n",
                             Start, Stop, i);
               }
             else
               {
                printf ("   No stop found in start frame\n");
                continue;
               }

           for  (i = Gene_Len - 6;  i >= 0;  i -= 3)
             if  (strncmp (Buffer + i, "taa", 3) == 0
                     || strncmp (Buffer + i, "tag", 3) == 0
                     || strncmp (Buffer + i, "tga", 3) == 0)
                 break;
           i += 3;
           Begin = Start + Direction * i;
           if  (Begin < 1)
               Begin += Len;
           else if  (Stop > Len)
               Begin -= Len;
           printf ("   Best suffix is %8ld %8ld   Len = %ld\n",
                        Begin, End, Gene_Len - i - 3);

          }
        else
          {
           for  (i = 0;  i < Gene_Len - 3;  i += 3)
             if  (strncmp (Buffer + i, "taa", 3) == 0
                     || strncmp (Buffer + i, "tag", 3) == 0
                     || strncmp (Buffer + i, "tga", 3) == 0)
                 printf ("%-10s has stop codon %.3s at offset %ld\n",
                                 Name, Buffer + i, i);
          }
     }

   return  0;
  }



